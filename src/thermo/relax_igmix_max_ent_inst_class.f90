!> Maximum-entropy sequential relaxation (relaxation_ent_max.tex): SG/NASG liquid + vapor/air ideal-gas mixture.
!> Transfers x~=[V~,E~,m~] (volume received, heat and mass given by the liquid) move x_l=[V_l,E_l,m_l] along D=[d_M,d_T,d_C].
!> Newton on F^=[p_l-p_g-p_j,T_l-T_g,g_l-g_v(p_v,T_g)] with J^=(Z_l+Z_g*A)*D, eliminated M->T->C (eqs. newton_rows-row_C);
!> relaxed states by the iteration of the stages (EOS update, D held at the start, damped steps; tex: Iteration of the stages).
module relax_igmix_max_ent_inst_class
   use precision,            only: WP
   use relax_igmix_sg_class, only: relax_igmix_sg,Prelax,PTrelax,PTgrelax
   use stiffened_gas_class,  only: stiffened_gas
   use igmix_class,          only: igmix
   use thermorelax_class,    only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_DEGENERATE,RELAX_SINGULAR,RELAX_NUC_FAILED
   implicit none
   private

   public :: relax_igmix_max_ent_inst,phasic_state,max_ent_apply
   public :: pI_liquid,pI_gas,pI_mean,pI_impedance,pI_end

   !> Interface pressure p_I of d_M (eq. pI_w): w_M=0 liquid, 1 gas, 1/2 mean, or Pelanti's impedance weighting
   !> pI_end: p_I=p*, the end pressure of the instantaneous mechanical stage (tex: Choice of p_I)
   integer, parameter :: pI_liquid=0,pI_gas=1,pI_mean=2,pI_impedance=3,pI_end=4

   !> Cell state (cell volume 1): extensive x_k, intensive z_k=[p,T,g], forces F^ and Z_k=dz_k/dx_k
   type :: phasic_state
      real(WP) :: Vl,El,ml,Vg,Eg,mv,ma
      real(WP) :: pl,Tl,gl,sl,cl
      real(WP) :: pg,Tg,gv,sv,xv,cg
      logical  :: dry_gas
      real(WP), dimension(3)   :: F
      real(WP), dimension(3,3) :: Zl,Zg
   end type phasic_state

   type, extends(relax_igmix_sg) :: relax_igmix_max_ent_inst
      !> Interface-pressure closure of d_M (eq. pI_w). Default pI_end, p_I=p*: the only constant p_I whose end state is the pressure
      !> equilibrium and the maximum of S on its path (tex: Choice of p_I), by Newton on eqs. mech_end (NR_stage_pI_end)
      !> Other modes hold p_I at the start state: of these only the impedance weighting keeps K^M>0 except at extreme forces
      !> (tex Algorithm step 2), and the end state is not the maximum of S on its path (eq. FM_end), so p-only does not
      !> guarantee dS>=0 when T_l and T_g differ. A finite mechanical step uses the impedance p_I for pI_end; pT, pTg do not depend on p_I
      integer  :: pI_mode=pI_end
      !> Stage tolerances: |F_M|<FM_tol*p, |F_T|<FT_tol*T, |F_C|<FC_tol*Rv*T (floored at round-off)
      real(WP) :: FM_tol=1.0e-10_WP,FT_tol=1.0e-10_WP,FC_tol=1.0e-10_WP
      integer  :: itmax=50
      !> Maximum halvings of a damped step until both phases have a state
      integer  :: itmax_half=60
      !> Vapor mole fraction seeded into a dry gas before the chemical stage
      real(WP) :: xv_dry_seed=1.0e-10_WP
      !> Diagnostics: iterations of the last stages, and pivots K^M,K^T,K^C at their start (eq. stiff_def)
      integer  :: it_M=0,it_T=0,it_C=0
      real(WP) :: KM=0.0_WP,KT=0.0_WP,KC=0.0_WP
   contains
      procedure :: initialize
      procedure :: apply=>max_ent_apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: require_unclamped
      procedure :: admissible
      procedure :: iVQ
      procedure :: get_state
      procedure :: liq_state
      procedure :: gas_state
      procedure :: get_dir
      procedure :: get_jac
      procedure :: eliminate
      procedure :: transfer
      procedure :: force_tol
      procedure :: converged
      procedure :: dampedNR
      procedure :: NR_stage_pI_end
      procedure :: NR_stage_block
      procedure :: relax_M
      procedure :: relax_MT
      procedure :: relax_C
      procedure :: liq_v_from_p_T
      procedure :: nucleate
      procedure :: get_entropy
   end type relax_igmix_max_ent_inst

contains

   !> Parent initialization (EOS pointers, saturation coefficients), then the NASG saturation coefficient
   subroutine initialize(this,liq,gas,indV,indA)
      use nasg_class, only: nasg
      implicit none
      class(relax_igmix_max_ent_inst),    intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      call this%relax_igmix_sg%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      select type (liq)
      type is (nasg)
         this%ES=liq%b/gas%R(indV)
      end select
   end subroutine initialize

   !> Index of the vapor partial density in Q
   integer function iVQ(this)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      iVQ=7+this%liq%ns+this%indV-1
   end function iVQ

   !> Model dispatch on mixture cells; pTg first seeds the opposite phase in metastable near-pure cells.
   !> A cell that is not relaxed is returned untouched with a non-OK code.
   !> Exception: a pure cell that pTg does not seed is already the maximum of S and returns RELAX_OK.
   subroutine max_ent_apply(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      use messager,     only: die
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      real(WP), parameter :: VF_nuc=0.01_WP
      real(WP), dimension(size(Q)) :: Q0
      real(WP) :: VF0
      logical  :: mixture
      integer  :: ier
      this%nucleated=.false.
      this%CRR_loc=0.0_WP; this%tauC_loc=0.0_WP
      this%it_M=0; this%it_T=0; this%it_C=0
      this%KM=0.0_WP; this%KT=0.0_WP; this%KC=0.0_WP
      VF0=VF; Q0=Q
      mixture=(VF.ge.VFlo.and.VF.le.VFhi.and.Q(1).gt.0.0_WP.and.Q(2).gt.0.0_WP)
      select case (this%model)
      case (Prelax,PTrelax)
         if (.not.mixture) then
            if (present(ierr)) ierr=RELAX_DEGENERATE
            return
         end if
         ! Only a start state of the EOS has an entropy to maximize
         call this%admissible(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) then
            if (this%model.eq.Prelax) then
               call this%p_relax(dt,VF,Q,Pjump,ier)
            else
               call this%pT_relax(dt,VF,Q,Pjump,ier)
            end if
         end if
      case (PTgrelax)
         ! Near-pure cells are relaxed here only if seeded; the others are returned unchanged for the split-remix reactor
         if (VF.lt.VF_nuc.or.VF.gt.1.0_WP-VF_nuc) then
            if (this%do_nucleate) call this%nucleate(VF,Q,Pjump,this%nucleated)
            if (.not.this%nucleated) then
               ! A pure cell that is not seeded is the maximum of S (eq. nuc_cond)
               ier=RELAX_DEGENERATE
               if (VF.lt.VFlo.or.VF.gt.VFhi) ier=RELAX_OK
               if (present(ierr)) ierr=ier
               return
            end if
         else if (.not.mixture) then
            if (present(ierr)) ierr=RELAX_DEGENERATE
            return
         end if
         call this%admissible(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%pTg_relax(dt,VF,Q,Pjump,ier)
      case default;    call die('[relax_igmix_max_ent_inst apply] model must be Prelax, PTrelax or PTgrelax')
      end select
      if (ier.ne.RELAX_OK) then
         VF=VF0; Q=Q0
         if (this%nucleated) then
            this%nucleated=.false.
            ier=RELAX_NUC_FAILED
         end if
      end if
      if (present(ierr)) ierr=ier
   end subroutine max_ent_apply

   !> Mechanical stage: maximum of S along x^0+span{d_M}
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer :: ier
      call this%relax_M(VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine p_relax

   !> Direct iteration of the thermo-mechanical stage: maximum of S at fixed masses (pT-relaxed state)
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer :: ier
      call this%relax_MT(VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine pT_relax

   !> Maximum of S over all exchanges (pTg-relaxed state, F=0), tex: Iteration of the stages, Order: the pT-relaxed state at
   !> the start masses, then the chemical stage on the pT-relaxed states
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer  :: ier
      call this%require_unclamped()
      call this%relax_MT(VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) call this%relax_C(VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine pTg_relax

   !> pTg works on the NASG liquid without the brhomax clamp of nasg_class: die unless brhomax>=1
   subroutine require_unclamped(this)
      use messager,   only: die
      use nasg_class, only: nasg
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      select type (liq=>this%liq)
      type is (nasg)
         if (liq%brhomax.lt.1.0_WP) call die('[relax_igmix_max_ent pTg_relax] requires an unclamped NASG liquid: set brhomax=1')
      end select
   end subroutine require_unclamped

   !> Start state of a relaxation call: both phases have a state of their EOS (get_state)
   subroutine admissible(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)  :: this
      real(WP),                        intent(in)  :: VF
      real(WP), dimension(:),          intent(in)  :: Q
      real(WP),                        intent(in)  :: Pjump
      integer,                         intent(out) :: ierr
      type(phasic_state) :: st
      call this%get_state(VF,Q,Pjump,st,ierr)
   end subroutine admissible

   !> Phasic states, forces F^ (eq. Fhat_def) and Z_k=dz_k/dx_k from (VF,Q); x_l=[V_l,E_l,m_l], x_g=[V_g,E_g,m_v]
   subroutine get_state(this,VF,Q,Pjump,st,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)  :: this
      real(WP),                   intent(in)  :: VF
      real(WP), dimension(:),     intent(in)  :: Q
      real(WP),                   intent(in)  :: Pjump
      type(phasic_state),              intent(out) :: st
      integer,                    intent(out) :: ierr
      ierr=RELAX_DEGENERATE
      st%Vl=VF; st%Vg=1.0_WP-VF
      st%ml=Q(1); st%El=Q(3); st%Eg=Q(4)
      st%mv=min(Q(this%iVQ()),Q(2)); st%ma=Q(2)-st%mv
      if (st%Vl.le.0.0_WP.or.st%Vg.le.0.0_WP.or.st%ml.le.0.0_WP.or.Q(2).le.0.0_WP.or.st%mv.lt.0.0_WP) return
      call this%liq_state(st,ierr)
      if (ierr.ne.RELAX_OK) return
      call this%gas_state(st,ierr)
      if (ierr.ne.RELAX_OK) return
      ! Forces (eq. Fhat_def)
      st%F=[st%pl-st%pg-Pjump,st%Tl-st%Tg,st%gl-st%gv]
      if (st%dry_gas) st%F(3)=0.0_WP
   end subroutine get_state

   !> Liquid p,T,s,g,c and Z_l from st%Vl,El,ml through the material interface
   subroutine liq_state(this,st,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)    :: this
      type(phasic_state),              intent(inout) :: st
      integer,                    intent(out)   :: ierr
      real(WP) :: rho,el,Gam,Ks,Cv
      ierr=RELAX_BAD_LIQUID
      rho=st%ml/st%Vl; el=st%El/st%ml
      ! No state past the EOS density limit (NASG co-volume); a finite test, valid under -ffast-math
      if (rho.ge.this%liq%get_rho_max(y=[1.0_WP])) return
      st%pl=this%liq%get_p_from_rho_e(rho=rho,e=el,y=[1.0_WP])
      st%Tl=this%liq%get_T_from_rho_e(rho=rho,e=el,y=[1.0_WP])
      st%cl=this%liq%get_c_from_rho_e(rho=rho,e=el,y=[1.0_WP])
      if (.not.(st%Tl.gt.0.0_WP.and.st%cl.gt.0.0_WP)) return
      st%sl=this%liq%get_s_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%gl=this%liq%get_g_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      ! Z_l from dp=((Gamma*p-K_s)dV+Gamma*dE)/V, dT=(p/C_v-Gamma*T/V)dV+dE/C_v (tex: Jacobian in the transfer coordinates),
      ! the mass column from homogeneity of degree 0, dg=v*dp-s*dT
      Gam=this%liq%get_gruneisen_from_rho_e(rho=rho,e=el,y=[1.0_WP])
      Ks=rho*st%cl**2
      Cv=st%ml*this%liq%get_cv_from_rho_e(rho=rho,e=el,y=[1.0_WP])
      st%Zl(1,1:2)=[(Gam*st%pl-Ks)/st%Vl,Gam/st%Vl]
      st%Zl(2,1:2)=[st%pl/Cv-Gam*st%Tl/st%Vl,1.0_WP/Cv]
      st%Zl(1:2,3)=-(st%Vl*st%Zl(1:2,1)+st%El*st%Zl(1:2,2))/st%ml
      st%Zl(3,:)=st%Vl/st%ml*st%Zl(1,:)-st%sl*st%Zl(2,:)
      ierr=RELAX_OK
   end subroutine liq_state

   !> Liquid specific volume at (p,T) through the material interface; negative where the liquid has no state
   real(WP) function liq_v_from_p_T(this,p,T) result(v)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      real(WP),                        intent(in) :: p,T
      real(WP) :: rho
      v=-1.0_WP
      if (T.le.0.0_WP) return
      rho=this%liq%get_rho_from_p_T(p=p,T=T,y=[1.0_WP])
      if (rho.gt.0.0_WP) v=1.0_WP/rho
   end function liq_v_from_p_T

   !> Gas p,T,c, vapor x_v,s_v,g_v at its partial pressure and Z_g from st%Vg,Eg,mv,ma
   subroutine gas_state(this,st,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)    :: this
      type(phasic_state),              intent(inout) :: st
      integer,                    intent(out)   :: ierr
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: mg,Rv,Cvg,mR,dTdm
      mg=st%mv+st%ma
      y=0.0_WP; y(this%indV)=st%mv/mg; y(this%indA)=st%ma/mg
      st%pg=this%gas%get_p_from_rho_e(rho=mg/st%Vg,e=st%Eg/mg,y=y)
      st%Tg=this%gas%get_T_from_rho_e(rho=mg/st%Vg,e=st%Eg/mg,y=y)
      if (st%pg.le.0.0_WP.or.st%Tg.le.0.0_WP) then; ierr=RELAX_BAD_GAS; return; end if
      st%cg=this%gas%get_c_from_p_rho(p=st%pg,rho=mg/st%Vg,y=y)
      Rv=this%gas%R(this%indV)
      Cvg=st%mv*this%gas%cv(this%indV)+st%ma*this%gas%cv(this%indA)
      mR=st%mv*Rv+st%ma*this%gas%R(this%indA)
      st%xv=st%mv*Rv/mR
      dTdm=-(this%gas%cv(this%indV)*st%Tg+this%gas%q(this%indV))/Cvg
      st%Zg(1,:)=[-st%pg/st%Vg,mR/(st%Vg*Cvg),(Rv*st%Tg+mR*dTdm)/st%Vg]
      st%Zg(2,:)=[0.0_WP,1.0_WP/Cvg,dTdm]
      ! Vapor chemical potential at its partial pressure, dg_v=Rv*T*(dm_v/m_v-dV_g/V_g)+(Rv-s_v)*dT_g
      st%dry_gas=(st%mv.le.0.0_WP)
      if (st%dry_gas) then
         st%sv=0.0_WP; st%gv=0.0_WP; st%Zg(3,:)=0.0_WP
      else
         st%sv=this%gas%cp(this%indV)*log(st%Tg)-Rv*log(st%xv*st%pg)+this%gas%qp(this%indV)
         st%gv=this%gas%cp(this%indV)*st%Tg+this%gas%q(this%indV)-st%Tg*st%sv
         st%Zg(3,:)=[-Rv*st%Tg/st%Vg,(Rv-st%sv)/Cvg,Rv*st%Tg/st%mv+(Rv-st%sv)*dTdm]
      end if
      ierr=RELAX_OK
   end subroutine gas_state

   !> Transfer directions at the reference state (tex: Matrix decomposition): columns d_M=[1,-p_I,0], d_T=[0,-1,0], d_C=[-v_l,-e_l,-1].
   !> With ns>=2 p_I=p_l (eqs. alg_T, alg_C): only the chemical stage on the pT-relaxed states (relax_C, ns=3) uses it;
   !> the direct iterations (NR_stage_block) hold the p_I of pI_mode at their start state (tex: Iteration of the stages).
   function get_dir(this,st,Pjump,ns) result(D)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      type(phasic_state),              intent(in) :: st
      real(WP),                   intent(in) :: Pjump
      integer, optional,          intent(in) :: ns
      real(WP), dimension(3,3) :: D
      real(WP) :: pI,impl,impg
      select case (this%pI_mode)
      case (pI_liquid); pI=st%pl
      case (pI_gas);    pI=st%pg+Pjump
      case (pI_mean);   pI=0.5_WP*(st%pl+st%pg+Pjump)
      case default
         ! Impedance weighting, also for pI_end: p* is not known at the reference state
         impl=st%ml/st%Vl*st%cl; impg=(st%mv+st%ma)/st%Vg*st%cg
         pI=(impg*st%pl+impl*(st%pg+Pjump))/(impl+impg)
      end select
      if (present(ns)) then
         if (ns.ge.2) pI=st%pl
      end if
      D(:,1)=[1.0_WP,-pI,0.0_WP]
      D(:,2)=[0.0_WP,-1.0_WP,0.0_WP]
      D(:,3)=[-st%Vl/st%ml,-st%El/st%ml,-1.0_WP]
   end function get_dir

   !> Jacobian of F^ in the transfer coordinates, J^=(Z_l+Z_g*A)*D (eq. dFhat), A=[[1,0,0],[p_j,1,0],[0,0,1]]
   function get_jac(this,st,D,Pjump) result(J)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      type(phasic_state),              intent(in) :: st
      real(WP), dimension(3,3),   intent(in) :: D
      real(WP),                   intent(in) :: Pjump
      real(WP), dimension(3,3) :: J,ZA
      ZA=st%Zl+st%Zg
      ZA(:,1)=ZA(:,1)+Pjump*st%Zg(:,2)
      J=matmul(ZA,D)
   end function get_jac

   !> Leading ns x ns block of J^*dx~=-F^ by the elimination M->T->C (eqs. elim_V, schur_1, elim_Q, row_C), pivots -J_MM,
   !> -J'_TT, -J''_CC in K. The last transfer is omega times its Newton value (or dxs), the faster ones hold their forces.
   subroutine eliminate(this,ns,J,F,omega,dx,K,ierr,dxs)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)  :: this
      integer,                    intent(in)  :: ns
      real(WP), dimension(3,3),   intent(in)  :: J
      real(WP), dimension(3),     intent(in)  :: F
      real(WP),                   intent(in)  :: omega
      real(WP), dimension(3),     intent(out) :: dx,K
      integer,                    intent(out) :: ierr
      real(WP), optional,         intent(in)  :: dxs
      real(WP), dimension(2:3,2:3) :: J1
      real(WP), dimension(2:3) :: F1
      real(WP) :: J2,F2
      integer  :: a,b
      ierr=RELAX_SINGULAR
      dx=0.0_WP; K=0.0_WP
      ! First elimination (M)
      K(1)=-J(1,1)
      if (K(1).le.0.0_WP) return
      do a=2,3
         do b=2,3
            J1(a,b)=J(a,b)-J(a,1)*J(1,b)/J(1,1)
         end do
         F1(a)=F(a)-J(a,1)/J(1,1)*F(1)
      end do
      select case (ns)
      case (1)
         dx(1)=-omega*F(1)/J(1,1)
         if (present(dxs)) dx(1)=dxs
      case (2)
         ! Second elimination (T)
         K(2)=-J1(2,2)
         if (K(2).le.0.0_WP) return
         dx(2)=-omega*F1(2)/J1(2,2)
         if (present(dxs)) dx(2)=dxs
         dx(1)=-(F(1)+J(1,2)*dx(2))/J(1,1)
      case (3)
         K(2)=-J1(2,2)
         if (K(2).le.0.0_WP) return
         if (present(dxs)) then
            dx(3)=dxs
         else
            ! Last pivot (C), eq. row_C
            J2=J1(3,3)-J1(3,2)*J1(2,3)/J1(2,2)
            F2=F1(3)-J1(3,2)/J1(2,2)*F1(2)
            K(3)=-J2
            if (K(3).le.0.0_WP) return
            dx(3)=-omega*F2/J2
         end if
         dx(2)=-(F1(2)+J1(2,3)*dx(3))/J1(2,2)
         dx(1)=-(F(1)+J(1,2)*dx(2)+J(1,3)*dx(3))/J(1,1)
      end select
      ierr=RELAX_OK
   end subroutine eliminate

   !> x_l+=D*dx~, x_g-=A*D*dx~ (eq. update_X): exact conservation of V, E_l+E_g+p_j*V_l and mass.
   !> If mv is given, the vapor mass is set to it with the air mass held (tex: Dilute vapor); dx(3) must be mv-m_v.
   subroutine transfer(this,D,dx,Pjump,VF,Q,mv)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)    :: this
      real(WP), dimension(3,3),   intent(in)    :: D
      real(WP), dimension(3),     intent(in)    :: dx
      real(WP),                   intent(in)    :: Pjump
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP), optional,         intent(in)    :: mv
      real(WP), dimension(3) :: dxl
      real(WP) :: ma
      integer :: iV
      iV=this%iVQ()
      ! Air mass as in get_state; max also keeps Q(2)=ma+mv from being reassociated with the old Q(2) (fast-math)
      ma=max(Q(2)-Q(iV),0.0_WP)
      dxl=matmul(D,dx)
      VF=VF+dxl(1)
      Q(3)=Q(3)+dxl(2)
      Q(4)=Q(4)-dxl(2)-Pjump*dxl(1)
      Q(1)=Q(1)+dxl(3)
      Q(2)=Q(2)-dxl(3)
      Q(iV)=Q(iV)-dxl(3)
      if (present(mv)) then
         Q(iV)=mv
         Q(2)=ma+mv
      end if
   end subroutine transfer

   !> Tolerance of force k (eq. conv_tol), floored at the round-off of the difference it measures
   real(WP) function force_tol(this,st,k,Pjump) result(tol)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      type(phasic_state),              intent(in) :: st
      integer,                         intent(in) :: k
      real(WP),                        intent(in) :: Pjump
      real(WP) :: rs
      ! Round-off of each force: V_g=1-VF and E_l (energy conservation) carry the round-off of the cell totals
      rs=10.0_WP*epsilon(1.0_WP)*((st%Vl+st%Vg)*(abs(st%Zl(k,1))+abs(st%Zg(k,1)))+(abs(st%El)+abs(st%Eg))*abs(st%Zl(k,2)) &
      &  +abs(st%Eg)*abs(st%Zg(k,2)))
      select case (k)
      case (1); tol=max(this%FM_tol*max(abs(st%pl),abs(st%pg+Pjump)),10.0_WP*epsilon(1.0_WP)*abs(st%pl)+rs)
      case (2); tol=max(max(this%FT_tol,10.0_WP*epsilon(1.0_WP))*max(st%Tl,st%Tg),10.0_WP*epsilon(1.0_WP)*max(st%Tl,st%Tg)+rs)
      case default; tol=max(this%FC_tol*this%gas%R(this%indV)*st%Tg,10.0_WP*epsilon(1.0_WP)*max(abs(st%gl),abs(st%gv))+rs)
      end select
   end function force_tol

   !> First ns forces below tolerance, each floored at the round-off of the difference it measures
   logical function converged(this,st,ns,Pjump)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      type(phasic_state),              intent(in) :: st
      integer,                    intent(in) :: ns
      real(WP),                   intent(in) :: Pjump
      integer  :: k
      converged=.false.
      do k=1,ns
         if (abs(st%F(k)).gt.this%force_tol(st,k,Pjump)) return
      end do
      converged=.true.
   end function converged

   !> Damped Newton step on the leading ns block in the unknowns V~,E~,u=ln m_v (tex: Iteration of the stages): lambda=1,1/2,...
   !> until both phases have a state. On return st is the state at the new (VF,Q) and K the pivots at the old one.
   subroutine dampedNR(this,ns,st,D,J,Pjump,VF,Q,K,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)    :: this
      integer,                    intent(in)    :: ns
      type(phasic_state),              intent(inout) :: st
      real(WP), dimension(3,3),   intent(in)    :: D,J
      real(WP),                   intent(in)    :: Pjump
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP), dimension(3),     intent(out)   :: K
      integer,                    intent(out)   :: ierr
      real(WP), dimension(size(Q)) :: Qt
      real(WP), dimension(3) :: F,dxN,dx,K2
      real(WP) :: du,dumax,lam,mv0,VFt
      integer :: ih
      F=st%F; mv0=st%mv
      call this%eliminate(ns,J,F,1.0_WP,dxN,K,ierr)
      if (ierr.ne.RELAX_OK) return
      ! Mass step in u=ln m_v; the liquid can give at most m_l
      du=0.0_WP; dumax=0.0_WP
      if (ns.eq.3) then
         du=dxN(3)/mv0
         dumax=log(1.0_WP+st%ml/mv0)
      end if
      lam=1.0_WP
      do ih=0,this%itmax_half
         ierr=RELAX_FAILED
         if (ns.lt.3.or.(lam*du.lt.dumax.and.mv0*exp(lam*du).gt.0.0_WP)) then
            VFt=VF; Qt=Q
            if (ns.eq.3) then
               call this%eliminate(3,J,lam*F,1.0_WP,dx,K2,ierr,dxs=mv0*(exp(lam*du)-1.0_WP))
               call this%transfer(D,dx,Pjump,VFt,Qt,mv=mv0*exp(lam*du))
            else
               call this%transfer(D,lam*dxN,Pjump,VFt,Qt)
            end if
            call this%get_state(VFt,Qt,Pjump,st,ierr)
            if (ierr.eq.RELAX_OK) exit
         end if
         lam=0.5_WP*lam
      end do
      if (ierr.ne.RELAX_OK) return
      VF=VFt; Q=Qt
   end subroutine dampedNR

   !> Direct iteration of the leading ns block (1: M, 2: M,T, 3: M,T,C) from the current state: damped Newton steps with
   !> D held at the start (p_I of the mechanical step), EOS update at every iterate, until the first ns forces vanish
   subroutine NR_stage_block(this,ns,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      integer,                    intent(in)    :: ns
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,                    intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: K
      integer :: it
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      ! A dry gas has g_v=-inf; seed x_v=xv_dry_seed by a transfer along d_C before the chemical stage (tex: Dilute vapor)
      if (ns.eq.3.and.st%dry_gas) then
         D=this%get_dir(st,Pjump)
         call this%transfer(D,[0.0_WP,0.0_WP,this%xv_dry_seed*st%ma*this%gas%R(this%indA)/this%gas%R(this%indV)],Pjump,VF,Q)
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
      end if
      ! Directions held at the start
      D=this%get_dir(st,Pjump)
      do it=1,this%itmax
         if (this%converged(st,ns,Pjump)) then
            select case (ns)
            case (1); this%it_M=it-1
            case (2); this%it_T=it-1
            case (3); this%it_C=it-1
            end select
            ierr=RELAX_OK
            return
         end if
         J=this%get_jac(st,D,Pjump)
         call this%dampedNR(ns,st,D,J,Pjump,VF,Q,K,ierr)
         if (ierr.ne.RELAX_OK) return
         if (it.eq.1) then
            this%KM=K(1)
            if (ns.ge.2) this%KT=K(2)
            if (ns.eq.3) this%KC=K(3)
         end if
      end do
      ierr=RELAX_FAILED
   end subroutine NR_stage_block

   !> Mechanical stage: p_I=p* (pI_end) by Newton on eqs. mech_end, or Newton on the leading 1x1 block with p_I held at the
   !> start state (other pI_mode)
   subroutine relax_M(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                        intent(inout) :: VF
      real(WP), dimension(:),          intent(inout) :: Q
      real(WP),                        intent(in)    :: Pjump
      integer,                         intent(out)   :: ierr
      if (this%pI_mode.eq.pI_end) then
         call this%NR_stage_pI_end(VF,Q,Pjump,ierr)
      else
         call this%NR_stage_block(1,VF,Q,Pjump,ierr)
      end if
   end subroutine relax_M

   !> Mechanical stage with p_I=p*, its own end pressure (tex: Mechanical stage with p_I=p* for a general equation of state):
   !> damped Newton in (delta,p*) on R_l=p_l-p*, R_g=p_g+p_j-p* (eqs. mech_end, mech_newton), every iterate on the path
   !> eq. mech_path from the start state, so that E_l+E_g+p_j*V_l is conserved exactly
   subroutine NR_stage_pI_end(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                        intent(inout) :: VF
      real(WP), dimension(:),          intent(inout) :: Q
      real(WP),                        intent(in)    :: Pjump
      integer,                         intent(out)   :: ierr
      type(phasic_state) :: st0,st
      real(WP), dimension(size(Q)) :: Q0,Qt
      real(WP), dimension(3,3) :: D,J
      real(WP) :: VF0,VFt,dl,ps,Rl,Rg,a11,a12,a21,a22,det,ddl,dps,lam,tol
      integer  :: it,ih
      VF0=VF; Q0=Q
      call this%get_state(VF0,Q0,Pjump,st0,ierr)
      if (ierr.ne.RELAX_OK) return
      ! Start at delta=0 with p* the impedance-weighted p_I of the start state (Algorithm step 2)
      D=this%get_dir(st0,Pjump)
      dl=0.0_WP; ps=-D(2,1)
      st=st0
      do it=1,this%itmax
         Rl=st%pl-ps; Rg=st%pg+Pjump-ps
         tol=this%force_tol(st,1,Pjump)
         if (abs(st%F(1)).le.tol.and.abs(Rl).le.tol) then
            this%it_M=it-1
            ! Pivot K^M at the start with p_I=p* (eq. K_M)
            D(:,1)=[1.0_WP,-ps,0.0_WP]
            J=this%get_jac(st0,D,Pjump)
            this%KM=-J(1,1)
            ierr=RELAX_OK
            return
         end if
         ! Eq. mech_newton: columns d/d(delta) (the terms of K^M at p_I=p*, eq. dM_phase) and d/dp*
         a11=st%Zl(1,1)-ps*st%Zl(1,2);               a12=-1.0_WP-dl*st%Zl(1,2)
         a21=-st%Zg(1,1)+(ps-Pjump)*st%Zg(1,2);      a22=-1.0_WP+dl*st%Zg(1,2)
         det=a11*a22-a12*a21
         ierr=RELAX_SINGULAR
         if (det.eq.0.0_WP) return
         ddl=(a12*Rg-a22*Rl)/det
         dps=(a21*Rl-a11*Rg)/det
         ! Halve the step until both phases have a state
         lam=1.0_WP
         do ih=0,this%itmax_half
            VFt=VF0+dl+lam*ddl; Qt=Q0
            Qt(3)=Q0(3)-(ps+lam*dps)*(dl+lam*ddl)
            Qt(4)=Q0(4)+(ps+lam*dps-Pjump)*(dl+lam*ddl)
            call this%get_state(VFt,Qt,Pjump,st,ierr)
            if (ierr.eq.RELAX_OK) exit
            lam=0.5_WP*lam
         end do
         if (ierr.ne.RELAX_OK) return
         dl=dl+lam*ddl; ps=ps+lam*dps
         VF=VFt; Q=Qt
      end do
      ierr=RELAX_FAILED
   end subroutine NR_stage_pI_end

   !> pT-relaxed state at fixed masses: direct damped Newton on the leading 2x2 block (tex: Iteration of the stages)
   subroutine relax_MT(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                        intent(inout) :: VF
      real(WP), dimension(:),          intent(inout) :: Q
      real(WP),                        intent(in)    :: Pjump
      integer,                         intent(out)   :: ierr
      call this%NR_stage_block(2,VF,Q,Pjump,ierr)
   end subroutine relax_MT

   !> Chemical stage on the pT-relaxed states from a pT-relaxed state (tex: Iteration of the stages): Newton in u=ln m_v on F^_C,
   !> bounded by the sign of F''_C; a trial moves the mass with E~,V~ of eqs. elim_Q, elim_V, then restores F_M=F_T=0 by relax_MT,
   !> and is halved in u if that fails
   subroutine relax_C(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,                    intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP), dimension(size(Q)) :: Qt
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: dx,K
      real(WP) :: u,ut,lo,hi,VFt
      integer :: it,ih
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      ! A dry gas is seeded along d_C (tex: Dilute vapor), then F_M=F_T=0 restored
      if (st%dry_gas) then
         D=this%get_dir(st,Pjump)
         call this%transfer(D,[0.0_WP,0.0_WP,this%xv_dry_seed*st%ma*this%gas%R(this%indA)/this%gas%R(this%indV)],Pjump,VF,Q)
         call this%relax_MT(VF,Q,Pjump,ierr)
         if (ierr.ne.RELAX_OK) return
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
      end if
      lo=-huge(1.0_WP); hi=huge(1.0_WP)
      do it=1,this%itmax
         if (this%converged(st,3,Pjump)) then
            this%it_C=it-1
            ierr=RELAX_OK
            return
         end if
         D=this%get_dir(st,Pjump,3)
         J=this%get_jac(st,D,Pjump)
         call this%eliminate(3,J,st%F,1.0_WP,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
         if (it.eq.1) this%KC=K(3)
         ! Bounds of the root from the sign of F''_C=K^C*dm~ (eq. row_C); a Newton point outside them goes to their midpoint
         u=log(st%mv)
         if (dx(3).gt.0.0_WP) then
            lo=u
         else
            hi=u
         end if
         ut=u+dx(3)/st%mv
         if (ut.le.lo.or.ut.ge.hi) ut=0.5_WP*(lo+hi)
         ! Trial: mass with E~,V~ holding F_M,F_T to first order, then F_M=F_T=0 restored; the liquid gives at most m_l
         do ih=0,this%itmax_half
            ierr=RELAX_FAILED
            if (ut-u.lt.log(1.0_WP+st%ml/st%mv).and.exp(ut).gt.0.0_WP) then
               VFt=VF; Qt=Q
               call this%eliminate(3,J,st%F,1.0_WP,dx,K,ierr,dxs=st%mv*(exp(ut-u)-1.0_WP))
               call this%transfer(D,dx,Pjump,VFt,Qt,mv=exp(ut))
               call this%relax_MT(VFt,Qt,Pjump,ierr)
               if (ierr.eq.RELAX_OK) exit
            end if
            ut=u+0.5_WP*(ut-u)
         end do
         if (ierr.ne.RELAX_OK) return
         VF=VFt; Q=Qt
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
      end do
      ierr=RELAX_FAILED
   end subroutine relax_C

   !> Seed the opposite phase in a metastable near-pure cell (tex Nucleation): superheated liquid (p_vsat(p_l,T_l)>p_l-p_j,
   !> p_l<p_cav) gets vapor at (p_vsat,T_l); supersaturated gas (T_g<T_sat(p_g+p_j,p_v)-Tctol) gets liquid at (p_g+p_j,T_g).
   !> Each seed occupies its own volume and is limited by eq. nuc_cap. Conserves mass and Q(3)+Q(4)+p_j*VF (eq. nuc_cons).
   subroutine nucleate(this,VF,Q,Pjump,seeded)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      logical,                    intent(out)   :: seeded
      real(WP), parameter :: VF_nuc=0.01_WP
      type(phasic_state) :: st
      real(WP), dimension(3) :: sig
      real(WP) :: pvs,rhov,ev,sv,dm,pv,Tsat,rhol,el,sl,A0,dA
      logical  :: conv
      integer  :: iV,it,ier
      seeded=.false.
      iV=this%iVQ()
      if (VF.gt.0.5_WP) then
         ! Cavitation
         if (Q(1).le.0.0_WP) return
         st%Vl=VF; st%El=Q(3); st%ml=Q(1)
         call this%liq_state(st,ier)
         if (ier.ne.RELAX_OK) return
         pvs=this%get_pvsat(st%pl,st%Tl)
         if (pvs.le.st%pl-Pjump.or.st%pl.gt.this%p_cav) return
         rhov=this%gas%get_rho_from_p_T(p=pvs,T=st%Tl,y=vapor())
         ev  =this%gas%get_e_from_p_T  (p=pvs,T=st%Tl,y=vapor())
         sv  =this%gas%get_s_from_p_T  (p=pvs,T=st%Tl,y=vapor())
         ! Seed limit (eqs. nuc_a, nuc_dA, nuc_cap); the liquid gives A^-1*sigma=[v_v,e_v-p_j*v_v,1] per unit mass
         A0=st%gl-(ev+(st%pl-Pjump)/rhov-st%Tl*sv)
         sig=[1.0_WP/rhov,ev-Pjump/rhov,1.0_WP]
         dA=-dot_product(st%Zl(3,:)-st%Zl(1,:)/rhov+sv*st%Zl(2,:),sig)
         dm=min(VF_nuc*rhov,0.5_WP*Q(1))
         if (dA.lt.0.0_WP) dm=min(dm,A0/abs(dA))
         if (A0.le.0.0_WP.or.dm.le.0.0_WP) return
         Q(1)=Q(1)-dm;    Q(2)=Q(2)+dm;    Q(iV)=Q(iV)+dm
         Q(3)=Q(3)-dm*(ev-Pjump/rhov); Q(4)=Q(4)+dm*ev
         VF=VF-dm/rhov
      else
         ! Condensation
         if (Q(2).le.0.0_WP.or.Q(iV).le.0.0_WP.or.VF.ge.1.0_WP) return
         st%Vg=1.0_WP-VF; st%Eg=Q(4); st%mv=min(Q(iV),Q(2)); st%ma=Q(2)-st%mv
         call this%gas_state(st,ier)
         if (ier.ne.RELAX_OK) return
         pv=st%xv*st%pg
         call this%get_Tsat(st%pg+Pjump,pv,st%Tg,Tsat,conv,it)
         if (.not.conv.or.st%Tg.gt.Tsat-this%Tctol) return
         rhol=1.0_WP/this%liq_v_from_p_T(st%pg+Pjump,st%Tg)
         el  =this%liq%get_e_from_p_T  (p=st%pg+Pjump,T=st%Tg,y=[1.0_WP])
         sl  =this%liq%get_s_from_p_T  (p=st%pg+Pjump,T=st%Tg,y=[1.0_WP])
         ! Seed limit (eqs. nuc_a, nuc_dA, nuc_cap); the gas gives A*sigma=[v_l,e_l+p_j*v_l,1] per unit mass
         A0=st%gv-(el+(st%pg+Pjump)/rhol-st%Tg*sl)
         sig=[1.0_WP/rhol,el+Pjump/rhol,1.0_WP]
         dA=-dot_product(st%Zg(3,:)-st%Zg(1,:)/rhol+sl*st%Zg(2,:),sig)
         dm=min(VF_nuc*rhol,0.5_WP*st%mv,0.5_WP*Q(2))
         if (dA.lt.0.0_WP) dm=min(dm,A0/abs(dA))
         if (A0.le.0.0_WP.or.dm.le.0.0_WP) return
         Q(1)=Q(1)+dm;    Q(2)=Q(2)-dm;    Q(iV)=Q(iV)-dm
         Q(3)=Q(3)+dm*el; Q(4)=Q(4)-dm*(el+Pjump/rhol)
         VF=VF+dm/rhol
      end if
      seeded=.true.
   contains
      !> Pure-vapor composition
      function vapor() result(y)
         real(WP), dimension(this%gas%ns) :: y
         y=0.0_WP; y(this%indV)=1.0_WP
      end function vapor
   end subroutine nucleate

   !> Total entropy S=m_l*s_l+m_g*s_g (the igmix entropy includes mixing); 0 if the state is not two-phase
   real(WP) function get_entropy(this,VF,Q) result(S)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      real(WP),                   intent(in) :: VF
      real(WP), dimension(:),     intent(in) :: Q
      type(phasic_state) :: st
      real(WP), dimension(this%gas%ns) :: y
      integer :: ier
      S=0.0_WP
      call this%get_state(VF,Q,0.0_WP,st,ier)
      if (ier.ne.RELAX_OK) return
      y=0.0_WP; y(this%indV)=st%mv/(st%mv+st%ma); y(this%indA)=st%ma/(st%mv+st%ma)
      S=st%ml*st%sl+(st%mv+st%ma)*this%gas%get_s_from_p_T(p=st%pg,T=st%Tg,y=y)
   end function get_entropy

end module relax_igmix_max_ent_inst_class
