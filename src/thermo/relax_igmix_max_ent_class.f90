!> Maximum-entropy sequential relaxation (relaxation_ent_max.tex): SG/NASG liquid + vapor/air ideal-gas mixture.
!> Transfers x~=[V~,E~,m~] (volume received, heat and mass given by the liquid) move x_l=[V_l,E_l,m_l] along D=[d_M,d_T,d_C].
!> Newton on F^=[p_l-p_g-p_j,T_l-T_g,g_l-g_v(p_v,T_g)] with J^=(Z_l+Z_g*A)*D, eliminated M->T->C (eqs. newton_rows-row_C);
!> each stage is iterated with the EOS update (option 2) and D held at the stage start.
module relax_igmix_max_ent_class
   use precision,            only: WP
   use relax_igmix_sg_class, only: relax_igmix_sg,Prelax,PTrelax,PTgrelax
   use stiffened_gas_class,  only: stiffened_gas
   use igmix_class,          only: igmix
   use thermorelax_class,    only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_DEGENERATE,RELAX_SINGULAR,RELAX_NUC_FAILED
   implicit none
   private

   public :: relax_igmix_max_ent,mestate,max_ent_apply
   public :: pI_liquid,pI_gas,pI_mean,pI_impedance

   !> Interface pressure p_I of d_M (eq. pI_w): w_M=0 liquid, 1 gas, 1/2 mean, or Pelanti's impedance weighting
   integer, parameter :: pI_liquid=0,pI_gas=1,pI_mean=2,pI_impedance=3

   !> Cell state (cell volume 1): extensive x_k, intensive z_k=[p,T,g], forces F^ and Z_k=dz_k/dx_k
   type :: mestate
      real(WP) :: Vl,El,ml,Vg,Eg,mv,ma
      real(WP) :: pl,Tl,gl,sl,cl
      real(WP) :: pg,Tg,gv,sv,xv,cg
      logical  :: dry
      real(WP), dimension(3)   :: F
      real(WP), dimension(3,3) :: Zl,Zg
   end type mestate

   type, extends(relax_igmix_sg) :: relax_igmix_max_ent
      !> Liquid co-volume and packing clamp (0 and 1 for a plain stiffened gas)
      real(WP) :: bL=0.0_WP,brhomax=1.0_WP
      !> Interface-pressure closure of d_M
      integer  :: pI_mode=pI_impedance
      !> Stage tolerances: |F_M|<FM_tol*p, |F_T|<FT_tol*T, |F_C|<FC_tol*Rv*T (floored at round-off)
      real(WP) :: FM_tol=1.0e-10_WP,FT_tol=1.0e-10_WP,FC_tol=1.0e-10_WP
      integer  :: itmax=50
      !> Vapor mole fraction seeded into a dry gas before the chemical stage
      real(WP) :: xv_seed=1.0e-10_WP
      !> Diagnostics: iterations of the last stages, and pivots K^M,K^T,K^C at their start (eq. stiff_def)
      integer  :: nit_M=0,nit_T=0,nit_C=0
      real(WP) :: KM=0.0_WP,KT=0.0_WP,KC=0.0_WP
   contains
      procedure :: initialize
      procedure :: apply=>max_ent_apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: iVQ
      procedure :: get_state
      procedure :: get_D
      procedure :: get_jac
      procedure :: eliminate
      procedure :: transfer
      procedure :: stage
      procedure :: nucleate
      procedure :: get_entropy
   end type relax_igmix_max_ent

contains

   !> Parent initialization (EOS pointers, saturation coefficients), then the NASG co-volume
   subroutine initialize(this,liq,gas,indV,indA)
      use nasg_class, only: nasg
      implicit none
      class(relax_igmix_max_ent),    intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      call this%relax_igmix_sg%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      this%bL=0.0_WP; this%brhomax=1.0_WP
      select type (liq)
      type is (nasg)
         this%bL=liq%b
         this%brhomax=liq%brhomax
         this%ES=liq%b/gas%R(indV)
      end select
   end subroutine initialize

   !> Index of the vapor partial density in Q
   integer function iVQ(this)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      iVQ=7+this%liq%ns+this%indV-1
   end function iVQ

   !> Model dispatch on mixture cells; pTg first seeds the opposite phase in metastable near-pure cells.
   !> A cell that is not relaxed is returned untouched with a non-OK code.
   subroutine max_ent_apply(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      use messager,     only: die
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
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
      this%nit_M=0; this%nit_T=0; this%nit_C=0
      this%KM=0.0_WP; this%KT=0.0_WP; this%KC=0.0_WP
      VF0=VF; Q0=Q
      mixture=(VF.ge.VFlo.and.VF.le.VFhi.and.Q(1).gt.0.0_WP.and.Q(2).gt.0.0_WP)
      ! Nucleation in near-pure cells
      if (this%model.eq.PTgrelax.and.this%do_nucleate.and.(VF.lt.VF_nuc.or.VF.gt.1.0_WP-VF_nuc)) call this%nucleate(VF,Q,this%nucleated)
      if (.not.(mixture.or.this%nucleated)) then
         if (present(ierr)) ierr=merge(RELAX_NUC_FAILED,RELAX_DEGENERATE,this%model.eq.PTgrelax)
         return
      end if
      select case (this%model)
      case (Prelax);   call this%p_relax  (dt,VF,Q,Pjump,ier)
      case (PTrelax);  call this%pT_relax (dt,VF,Q,Pjump,ier)
      case (PTgrelax); call this%pTg_relax(dt,VF,Q,Pjump,ier)
      case default;    call die('[relax_igmix_max_ent apply] model must be Prelax, PTrelax or PTgrelax')
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
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer :: ier
      call this%stage(1,VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine p_relax

   !> Mechanical then thermal stage: maximum of S at fixed masses (pT-relaxed state)
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer :: ier
      call this%stage(1,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) call this%stage(2,VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine pT_relax

   !> Mechanical, thermal then chemical stage: maximum of S over all exchanges (pTg-relaxed state, F=0)
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                   intent(in)    :: dt
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,  optional,         intent(out)   :: ierr
      integer :: ier
      call this%stage(1,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) call this%stage(2,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) call this%stage(3,VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine pTg_relax

   !> Phasic states, forces F^ (eq. Fhat_def) and Z_k=dz_k/dx_k from (VF,Q); x_l=[V_l,E_l,m_l], x_g=[V_g,E_g,m_v]
   subroutine get_state(this,VF,Q,Pjump,st,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(in)  :: this
      real(WP),                   intent(in)  :: VF
      real(WP), dimension(:),     intent(in)  :: Q
      real(WP),                   intent(in)  :: Pjump
      type(mestate),              intent(out) :: st
      integer,                    intent(out) :: ierr
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: vl,el,vb,dpdv,dpde,dTdv,dTde,Rv,Cvg,mR,dTdm
      ierr=RELAX_DEGENERATE
      st%Vl=VF; st%Vg=1.0_WP-VF
      st%ml=Q(1); st%El=Q(3); st%Eg=Q(4)
      st%mv=min(Q(this%iVQ()),Q(2)); st%ma=Q(2)-st%mv
      if (st%Vl.le.0.0_WP.or.st%Vg.le.0.0_WP.or.st%ml.le.0.0_WP.or.Q(2).le.0.0_WP.or.st%mv.lt.0.0_WP) return
      ! Liquid
      vl=st%Vl/st%ml; el=st%El/st%ml
      st%pl=this%liq%get_p_from_rho_e(rho=1.0_WP/vl,e=el,y=[1.0_WP])
      st%Tl=this%liq%get_T_from_rho_e(rho=1.0_WP/vl,e=el,y=[1.0_WP])
      if (st%pl+this%liq%pinf.le.0.0_WP.or.st%Tl.le.0.0_WP) then; ierr=RELAX_BAD_LIQUID; return; end if
      st%sl=this%liq%get_s_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%gl=this%liq%get_g_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%cl=this%liq%get_c_from_p_rho(p=st%pl,rho=1.0_WP/vl,y=[1.0_WP])
      ! Liquid Z_l from p=(gamma-1)(e-q)/(v-b)-gamma*pinf, T=(e-q-pinf*(v-b))/cv, dg=v*dp-s*dT
      vb=max(vl-this%bL,(1.0_WP-this%brhomax)*vl)
      dpdv=-(st%pl+this%liq%gamma*this%liq%pinf)/vb
      dpde=(this%liq%gamma-1.0_WP)/vb
      dTdv=-this%liq%pinf/this%liq%cv
      dTde=1.0_WP/this%liq%cv
      st%Zl(1,:)=[dpdv,dpde,-(vl*dpdv+el*dpde)]/st%ml
      st%Zl(2,:)=[dTdv,dTde,-(vl*dTdv+el*dTde)]/st%ml
      st%Zl(3,:)=vl*st%Zl(1,:)-st%sl*st%Zl(2,:)
      ! Gas mixture
      y=0.0_WP; y(this%indV)=st%mv/Q(2); y(this%indA)=st%ma/Q(2)
      st%pg=this%gas%get_p_from_rho_e(rho=Q(2)/st%Vg,e=st%Eg/Q(2),y=y)
      st%Tg=this%gas%get_T_from_rho_e(rho=Q(2)/st%Vg,e=st%Eg/Q(2),y=y)
      if (st%pg.le.0.0_WP.or.st%Tg.le.0.0_WP) then; ierr=RELAX_BAD_GAS; return; end if
      st%cg=this%gas%get_c_from_p_rho(p=st%pg,rho=Q(2)/st%Vg,y=y)
      Rv=this%gas%R(this%indV)
      Cvg=st%mv*this%gas%cv(this%indV)+st%ma*this%gas%cv(this%indA)
      mR=st%mv*Rv+st%ma*this%gas%R(this%indA)
      st%xv=st%mv*Rv/mR
      dTdm=-(this%gas%cv(this%indV)*st%Tg+this%gas%q(this%indV))/Cvg
      st%Zg(1,:)=[-st%pg/st%Vg,mR/(st%Vg*Cvg),(Rv*st%Tg+mR*dTdm)/st%Vg]
      st%Zg(2,:)=[0.0_WP,1.0_WP/Cvg,dTdm]
      ! Vapor chemical potential at its partial pressure, dg_v=Rv*T*(dm_v/m_v-dV_g/V_g)+(Rv-s_v)*dT_g
      st%dry=(st%mv.le.0.0_WP)
      if (st%dry) then
         st%sv=0.0_WP; st%gv=0.0_WP; st%Zg(3,:)=0.0_WP
      else
         st%sv=this%gas%cp(this%indV)*log(st%Tg)-Rv*log(st%xv*st%pg)+this%gas%qp(this%indV)
         st%gv=this%gas%cp(this%indV)*st%Tg+this%gas%q(this%indV)-st%Tg*st%sv
         st%Zg(3,:)=[-Rv*st%Tg/st%Vg,(Rv-st%sv)/Cvg,Rv*st%Tg/st%mv+(Rv-st%sv)*dTdm]
      end if
      ! Forces (eq. Fhat_def)
      st%F=[st%pl-st%pg-Pjump,st%Tl-st%Tg,st%gl-st%gv]
      if (st%dry) st%F(3)=0.0_WP
      ierr=RELAX_OK
   end subroutine get_state

   !> Transfer directions at the reference state (eqs. d_M, d_T, d_C): columns [1,-p_I,0], [0,-1,0], [-v_l,-e_l,-1]
   function get_D(this,st,Pjump) result(D)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      type(mestate),              intent(in) :: st
      real(WP),                   intent(in) :: Pjump
      real(WP), dimension(3,3) :: D
      real(WP) :: pI,impl,impg
      select case (this%pI_mode)
      case (pI_liquid); pI=st%pl
      case (pI_gas);    pI=st%pg+Pjump
      case (pI_mean);   pI=0.5_WP*(st%pl+st%pg+Pjump)
      case default
         impl=st%ml/st%Vl*st%cl; impg=(st%mv+st%ma)/st%Vg*st%cg
         pI=(impg*st%pl+impl*(st%pg+Pjump))/(impl+impg)
      end select
      D(:,1)=[1.0_WP,-pI,0.0_WP]
      D(:,2)=[0.0_WP,-1.0_WP,0.0_WP]
      D(:,3)=[-st%Vl/st%ml,-st%El/st%ml,-1.0_WP]
   end function get_D

   !> Jacobian of F^ in the transfer coordinates, J^=(Z_l+Z_g*A)*D (eq. dFhat), A=[[1,0,0],[p_j,1,0],[0,0,1]]
   function get_jac(this,st,D,Pjump) result(J)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      type(mestate),              intent(in) :: st
      real(WP), dimension(3,3),   intent(in) :: D
      real(WP),                   intent(in) :: Pjump
      real(WP), dimension(3,3) :: J,ZA
      ZA=st%Zl+st%Zg
      ZA(:,1)=ZA(:,1)+Pjump*st%Zg(:,2)
      J=matmul(ZA,D)
   end function get_jac

   !> Leading ns x ns block of J^*dx~=-F^ by the elimination M->T->C (eqs. elim_V, schur_1, elim_Q, row_C), pivots -J_MM,
   !> -J'_TT, -J''_CC in K. The last transfer is omega times its Newton value (or dxs), the faster ones hold their forces.
   subroutine eliminate(this,ns,J,F,mv,logm,omega,dx,K,ierr,dxs)
      implicit none
      class(relax_igmix_max_ent), intent(in)  :: this
      integer,                    intent(in)  :: ns
      real(WP), dimension(3,3),   intent(in)  :: J
      real(WP), dimension(3),     intent(in)  :: F
      real(WP),                   intent(in)  :: mv
      logical,                    intent(in)  :: logm
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
            ! Departure: Newton unknown ln m_v for the mass (K_CC~1/m_v)
            if (logm) dx(3)=mv*(exp(dx(3)/mv)-1.0_WP)
         end if
         dx(2)=-(F1(2)+J1(2,3)*dx(3))/J1(2,2)
         dx(1)=-(F(1)+J(1,2)*dx(2)+J(1,3)*dx(3))/J(1,1)
      end select
      ierr=RELAX_OK
   end subroutine eliminate

   !> x_l+=D*dx~, x_g-=A*D*dx~ (eq. update_X): exact conservation of V, E_l+E_g+p_j*V_l and mass
   subroutine transfer(this,D,dx,Pjump,VF,Q)
      implicit none
      class(relax_igmix_max_ent), intent(in)    :: this
      real(WP), dimension(3,3),   intent(in)    :: D
      real(WP), dimension(3),     intent(in)    :: dx
      real(WP),                   intent(in)    :: Pjump
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP), dimension(3) :: dxl
      integer :: iV
      iV=this%iVQ()
      dxl=matmul(D,dx)
      VF=VF+dxl(1)
      Q(3)=Q(3)+dxl(2)
      Q(4)=Q(4)-dxl(2)-Pjump*dxl(1)
      Q(1)=Q(1)+dxl(3)
      Q(2)=Q(2)-dxl(3)
      Q(iV)=Q(iV)-dxl(3)
   end subroutine transfer

   !> Stage ns=1 (M), 2 (M,T) or 3 (M,T,C): Newton on the leading block of J^ with D held at the stage start,
   !> EOS update at every iterate, until the first ns forces vanish
   subroutine stage(this,ns,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      integer,                    intent(in)    :: ns
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      real(WP),                   intent(in)    :: Pjump
      integer,                    intent(out)   :: ierr
      type(mestate) :: st
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: dx,K
      integer :: it
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      ! Departure: a dry gas has g_v=-inf; seed x_v=xv_seed by a transfer along d_C before the chemical stage
      if (ns.eq.3.and.st%dry) then
         D=this%get_D(st,Pjump)
         call this%transfer(D,[0.0_WP,0.0_WP,this%xv_seed*st%ma*this%gas%R(this%indA)/this%gas%R(this%indV)],Pjump,VF,Q)
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
      end if
      ! Directions held at the stage start
      D=this%get_D(st,Pjump)
      do it=1,this%itmax
         if (converged()) then
            select case (ns)
            case (1); this%nit_M=it-1
            case (2); this%nit_T=it-1
            case (3); this%nit_C=it-1
            end select
            ierr=RELAX_OK
            return
         end if
         J=this%get_jac(st,D,Pjump)
         call this%eliminate(ns,J,st%F,st%mv,.true.,1.0_WP,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
         if (it.eq.1) then
            select case (ns)
            case (1); this%KM=K(1)
            case (2); this%KT=K(2)
            case (3); this%KC=K(3)
            end select
         end if
         call this%transfer(D,dx,Pjump,VF,Q)
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
      end do
      ierr=RELAX_FAILED
   contains
      !> First ns forces below tolerance, each floored at the round-off of the difference it measures
      logical function converged()
         real(WP) :: tol
         converged=.false.
         tol=max(this%FM_tol*max(abs(st%pl),abs(st%pg+Pjump)),10.0_WP*epsilon(1.0_WP)*this%liq%gamma*this%liq%pinf)
         if (abs(st%F(1)).gt.tol) return
         if (ns.ge.2) then
            tol=max(this%FT_tol,10.0_WP*epsilon(1.0_WP))*max(st%Tl,st%Tg)
            if (abs(st%F(2)).gt.tol) return
         end if
         if (ns.ge.3) then
            tol=max(this%FC_tol*this%gas%R(this%indV)*st%Tg,10.0_WP*epsilon(1.0_WP)*max(abs(st%gl),abs(st%gv)))
            if (abs(st%F(3)).gt.tol) return
         end if
         converged=.true.
      end function converged
   end subroutine stage

   !> Seed the opposite phase in a metastable near-pure cell: superheated liquid (p_vsat(p_l,T_l)>p_l, p_l<p_cav) gets vapor
   !> at (p_vsat,T_l); supersaturated gas (T_g<T_sat(p_g,p_v)-Tctol) gets liquid at (p_g,T_g). Each seed occupies its own volume.
   !> Conserves mass and Q(3)+Q(4); the pTg end state does not depend on the seed, which only sets the Newton start.
   subroutine nucleate(this,VF,Q,seeded)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                   intent(inout) :: VF
      real(WP), dimension(:),     intent(inout) :: Q
      logical,                    intent(out)   :: seeded
      real(WP), parameter :: VF_nuc=0.01_WP,dT_nuc_max=0.05_WP
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: rhol,pl,Tl,pvs,rhov,ev,dm,rhog,pg,Tg,pv,Tsat,el,eg,cvg
      logical  :: conv
      integer  :: iV,it
      seeded=.false.
      iV=this%iVQ()
      if (VF.gt.0.5_WP) then
         ! Cavitation
         if (Q(1).le.0.0_WP) return
         rhol=Q(1)/VF
         pl=this%liq%get_p_from_rho_e(rho=rhol,e=Q(3)/Q(1),y=[1.0_WP])
         Tl=this%liq%get_T_from_p_rho(p=pl,rho=rhol,y=[1.0_WP])
         if (pl.le.-this%liq%pinf.or.Tl.le.0.0_WP) return
         pvs=this%get_pvsat(pl,Tl)
         if (pvs.le.pl.or.pl.gt.this%p_cav) return
         y=0.0_WP; y(this%indV)=1.0_WP
         rhov=this%gas%get_rho_from_p_T(p=pvs,T=Tl,y=y)
         ev  =this%gas%get_e_from_p_T  (p=pvs,T=Tl,y=y)
         ! Seed volume dm/rhov taken from the liquid, capped so the liquid is compressed at most to p_vsat
         dm=min(VF_nuc*rhov,0.5_WP*Q(1),rhov*VF*(pvs-pl)/(rhol*this%liq%get_c_from_p_rho(p=pl,rho=rhol,y=[1.0_WP])**2))
         ! Seed energy limited to a dT_nuc_max relative heating of the liquid
         if (Q(3)/Q(1).gt.ev) dm=min(dm,dT_nuc_max*Tl*this%liq%cv*Q(1)/(Q(3)/Q(1)-ev))
         if (dm.le.0.0_WP) return
         Q(1)=Q(1)-dm;    Q(2)=Q(2)+dm;    Q(iV)=Q(iV)+dm
         Q(3)=Q(3)-dm*ev; Q(4)=Q(4)+dm*ev
         VF=VF-dm/rhov
      else
         ! Condensation
         if (Q(2).le.0.0_WP.or.Q(iV).le.0.0_WP.or.VF.ge.1.0_WP) return
         y=0.0_WP; y(this%indV)=min(Q(iV)/Q(2),1.0_WP); y(this%indA)=1.0_WP-y(this%indV)
         rhog=Q(2)/(1.0_WP-VF)
         pg=this%gas%get_p_from_rho_e(rho=rhog,e=Q(4)/Q(2),y=y)
         Tg=this%gas%get_T_from_p_rho(p=pg,rho=rhog,y=y)
         if (pg.le.0.0_WP.or.Tg.le.0.0_WP) return
         pv=y(this%indV)*this%gas%R(this%indV)/sum(y*this%gas%R(1:this%gas%ns))*pg
         call this%get_Tsat(pg,pv,Tg,Tsat,conv,it)
         if (.not.conv.or.Tg.gt.Tsat-this%Tctol) return
         rhol=this%liq%get_rho_from_p_T(p=pg,T=Tg,y=[1.0_WP])
         el  =this%liq%get_e_from_p_T  (p=pg,T=Tg,y=[1.0_WP])
         eg  =this%gas%get_e_from_p_T  (p=pg,T=Tg,y=y)
         cvg =sum(y*this%gas%cv(1:this%gas%ns))
         dm=min(VF_nuc*rhol,0.5_WP*Q(iV),0.5_WP*Q(2))
         ! Seed energy limited to a dT_nuc_max relative heating of the gas
         if (eg.gt.el) dm=min(dm,dT_nuc_max*Tg*cvg*Q(2)/(eg-el))
         if (dm.le.0.0_WP) return
         Q(1)=Q(1)+dm;    Q(2)=Q(2)-dm;    Q(iV)=Q(iV)-dm
         Q(3)=Q(3)+dm*el; Q(4)=Q(4)-dm*el
         VF=VF+dm/rhol
      end if
      seeded=.true.
   end subroutine nucleate

   !> Total entropy S=m_l*s_l+m_g*s_g (the igmix entropy includes mixing); 0 if the state is not two-phase
   real(WP) function get_entropy(this,VF,Q) result(S)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      real(WP),                   intent(in) :: VF
      real(WP), dimension(:),     intent(in) :: Q
      type(mestate) :: st
      real(WP), dimension(this%gas%ns) :: y
      integer :: ier
      S=0.0_WP
      call this%get_state(VF,Q,0.0_WP,st,ier)
      if (ier.ne.RELAX_OK) return
      y=0.0_WP; y(this%indV)=st%mv/(st%mv+st%ma); y(this%indA)=st%ma/(st%mv+st%ma)
      S=st%ml*st%sl+(st%mv+st%ma)*this%gas%get_s_from_p_T(p=st%pg,T=st%Tg,y=y)
   end function get_entropy

end module relax_igmix_max_ent_class
