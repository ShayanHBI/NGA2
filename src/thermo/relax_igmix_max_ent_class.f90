!> Maximum-entropy (capacity-split) relaxation: SG/NASG liquid + ideal-gas mixture.
!> Implements relaxation as constrained entropy maximization. The extensive state
!> x_k=[V,E,m] is exchanged between the phases until the potentials
!> y_k=[p/T,1/T,-mu_hat/T] equilibrate; linearizing y_k about the current state gives
!> (K_l+K_g)*xi=F, so each phase's potential moves by the OTHER phase's share of the
!> total capacity. Solved sequentially: mechanical (volume), thermal (heat), chemical
!> (mass).
!>
!> Each step is a closed-form solve in response coefficients (compliance W, heat
!> capacity C_p, expansivity alpha_p, chemical stiffness K_CC), and commits
!> (dV_l,dE_l,dm_l) to (VF,Q) so that total mass and total energy (surface energy
!> included) are conserved EXACTLY, while equilibrium is reached to first order only.
!> Iterating a step to convergence is a Newton iteration on the entropy maximum, and
!> lands on the same equilibrium as the exact quadratics of relax_igmix_sg/nasg; a
!> single step is the a-priori predictor the starvation criterion needs.
!>
!> Liquid: class(stiffened_gas) (SG or NASG by inheritance). Gas: class(igmix), vapor
!> at indV, carrier (air) at indA. Vapor mass sits at Q(7+liq%ns+indV-1).
module relax_igmix_max_ent_class
   use precision,             only: WP
   use thermorelax_class,     only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS, &
   &                                RELAX_DEGENERATE,RELAX_SINGULAR,RELAX_VACUUM_VAPOR
   use stiffened_gas_class,   only: stiffened_gas
   use nasg_class,            only: nasg
   use igmix_class,           only: igmix
   use relax_igmix_sg_class,  only: Prelax,PTrelax,PTgrelax
   implicit none
   private

   public :: relax_igmix_max_ent
   public :: cellstate

   !> Cell state and the response coefficients both phases need. Everything here is
   !> extensive per unit cell volume (cell volume 1), matching the Q layout: V_l=VF,
   !> m_l=Q(1), E_l=Q(3), and so on.
   type :: cellstate
      !> Extensive state
      real(WP) :: Vl,Vg,ml,mg,mv,ma,El,Eg
      !> Primitives
      real(WP) :: rhol,rhog,pl,pg,Tl,Tg
      !> Gas composition
      real(WP) :: Yv,xv,pv,Rm,cvm,cpm,qm
      !> Liquid response coefficients
      real(WP) :: Ksl,KTl,apl,Gaml,vspl,hl,sl,gl,Cpl,Cvl,ombm
      !> Gas response coefficients
      real(WP) :: Ksg,KTg,apg,Gamg,vspg,Cpg,Cvg
      !> Vapor partial specific quantities (dX_g/dm_v at fixed p_g,T_g,m_a)
      real(WP) :: vbarv,hbarv,sbarv,mubarv,KCC
      !> Aggregates shared by the thermal and chemical steps
      real(WP) :: Ap,ApT,WT,Cp,D
   end type cellstate

   type, extends(thermorelax) :: relax_igmix_max_ent
      class(stiffened_gas), pointer :: liq=>null()
      class(igmix),         pointer :: gas=>null()
      !> Liquid co-volume, extracted once at init (0 for a plain stiffened gas)
      real(WP) :: bL=0.0_WP
      !> Species indices for vapor and carrier in the gas mixture
      integer  :: indV=0,indA=0
      !> Dispatch
      integer  :: model=Prelax
      !> Outer passes. One pass is M then T then C, which is the block elimination of
      !> the full linearized system; repeating it is Newton on the entropy maximum.
      integer  :: nouter=1
      !> Relative force tolerances: |F_M|<FM_tol*max(|p_l|,|p_g|), |F_T|<FT_tol*T,
      !> |F_C|<FC_tol*Rv*T. Each test is additionally floored at the arithmetic noise of
      !> the difference it measures -- eps*gamma*pinf for p_l-p_g, eps*T for T_l-T_g,
      !> eps*|g| for g_l-mu_v -- because each is formed by cancelling larger numbers and
      !> cannot be resolved below that. Zero therefore means "iterate until each force is
      !> indistinguishable from zero in double precision", which is the default; raise a
      !> tolerance to stop earlier and more cheaply.
      !>
      !> Note what the mechanical scale is NOT: pinf. p_l comes out of the EOS as a
      !> difference against gamma*pinf, so scaling the test on pinf would pin it at
      !> FM_tol*pinf -- 5e-4 Pa here, at every pressure, thousands of times the actual
      !> noise floor in a 1 bar cell. pinf belongs in the floor, not in the scale.
      real(WP) :: FM_tol=0.0_WP,FT_tol=0.0_WP,FC_tol=0.0_WP
      !> Interfacial-pressure closure. The derivation fixes only that p_I lies between
      !> p_g+Pjump and p_l, and that nothing depends on it at first order; it enters
      !> solely through the bracket of W^eff, which differs from 1 by
      !> Gamma_k*(p_k-p_I,k)/K_s,k. 2 (mean) is the default because it is the choice
      !> that makes that bracket smallest for BOTH phases at once -- it is the midpoint
      !> of the admissible interval, so neither |p_l-p_I| nor |p_g+Pjump-p_I| can exceed
      !> half the force. 0/1 reproduce the two one-sided closures of the reference
      !> solver, 3 the post-step liquid pressure the deployed matm-style p_relax uses.
      !>   0 = liquid, p_I=p_l           1 = gas, p_I=p_g+Pjump
      !>   2 = mean (default)            3 = post-step liquid, p_I=p_l-phi_l^M*F_M
      integer  :: pI_mode=2
      !> .true.: conserve E_l+E_g+sigma*A (the derivation). .false.: conserve E_l+E_g
      !> only, i.e. one shared p_I, which is what relax_igmix_sg/nasg do. Identical when
      !> Pjump=0; set .false. to compare against those at Pjump/=0.
      logical  :: conserve_surface=.true.
      !> Volume-fraction bounds the steps are allowed to approach
      real(WP) :: VFmin=1.0e-12_WP,VFmax=1.0_WP-1.0e-12_WP
      !> Liquid packing guard: step rejected if 1-bL*rhol falls below this
      real(WP) :: ombm_min=1.0e-3_WP
      !> Diagnostics from the last apply: iterations taken per channel
      integer  :: nit_M=0,nit_T=0,nit_C=0
      !> Diagnostics: the liquid's share of the capacity in each channel, i.e. the
      !> fraction of the gap the liquid's own potential crossed. These ARE the
      !> starvation shares -- phi->1 means the liquid moved the whole way.
      real(WP) :: phi_M=0.0_WP,phi_T=0.0_WP
      !> Chemical step: reciprocal chemical stiffness (not a share; carries units)
      real(WP) :: phi_C=0.0_WP
      !> Diagnostics: the two halves of the chemical stiffness, K_C=K_pT+K_CC. K_CC is
      !> the composition part (exact in closed form, see the logit root); K_pT is the
      !> thermo-mechanical back-reaction, which no closed form captures. Their ratio
      !> says how much of the chemical step is left approximate once composition is exact.
      real(WP) :: K_pT=0.0_WP,K_CC=0.0_WP
      !> Diagnostics: forces at entry
      real(WP) :: FM0=0.0_WP,FT0=0.0_WP,FC0=0.0_WP
   contains
      procedure :: initialize
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: get_state
      procedure :: get_entropy
      procedure, private :: pass_M
      procedure, private :: pass_T
      procedure, private :: pass_C
      procedure, private :: step_M
      procedure, private :: step_T
      procedure, private :: step_C
      procedure, private :: iVQ
   end type relax_igmix_max_ent

contains

   !> Store EOS pointers, species indices, and extract the liquid co-volume
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_max_ent),     intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      this%liq=>liq
      this%gas=>gas
      this%indV=indV
      this%indA=indA
      ! Co-volume is NASG-only; a plain stiffened gas keeps bL=0
      this%bL=0.0_WP
      select type (liq)
      type is (nasg)
         this%bL=liq%b
      end select
   end subroutine initialize

   !> Index of the transported vapor mass in Q
   integer function iVQ(this)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      iVQ=7+this%liq%ns+this%indV-1
   end function iVQ

   !> Extract the cell state and every response coefficient from (VF,Q).
   !>
   !> Liquid (NASG, v=b+R*T/(p+pinf)):
   !>   K_s=gamma*(p+pinf)/(1-b*rho), K_T=(p+pinf)/(1-b*rho), alpha_p=(1-b*rho)/T,
   !>   Gamma=(gamma-1)/(1-b*rho)
   !> Gas (ideal mixture): K_s=gamma_m*p, K_T=p, alpha_p=1/T, Gamma=R_m/cv_m
   !>
   !> Vapor mole fraction uses the EOS-implied x_v=y_v*R_v/R_m, the same convention
   !> igmix_get_s_from_p_T uses to build partial pressures, so s_bar_v and mu_hat_v stay
   !> consistent with the mixture entropy the EOS reports.
   !>
   !> need_chem=.true. additionally requires a usable vapor partial pressure.
   subroutine get_state(this,VF,Q,st,need_chem,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(in)  :: this
      real(WP),                  intent(in)  :: VF
      real(WP), dimension(:),    intent(in)  :: Q
      type(cellstate),           intent(out) :: st
      logical,                   intent(in)  :: need_chem
      integer,                   intent(out) :: ierr
      real(WP) :: y(this%gas%ns),Rv,cpv,qv,qpv,gammam
      ierr=RELAX_DEGENERATE
      ! Degenerate cells: both phases must be present with positive mass and energy
      if (VF.le.this%VFmin.or.VF.ge.this%VFmax) return
      if (any(Q(1:4).le.0.0_WP)) return
      ! Extensive state (cell volume 1)
      st%Vl=VF; st%Vg=1.0_WP-VF
      st%ml=Q(1); st%mg=Q(2)
      st%El=Q(3); st%Eg=Q(4)
      st%mv=max(0.0_WP,min(Q(this%iVQ()),Q(2)))
      st%ma=st%mg-st%mv
      ! Gas composition and mixture coefficients
      st%Yv=st%mv/st%mg
      y=0.0_WP; y(this%indV)=st%Yv; y(this%indA)=1.0_WP-st%Yv
      st%Rm =sum(y(1:this%gas%ns)*this%gas%R (1:this%gas%ns))
      st%cvm=sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      st%cpm=sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      st%qm =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      ! Primitives, through the EOS accessors so the coefficients below are consistent
      st%rhol=st%ml/st%Vl
      st%rhog=st%mg/st%Vg
      st%pl=this%liq%get_p_from_rho_e(rho=st%rhol,e=st%El/st%ml,y=[1.0_WP])
      st%Tl=this%liq%get_T_from_rho_e(rho=st%rhol,e=st%El/st%ml,y=[1.0_WP])
      st%pg=this%gas%get_p_from_rho_e(rho=st%rhog,e=st%Eg/st%mg,y=y)
      st%Tg=this%gas%get_T_from_rho_e(rho=st%rhog,e=st%Eg/st%mg,y=y)
      ! Liquid soundness: past the co-volume packing limit or below vacuum
      st%ombm=1.0_WP-this%bL*st%rhol
      if (st%ombm.le.this%ombm_min) then; ierr=RELAX_BAD_LIQUID; return; end if
      if (st%pl+this%liq%pinf.le.0.0_WP.or.st%Tl.le.0.0_WP) then; ierr=RELAX_BAD_LIQUID; return; end if
      if (st%pg.le.0.0_WP.or.st%Tg.le.0.0_WP) then; ierr=RELAX_BAD_GAS; return; end if
      ! Liquid response coefficients
      st%Ksl=this%liq%gamma*(st%pl+this%liq%pinf)/st%ombm
      st%KTl=(st%pl+this%liq%pinf)/st%ombm
      st%apl=st%ombm/st%Tl
      st%Gaml=(this%liq%gamma-1.0_WP)/st%ombm
      st%vspl=1.0_WP/st%rhol
      st%Cvl=st%ml*this%liq%cv
      st%Cpl=st%ml*this%liq%cp
      st%hl=this%liq%get_h_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%sl=this%liq%get_s_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%gl=this%liq%get_g_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      ! Gas response coefficients
      gammam=st%cpm/st%cvm
      st%Ksg=gammam*st%pg
      st%KTg=st%pg
      st%apg=1.0_WP/st%Tg
      st%Gamg=st%Rm/st%cvm
      st%vspg=1.0_WP/st%rhog
      st%Cvg=st%mg*st%cvm
      st%Cpg=st%mg*st%cpm
      ! Aggregates. The NASG closed forms collapse these: V_l*T_l*alpha_p,l=V_l-b*m_l
      ! and V_g*T_g*alpha_p,g=V_g, so no small-denominator division appears.
      st%Ap =st%Vl*st%apl+st%Vg*st%apg
      st%ApT=st%Vl*st%Tl*st%apl+st%Vg*st%Tg*st%apg
      st%WT =st%Vl/st%KTl+st%Vg/st%KTg
      st%Cp =st%Cpl+st%Cpg
      st%D  =st%Cp*st%WT-st%Ap*st%ApT
      ! Vapor partial specific quantities. For a single phase D=V*C_v/K_T>0 identically;
      ! the two-phase cross term is what needs checking at run time.
      Rv =this%gas%R (this%indV)
      cpv=this%gas%cp(this%indV)
      qv =this%gas%q (this%indV)
      qpv=this%gas%qp(this%indV)
      st%xv=st%Yv*Rv/st%Rm
      st%pv=st%xv*st%pg
      st%vbarv=Rv*st%Tg/st%pg
      st%hbarv=cpv*st%Tg+qv
      if (need_chem) then
         if (st%mv.le.0.0_WP.or.st%pv.le.0.0_WP) then; ierr=RELAX_VACUUM_VAPOR; return; end if
         st%sbarv=cpv*log(st%Tg)-Rv*log(st%pv)+qpv
         st%mubarv=st%hbarv-st%Tg*st%sbarv
         st%KCC=Rv*st%Tg*(1.0_WP-st%xv)/st%mv
      else
         st%sbarv=0.0_WP; st%mubarv=0.0_WP; st%KCC=0.0_WP
      end if
      ierr=RELAX_OK
   end subroutine get_state

   !> Total entropy of the cell, S=m_l*s_l+m_g*s_g, through the EOS accessors
   !> (igmix_get_s_from_p_T already carries the entropy of mixing). Used to check that
   !> each step produces entropy.
   real(WP) function get_entropy(this,VF,Q) result(S)
      implicit none
      class(relax_igmix_max_ent), intent(in) :: this
      real(WP),                  intent(in) :: VF
      real(WP), dimension(:),    intent(in) :: Q
      type(cellstate) :: st
      real(WP) :: y(this%gas%ns)
      integer  :: ier
      S=0.0_WP
      call this%get_state(VF,Q,st,.false.,ier)
      if (ier.ne.RELAX_OK) return
      y=0.0_WP; y(this%indV)=st%Yv; y(this%indA)=1.0_WP-st%Yv
      S=st%ml*this%liq%get_s_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP]) &
      & +st%mg*this%gas%get_s_from_p_T(p=st%pg,T=st%Tg,y=y)
   end function get_entropy

   !> Dispatch on model. Leaves the cell untouched and reports a RELAX_* code on any
   !> non-OK verdict, the same contract relax_igmix_sg/nasg honour.
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                  intent(in)    :: dt
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      real(WP),                  intent(in)    :: Pjump
      integer,  optional,        intent(out)   :: ierr
      integer :: ier
      ier=RELAX_OK
      select case (this%model)
      case (Prelax);   call this%p_relax  (dt,VF,Q,Pjump,ier)
      case (PTrelax);  call this%pT_relax (dt,VF,Q,Pjump,ier)
      case (PTgrelax); call this%pTg_relax(dt,VF,Q,Pjump,ier)
      case default;    ier=RELAX_FAILED
      end select
      if (present(ierr)) ierr=ier
   end subroutine apply

   !> One mechanical pass: evaluate F_M at the current state and take a single step_M,
   !> unless F_M is already below tolerance, in which case done=.true. and nothing moves.
   subroutine pass_M(this,it,Pjump,VF,Q,done,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      integer,                   intent(in)    :: it
      real(WP),                  intent(in)    :: Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      logical,                   intent(out)   :: done
      integer,                   intent(out)   :: ierr
      type(cellstate) :: st
      real(WP) :: FM,tol
      done=.false.
      call this%get_state(VF,Q,st,.false.,ierr)
      if (ierr.ne.RELAX_OK) return
      FM=st%pl-st%pg-Pjump
      if (it.eq.1) this%FM0=FM
      ! Relative to the cell's own pressure, but never below the arithmetic noise of p_l:
      ! the EOS forms p_l as a difference against gamma*pinf, so eps*gamma*pinf is the
      ! smallest pressure difference double precision can resolve here.
      tol=max(this%FM_tol*max(abs(st%pl),abs(st%pg)),epsilon(1.0_WP)*this%liq%gamma*this%liq%pinf)
      if (abs(FM).le.tol) then; done=.true.; return; end if
      call this%step_M(st,FM,Pjump,VF,Q,ierr)
      if (ierr.ne.RELAX_OK) return
      this%nit_M=this%nit_M+1
   end subroutine pass_M

   !> One thermal pass. Assumes a mechanical pass has just run: step_T is derived at
   !> mechanical equilibrium, which is what lets both phases share one dp.
   subroutine pass_T(this,it,Pjump,VF,Q,done,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      integer,                   intent(in)    :: it
      real(WP),                  intent(in)    :: Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      logical,                   intent(out)   :: done
      integer,                   intent(out)   :: ierr
      type(cellstate) :: st
      real(WP) :: FT,tol
      done=.false.
      call this%get_state(VF,Q,st,.false.,ierr)
      if (ierr.ne.RELAX_OK) return
      FT=st%Tl-st%Tg
      if (it.eq.1) this%FT0=FT
      tol=max(this%FT_tol,epsilon(1.0_WP))*max(st%Tl,st%Tg)
      if (abs(FT).le.tol) then; done=.true.; return; end if
      call this%step_T(st,FT,Pjump,VF,Q,ierr)
      if (ierr.ne.RELAX_OK) return
      this%nit_T=this%nit_T+1
   end subroutine pass_T

   !> One chemical pass. Assumes mechanical and thermal passes have just run: step_C is
   !> derived at thermo-mechanical equilibrium, which is what lets both phases share one
   !> dp and one dT.
   subroutine pass_C(this,it,Pjump,VF,Q,done,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      integer,                   intent(in)    :: it
      real(WP),                  intent(in)    :: Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      logical,                   intent(out)   :: done
      integer,                   intent(out)   :: ierr
      type(cellstate) :: st
      real(WP) :: FC,tol
      done=.false.
      call this%get_state(VF,Q,st,.true.,ierr)
      if (ierr.ne.RELAX_OK) return
      FC=st%gl-st%mubarv
      if (it.eq.1) this%FC0=FC
      ! g_l and mu_v are O(q) and nearly cancel, so eps*|g| is the floor on their difference
      tol=max(this%FC_tol*this%gas%R(this%indV)*st%Tg, &
      &       epsilon(1.0_WP)*max(abs(st%gl),abs(st%mubarv)))
      if (abs(FC).le.tol) then; done=.true.; return; end if
      call this%step_C(st,FC,Pjump,VF,Q,ierr)
      if (ierr.ne.RELAX_OK) return
      this%nit_C=this%nit_C+1
   end subroutine pass_C

   !> Mechanical relaxation: exchange volume until p_l-p_g=Pjump.
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                  intent(in)    :: dt
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      real(WP),                  intent(in)    :: Pjump
      integer,  optional,        intent(out)   :: ierr
      real(WP), dimension(size(Q)) :: Q0
      real(WP) :: VF0
      logical  :: okM
      integer  :: it,ier
      VF0=VF; Q0=Q
      this%nit_M=0
      do it=1,this%nouter
         call this%pass_M(it,Pjump,VF,Q,okM,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         if (okM) exit
      end do
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_relax

   !> Mechanical + thermal. One pass is one step_M followed by one step_T; the pass is
   !> then repeated. step_T is built to leave the Laplace relation intact, so the
   !> mechanical equilibrium step_M just established survives it to second order and
   !> there is nothing to restore in between -- re-running step_M before step_T within a
   !> pass would only redo work the next pass does anyway.
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                  intent(in)    :: dt
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      real(WP),                  intent(in)    :: Pjump
      integer,  optional,        intent(out)   :: ierr
      real(WP), dimension(size(Q)) :: Q0
      real(WP) :: VF0
      logical  :: okM,okT
      integer  :: it,ier
      VF0=VF; Q0=Q
      this%nit_M=0; this%nit_T=0
      do it=1,this%nouter
         call this%pass_M(it,Pjump,VF,Q,okM,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         call this%pass_T(it,Pjump,VF,Q,okT,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         if (okM.and.okT) exit
      end do
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_relax

   !> The full sequential pass: mechanical, then thermal, then chemical.
   !>   step_M -- volume crosses, p_l-p_g -> Pjump
   !>   step_T -- heat crosses, T_l -> T_g, sharing one dp so the Laplace relation holds
   !>   step_C -- mass crosses, g_l -> mu_v, sharing one dp and one dT so both hold
   !> The three directions are conjugate with respect to the stiffness, so one pass IS
   !> the block elimination of the full linearized 3x3 -- nothing is gained by restoring
   !> an earlier channel in the middle of a pass. Repeating the pass is Newton, and the
   !> residual each pass leaves is O(F^2) from the variation of the stiffness alone.
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      real(WP),                  intent(in)    :: dt
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      real(WP),                  intent(in)    :: Pjump
      integer,  optional,        intent(out)   :: ierr
      real(WP), dimension(size(Q)) :: Q0
      real(WP) :: VF0
      logical  :: okM,okT,okC
      integer  :: it,ier
      VF0=VF; Q0=Q
      this%nit_M=0; this%nit_T=0; this%nit_C=0
      do it=1,this%nouter
         call this%pass_M(it,Pjump,VF,Q,okM,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         call this%pass_T(it,Pjump,VF,Q,okT,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         call this%pass_C(it,Pjump,VF,Q,okC,ier)
         if (ier.ne.RELAX_OK) then; VF=VF0; Q=Q0; if (present(ierr)) ierr=ier; return; end if
         if (okM.and.okT.and.okC) exit
      end do
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pTg_relax

   !> One linear mechanical step.
   !>   dp_k=-dV_k/W_k^eff,  W_k^eff=W_k/(1-Gamma_k*(p_k-p_I,k)/K_s,k),  W_k=V_k/K_s,k
   !>   dV_l=W_l^eff*W_g^eff/(W_l^eff+W_g^eff)*F_M,  dp_l=-phi_l*F_M,  dp_g=+phi_g*F_M
   !> Surface energy forces a per-phase working pressure, p_I,l=p_I and p_I,g=p_I-Pjump,
   !> so that dE_l+dE_g=-Pjump*dV_l=-sigma*dA. p_I itself cancels at first order; it
   !> enters only through the small bracket in W^eff, which is why a single
   !> predictor-corrector pass on it is enough.
   subroutine step_M(this,st,FM,Pjump,VF,Q,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      type(cellstate),           intent(in)    :: st
      real(WP),                  intent(in)    :: FM,Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      integer,                   intent(out)   :: ierr
      real(WP) :: Wl,Wg,Wle,Wge,phil,pI,pIl,pIg,bl,bg,dVl,dEl,dEg
      ierr=RELAX_FAILED
      ! Plain compliances and the uncorrected share, used only to place p_I
      Wl=st%Vl/st%Ksl
      Wg=st%Vg/st%Ksg
      if (Wl+Wg.le.0.0_WP) return
      phil=Wg/(Wl+Wg)
      select case (this%pI_mode)
      case (0); pI=st%pl                                  ! liquid side
      case (1); pI=st%pg+Pjump                            ! gas side
      case (3); pI=st%pl-phil*FM                          ! post-step liquid (matm-like)
      case default; pI=0.5_WP*(st%pl+st%pg+Pjump)         ! mean of the admissible interval
      end select
      pIl=pI
      pIg=pI-Pjump
      if (.not.this%conserve_surface) pIg=pI
      ! Effective compliances
      bl=1.0_WP-st%Gaml*(st%pl-pIl)/st%Ksl
      bg=1.0_WP-st%Gamg*(st%pg-pIg)/st%Ksg
      if (bl.le.0.0_WP.or.bg.le.0.0_WP) return
      Wle=Wl/bl
      Wge=Wg/bg
      if (Wle+Wge.le.0.0_WP) return
      this%phi_M=Wge/(Wle+Wge)
      ! Volume the liquid receives, and the p*dV work each phase does on the interface
      dVl=Wle*Wge/(Wle+Wge)*FM
      if (VF+dVl.le.this%VFmin.or.VF+dVl.ge.this%VFmax) return
      dEl=-pIl*dVl
      dEg=+pIg*dVl
      ! Commit; gas energy enforced from conservation, not from its own linearization
      VF=VF+dVl
      Q(3)=Q(3)+dEl
      if (this%conserve_surface) then
         Q(4)=Q(4)-(dEl+Pjump*dVl)
      else
         Q(4)=Q(4)-dEl
      end if
      ierr=RELAX_OK
   end subroutine step_M

   !> One linear thermal step at frozen mass, holding the Laplace relation so both
   !> phases share one pressure change. Solves
   !>   V_l*ap_l*dT_l + V_g*ap_g*dT_g - W_T*dp = 0            (total volume fixed)
   !>   C_p,l*dT_l + C_p,g*dT_g - A_p^T*dp = 0                (total enthalpy - V*dp)
   !>   dT_g - dT_l = F_T                                      (thermal equilibrium)
   !> giving dT_l=-phi_l^T*F_T, dT_g=+phi_g^T*F_T with phi_l^T+phi_g^T=1.
   subroutine step_T(this,st,FT,Pjump,VF,Q,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      type(cellstate),           intent(in)    :: st
      real(WP),                  intent(in)    :: FT,Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      integer,                   intent(out)   :: ierr
      real(WP) :: phil,phig,dTl,dTg,dp,dVl,dHl,dEl
      ierr=RELAX_FAILED
      ! D=C_p*W_T-A_p*A_p^T must be positive for the 3x3 to be solvable
      if (st%D.le.0.0_WP) then; ierr=RELAX_SINGULAR; return; end if
      phil=(st%WT*st%Cpg-st%ApT*st%Vg*st%apg)/st%D
      phig=(st%WT*st%Cpl-st%ApT*st%Vl*st%apl)/st%D
      this%phi_T=phil
      dTl=-phil*FT
      dTg=+phig*FT
      dp=(st%Vl*st%apl*dTl+st%Vg*st%apg*dTg)/st%WT
      ! Liquid volume and energy change. V_l*(1-T_l*ap_l)=b*m_l for NASG, i.e. the dp
      ! coefficient of dH_l is exactly the co-volume term of h=cp*T+b*p+q.
      dVl=st%Vl*st%apl*dTl-st%Vl/st%KTl*dp
      dHl=st%Cpl*dTl+st%Vl*(1.0_WP-st%Tl*st%apl)*dp
      dEl=dHl-st%pl*dVl-st%Vl*dp
      if (VF+dVl.le.this%VFmin.or.VF+dVl.ge.this%VFmax) return
      VF=VF+dVl
      Q(3)=Q(3)+dEl
      if (this%conserve_surface) then
         Q(4)=Q(4)-(dEl+Pjump*dVl)
      else
         Q(4)=Q(4)-dEl
      end if
      ierr=RELAX_OK
   end subroutine step_T

   !> One linear chemical step at thermo-mechanical equilibrium, so both phases share
   !> dp^C and dT^C. Solves
   !>   -W_T*dp + A_p*dT + Dv*dm = 0                           (total volume fixed)
   !>   -A_p^T*dp + C_p*dT + Dh*dm = 0                          (total energy fixed)
   !>   -Dv*dp + Ds*dT - K_CC*dm = F_C                          (chemical equilibrium)
   !> with Dv=v_l-v_bar_v, Dh=h_l-h_bar_v, Ds=s_l-s_bar_v, and the vapor partial
   !> quantities taken at fixed (p_g,T_g,m_a). Eliminating the first two rows,
   !>   dp=Pi*dm, dT=Theta*dm, dm=-F_C/K_C,  K_C=Dv*Pi-Ds*Theta+K_CC.
   subroutine step_C(this,st,FC,Pjump,VF,Q,ierr)
      implicit none
      class(relax_igmix_max_ent), intent(inout) :: this
      type(cellstate),           intent(in)    :: st
      real(WP),                  intent(in)    :: FC,Pjump
      real(WP),                  intent(inout) :: VF
      real(WP), dimension(:),    intent(inout) :: Q
      integer,                   intent(out)   :: ierr
      real(WP) :: Dv,Dh,Ds,Pi,Theta,KC,dml,dp,dT,dVl,dHl,dEl,du
      ierr=RELAX_FAILED
      if (st%D.le.0.0_WP) then; ierr=RELAX_SINGULAR; return; end if
      Dv=st%vspl-st%vbarv
      Dh=st%hl-st%hbarv
      Ds=st%sl-st%sbarv
      Pi   =(st%Cp*Dv-st%Ap*Dh)/st%D
      Theta=(st%ApT*Dv-st%WT*Dh)/st%D
      KC=Dv*Pi-Ds*Theta+st%KCC
      if (KC.le.0.0_WP) then; ierr=RELAX_SINGULAR; return; end if
      this%K_pT=Dv*Pi-Ds*Theta
      this%K_CC=st%KCC
      this%phi_C=1.0_WP/KC
      dml=-FC/KC
      ! mu_v depends on the vapor mass only through ln x_v, so d(mu_v)/d(ln m_v) is the
      ! bounded Rv*T*(1-x_v) while d(mu_v)/dm_v carries a 1/m_v singularity. The solved
      ! increment is therefore the logarithmic one: apply it as such. This keeps both
      ! masses positive identically and reduces to dml when the transfer is small.
      du=-dml/st%mv                                    ! d(ln m_v)
      du=min(du,log(1.0_WP+st%ml/st%mv))               ! cannot evaporate more liquid than there is
      dml=-st%mv*(exp(du)-1.0_WP)
      dp=Pi*dml
      dT=Theta*dml
      dVl=st%Vl*st%apl*dT-st%Vl/st%KTl*dp+st%vspl*dml
      if (VF+dVl.le.this%VFmin.or.VF+dVl.ge.this%VFmax) then; ierr=RELAX_FAILED; return; end if
      dHl=st%Cpl*dT+st%Vl*(1.0_WP-st%Tl*st%apl)*dp+st%hl*dml
      dEl=dHl-st%pl*dVl-st%Vl*dp
      ! Commit: mass moves liquid<->vapor, air is untouched, gas energy enforced
      VF=VF+dVl
      Q(1)=Q(1)+dml
      Q(2)=Q(2)-dml
      Q(this%iVQ())=Q(this%iVQ())-dml
      Q(3)=Q(3)+dEl
      if (this%conserve_surface) then
         Q(4)=Q(4)-(dEl+Pjump*dVl)
      else
         Q(4)=Q(4)-dEl
      end if
      ierr=RELAX_OK
   end subroutine step_C

end module relax_igmix_max_ent_class
