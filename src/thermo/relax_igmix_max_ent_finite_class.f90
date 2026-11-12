!> Finite-rate maximum-entropy relaxation (relaxation_ent_max.tex, Finite-rate relaxation): rate law dx~/dt=R*F^,
!> R=diag(lambda,vartheta,nu). The transfers faster than the finite one are instantaneous, the slower ones inactive;
!> a finite mechanical or thermal stage moves dx~=omega*K^-1*F^0, omega=1-exp(-r*K*dt) (eq. finite_solution), with K its pivot
!> at the step start; the finite chemical stage takes the same exponential where it moves m_v by at most dmv_lin*m_v**, and
!> otherwise integrates the rate law along its frozen path (eq. dm_dt_exact).
module relax_igmix_max_ent_finite_class
   use precision,                 only: WP
   use relax_igmix_max_ent_inst_class, only: relax_igmix_max_ent_inst,phasic_state,max_ent_apply
   use relax_igmix_sg_class,      only: Prelax,PTrelax,PTgrelax
   use thermorelax_class,         only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS
   implicit none
   private

   public :: relax_igmix_max_ent_finite
   public :: restore_faster_eq,energy_split

   !> Closure after the finite step: restore the faster equilibria by their stages (Activation), or Pelanti's
   !> mixture-energy relation at the new V_l and masses (Pelanti 2022 eqs. 59, 70, 82)
   integer, parameter :: restore_faster_eq=1,energy_split=2

   type, extends(relax_igmix_max_ent_inst) :: relax_igmix_max_ent_finite
      real(WP) :: MRR =huge(1.0_WP)   !< Mechanical rate lambda [1/(Pa s)]; huge = instantaneous
      real(WP) :: TRR =huge(1.0_WP)   !< Thermal rate vartheta [W/(m^3 K)]; huge = instantaneous
      real(WP) :: CRR =huge(1.0_WP)   !< Chemical rate nu [kg s/m^5]; huge = instantaneous
      real(WP) :: tauC=-1.0_WP        !< Chemical relaxation time; >0 sets nu=1/(K^C*tauC), K^C at the end state (tex: Chemical step with exact vapor potential)
      !> Default restore_faster_eq: the faster stages end exactly at their equilibria. energy_split is Pelanti's one-pass closure
      !> (F^_M=0 only), whose second-order energy error stays in the phase of smaller heat capacity: T_l/=T_g grows with the transfer
      integer  :: closure=restore_faster_eq
      !> p_I: inherited default pI_end, p* for the instantaneous mechanical stage (NR_stage_pI_end); the finite mechanical step
      !> holds the impedance p_I of its start state, between the pressures (tex: Choice of p_I)
      !> Iteration limits and tolerances of the chemical step (eqs. FC_exact, dm_dt_exact); not converged gives RELAX_FAILED
      integer  :: itmax_bracket=2000      !< Doublings of the evaporation bracket of the root m_v^inf
      integer  :: itmax_mv_eq=200          !< Newton iterations in ln m_v for the root m_v^inf
      real(WP) :: tol_mv_eq=1.0e-14_WP     !< Tolerance on the ln m_v step or bracket width, relative to max(1,|ln m_v|)
      integer  :: itmax_time=60           !< Newton iterations for the upper limit of eq. dm_dt_exact
      real(WP) :: time_tol=1.0e-13_WP     !< Tolerance on the theta step or bracket width of the upper limit
      !> Chemical step: the frozen exponential (eq. finite_solution) where it moves |dm|<=dmv_lin*m_v** (error <~dmv_lin/2, K^C~R_v*T/m_v
      !> nearly constant); the exact integral (eq. dm_dt_exact) is a safeguard for larger steps and a dry gas, not meant to be reached
      real(WP) :: dmv_lin=0.1_WP
      logical  :: chem_exact=.false.      !< Diagnostic: the last chemical step needed the exact integral
   contains
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: finite_step
      procedure :: finite_stage_C
      procedure :: mixture_energy_closure
   end type relax_igmix_max_ent_finite

contains

   !> Rate-separation checks, then the instantaneous class's apply on this, so that it dispatches to the finite workers below
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use messager, only: die
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      logical :: nuc0
      if (this%model.ne.Prelax.and.this%MRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pT and pTg require MRR=huge')
      if (this%model.eq.PTgrelax.and.this%TRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pTg requires TRR=huge')
      ! No nucleation at zero chemical rate: the seed is itself a phase change
      nuc0=this%do_nucleate
      if (this%CRR.le.0.0_WP.and.this%tauC.le.0.0_WP) this%do_nucleate=.false.
      call max_ent_apply(this,dt,VF,Q,Pjump,ierr)
      this%do_nucleate=nuc0
   end subroutine apply

   !> Mechanical relaxation at rate lambda
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      if (this%MRR.ge.huge(1.0_WP)) then
         call this%relax_M(VF,Q,Pjump,ier)
      else
         call this%finite_step(1,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine p_relax

   !> Instantaneous mechanical stage, then thermal relaxation at rate vartheta
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      if (this%TRR.ge.huge(1.0_WP)) then
         call this%relax_MT(VF,Q,Pjump,ier)
      else
         call this%relax_M(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%finite_step(2,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine pT_relax

   !> Instantaneous mechanical and thermal stages, then chemical relaxation at rate nu
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      call this%require_unclamped()
      if (this%CRR.ge.huge(1.0_WP).and.this%tauC.le.0.0_WP) then
         ! Instantaneous chemical stage as in the inst class, through the stages of this type
         call this%relax_MT(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%relax_C(VF,Q,Pjump,ier)
      else
         call this%relax_MT(VF,Q,Pjump,ier)
         if (ier.eq.RELAX_OK) call this%finite_step(3,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine pTg_relax

   !> Finite-rate step of stage ns from a state where the faster forces vanish (eq. finite_solution), then the closure
   subroutine finite_step(this,ns,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      integer,                           intent(in)    :: ns
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: dx,K
      real(WP) :: rate,omega,dm
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      D=this%get_dir(st,Pjump,ns)
      J=this%get_jac(st,D,Pjump)
      if (ns.eq.3) then
         ! Chemical transfer with the composition terms of F^_C exact (eq. dm_dt_exact)
         call this%finite_stage_C(st,D,J,dt,dm,ierr)
         if (ierr.ne.RELAX_OK) return
         call this%eliminate(ns,J,st%F,1.0_WP,dx,K,ierr,dxs=dm)
         if (ierr.ne.RELAX_OK) return
      else
         ! Pivot K of the stage at the step start (eq. stiff_def)
         call this%eliminate(ns,J,st%F,1.0_WP,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
         if (ns.eq.1) then
            rate=this%MRR; this%KM=K(1)
         else
            rate=this%TRR; this%KT=K(2)
         end if
         ! omega=1-exp(-r*K*dt), 1 once r*K*dt>=50
         omega=0.0_WP
         if (rate.gt.0.0_WP) then
            omega=1.0_WP
            if (log(rate)+log(K(ns)*dt).lt.log(50.0_WP)) omega=1.0_WP-exp(-rate*K(ns)*dt)
         end if
         call this%eliminate(ns,J,st%F,omega,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
      end if
      call this%transfer(D,dx,Pjump,VF,Q)
      ! Closure of the faster equilibria
      if (ns.eq.1) return
      select case (this%closure)
      case (energy_split)
         call this%mixture_energy_closure(VF,Q,Pjump,ierr)
      case default
         if (ns.eq.2) then
            call this%relax_M(VF,Q,Pjump,ierr)
         else
            call this%relax_MT(VF,Q,Pjump,ierr)
         end if
      end select
   end subroutine finite_step

   !> Vapor mass dm received by the gas over dt, from a pT-relaxed state (tex: Chemical step with exact vapor potential).
   !> p_l, T, V_g are linear in the transfer with Pi, Theta, chi_V frozen (eq. FC_path); F^_C=g_l(p_l,T)-g_v(p_v,T) exactly
   subroutine finite_stage_C(this,st,D,J,dt,dm,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      type(phasic_state),                     intent(in)    :: st
      real(WP), dimension(3,3),          intent(in)    :: D,J
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(out)   :: dm
      integer,                           intent(out)   :: ierr
      real(WP), dimension(8), parameter :: xg=[-0.9602898564975363_WP,-0.7966664774136267_WP,-0.5255324099163290_WP,-0.1834346424956498_WP, &
      &                                         0.1834346424956498_WP, 0.5255324099163290_WP, 0.7966664774136267_WP, 0.9602898564975363_WP]
      real(WP), dimension(8), parameter :: wg=[ 0.1012285362903763_WP, 0.2223810344533745_WP, 0.3137066458778873_WP, 0.3626837833783620_WP, &
      &                                         0.3626837833783620_WP, 0.3137066458778873_WP, 0.2223810344533745_WP, 0.1012285362903763_WP]
      real(WP), parameter :: thmin=log(1.0e-16_WP),thmax=20.0_WP
      real(WP), dimension(3) :: chi,K
      real(WP) :: Rv,c,m0,Pi,Theta,dVl,minf,rate,tau,tp,tha,thb,thet,lo,hi,u,du,fv,dd,F0,K0
      logical  :: evap,found,conv
      integer  :: it
      dm=0.0_WP; F0=0.0_WP; K0=0.0_WP
      this%chem_exact=.false.
      Rv=this%gas%R(this%indV)
      c=st%ma*this%gas%R(this%indA)/Rv
      m0=st%mv
      ! Pi, Theta, dV_l/dm~: changes of p, T, V_l per unit mass given by the liquid with F_M, F_T held (eqs. VE_C, PiTh)
      call this%eliminate(3,J,[0.0_WP,0.0_WP,0.0_WP],1.0_WP,chi,K,ierr,dxs=1.0_WP)
      if (ierr.ne.RELAX_OK) return
      Pi=-dot_product(st%Zl(1,:),matmul(D,chi))
      Theta=-dot_product(st%Zl(2,:),matmul(D,chi))
      dVl=dot_product(D(1,:),chi)
      ! No transfer at zero rate
      if (this%CRR.le.0.0_WP.and.this%tauC.le.0.0_WP) then
         call chemical_rate(1.0_WP)
         return
      end if
      ! Frozen exponential (eq. finite_solution) at the rate CRR, kept where |dm|<=dmv_lin*m_v** and the liquid keeps mass
      ! (tex: Choice of the solution)
      if (m0.gt.0.0_WP) then
         F0=chemical_force(m0)
         if (F0.eq.0.0_WP) return
         K0=chemical_stiffness(m0)
         if (this%tauC.le.0.0_WP.and.K0.gt.0.0_WP) then
            call chemical_rate(K0)
            this%KC=K0
            dm=frozen_step(K0)
            if (abs(dm).le.this%dmv_lin*m0.and.st%ml-dm.gt.0.0_WP) return
            dm=0.0_WP
         end if
      end if
      ! Root m_v^inf of F^_C: Newton in u=ln m_v inside a bracket where F^_C changes sign
      evap=(m0.le.0.0_WP)
      if (.not.evap) evap=(F0.gt.0.0_WP)
      ierr=RELAX_FAILED
      if (evap) then
         ! Evaporation: double m_v from m** until F^_C<0, on the admissible path, while the liquid has mass
         lo=log(max(m0,1.0e-300_WP))
         hi=log(max(2.0_WP*m0,c*1.0e-12_WP))
         conv=.false.
         do it=1,this%itmax_bracket
            if (.not.path_admissible(exp(hi))) return
            if (chemical_force(exp(hi)).lt.0.0_WP) then; conv=.true.; exit; end if
            if (exp(hi).ge.m0+st%ml) return
            lo=hi; hi=hi+log(2.0_WP)
         end do
         if (.not.conv) return
      else
         ! Condensation: halve m_v from m** until F^_C>0, on the admissible path, while the vapor has mass
         hi=log(m0)
         lo=hi-log(2.0_WP)
         conv=.false.
         do it=1,this%itmax_bracket
            if (.not.path_admissible(exp(lo))) return
            if (chemical_force(exp(lo)).gt.0.0_WP) then; conv=.true.; exit; end if
            if (lo.le.log(1.0e-300_WP)) return
            hi=lo; lo=lo-log(2.0_WP)
         end do
         if (.not.conv) return
      end if
      u=0.5_WP*(lo+hi)
      conv=.false.
      do it=1,this%itmax_mv_eq
         fv=chemical_force(exp(u))
         if (fv.gt.0.0_WP) then; lo=u; else; hi=u; end if
         if (hi-lo.lt.this%tol_mv_eq*max(1.0_WP,abs(u))) then; u=0.5_WP*(lo+hi); conv=.true.; exit; end if
         du=fv/(exp(u)*chemical_stiffness(exp(u)))
         if (u+du.le.lo.or.u+du.ge.hi) du=0.5_WP*(lo+hi)-u
         u=u+du
         if (abs(du).lt.this%tol_mv_eq*max(1.0_WP,abs(u))) then; conv=.true.; exit; end if
      end do
      if (.not.conv) return
      minf=exp(u)
      this%KC=chemical_stiffness(minf)
      call chemical_rate(this%KC)
      ierr=RELAX_OK
      if (rate.le.0.0_WP) return
      ! Relaxation time tau_C: nu needed the root; the frozen exponential at this nu still saves the quadrature
      if (this%tauC.gt.0.0_WP.and.m0.gt.0.0_WP.and.K0.gt.0.0_WP) then
         dm=frozen_step(K0)
         if (abs(dm).le.this%dmv_lin*m0.and.st%ml-dm.gt.0.0_WP) return
         dm=0.0_WP
      end if
      ! Safeguard, not meant to be reached in practice: the frozen step exceeds dmv_lin*m_v** (K^C~R_v*T/m_v not constant over it:
      ! over-condensation, or under-evaporation into a nearly dry gas) or the gas is dry (seeded vapor would set the evaporation, not nu)
      this%chem_exact=.true.
      ! Rate law dt=dm_v/(nu*F^_C) integrated in theta, m_v=minf/(1+exp(-theta)) (evaporation) or
      ! m_v=minf+(m0-minf)*exp(-theta) (condensation), on unit Gauss panels until nu*dt is reached
      dd=m0-minf
      if (evap) then
         thet=thmin
         if (m0.gt.0.0_WP) thet=max(thmin,log(m0/(minf-m0)))
      else
         thet=0.0_WP
      end if
      tau=0.0_WP; found=.false.
      tha=thet
      do while (tha.lt.thmax)
         thb=min(tha+1.0_WP,thmax)
         tp=gauss_panel(tha,thb)
         if (tau+tp.ge.rate*dt) then
            ! Upper limit within the panel: Newton on tau(theta)=nu*dt, bisection-safeguarded
            lo=tha; hi=thb; thet=tha+(rate*dt-tau)/dt_dtheta(tha)
            conv=.false.
            do it=1,this%itmax_time
               if (thet.le.lo.or.thet.ge.hi) thet=0.5_WP*(lo+hi)
               fv=tau+gauss_panel(tha,thet)-rate*dt
               if (fv.gt.0.0_WP) then; hi=thet; else; lo=thet; end if
               if (hi-lo.lt.this%time_tol) then; thet=0.5_WP*(lo+hi); conv=.true.; exit; end if
               du=-fv/dt_dtheta(thet)
               thet=thet+du
               if (abs(du).lt.this%time_tol) then; conv=.true.; exit; end if
            end do
            if (.not.conv) then; ierr=RELAX_FAILED; return; end if
            found=.true.
            exit
         end if
         tau=tau+tp
         tha=thb
      end do
      if (found) then
         dm=mv_of_theta(thet)-m0
      else
         dm=minf-m0
      end if
   contains
      !> Rate nu, or 1/(K*tauC) with K the stiffness at the end state (tex: relaxation time tau_C); diagnostics
      subroutine chemical_rate(Kend)
         real(WP), intent(in) :: Kend
         rate=this%CRR
         if (this%tauC.gt.0.0_WP) rate=1.0_WP/(Kend*this%tauC)
         this%CRR_loc=rate
         if (rate.gt.0.0_WP) this%tauC_loc=1.0_WP/(rate*Kend)
      end subroutine chemical_rate
      !> Frozen exponential dm=omega*F^_C**/K, omega=1-exp(-nu*K*dt) (1 once nu*K*dt>=50)
      real(WP) function frozen_step(Kc)
         real(WP), intent(in) :: Kc
         real(WP) :: om
         om=1.0_WP
         if (log(rate)+log(Kc*dt).lt.log(50.0_WP)) om=1.0_WP-exp(-rate*Kc*dt)
         frozen_step=om*F0/Kc
      end function frozen_step
      !> Liquid pressure, temperature and gas volume on the path (eq. FC_path)
      subroutine fc_path(m,p,T,Vg)
         real(WP), intent(in)  :: m
         real(WP), intent(out) :: p,T,Vg
         p=st%pl-Pi*(m-m0)
         T=st%Tg-Theta*(m-m0)
         Vg=st%Vg-dVl*(m-m0)
      end subroutine fc_path
      !> The path stays where both equations of state hold
      logical function path_admissible(m)
         real(WP), intent(in) :: m
         real(WP) :: p,T,Vg
         call fc_path(m,p,T,Vg)
         path_admissible=(Vg.gt.0.0_WP.and.this%liq_v_from_p_T(p,T).gt.0.0_WP)
      end function path_admissible
      !> Vapor entropy at its partial pressure p_v=m_v*Rv*T/V_g
      real(WP) function vapor_entropy(m,T,Vg)
         real(WP), intent(in) :: m,T,Vg
         vapor_entropy=this%gas%cp(this%indV)*log(T)-Rv*log(m*Rv*T/Vg)+this%gas%qp(this%indV)
      end function vapor_entropy
      !> F^_C(m_v)=g_l(p_l,T)-g_v(p_v,T) on the path (eq. FC_exact)
      real(WP) function chemical_force(m)
         real(WP), intent(in) :: m
         real(WP) :: p,T,Vg
         call fc_path(m,p,T,Vg)
         chemical_force=this%liq%get_g_from_p_T(p=p,T=T,y=[1.0_WP])-(this%gas%cp(this%indV)*T+this%gas%q(this%indV)-T*vapor_entropy(m,T,Vg))
      end function chemical_force
      !> K^C(m_v)=-dF^_C/dm_v (eq. KC_exact)
      real(WP) function chemical_stiffness(m)
         real(WP), intent(in) :: m
         real(WP) :: p,T,Vg
         call fc_path(m,p,T,Vg)
         chemical_stiffness=Pi*this%liq_v_from_p_T(p,T)-this%liq%get_s_from_p_T(p=p,T=T,y=[1.0_WP])*Theta &
         &   +Rv*T*(1.0_WP/m-Theta/T+dVl/Vg)+vapor_entropy(m,T,Vg)*Theta
      end function chemical_stiffness
      !> Vapor mass at theta
      real(WP) function mv_of_theta(t)
         real(WP), intent(in) :: t
         if (evap) then
            mv_of_theta=minf/(1.0_WP+exp(-t))
         else
            mv_of_theta=minf+dd*exp(-t)
         end if
      end function mv_of_theta
      !> d(nu*t)/dtheta=|dm_v/dtheta|/|F^_C|
      real(WP) function dt_dtheta(t)
         real(WP), intent(in) :: t
         real(WP) :: m,s
         if (evap) then
            s=1.0_WP/(1.0_WP+exp(-t))
            m=minf*s
            dt_dtheta=minf*s*(1.0_WP-s)/chemical_force(m)
         else
            m=minf+dd*exp(-t)
            dt_dtheta=dd*exp(-t)/(-chemical_force(m))
         end if
      end function dt_dtheta
      !> 8-point Gauss-Legendre integral of dt_dtheta over [a,b]
      real(WP) function gauss_panel(a,b)
         real(WP), intent(in) :: a,b
         integer :: q
         gauss_panel=0.0_WP
         do q=1,8
            gauss_panel=gauss_panel+wg(q)*dt_dtheta(0.5_WP*(a+b)+0.5_WP*(b-a)*xg(q))
         end do
         gauss_panel=0.5_WP*(b-a)*gauss_panel
      end function gauss_panel
   end subroutine finite_stage_C

   !> Pelanti 2022 eqs. 59, 70, 82: V_l and the masses kept, energies split so that p_l-p_g=p_j with E_l+E_g fixed
   !> (tex: Update of pressure and temperature, option 3: heat dE~ along d_T until F^_M=0)
   subroutine mixture_energy_closure(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      type(phasic_state) :: st
      real(WP) :: dE
      integer  :: it
      ! Newton in dE~ for a general equation of state on the root that defines option 3; one step for a stiffened-gas/NASG liquid and an ideal gas
      do it=1,this%itmax
         call this%get_state(VF,Q,Pjump,st,ierr)
         if (ierr.ne.RELAX_OK) return
         if (this%converged(st,1,Pjump)) return
         dE=st%F(1)/(st%Zl(1,2)+st%Zg(1,2))
         Q(3)=Q(3)-dE
         Q(4)=Q(4)+dE
      end do
      ierr=RELAX_FAILED
   end subroutine mixture_energy_closure

end module relax_igmix_max_ent_finite_class
