!> SG liquid and ideal gas relaxation model
module relax_sg_ig_class
   use precision,   only: WP
   use relax_class, only: relax
   use sg_class,    only: sg
   use igmix_class, only: igmix
   implicit none
   private

   public :: relax_sg_ig,Mv,Ma

   !> Molar mass of vapor and air [kg/mol]
   real(WP), parameter :: Mv=0.0180153_WP
   real(WP), parameter :: Ma=0.02897_WP

   type, extends(relax) :: relax_sg_ig
      !> Typed EOS pointers (liq is class(sg) so nasg IS-A sg works for child class)
      class(sg),    pointer :: liq => null()
      class(igmix), pointer :: gas => null()
      !> Species indices for vapor and air in the gas mixture
      integer :: indV=0,indA=0
      !> Saturation curve coefficients
      real(WP) :: AS=0.0_WP,BS=0.0_WP,CS=0.0_WP,DS=0.0_WP,ES=0.0_WP
      !> Convergence tolerances
      real(WP) :: p_tol   =1.0e-4_WP
      real(WP) :: Yv_tol  =1.0e-4_WP
      real(WP) :: Tsat_tol=1.0e-5_WP
      real(WP) :: rho_tol =1.0e-4_WP
      real(WP) :: rhoe_tol=1.0e-4_WP
      real(WP) :: F1_tol  =1.0e-4_WP
      real(WP) :: F2_tol  =1.0e-4_WP
      !> Iteration limits
      integer :: Tsat_itmax=40
      integer :: NR_itmax  =40
   contains
      procedure :: initialize   =>relax_sg_ig_initialize
      procedure :: relax_p      =>relax_sg_ig_relax_p
      procedure :: relax_pT     =>relax_sg_ig_relax_pT
      procedure :: relax_pTg    =>relax_sg_ig_relax_pTg
      procedure :: get_T_lvg    =>relax_sg_ig_get_T_lvg
      procedure :: get_coeffs_lv=>relax_sg_ig_get_coeffs_lv
      procedure :: pTsat        =>relax_sg_ig_pTsat
      procedure :: dpTsatdT     =>relax_sg_ig_dpTsatdT
      procedure :: dpTsatdp_lv  =>relax_sg_ig_dpTsatdp_lv
      procedure :: get_Tsat     =>relax_sg_ig_get_Tsat
      procedure :: get_xv       =>relax_sg_ig_get_xv
   end type relax_sg_ig

contains

   !> Initialize: store EOS pointers and compute saturation curve coefficients
   subroutine relax_sg_ig_initialize(this,liq,gas,indV,indA)
      class(relax_sg_ig), intent(inout)    :: this
      class(sg),    target, intent(in)     :: liq
      class(igmix), target, intent(in)     :: gas
      integer, intent(in) :: indV,indA
      real(WP) :: cpV,cvV,RV
      this%liq=>liq
      this%gas=>gas
      this%indV=indV
      this%indA=indA
      cpV=gas%get_species_cp(indV)
      cvV=gas%get_species_cv(indV)
      RV=cpV-cvV
      this%AS=(liq%cp-cpV+gas%get_species_qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%get_species_q(indV))            /RV
      this%CS=(cpV-liq%cp)                                /RV
      this%DS=(liq%cp-liq%cv)                             /RV
      ! ES stays 0.0_WP (SG has no co-volume b); NASG overrides this after calling parent init
   end subroutine relax_sg_ig_initialize

   !> Mechanical relaxation
   subroutine relax_sg_ig_relax_p(this,VF,Q,Pjump)
      class(relax_sg_ig),      intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: y
      real(WP) :: PL,PG,ZL,ZG,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Get phasic pressures
      PL=this%liq%get_p_from_rho_e(Q(1)/(       VF),Q(3)/Q(1))
      PG=this%gas%get_p_from_rho_e(Q(2)/(1.0_WP-VF),Q(4)/Q(2),y)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(8)=Yv*Q(2)
         deallocate(y)
         return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP
         deallocate(y)
         return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(PL,Q(1)/(       VF))**2
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(PG,Q(2)/(1.0_WP-VF),y)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(this%liq%gamma-1.0_WP)*Pint+2.0_WP*this%liq%gamma*this%liq%pinf
      coeffG=(gammaG        -1.0_WP)*Pint
      a=1.0_WP+gammaG*VF+this%liq%gamma*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+gammaG)*VF*PL-(1.0_WP+this%liq%gamma)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=VF*((this%liq%gamma-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+this%liq%gamma)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      deallocate(y)
   end subroutine relax_sg_ig_relax_p

   !> Mechanical and thermal relaxation
   subroutine relax_sg_ig_relax_pT(this,VF,Q,Pjump)
      class(relax_sg_ig), intent(inout)      :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: y
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! ================ First step for mechanical relaxation ================
      call this%relax_p(VF,Q,Pjump)
      ! ================= Second step for thermal relaxation =================
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Setup quadratic problem
      a=Q(1)*this%liq%cv+Q(2)*cvG
      b=this%liq%q*this%liq%cv*(this%liq%gamma-1.0_WP)*Q(1)**2+qG*cvG*(gammaG-1.0_WP)*Q(2)**2+&
      &  Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf+Q(2)*cvG*this%liq%pinf                +&
      &  Q(1)*Q(2)*(this%liq%q*cvG*(gammaG-1.0_WP)+qG*this%liq%cv*(this%liq%gamma-1.0_WP))   -&
      &  sum(Q(3:4))*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*cvG*(gammaG-1.0_WP))
      d=cvG*(gammaG-1.0_WP)*this%liq%pinf*(qG*Q(2)**2+this%liq%q*Q(1)*Q(2)-sum(Q(3:4))*Q(2))
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y)
         return
      end if
      ! Get equilibrium volume fraction
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq/(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*cvG*(gammaG-1.0_WP)*(Peq+this%liq%pinf))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then
         VFeq=0.0_WP
         Peq=max(Peq,-this%liq%pinf)
      end if
      if (VFeq.gt.1.0_WP) then
         VFeq=1.0_WP
         Peq=max(Peq,0.0_WP)
      end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(Peq,Q(1)/max(VFeq,tiny(1.0_WP)))
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(Peq,Q(2)/max(1.0_WP-VFeq,tiny(1.0_WP)),y)
      VF=VFeq
      deallocate(y)
   end subroutine relax_sg_ig_relax_pT

   !> Mechanical, thermal, and chemical relaxation
   subroutine relax_sg_ig_relax_pTg(this,VF,Q,Pjump)
      class(relax_sg_ig),      intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,y
      real(WP) :: VF0,p,T,Yv
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: RHOL,RHOG
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP), parameter :: p_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,pv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
      logical :: chem_relax
      ! Cavitation nucleation: Conservatively move a little mass and energy from liquid to vapor so the chemical relaxation starts
      ! from a non-stiff initial condition.
      nucleation: block
         real(WP), parameter :: VF_cav_seed=1.0e-7_WP
         real(WP) :: rhoL_cav,pL_cav,TL_cav,pv_sat,rhoV_cav,eV_cav,drho,de
         real(WP) :: y_nuc(this%gas%ns)
         if (VF.ge.1.0_WP-VF_cav_seed) then
            rhoL_cav=Q(1)/VF
            pL_cav=this%liq%get_p_from_rho_e(rhoL_cav,Q(3)/Q(1))
            TL_cav=this%liq%get_T_from_p_rho(pL_cav,rhoL_cav)
            if (pL_cav.le.-this%liq%pinf.or.TL_cav.le.0.0_WP) exit nucleation
            ! Saturation vapor pressure at current liquid state
            pv_sat=exp(this%AS+(this%BS+this%ES*pL_cav)/TL_cav)*TL_cav**this%CS*(pL_cav+this%liq%pinf)**this%DS
            if (pv_sat.le.pL_cav) exit nucleation
            ! Superheated liquid: seed a tiny vapor phase.
            ! Transfer volumetric mass drho and its energy de from the liquid to the vapor; total rho and rhoe conserved.
            y_nuc           =0.0_WP
            y_nuc(this%indV)=1.0_WP
            rhoV_cav=this%gas%get_rho_from_p_T(pL_cav,TL_cav,y_nuc)
            eV_cav=this%gas%get_e_from_p_T(pL_cav,TL_cav,y_nuc)
            drho=VF_cav_seed*rhoV_cav
            de=drho*eV_cav
            Q(1)=Q(1)-drho; Q(2)=drho
            Q(3)=Q(3)-de;   Q(4)=de
            Q(8)=drho
            VF=1.0_WP-VF_cav_seed
         end if
      end block nucleation
      ! ================ First and second steps for mechanical and thermal relaxation ================
      call this%relax_pT(VF,Q,Pjump)
      ! ================= Third step for chemical relaxation =================
      ! Set gas EOS coefficients from vapor mass fraction (Q(8), Q(2) unchanged by relax_pT)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Recover p and T from the dominant phase for better numerical stability
      if (VF.gt.0.5_WP) then
         p=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1))
         T=this%liq%get_T_from_p_rho(p=p,rho=Q(1)/VF)
      else
         p=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
         T=this%gas%get_T_from_p_rho(p=p,rho=Q(2)/(1.0_WP-VF),y=y)
      end if
      ! Store input state to the chemical relaxation algorithm
      VF0=VF
      allocate(Q0(size(Q)))
      Q0=Q
      rho0 =sum(Q0(1:2))
      rhoe0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! print '(A)',       '============ PT_relax ============='
      ! print '(A,ES15.7)','p   =',p
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A,ES15.7)','T   =',T
      ! print '(A)',       '==================================='
      ! Check if chemical relaxation should be activated
      chem_relax=activate_chem(p,T,Yv)
      if (.not.chem_relax) then
         call restore()
         call dealloc()
         return
      end if
      ! Solve chemical equilibrium without touching the conserved variables.
      ! If the inert gas content is only a numerical trace, the LVG equations
      ! become ill-conditioned; use the pure liquid-vapor branch instead.
      if (Yv.gt.Yv_pure) then
         ! print*,'****************** Using pure LV chemical relaxation!'
         Yv=Yvmax
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
         call solve_lv(p,T,chem_relax)
      else
         call solve_lvg(p,T,Yv,chem_relax)
      end if
      ! Skip if not converged
      if (.not.chem_relax) then
         call restore()
         call dealloc()
         return
      end if
      ! print '(A)',       '=========== PTg_relax ============'
      ! print '(A,ES15.7)','p   =',p
      ! print '(A,ES15.7)','T   =',T
      ! print '(A,ES15.7)','Yv  =',Yv
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A)',       '=================================='
      ! Adjust vapor mass fraction and update gas EOS parameters
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Adjust densities
      RHOL=this%liq%get_rho_from_p_T(p=p,T=T)
      RHOG=this%gas%get_rho_from_p_T(p=p,T=T,y=y)
      ! Adjust VF
      VF=(rho0-RHOG)/(RHOL-RHOG)
      ! Clean up solution
      if (VF.lt.0.0_WP) then
         VF=0.0_WP
         p=max(p,-this%liq%pinf)
      end if
      if (VF.gt.1.0_WP) then
         VF=1.0_WP
         p=max(p,0.0_WP)
      end if
      ! Adjust conserved quantities (velocities remain unchanged)
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*this%liq%get_e_from_p_T(p=p,T=T)
      Q(4)=Q(2)*this%gas%get_e_from_p_T(p=p,T=T,y=y)
      Q(8)=Q(2)*Yv
      ! Evaluate conservation
      if (.not.check_cons()) then
         call restore()
         call dealloc()
         ! print*,"****************** Conservation is violated. Skipping the cell!!"
         return
      end if
      ! Check VOF
      if (VF.lt.VFmin) then
         call restore()
         call dealloc()
         ! print*,"****************** Not enough liquid to vaporize. Skipping the cell!!"
         return
      end if
      ! Release memory
      call dealloc()
   contains
      !> Reset the output to the initial values fed into chemical relaxation algorithm
      subroutine restore()
         VF=VF0
         Q=Q0
      end subroutine restore
      !> Release memory allocated for Q0 and y
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
      end subroutine dealloc
      !> Check sanity of a vapor partial pressure
      logical function check_pv(pv_)
         real(WP), intent(in) :: pv_
         check_pv=(pv_.gt.p_eps)
      end function check_pv
      !> Check mass and energy conservation
      logical function check_cons()
         real(WP) :: re,ee
         re=(sum(Q(1:2))-rho0)/rho0
         ee=(sum(Q(3:4))-rhoe0)/rhoe0
         check_cons=(abs(re).le.this%rho_tol).and.(abs(ee).le.this%rhoe_tol)
      end function check_cons
      !> Decide if chemical relaxation needs to be activated and adjust Yv if needed
      logical function activate_chem(p_,T_,Yv_)
         real(WP), intent(in)    :: p_,T_
         real(WP), intent(inout) :: Yv_
         real(WP) :: xv,pv_,Tsat
         integer  :: Tsat_it
         logical  :: conv
         activate_chem=.false.
         ! Get vapor mole fraction and partial pressure
         xv=this%get_xv(Yv_)
         pv_=xv*p_
         if ((Yv_.le.Yv_dry).or.(pv_.le.pv_dry).or.(.not.check_pv(pv_))) then
            ! Dry/nearly-dry air edge case: pv_ is zero or so tiny that solving Tsat(pv_) is log-singular/ill-conditioned.
            ! Seed Yv from saturation at the current thermally-relaxed state and then continue with the ordinary LVG Newton solve.
            pv_=exp(this%AS+(this%BS+this%ES*p_)/T_)*T_**this%CS*(p_+this%liq%pinf)**this%DS
            if (.not.check_pv(pv_)) then
               ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
               return
            end if
            if (pv_.ge.p_) then
               ! pv_sat >= p: flash regime. Seed with a small Yv so that get_T stays well-conditioned.
               Yv_=sqrt(0.0001_WP)
            else
               xv=pv_/p_
               Yv_=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
            end if
            Yv_=max(Yvmin,min(Yvmax,Yv_))
            ! print '(A,ES15.7)','Seeded Yv from saturation at Teq=',Yv_
         else
            ! Get saturation temperature at current pressure
            call this%get_Tsat(p_,pv_,T_,Tsat,conv,Tsat_it)
            if (.not.conv) then
               ! print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
               return
            end if
            ! print '(A)',       '========== Finding Tsat ==========='
            ! print '(A,I2)','Tsat it= ',Tsat_it
            ! print '(A,ES15.7)','Tsat   =',Tsat
            ! print '(A)',       '==================================='
            ! Activate chemical relaxation only for vaporizing metastable states
            if (T_.le.Tsat) return
         end if
         activate_chem=.true.
      end function activate_chem
      !> Equilibrium T for LV
      real(WP) function get_T_lv(ap,bp,dp)
         real(WP), intent(in) :: ap,bp,dp
         get_T_lv=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
      end function get_T_lv
      !> dT/dp from the quadratic coefficients for LV
      real(WP) function get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         get_dTdp_lv=(ap*(-dbpdp+(bp*dbpdp-2.0_WP*(dapdp*dp+ap*ddpdp))/sqrt(bp**2-4.0_WP*ap*dp))-dapdp*(-bp+sqrt(bp**2-4.0_WP*ap*dp)))/(2.0_WP*ap**2)
      end function get_dTdp_lv
      !> Energy residual for LVG
      real(WP) function rhoe_res_lvg(p_,T_,Yv_)
         real(WP), intent(in) :: p_,T_,Yv_
         rhoe_res_lvg=(rho0*(1.0_WP-Yv_)-rhoA0)*this%liq%get_e_from_p_T(p_,T_)+&
         &             rhoA0*this%gas%get_e_from_p_T(p_,T_,[Yv_,1.0_WP-Yv_])-  &
         &             rhoe0*(1.0_WP-Yv_)
      end function rhoe_res_lvg
      !> Energy residual for LV
      real(WP) function rhoe_res_lv(p_,T_)
         real(WP), intent(in) :: p_,T_
         real(WP) :: rho_l,rho_g
         rho_l=this%liq%get_rho_from_p_T(p=p_,T=T_)
         rho_g=this%gas%get_rho_from_p_T(p=p_,T=T_,y=[1.0_WP,0.0_WP])
         rhoe_res_lv=(rho0-rho_g)/(rho_l-rho_g)*this%liq%get_e_from_p_T(p=p_,T=T_)+                    &
         &           (rho_l-rho0)/(rho_l-rho_g)*this%gas%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP,0.0_WP])-&
         &            rhoe0
      end function rhoe_res_lv
      !> Pure liquid and vapor chemical relaxation
      subroutine solve_lv(p_eq,T_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq
         logical,  intent(out)   :: conv
         real(WP) :: pOld,ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: dTdp
         integer  :: it
         ! Iteratively solve for the equilibrium pressure in pure vapor case
         conv=.false.
         do it=1,this%NR_itmax
            ! Get the coefficients
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Get temperature
            T_eq=get_T_lv(ap,bp,dp)
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Newton-Raphson iteration: Pure-vapor branch: Y_v=1 so the vapor mole fraction x_v=1.
            pOld=p_eq
            p_eq=pOld-this%pTsat(pOld,pOld,T_eq)/this%dpTsatdp_lv(pOld,T_eq,dTdp)
            ! Evaluate the error
            if (abs((p_eq-pOld)/pOld).lt.this%p_tol) then
               conv=.true.
               exit
            end if
         end do
         ! Check convergence
         if (.not.conv) then
            ! print*,"****************** p iterations blew up. Skipping the cell!!"
            return
         end if
         ! Update equilibrium temperature
         call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
         T_eq=get_T_lv(ap,bp,dp)
      end subroutine solve_lv
      !> Pure liquid and vapor-gas mixture chemical relaxation
      subroutine solve_lvg(p_eq,T_eq,Yv_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: conv
         real(WP) :: xv,pv
         real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
         real(WP) :: p_pert,Yv_pert,T_pert,xv_pert,ppv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
         real(WP) :: pOld,YvOld,p_err,Yv_err
         real(WP) :: alpha,res0,res_try
         real(WP) :: p_try,Yv_try,T_try,xv_try,ppv_try,F1_try,F2_try
         real(WP) :: Yv_max_phys
         integer  :: it,lsit
         logical  :: accepted
         ! Calculate the absolute physical ceiling for Yv based on available liquid
         Yv_max_phys=1.0_WP-(rhoA0/rho0)
         ! Iteratively solve for the equilibrium pressure and vapor mass fraction
         conv=.false.
         p_err=10.0_WP*this%p_tol
         Yv_err=10.0_WP*this%Yv_tol
         do it=1,this%NR_itmax
            ! Get temperature
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0)
            ! Get vapor partial pressure and mole fraction
            xv=this%get_xv(Yv_eq)
            pv=xv*p_eq
            if (.not.check_pv(pv)) return
            ! Get residuals at current state
            F1=this%pTsat(p_eq,pv,T_eq)/p_eq
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            res0=sqrt(F1**2+F2**2)
            ! Perturbation in p
            p_pert=p_eq*(1.0_WP+fd_eps)
            ppv_pert=xv*p_pert
            if (.not.check_pv(ppv_pert)) return
            ! Get the corresponding T
            T_pert=this%get_T_lvg(p_pert,Yv_eq,rho0,rhoA0)
            ! Get residuals
            F1p=this%pTsat(p_pert,ppv_pert,T_pert)/p_eq
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)/rhoe0
            dF1dp=(F1p-F1)/(p_pert-p_eq)
            dF2dp=(F2p-F2)/(p_pert-p_eq)
            ! Perturbation in Yv
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yvmax-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            xv_pert=this%get_xv(Yv_pert)
            ppv_pert=xv_pert*p_eq
            if (.not.check_pv(ppv_pert)) return
            ! Get the corresponding T
            T_pert=this%get_T_lvg(p_eq,Yv_pert,rho0,rhoA0)
            ! Get residuals
            F1Y=this%pTsat(p_eq,ppv_pert,T_pert)/p_eq
            F2Y=rhoe_res_lvg(p_eq,T_pert,Yv_pert)/rhoe0
            dF1dYv=(F1Y-F1)/(Yv_pert-Yv_eq)
            dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
            ! Solve 2x2 system: J*[dp; dYv]=-[F1; F2]
            detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
            if (abs(detJ).lt.1.0e-30_WP) exit
            dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
            dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
            ! Direction preserving step limiter
            step_limit: block
               real(WP) :: ms
               ms=1.0_WP
               ! 1. Prevent pressure from changing by more than 50% in a single step
               if (abs(dp_nr).gt.0.5_WP*p_eq) ms=min(ms,0.5_WP*p_eq/abs(dp_nr))
               ! 2. Prevent Yv from crossing physical boundaries
               if (dYv_nr.gt.0.0_WP) then
                  if (Yv_eq+dYv_nr.ge.Yv_max_phys)  ms=min(ms,0.9_WP*(Yv_max_phys-Yv_eq)/dYv_nr)
                  if (Yv_eq+dYv_nr.ge.Yvmax-fd_eps) ms=min(ms,0.9_WP*(Yvmax-fd_eps-Yv_eq)/dYv_nr)
               else if (dYv_nr.lt.0.0_WP) then
                  if (Yv_eq+dYv_nr.le.Yvmin+fd_eps) ms=min(ms,0.9_WP*(Yv_eq-Yvmin-fd_eps)/abs(dYv_nr))
               end if
               ! Scale the Newton-Raphson step uniformly to keep pointing directly at the root
               dp_nr=dp_nr*ms
               dYv_nr=dYv_nr*ms
            end block step_limit
            ! Damped Newton-Raphson update
            pOld=p_eq
            YvOld=Yv_eq
            alpha=1.0_WP
            lsit=0
            if ((abs(F1).lt.F_line_search_tol).and.(abs(F2).lt.F_line_search_tol)) then
               p_eq=pOld+dp_nr
               Yv_eq=YvOld+dYv_nr
            else
               accepted=.false.
               do while (alpha.gt.1.0e-8_WP)
                  lsit=lsit+1
                  p_try=p_eq+alpha*dp_nr
                  Yv_try=Yv_eq+alpha*dYv_nr
                  ! Keep the trial state inside the physical/log-safe domain
                  if (p_try.le.p_eps) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  if ((Yv_try.le.Yvmin+fd_eps).or.(Yv_try.ge.Yvmax-fd_eps)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  if ((rho0*(1.0_WP-Yv_try)-rhoA0).le.0.0_WP) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  T_try=this%get_T_lvg(p_try,Yv_try,rho0,rhoA0)
                  if (T_try.le.0.0_WP) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  xv_try=this%get_xv(Yv_try)
                  ppv_try=xv_try*p_try
                  if (.not.check_pv(ppv_try)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  F1_try=this%pTsat(p_try,ppv_try,T_try)/p_try
                  F2_try=rhoe_res_lvg(p_try,T_try,Yv_try)/rhoe0
                  res_try=sqrt(F1_try**2+F2_try**2)
                  if (res_try.lt.res0) then
                     p_eq=p_try
                     Yv_eq=Yv_try
                     T_eq=T_try
                     accepted=.true.
                     exit
                  end if
                  alpha=0.5_WP*alpha
               end do
               if (.not.accepted) then
                  ! print*,"****************** Line search failed at it=",it," with alpha=",alpha
                  exit
               end if
            end if
            ! Evaluate errors
            p_err=abs((p_eq-pOld)/pOld)
            Yv_err=abs((Yv_eq-YvOld)/(YvOld+1.0e-30_WP))
            ! Refresh xv and pv from the accepted solution
            xv=this%get_xv(Yv_eq)
            pv=xv*p_eq
            ! Get temperature
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0)
            ! Evaluate residuals
            F1=this%pTsat(p_eq,pv,T_eq)/p_eq
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            if ((p_err.lt.this%p_tol).and.(Yv_err.lt.this%Yv_tol).and.(abs(F1).lt.this%F1_tol).and.(abs(F2).lt.this%F2_tol)) then
               conv=.true.
               exit
            end if
         end do
         ! Check convergence
         if (.not.conv) then
            ! print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
            return
         end if
      end subroutine solve_lvg
   end subroutine relax_sg_ig_relax_pTg

   !> p-T saturation curve (general form: ES=0 for SG reduces exactly to the SG formula)
   real(WP) function relax_sg_ig_pTsat(this,pl_,pv_,T_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,pv_,T_
      relax_sg_ig_pTsat=this%AS+(this%BS+this%ES*pl_)/T_+this%CS*log(T_)+this%DS*log(pl_+this%liq%pinf)-log(pv_)
   end function relax_sg_ig_pTsat

   !> d(pTsat)/dT
   real(WP) function relax_sg_ig_dpTsatdT(this,pl_,T_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,T_
      relax_sg_ig_dpTsatdT=-(this%BS+this%ES*pl_)/T_**2+this%CS/T_
   end function relax_sg_ig_dpTsatdT

   !> d(pTsat)/dp for the LV solve (pv_ = pl_)
   real(WP) function relax_sg_ig_dpTsatdp_lv(this,p_,T_,dTdp)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdp
      relax_sg_ig_dpTsatdp_lv=this%dpTsatdT(p_,T_)*dTdp+this%ES/T_+this%DS/(p_+this%liq%pinf)-1.0_WP/p_
   end function relax_sg_ig_dpTsatdp_lv

   !> Safeguarded Newton-Raphson for saturation temperature at fixed (pl_, pv_)
   subroutine relax_sg_ig_get_Tsat(this,pl_,pv_,Tguess,Tsat,conv,Tsat_it)
      class(relax_sg_ig), intent(inout) :: this
      real(WP), intent(in)  :: pl_,pv_,Tguess
      real(WP), intent(out) :: Tsat
      logical,  intent(out) :: conv
      integer,  intent(out) :: Tsat_it
      real(WP) :: Tlo,Thi,Told,Tnew,Flo,Fhi,Fold,Fnew,dFold
      integer  :: it,expand_it
      conv=.false.
      Tsat_it=0
      ! Broad physical bracket for the current water EOS.  Newton is accepted only when it remains inside this bracket;
      ! otherwise the update falls back to bisection.
      Tlo=250.0_WP
      Thi=900.0_WP
      Flo=this%pTsat(pl_,pv_,Tlo)
      Fhi=this%pTsat(pl_,pv_,Thi)
      ! Expand the bracket if needed.  In the normal liquid-vapor range pTsat is monotone in T, so these one-sided expansions are enough.
      expand_it=0
      do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
         if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
            Tlo=max(1.0_WP,0.8_WP*Tlo)
            Flo=this%pTsat(pl_,pv_,Tlo)
         else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
            Thi=1.2_WP*Thi
            Fhi=this%pTsat(pl_,pv_,Thi)
         else
            exit
         end if
         expand_it=expand_it+1
      end do
      if (Flo*Fhi.gt.0.0_WP) then
         ! print*,'****************** Could not bracket Tsat!',pl_,pv_,Flo,Fhi,Tlo,Thi
         return
      end if
      ! Use the caller's guess only after clamping it to the safe bracket.
      Tsat=max(Tlo,min(Thi,Tguess))
      do it=1,this%Tsat_itmax
         Told=Tsat
         Fold=this%pTsat(pl_,pv_,Told)
         dFold=this%dpTsatdT(pl_,Told)
         if (abs(dFold).gt.tiny(1.0_WP)) then
            Tnew=Told-Fold/dFold
         else
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         ! Safeguard: reject Newton steps that leave the bracket or are NaN.
         if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) then
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         Fnew=this%pTsat(pl_,pv_,Tnew)
         if (Fnew.ne.Fnew) then
            Tnew=0.5_WP*(Tlo+Thi)
            Fnew=this%pTsat(pl_,pv_,Tnew)
         end if
         ! Update the bracket around the root.
         if (Flo*Fnew.le.0.0_WP) then
            Thi=Tnew
            Fhi=Fnew
         else
            Tlo=Tnew
            Flo=Fnew
         end if
         Tsat_it=it
         Tsat=Tnew
         if ((abs((Tnew-Told)/max(abs(Told),tiny(1.0_WP))).lt.this%Tsat_tol).or.(abs(Fnew).lt.this%F1_tol)) then
            conv=.true.
            return
         end if
      end do
   end subroutine relax_sg_ig_get_Tsat

   !> Vapor mole fraction from mass fraction
   real(WP) function relax_sg_ig_get_xv(this,Yv_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: Yv_
      relax_sg_ig_get_xv=Yv_*Ma/(Yv_*Ma+(1.0_WP-Yv_)*Mv)
   end function relax_sg_ig_get_xv

   !> Equilibrium T from energy conservation (SG form: no co-volume b correction)
   real(WP) function relax_sg_ig_get_T_lvg(this,p_,Yv_,rho0,rhoA0)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,Yv_,rho0,rhoA0
      relax_sg_ig_get_T_lvg=(1.0_WP-Yv_)/                                                                                                      &
      &         ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)+                                             &
      &           rhoA0*((this%gas%get_species_gamma(this%indV)-1.0_WP)*this%gas%get_species_cv(this%indV)*Yv_+                                  &
      &                  (this%gas%get_species_gamma(this%indA)-1.0_WP)*this%gas%get_species_cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function relax_sg_ig_get_T_lvg

   !> Update the quadratic coefficients of equilibrium temperature equation (SG form: PinfG=0)
   subroutine relax_sg_ig_get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      class(relax_sg_ig), intent(in)  :: this
      real(WP), intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      ! Coefficients (rho0 = sum of phase masses; rhoe0 = sum of phase energies — both invariant during steps 1-2)
      ap=rho0*this%liq%cv*cvG*((GammaG-1.0_WP)*(p_eq+this%liq%gamma*this%liq%pinf)-(this%liq%gamma-1.0_WP)*p_eq)
      bp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv*p_eq-(GammaG-1.0_WP)*cvG*(p_eq+this%liq%pinf))              +&
      &  rho0*((GammaG-1.0_WP)*cvG*this%liq%q*(p_eq+this%liq%pinf)-(this%liq%gamma-1.0_WP)*this%liq%cv*qG*p_eq) +&
      &  cvG*p_eq*(p_eq+this%liq%pinf)-this%liq%cv*p_eq*(p_eq+this%liq%gamma*this%liq%pinf)
      dp=(qG-this%liq%q)*(p_eq+this%liq%pinf)*p_eq
      ! Pressure derivative of the coefficients
      dapdp=rho0*this%liq%cv*cvG*(GammaG-this%liq%gamma)
      dbpdp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(GammaG-1.0_WP)*cvG)                                     +&
      &     rho0*((GammaG-1.0_WP)*cvG*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qG)                        +&
      &     cvG*(2.0_WP*p_eq+this%liq%pinf)-this%liq%cv*(2.0_WP*p_eq+this%liq%gamma*this%liq%pinf)
      ddpdp=(qG-this%liq%q)*(2.0_WP*p_eq+this%liq%pinf)
   end subroutine relax_sg_ig_get_coeffs_lv

end module relax_sg_ig_class
