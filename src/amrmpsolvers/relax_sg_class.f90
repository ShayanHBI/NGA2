!> SG liquid EOS relaxation model.
!> Implements the 3-step pressure-temperature-chemical relaxation
!> algorithm for SG liquid + ideal-gas mixture (PinfG=0).
module relax_sg_class
   use precision,   only: WP
   use relax_class, only: relax
   use sg_class,    only: sg
   use igmix_class, only: igmix
   implicit none
   private

   public :: relax_sg

   real(WP), parameter :: Mv=0.0180153_WP   !< Molar mass of vapor [kg/mol]
   real(WP), parameter :: Ma=0.02897_WP     !< Molar mass of air [kg/mol]

   type, extends(relax) :: relax_sg
      !> Typed EOS pointers
      type(sg),     pointer :: liq => null()   !< Liquid EOS (SG)
      class(igmix), pointer :: gas => null()   !< Gas mixture (igmix)
      !> Species indices for vapor and air/carrier in the gas mixture
      integer :: indV=0,indA=0
      !> Saturation curve coefficients (Clausius-Clapeyron, SG form: ES=0)
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
      procedure :: initialize => relax_sg_initialize
      procedure :: apply      => relax_sg_apply
   end type relax_sg

contains

   !> Initialize: store typed EOS pointers and compute saturation curve coefficients
   subroutine relax_sg_initialize(this,liq,gas,indV,indA)
      class(relax_sg), intent(inout)   :: this
      type(sg),     target, intent(in) :: liq
      class(igmix), target, intent(in) :: gas
      integer, intent(in) :: indV,indA
      real(WP) :: CpV,CvV,RV
      this%liq  => liq
      this%gas  => gas
      this%indV =  indV
      this%indA =  indA
      CpV=gas%get_species_cp(indV)
      CvV=gas%get_species_cv(indV)
      RV =CpV-CvV
      this%AS=(liq%cp-CpV+gas%get_species_qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%get_species_q(indV))            /RV
      this%CS=(CpV-liq%cp)                                /RV
      this%DS=(liq%cp-liq%cv)                             /RV
      this%ES=0.0_WP
   end subroutine relax_sg_initialize

   !> Apply the 3-step SG relaxation to a single mixture cell
   subroutine relax_sg_apply(this,VF,Q,Pjump)
      class(relax_sg), intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: Q0,y
      real(WP) :: PL,PG,ZL,ZG,Pint
      real(WP) :: a,b,d,coeffL,coeffG
      real(WP) :: VFeq,VF0,Peq,p,T,Teq,Yv,RHOL,RHOG
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: GammaG,CvG,qG,CpG
      real(WP) :: rho_err,rhoe_err
      real(WP), parameter :: p_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,ppv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
      logical :: chem_relax
      ! Set gas EoS coefficients from vapor mass fraction
      allocate(y(this%gas%ns))
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=CvG,cp=CpG,q=qG,gamma=GammaG)
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=this%liq%get_p_from_rho_e(Q(1)/(       VF),Q(3)/Q(1))
      PG=this%gas%get_p_from_rho_e(Q(2)/(1.0_WP-VF),Q(4)/Q(2),y)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         return
      end if
      ! Get phasic impedances (use C^2 for the SG quadratic formulation)
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(PL,Q(1)/(       VF))**2
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(PG,Q(2)/(1.0_WP-VF),y)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(this%liq%gamma-1.0_WP)*Pint+2.0_WP*this%liq%gamma*this%liq%pinf
      coeffG=(GammaG        -1.0_WP)*Pint   ! PinfG=0 for ideal gas
      a=1.0_WP+GammaG*VF+this%liq%gamma*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+GammaG)*VF*PL-(1.0_WP+this%liq%gamma)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=VF*((this%liq%gamma-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+this%liq%gamma)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      ! ================= Second step for thermal relaxation =================
      ! etaL=qL, etaG=qG, PinfG=0 (ideal gas)
      ! Setup quadratic problem
      a=Q(1)*this%liq%cv+Q(2)*CvG
      b=this%liq%q*this%liq%cv*(this%liq%gamma-1.0_WP)*Q(1)**2 &
      &+qG*CvG*(GammaG-1.0_WP)*Q(2)**2 &
      &+Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf &
      &+Q(1)*Q(2)*(this%liq%q*CvG*(GammaG-1.0_WP)+qG*this%liq%cv*(this%liq%gamma-1.0_WP)) &
      &-sum(Q(3:4))*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*CvG*(GammaG-1.0_WP))
      d=CvG*(GammaG-1.0_WP)*this%liq%pinf &
      & *(qG*Q(2)**2+this%liq%q*Q(1)*Q(2)-sum(Q(3:4))*Q(2))
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) return
      ! Get equilibrium volume fraction
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq &
          /(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*CvG*(GammaG-1.0_WP)*(Peq+this%liq%pinf))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then
         VFeq=0.0_WP
         Peq=max(Peq,-this%liq%pinf)
      end if
      if (VFeq.gt.1.0_WP) then
         VFeq=1.0_WP
         Peq=max(Peq,0.0_WP)
      end if
      ! Get the thermo-mechanically relaxed temperature
      if (VFeq.gt.0.0_WP) then
         Teq=this%liq%get_T_from_p_rho(Peq,Q(1)/VFeq)
      else
         Teq=this%gas%get_T_from_p_rho(Peq,Q(2)/(1.0_WP-VFeq),y)
      end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(Peq,Q(1)/max(VFeq,tiny(1.0_WP)))
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(Peq,Q(2)/max(1.0_WP-VFeq,tiny(1.0_WP)),y)
      VF=VFeq
      ! ================= Third step for chemical relaxation =================
      ! Store input state to the chemical relaxation algorithm
      allocate(Q0(size(Q)))
      Q0=Q
      VF0=VF
      p=Peq
      T=Teq
      rho0=sum(Q0(1:2))
      rhoe0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! print '(A)',       '============ PT_relax ============='
      ! print '(A,ES15.7)','p    = ',p
      ! print '(A,ES15.7)','VF   = ',VF
      ! print '(A,ES15.7)','T    = ',T
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
         call this%gas%get_mix_coeffs(y=y,cv=CvG,cp=CpG,q=qG,gamma=GammaG)
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
      ! print '(A,ES15.7)','p    = ',p
      ! print '(A,ES15.7)','T    = ',T
      ! print '(A,ES15.7)','Yv   = ',Yv
      ! print '(A,ES15.7)','VF   = ',VF
      ! print '(A)',       '=================================='
      ! Apply the converged state
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=CvG,cp=CpG,q=qG,gamma=GammaG)
      RHOL=this%liq%get_rho_from_p_T(p,T)
      RHOG=this%gas%get_rho_from_p_T(p,T,y)
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
      ! Adjust conserved quantities
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*this%liq%get_e_from_p_T(p,T)
      Q(4)=Q(2)*this%gas%get_e_from_p_T(p,T,y)
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
      !> Sanity check pressure value
      logical function check_p(pv)
         real(WP), intent(in) :: pv
         check_p=(pv.gt.p_eps).and.((pv+this%liq%pinf).gt.p_eps)
      end function check_p
      !> Vapor mole fraction (Ideal gas mixture)
      real(WP) function get_xv(Yv_loc)
         real(WP), intent(in) :: Yv_loc
         get_xv=Yv_loc*Ma/(Yv_loc*Ma+(1.0_WP-Yv_loc)*Mv)
      end function get_xv
      !> Check conservation
      logical function check_cons()
         real(WP) :: re,ee
         re=(sum(Q(1:2))-rho0)/rho0
         ee=(sum(Q(3:4))-rhoe0)/rhoe0
         check_cons=(abs(re).le.this%rho_tol).and.(abs(ee).le.this%rhoe_tol)
      end function check_cons
      !> Decide if chemical relaxation needs to be activated and adjust Yv if needed
      logical function activate_chem(p_eq,T_eq,Yv_eq)
         real(WP), intent(in)    :: p_eq,T_eq
         real(WP), intent(inout) :: Yv_eq
         real(WP) :: xv,ppv,Tsat
         integer  :: Tsat_it
         logical  :: conv
         activate_chem=.false.
         ! Get vapor mole fraction and partial pressure
         xv=get_xv(Yv_eq)
         ppv=xv*p_eq
         if ((Yv_eq.le.Yv_dry).or.(ppv.le.ppv_dry).or.(.not.check_p(ppv))) then
            ! Dry/nearly-dry air edge case: ppv is zero or so tiny that
            ! solving Tsat(ppv) is log-singular/ill-conditioned.  Seed Yv
            ! from saturation at the current thermally-relaxed state and
            ! then continue with the ordinary LVG Newton solve.
            ! SG saturation seed: no ES*p term
            ppv=exp(this%AS+this%BS/T_eq+this%CS*log(T_eq)+this%DS*log(p_eq+this%liq%pinf))
            if (.not.check_p(ppv)) then
               ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
               return
            end if
            if (ppv.ge.p_eq) then
               ! ppv_sat >= p: flash regime. Seed with a small Yv so that
               ! get_T stays well-conditioned (energy balance requires Yv << 1).
               Yv_eq=sqrt(0.0001_WP)
            else
               xv=ppv/p_eq
               Yv_eq=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
            end if
            Yv_eq=max(Yvmin,min(Yvmax,Yv_eq))
            ! print '(A,ES15.7)','Seeded Yv from saturation at Teq = ',Yv_eq
         else
            ! Get saturation temperature from total pressure and vapor partial
            ! pressure.  Use a safeguarded Newton solve rather than trusting T_eq
            ! as the initial guess.  At high-T metastable states, T_eq can be far
            ! above the saturation temperature corresponding to (p_eq,ppv).
            call get_Tsat(p_eq,ppv,T_eq,Tsat,conv,Tsat_it)
            if (.not.conv) then
               ! print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
               return
            end if
            ! print '(A)',       '========== Finding Tsat ==========='
            ! print '(A,I2)','Tsat it= ',Tsat_it
            ! print '(A,ES15.7)','Tsat    = ',Tsat
            ! print '(A)',       '==================================='
            ! Activate chemical relaxation only for metastable states
            if (T_eq.le.Tsat) return
         end if
         activate_chem=.true.
      end function activate_chem
      !> Safeguarded Newton solve for saturation temperature at fixed (p_l, p_v)
      !> The full multispecies saturation relation g_l(p_l,T) = g_v(p_v,T) is a
      !> function of two pressures, not one.  We hold both p_l and p_v fixed and
      !> solve for T.  In the pure-vapor branch p_v = p_l and the result reduces
      !> to the classical Tsat(p) curve.
      subroutine get_Tsat(p_l,p_v,Tguess,Tsat,conv,Tsat_it)
         real(WP), intent(in)  :: p_l,p_v,Tguess
         real(WP), intent(out) :: Tsat
         logical,  intent(out) :: conv
         integer,  intent(out) :: Tsat_it
         real(WP) :: Tlo,Thi,Told,Tnew,Flo,Fhi,Fold,Fnew,dFold
         integer  :: it,expand_it
         conv=.false.
         Tsat_it=0
         ! Broad physical bracket for the current water EOS.  Newton is
         ! accepted only when it remains inside this bracket; otherwise the
         ! update falls back to bisection.
         Tlo=250.0_WP
         Thi=900.0_WP
         Flo=PTsat(p_l,p_v,Tlo)
         Fhi=PTsat(p_l,p_v,Thi)
         ! Expand the bracket if needed.  In the normal liquid-vapor range
         ! PTsat is monotone in T, so these one-sided expansions are enough.
         expand_it=0
         do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
            if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
               Tlo=max(1.0_WP,0.8_WP*Tlo)
               Flo=PTsat(p_l,p_v,Tlo)
            else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
               Thi=1.2_WP*Thi
               Fhi=PTsat(p_l,p_v,Thi)
            else
               exit
            end if
            expand_it=expand_it+1
         end do
         if (Flo*Fhi.gt.0.0_WP) then
            ! print*,'****************** Could not bracket Tsat!',p_l,p_v,Flo,Fhi,Tlo,Thi
            return
         end if
         ! Use the caller's guess only after clamping it to the safe bracket.
         Tsat=max(Tlo,min(Thi,Tguess))
         do it=1,this%Tsat_itmax
            Told=Tsat
            Fold=PTsat(p_l,p_v,Told)
            dFold=dPTsatdT(Told)
            if (abs(dFold).gt.tiny(1.0_WP)) then
               Tnew=Told-Fold/dFold
            else
               Tnew=0.5_WP*(Tlo+Thi)
            end if
            ! Safeguard: reject Newton steps that leave the bracket or are NaN.
            if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) then
               Tnew=0.5_WP*(Tlo+Thi)
            end if
            Fnew=PTsat(p_l,p_v,Tnew)
            if (Fnew.ne.Fnew) then
               Tnew=0.5_WP*(Tlo+Thi)
               Fnew=PTsat(p_l,p_v,Tnew)
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
      end subroutine get_Tsat
      !> Function that defines p-T saturation curve (SG: no ES*p/T term)
      real(WP) function PTsat(p_l,p_v,T_eq)
         real(WP), intent(in) :: p_l,p_v,T_eq
         PTsat=this%AS+this%BS/T_eq+this%CS*log(T_eq)+this%DS*log(p_l+this%liq%pinf)-log(p_v)
      end function PTsat
      !> Temperature derivative of p-T saturation curve function
      real(WP) function dPTsatdT(T_eq)
         real(WP), intent(in) :: T_eq
         dPTsatdT=-this%BS/T_eq**2+this%CS/T_eq
      end function dPTsatdT
      !> Pressure derivative of p-T saturation curve function for pure vapor p iteration
      real(WP) function dPTsatdp_lv(p_eq,T_eq,dTdp)
         real(WP), intent(in) :: p_eq,T_eq,dTdp
         dPTsatdp_lv=dPTsatdT(T_eq)*dTdp+this%DS/(p_eq+this%liq%pinf)-1.0_WP/p_eq
      end function dPTsatdp_lv
      !> Equilibrium temperature as a function of pressure and vapor mass fraction (SG form, no bL correction)
      real(WP) function get_T_lvg(p_eq,Yv_eq)
         real(WP), intent(in) :: p_eq,Yv_eq
         get_T_lvg=(1.0_WP-Yv_eq) &
         &        /((rho0*(1.0_WP-Yv_eq)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_eq+this%liq%pinf) &
         &         +rhoA0*((this%gas%get_species_gamma(this%indV)-1.0_WP)*this%gas%get_species_cv(this%indV)*Yv_eq+   &
         &                 (this%gas%get_species_gamma(this%indA)-1.0_WP)*this%gas%get_species_cv(this%indA)*(1.0_WP-Yv_eq))/p_eq)
      end function get_T_lvg
      !> Equilibrium temperature as a function of equilibrium pressure
      real(WP) function get_T_lv(ap,bp,dp)
         real(WP), intent(in) :: ap,bp,dp
         get_T_lv=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
      end function get_T_lv
      !> Pressure derivative of the equilibrium temperature as a function of equilibrium pressure
      real(WP) function get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         get_dTdp_lv=(ap*(-dbpdp+(bp*dbpdp-2.0_WP*(dapdp*dp+ap*ddpdp))/sqrt(bp**2-4.0_WP*ap*dp)) &
         &           -dapdp*(-bp+sqrt(bp**2-4.0_WP*ap*dp)))/(2.0_WP*ap**2)
      end function get_dTdp_lv
      !> Residual of the internal energy conservation equation for the LVG case
      real(WP) function rhoe_res_lvg(p_eq,T_eq,Yv_eq)
         real(WP), intent(in) :: p_eq,T_eq,Yv_eq
         rhoe_res_lvg=(rho0*(1.0_WP-Yv_eq)-rhoA0)*this%liq%get_e_from_p_T(p_eq,T_eq) &
         &           +rhoA0*this%gas%get_e_from_p_T(p_eq,T_eq,[Yv_eq,1.0_WP-Yv_eq]) &
         &           -rhoe0*(1.0_WP-Yv_eq)
      end function rhoe_res_lvg
      !> Residual of the internal energy conservation equation for the LV case
      real(WP) function rhoe_res_lv(p_eq,T_eq)
         real(WP), intent(in) :: p_eq,T_eq
         real(WP) :: rho_l,rho_g
         rho_l=this%liq%get_rho_from_p_T(p_eq,T_eq)
         rho_g=this%gas%get_rho_from_p_T(p_eq,T_eq,[1.0_WP,0.0_WP])
         rhoe_res_lv=(rho0-rho_g)/(rho_l-rho_g)*this%liq%get_e_from_p_T(p_eq,T_eq) &
         &         +(rho_l-rho0)/(rho_l-rho_g)*this%gas%get_e_from_p_T(p_eq,T_eq,[1.0_WP,0.0_WP]) &
         &         -rhoe0
      end function rhoe_res_lv
      !> (SG) Subroutine that updates the coefficients of the quadratic equilibrium temperature equation as functions of equilibrium pressure
      subroutine get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in)  :: p_eq
         real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         ! Coefficients
         ap=sum(Q(1:2))*this%liq%cv*CvG*((GammaG-1.0_WP)*(p_eq+this%liq%gamma*this%liq%pinf)-(this%liq%gamma-1.0_WP)*p_eq)
         bp=sum(Q(3:4))*((this%liq%gamma-1.0_WP)*this%liq%cv*p_eq-(GammaG-1.0_WP)*CvG*(p_eq+this%liq%pinf))  &
         &  +sum(Q(1:2))*((GammaG-1.0_WP)*CvG*this%liq%q*(p_eq+this%liq%pinf)-(this%liq%gamma-1.0_WP)*this%liq%cv*qG*p_eq) &
         &  +CvG*p_eq*(p_eq+this%liq%pinf)-this%liq%cv*p_eq*(p_eq+this%liq%gamma*this%liq%pinf)    ! PinfG=0
         dp=(qG-this%liq%q)*(p_eq+this%liq%pinf)*p_eq                                               ! PinfG=0
         ! Pressure derivative of the coefficients
         dapdp=sum(Q(1:2))*this%liq%cv*CvG*(GammaG-this%liq%gamma)
         dbpdp=sum(Q(3:4))*((this%liq%gamma-1.0_WP)*this%liq%cv-(GammaG-1.0_WP)*CvG)                           &
         &     +sum(Q(1:2))*((GammaG-1.0_WP)*CvG*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qG)            &
         &     +CvG*(2.0_WP*p_eq+this%liq%pinf)-this%liq%cv*(2.0_WP*p_eq+this%liq%gamma*this%liq%pinf)  ! PinfG=0
         ddpdp=(qG-this%liq%q)*(2.0_WP*p_eq+this%liq%pinf)                                           ! PinfG=0
      end subroutine get_coeffs_lv
      !> Pure liquid and vapor chemical relaxation (Pelanti and Shyue 2014)
      !> In the pure-vapor branch, Yv=1 so x_v=1 and p_v=p; pass p_eq for both
      !> the total-pressure and vapor-partial-pressure arguments of PTsat.
      subroutine solve_lv(p_eq,T_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq
         logical,  intent(out)   :: conv
         real(WP) :: pOld,ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: dTdp
         integer  :: it
         ! Iteratively solve for the equilibrium pressure in the pure vapor case
         conv=.false.
         do it=1,this%NR_itmax
            ! Get the coefficients
            call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Get temperature
            T_eq=get_T_lv(ap,bp,dp)
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Newton-Raphson iteration
            ! Pure-vapor branch: Y_v = 1 so the vapor mole fraction x_v = 1.
            pOld=p_eq
            p_eq=pOld-PTsat(pOld,pOld,T_eq)/dPTsatdp_lv(pOld,T_eq,dTdp)
            ! Evaluate the error
            if (abs((p_eq-pOld)/pOld).lt.this%p_tol) then
               conv=.true.
               exit
            end if
         end do
         if (.not.conv) then
            ! print*,"****************** p iterations blew up. Skipping the cell!!"
            return
         end if
         ! Update equilibrium temperature
         call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
         T_eq=get_T_lv(ap,bp,dp)
      end subroutine solve_lv
      !> Pure liquid and vapor-gas mixture chemical relaxation
      subroutine solve_lvg(p_eq,T_eq,Yv_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: conv
         real(WP) :: xv,ppv
         real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
         real(WP) :: p_pert,Yv_pert,T_pert,xv_pert,ppv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
         real(WP) :: pOld,YvOld,p_err,Yv_err
         real(WP) :: alpha,res0,res_try
         real(WP) :: p_try,Yv_try,T_try,xv_try,ppv_try,F1_try,F2_try
         integer  :: it
         logical  :: accepted
         ! Iteratively solve for the equilibrium pressure and vapor mass fraction
         conv=.false.
         p_err=10.0_WP*this%p_tol
         Yv_err=10.0_WP*this%Yv_tol
         do it=1,this%NR_itmax
            ! Evaluate residuals at current state
            T_eq=get_T_lvg(p_eq,Yv_eq)
            xv=get_xv(Yv_eq)
            ppv=xv*p_eq
            if (.not.check_p(ppv)) then
               ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
               return
            end if
            F1=PTsat(p_eq,ppv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)
            res0=sqrt(F1**2+(F2/rhoe0)**2)
            ! Compute Jacobian via finite difference
            ! --- Perturbation in p ---
            p_pert=p_eq*(1.0_WP+fd_eps)
            ppv_pert=xv*p_pert
            if (.not.check_p(ppv_pert)) then
               ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
               return
            end if
            T_pert=get_T_lvg(p_pert,Yv_eq)
            F1p=PTsat(p_pert,ppv_pert,T_pert)
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)
            dF1dp=(F1p-F1)/(p_pert-p_eq)
            dF2dp=(F2p-F2)/(p_pert-p_eq)
            ! --- Perturbation in Yv ---
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yvmax-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            xv_pert=get_xv(Yv_pert)
            ppv_pert=xv_pert*p_eq
            if (.not.check_p(ppv_pert)) then
               ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
               return
            end if
            T_pert=get_T_lvg(p_eq,Yv_pert)
            F1Y=PTsat(p_eq,ppv_pert,T_pert)
            F2Y=rhoe_res_lvg(p_eq,T_pert,Yv_pert)
            dF1dYv=(F1Y-F1)/(Yv_pert-Yv_eq)
            dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
            ! Solve 2x2 system: J * [dp; dYv] = -[F1; F2]
            detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
            if (abs(detJ).lt.1.0e-30_WP) exit
            dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
            dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
            ! Damped Newton-Raphson update
            pOld=p_eq
            YvOld=Yv_eq
            alpha=1.0_WP
            if ((abs(F1).lt.F_line_search_tol).and.(abs(F2/rhoe0).lt.F_line_search_tol)) then
               p_eq=pOld+dp_nr
               Yv_eq=YvOld+dYv_nr
            else
               accepted=.false.
               do while (alpha.gt.1.0e-8_WP)
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
                  T_try=get_T_lvg(p_try,Yv_try)
                  if ((T_try.ne.T_try).or.(T_try.le.0.0_WP)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  xv_try=get_xv(Yv_try)
                  ppv_try=xv_try*p_try
                  if (.not.check_p(ppv_try)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  F1_try=PTsat(p_try,ppv_try,T_try)
                  F2_try=rhoe_res_lvg(p_try,T_try,Yv_try)
                  res_try=sqrt(F1_try**2+(F2_try/rhoe0)**2)
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
            ! Refresh xv and ppv from the accepted (p_eq, Yv_eq) before re-evaluating
            ! the residuals.  The previous ppv was stale from before the line search.
            xv=get_xv(Yv_eq)
            ppv=xv*p_eq
            F1=PTsat(p_eq,ppv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)
            ! Per-iteration diagnostic (fires on EVERY iteration, including the
            ! one that converges, so we can always see the full trace).
            ! write(*,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
            !    '  it=',it,                                  &
            !    '  p=',p_eq,                                 &
            !    '  Yv=',Yv_eq,                               &
            !    '  T=',T_eq,                                 &
            !    '  alpha=',alpha,                            &
            !    '  |F1|=',abs(F1),                           &
            !    '  |F2/E0|=',abs(F2/rhoe0),                  &
            !    '  perr=',p_err,                             &
            !    '  Yverr=',Yv_err
            if ((p_err.lt.this%p_tol).and.(Yv_err.lt.this%Yv_tol) &
            &   .and.(abs(F1).lt.this%F1_tol).and.(abs(F2/rhoe0).lt.this%F2_tol)) then
               conv=.true.
               exit
            end if
         end do
         if (.not.conv) then
            ! print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
            return
         end if
         ! Update equilibrium temperature
         T_eq=get_T_lvg(p_eq,Yv_eq)
      end subroutine solve_lvg
   end subroutine relax_sg_apply

end module relax_sg_class
