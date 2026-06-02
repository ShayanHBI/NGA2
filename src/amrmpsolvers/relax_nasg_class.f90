!> NASG liquid EOS relaxation model.
!> relaxation algorithm for NASG liquid + ideal gas mixture.
module relax_nasg_class
   use precision,   only: WP
   use relax_class, only: relax
   use nasg_class,  only: nasg
   use igmix_class, only: igmix
   implicit none
   private

   public :: relax_nasg

   real(WP), parameter :: Mv=0.0180153_WP   !< Molar mass of vapor [kg/mol]
   real(WP), parameter :: Ma=0.02897_WP     !< Molar mass of air [kg/mol]

   type, extends(relax) :: relax_nasg
      !> Typed EOS pointers
      type(nasg),   pointer :: liq => null()
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
      procedure :: initialize=>relax_nasg_initialize
      procedure :: apply     =>relax_nasg_apply
   end type relax_nasg

contains

   !> Initialize: store typed EOS pointers and compute saturation curve coefficients
   subroutine relax_nasg_initialize(this,liq,gas,indV,indA)
      class(relax_nasg), intent(inout)  :: this
      type(nasg),   target, intent(in)  :: liq
      class(igmix), target, intent(in)  :: gas
      integer, intent(in) :: indV,indA
      real(WP) :: cpV,cvV,RV
      this%liq =>liq
      this%gas =>gas
      this%indV= indV
      this%indA= indA
      cpV=gas%get_species_cp(indV)
      cvV=gas%get_species_cv(indV)
      RV =cpV-cvV
      this%AS=(liq%cp-cpV+gas%get_species_qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%get_species_q(indV))            /RV
      this%CS=(cpV-liq%cp)                                /RV
      this%DS=(liq%cp-liq%cv)                             /RV
      this%ES=liq%b                                       /RV
   end subroutine relax_nasg_initialize

   subroutine relax_nasg_apply(this,VF,Q,Pjump)
      class(relax_nasg), intent(inout)       :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable    :: Q0,y
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,PL,PG,TL,TG,IL,IG,ZL,ZG,PHIL,PHIG,Pint
      real(WP) :: xiL,xiG,xiLinv,xiGinv
      real(WP) :: xiTL,xiTG,xiTLinv,xiTGinv,zetaL,zetaG
      real(WP) :: Z,D,COF
      real(WP) :: VFeq,VF0,Peq,p,T,Teq,Yv
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: gammaG,cvG,qG,cpG
      real(WP) :: rho_err,rhoe_err
      real(WP), parameter :: p_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,pv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
      logical :: chem_relax
      ! Store the input state
      allocate(Q0(size(Q)))
      VF0=VF
      Q0=Q
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! ================ First step for mechanical relaxation ================
      ! Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
      ! Get phasic thermodynamic quantities
      RHOL=Q(1)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF)
      IL=Q(3)/Q(1)
      IG=Q(4)/Q(2)
      PL=this%liq%get_p_from_rho_e(rho=RHOL,e=IL)
      PG=this%gas%get_p_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=PL,rho=RHOL)
      CG=this%gas%get_c_from_p_rho(p=PG,rho=RHOG,y=y)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(8)=Yv*Q(2)
         call dealloc()
         return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP
         call dealloc()
         return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*CL
      ZG=Q(2)/(1.0_WP-VF)*CG
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup the ODE coefficients
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL)
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      xiL=         VF/(GL*(Pint-PL)+RHOL*CL**2)
      xiG=(1.0_WP-VF)/(GG*(Pint-PG)+RHOG*CG**2)
      xiLinv=1.0_WP/xiL
      xiGinv=1.0_WP/xiG
      ! Get equilibrium volume fraction
      VFeq=VF-(PG-PL)/(xiLinv+xiGinv)
      if ((VFeq.lt.0.0_WP).or.(VFeq.gt.1.0_WP)) then
         call restore()
         call dealloc()
         return
      end if
      ! Get equilibrium pressure
      Peq=get_p_eq(VFeq)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         call restore()
         call dealloc()
         return
      end if
      ! Adjust densities
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL)
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! ================= Second step for thermal relaxation =================
      ! Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
      ! Update thermodynamic state
      TL=this%liq%get_T_from_p_rho(p=Peq,rho=RHOL)
      TG=this%gas%get_T_from_p_rho(p=Peq,rho=RHOG,y=y)
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL)
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=Peq,rho=RHOL)
      CG=this%gas%get_c_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! print '(A)',       '============ P_relax ============='
      ! print '(A,ES15.7)','p   =',Peq
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A,ES15.7)','TL  =',TL
      ! print '(A,ES15.7)','TG  =',TG
      ! print '(A)',       '==================================='
      ! Setup the ODE coefficients
      Z=(1.0_WP-VF)*GL+VF*GG
      D=VF*RHOG*CG**2+(1.0_WP-VF)*RHOL*CL**2
      PHIL=-(this%liq%gamma-1.0_WP)*this%liq%cv*RHOL**2/(Peq+this%liq%pinf)
      PHIG=-(gammaG-1.0_WP)*cvG*RHOG**2/Peq
      zetaL=RHOL*(1.0_WP-this%liq%b*RHOL)/(Peq+this%liq%pinf)
      zetaG=RHOG/Peq
      COF=GL*RHOG*CG**2-GG*RHOL*CL**2
      xiTL=-PHIL*D/(RHOL/(       VF)*Z+zetaL*COF)
      xiTG=-PHIG*D/(RHOG/(1.0_WP-VF)*Z-zetaG*COF)
      xiTLinv=1.0_WP/xiTL
      xiTGinv=1.0_WP/xiTG
      ! Get equilibrium VF and T
      VFeq=VF+Z/D*(TG-TL)/(xiTLinv+xiTGinv)
      Teq=(xiTL*TL+xiTG*TG)/(xiTL+xiTG)
      ! Get equilibrium pressure
      Peq=get_p_eq(VFeq)
      ! Check if pressure is sound
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         call dealloc()
         return
      end if
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then
         VFeq=0.0_WP
         Peq=max(Peq,-this%liq%pinf)
      end if
      if (VFeq.gt.1.0_WP) then
         VFeq=1.0_WP
         Peq=max(Peq,0.0_WP)
      end if
      ! Adjust densities
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL)
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! ================= Third step for chemical relaxation =================
      ! Store input state to the chemical relaxation step
      VF0=VF
      Q0=Q
      p=Peq
      T=Teq
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
      ! Solve chemical equilibrium
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
      ! Apply the converged state
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
      !> Energy-conserving equilibrium pressure at given VF_ (Used for p and pT relaxation steps only)
      real(WP) function get_p_eq(VF_)
         real(WP), intent(in) :: VF_
         real(WP) :: one_brho
         one_brho=1.0_WP-this%liq%b*Q0(1)/VF_
         get_p_eq=(sum(Q0(3:4))-Q0(1)*this%liq%q-one_brho*VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0(2)*qG)/(one_brho*VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG-1.0_WP))
      end function get_p_eq
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
      !> Vapor mole fraction from mass fraction
      real(WP) function get_xv(Yv_)
         real(WP), intent(in) :: Yv_
         get_xv=Yv_*Ma/(Yv_*Ma+(1.0_WP-Yv_)*Mv)
      end function get_xv
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
         xv=get_xv(Yv_)
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
               ! pv_sat >= p: flash regime. Seed with a small Yv so that get_T stays well-conditioned (energy balance requires Yv << 1).
               Yv_=sqrt(0.0001_WP)
            else
               xv=pv_/p_
               Yv_=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
            end if
            Yv_=max(Yvmin,min(Yvmax,Yv_))
            ! print '(A,ES15.7)','Seeded Yv from saturation at Teq=',Yv_
         else
            ! Get saturation temperature from total pressure and vapor partial pressure.
            ! Use a safeguarded Newton solve rather than trusting T_ as the initial guess.
            ! At high-T metastable states, T_ can be far above the saturation temperature corresponding to (p_,pv_).
            call get_Tsat(p_,pv_,T_,Tsat,conv,Tsat_it)
            if (.not.conv) then
               ! print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
               return
            end if
            ! print '(A)',       '========== Finding Tsat ==========='
            ! print '(A,I2)','Tsat it= ',Tsat_it
            ! print '(A,ES15.7)','Tsat   =',Tsat
            ! print '(A)',       '==================================='
            ! Activate chemical relaxation only for metastable states
            if (T_.le.Tsat) return
         end if
         activate_chem=.true.
      end function activate_chem
      !> Safeguarded Newton-Raphson for saturation temperature at fixed (pl_, pv_)
      subroutine get_Tsat(pl_,pv_,Tguess,Tsat,conv,Tsat_it)
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
         Flo=PTsat(pl_,pv_,Tlo)
         Fhi=PTsat(pl_,pv_,Thi)
         ! Expand the bracket if needed.  In the normal liquid-vapor range PTsat is monotone in T, so these one-sided expansions are enough.
         expand_it=0
         do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
            if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
               Tlo=max(1.0_WP,0.8_WP*Tlo)
               Flo=PTsat(pl_,pv_,Tlo)
            else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
               Thi=1.2_WP*Thi
               Fhi=PTsat(pl_,pv_,Thi)
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
            Fold=PTsat(pl_,pv_,Told)
            dFold=dPTsatdT(pl_,Told)
            if (abs(dFold).gt.tiny(1.0_WP)) then
               Tnew=Told-Fold/dFold
            else
               Tnew=0.5_WP*(Tlo+Thi)
            end if
            ! Safeguard: reject Newton steps that leave the bracket or are NaN.
            if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) then
               Tnew=0.5_WP*(Tlo+Thi)
            end if
            Fnew=PTsat(pl_,pv_,Tnew)
            if (Fnew.ne.Fnew) then
               Tnew=0.5_WP*(Tlo+Thi)
               Fnew=PTsat(pl_,pv_,Tnew)
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
      !> p-T saturation curve
      real(WP) function PTsat(pl_,pv_,T_)
         real(WP), intent(in) :: pl_,pv_,T_
         PTsat=this%AS+(this%BS+this%ES*pl_)/T_+this%CS*log(T_)+this%DS*log(pl_+this%liq%pinf)-log(pv_)
      end function PTsat
      !> d(PTsat)/dT
      real(WP) function dPTsatdT(pl_,T_)
         real(WP), intent(in) :: pl_,T_
         dPTsatdT=-(this%BS+this%ES*pl_)/T_**2+this%CS/T_
      end function dPTsatdT
      !> d(PTsat)/dp for the LV solve (pv_ = pl_)
      real(WP) function dPTsatdp_lv(p_,T_,dTdp)
         real(WP), intent(in) :: p_,T_,dTdp
         dPTsatdp_lv=dPTsatdT(p_,T_)*dTdp+this%ES/T_+this%DS/(p_+this%liq%pinf)-1.0_WP/p_
      end function dPTsatdp_lv
      !> Equilibrium T from energy conservation
      real(WP) function get_T_lvg(p_,Yv_)
         real(WP), intent(in) :: p_,Yv_
         get_T_lvg=(1.0_WP-Yv_-this%liq%b*(rho0*(1.0_WP-Yv_)-rhoA0))/                                                                                                      &
         &         ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)+                                                                       &
         &           rhoA0*((this%gas%get_species_gamma(this%indV)-1.0_WP)*this%gas%get_species_cv(this%indV)*Yv_+                                                            &
         &                  (this%gas%get_species_gamma(this%indA)-1.0_WP)*this%gas%get_species_cv(this%indA)*(1.0_WP-Yv_))/p_)
      end function get_T_lvg
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
         rhoe_res_lv=(rho0-rho_g)/(rho_l-rho_g)*this%liq%get_e_from_p_T(p=p_,T=T_)+                  &
         &           (rho_l-rho0)/(rho_l-rho_g)*this%gas%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP,0.0_WP])-&
         &            rhoe0
      end function rhoe_res_lv
      !> Update the quadratic coefficients of equilibrium temperature equation
      subroutine get_coeffs_lv(p_,ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in)  :: p_
         real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: cvV_,gammaV_,qV_
         cvV_   =this%gas%get_species_cv(this%indV)
         gammaV_=this%gas%get_species_gamma(this%indV)
         qV_    =this%gas%get_species_q(this%indV)
         ! Coefficients
         ap=sum(Q(1:2))*this%liq%cv*cvV_*((gammaV_-this%liq%gamma)*p_+this%liq%gamma*(gammaV_-1.0_WP)*this%liq%pinf)
         bp=(cvV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-this%liq%cv)*p_**2+                                                      &
         &  (this%liq%pinf*(cvV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-                                                          &
         &   this%liq%gamma*this%liq%cv)+sum(Q(1:2))*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
         &   sum(Q(3:4))*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_+                                  &
         &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(sum(Q(1:2))*this%liq%q-sum(Q(3:4)))
         dp=p_*(p_+this%liq%pinf)*(qV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-this%liq%q+this%liq%b*sum(Q(3:4)))
         ! Pressure derivative of the coefficients
         dapdp=sum(Q(1:2))*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
         dbpdp=2.0_WP*(cvV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-this%liq%cv)*p_+                               &
         &     this%liq%pinf*(cvV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-this%liq%gamma*this%liq%cv)+            &
         &     sum(Q(1:2))*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+    &
         &     sum(Q(3:4))*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)
         ddpdp=(2.0_WP*p_+this%liq%pinf)*(qV_*(1.0_WP-sum(Q(1:2))*this%liq%b)-this%liq%q+this%liq%b*sum(Q(3:4)))
      end subroutine get_coeffs_lv
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
            call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Get temperature
            T_eq=get_T_lv(ap,bp,dp)
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Newton-Raphson iteration: Pure-vapor branch: Y_v=1 so the vapor mole fraction x_v=1.
            pOld=p_eq
            p_eq=pOld-PTsat(pOld,pOld,T_eq)/dPTsatdp_lv(pOld,T_eq,dTdp)
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
         call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
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
            T_eq=get_T_lvg(p_eq,Yv_eq)
            ! Get vapor partial pressure and mole fraction
            xv=get_xv(Yv_eq)
            pv=xv*p_eq
            if (.not.check_pv(pv)) return
            ! Get residuals at current state
            F1=PTsat(p_eq,pv,T_eq)/p_eq
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            res0=sqrt(F1**2+F2**2)
            ! Perturbation in p
            p_pert=p_eq*(1.0_WP+fd_eps)
            ppv_pert=xv*p_pert
            if (.not.check_pv(ppv_pert)) return
            ! Get the corresponding T
            T_pert=get_T_lvg(p_pert,Yv_eq)
            ! Get residuals
            F1p=PTsat(p_pert,ppv_pert,T_pert)/p_eq
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)/rhoe0
            dF1dp=(F1p-F1)/(p_pert-p_eq)
            dF2dp=(F2p-F2)/(p_pert-p_eq)
            ! Perturbation in Yv
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yvmax-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            xv_pert=get_xv(Yv_pert)
            ppv_pert=xv_pert*p_eq
            if (.not.check_pv(ppv_pert)) return
            ! Get the corresponding T
            T_pert=get_T_lvg(p_eq,Yv_pert)
            ! Get residuals
            F1Y=PTsat(p_eq,ppv_pert,T_pert)/p_eq
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
                  T_try=get_T_lvg(p_try,Yv_try)
                  if (T_try.le.0.0_WP) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  xv_try=get_xv(Yv_try)
                  ppv_try=xv_try*p_try
                  if (.not.check_pv(ppv_try)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  F1_try=PTsat(p_try,ppv_try,T_try)/p_try
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
            xv=get_xv(Yv_eq)
            pv=xv*p_eq
            ! Get temperature
            T_eq=get_T_lvg(p_eq,Yv_eq)
            ! Evaluate residuals
            F1=PTsat(p_eq,pv,T_eq)/p_eq
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
   end subroutine relax_nasg_apply

end module relax_nasg_class
