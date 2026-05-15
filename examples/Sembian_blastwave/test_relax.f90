! Standalone single-cell test for PTg_relax
! Compile: gfortran -O2 -o test_relax test_relax.f90
! Run:     ./test_relax
module eos_mod
   implicit none
   integer, parameter :: WP = selected_real_kind(15,307)
   real(WP) :: GammaL,PinfL,etaL,etapL,CvL,CpL
   real(WP) :: GammaG,PinfG,etaG,etapG,CvG,CpG
   real(WP) :: GammaA,PinfA,etaA,etapA,CvA,CpA
   real(WP) :: GammaV,PinfV,etaV,etapV,CvV,CpV
   real(WP) :: AS,BS,CS,DS
   real(WP), parameter :: Mv=0.0180153_WP, Ma=0.02897_WP
contains
   !> Liquid EOS: P=f(RHO,I) - Stiffened gas
   real(WP) function get_PL(RHO,I,Yv)
      implicit none
      real(WP), intent(in) :: RHO,I
      real(WP), intent(in), optional :: Yv
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL-(GammaL-1.0_WP)*etaL*RHO
   end function get_PL
   !> Liquid EOS: T=f(RHO,P)
   real(WP) function get_TL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> Liquid EOS: C=f(RHO,P)
   real(WP) function get_CL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_CL=sqrt(max(0.0_WP,GammaL*(P+PinfL)/RHO))
   end function get_CL
   !> Liquid EOS: I=f(RHO,P) (used for initialization)
   real(WP) function get_IL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_IL=(P+GammaL*PinfL)/(RHO*(GammaL-1.0_WP))+etaL
   end function get_IL
   !> Gas EOS: I=f(P,T)
   real(WP) function get_IL_PT(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_IL_PT=CvL*T*(P+GammaL*PinfL)/(P+PinfL)+etaL
   end function get_IL_PT
   !> Liquid EOS: RHOLIL=f(P,T)
   real(WP) function get_RHOLIL(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_RHOLIL=(P+GammaL*PinfL+etaL*(P+PinfL)/(CvL*T))/(GammaL-1.0_WP)
   end function get_RHOLIL
   !> Liquid EOS: RHOL=f(P,T)
   real(WP) function get_RHOL(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_RHOL=(P+PinfL)/((GammaL-1.0_WP)*CvL*T)
   end function get_RHOL
   !> Gas EOS: P=f(RHO,I)
   real(WP) function get_PG(RHO,I,Yv)
      implicit none
      real(WP), intent(in) :: RHO,I
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG-(GammaG-1.0_WP)*etaG*RHO
   end function get_PG
   !> Gas EOS: T=f(RHO,P)
   real(WP) function get_TG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> Gas EOS: C=f(RHO,P)
   real(WP) function get_CG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_CG=sqrt(max(0.0_WP,GammaG*(P+PinfG)/RHO))
   end function get_CG
   !> Gas EOS: I=f(RHO,P) (used for initialization)
   real(WP) function get_IG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_IG=(P+GammaG*PinfG)/(RHO*(GammaG-1.0_WP))+etaG
   end function get_IG
   !> Gas EOS: I=f(P,T,Yv)
   real(WP) function get_IG_PT(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_IG_PT=CvG*T*(P+GammaG*PinfG)/(P+PinfG)+etaG
   end function get_IG_PT
   !> Gas EOS: RHOGIG=f(P,T,Yv)
   real(WP) function get_RHOGIG(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_RHOGIG=(P+GammaG*PinfG+etaG*(P+PinfG)/(CvG*T))/(GammaG-1.0_WP)
   end function get_RHOGIG
   !> Liquid EOS: RHOG=f(P,T,Yv)
   real(WP) function get_RHOG(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_RHOG=(P+PinfG)/((GammaG-1.0_WP)*CvG*T)
   end function get_RHOG
   !> Set gas EOS coefficients from vapor mass fraction (Ideal gas mixture of vapor and air)
   subroutine set_gas_eos_cof(Yv)
      implicit none
      real(WP), intent(in) :: Yv
      real(WP) :: Ya
      Ya    =1.0_WP-Yv
      CvG   =Yv*CvV+Ya*CvA
      CpG   =Yv*CpV+Ya*CpA
      etaG  =Yv*etaV +Ya*etaA
      etapG =Yv*etapV+Ya*etapA
      PinfG =0.0_WP
      GammaG=CpG/CvG
   end subroutine set_gas_eos_cof
   !> Vapor mole fraction (Ideal gas mixture)
   real(WP) function get_xv(Yv)
      implicit none
      real(WP), intent(in) :: Yv
      get_xv=Yv*Ma/(Yv*Ma+(1.0_WP-Yv)*Mv)
   end function get_xv
end module eos_mod

module relax_mod
   use eos_mod
   implicit none
contains
   subroutine PTg_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP), dimension(:), allocatable    :: Q0
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG
      real(WP) :: VFeq,VF0,Peq,p,T,Teq,Yv,RHOL,RHOG
      real(WP) :: rho0,Eps0,rhoA0
      real(WP) :: rho_err,Eps_err
      real(WP), parameter :: lnP_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,ppv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: p_tol=1.0e-4_WP,Yv_tol=1.0e-4_WP,Tsat_tol=1.0e-5_WP,rho_tol=1.0e-4_WP,Eps_tol=1.0e-4_WP,F1_tol=1e-7_WP,F2_tol=1e-7_WP
      integer,  parameter :: Tsat_itmax=40,NR_itmax=40
      logical  :: chem_relax
      ! Set gas EoS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      call set_gas_eos_cof(Yv)
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2),Yv=Yv)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; return
      end if
      if (PG.le.-PinfG) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)**2
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG,Yv=Yv)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(GammaL-1.0_WP)*Pint+2.0_WP*GammaL*PinfL
      coeffG=(GammaG-1.0_WP)*Pint+2.0_WP*GammaG*PinfG
      a=1.0_WP+GammaG*VF+GammaL*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+GammaG)*VF*PL-(1.0_WP+GammaL)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      ! ================= Second step for thermal relaxation =================
      ! Setup quadratic problem
      a=Q(1)*CvL+Q(2)*CvG
      b=etaL*CvL*(GammaL-1.0_WP)*Q(1)**2+etaG*CvG*(GammaG-1.0_WP)*Q(2)**2+&
      & Q(1)*CvL*(GammaL*PinfL+PinfG)+Q(2)*CvG*(GammaG*PinfG+PinfL)      +&
      & Q(1)*Q(2)*(etaL*CvG*(GammaG-1.0_WP)+etaG*CvL*(GammaL-1.0_WP))    -&
      & sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)+Q(2)*CvG*(GammaG-1.0_WP))
      d=etaL*CvL*(GammaL-1.0_WP)*PinfG*Q(1)**2+etaG*CvG*(GammaG-1.0_WP)*PinfL*Q(2)**2+&
      & (Q(1)*CvL*GammaL+Q(2)*CvG*GammaG)*PinfL*PinfG                                +&
      Q(1)*Q(2)*(etaL*CvG*(GammaG-1.0_WP)*PinfL+etaG*CvL*(GammaL-1.0_WP)*PinfG)      -&
      & sum(Q(3:4))*(Q(1)*CvL*(GammaL-1.0_WP)*PinfG+Q(2)*CvG*(GammaG-1.0_WP)*PinfL)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)/(Q(1)*CvL*(GammaL-1.0_WP)*(Peq+PinfG)+Q(2)*CvG*(GammaG-1.0_WP)*(Peq+PinfL))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-PinfL); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,-PinfG); end if
      ! Get the thermo-mechanically relaxed temperature
      if (VFeq.gt.0.0_WP) then
         Teq=get_TL(RHO=Q(1)/VFeq,P=Peq)
      else
         Teq=get_TG(RHO=Q(2)/(1.0_WP-VFeq),P=Peq)
      end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*get_RHOLIL(Peq,Teq)
      Q(4)=(1.0_WP-VFeq)*get_RHOGIG(Peq,Teq,Yv)
      VF=VFeq
      ! ================= Third step for chemical relaxation =================
      ! Store input state to the chemical relaxation algorithm
      allocate(Q0(size(Q)))
      Q0=Q
      VF0=VF
      p=Peq
      T=Teq
      rho0=sum(Q0(1:2))
      Eps0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! print '(A)',       '============ PT_relax ============='
      ! print '(A,ES15.7)','p    = ',p
      ! print '(A,ES15.7)','VF   = ',VF
      ! print '(A,ES15.7)','T    = ',T
      ! print '(A)',       '==================================='
      ! Check if chemical relaxation should be activated
      chem_relax=activate_chem_relax(p,T,Yv)
      if (.not.chem_relax) then
         call restore_VFQ()
         call deallocate_Q0()
         return
      end if
      ! Solve chemical equilibrium without touching the conserved variables.
      ! If the inert gas content is only a numerical trace, the LVG equations
      ! become ill-conditioned; use the pure liquid-vapor branch instead.
      if (Yv.gt.Yv_pure) then
         ! print*,'****************** Using pure LV chemical relaxation!'
         Yv=Yvmax
         call set_gas_eos_cof(Yv)
         call solve_PTg_relax_lv(p,T,chem_relax)
      else
         call solve_PTg_relax_lvg(p,T,Yv,chem_relax)
      end if
      ! Skip if not converged
      if (.not.chem_relax) then
         call restore_VFQ()
         call deallocate_Q0()
         return
      end if
      ! print '(A)',       '=========== PTg_relax ============'
      ! print '(A,ES15.7)','p    = ',p
      ! print '(A,ES15.7)','T    = ',T
      ! print '(A,ES15.7)','Yv   = ',Yv
      ! print '(A,ES15.7)','VF   = ',VF
      ! print '(A)',       '=================================='
      ! Apply the converged state
      call set_gas_eos_cof(Yv)
      RHOL=get_RHOL(p,T)
      RHOG=get_RHOG(p,T,Yv)
      VF=(rho0-RHOG)/(RHOL-RHOG)
      ! Clean up solution
      if (VF.lt.0.0_WP) then; VF=0.0_WP; p=max(p,-PinfL); end if
      if (VF.gt.1.0_WP) then; VF=1.0_WP; p=max(p,-PinfG); end if
      ! Adjust conserved quantities
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*get_IL_PT(p,T)
      Q(4)=Q(2)*get_IG_PT(p,T,Yv)
      Q(8)=Q(2)*Yv
      ! Evaluate conservation
      if (.not.check_cons()) then
         call restore_VFQ()
         call deallocate_Q0()
         ! print*,"****************** Conservation is violated. Skipping the cell!!"
         return
      end if
      ! Check VOF
      if (VF.lt.VFmin) then
         call restore_VFQ()
         call deallocate_Q0()
         ! print*,"****************** Not enough liquid to vaporize. Skipping the cell!!"
         return
      end if
      ! Release memory
      call deallocate_Q0()
      contains
         !> Reset the output to the initial values fed into chemical relaxation algorithm
         subroutine restore_VFQ()
            VF=VF0
            Q=Q0
         end subroutine restore_VFQ
         !> Release memory allocated for Q0
         subroutine deallocate_Q0()
            if (allocated(Q0)) deallocate(Q0)
         end subroutine deallocate_Q0
         !> Decide if chemical relaxation needs to be activated and adjust Yv if needed
         logical function activate_chem_relax(p_eq,T_eq,Yv_eq)
            real(WP), intent(in)    :: p_eq,T_eq
            real(WP), intent(inout) :: Yv_eq
            real(WP) :: xv,ppv,Tsat
            integer  :: Tsat_it
            logical  :: converge
            activate_chem_relax=.false.
            ! Get vapor mole fraction and partial pressure
            xv=get_xv(Yv_eq)
            ppv=xv*p_eq
            if ((Yv_eq.le.Yv_dry).or.(ppv.le.ppv_dry).or.(.not.check_p(ppv))) then
               if (PinfV.lt.lnP_eps) then
                  ! Dry/nearly-dry air edge case: ppv is zero or so tiny that
                  ! solving Tsat(ppv) is log-singular/ill-conditioned.  Seed Yv
                  ! from saturation at the current thermally-relaxed state and
                  ! then continue with the ordinary LVG Newton solve.
                  ppv=exp(AS+BS/T_eq+CS*log(T_eq)+DS*log(p_eq+PinfL))
                  if (.not.check_p(ppv)) then
                     ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                     return
                  end if
                  if (ppv.ge.p_eq) then
                     ! ppv_sat >= p: flash regime. Seed with a small Yv so that
                     ! get_T stays well-conditioned (energy balance requires Yv << 1).
                     Yv_eq=sqrt(0.0001)
                  else
                     xv=ppv/p_eq
                     Yv_eq=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
                  end if
                  Yv_eq=max(Yvmin,min(Yvmax,Yv_eq))
                  ! print '(A,ES15.7)','Seeded Yv from saturation at Teq = ',Yv_eq
               else
                  ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
            else
               ! Get saturation temperature from total pressure and vapor partial
               ! pressure.  Use a safeguarded Newton solve rather than trusting T_eq
               ! as the initial guess.  At high-T metastable states, T_eq can be far
               ! above the saturation temperature corresponding to (p_eq,ppv).
               call get_Tsat(p_eq,ppv,T_eq,Tsat,converge,Tsat_it)
               if (.not.converge) then
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
            activate_chem_relax=.true.
         end function activate_chem_relax
         !> Safeguarded Newton solve for saturation temperature at fixed (p_l, p_v)
         !> The full multispecies saturation relation g_l(p_l,T) = g_v(p_v,T) is a
         !> function of two pressures, not one.  We hold both p_l and p_v fixed and
         !> solve for T.  In the pure-vapor branch p_v = p_l and the result reduces
         !> to the classical Tsat(p) curve.
         subroutine get_Tsat(p_l,p_v,Tguess,Tsat,converged,Tsat_it)
            real(WP), intent(in)  :: p_l,p_v,Tguess
            real(WP), intent(out) :: Tsat
            logical,  intent(out) :: converged
            integer,  intent(out) :: Tsat_it
            real(WP) :: Tlo,Thi,Told,Tnew
            real(WP) :: Flo,Fhi,Fold,Fnew,dFold
            real(WP) :: Tsat_err
            integer  :: it,expand_it
            converged=.false.
            Tsat_it=0
            ! if (.not.check_p(ppv)) return
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
            do it=1,Tsat_itmax
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
               Tsat_err=abs((Tnew-Told)/max(abs(Told),tiny(1.0_WP)))
               Tsat=Tnew
               Tsat_it=it
               if ((Tsat_err.lt.Tsat_tol).or.(abs(Fnew).lt.F1_tol)) then
                  converged=.true.
                  return
               end if
            end do
         end subroutine get_Tsat
         !> Pure liquid and vapor chemical relaxation (Pelanti and Shyue 2014)
         !> In the pure-vapor branch, Yv=1 so x_v=1 and p_v=p; pass p_eq for both
         !> the total-pressure and vapor-partial-pressure arguments of PTsat.
         subroutine solve_PTg_relax_lv(p_eq,T_eq,converged)
            real(WP), intent(inout) :: p_eq,T_eq
            logical,  intent(out)   :: converged
            real(WP) :: pOld,p_err,dTdp
            real(WP) :: ap,bp,dp,dapdp,dbpdp,ddpdp
            integer  :: it
            ! Iteratively solve for the equilibrium pressure in the pure vapor case
            converged=.false.
            do it=1,NR_itmax
               ! Get the coefficients
               call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
               ! Get temperature
               T_eq=get_T_lv(ap,bp,dp)
               dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
               ! Newton-Raphson iteration
               ! Pure-vapor branch: Y_v = 1 so the vapor mole fraction x_v = 1.
               ! dPTsatdp takes (p_eq, x_v, T_eq, dTdp); pass 1 for x_v here.
               pOld=p_eq
               p_eq=pOld-PTsat(pOld,pOld,T_eq)/dPTsatdp(pOld,1.0_WP,T_eq,dTdp)
               ! Evaluate the error
               p_err=abs((p_eq-pOld)/pOld)
               if (p_err.lt.p_tol) then
                  converged=.true.
                  exit
               end if
            end do
            print*,'NR iterations = ',it
            if (.not.converged) then
               ! print*,"****************** p iterations blew up. Skipping the cell!!"
               return
            end if
            ! Update equilibrium temperature
            call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp)
         end subroutine solve_PTg_relax_lv
         !> Pure liquid and vapor-gas mixture chemical relaxation
         subroutine solve_PTg_relax_lvg(p_eq,T_eq,Yv_eq,converged)
            real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
            logical,  intent(out)   :: converged
            real(WP) :: xv,ppv
            real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
            real(WP) :: p_pert,Yv_pert,T_pert,xv_pert,ppv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
            real(WP) :: pOld,YvOld,p_err,Yv_err
            real(WP) :: alpha,res0,res_try
            real(WP) :: p_try,Yv_try,T_try,xv_try,ppv_try,F1_try,F2_try
            real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
            integer  :: it
            logical  :: accepted
            ! Iteratively solve for the equilibrium pressure and vapor mass fraction
            converged=.false.
            p_err=10.0_WP*p_tol
            Yv_err=10.0_WP*Yv_tol
            do it=1,NR_itmax
               ! Evaluate residuals at current state
               T_eq=get_T_lvg(p_eq,Yv_eq)
               xv=get_xv(Yv_eq)
               ppv=xv*p_eq
               if (.not.check_p(ppv)) then
                  ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
               F1=PTsat(p_eq,ppv,T_eq)
               F2=Eps_res_lvg(p_eq,Yv_eq,T_eq)
               res0=sqrt(F1**2+(F2/Eps0)**2)
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
               F2p=Eps_res_lvg(p_pert,Yv_eq,T_pert)
               dF1dp=(F1p-F1)/(p_pert-p_eq)
               dF2dp=(F2p-F2)/(p_pert-p_eq)
               ! --- Perturbation in Yv ---
               Yv_pert=Yv_eq+fd_eps
               if (Yv_pert.gt.Yvmax-fd_eps) then
                  Yv_pert=Yv_eq-fd_eps
               end if
               if (Yv_pert.lt.Yvmin+fd_eps) then
                  Yv_pert=Yv_eq+fd_eps
               end if
               xv_pert=get_xv(Yv_pert)
               ppv_pert=xv_pert*p_eq
               if (.not.check_p(ppv_pert)) then
                  ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
               T_pert=get_T_lvg(p_eq,Yv_pert)
               F1Y=PTsat(p_eq,ppv_pert,T_pert)
               F2Y=Eps_res_lvg(p_eq,Yv_pert,T_pert)
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
               if ((abs(F1).lt.F_line_search_tol).and.(abs(F2/Eps0).lt.F_line_search_tol)) then
                  p_eq=pOld+dp_nr
                  Yv_eq=YvOld+dYv_nr
               else
                  accepted=.false.
                  do while (alpha.gt.1.0e-8_WP)
                     p_try=p_eq+alpha*dp_nr
                     Yv_try=Yv_eq+alpha*dYv_nr
                     ! Keep the trial state inside the physical/log-safe domain
                     if (p_try.le.lnP_eps) then
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
                     F2_try=Eps_res_lvg(p_try,Yv_try,T_try)
                     res_try=sqrt(F1_try**2+(F2_try/Eps0)**2)
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
               F2=Eps_res_lvg(p_eq,Yv_eq,T_eq)
               ! Per-iteration diagnostic (fires on EVERY iteration, including the
               ! one that converges, so we can always see the full trace).
               ! write(*,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
               !    '  it=',it,                                  &
               !    '  p=',p_eq,                                 &
               !    '  Yv=',Yv_eq,                               &
               !    '  T=',T_eq,                                 &
               !    '  alpha=',alpha,                            &
               !    '  |F1|=',abs(F1),                           &
               !    '  |F2/E0|=',abs(F2/Eps0),                   &
               !    '  perr=',p_err,                             &
               !    '  Yverr=',Yv_err
               if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol).and.(abs(F1).lt.F1_tol).and.(abs(F2/Eps0).lt.F2_tol)) then
                  converged=.true.
                  exit
               end if
            end do
            print*,'NR iterations = ',it
            if (.not.converged) then
               ! print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
               return
            end if
            ! Update equilibrium temperature
            T_eq=get_T_lvg(p_eq,Yv_eq)
         end subroutine solve_PTg_relax_lvg
         !> Check conservation
         logical function check_cons()
            rho_err=(sum(Q(1:2))-rho0)/rho0
            Eps_err=(sum(Q(3:4))-Eps0)/Eps0
            check_cons=(abs(rho_err).le.rho_tol).and.(abs(Eps_err).le.Eps_tol)
         end function check_cons
         !> Equilibrium temperature as a function of pressure and vapor mass fraction
         function get_T_lvg(p_eq,Yv_eq)
            real(WP), intent(in) :: p_eq,Yv_eq
            real(WP) :: get_T_lvg
            get_T_lvg=(1.0_WP-Yv_eq)/((rho0*(1.0_WP-Yv_eq)-rhoA0)*(GammaL-1.0_WP)*CvL/(p_eq+PinfL)+rhoA0*((GammaV-1.0_WP)*CvV*Yv_eq/(p_eq+PinfV)+(GammaA-1.0_WP)*CvA*(1.0_WP-Yv_eq)/(p_eq+PinfA)))
         end function get_T_lvg
         !> Function that defines p-T saturation curve
         function PTsat(p_l,p_v,T_eq)
            real(WP), intent(in) :: p_l,p_v,T_eq
            real(WP) :: PTsat
            PTsat=AS+BS/T_eq+CS*log(T_eq)+DS*log(p_l+PinfL)-log(p_v+PinfV)
         end function PTsat
         !> Temperature derivative of p-T saturation curve function
         function dPTsatdT(T_eq)
            real(WP), intent(in) :: T_eq
            real(WP) :: dPTsatdT
            dPTsatdT=-BS/T_eq**2+CS/T_eq
         end function dPTsatdT
         !> Pressure derivative of p-T saturation curve function for pure vapor p iteration
         function dPTsatdp(p_eq,x_v,T_eq,dTdp)
            real(WP), intent(in) :: p_eq,x_v,T_eq,dTdp
            real(WP) :: dPTsatdp
            dPTsatdp=dPTsatdT(T_eq)*dTdp+DS/(p_eq+PinfL)-x_v/(x_v*p_eq+PinfV)
         end function dPTsatdp
         !> Residual of the internal energy conservation equation
         function Eps_res_lvg(p_eq,Yv_eq,T_eq)
            real(WP), intent(in) :: p_eq,Yv_eq,T_eq
            real(WP) :: Eps_res_lvg
            Eps_res_lvg=(rho0*(1.0_WP-Yv_eq)-rhoA0)*(CvL*T_eq*(p_eq+GammaL*PinfL)/(p_eq+PinfL)+etaL)+rhoA0*((CvV*T_eq*(p_eq+GammaV*PinfV)/(p_eq+PinfV)+etaV)*Yv_eq+(CvA*T_eq*(p_eq+GammaA*PinfA)/(p_eq+PinfA)+etaA)*(1.0_WP-Yv_eq))-Eps0*(1.0_WP-Yv_eq)
         end function Eps_res_lvg
         !> Equilibrium temperature as a function of equilibrium pressure
         function get_T_lv(ap,bp,dp)
            real(WP), intent(in) :: ap,bp,dp
            real(WP) :: get_T_lv
            get_T_lv=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
         end function get_T_lv
         !> Pressure derivative of the equilibrium temperature as a function of equilibrium pressure
         function get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
            real(WP), intent(in) :: ap,bp,dp,dapdp,dbpdp,ddpdp
            real(WP) :: get_dTdp_lv
            get_dTdp_lv=(ap*(-dbpdp+(bp*dbpdp-2.0_WP*(dapdp*dp+ap*ddpdp))/sqrt(bp**2-4.0_WP*ap*dp))-dapdp*(-bp+sqrt(bp**2-4.0_WP*ap*dp)))/(2.0_WP*ap**2)
         end function get_dTdp_lv
         !> Subroutine that updates the coefficients of the quadradic equilibrium temperature equation as functions of equilibrium pressure
         subroutine get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            real(WP), intent(in)  :: p_eq
            real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
            ! Coefficients
            ap=sum(Q(1:2))*CvL*CvG*((GammaG-1.0_WP)*(p_eq+GammaL*PinfL)-(GammaL-1.0_WP)*(p_eq+GammaG*PinfG))
            bp=sum(Q(3:4))*((GammaL-1.0_WP)*CvL*(p_eq+PinfG)-(GammaG-1.0_WP)*CvG*(p_eq+PinfL))+&
            &  sum(Q(1:2))*((GammaG-1.0_WP)*CvG*etaL*(p_eq+PinfL)-(GammaL-1.0_WP)*CvL*etaG*(p_eq+PinfG))+&
            &  CvG*(p_eq+PinfL)*(p_eq+GammaG*PinfG)-CvL*(p_eq+PinfG)*(p_eq+GammaL*PinfL)
            dp=(etaG-etaL)*(p_eq+PinfL)*(p_eq+PinfG)
            ! Pressure derivative of the coefficients
            dapdp=sum(Q(1:2))*CvL*CvG*(GammaG-GammaL)
            dbpdp=sum(Q(3:4))*((GammaL-1.0_WP)*CvL-(GammaG-1.0_WP)*CvG)+&
            &     sum(Q(1:2))*((GammaG-1.0_WP)*CvG*etaL-(GammaL-1.0_WP)*CvL*etaG)+&
            &     CvG*(2.0_WP*p_eq+PinfL+GammaG*PinfG)-CvL*(2.0_WP*p_eq+PinfG+GammaL*PinfL)
            ddpdp=(etaG-etaL)*(2.0_WP*p_eq+PinfL+PinfG)
         end subroutine get_coeffs_lv
         !> Sanity check pressure value
         logical function check_p(p_v)
            real(WP), intent(in) :: p_v
            check_p=((p_v+PinfV).gt.lnP_eps).and.((p_v+PinfL).gt.lnP_eps)
         end function check_p
   end subroutine PTg_relax
end module relax_mod


program test_relax
   use eos_mod
   use relax_mod
   implicit none
   real(WP) :: VF,Q(8),Q0(8),VF0
   real(WP) :: T0,p0,rhoL,rhoG,Yv0

   ! Set EOS parameters (Pelanti & Shyue 2014 Table 1 + air)
   GammaL=2.35_WP;  PinfL=1.0e9_WP;  etaL=-1.167e6_WP;  etapL=0.0_WP;      CvL=1816.0_WP
   GammaV=1.43_WP;  PinfV=0.0_WP;    etaV=2030.0e3_WP;  etapV=-23.4e3_WP;  CvV=1040.0_WP
   GammaA=1.40_WP;  PinfA=0.0_WP;    etaA=0.0_WP;       etapA=0.0_WP;      CvA=718.0_WP
   CpL=GammaL*CvL;  CpA=GammaA*CvA;  CpV=GammaV*CvV

   ! Saturation curve coefficients (Pelanti & Shyue 2014, Eq. 11b)
   AS=(CpL-CpV+etapV-etapL)/(CpV-CvV)
   BS=(etaL-etaV)/(CpV-CvV)
   CS=(CpV-CpL)/(CpV-CvV)
   DS=(CpL-CvL)/(CpV-CvV)
   print '(A,4ES15.7)', 'Sat coeffs AS,BS,CS,DS = ',AS,BS,CS,DS

   ! --- Setup metastable liquid ---
   ! Very low VOF means there is not enough liquid to vaporize. The solver skips those cells for now. Later, we can think of clustering.
   ! Very low Yv is replaced by 1e-4 or from vapor partial pressure taken from saturation curve depending on which outcome is more physical.
   ! Very high Yv is treated as pure liquid pure vapor case
   ! p0=1e5_WP; T0=480.0_WP; Yv0=0.5_WP; VF0=0.5_WP
   ! call set_gas_eos_cof(Yv0)
   ! rhoL=get_RHOL(p0,T0)
   ! rhoG=get_RHOG(p0,T0,Yv0)
   ! Q0(1)=VF0*rhoL
   ! Q0(2)=(1.0_WP-VF0)*rhoG
   ! Q0(3)=Q0(1)*get_IL(rhoL,p0,Yv0)
   ! Q0(4)=Q0(2)*get_IG(rhoG,p0,Yv0)
   ! Q0(5:7)=0.0_WP
   ! Q0(8)=Q0(2)*Yv0
   ! VF=VF0; Q=Q0


   ! Test (this is the inputs for one of te cells that at the end of the time step, has a big value for Yv (close to 1))
   ! VF=0.90572023519589784_WP
   ! Q(1)=999.75432068566556_WP
   ! Q(2)=3.4482498824993049E-004_WP
   ! Q(3)=409911095.80795026_WP
   ! Q(4)=165.59270259183819_WP
   ! Q(5)=2265.3228035044976_WP
   ! Q(6)=-2295.9408590195899_WP
   ! Q(7)=-5.6361630364757843E-014_WP
   ! Q(8)=3.4482498824993050E-005_WP
   ! VF0=VF
   ! Q0=Q
   ! T0=369.53139902482570_WP
   ! Yv0=0.10000000000000001_WP

   ! print '(A)',       '========== INITIAL STATE =========='
   ! print '(A,ES15.7)','VF    = ',VF0
   ! print '(A,ES15.7)','rhoL  = ',Q0(1)/VF0
   ! print '(A,ES15.7)','rhoG  = ',Q0(2)/(1.0_WP-VF0)
   ! print '(A,ES15.7)','p_L   = ',get_PL(Q0(1)/VF0,Q0(3)/Q0(1))
   ! print '(A,ES15.7)','p_G   = ',get_PG(Q0(2)/(1.0_WP-VF0),Q0(4)/Q0(2),Yv0)
   ! print '(A,ES15.7)','T     = ',T0
   ! print '(A,ES15.7)','Yv    = ',Yv0
   ! print '(A,ES15.7)','Eps0 = ',sum(Q0(3:4))
   ! print '(A)',       '==================================='
   ! call PTg_relax(VF,Q)
   ! print '(A)',       '==================================='
   ! print '(A)',       '=========== FINAL STATE ==========='
   ! print '(A,ES15.7)','VF    = ',VF
   ! if (Q(2).gt.0.0_WP) print '(A,ES15.7)','Yv    = ',Q(8)/Q(2)
   ! print '(A)',       '==================================='
   ! print '(A)',       '========== CONSERVATION ==========='
   ! print '(A,ES15.7)','d(rho)/rho0  = ',(sum(Q(1:2))-sum(Q0(1:2)))/sum(Q0(1:2))
   ! print '(A,ES15.7)','d(Eps)/Eps0 = ',(sum(Q(3:4))-sum(Q0(3:4)))/sum(Q0(3:4))
   ! print '(A)',       '==================================='

   ! Write PTg_relax p-T data for comparison with reference saturation data.
   ! For the pure-water case, p_v = p_total because Yv0 = 1.
   ! For the air-water case, compare p_v = x_v*p_total against reference p_sat(T).
   call write_PTg_pT_curve('ptg_relax_pT_pure_water.csv', 1.0e5_WP, 0.5_WP, 1.0_WP, 300.0_WP, 500.0_WP, 10.0_WP)
   call write_PTg_pT_curve('ptg_relax_pT_Yv0_0p5.csv', 1.0e5_WP, 0.5_WP, 0.5_WP, 300.0_WP, 500.0_WP, 10.0_WP)

contains

   subroutine initialize_PTg_case(p0,T0,VF0,Yv0,VF,Q)
      implicit none
      real(WP), intent(in)  :: p0,T0,VF0,Yv0
      real(WP), intent(out) :: VF
      real(WP), intent(out) :: Q(8)
      real(WP) :: rhoL0,rhoG0

      call set_gas_eos_cof(Yv0)
      rhoL0=get_RHOL(p0,T0)
      rhoG0=get_RHOG(p0,T0,Yv0)

      Q(1)=VF0*rhoL0
      Q(2)=(1.0_WP-VF0)*rhoG0
      Q(3)=Q(1)*get_IL(rhoL0,p0,Yv0)
      Q(4)=Q(2)*get_IG(rhoG0,p0,Yv0)
      Q(8)=Q(2)*Yv0
      VF=VF0
   end subroutine initialize_PTg_case

   subroutine get_relaxed_PT_state(VF,Q,p_final,T_final,pv_final,Yv_final,xv_final,rhoL_final,rhoG_final,valid_state)
      implicit none
      real(WP), intent(in)  :: VF
      real(WP), intent(in)  :: Q(8)
      real(WP), intent(out) :: p_final,T_final,pv_final,Yv_final,xv_final,rhoL_final,rhoG_final
      integer,  intent(out) :: valid_state
      real(WP) :: pL,pG,TL,TG

      p_final=-1.0_WP
      T_final=-1.0_WP
      pv_final=-1.0_WP
      Yv_final=0.0_WP
      xv_final=0.0_WP
      rhoL_final=0.0_WP
      rhoG_final=0.0_WP
      valid_state=0

      if (Q(2).gt.0.0_WP) then
         Yv_final=max(0.0_WP,min(1.0_WP,Q(8)/Q(2)))
      else
         Yv_final=0.0_WP
      end if
      call set_gas_eos_cof(Yv_final)
      xv_final=get_xv(Yv_final)

      if ((VF.gt.1.0e-12_WP).and.(Q(1).gt.0.0_WP)) then
         rhoL_final=Q(1)/VF
         pL=get_PL(rhoL_final,Q(3)/Q(1),Yv_final)
         TL=get_TL(rhoL_final,pL,Yv_final)
         p_final=pL
         T_final=TL
      else if (((1.0_WP-VF).gt.1.0e-12_WP).and.(Q(2).gt.0.0_WP)) then
         rhoG_final=Q(2)/(1.0_WP-VF)
         pG=get_PG(rhoG_final,Q(4)/Q(2),Yv_final)
         TG=get_TG(rhoG_final,pG,Yv_final)
         if (p_final.le.0.0_WP) then
            p_final=pG
            T_final=TG
         end if
      end if

      if ((p_final.gt.0.0_WP).and.(T_final.gt.0.0_WP)) then
         pv_final=xv_final*p_final
         valid_state=1
      end if
   end subroutine get_relaxed_PT_state

   real(WP) function PTsat_residual_for_output(p_l,p_v,T)
      implicit none
      real(WP), intent(in) :: p_l,p_v,T
      if (((p_v+PinfV).gt.1.0e-10_WP).and.((p_l+PinfL).gt.1.0e-10_WP).and.(T.gt.0.0_WP)) then
         PTsat_residual_for_output=AS+BS/T+CS*log(T)+DS*log(p_l+PinfL)-log(p_v+PinfV)
      else
         PTsat_residual_for_output=huge(1.0_WP)
      end if
   end function PTsat_residual_for_output

   subroutine write_PTg_pT_curve(file_name,p0,VF0,Yv0,Tmin,Tmax,dT)
      implicit none
      character(len=*), intent(in) :: file_name
      real(WP),         intent(in) :: p0,VF0,Yv0,Tmin,Tmax,dT
      real(WP) :: VF,Q(8),Qinit(8)
      real(WP) :: T0,p_final,T_final,pv_final,Yv_final,xv_final,rhoL_final,rhoG_final
      real(WP) :: rho_err,Eps_err,rhoA0,rhoA_err,PTsat_res
      integer  :: unit,case_id,valid_state

      open(newunit=unit,file=trim(file_name),status='replace',action='write')
      write(unit,'(A)') 'case_id,p0,T0,VF0,Yv0,p_total_final,T_final,pv_final,Yv_final,xv_final,VF_final,rhoL_final,rhoG_final,rho_err,Eps_err,rhoA_err,PTsat_res,valid_state'

      case_id=0
      T0=Tmin
      do while (T0.le.Tmax+0.5_WP*dT)
         case_id=case_id+1
         call initialize_PTg_case(p0,T0,VF0,Yv0,VF,Q)
         Qinit=Q
         rhoA0=(1.0_WP-Yv0)*Qinit(2)

         call PTg_relax(VF,Q)
         call get_relaxed_PT_state(VF,Q,p_final,T_final,pv_final,Yv_final,xv_final,rhoL_final,rhoG_final,valid_state)

         rho_err=(sum(Q(1:2))-sum(Qinit(1:2)))/max(sum(Qinit(1:2)),tiny(1.0_WP))
         Eps_err=(sum(Q(3:4))-sum(Qinit(3:4)))/max(sum(Qinit(3:4)),tiny(1.0_WP))
         if (Q(2).gt.0.0_WP) then
            rhoA_err=((1.0_WP-Yv_final)*Q(2)-rhoA0)/max(sum(Qinit(1:2)),tiny(1.0_WP))
         else
            rhoA_err=-rhoA0/max(sum(Qinit(1:2)),tiny(1.0_WP))
         end if
         PTsat_res=PTsat_residual_for_output(p_final,pv_final,T_final)

         write(unit,'(*(G0.16,:,","))') case_id,p0,T0,VF0,Yv0,p_final,T_final,pv_final,Yv_final,xv_final,VF, &
              rhoL_final,rhoG_final,rho_err,Eps_err,rhoA_err,PTsat_res,valid_state

         T0=T0+dT
      end do

      close(unit)
      print '(A,A)', 'Wrote PTg_relax p-T curve to ', trim(file_name)
   end subroutine write_PTg_pT_curve

end program test_relax
