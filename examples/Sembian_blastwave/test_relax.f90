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

   !> Thermo-chemical relaxation model
   subroutine PTg_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP), dimension(:), allocatable    :: Q0
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG
      real(WP) :: VFeq,VF0,Peq,p,ppv,Yv,xv,T,Teq,Tsat,pOld,YvOld
      real(WP) :: RHOL,RHOG
      real(WP) :: rho0,Eps0,rhoA0
      real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
      real(WP) :: p_pert,Yv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
      real(WP), parameter :: fd_eps=1.0e-8_WP,lnP_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=1.0e-5_WP,Yvmax=1.0_WP-Yvmin
      real(WP), parameter :: p_tol=1.0e-4_WP,Yv_tol=1.0e-4_WP,Tsat_tol=1.0e-5_WP,rho_tol=1.0e-4_WP,Eps_tol=1.0e-4_WP,F1_tol=1e-4_WP,F2_tol=1e-4_WP
      integer  :: it,itmax
      real(WP) :: p_err,Yv_err,Tsat_err,rho_err,Eps_err
      logical  :: converge
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
      ! Get the thermo-mechanically relaxed temperature
      Teq=get_TL(RHO=Q(1)/VF,P=p)
      ! Get vapor mole fraction and partial pressure
      xv=get_xv(Yv)
      ppv=xv*p
      if (.not.check_p(ppv)) then
         print*,"****************** Vapor partial pressure too low. Skipping the cell!"
         deallocate(Q0)
         return
      end if
      ! Get saturation temperature from vapor partial pressure
      itmax=20
      converge=.false.
      Tsat=Teq
      do it=1,itmax
         T=Tsat
         Tsat=T-PTsat(ppv,T)/dPTsatdT(T)
         Tsat_err=abs((Tsat-T)/T)
         if (Tsat_err.lt.Tsat_tol) then
            converge=.true.
            exit
         end if
      end do
      if (.not.converge) then
         print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
         deallocate(Q0)
         return
      end if
      print '(A,I2)','Tsat it= ',it
      print '(A,ES15.7)','T from PT_relax    = ',Teq
      print '(A,ES15.7)','Tsat    = ',Tsat
      ! Activate chemical relaxation only for metastable states
      if (Teq.le.Tsat) then
         deallocate(Q0)
         return
      end if
      ! Get the initial quantities
      rho0=sum(Q(1:2))
      Eps0=sum(Q(3:4))
      rhoA0=(1.0_WP-Yv)*Q(2)
      ! Iteratively solve for the equilibrium pressure and vapor mass fraction
      itmax=20
      converge=.false.
      p_err=10.0_WP*p_tol
      Yv_err=10.0_WP*Yv_tol
      do it=1,itmax
         print '(A,I2)','it = ',it
         print '(A,ES15.7)','p    = ',p
         print '(A,ES15.7)','Yv    = ',Yv
         ! Evaluate residuals at current state
         call get_T(p,Yv)
         print '(A,ES15.7)','T    = ',T
         xv=get_xv(Yv)
         ppv=xv*p
         print '(A,ES15.7)','xv    = ',xv
         print '(A,ES15.7)','ppv    = ',ppv
         if (.not.check_p(ppv)) then
            print*,"****************** Vapor partial pressure too low. Skipping the cell!"
            deallocate(Q0)
            return
         end if
         F1=PTsat(ppv,T)
         F2=Eps_res(p,Yv)
         print '(A,ES15.7)','F1 error    = ',abs(F1)
         print '(A,ES15.7)','F2 error    = ',abs(F2/Eps0)
         if ((abs(F1).lt.F1_tol).and.(abs(F2/Eps0).lt.F2_tol)) then
            converge=.true.
            exit
         end if
         ! Compute Jacobian via finite difference
         ! --- Perturbation in p ---
         p_pert=p*(1.0_WP+fd_eps)
         ppv=xv*p_pert
         if (.not.check_p(ppv)) then
            print*,"****************** Vapor partial pressure too low. Skipping the cell!"
            deallocate(Q0)
            return
         end if
         call get_T(p_pert,Yv)
         F1p=PTsat(ppv,T)
         F2p=Eps_res(p_pert,Yv)
         dF1dp=(F1p-F1)/(p_pert-p)
         dF2dp=(F2p-F2)/(p_pert-p)
         ! --- Perturbation in Yv ---
         Yv_pert=Yv+fd_eps
         if (Yv_pert.gt.Yvmax) then
            Yv_pert=Yv-fd_eps
         end if
         if (Yv_pert.lt.Yvmin) then
            Yv_pert=Yv+fd_eps
         end if
         xv=get_xv(Yv_pert)
         ppv=xv*p
         if (.not.check_p(ppv)) then
            print*,"****************** Vapor partial pressure too low. Skipping the cell!"
            deallocate(Q0)
            return
         end if
         call get_T(p,Yv_pert)
         F1Y=PTsat(ppv,T)
         F2Y=Eps_res(p,Yv_pert)
         dF1dYv=(F1Y-F1)/(Yv_pert-Yv)
         dF2dYv=(F2Y-F2)/(Yv_pert-Yv)
         ! Solve 2x2 system: J * [dp; dYv] = -[F1; F2]
         detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
         if (abs(detJ).lt.1.0e-30_WP) exit
         dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
         dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
         ! Newton-Raphson update
         pOld=p
         YvOld=Yv
         p=p+dp_nr
         print '(A,ES15.7)','Yv new unclipped = ',Yv+dYv_nr
         Yv=max(Yvmin,min(Yvmax,Yv+dYv_nr))
         print '(A,ES15.7)','p  new = ',p
         print '(A,ES15.7)','Yv new = ',Yv
         ! Evaluate p and T errors
         p_err=abs((p-pOld)/pOld)
         Yv_err=abs((Yv-YvOld)/(YvOld+1.0e-30_WP))
         print '(A,ES15.7)','p  error = ',p_err
         print '(A,ES15.7)','Yv error = ',Yv_err
         if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol)) then
            converge=.true.
            exit
         end if
      end do
      print '(A,I2)','p-Yv iterations = ',it
      print '(A,ES15.7)','p    = ',p
      print '(A,ES15.7)','Yv    = ',Yv
      if (.not.converge) then
         print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
         deallocate(Q0)
         return
      end if
      ! Get equilibrium quantities at converged (p, Yv)
      call get_T(p,Yv)
      RHOL=get_RHOL(p,T)
      call set_gas_eos_cof(Yv)
      RHOG=get_RHOG(p,T,Yv)
      VFeq=(rho0-RHOG)/(RHOL-RHOG)
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; p=max(p,-PinfL); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; p=max(p,-PinfG); end if
      ! Adjust conserved quantities
      Q(1)=(       VFeq)*RHOL
      Q(2)=(1.0_WP-VFeq)*RHOG
      Q(3)=Q(1)*get_IL_PT(p,T)
      Q(4)=Q(2)*get_IG_PT(p,T,Yv)
      Q(8)=Q(2)*Yv
      VF=VFeq
      ! Evaluate conservation
      rho_err=(sum(Q(1:2))-rho0)/rho0
      Eps_err=(sum(Q(3:4))-Eps0)/Eps0
      print '(A,ES15.7)','rho_err    = ',rho_err
      print '(A,ES15.7)','rho_tol    = ',rho_tol
      print '(A,ES15.7)','Eps_err    = ',Eps_err
      print '(A,ES15.7)','Eps_tol    = ',Eps_tol
      if ((abs(rho_err).gt.rho_tol).or.(abs(Eps_err).gt.Eps_tol)) then
         ! Reset the output to the initial values fed into chemical relaxation algorithm
         Q=Q0
         VF=VF0
         deallocate(Q0)
         print*,"****************** Conservation is violated. Skipping the cell!!"
         return
      end if
      if (VF.lt.VFmin) then
         ! Reset the output to the initial values fed into chemical relaxation algorithm
         Q=Q0
         VF=VF0
         deallocate(Q0)
         print*,"****************** Not enough liquid to vaporize. Skipping the cell!!"
         return
      end if
      if (allocated(Q0)) deallocate(Q0)
      contains
      ! Equilibrium temperature as a function of pressure and vapor mass fraction
      subroutine get_T(p_in,Yv_in)
         real(WP), intent(in) :: p_in,Yv_in
         T=(1.0_WP-Yv_in)/((rho0*(1.0_WP-Yv_in)-rhoA0)*(GammaL-1.0_WP)*CvL/(p_in+PinfL)+rhoA0*((GammaV-1.0_WP)*CvV*Yv_in/(p_in+PinfV)+(GammaA-1.0_WP)*CvA*(1.0_WP-Yv_in)/(p_in+PinfA)))
      end subroutine get_T
      ! Function that defines p-T saturation curve
      function PTsat(p_in,T_in)
         real(WP), intent(in) :: p_in,T_in
         real(WP) :: PTsat
         PTsat=AS+BS/T_in+CS*log(T_in)+DS*log(p_in+PinfL)-log(p_in+PinfV)
      end function PTsat
      ! Temperature derivative of p-T saturation curve function
      function dPTsatdT(T_in)
         real(WP), intent(in) :: T_in
         real(WP) :: dPTsatdT
         dPTsatdT=-BS/T_in**2+CS/T_in
      end function dPTsatdT
      ! Residual of the internal energy conservation equation
      function Eps_res(p_in,Yv_in)
         real(WP), intent(in) :: p_in,Yv_in
         real(WP) :: Eps_res
         Eps_res=(rho0*(1.0_WP-Yv_in)-rhoA0)*(CvL*T*(p_in+GammaL*PinfL)/(p_in+PinfL)+etaL)+rhoA0*((CvV*T*(p_in+GammaV*PinfV)/(p_in+PinfV)+etaV)*Yv_in+(CvA*T*(p_in+GammaA*PinfA)/(p_in+PinfA)+etaA)*(1.0_WP-Yv_in))-Eps0*(1.0_WP-Yv_in)
      end function Eps_res
      ! Sanity check pressure value
      logical function check_p(p_in)
         real(WP), intent(in) :: p_in
         check_p=((p_in+PinfV).gt.lnP_eps).and.((p_in+PinfL).gt.lnP_eps)
      end function check_p
   end subroutine PTg_relax

end module relax_mod

program test_relax
   use eos_mod
   use relax_mod
   implicit none
   real(WP) :: VF,Q(8),Q0(8),VF0
   real(WP) :: T0,p0,rhoL,rhoG,Yv0,VFmin,VFmax,Yvmin,Yvmax

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

   ! Outside the following range either blows up or violates mass/energy conservation
   VFmin=1e-5_WP
   VFmax=1.0_WP-VFmin
   Yvmin=1e-5_WP
   Yvmax=1.0_WP-Yvmin

   ! --- Setup metastable liquid --- 
   ! (Yv0=Yvmin; VF0=VFmin) fails (conservation issue: This I think is where we need the cluster method).
   ! Also, (Yv0=Yvmax; VF0=VFmin) fails.
   ! (Yv0=Yvmax; VF0=VFmax) converges with some energy error.
   p0=1.0e5_WP; T0=350.0_WP; Yv0=0.5_WP; VF0=0.0000001_WP
   call set_gas_eos_cof(Yv0)
   rhoL=(p0+PinfL)/((GammaL-1.0_WP)*CvL*T0)
   rhoG=(p0+PinfG)/((GammaG-1.0_WP)*CvG*T0)
   Q0(1)=VF0*rhoL
   Q0(2)=(1.0_WP-VF0)*rhoG
   Q0(3)=Q0(1)*get_IL(rhoL,p0,Yv0)
   Q0(4)=Q0(2)*get_IG(rhoG,p0,Yv0)
   Q0(5:7)=0.0_WP
   Q0(8)=Q0(2)*Yv0
   VF=VF0; Q=Q0

   print '(A)',       '========== INITIAL STATE =========='
   print '(A,ES15.7)','VF    = ',VF0
   print '(A,ES15.7)','rhoL  = ',Q0(1)/VF0
   print '(A,ES15.7)','rhoG  = ',Q0(2)/(1.0_WP-VF0)
   print '(A,ES15.7)','p_L   = ',get_PL(Q0(1)/VF0,Q0(3)/Q0(1))
   print '(A,ES15.7)','p_G   = ',get_PG(Q0(2)/(1.0_WP-VF0),Q0(4)/Q0(2),Yv0)
   print '(A,ES15.7)','T     = ',T0
   print '(A,ES15.7)','Yv    = ',Yv0
   print '(A,ES15.7)','Eps0 = ',sum(Q0(3:4))

   call PTg_relax(VF,Q)

   print '(A)',       '========== FINAL STATE =========='
   print '(A,ES15.7)','VF    = ',VF
   if (Q(2).gt.0.0_WP) print '(A,ES15.7)','Yv    = ',Q(8)/Q(2)
   print '(A)',       '========== CONSERVATION =========='
   print '(A,ES15.7)','d(rho)/rho0  = ',(sum(Q(1:2))-sum(Q0(1:2)))/sum(Q0(1:2))
   print '(A,ES15.7)','d(Eps)/Eps0 = ',(sum(Q(3:4))-sum(Q0(3:4)))/sum(Q0(3:4))

contains

   subroutine run_captured_case(label,VFin,Qin)
      character(len=*), intent(in) :: label
      real(WP), intent(in) :: VFin
      real(WP), dimension(8), intent(in) :: Qin
      real(WP) :: VFcase,Qcase(8),Qcase0(8),Yvin

      VFcase=VFin
      Qcase=Qin
      Qcase0=Qcase
      if (Qcase(2).gt.0.0_WP) then
         Yvin=Qcase(8)/Qcase(2)
      else
         Yvin=0.0_WP
      end if

      print '(A)',       '========== CAPTURED CASE =========='
      print '(A)',       trim(label)
      print '(A,ES15.7)','VF in = ',VFcase
      print '(A,ES15.7)','Yv in = ',Yvin
      print '(A,8ES15.7)','Q in  = ',Qcase
      if (VFcase.gt.0.0_WP.and.Qcase(1).gt.0.0_WP) then
         print '(A,ES15.7)','rhoL in = ',Qcase(1)/VFcase
         print '(A,ES15.7)','PL in   = ',get_PL(Qcase(1)/VFcase,Qcase(3)/Qcase(1))
         print '(A,ES15.7)','TL in   = ',get_TL(Qcase(1)/VFcase,get_PL(Qcase(1)/VFcase,Qcase(3)/Qcase(1)))
      end if
      if ((1.0_WP-VFcase).gt.0.0_WP.and.Qcase(2).gt.0.0_WP) then
         print '(A,ES15.7)','rhoG in = ',Qcase(2)/(1.0_WP-VFcase)
         print '(A,ES15.7)','PG in   = ',get_PG(Qcase(2)/(1.0_WP-VFcase),Qcase(4)/Qcase(2),Yvin)
         print '(A,ES15.7)','TG in   = ',get_TG(Qcase(2)/(1.0_WP-VFcase), &
         & get_PG(Qcase(2)/(1.0_WP-VFcase),Qcase(4)/Qcase(2),Yvin),Yvin)
      end if

      call PTg_relax(VFcase,Qcase)

      print '(A,ES15.7)','VF out = ',VFcase
      if (Qcase(2).gt.0.0_WP) print '(A,ES15.7)','Yv out = ',Qcase(8)/Qcase(2)
      print '(A,8ES15.7)','Q out = ',Qcase
      print '(A,ES15.7)','d(rho)/rho0 = ',(sum(Qcase(1:2))-sum(Qcase0(1:2)))/sum(Qcase0(1:2))
      print '(A,ES15.7)','d(Eps)/Eps0 = ',(sum(Qcase(3:4))-sum(Qcase0(3:4)))/sum(Qcase0(3:4))
   end subroutine run_captured_case
end program test_relax
