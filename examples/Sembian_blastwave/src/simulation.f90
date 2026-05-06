!> Sembian blastwave case – high-pressure region initialization (dimensional)
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use amrio_class,       only: amrio
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> Stiffened gas EOS parameters (liquid and gas)
   real(WP) :: GammaL,PinfL,etaL,etapL,CvL,CpL
   real(WP) :: GammaG,PinfG,etaG,etapG,CvG,CpG

   !> Gas EOS parameters (air + vapor)
   real(WP) :: GammaA,PinfA,etaA,etapA,CvA,CpA
   real(WP) :: GammaV,PinfV,etaV,etapV,CvV,CpV
   real(WP) :: AS,BS,CS,DS

   !> Molar mass of the gas
   real(WP), parameter :: Mv=0.0180153_WP,Ma=0.02897_WP

   !> Flow parameters
   real(WP) :: Grho0,GP0          !< Pre-shock gas state
   real(WP) :: Grho1,GP1,u1       !< Post-shock gas state (from Rankine-Hugoniot)
   real(WP) :: M1                 !< Post-shock Mach number
   real(WP) :: relshockvel        !< Velocity at which shock moves
   real(WP) :: Lrho0              !< Liquid density
   real(WP) :: Xs                 !< Shock location [m]
   real(WP) :: Ms                 !< Shock Mach number
   real(WP) :: muG,muL            !< Dynamic viscosities

   !> High-pressure region parameters
   real(WP) :: HP_thickness       !< Thickness of HP region [m]
   real(WP) :: HP_center          !< Center of HP region [m]
   real(WP) :: HP_density         !< Density in HP region [kg/m^3]
   real(WP) :: HP_pressure        !< Pressure in HP region [Pa]

   !> Cylinder parameters
   real(WP) :: dcyl               !< Cylinder diameter [m]
   real(WP) :: xcyl               !< Cylinder center x location [m]

   !> Domain dimensions
   real(WP) :: Lx,Ly              !< Domain lengths [m]

   !> Tagging parameters
   real(WP) :: vorticity_tag=huge(1.0_WP)
   real(WP) :: rho_ratio_tag=huge(1.0_WP)

contains

   !> Levelset function for 2D cylinder centered at (xcyl, 0)
   function levelset_cyl(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function levelset_cyl

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
   !> Liquid EOS: I=f(P,T)
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
   !> Gas EOS: RHOG=f(P,T,Yv)
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

   !> Mechanical relaxation model
   subroutine P_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq,Yv
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! Set gas EoS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      call set_gas_eos_cof(Yv)
      ! ================ Handle gas flotsams ================
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! ================ First step for mechanical relaxation ================
      ! Get phasic pressures
      PL=get_PL(RHO=Q(1)/(       VF),I=Q(3)/Q(1))
      PG=get_PG(RHO=Q(2)/(1.0_WP-VF),I=Q(4)/Q(2),Yv=Yv)
      ! Handle limit cases - should mass/energy be tranasfered or lost? - this should probably never happen...
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
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=VF*((gammaL-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+gammaL)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax

   !> Implicit mechanical relaxation for stiffened gas EOS pair
   subroutine P_relax_implicit(VF,Q)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP) :: a,b,d,Peq,VFeq,Yv
      real(WP) :: invG1G,invG1L,d0,d1,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      ! Set gas EoS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      call set_gas_eos_cof(Yv)
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Precompute EOS constants
      invG1G=1.0_WP/(GammaG-1.0_WP)
      invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=GammaL*PinfL*invG1L
      d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G
      facL=invG1G+VF
      ! Quadratic coefficients: a*Peq^2 + b*Peq + d = 0
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Solve for equilibrium pressure (positive root)
      if (b**2-4.0_WP*a*d.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Bail if pressure is unphysical
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Equilibrium volume fraction from liquid energy constraint
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Update internal energies via p*dV work exchange
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit

   !> Thermo-mechanical relaxation model
   subroutine PT_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq,Teq,Yv
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
      Teq=get_TL(RHO=Q(1)/VFeq,P=Peq)
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*get_RHOLIL(Peq,Teq)
      Q(4)=(1.0_WP-VFeq)*get_RHOGIG(Peq,Teq,Yv)
      VF=VFeq
   end subroutine PT_relax

   !> Thermo-chemical relaxation model
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
      real(WP), parameter :: fd_eps=1.0e-8_WP,lnP_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,ppv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: p_tol=1.0e-4_WP,Yv_tol=1.0e-4_WP,Tsat_tol=1.0e-5_WP,rho_tol=1.0e-4_WP,Eps_tol=1.0e-4_WP,F1_tol=1e-12_WP,F2_tol=1e-12_WP
      real(WP), parameter :: p_tol_lv=1.0e-7_WP
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
         ! print*,'****************** Trace inert gas mass fraction. Using pure LV chemical relaxation!'
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
         print*,"****************** Conservation is violated. Skipping the cell!!"
         return
      end if
      ! Check VOF
      if (VF.lt.VFmin) then
         call restore_VFQ()
         call deallocate_Q0()
         print*,"****************** Not enough liquid to vaporize. Skipping the cell!!"
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
                     print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                     return
                  end if
                  if (ppv.ge.p_eq) then
                     ! ppv_sat >= p: flash regime. Seed with a small Yv so that
                     ! get_T stays well-conditioned (energy balance requires Yv << 1).
                     Yv_eq=sqrt(fd_eps)
                  else
                     xv=ppv/p_eq
                     Yv_eq=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
                  end if
                  ! Yv_eq=max(Yvmin+fd_eps,min(Yvmax-fd_eps,Yv_eq))
                  Yv_eq=max(Yvmin,min(Yvmax,Yv_eq))
                  print '(A,ES15.7)','Seeded Yv from saturation at Teq = ',Yv_eq
               else
                  print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
            else
               ! Get saturation temperature from
               call get_Tsat(ppv,T_eq,Tsat,converge,Tsat_it)
               if (.not.converge) then
                  print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
                  return
               end if
               ! Activate chemical relaxation only for metastable states
               if (T_eq.le.Tsat) return
            end if
            activate_chem_relax=.true.
         end function activate_chem_relax
         !> Safeguarded Newton solve for saturation temperature at fixed vapor partial pressure
         subroutine get_Tsat(ppv,Tguess,Tsat,converged,Tsat_it)
            real(WP), intent(in)  :: ppv,Tguess
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
            Flo=PTsat(ppv,Tlo)
            Fhi=PTsat(ppv,Thi)
            ! Expand the bracket if needed.  In the normal liquid-vapor range
            ! PTsat is monotone in T, so these one-sided expansions are enough.
            expand_it=0
            do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
               if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
                  Tlo=max(1.0_WP,0.8_WP*Tlo)
                  Flo=PTsat(ppv,Tlo)
               else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
                  Thi=1.2_WP*Thi
                  Fhi=PTsat(ppv,Thi)
               else
                  exit
               end if
               expand_it=expand_it+1
            end do
            if (Flo*Fhi.gt.0.0_WP) then
               print*,'****************** Could not bracket Tsat!',ppv,Flo,Fhi,Tlo,Thi
               return
            end if
            ! Use the caller's guess only after clamping it to the safe bracket.
            Tsat=max(Tlo,min(Thi,Tguess))
            do it=1,Tsat_itmax
               Told=Tsat
               Fold=PTsat(ppv,Told)
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
               Fnew=PTsat(ppv,Tnew)
               if (Fnew.ne.Fnew) then
                  Tnew=0.5_WP*(Tlo+Thi)
                  Fnew=PTsat(ppv,Tnew)
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
               pOld=p_eq
               p_eq=pOld-PTsat(pOld,T_eq)/dPTsatdp(pOld,T_eq,dTdp)
               ! Evaluate the error
               p_err=abs((p_eq-pOld)/pOld)
               if (p_err.lt.p_tol) then
                  converged=.true.
                  exit
               end if
            end do
            if (.not.converged) then
               print*,"****************** p iterations blew up. Skipping the cell!!"
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
                  print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
               F1=PTsat(ppv,T_eq)
               F2=Eps_res_lvg(p_eq,Yv_eq,T_eq)
               res0=sqrt(F1**2+(F2/Eps0)**2)
               if ((abs(F1).lt.F1_tol).and.(abs(F2/Eps0).lt.F2_tol)) then
                  converged=.true.
                  exit
               end if
               ! Compute Jacobian via finite difference
               ! --- Perturbation in p ---
               p_pert=p_eq*(1.0_WP+fd_eps)
               ppv_pert=xv*p_pert
               if (.not.check_p(ppv_pert)) then
                  print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
               T_pert=get_T_lvg(p_pert,Yv_eq)
               F1p=PTsat(ppv_pert,T_pert)
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
                  print*,"****************** Vapor partial pressure too low. Skipping the cell!"
                  return
               end if
               T_pert=get_T_lvg(p_eq,Yv_pert)
               F1Y=PTsat(ppv_pert,T_pert)
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
                  F1_try=PTsat(ppv_try,T_try)
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
               if (.not.accepted) exit
               ! Evaluate p and Yv errors
               p_err=abs((p_eq-pOld)/pOld)
               Yv_err=abs((Yv_eq-YvOld)/(YvOld+1.0e-30_WP))
               if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol).and.(abs(F1_try).lt.F1_tol).and.(abs(F2_try/Eps0).lt.F2_tol)) then
                  converged=.true.
                  exit
               end if
            end do
            if (.not.converged) then
               print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
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
         function PTsat(pv,T_eq)
            real(WP), intent(in) :: pv,T_eq
            real(WP) :: PTsat
            PTsat=AS+BS/T_eq+CS*log(T_eq)+DS*log(pv+PinfL)-log(pv+PinfV)
         end function PTsat
         !> Temperature derivative of p-T saturation curve function
         function dPTsatdT(T_eq)
            real(WP), intent(in) :: T_eq
            real(WP) :: dPTsatdT
            dPTsatdT=-BS/T_eq**2+CS/T_eq
         end function dPTsatdT
         !> Pressure derivative of p-T saturation curve function for pure vapor p iteration
         function dPTsatdp(pv,T_eq,dTdp)
            real(WP), intent(in) :: pv,T_eq,dTdp
            real(WP) :: dPTsatdp
            dPTsatdp=dPTsatdT(T_eq)*dTdp+DS/(pv+PinfL)-1.0_WP/(pv+PinfV)
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
         logical function check_p(pv)
            real(WP), intent(in) :: pv
            check_p=((pv+PinfV).gt.lnP_eps).and.((pv+PinfL).gt.lnP_eps)
         end function check_p
   end subroutine PTg_relax

   !> Compute viscosity: constant gas and liquid, VF-weighted blend
   !> Contains commented-out Sutherland law for variable gas viscosity (dimensional form)
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pVisc,pBeta,pCond,pDiff
      real(WP) :: mu_g,mu_l
      real(WP), parameter :: myeps=1.0e-15_WP
      !> Sutherland's law parameters (dimensional, SI units)
      !> mu_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K
      ! real(WP), parameter :: mu_ref=1.716e-5_WP  !< Reference viscosity [Pa·s]
      ! real(WP), parameter :: T_ref=273.15_WP      !< Reference temperature [K]
      ! real(WP), parameter :: S_suth=110.4_WP      !< Sutherland constant [K]
      ! real(WP) :: T_gas                           !< Local gas temperature [K]
      ! Loop over levels
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pCond=>fs%cond%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! --- Gas viscosity ---
               ! Constant gas viscosity (current default)
               mu_g=muG
               ! Sutherland's law for variable gas viscosity (dimensional):
               ! T_gas=pTG(i,j,k,1)  ! Gas temperature from EOS [K]
               ! mu_g=mu_ref*(T_gas/T_ref)**1.5_WP*(T_ref+S_suth)/(T_gas+S_suth)
               ! --- Liquid viscosity ---
               mu_l=muL
               ! Mixture viscosity (harmonic averaging)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps))
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Zero thermal conductivity
               pCond(i,j,k,1)=0.0_WP
               ! Zero mass diffusivity
               pDiff(i,j,k,1)=0.0_WP
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback – set Q and VF/barycenters for blastwave + cylinder
   subroutine blastwave_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      use param, only: param_read
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG_local,pG_local,IG_local,Yv0
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energy of liquid at ambient
      IEL=get_IL(Lrho0,GP0)
      ! Get initial vapor mass fraction
      call param_read('Initial vapor mass fraction',Yv0,default=1.0e-1_WP)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset (cylinder)
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=levelset_cyl,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Cell-center x coordinate
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            ! Determine gas state: high pressure region or ambient
            if ((x_cc.gt.(HP_center-0.5_WP*HP_thickness)).and.(x_cc.lt.(HP_center+0.5_WP*HP_thickness))) then
               ! High pressure region
               rhoG_local=HP_density
               pG_local  =HP_pressure
            else
               ! Ambient region
               rhoG_local=Grho0
               pG_local  =GP0
            end if
            IG_local=get_IG(rhoG_local,pG_local,Yv0)
            ! Set conserved variables
            pQ(i,j,k,1)=(       myVF)*Lrho0
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG_local
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*IG_local
            pQ(i,j,k,5)=0.0_WP
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
            pQ(i,j,k,8)=pQ(i,j,k,2)*Yv0
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine blastwave_init

   !> Tagger based on vorticity and density ratio (from amrcomp_drop)
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP) :: irho_cc,irho_xp,irho_xm,irho_yp,irho_ym,irho_zp,irho_zm
      real(WP) :: vort_x,vort_y,vort_z,vort_mag
      real(WP) :: rho_max,rho_min,rho_nb,rho_ratio
      integer :: i,j,k,ii,jj,kk
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                  ! Get local inverse densities
                  irho_cc=1.0_WP/max(sum(pQ(i  ,j,  k,  1:2)),solver%rho_floor)
                  irho_xp=1.0_WP/max(sum(pQ(i+1,j,  k,  1:2)),solver%rho_floor)
                  irho_xm=1.0_WP/max(sum(pQ(i-1,j,  k,  1:2)),solver%rho_floor)
                  irho_yp=1.0_WP/max(sum(pQ(i,  j+1,k,  1:2)),solver%rho_floor)
                  irho_ym=1.0_WP/max(sum(pQ(i,  j-1,k,  1:2)),solver%rho_floor)
                  irho_zp=1.0_WP/max(sum(pQ(i,  j,  k+1,1:2)),solver%rho_floor)
                  irho_zm=1.0_WP/max(sum(pQ(i,  j,  k-1,1:2)),solver%rho_floor)
                  ! Compute vorticity and tag based on it
                  vort_x=(pQ(i,j+1,k,7)*irho_yp-pQ(i,j-1,k,7)*irho_ym)*0.5_WP*dyi-(pQ(i,j,k+1,6)*irho_zp-pQ(i,j,k-1,6)*irho_zm)*0.5_WP*dzi
                  vort_y=(pQ(i,j,k+1,5)*irho_zp-pQ(i,j,k-1,5)*irho_zm)*0.5_WP*dzi-(pQ(i+1,j,k,7)*irho_xp-pQ(i-1,j,k,7)*irho_xm)*0.5_WP*dxi
                  vort_z=(pQ(i+1,j,k,6)*irho_xp-pQ(i-1,j,k,6)*irho_xm)*0.5_WP*dxi-(pQ(i,j+1,k,5)*irho_yp-pQ(i,j-1,k,5)*irho_ym)*0.5_WP*dyi
                  vort_mag=sqrt(vort_x**2+vort_y**2+vort_z**2)
                  if (vort_mag.gt.vorticity_tag) tagarr(i,j,k,1)=SETtag
                  ! Compute density ratio in 3x3x3 stencil and tag based on it
                  rho_max=solver%rho_floor; rho_min=huge(1.0_WP)
                  do kk=-1,1; do jj=-1,1; do ii=-1,1
                           rho_nb=sum(pQ(i+ii,j+jj,k+kk,1:2))
                           rho_max=max(rho_max,rho_nb)
                           rho_min=min(rho_min,max(rho_nb,solver%rho_floor))
                        end do; end do; end do
                  rho_ratio=rho_max/rho_min
                  if (rho_ratio.gt.rho_ratio_tag) tagarr(i,j,k,1)=SETtag
               end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         ! Get fluid parameters
         call param_read('Liquid gamma',GammaL)
         !call param_read('Gas gamma',GammaG)
         call param_read('Liquid Pinf',PinfL)
         !call param_read('Gas Pinf',PinfG)
         call param_read('Liquid eta',etaL)
         !call param_read('Gas eta',etaG)
         call param_read('Liquid etap',etapL)
         !call param_read('Gas etap',etapG)
         call param_read('Liquid specific heat',CvL)
         !call param_read('Gas specific heat',CvG)
         ! Read gas species parameters (air + vapor)
         call param_read('Air gamma',GammaA)
         call param_read('Vapor gamma',GammaV)
         call param_read('Air specific heat',CvA)
         call param_read('Vapor specific heat',CvV)
         call param_read('Air Pinf',PinfA)
         call param_read('Vapor Pinf',PinfV)
         call param_read('Air eta',etaA)
         call param_read('Vapor eta',etaV)
         call param_read('Air etap',etapA)
         call param_read('Vapor etap',etapV)
         ! Derive Cp from Gamma and Cv
         CpL=GammaL*CvL
         CpA=GammaA*CvA
         CpV=GammaV*CvV
         call param_read('Liquid density',Lrho0)
         call param_read('Pre-shock density',Grho0)
         call param_read('Pre-shock pressure',GP0)
         call param_read('Mach number of shock',Ms)
         call param_read('Shock location',Xs)
         call param_read('Liquid viscosity',muL)
         call param_read('Gas viscosity',muG)
         ! High-pressure region parameters
         call param_read('HP region thickness',HP_thickness)
         call param_read('HP center location',HP_center)
         call param_read('HP density',HP_density)
         call param_read('HP pressure',HP_pressure)
         ! Cylinder geometry
         call param_read('Cylinder diameter',dcyl)
         call param_read('Cylinder location',xcyl)
         ! Domain dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         ! Saturation curve coefficients
         AS=(CpL-CpV+etapV-etapL)/(CpV-CvV)
         BS=(etaL-etaV)/(CpV-CvV)
         CS=(CpV-CpL)/(CpV-CvV)
         DS=(CpL-CvL)/(CpV-CvV)
         ! Set gas EOS coefficients using air properties
         call set_gas_eos_cof(0.0_WP)
         ! Use shock relations to get post-shock numbers (for informational / tagging purposes)
         GP1 = GP0 * (2.0_WP*GammaG*Ms**2 - (GammaG-1.0_WP)) / (GammaG+1.0_WP)
         Grho1 = Grho0 * (Ms**2 * (GammaG+1.0_WP) / ((GammaG-1.0_WP)*Ms**2 + 2.0_WP))
         ! Calculate post-shock Mach number
         M1 = sqrt(((GammaG-1.0_WP)*(Ms**2)+2.0_WP)/(2.0_WP*GammaG*(Ms**2)-(GammaG-1.0_WP)))
         ! Calculate post-shock velocity
         u1 = -M1 * sqrt(GammaG*GP1/Grho1) + Ms*sqrt(GammaG*GP0/Grho0)
         ! Velocity at which shock moves
         relshockvel = -Grho1*u1/(Grho0-Grho1)
         ! Log setup
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  Grho0=",es12.5," GP0=",es12.5)') Grho0,GP0; call log(message)
         write(message,'("[Post-shock] Grho1=",es12.5," GP1=",es12.5," u1=",es12.5)') Grho1,GP1,u1; call log(message)
         write(message,'("[Shock velocity] =",es12.5)') relshockvel; call log(message)
         write(message,'("[HP region] center=",es12.5," thickness=",es12.5)') HP_center,HP_thickness; call log(message)
         write(message,'("[HP region] density=",es12.5," pressure=",es12.5)') HP_density,HP_pressure; call log(message)
         write(message,'("[Liquid] rhoL=",es12.5," GammaL=",es12.5," PinfL=",es12.5)') Lrho0,GammaL,PinfL; call log(message)
         write(message,'("[Gas]    GammaG=",es12.5," PinfG=",es12.5)') GammaG,PinfG; call log(message)
         write(message,'("[Visc]   muG=",es12.5," muL=",es12.5)') muG,muL; call log(message)
         write(message,'("[Cylinder] diameter=",es12.5," location=",es12.5)') dcyl,xcyl; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='Sembian_blastwave'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('AMReX blocking factor',amr%nbloc,default=8)
         amr%xlo=0.0_WP;      amr%xhi=Lx
         amr%ylo=-0.5_WP*Ly;  amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Ly;  amr%zhi=+0.5_WP*Ly
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zper=.true.
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap, amrex_bc_reflect_even, amrex_bc_reflect_odd
         use amrmpcomp_class, only: BC_REFLECT
         ! Create flow solver
         call fs%initialize(amr=amr,name='Sembian_blastwave')
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide pressure relaxation model
         fs%relax=>PTg_relax
         ! Set initial conditions via blastwave callback
         fs%user_mpcomp_init=>blastwave_init
         ! Set BCs
         if (.not.amr%xper) then
            ! x-low is a wall
            fs%lo_bc(1)=BC_REFLECT
            fs%Q%lo_bc(1,:)=amrex_bc_reflect_even
            fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd
            ! x-hi is extrapolation
            fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         end if
         if (.not.amr%yper) then
            fs%Q%lo_bc(2,:)=amrex_bc_foextrap
            fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         end if
         if (.not.amr%zper) then
            fs%Q%lo_bc(3,:)=amrex_bc_foextrap
            fs%Q%hi_bc(3,:)=amrex_bc_foextrap
         end if
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=8,ng=0,interp=amrex_interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=amrex_interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_mpcomp_tagging=>my_tagger
         call param_read('Tagging vorticity',vorticity_tag)
         call param_read('Tagging rho ratio',rho_ratio_tag)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
            ! Force a regrid immediately after restoring checkpoint
            call amr%regrid(baselvl=0,time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC and reset moments
            call fs%build_plic(time%t)
            call fs%reset_moments()
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'Sembian_blastwave',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_scalar(fs%Yv,1,'Yv')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs().or.restarted) call viz%write(time=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%Yvmin,'Yvmin')
         call mfile%add_column(fs%Yvmax,'Yvmax')
         call mfile%add_column(fs%Qint(8),'Vapor mass')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%Qint(8),'Vapor Mass')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_sl,'sl_max')
         call tfile%add_column(fs%wtmin_sl,'sl_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         if (time%n.le.20) then
            time%dt=1e-9_WP
         else
            call time%adjust_dt()
         end if
         call time%increment()

         ! Remember old state
         call fs%Qold%copy(src=fs%Q)
         call fs%store_old()

         ! ===== RK2 Stage 1: dQdt = f(t, Q) =====
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=0.5_WP*time%dt,time=time%t)

         ! ===== RK2 Stage 2: Q* = Qold + dt/2*dQdt, dQdt* = f(t+dt/2, Q*) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%apply_relax(time=time%t+0.5_WP*time%dt)
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=time%dt,time=time%t+0.5_WP*time%dt)

         ! ===== RK2 Final: Q = Qold + dt*dQdt* =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%apply_relax(time=time%t)

         ! Rebuild PLIC
         call fs%build_plic(time%t)

         ! Recompute primitive variables
         call fs%get_primitive(fs%Q)

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         if (fs%amr%amRoot) print*, 'time: ', time%t, ' dt: ', time%dt,' occurs: ',viz_evt%occurs()
         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time%t)
         end if
         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/Sembian_blastwave_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
