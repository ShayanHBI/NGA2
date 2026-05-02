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
   real(WP) :: Yvmin_half_transport=0.0_WP,Yvmax_half_transport=0.0_WP
   real(WP) :: Yvmin_half_relax=0.0_WP,Yvmax_half_relax=0.0_WP
   real(WP) :: Yvmin_full_transport=0.0_WP,Yvmax_full_transport=0.0_WP
   real(WP) :: Yvmin_full_relax=0.0_WP,Yvmax_full_relax=0.0_WP
   real(WP) :: dYvmax_half_relax=0.0_WP,Yvpre_half_relax=0.0_WP,Yvpost_half_relax=0.0_WP
   real(WP) :: VF_dYv_half_relax=0.0_WP,p_dYv_half_relax=0.0_WP,superheat_dYv_half_relax=0.0_WP
   real(WP) :: Q2_dYv_half_relax=0.0_WP,Q8_dYv_half_relax=0.0_WP
   real(WP) :: dYvmax_full_relax=0.0_WP,Yvpre_full_relax=0.0_WP,Yvpost_full_relax=0.0_WP
   real(WP) :: VF_dYv_full_relax=0.0_WP,p_dYv_full_relax=0.0_WP,superheat_dYv_full_relax=0.0_WP
   real(WP) :: Q2_dYv_full_relax=0.0_WP,Q8_dYv_full_relax=0.0_WP

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

   !> Rate-limit PTg relaxation diagnostics in parallel runs
   integer :: PTg_fail_prints=0
   integer, parameter :: PTg_fail_print_max=40
   integer :: PTg_large_dYv_prints=0
   integer, parameter :: PTg_large_dYv_print_max=80
   real(WP), parameter :: PTg_large_dYv_threshold=2.0e-2_WP
   real(WP) :: PTg_diag_dYvmax=0.0_WP,PTg_diag_Yvpre=0.0_WP,PTg_diag_Yvpost=0.0_WP
   real(WP) :: PTg_diag_VF=0.0_WP,PTg_diag_p=0.0_WP,PTg_diag_superheat=0.0_WP
   real(WP) :: PTg_diag_Q2=0.0_WP,PTg_diag_Q8=0.0_WP

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

   !> Thermo-chemical relaxation model (Pure liquid and gas mixture)
   subroutine PTg_relax(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: PG,PL,ZG,ZL,Pint
      real(WP) :: a,b,d,coeffL,coeffG
      real(WP) :: VFeq,Peq,p,ppv,Yv,xv,T,Teq,Tsat,pOld,YvOld,Yv_chem0
      real(WP) :: RHOL,RHOG
      real(WP) :: rho0,Eps0,rhoA0
      real(WP) :: Yv_tol,Yv_err
      real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
      real(WP) :: p_pert,Yv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
      real(WP) :: step_fac,p_try,Yv_try,xv_try,ppv_try
      real(WP) :: F1_try,F2_try,res_norm,res_norm_try,Eps_scale
      real(WP) :: Tlo,Thi,Tmid,Flo,Fhi,Fmid
      real(WP) :: Yv_hi
      real(WP) :: VF_input
      real(WP), dimension(8) :: Q_input
      character(len=64) :: fail_reason
      real(WP), parameter :: fd_eps=1.0e-8_WP,lnP_eps=1.0e-12_WP,Yv_min=1.0e-12_WP,Yv_max=1.0_WP-Yv_min
      real(WP), parameter :: chem_VF_min=1.0e-5_WP
      real(WP), parameter :: mech_VF_min=1.0e-4_WP
      real(WP), parameter :: chem_Ya_min=1.0e-8_WP
      real(WP), parameter :: chem_air_mass_frac_min=1.0e-6_WP
      real(WP), parameter :: chem_mass_margin=1.0e-10_WP
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      integer  :: it,itmax
      real(WP) :: p_tol,p_err,T_tol,T_err
      real(WP), parameter :: F1_tol=1.0e-3_WP,F2_rel_tol=1.0e-4_WP
      logical  :: converge
      VF_input=VF
      Q_input=Q(1:8)
      ! Set gas EoS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      call set_gas_eos_cof(Yv)
      ! Remove numerically meaningless gas slivers before evaluating the gas EOS.
      ! These states usually come from tiny transport undershoots near pure liquid.
      if ((1.0_WP-VF).le.chem_VF_min.or.Q(2).le.0.0_WP.or.Q(4).le.0.0_WP.or.Q(2)/(1.0_WP-VF).lt.RHOGmin) then
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP
         return
      end if
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
         if ((1.0_WP-VF).le.mech_VF_min) then
            VF=1.0_WP
            Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
            Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
            Q(8)=0.0_WP
            return
         end if
         print*,"****************** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; Q(8)=0.0_WP; return
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
      ! ================= Third step for chemical relaxation =================
      ! Return if there is no interface left after thermal relaxation
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP) return
      ! The p-Yv equilibrium solve is ill-conditioned when either phase is a
      ! volume sliver; keep pressure/thermal relaxation, but skip chemistry.
      if (VF.le.chem_VF_min.or.(1.0_WP-VF).le.chem_VF_min) return
      p=Peq
      ! Get the thermo-mechanically relaxed temperature
      Teq=get_TL(RHO=Q(1)/VF,P=p)
      ! Get vapor mole fraction and partial pressure
      xv=get_xv(Yv)
      ppv=xv*p
      if (.not.valid_sat_state(ppv,Teq)) then
         print*,"****************** Vapor partial pressure too low. Skipping the cell!"
         return
      end if
      ! Get saturation temperature from vapor partial pressure
      itmax=20
      T_tol=1e-7_WP
      converge=.false.
      Tsat=Teq
      do it=1,itmax
         T=Tsat
         if (.not.valid_sat_state(ppv,T)) return
         Tsat=T-PTsat(ppv,T)/dPTsatdT(T)
         T_err=abs((Tsat-T)/T)
         if (T_err.lt.T_tol) then
            converge=.true.
            exit
         end if
      end do
      ! If Newton fails, fall back to a bracketed solve in temperature.
      if (.not.converge) then
         Tlo=max(lnP_eps,100.0_WP)
         Thi=max(2000.0_WP,2.0_WP*abs(Teq)+100.0_WP)
         Flo=PTsat(ppv,Tlo)
         Fhi=PTsat(ppv,Thi)
         do it=1,20
            if (Flo*Fhi.le.0.0_WP) exit
            Tlo=max(lnP_eps,0.5_WP*Tlo)
            Thi=2.0_WP*Thi
            Flo=PTsat(ppv,Tlo)
            Fhi=PTsat(ppv,Thi)
         end do
         if (Flo*Fhi.le.0.0_WP) then
            do it=1,80
               Tmid=0.5_WP*(Tlo+Thi)
               Fmid=PTsat(ppv,Tmid)
               if (abs(Thi-Tlo)/max(abs(Tmid),1.0_WP).lt.T_tol) then
                  Tsat=Tmid
                  converge=.true.
                  exit
               end if
               if (Flo*Fmid.le.0.0_WP) then
                  Thi=Tmid
                  Fhi=Fmid
               else
                  Tlo=Tmid
                  Flo=Fmid
               end if
            end do
         end if
      end if
      if (.not.converge) then
         print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
         return
      end if
      ! Activate chemical relaxation only for metastable states
      if (Teq.le.Tsat) return
      ! Get the initial quantities
      rho0=sum(Q(1:2))
      Eps0=sum(Q(3:4))
      rhoA0=(1.0_WP-Yv)*Q(2)
      ! This mixture formulation uses conserved air mass through
      ! rhoA0/(1-Yv); when the gas sliver is essentially pure vapor, the
      ! p-Yv solve becomes singular and there is no air reservoir left.
      if ((1.0_WP-Yv).le.chem_Ya_min.or.rhoA0.le.chem_air_mass_frac_min*max(rho0,1.0_WP)) return
      ! Conserved air mass also imposes rhoA0/(1-Yv) < rho0.  Enforce that
      ! physical upper bound directly instead of chasing it with thresholds.
      Yv_hi=min(Yv_max,1.0_WP-rhoA0/(rho0*(1.0_WP-chem_mass_margin)))
      if (Yv_hi.le.Yv_min.or.Yv.ge.Yv_hi) return
      Yv_chem0=Yv
      ! Iteratively solve for the equilibrium pressure and vapor mass fraction
      itmax=60
      p_tol=1e-7_WP
      Yv_tol=1e-7_WP
      Eps_scale=max(abs(Eps0),1.0_WP)
      converge=.false.
      fail_reason='iteration limit'
      do it=1,itmax
         ! Evaluate residuals at current state
         call get_T(p,Yv)
         xv=get_xv(Yv)
         ppv=xv*p
         if (.not.valid_sat_state(ppv,T)) then
            fail_reason='invalid current saturation state'
            exit
         end if
         F1=PTsat(ppv,T)
         F2=Eps_res(p,Yv)
         res_norm=F1**2+(F2/Eps_scale)**2
         ! Compute Jacobian via finite differences
         ! --- Perturbation in p ---
         p_pert=p+fd_eps*max(abs(p),1.0_WP)
         ppv=xv*p_pert
         call get_T(p_pert,Yv)
         if (.not.valid_sat_state(ppv,T)) then
            fail_reason='invalid pressure perturbation'
            exit
         end if
         F1p=PTsat(ppv,T)
         F2p=Eps_res(p_pert,Yv)
         dF1dp=(F1p-F1)/(p_pert-p)
         dF2dp=(F2p-F2)/(p_pert-p)
         ! --- Perturbation in Yv ---
         Yv_pert=min(Yv_hi,Yv+fd_eps)
         if (abs(Yv_pert-Yv).lt.tiny(1.0_WP)) Yv_pert=max(Yv_min,Yv-fd_eps)
         xv=get_xv(Yv_pert)
         ppv=xv*p
         call get_T(p,Yv_pert)
         if (.not.valid_sat_state(ppv,T)) then
            fail_reason='invalid Yv perturbation'
            exit
         end if
         F1Y=PTsat(ppv,T)
         F2Y=Eps_res(p,Yv_pert)
         dF1dYv=(F1Y-F1)/(Yv_pert-Yv)
         dF2dYv=(F2Y-F2)/(Yv_pert-Yv)
         ! Solve 2x2 system: J * [dp; dYv] = -[F1; F2]
         detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
         if (abs(detJ).lt.1.0e-30_WP) then
            fail_reason='singular Jacobian'
            exit
         end if
         dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
         dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
         ! Newton-Raphson update
         pOld=p
         YvOld=Yv
         step_fac=1.0_WP
         do while (step_fac.gt.1.0e-6_WP)
            p_try=pOld+step_fac*dp_nr
            Yv_try=max(Yv_min,min(Yv_hi,YvOld+step_fac*dYv_nr))
            call get_T(p_try,Yv_try)
            xv_try=get_xv(Yv_try)
            ppv_try=xv_try*p_try
            if (valid_sat_state(ppv_try,T)) then
               F1_try=PTsat(ppv_try,T)
               F2_try=Eps_res(p_try,Yv_try)
               res_norm_try=F1_try**2+(F2_try/Eps_scale)**2
               if (res_norm_try.lt.res_norm) exit
            end if
            step_fac=0.5_WP*step_fac
         end do
         if (step_fac.le.1.0e-6_WP) then
            fail_reason='line search failed'
            exit
         end if
         p=p_try
         Yv=Yv_try
         ! Evaluate convergence
         p_err=abs((p-pOld)/max(abs(pOld),1.0_WP))
         Yv_err=abs((Yv-YvOld)/max(abs(YvOld),Yv_min))
         if (((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol)).or. &
         &   (abs(F1_try).lt.F1_tol.and.abs(F2_try).lt.F2_rel_tol*Eps_scale)) then
            converge=.true.
            exit
         end if
      end do
      if (.not.converge.and.abs(F1).lt.F1_tol.and.abs(F2).lt.F2_rel_tol*Eps_scale) converge=.true.
      if (.not.converge) then
         call print_PTG_failure(trim(fail_reason),it,VF,Q,p,Yv,T,ppv,F1,F2,detJ,res_norm,dp_nr,dYv_nr, &
         & step_fac,Teq,Tsat,rho0,Eps0,rhoA0)
         print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
         return
      end if
      ! Get equilibrium quantities at converged (p, Yv)
      call get_T(p,Yv)
      RHOL=get_RHOL(p,T)
      call set_gas_eos_cof(Yv)
      RHOG=get_RHOG(p,T,Yv)
      VFeq=(rho0-RHOG)/(RHOL-RHOG)
      ! Do not apply clipped chemical states; the pure-phase limit needs a
      ! separate conservative treatment.
      if (VFeq.le.chem_VF_min.or.VFeq.ge.1.0_WP-chem_VF_min) return
      ! Adjust conserved quantities
      call print_PTG_large_dYv(VF_input,Q_input,Yv_chem0,Yv,VFeq,p,T,Tsat, &
      & (1.0_WP-VFeq)*RHOG,(1.0_WP-VFeq)*RHOG*Yv)
      call update_PTG_diag(Yv_chem0,Yv,VFeq,p,T,Tsat,(1.0_WP-VFeq)*RHOG,(1.0_WP-VFeq)*RHOG*Yv)
      Q(1)=(       VFeq)*RHOL
      Q(2)=(1.0_WP-VFeq)*RHOG
      Q(3)=Q(1)*get_IL_PT(p,T)
      Q(4)=Q(2)*get_IG_PT(p,T,Yv)
      Q(8)=Q(2)*Yv
      VF=VFeq
      contains
      ! Equilibrium temperature as a function of pressure and vapor mass fraction
      subroutine get_T(p_in,Yv_in)
         real(WP), intent(in) :: p_in,Yv_in
         T=1.0_WP/((rho0-rhoA0/(1.0_WP-Yv_in))*(GammaL-1.0_WP)*CvL/(p_in+PinfL)+rhoA0/(1.0_WP-Yv_in)*((GammaV-1.0_WP)*CvV*Yv_in/(p_in+PinfV)+(GammaA-1.0_WP)*CvA*(1.0_WP-Yv_in)/(p_in+PinfA)))
      end subroutine get_T
      ! Function that defines p-T saturation curve
      function PTsat(p_in,T_in)
         real(WP), intent(in) :: p_in,T_in
         real(WP) :: PTsat
         PTsat=AS+BS/T_in+CS*log(T_in)+DS*log(p_in+PinfL)-log(p_in+PinfV)
      end function PTsat
      ! Saturation-law domain guard; this rejects invalid Newton candidates
      ! without imposing a pressure floor on the mechanical relaxation state.
      logical function valid_sat_state(p_in,T_in)
         real(WP), intent(in) :: p_in,T_in
         valid_sat_state=((p_in+PinfV).gt.lnP_eps).and.((p_in+PinfL).gt.lnP_eps)
      end function valid_sat_state
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
         Eps_res=(rho0-rhoA0/(1.0_WP-Yv_in))*(CvL*T*(p_in+GammaL*PinfL)/(p_in+PinfL)+etaL)+rhoA0/(1.0_WP-Yv_in)*((CvV*T*(p_in+GammaV*PinfV)/(p_in+PinfV)+etaV)*Yv_in+(CvA*T*(p_in+GammaA*PinfA)/(p_in+PinfA)+etaA)*(1.0_WP-Yv_in))-Eps0
      end function Eps_res
   end subroutine PTg_relax

   subroutine print_PTG_failure(reason,it,VF,Q,p,Yv,T,ppv,F1,F2,detJ,res_norm,dp_nr,dYv_nr, &
   & step_fac,Teq,Tsat,rho0,Eps0,rhoA0)
      implicit none
      character(len=*), intent(in) :: reason
      integer, intent(in) :: it
      real(WP), intent(in) :: VF,p,Yv,T,ppv,F1,F2,detJ,res_norm,dp_nr,dYv_nr
      real(WP), intent(in) :: step_fac,Teq,Tsat,rho0,Eps0,rhoA0
      real(WP), dimension(1:), intent(in) :: Q
      real(WP) :: PL_dbg,PG_dbg,rhoL_dbg,rhoG_dbg,IL_dbg,IG_dbg
      if (PTg_fail_prints.ge.PTg_fail_print_max) return
      PTg_fail_prints=PTg_fail_prints+1
      rhoL_dbg=0.0_WP; rhoG_dbg=0.0_WP; IL_dbg=0.0_WP; IG_dbg=0.0_WP; PL_dbg=0.0_WP; PG_dbg=0.0_WP
      if (VF.gt.0.0_WP.and.Q(1).gt.0.0_WP) then
         rhoL_dbg=Q(1)/VF
         IL_dbg=Q(3)/Q(1)
         PL_dbg=get_PL(rhoL_dbg,IL_dbg)
      end if
      if ((1.0_WP-VF).gt.0.0_WP.and.Q(2).gt.0.0_WP) then
         rhoG_dbg=Q(2)/(1.0_WP-VF)
         IG_dbg=Q(4)/Q(2)
         PG_dbg=get_PG(rhoG_dbg,IG_dbg,Yv)
      end if
      write(*,'(A,I0,3A,I0)') 'PTg failure #',PTg_fail_prints,': ',trim(reason),' it=',it
      write(*,'(A,12ES15.7)') '  VF p Yv T ppv Teq Tsat rho0 Eps0 rhoA0 detJ res = ', &
      & VF,p,Yv,T,ppv,Teq,Tsat,rho0,Eps0,rhoA0,detJ,res_norm
      write(*,'(A,10ES15.7)') '  F1 F2 dp dYv step Q1 Q2 Q3 Q4 Q8 = ', &
      & F1,F2,dp_nr,dYv_nr,step_fac,Q(1),Q(2),Q(3),Q(4),Q(8)
      write(*,'(A,6ES15.7)') '  rhoL rhoG IL IG PL PG = ',rhoL_dbg,rhoG_dbg,IL_dbg,IG_dbg,PL_dbg,PG_dbg
   end subroutine print_PTG_failure

   subroutine print_PTG_large_dYv(VF_input,Q_input,Yvpre,Yvpost,VFout,p,T,Tsat,Q2out,Q8out)
      implicit none
      real(WP), intent(in) :: VF_input,Yvpre,Yvpost,VFout,p,T,Tsat,Q2out,Q8out
      real(WP), dimension(8), intent(in) :: Q_input
      real(WP) :: dYv,Yv_input
      dYv=Yvpost-Yvpre
      if (dYv.lt.PTg_large_dYv_threshold) return
      if (PTg_large_dYv_prints.ge.PTg_large_dYv_print_max) return
      PTg_large_dYv_prints=PTg_large_dYv_prints+1
      if (Q_input(2).gt.0.0_WP) then
         Yv_input=Q_input(8)/Q_input(2)
      else
         Yv_input=0.0_WP
      end if
      write(*,'(A,I0,A)') 'PTg large dYv #',PTg_large_dYv_prints,': successful chemical relaxation'
      write(*,'(A,6ES15.7)') '  input VF Yv dYv Yvpre Yvpost T-Tsat = ', &
      & VF_input,Yv_input,dYv,Yvpre,Yvpost,T-Tsat
      write(*,'(A,8ES15.7)') '  input Q(1:8) = ',Q_input
      write(*,'(A,5ES15.7)') '  output VF p T Tsat Q2 = ',VFout,p,T,Tsat,Q2out
      write(*,'(A,1ES15.7)') '  output Q8 = ',Q8out
   end subroutine print_PTG_large_dYv

   subroutine reset_PTG_diag()
      implicit none
      PTg_diag_dYvmax=0.0_WP
      PTg_diag_Yvpre=0.0_WP
      PTg_diag_Yvpost=0.0_WP
      PTg_diag_VF=0.0_WP
      PTg_diag_p=0.0_WP
      PTg_diag_superheat=0.0_WP
      PTg_diag_Q2=0.0_WP
      PTg_diag_Q8=0.0_WP
   end subroutine reset_PTG_diag

   subroutine update_PTG_diag(Yvpre,Yvpost,VF,p,T,Tsat,Q2,Q8)
      implicit none
      real(WP), intent(in) :: Yvpre,Yvpost,VF,p,T,Tsat,Q2,Q8
      real(WP) :: dYv
      dYv=Yvpost-Yvpre
      if (dYv.le.PTg_diag_dYvmax) return
      PTg_diag_dYvmax=dYv
      PTg_diag_Yvpre=Yvpre
      PTg_diag_Yvpost=Yvpost
      PTg_diag_VF=VF
      PTg_diag_p=p
      PTg_diag_superheat=T-Tsat
      PTg_diag_Q2=Q2
      PTg_diag_Q8=Q8
   end subroutine update_PTG_diag

   subroutine collect_PTG_diag(dYvmax,Yvpre,Yvpost,VF,p,superheat,Q2,Q8)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_ALLGATHER,MPI_COMM_SIZE
      implicit none
      real(WP), intent(out) :: dYvmax,Yvpre,Yvpost,VF,p,superheat,Q2,Q8
      real(WP), dimension(8) :: sendbuf
      real(WP), dimension(:,:), allocatable :: recvbuf
      integer :: nproc,ierr,owner
      sendbuf=[PTg_diag_dYvmax,PTg_diag_Yvpre,PTg_diag_Yvpost,PTg_diag_VF, &
      &        PTg_diag_p,PTg_diag_superheat,PTg_diag_Q2,PTg_diag_Q8]
      call MPI_COMM_SIZE(amr%comm,nproc,ierr)
      allocate(recvbuf(8,nproc))
      call MPI_ALLGATHER(sendbuf,8,MPI_REAL_WP,recvbuf,8,MPI_REAL_WP,amr%comm,ierr)
      owner=maxloc(recvbuf(1,:),dim=1)
      dYvmax  =recvbuf(1,owner)
      Yvpre   =recvbuf(2,owner)
      Yvpost  =recvbuf(3,owner)
      VF      =recvbuf(4,owner)
      p       =recvbuf(5,owner)
      superheat=recvbuf(6,owner)
      Q2      =recvbuf(7,owner)
      Q8      =recvbuf(8,owner)
      deallocate(recvbuf)
   end subroutine collect_PTG_diag

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
      call param_read('Initial vapor mass fraction',Yv0,default=6.3e-1_WP)
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
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
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
         Yvmin_half_transport=fs%Yvmin; Yvmax_half_transport=fs%Yvmax
         Yvmin_half_relax    =fs%Yvmin; Yvmax_half_relax    =fs%Yvmax
         Yvmin_full_transport=fs%Yvmin; Yvmax_full_transport=fs%Yvmax
         Yvmin_full_relax    =fs%Yvmin; Yvmax_full_relax    =fs%Yvmax
         dYvmax_half_relax=0.0_WP; Yvpre_half_relax=fs%Yvmax; Yvpost_half_relax=fs%Yvmax
         VF_dYv_half_relax=0.0_WP; p_dYv_half_relax=0.0_WP; superheat_dYv_half_relax=0.0_WP
         Q2_dYv_half_relax=0.0_WP; Q8_dYv_half_relax=0.0_WP
         dYvmax_full_relax=0.0_WP; Yvpre_full_relax=fs%Yvmax; Yvpost_full_relax=fs%Yvmax
         VF_dYv_full_relax=0.0_WP; p_dYv_full_relax=0.0_WP; superheat_dYv_full_relax=0.0_WP
         Q2_dYv_full_relax=0.0_WP; Q8_dYv_full_relax=0.0_WP
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
         call mfile%add_column(fs%Yvmin,'Yvmin')
         call mfile%add_column(fs%Yvmax,'Yvmax')
         call mfile%add_column(Yvmin_half_transport,'Yvmin half transport')
         call mfile%add_column(Yvmax_half_transport,'Yvmax half transport')
         call mfile%add_column(Yvmin_half_relax,'Yvmin half relax')
         call mfile%add_column(Yvmax_half_relax,'Yvmax half relax')
         call mfile%add_column(Yvmin_full_transport,'Yvmin full transport')
         call mfile%add_column(Yvmax_full_transport,'Yvmax full transport')
         call mfile%add_column(Yvmin_full_relax,'Yvmin full relax')
         call mfile%add_column(Yvmax_full_relax,'Yvmax full relax')
         call mfile%add_column(dYvmax_half_relax,'max dYv half relax')
         call mfile%add_column(Yvpre_half_relax,'Yv before max dYv half')
         call mfile%add_column(Yvpost_half_relax,'Yv after max dYv half')
         call mfile%add_column(VF_dYv_half_relax,'VF max dYv half')
         call mfile%add_column(p_dYv_half_relax,'p max dYv half')
         call mfile%add_column(superheat_dYv_half_relax,'T-Tsat max dYv half')
         call mfile%add_column(Q2_dYv_half_relax,'Q2 max dYv half')
         call mfile%add_column(Q8_dYv_half_relax,'Q8 max dYv half')
         call mfile%add_column(dYvmax_full_relax,'max dYv full relax')
         call mfile%add_column(Yvpre_full_relax,'Yv before max dYv full')
         call mfile%add_column(Yvpost_full_relax,'Yv after max dYv full')
         call mfile%add_column(VF_dYv_full_relax,'VF max dYv full')
         call mfile%add_column(p_dYv_full_relax,'p max dYv full')
         call mfile%add_column(superheat_dYv_full_relax,'T-Tsat max dYv full')
         call mfile%add_column(Q2_dYv_full_relax,'Q2 max dYv full')
         call mfile%add_column(Q8_dYv_full_relax,'Q8 max dYv full')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
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
         call record_Yv_stage(Yvmin_half_transport,Yvmax_half_transport)
         call reset_PTG_diag()
         call fs%apply_relax(time=time%t+0.5_WP*time%dt)
         call collect_PTG_diag(dYvmax_half_relax,Yvpre_half_relax,Yvpost_half_relax,VF_dYv_half_relax, &
         & p_dYv_half_relax,superheat_dYv_half_relax,Q2_dYv_half_relax,Q8_dYv_half_relax)
         call record_Yv_stage(Yvmin_half_relax,Yvmax_half_relax)
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=time%dt,time=time%t+0.5_WP*time%dt)

         ! ===== RK2 Final: Q = Qold + dt*dQdt* =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call record_Yv_stage(Yvmin_full_transport,Yvmax_full_transport)
         call reset_PTG_diag()
         call fs%apply_relax(time=time%t)
         call collect_PTG_diag(dYvmax_full_relax,Yvpre_full_relax,Yvpost_full_relax,VF_dYv_full_relax, &
         & p_dYv_full_relax,superheat_dYv_full_relax,Q2_dYv_full_relax,Q8_dYv_full_relax)
         call record_Yv_stage(Yvmin_full_relax,Yvmax_full_relax)

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

   !> Record Yv extrema for the current conserved state.
   subroutine record_Yv_stage(Yvmin,Yvmax)
      implicit none
      real(WP), intent(out) :: Yvmin,Yvmax
      call fs%get_primitive(fs%Q)
      call fs%get_info()
      Yvmin=fs%Yvmin
      Yvmax=fs%Yvmax
   end subroutine record_Yv_stage

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
