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
   integer  :: restart_step

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> EOS parameters (NASG liquid and ideal gas)
   real(WP) :: GammaL,PinfL,qL,qpL,CvL,CpL,bL
   real(WP) :: GammaG,qG,qpG,CvG,CpG
   real(WP) :: GammaA,qA,qpA,CvA,CpA
   real(WP) :: GammaV,qV,qpV,CvV,CpV

   !> Saturation curve coefficients
   real(WP) :: AS,BS,CS,DS,ES

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
   real(WP) :: PrL,PrG,ScV        !< Prandtle and Schmidt numbers

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

   !> Time stepping
   real(WP) :: dt_init

contains

   !> Levelset function for 2D cylinder centered at (xcyl, 0)
   function levelset_cyl(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2+xyz(3)**2)
      if (amr%nz.eq.1) G=0.5_WP*dcyl-sqrt((xyz(1)-xcyl)**2+xyz(2)**2) ! Enable quasi-2D runs
   end function levelset_cyl

   !> Liquid EOS: NASG
   !> P=f(RHO,I)
   real(WP) function get_PL(RHO,I,Yv)
      implicit none
      real(WP), intent(in) :: RHO,I
      real(WP), intent(in), optional :: Yv
      get_PL=(GammaL-1.0_WP)*RHO*(I-qL)/(1.0_WP-bL*RHO)-GammaL*PinfL
   end function get_PL
   !> T=f(RHO,P)
   real(WP) function get_TL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_TL=(P+PinfL)*(1.0_WP-bL*RHO)/((GammaL-1.0_WP)*CvL*RHO)
   end function get_TL
   !> C=f(RHO,P)
   real(WP) function get_CL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_CL=sqrt(max(0.0_WP,GammaL*(P+PinfL)/(RHO*(1.0_WP-bL*RHO))))
   end function get_CL
   !> I=f(RHO,P)
   real(WP) function get_IL(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      get_IL=(1.0_WP-bL*RHO)*(P+GammaL*PinfL)/((GammaL-1.0_WP)*RHO)+qL
   end function get_IL
   !> I=f(P,T)
   real(WP) function get_IL_PT(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_IL_PT=CvL*T*(P+GammaL*PinfL)/(P+PinfL)+qL
   end function get_IL_PT
   !> RHOLIL=f(P,RHO)
   real(WP) function get_RHOLIL(P,RHO)
      implicit none
      real(WP), intent(in) :: P,RHO
      get_RHOLIL=(1.0_WP-bL*RHO)*(P+GammaL*PinfL)/(GammaL-1.0_WP)+RHO*qL
   end function get_RHOLIL
   !> RHOLIL=f(P,T)
   real(WP) function get_RHOLIL_PT(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_RHOLIL_PT=(P+GammaL*PinfL)*CvL*T+qL*(P+PinfL)/((GammaL-1.0_WP)*CvL*T+bL*(P+PinfL))
   end function get_RHOLIL_PT
   !> RHOL=f(P,T)
   real(WP) function get_RHOL(P,T)
      implicit none
      real(WP), intent(in) :: P,T
      get_RHOL=(P+PinfL)/((GammaL-1.0_WP)*CvL*T+bL*(P+PinfL))
   end function get_RHOL
   !> GL=f(RHO)
   real(WP) function get_GL(RHO)
      implicit none
      real(WP), intent(in) :: RHO
      get_GL=(GammaL-1.0_WP)/(1.0_WP-bL*RHO)
   end function get_GL

   !> Gas EOS: Ideal Gas
   !> P=f(RHO,I)
   real(WP) function get_PG(RHO,I,Yv)
      implicit none
      real(WP), intent(in) :: RHO,I
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_PG=RHO*I*(GammaG-1.0_WP)-(GammaG-1.0_WP)*qG*RHO
   end function get_PG
   !> T=f(RHO,P)
   real(WP) function get_TG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_TG=P/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> C=f(RHO,P)
   real(WP) function get_CG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_CG=sqrt(max(0.0_WP,GammaG*P/RHO))
   end function get_CG
   !> I=f(RHO,P)
   real(WP) function get_IG(RHO,P,Yv)
      implicit none
      real(WP), intent(in) :: RHO,P
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_IG=P/(RHO*(GammaG-1.0_WP))+qG
   end function get_IG
   !> I=f(P,T,Yv)
   real(WP) function get_IG_PT(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_IG_PT=CvG*T+qG
   end function get_IG_PT
   !> RHOGIG=f(P,RHO,Yv)
   real(WP) function get_RHOGIG(P,RHO,Yv)
      implicit none
      real(WP), intent(in) :: P,RHO
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_RHOGIG=P/(GammaG-1.0_WP)+RHO*qG
   end function get_RHOGIG
   !> RHOGIG=f(P,T,Yv)
   real(WP) function get_RHOGIG_PT(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_RHOGIG_PT=P*(CvG*T+qG)/((GammaG-1.0_WP)*CvG*T)
   end function get_RHOGIG_PT
   !> RHOG=f(P,T,Yv)
   real(WP) function get_RHOG(P,T,Yv)
      implicit none
      real(WP), intent(in) :: P,T
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_RHOG=P/((GammaG-1.0_WP)*CvG*T)
   end function get_RHOG
   !> GG=f(Yv)
   real(WP) function get_GG(Yv)
      implicit none
      real(WP), intent(in), optional :: Yv
      if (present(Yv)) then
         call set_gas_eos_cof(Yv)
      else
         call set_gas_eos_cof(0.0_WP)
      end if
      get_GG=GammaG-1.0_WP
   end function get_GG
   !> Set gas EOS coefficients from vapor mass fraction (Ideal gas mixture of vapor and air)
   subroutine set_gas_eos_cof(Yv)
      implicit none
      real(WP), intent(in) :: Yv
      real(WP) :: Ya
      Ya    =1.0_WP-Yv
      CvG   =Yv*CvV+Ya*CvA
      CpG   =Yv*CpV+Ya*CpA
      qG    =Yv*qV +Ya*qA
      qpG   =Yv*qpV+Ya*qpA
      GammaG=CpG/CvG
   end subroutine set_gas_eos_cof
   !> Vapor mole fraction (Ideal gas mixture)
   real(WP) function get_xv(Yv)
      implicit none
      real(WP), intent(in) :: Yv
      get_xv=Yv*Ma/(Yv*Ma+(1.0_WP-Yv)*Mv)
   end function get_xv

   !> Chemical relaxation
   subroutine PTg_relax(VF,Q,Pjump)
      use messager, only: die
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable    :: Q0
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,PL,PG,TL,TG,ZL,ZG,PHIL,PHIG,Pint
      real(WP) :: xiL,xiG,xiLinv,xiGinv
      real(WP) :: xiTL,xiTG,xiTLinv,xiTGinv,zetaL,zetaG
      real(WP) :: Z,D,COF
      real(WP) :: VFeq,VF0,Peq,p,T,Teq,Yv
      real(WP) :: rho0,Eps0,rhoA0
      real(WP) :: rho_err,Eps_err
      real(WP), parameter :: lnP_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,ppv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: p_tol=1.0e-4_WP,Yv_tol=1.0e-4_WP,Tsat_tol=1.0e-5_WP,rho_tol=1.0e-4_WP,Eps_tol=1.0e-4_WP,F1_tol=1e-4_WP,F2_tol=1e-4_WP
      integer,  parameter :: Tsat_itmax=40,NR_itmax=40
      logical  :: chem_relax
      ! Store the input state
      allocate(Q0(size(Q)))
      VF0=VF
      Q0=Q
      ! Set gas EoS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      call set_gas_eos_cof(Yv)
      ! ================ First step for mechanical relaxation ================
      ! Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
      ! Get phasic quantities
      RHOL=Q(1)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF)
      PL=get_PL(RHO=RHOL,I=Q(3)/Q(1))
      PG=get_PG(RHO=RHOG,I=Q(4)/Q(2),Yv=Yv)
      CL=get_CL(RHO=RHOL,P=PL)
      CG=get_CG(RHO=RHOG,P=PG,Yv=Yv)
      ! Handle limit cases-should mass/energy be transfered or lost?-this should probably never happen...
      if (PL.le.-PinfL) then
         print*,"****************** LIQUID CLIPPED!",PL,VF,Q,Yv
         VF=0.0_WP; Q(2)=sum(Q(1:2)); Q(1)=0.0_WP; Q(4)=sum(Q(3:4)); Q(3)=0.0_WP; Q(8)=Yv*Q(2); call deallocate_Q0(); return
      end if
      if (PG.le.0.0_WP) then
         print*,"****************** GAS CLIPPED!",PG,VF,Q,Yv
         VF=1.0_WP; Q(1)=sum(Q(1:2)); Q(2)=0.0_WP; Q(3)=sum(Q(3:4)); Q(4)=0.0_WP; Q(8)=0.0_WP;  call deallocate_Q0(); return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*get_CL(RHO=Q(1)/(       VF),P=PL)
      ZG=Q(2)/(1.0_WP-VF)*get_CG(RHO=Q(2)/(1.0_WP-VF),P=PG,Yv=Yv)
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup the ODE coefficients
      GL=get_GL(RHO=RHOL)
      GG=get_GG(Yv=Yv)
      xiL=         VF/(GL*(Pint-PL)+RHOL*CL**2)
      xiG=(1.0_WP-VF)/(GG*(Pint-PG)+RHOG*CG**2)
      xiLinv=1.0_WP/xiL
      xiGinv=1.0_WP/xiG
      ! Get equilibrium volume fraction
      VFeq=VF-(PG-PL)/(xiLinv+xiGinv)
      ! Guard: if the linear formula drives VFeq outside [0,1] (strong disequilibrium,
      ! get_Peq_from_Econs (one_brho diverges as VFeq→0 or 1-VFeq→0).
      if ((VFeq.lt.0.0_WP).or.(VFeq.gt.1.0_WP)) then
         call restore_VFQ(); call deallocate_Q0(); return
      end if
      ! Get the equilibrium pressure
      Peq=get_Peq_from_Econs(VFeq)
      ! Guard: if energy-conserved Peq is unphysical (can happen when Q is already
      ! slightly corrupted by prior transport), skip rather than overwriting Q(3)/Q(4)
      ! with garbage.  Q has not been touched yet so no restore needed.
      if (Peq.le.max(0.0_WP,-PinfL)) then
         call restore_VFQ(); call deallocate_Q0(); return
      end if
      ! Adjust densities
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*get_RHOLIL(P=Peq,RHO=RHOL)
      Q(4)=(1.0_WP-VFeq)*get_RHOGIG(P=Peq,RHO=RHOG,Yv=Yv)
      ! ================= Second step for thermal relaxation =================
      ! Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
      ! Update the thermodynamic state
      TL=get_TL(P=Peq,RHO=RHOL)
      TG=get_TG(P=Peq,RHO=RHOG,Yv=Yv)
      GL=get_GL(RHOL)
      GG=get_GG(Yv)
      CL=get_CL(P=Peq,RHO=RHOL)
      CG=get_CG(P=Peq,RHO=RHOG,Yv=Yv)
      ! print '(A)',       '============ P_relax ============='
      ! print '(A,ES15.7)','p   =',Peq
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A,ES15.7)','TL   =',TL
      ! print '(A,ES15.7)','TG   =',TG
      ! print '(A)',       '==================================='
      ! Setup the ODE coefficients
      Z=(1.0_WP-VF)*GL+VF*GG
      D=VF*RHOG*CG**2+(1.0_WP-VF)*RHOL*CL**2
      PHIL=-(GammaL-1.0_WP)*CvL*RHOL**2/(Peq+PinfL)
      PHIG=-(GammaG-1.0_WP)*CvG*RHOG**2/Peq
      zetaL=RHOL*(1.0_WP-bL*RHOL)/(Peq+PinfL)
      zetaG=RHOG/Peq
      COF=GL*RHOG*CG**2-GG*RHOL*CL**2
      xiTL=-PHIL*D/(RHOL/(       VF)*Z+zetaL*COF)
      xiTG=-PHIG*D/(RHOG/(1.0_WP-VF)*Z-zetaG*COF)
      xiTLinv=1.0_WP/xiTL
      xiTGinv=1.0_WP/xiTG
      ! Get equilibrium VF and T
      VFeq=VF+Z/D*(TG-TL)/(xiTLinv+xiTGinv)
      Teq=(xiTL*TL+xiTG*TG)/(xiTL+xiTG)
      ! Get the equilibrium pressure
      Peq=get_Peq_from_Econs(VFeq)
      ! Check if pressure is sound; restore to pre-relaxation state if bad
      ! (Q0 still holds state A here — the Q0=Q overwrite at the start of step 3
      ! has not happened yet, so restore_VFQ() correctly undoes both step 1 and 2).
      if (Peq.le.max(0.0_WP,-PinfL)) then
         call deallocate_Q0(); return
      end if
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-PinfL); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,0.0_WP); end if
      ! Adjust densities
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*get_RHOLIL(P=Peq,RHO=RHOL)
      Q(4)=(1.0_WP-VFeq)*get_RHOGIG(P=Peq,RHO=RHOG,Yv=Yv)
      ! ================= Third step for chemical relaxation =================
      ! Store input state to the chemical relaxation algorithm
      Q0=Q
      VF0=VF
      p=Peq
      T=Teq
      rho0=sum(Q0(1:2))
      Eps0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! print '(A)',       '============ PT_relax ============='
      ! print '(A,ES15.7)','p   =',p
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A,ES15.7)','T   =',T
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
      ! print '(A,ES15.7)','p   =',p
      ! print '(A,ES15.7)','T   =',T
      ! print '(A,ES15.7)','Yv  =',Yv
      ! print '(A,ES15.7)','VF  =',VF
      ! print '(A)',       '=================================='
      ! Apply the converged state
      call set_gas_eos_cof(Yv)
      RHOL=get_RHOL(p,T)
      RHOG=get_RHOG(p,T,Yv)
      VF=(rho0-RHOG)/(RHOL-RHOG)
      ! Clean up solution
      if (VF.lt.0.0_WP) then; VF=0.0_WP; p=max(p,-PinfL); end if
      if (VF.gt.1.0_WP) then; VF=1.0_WP; p=max(p,0.0_WP); end if
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
         real(WP) function get_Peq_from_Econs(VF_eq)
            real(WP), intent(in) :: VF_eq
            real(WP) :: one_brho
            one_brho=1.0_WP-bL*Q0(1)/VF_eq
            get_Peq_from_Econs=(sum(Q0(3:4))-Q0(1)*qL-one_brho*VF_eq*GammaL*PinfL/(GammaL-1.0_WP)-Q0(2)*qG)/(one_brho*VF_eq/(GammaL-1.0_WP)+(1.0_WP-VF_eq)/(GammaG-1.0_WP))
         end function get_Peq_from_Econs
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
               ! Dry/nearly-dry air edge case: ppv is zero or so tiny that
               ! solving Tsat(ppv) is log-singular/ill-conditioned.  Seed Yv
               ! from saturation at the current thermally-relaxed state and
               ! then continue with the ordinary LVG Newton solve.
               ppv=exp(AS+(BS+ES*p_eq)/T_eq)*T_eq**CS*(p_eq+PinfL)**DS
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
               ! print '(A,ES15.7)','Seeded Yv from saturation at Teq=',Yv_eq
            else
               ! Get saturation temperature from total pressure and vapor partial
               ! pressure.  Use a safeguarded Newton solve rather than trusting T_eq
               ! as the initial guess.  At high-T mqstable states, T_eq can be far
               ! above the saturation temperature corresponding to (p_eq,ppv).
               call get_Tsat(p_eq,ppv,T_eq,Tsat,converge,Tsat_it)
               if (.not.converge) then
                  ! print*,"****************** Saturation temperature iterations blew up. Skipping the cell!!"
                  return
               end if
               ! print '(A)',       '========== Finding Tsat ==========='
               ! print '(A,I2)','Tsat it= ',Tsat_it
               ! print '(A,ES15.7)','Tsat   =',Tsat
               ! print '(A)',       '==================================='
               ! Activate chemical relaxation only for mqstable states
               ! print*,'T   =',T_eq
               ! print*,'Tsat=',Tsat
               if (T_eq.le.Tsat) return
            end if
            activate_chem_relax=.true.
         end function activate_chem_relax
         !> Safeguarded Newton solve for saturation temperature at fixed (p_l, p_v)
         !> The full multispecies saturation relation g_l(p_l,T)=g_v(p_v,T) is a
         !> function of two pressures, not one.  We hold both p_l and p_v fixed and
         !> solve for T.  In the pure-vapor branch p_v=p_l and the result reduces
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
               dFold=dPTsatdT(p_l,Told)
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
         !> Pure liquid and vapor chemical relaxation
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
               ! Pure-vapor branch: Y_v=1 so the vapor mole fraction x_v=1.
               pOld=p_eq
               p_eq=pOld-PTsat(pOld,pOld,T_eq)/dPTsatdp_lv(pOld,T_eq,dTdp)
               ! Evaluate the error
               p_err=abs((p_eq-pOld)/pOld)
               if (p_err.lt.p_tol) then
                  converged=.true.
                  exit
               end if
            end do
            ! print*,'NR iterations=',it
            if (.not.converged) then
               ! print*,"****************** p iterations blew up. Skipping the cell!!"
               return
            end if
            ! Update equilibrium temperature
            call get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp)
         end subroutine solve_PTg_relax_lv
         !> Pure liquid and vapor-gas mixture chemical relaxation
         ! subroutine solve_PTg_relax_lvg(p_eq,T_eq,Yv_eq,converged)
         !    real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         !    logical,  intent(out)   :: converged
         !    real(WP) :: xv,ppv
         !    real(WP) :: F1,F2,dF1dp,dF1dYv,dF2dp,dF2dYv,detJ
         !    real(WP) :: p_pert,Yv_pert,T_pert,xv_pert,ppv_pert,F1p,F2p,F1Y,F2Y,dp_nr,dYv_nr
         !    real(WP) :: pOld,YvOld,p_err,Yv_err
         !    real(WP) :: alpha,res0,res_try
         !    real(WP) :: p_try,Yv_try,T_try,xv_try,ppv_try,F1_try,F2_try
         !    real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
         !    integer  :: it
         !    logical  :: accepted
         !    ! Iteratively solve for the equilibrium pressure and vapor mass fraction
         !    converged=.false.
         !    p_err=10.0_WP*p_tol
         !    Yv_err=10.0_WP*Yv_tol
         !    do it=1,NR_itmax
         !       ! Evaluate residuals at current state
         !       T_eq=get_T_lvg(p_eq,Yv_eq)
         !       xv=get_xv(Yv_eq)
         !       ppv=xv*p_eq
         !       if (.not.check_p(ppv)) then
         !          ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
         !          return
         !       end if
         !       F1=PTsat(p_eq,ppv,T_eq)
         !       ! F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)
         !       ! res0=sqrt(F1**2+(F2/Eps0)**2)
         !       F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)/Eps0
         !       res0=sqrt(F1**2+F2**2)
         !       ! Compute Jacobian via finite difference
         !       ! --- Perturbation in p ---
         !       p_pert=p_eq*(1.0_WP+fd_eps)
         !       ppv_pert=xv*p_pert
         !       if (.not.check_p(ppv_pert)) then
         !          ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
         !          return
         !       end if
         !       T_pert=get_T_lvg(p_pert,Yv_eq)
         !       F1p=PTsat(p_pert,ppv_pert,T_pert)
         !       ! F2p=Eps_res_lvg(p_pert,T_pert,Yv_eq)
         !       F2p=Eps_res_lvg(p_pert,T_pert,Yv_eq)/Eps0
         !       dF1dp=(F1p-F1)/(p_pert-p_eq)
         !       dF2dp=(F2p-F2)/(p_pert-p_eq)
         !       ! --- Perturbation in Yv ---
         !       Yv_pert=Yv_eq+fd_eps
         !       if (Yv_pert.gt.Yvmax-fd_eps) then
         !          Yv_pert=Yv_eq-fd_eps
         !       end if
         !       if (Yv_pert.lt.Yvmin+fd_eps) then
         !          Yv_pert=Yv_eq+fd_eps
         !       end if
         !       xv_pert=get_xv(Yv_pert)
         !       ppv_pert=xv_pert*p_eq
         !       if (.not.check_p(ppv_pert)) then
         !          ! print*,"****************** Vapor partial pressure too low. Skipping the cell!"
         !          return
         !       end if
         !       T_pert=get_T_lvg(p_eq,Yv_pert)
         !       F1Y=PTsat(p_eq,ppv_pert,T_pert)
         !       ! F2Y=Eps_res_lvg(p_eq,T_pert,Yv_pert)
         !       F2Y=Eps_res_lvg(p_eq,T_pert,Yv_pert)/Eps0
         !       dF1dYv=(F1Y-F1)/(Yv_pert-Yv_eq)
         !       dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
         !       ! Solve 2x2 system: J*[dp; dYv]=-[F1; F2]
         !       detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
         !       if (abs(detJ).lt.1.0e-30_WP) exit
         !       dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
         !       dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
         !       ! Damped Newton-Raphson update
         !       pOld=p_eq
         !       YvOld=Yv_eq
         !       alpha=1.0_WP
         !       ! if ((abs(F1).lt.F_line_search_tol).and.(abs(F2/Eps0).lt.F_line_search_tol)) then
         !       if ((abs(F1).lt.F_line_search_tol).and.(abs(F2).lt.F_line_search_tol)) then
         !          p_eq=pOld+dp_nr
         !          Yv_eq=YvOld+dYv_nr
         !       else
         !          accepted=.false.
         !          do while (alpha.gt.1.0e-8_WP)
         !             p_try=p_eq+alpha*dp_nr
         !             Yv_try=Yv_eq+alpha*dYv_nr
         !             ! Keep the trial state inside the physical/log-safe domain
         !             if (p_try.le.lnP_eps) then
         !                alpha=0.5_WP*alpha
         !                cycle
         !             end if
         !             ! if ((Yv_try.le.Yvmin+fd_eps).or.(Yv_try.ge.Yvmax-fd_eps)) then
         !             !    alpha=0.5_WP*alpha
         !             !    cycle
         !             ! end if
         !             if ((rho0*(1.0_WP-Yv_try)-rhoA0).le.0.0_WP) then
         !                alpha=0.5_WP*alpha
         !                cycle
         !             end if
         !             T_try=get_T_lvg(p_try,Yv_try)
         !             if (T_try.le.0.0_WP) then
         !                alpha=0.5_WP*alpha
         !                cycle
         !             end if
         !             xv_try=get_xv(Yv_try)
         !             ppv_try=xv_try*p_try
         !             if (.not.check_p(ppv_try)) then
         !                alpha=0.5_WP*alpha
         !                cycle
         !             end if
         !             F1_try=PTsat(p_try,ppv_try,T_try)
         !             ! F2_try=Eps_res_lvg(p_try,T_try,Yv_try)
         !             ! res_try=sqrt(F1_try**2+(F2_try/Eps0)**2)
         !             F2_try=Eps_res_lvg(p_try,T_try,Yv_try)/Eps0
         !             res_try=sqrt(F1_try**2+F2_try**2)
         !             if (res_try.lt.res0) then
         !                p_eq=p_try
         !                Yv_eq=Yv_try
         !                T_eq=T_try
         !                accepted=.true.
         !                exit
         !             end if
         !             alpha=0.5_WP*alpha
         !          end do
         !          if (.not.accepted) then
         !             ! print*,"****************** Line search failed at it=",it," with alpha=",alpha
         !             exit
         !          end if
         !       end if
         !       ! Evaluate errors
         !       p_err=abs((p_eq-pOld)/pOld)
         !       Yv_err=abs((Yv_eq-YvOld)/(YvOld+1.0e-30_WP))
         !       ! Refresh xv and ppv from the accepted (p_eq, Yv_eq) before re-evaluating
         !       ! the residuals.  The previous ppv was stale from before the line search.
         !       xv=get_xv(Yv_eq)
         !       ppv=xv*p_eq
         !       F1=PTsat(p_eq,ppv,T_eq)
         !       ! F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)
         !       F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)/Eps0
         !       ! Per-iteration diagnostic (fires on EVERY iteration, including the
         !       ! one that converges, so we can always see the full trace).
         !       ! write(*,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') &
         !       !    '  it=',it,                                  &
         !       !    '  p=',p_eq,                                 &
         !       !    '  Yv=',Yv_eq,                               &
         !       !    '  T=',T_eq,                                 &
         !       !    '  alpha=',alpha,                            &
         !       !    '  |F1|=',abs(F1),                           &
         !       !    '  |F2/E0|=',abs(F2/Eps0),                   &
         !       !    '  perr=',p_err,                             &
         !       !    '  Yverr=',Yv_err
         !       ! if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol).and.(abs(F1).lt.F1_tol).and.(abs(F2/Eps0).lt.F2_tol)) then
         !       if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol).and.(abs(F1).lt.F1_tol).and.(abs(F2).lt.F2_tol)) then
         !          converged=.true.
         !          exit
         !       end if
         !    end do
         !    print*,'NR iterations=',it
         !    if (.not.converged) then
         !       ! print*,"****************** p-Yv iterations blew up. Skipping the cell!!"
         !       return
         !    end if
         !    ! Update equilibrium temperature
         !    T_eq=get_T_lvg(p_eq,Yv_eq)
         ! end subroutine solve_PTg_relax_lvg
         ! Gemini:
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
            real(WP) :: Yv_max_phys
            real(WP), parameter :: fd_eps=1.0e-8_WP,F_line_search_tol=0.1_WP
            integer  :: it,lsit
            logical  :: accepted
            ! Calculate the absolute physical ceiling for Yv based on available liquid
            Yv_max_phys=1.0_WP-(rhoA0/rho0)
            ! Iteratively solve for the equilibrium pressure and vapor mass fraction
            converged=.false.
            p_err=10.0_WP*p_tol
            Yv_err=10.0_WP*Yv_tol
            do it=1,NR_itmax
               ! Evaluate residuals at current state
               T_eq=get_T_lvg(p_eq,Yv_eq)
               xv=get_xv(Yv_eq)
               ppv=xv*p_eq
               if (.not.check_p(ppv)) return
               ! FULLY SCALED RESIDUALS (Makes Jacobian O(1))
               F1=PTsat(p_eq,ppv,T_eq)/p_eq
               F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)/Eps0
               res0=sqrt(F1**2+F2**2)
               ! --- Perturbation in p ---
               p_pert=p_eq*(1.0_WP+fd_eps)
               ppv_pert=xv*p_pert
               if (.not.check_p(ppv_pert)) return
               T_pert=get_T_lvg(p_pert,Yv_eq)
               F1p=PTsat(p_pert,ppv_pert,T_pert)/p_eq
               F2p=Eps_res_lvg(p_pert,T_pert,Yv_eq)/Eps0
               dF1dp=(F1p-F1)/(p_pert-p_eq)
               dF2dp=(F2p-F2)/(p_pert-p_eq)
               ! --- Perturbation in Yv ---
               Yv_pert=Yv_eq+fd_eps
               if (Yv_pert.gt.Yvmax-fd_eps) Yv_pert=Yv_eq-fd_eps
               if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
               xv_pert=get_xv(Yv_pert)
               ppv_pert=xv_pert*p_eq
               if (.not.check_p(ppv_pert)) return
               T_pert=get_T_lvg(p_eq,Yv_pert)
               F1Y=PTsat(p_eq,ppv_pert,T_pert)/p_eq
               F2Y=Eps_res_lvg(p_eq,T_pert,Yv_pert)/Eps0
               dF1dYv=(F1Y-F1)/(Yv_pert-Yv_eq)
               dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
               ! Solve 2x2 system: J*[dp; dYv]=-[F1; F2]
               detJ=dF1dp*dF2dYv-dF1dYv*dF2dp
               if (abs(detJ).lt.1.0e-30_WP) exit
               dp_nr =-(dF2dYv*F1-dF1dYv*F2)/detJ
               dYv_nr=-(dF1dp *F2-dF2dp *F1)/detJ
               ! --- DIRECTION-PRESERVING STEP LIMITER ---
               block
                  real(WP) :: max_scale
                  max_scale=1.0_WP
                  ! 1. Prevent pressure from changing by more than 50% in a single step
                  if (abs(dp_nr) .gt. 0.5_WP*p_eq) then
                     max_scale=min(max_scale, 0.5_WP*p_eq/abs(dp_nr))
                  end if
                  ! 2. Prevent Yv from crossing physical boundaries
                  if (dYv_nr .gt. 0.0_WP) then
                     if (Yv_eq + dYv_nr .ge. Yv_max_phys) then
                        max_scale=min(max_scale, 0.9_WP*(Yv_max_phys-Yv_eq)/dYv_nr)
                     end if
                     if (Yv_eq + dYv_nr .ge. Yvmax-fd_eps) then
                        max_scale=min(max_scale, 0.9_WP*(Yvmax-fd_eps-Yv_eq)/dYv_nr)
                     end if
                  else if (dYv_nr .lt. 0.0_WP) then
                     if (Yv_eq + dYv_nr .le. Yvmin + fd_eps) then
                        max_scale=min(max_scale, 0.9_WP*(Yv_eq-Yvmin-fd_eps)/abs(dYv_nr))
                     end if
                  end if
                  ! Scale the Newton step uniformly to keep pointing directly at the root
                  dp_nr =dp_nr *max_scale
                  dYv_nr=dYv_nr*max_scale
               end block
               ! -----------------------------------------
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
                     if (T_try.le.0.0_WP) then
                        alpha=0.5_WP*alpha
                        cycle
                     end if
                     xv_try=get_xv(Yv_try)
                     ppv_try=xv_try*p_try
                     if (.not.check_p(ppv_try)) then
                        alpha=0.5_WP*alpha
                        cycle
                     end if
                     F1_try=PTsat(p_try,ppv_try,T_try)/p_try
                     F2_try=Eps_res_lvg(p_try,T_try,Yv_try)/Eps0
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
                  if (.not.accepted) exit
               end if
               p_err=abs((p_eq-pOld)/pOld)
               Yv_err=abs((Yv_eq-YvOld)/(YvOld+1.0e-30_WP))
               ! print*,'p_err = ',p_err
               ! print*,'Yv_err = ',Yv_err
               xv=get_xv(Yv_eq)
               ppv=xv*p_eq
               T_eq=get_T_lvg(p_eq,Yv_eq)
               ! Evaluate the errors
               F1=PTsat(p_eq,ppv,T_eq)/p_eq
               F2=Eps_res_lvg(p_eq,T_eq,Yv_eq)/Eps0
               if ((p_err.lt.p_tol).and.(Yv_err.lt.Yv_tol).and.(abs(F1).lt.F1_tol).and.(abs(F2).lt.F2_tol)) then
                  converged=.true.
                  exit
               end if
            end do
            ! print*,'NR iterations=',it
            ! print*,'Line sreach iterations=',lsit
            ! print*,'alpha=',alpha
            if (.not.converged) return
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
            get_T_lvg=(1.0_WP-Yv_eq-bL*(rho0*(1.0_WP-Yv_eq)-rhoA0))/((rho0*(1.0_WP-Yv_eq)-rhoA0)*(GammaL-1.0_WP)*CvL/(p_eq+PinfL)+rhoA0*((GammaV-1.0_WP)*CvV*Yv_eq+(GammaA-1.0_WP)*CvA*(1.0_WP-Yv_eq))/p_eq)
         end function get_T_lvg
         !> Function that defines p-T saturation curve
         function PTsat(p_l,p_v,T_eq)
            real(WP), intent(in) :: p_l,p_v,T_eq
            real(WP) :: PTsat
            PTsat=AS+(BS+ES*p_l)/T_eq+CS*log(T_eq)+DS*log(p_l+PinfL)-log(p_v)
         end function PTsat
         !> Temperature derivative of p-T saturation curve function
         function dPTsatdT(p_l,T_eq)
            real(WP), intent(in) :: p_l,T_eq
            real(WP) :: dPTsatdT
            dPTsatdT=-(BS+ES*p_l)/T_eq**2+CS/T_eq
         end function dPTsatdT
         !> Pressure derivative of p-T saturation curve function for pure vapor p iteration
         function dPTsatdp_lv(p_eq,T_eq,dTdp)
            real(WP), intent(in) :: p_eq,T_eq,dTdp
            real(WP) :: dPTsatdp_lv
            dPTsatdp_lv=dPTsatdT(p_eq,T_eq)*dTdp+ES/T_eq+DS/(p_eq+PinfL)-1.0_WP/p_eq
         end function dPTsatdp_lv
         !> Residual of the internal energy conservation equation for the lvg case
         function Eps_res_lvg(p_eq,T_eq,Yv_eq)
            real(WP), intent(in) :: p_eq,T_eq,Yv_eq
            real(WP) :: Eps_res_lvg
            Eps_res_lvg=(rho0*(1.0_WP-Yv_eq)-rhoA0)*get_IL_PT(p_eq,T_eq)+rhoA0*get_IG_PT(P_eq,T_eq,Yv_eq)-Eps0*(1.0_WP-Yv_eq)
         end function Eps_res_lvg
         !> Residual of the internal energy conservation equation for the lv case
         function Eps_res_lv(p_eq,T_eq)
            real(WP), intent(in) :: p_eq,T_eq
            real(WP) :: Eps_res_lv,rho_l,rho_g
            rho_l=get_RHOL(P=p_eq,T=T_eq)
            rho_g=get_RHOG(P=p_eq,T=T_eq,Yv=1.0_WP)
            Eps_res_lv=(rho0-rho_g)/(rho_l-rho_g)*get_IL_PT(P=p_eq,T=T_eq)+(rho_l-rho0)/(rho_l-rho_g)*get_IG_PT(P=p_eq,T=T_eq,Yv=1.0_WP)-Eps0
         end function Eps_res_lv
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
         !> (SG) Subroutine that updates the coefficients of the quadradic equilibrium temperature equation as functions of equilibrium pressure
         ! subroutine get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
         !    real(WP), intent(in)  :: p_eq
         !    real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         !    ! Coefficients
         !    ap=sum(Q(1:2))*CvL*CvV*((GammaV-1.0_WP)*(p_eq+GammaL*PinfL)-(GammaL-1.0_WP)*p_eq)
         !    bp=sum(Q(3:4))*((GammaL-1.0_WP)*CvL*p_eq-(GammaV-1.0_WP)*CvV*(p_eq+PinfL))+&
         !    &  sum(Q(1:2))*((GammaV-1.0_WP)*CvV*qL*(p_eq+PinfL)-(GammaL-1.0_WP)*CvL*qV*p_eq)+&
         !    &  CvV*(p_eq+PinfL)*p_eq-CvL*p_eq*(p_eq+GammaL*PinfL)
         !    dp=(qV-qL)*(p_eq+PinfL)*p_eq
         !    ! Pressure derivative of the coefficients
         !    dapdp=sum(Q(1:2))*CvL*CvV*(GammaV-GammaL)
         !    dbpdp=sum(Q(3:4))*((GammaL-1.0_WP)*CvL-(GammaV-1.0_WP)*CvV)+&
         !    &     sum(Q(1:2))*((GammaV-1.0_WP)*CvV*qL-(GammaL-1.0_WP)*CvL*qV)+&
         !    &     CvV*(2.0_WP*p_eq+PinfL)-CvL*(2.0_WP*p_eq+GammaL*PinfL)
         !    ddpdp=(qV-qL)*(2.0_WP*p_eq+PinfL)
         ! end subroutine get_coeffs_lv
         !> (NASG) Subroutine that updates the coefficients of the quadradic equilibrium temperature equation as functions of equilibrium pressure
         subroutine get_coeffs_lv(p_eq,ap,bp,dp,dapdp,dbpdp,ddpdp)
            real(WP), intent(in)  :: p_eq
            real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
            ! Coefficients
            ap=sum(Q(1:2))*CvL*CvV*((GammaV-GammaL)*p_eq+GammaL*(GammaV-1.0_WP)*PinfL)
            bp=(CvV*(1.0_WP-sum(Q(1:2))*bL)-CvL)*p_eq**2+&
            &  (PinfL*(CvV*(1.0_WP-sum(Q(1:2))*bL)-GammaL*CvL)+&
            &  sum(Q(1:2))*((GammaV-1.0_WP)*CvV*qL-(GammaL-1.0_WP)*CvL*qV)+&
            &  sum(Q(3:4))*((GammaL-1.0_WP)*CvL-(GammaV-1.0_WP)*CvV))*p_eq+&
            &  (GammaV-1.0_WP)*CvV*PinfL*(sum(Q(1:2))*qL-sum(Q(3:4)))
            dp=p_eq*(p_eq+PinfL)*(qV*(1.0_WP-sum(Q(1:2))*bL)-qL+bL*sum(Q(3:4)))
            ! Pressure derivative of the coefficients
            dapdp=sum(Q(1:2))*CvL*CvV*(GammaV-GammaL)
            dbpdp=2.0_WP*(CvV*(1.0_WP-sum(Q(1:2))*bL)-CvL)*p_eq+&
            &     PinfL*(CvV*(1.0_WP-sum(Q(1:2))*bL)-GammaL*CvL)+&
            &     sum(Q(1:2))*((GammaV-1.0_WP)*CvV*qL-(GammaL-1.0_WP)*CvL*qV)+&
            &     sum(Q(3:4))*((GammaL-1.0_WP)*CvL-(GammaV-1.0_WP)*CvV)
            ddpdp=(2.0_WP*p_eq+PinfL)*(qV*(1.0_WP-sum(Q(1:2))*bL)-qL+bL*sum(Q(3:4)))
         end subroutine get_coeffs_lv
         !> Sanity check pressure value
         logical function check_p(p_v)
            real(WP), intent(in) :: p_v
            check_p=(p_v.gt.lnP_eps).and.((p_v+PinfL).gt.lnP_eps)
         end function check_p
   end subroutine PTg_relax

   !> Compute viscosity: constant gas and liquid, VF-weighted blend
   !> Contains commented-out Sutherland law for variable gas viscosity (dimensional form)
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pYv,pVisc,pBeta,pCond,pDiff
      real(WP) :: mu_g,mu_l,cp_g
      real(WP), parameter :: myeps=1.0e-15_WP
      !> Sutherland's law parameters (dimensional, SI units)
      !> mu_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K
      ! real(WP), parameter :: mu_ref=1.716e-5_WP   !< Reference viscosity [Pa·s]
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
            pYv=>fs%Yv%mf(lvl)%dataptr(mfi)
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
               ! Thermal conductivity (zero when Pr=0 i.e. no heat conduction)
               cp_g=pYv(i,j,k,1)*CpV+(1.0_WP-pYv(i,j,k,1))*CpA
               if (PrL.gt.0.0_WP) then
                  pCond(i,j,k,1)=mu_l*CpL/PrL
               else
                  pCond(i,j,k,1)=0.0_WP
               end if
               if (PrG.gt.0.0_WP) then
                  pCond(i,j,k,2)=mu_g*cp_g/PrG
               else
                  pCond(i,j,k,2)=0.0_WP
               end if
               ! Vapor mass diffusivity (zero when Sc=0 i.e. no mass diffusion)
               if (ScV.gt.0.0_WP) then
                  pDiff(i,j,k,1)=mu_g/ScV
               else
                  pDiff(i,j,k,1)=0.0_WP
               end if
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
         ! Get liquid EOS parameters
         call param_read('Liquid specific heat capacity at constant volume',CvL)
         call param_read('Liquid specific heat capacity ratio',GammaL)
         call param_read('Liquid reference energy shift',qL)
         call param_read('Liquid reference entropy shift',qpL)
         call param_read('Liquid stiffening pressure',PinfL)
         call param_read('Liquid co-volume',bL)
         ! Get vapor EOS parameters
         call param_read('Vapor specific heat capacity at constant volume',CvV)
         call param_read('Vapor specific heat capacity ratio',GammaV)
         call param_read('Vapor reference energy shift',qV)
         call param_read('Vapor reference entropy shift',qpV)
         ! Get air EOS parameters
         call param_read('Air specific heat capacity at constant volume',CvA)
         call param_read('Air specific heat capacity ratio',GammaA)
         call param_read('Air reference energy shift',qA)
         call param_read('Air reference entropy shift',qpA)
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
         ! Prandtl and Schmidt numbers
         call param_read('Liquid Prandtl number',PrL)
         call param_read('Gas Prandtl number',PrG)
         call param_read('Vapor Schmidt number',ScV)
         ! Cylinder geometry
         call param_read('Cylinder diameter',dcyl)
         call param_read('Cylinder location',xcyl)
         ! Domain dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         ! Saturation curve coefficients
         AS=(CpL-CpV+qpV-qpL)/(CpV-CvV)
         BS=(qL-qV)/(CpV-CvV)
         CS=(CpV-CpL)/(CpV-CvV)
         DS=(CpL-CvL)/(CpV-CvV)
         ES=bL/(CpV-CvV)
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
         write(message,'("[Gas]    GammaG=",es12.5)') GammaG; call log(message)
         write(message,'("[Visc]   muG=",es12.5," muL=",es12.5)') muG,muL; call log(message)
         write(message,'("[Cylinder] diameter=",es12.5," location=",es12.5)') dcyl,xcyl; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='Sembian_blastwave_NASG'
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
         call param_read('Initial dt',dt_init)
         call param_read('Max CFL',time%cflmax)
         time%dt=dt_init
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap, amrex_bc_reflect_even, amrex_bc_reflect_odd
         use amrmpcomp_class, only: BC_REFLECT
         ! Create flow solver
         call fs%initialize(amr=amr,name='Sembian_blastwave_NASG')
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide relaxation model
         fs%relax=>PTg_relax
         ! Set initial conditions via blastwave callback
         fs%user_init=>blastwave_init
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
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=8,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
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
            ! Initialize Yv
            init_Yv_restart: block
               use amrex_amr_module, only: amrex_mfiter,amrex_mfiter_build,amrex_mfiter_destroy,amrex_box
               type(amrex_mfiter) :: mfi
               type(amrex_box) :: bx
               real(WP), contiguous, pointer, dimension(:,:,:,:) :: pQ
               real(WP) :: Yv0
               integer :: lvl,i,j,k
               call param_read('Initial vapor mass fraction',Yv0)
               do lvl=0,amr%maxlvl
                  call amrex_mfiter_build(mfi,fs%Q%mf(lvl),tiling=.false.)
                  do while (mfi%next())
                     bx=mfi%tilebox()
                     pQ=>fs%Q%mf(lvl)%dataptr(mfi)
                     do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
                        pQ(i,j,k,8)=pQ(i,j,k,2)*Yv0
                     end do; end do; end do
                  end do
                  call amrex_mfiter_destroy(mfi)
               end do
               call fs%Q%fill(time=time%t)
               call fs%get_primitive(fs%Q)
            end block init_Yv_restart
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%UVW,fs%UVW,fs%UVW,compX=1,compY=2,compZ=3)
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
         call viz%initialize(amr,'Sembian_blastwave_NASG',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_scalar(fs%Yv,1,'Yv')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         call viz%add_scalar(fs%IL,1,'IL')
         call viz%add_scalar(fs%IG,1,'IG')
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
            time%dt=dt_init
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
         call Umag%get_magnitude(fs%UVW,fs%UVW,fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         ! if (fs%amr%amRoot) print*, 'time: ', time%t, ' dt: ', time%dt,' occurs: ',viz_evt%occurs()
         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time%t)
         end if
         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/Sembian_blastwave_NASG_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

      ! Save the final checkpoint
      save_final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/Sembian_blastwave_NASG_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block save_final_checkpoint

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
