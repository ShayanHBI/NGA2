!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,           only: WP
   use geometry,            only: cfg,Lx
   use ddadi_class,         only: ddadi
   ! use fft3d_class,         only: fft3d
   use hypre_str_class,     only: hypre_str
   use lowmach_class,       only: lowmach
   use finitechem_class,    only: finitechem
   use timetracker_class,   only: timetracker
   use ensight_class,       only: ensight
   use event_class,         only: event
   use monitor_class,       only: monitor
   use fcmech
   implicit none
   private

   !> Single low Mach flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: ps
   ! type(fft3d),       public :: ps
   type(ddadi),       public :: vs,ss
   type(lowmach),     public :: fs
   type(finitechem),  public :: fc
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile,fcfile

   !> Simulation subroutines
   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW,resRHO
   real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:),   allocatable :: resSC,SCtmp

   !> Flame definition
   logical  :: use_reactions
   real(WP) :: xFlame,Tu,Tb,Uin
   real(WP) :: Schmidt,Prandtl
   real(WP), dimension(:), allocatable :: Yu,Yb

contains


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator

   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator

   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator

   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator

   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator

   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator

   !> Function that localizes the x- boundary for SC
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_exists
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate (resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (resRHO(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resRHO=0.0_WP
         allocate (Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate (resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,nspec+1)); resSC=0.0_WP
         allocate (SCtmp (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,nspec+1)); SCtmp=0.0_WP
      end block allocate_work_arrays

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub iterations',time%itmax)
         time%dt=time%dtmax
      end block initialize_timetracker

      ! Create a scalar solver
      create_fc: block
         use multivdscalar_class, only: bcond,dirichlet,neumann,bquick
         use string,              only: str_medium
         type(bcond),pointer :: mybc
         integer :: nsc,i,j,k,n
         character(len=str_medium), dimension(:), allocatable :: spname
         ! Create finite chem object
         fc=finitechem(cfg=cfg,scheme=bquick,name='fc')
         fc%use_scheduler=.false.
         ! Define boundary conditions
         call fc%add_bcond(name='inflow', type=dirichlet,locator=xm_locator_sc)
         call fc%add_bcond(name='outflow',type=neumann,  locator=xp_locator,dir='+x')
         ! Assign constant diffusivity
         ! call param_read('Dynamic diffusivity',diffusivity)
         ! fc%diff=diffusivity
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=13)
         ! Setup the solver
         call fc%setup(implicit_solver=ss)
         ! Allocate memory
         allocate (spname(nspec))
         allocate (Yu(nspec)); Yu=0.0_WP
         allocate (Yb(nspec)); Yb=0.0_WP
         ! Read in info
         call param_read('Use reactions',use_reactions)
         call param_read('Flame location',xFlame)
         call param_read('Unburned temperature',Tu)
         call param_read('Burned temperature',Tb)
         call param_read('Pressure',fc%Pthermo)
         ! Get species names
         call fcmech_get_speciesnames(spname)
         ! Unburned composition
         do nsc=1,nspec
            if (param_exists('Unburned '//trim(spname(nsc)))) then
               call param_read('Unburned '//trim(spname(nsc)),Yu(nsc))
            end if
         end do
         Yu=Yu/sum(Yu)
         ! Burned composition
         do nsc=1,nspec
            if (param_exists('Burned '//trim(spname(nsc)))) then
               call param_read('Burned '//trim(spname(nsc)),Yb(nsc))
            end if
         end do
         Yb=Yb/sum(Yb)
         ! Initialize the flame
         do nsc=1,nspec
            do k=fc%cfg%kmin_,fc%cfg%kmax_
               do j=fc%cfg%jmin_,fc%cfg%jmax_
                  do i=fc%cfg%imin_,fc%cfg%imax_
                     fc%SC(i,j,k,nsc)=Yu(nsc)+(Yb(nsc)-Yu(nsc))*0.5_WP*(1.0_WP+tanh((fc%cfg%xm(i)-xFlame)/(Lx/20.0_WP)))
                  end do
               end do
            end do
         end do
         do k=fc%cfg%kmin_,fc%cfg%kmax_
            do j=fc%cfg%jmin_,fc%cfg%jmax_
               do i=fc%cfg%imin_,fc%cfg%imax_
                  fc%SC(i,j,k,nspec+1)=Tu+(Tb-Tu)*0.5_WP*(1.0_WP+tanh((fc%cfg%xm(i)-xFlame)/(Lx/20.0_WP)))
               end do
            end do
         end do
         ! Sync scalars
         do nsc=1,nspec+1
            call fc%cfg%sync(fc%SC(:,:,:,nsc))
         end do
         ! Apply boundary conditions
         call fc%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fc%SC(i,j,k,1:nspec)=2.0_WP*Yu-fc%SC(i+1,j,k,1:nspec)
            fc%SC(i,j,k,nspec+1)=2.0_WP*Tu-fc%SC(i+1,j,k,nspec+1)
         end do
         call fc%apply_bcond(time%t,time%dt)
         ! Get fluid properties
         call fc%get_density()
         call fc%get_viscosity()
         call fc%get_diffusivity()
         ! Get monitoring quantities
         call fc%get_max()
      end block create_fc

      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use lowmach_class,   only: bcond,dirichlet,clipped_neumann
         use hypre_str_class, only: pcg_pfmg2
         integer :: n,i,j,k
         type(bcond),pointer :: mybc
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define boundary conditions
         call fs%add_bcond(name='inflow', type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.True., locator=xp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=6
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Configure implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initialize velocity field
         call param_read('Flame speed',Uin)
         fs%U=Uin
         fs%V=0.0_WP
         fs%W=0.0_WP
         ! Set density from scalar
         fs%rho=fc%rho
         fs%rhoold=fs%rho
         fs%visc=fc%visc
         ! Form momentum
         call fs%rho_multiply()
         ! Apply boundary conditions
         call fs%get_bcond('inflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=Uin
         end do
         call fs%apply_bcond(time%t,time%dt)
         ! Get cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(drhodt=resRHO)
         ! Compute MFR through all boundaries
         call fs%get_mfr()
      end block create_velocity_solver

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='flame1D')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('density',fs%rho)
         call ens_out%add_scalar('viscosity',fs%visc)
         call ens_out%add_scalar('thermal_diff',fc%diff(:,:,:,nspec+1))
         ! call ens_out%add_scalar('YCH4',fc%SC(:,:,:,sCH4))
         call ens_out%add_scalar('YOH',fc%SC(:,:,:,sOH))
         call ens_out%add_scalar('YO2',fc%SC(:,:,:,sO2))
         call ens_out%add_scalar('YN2',fc%SC(:,:,:,sN2))
         call ens_out%add_scalar('YCO2',fc%SC(:,:,:,sCO2))
         call ens_out%add_scalar('YH2O',fc%SC(:,:,:,sH2O))
         call ens_out%add_scalar('YCO',fc%SC(:,:,:,sCO))
         call ens_out%add_scalar('YNC12H26',fc%SC(:,:,:,sXC12H26))
         call ens_out%add_scalar('YHMN',fc%SC(:,:,:,sHMN))
         call ens_out%add_scalar('resRHO',resRHO)
         call ens_out%add_scalar('T',fc%SC(:,:,:,nspec+1))
         ! call ens_out%add_scalar('SRC_T',fc%SRCchem(:,:,:,nspec+1))
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call fc%get_max()
         call fc%get_int()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create FC monitor
         fcfile=monitor(fs%cfg%amRoot,'fc')
         call fcfile%add_column(time%n,'Timestep number')
         call fcfile%add_column(time%t,'Time')
         call fcfile%add_column(fs%CFLc_x,'Min Temperature')
         call fcfile%add_column(fs%CFLc_y,'Max Temperature')
         call fcfile%add_column(fs%CFLc_z,'Min sumY')
         call fcfile%add_column(fs%CFLv_x,'Max sumY')
         call fcfile%add_column(fc%rhomax,'RHOmax')
         call fcfile%add_column(fc%rhomin,'RHOmin')
         call fcfile%add_column(fc%Pthermo,'Pthermo')
         call fcfile%add_column(fc%SCmin(nspec+1),'Min Temperature')
         call fcfile%add_column(fc%SCmax(nspec+1),'Max Temperature')
         call fcfile%write()
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         ! call consfile%add_column(fc%SCint,'fc integral')
         call consfile%add_column(fc%rhoint,'RHO integral')
         ! call consfile%add_column(fc%rhoSCint,'rhoSC integral')
         call consfile%write()
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Perform time integration
      do while (.not. time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old scalar
         fc%rhoold=fc%rho
         fc%SCold=fc%SC

         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         fc%SRCchem=0.0_WP
         if (use_reactions) then
            call fc%react(time%dt)
         end if

         ! call fc%diffusive_source(time%dt)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! ============ SCALAR SOLVER ======================
            scalar_solver: block
               use messager, only: die
               integer :: nsc
               integer :: i,j,k
               logical, dimension(:,:,:,:), allocatable :: flag

               allocate (flag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,fc%nscalar))
               call fc%metric_reset()

               fc%SRC=0.0_WP
               ! call fc%pressure_source()

               ! Build mid-time scalar
               fc%SC=0.5_WP*(fc%SC+fc%SCold)

               call fc%diffusive_source(time%dt)
               ! Explicit calculation of drhoSC/dt from scalar equation
               call fc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
               ! Assemble explicit residual
               do nsc=1,fc%nscalar
                  resSC(:,:,:,nsc)=time%dt*resSC(:,:,:,nsc)-2.0_WP*fc%rho*fc%SC(:,:,:,nsc)+(fc%rho+fc%rhoold)*fc%SCold(:,:,:,nsc)+fc%rho*fc%SRCchem(:,:,:,nsc)+fc%SRC(:,:,:,nsc)
                  SCtmp(:,:,:,nsc)=2.0_WP*fc%SC(:,:,:,nsc)-fc%SCold(:,:,:,nsc)+resSC(:,:,:,nsc)/fc%rho
               end do

               ! Apply it to get explicit scalar prediction
               do nsc=1,fc%nscalar
                  do k=fc%cfg%kmino_,fc%cfg%kmaxo_
                     do j=fc%cfg%jmino_,fc%cfg%jmaxo_
                        do i=fc%cfg%imino_,fc%cfg%imaxo_
                           if (nsc.eq.nspec+1) then
                              if (SCtmp(i,j,k,nsc).le.290.0_WP.or.SCtmp(i,j,k,nsc).ge.4000.0_WP) then
                                 flag(i,j,k,nsc)=.true.
                              else
                                 flag(i,j,k,nsc)=.false.
                              end if
                           else
                              if (SCtmp(i,j,k,nsc).le.0.0_WP.or.SCtmp(i,j,k,nsc).ge.1.0_WP) then
                                 flag(i,j,k,nsc)=.true.
                              else
                                 flag(i,j,k,nsc)=.false.
                              end if
                           end if
                        end do
                     end do
                  end do
               end do
               ! Adjust metrics
               call fc%metric_adjust(SCtmp,flag)
               ! Clean up
               deallocate(flag)

               ! Recompute drhoSC/dt
               call fc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)
               ! resSC=-2.0_WP*(fc%SC-fc%SCold)+time%dt*resSC
               ! Assemble explicit residual
               do nsc=1,fc%nscalar
                  resSC(:,:,:,nsc)=time%dt*resSC(:,:,:,nsc)-2.0_WP*fc%rho*fc%SC(:,:,:,nsc)+(fc%rho+fc%rhoold)*fc%SCold(:,:,:,nsc)+fc%rho*fc%SRCchem(:,:,:,nsc)+fc%SRC(:,:,:,nsc)
               end do
               !    resSC(:,:,:,nsc)=time%dt*resSC(:,:,:,nsc)-2.0_WP*fc%rho*fc%SC(:,:,:,nsc)+(fc%rho+fc%rhoold)*fc%SCold(:,:,:,nsc)+fc%rho*fc%SRCchem(:,:,:,nsc)
               ! Form implicit residual
               call fc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)

               ! Advance scalar field
               fc%SC=2.0_WP*fc%SC-fc%SCold+resSC

               ! Apply boundary conditions
               call fc%apply_bcond(time%t,time%dt)
               ! Re-apply Dirichlet BCs
               dirichlet_scalar: block
                  use multivdscalar_class, only: bcond
                  type(bcond),pointer :: mybc
                  integer :: n
                  call fc%get_bcond('inflow',mybc)
                  do n=1,mybc%itr%no_
                     i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                     fc%SC(i,j,k,1:nspec)=2.0_WP*Yu-fc%SC(i+1,j,k,1:nspec)
                     fc%SC(i,j,k,nspec+1)=2.0_WP*Tu-fc%SC(i+1,j,k,nspec+1)
                  end do
               end block dirichlet_scalar

            end block scalar_solver

            ! ============ UPDATE PROPERTIES ====================
            call fc%get_density()
            ! call fc%rescale_density()
            call fc%get_viscosity()
            call fc%get_diffusivity()
            ! call fc%update_pressure()
            fs%visc=fc%visc

            ! ============ VELOCITY SOLVER ======================

            ! Build n+1 density
            fs%rho=0.5_WP*(fc%rho+fc%rhoold)

            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond),pointer :: mybc
               integer :: n,i,j,k
               call fs%get_bcond('inflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=Uin
               end do
            end block dirichlet_velocity

            call fs%rho_multiply()

            ! Solve Poisson equation
            call fc%get_drhodt(dt=time%dt,drhodt=resRHO)
            call fs%correct_mfr(drhodt=resRHO)
            call fs%get_div(drhodt=resRHO)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide
            ! ===================================================

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fc%get_drhodt(dt=time%dt,drhodt=resRHO)
         call fs%get_div(drhodt=resRHO)

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Perform and output monitoring
         call fs%get_max()
         call fc%get_max()
         call mfile%write()
         call cflfile%write()
         call consfile%write()

      end do

      ! Output profiles
      post_process: block
         integer :: i
         ! Open the file for writing
         open(unit=10, file=trim('NGA2.dat'), status='replace', action='write')
         ! Write the arrays in three columns
         write(10,'(a15,a15,a15,a15,a15,a15,a15,a15)') 'x', 'XC12H26', 'HMN', 'N2', 'OH', 'CO', 'T', 'u'
         do i=cfg%imin_,cfg%imax_
            write(10,'(8E15.6)') cfg%xm(i),fc%SC(i,1,1,sXC12H26),fc%SC(i,1,1,sHMN),fc%SC(i,1,1,sN2),fc%SC(i,1,1,sOH),fc%SC(i,1,1,sCO),fc%SC(i,1,1,nspec+1),Ui(i,1,1)
         end do
         ! Close the file
         close(10)
      end block post_process

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none

      ! Get rid of all objects-need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate (resU,resV,resW,Ui,Vi,Wi,resSC,resRHO,SCtmp)

   end subroutine simulation_final


end module simulation
