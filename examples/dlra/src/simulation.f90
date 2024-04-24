!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use geometry,          only: cfg
   use ddadi_class,       only: ddadi
   use hypre_str_class,   only: hypre_str
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use sgsmodel_class,    only: sgsmodel
   use flamelet_class,    only: flamelet
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Single phase low Mach flow solver and scalar solver and corresponding time tracker
   type(hypre_str),   public :: vs,ps
   type(ddadi),       public :: ss
   type(lowmach),     public :: fs
   type(vdscalar),    public :: sc    
   type(timetracker), public :: time
   type(sgsmodel),    public :: sgs
   type(flamelet),    public :: flm

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,consfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:),     allocatable :: resU,resV,resW,resSC
   real(WP), dimension(:,:,:),     allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU
   real(WP), dimension(:,:,:,:),   allocatable :: SR

   !> Inlet conditions
   real(WP) :: Z_jet,Z_cof
   real(WP) :: D_jet,D_cof
   real(WP) :: U_jet,U_cof
   real(WP) :: rho_jet,rho_cof

   !> Integral of pressure residual
   real(WP) :: int_RP=0.0_WP

   !> Turbulence
   ! real(WP) :: SchmidtSGS

   !> Combustion
   integer  :: nfilter,ncells,n_Y,iY
   logical  :: rho_limiter
   real(WP) :: rho_min,rho_max
   real(WP), dimension(:,:,:),   allocatable :: drho,ZgradMagSq,T
   real(WP), dimension(:,:,:,:), allocatable :: Y
   character(len=str_medium), dimension(:), allocatable :: Y_name

   !> Time stepping
   character(len=str_medium) :: time_stepping


contains


   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator


   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator


   !> Function that localizes jet at -x
   function jet(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Jet in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.le.0.5_WP*D_jet.and.i.eq.pg%imin) isIn=.true.
   end function jet


   !> Function that localizes coflow at -x
   function coflow(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Coflow in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.gt.0.5_WP*D_jet.and.radius.le.0.5_WP*D_cof.and.i.eq.pg%imin) isIn=.true.
   end function coflow


   !> Function that localizes jet at -x
   function jetsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Jet in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.le.0.5_WP*D_jet.and.i.eq.pg%imin-1) isIn=.true.
   end function jetsc


   !> Function that localizes coflow at -x
   function coflowsc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: radius
      logical :: isIn
      isIn=.false.
      ! Coflow in yz plane
      radius=norm2([pg%ym(j),pg%zm(k)]-[0.0_WP,0.0_WP])
      if (radius.gt.0.5_WP*D_jet.and.radius.le.0.5_WP*D_cof.and.i.eq.pg%imin-1) isIn=.true.
   end function coflowsc


   !> Spatially filter density changes
   subroutine filter_drho()
      use mpi_f08 , only: MPI_SUM,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      integer  :: i,j,k,ierr,n
      real(WP) :: vol,my_vol,volume,my_drho,my_drho2,drho_mean,drho2_mean,sdt
      real(WP), dimension(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_) :: Fdrho

      ! Find the deviation from mean
      my_vol  =0.0_WP
      my_drho =0.0_WP
      my_drho2=0.0_WP
      do k=sc%cfg%kmin_,sc%cfg%kmax_
         do j=sc%cfg%jmin_,sc%cfg%jmax_
            do i=sc%cfg%imin_,sc%cfg%imax_
               vol=sc%cfg%dx(i)*sc%cfg%dy(j)*sc%cfg%dz(k)*sc%cfg%VF(i,j,k)
               my_vol=my_vol+vol
               my_drho=my_drho+vol*drho(i,j,k)
               my_drho2=my_drho2+vol*drho(i,j,k)**2
            end do
         end do
      end do
      call MPI_ALLREDUCE(my_vol,volume,1,MPI_REAL_WP,MPI_SUM,sc%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_drho,drho_mean,1,MPI_REAL_WP,MPI_SUM,sc%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_drho2,drho2_mean,1,MPI_REAL_WP,MPI_SUM,sc%cfg%comm,ierr)
      drho_mean =drho_mean /volume
      drho2_mean=drho2_mean/volume
      sdt=sqrt(abs(drho2_mean-drho_mean**2))

      ! Perform the filtering
      do n=1,nfilter
         Fdrho=drho
         do k=sc%cfg%kmin_,sc%cfg%kmax_
            do j=sc%cfg%jmin_,sc%cfg%jmax_
               do i=sc%cfg%imin_,sc%cfg%imax_
                  if ((fs%mask(i,j,k).eq.0).and.(abs(drho(i,j,k)-drho_mean).gt.5.0_WP*sdt)) then
                     Fdrho(i,j,k)=sum(sgs%filtern(:,:,:,i,j,k)*drho(i-1:i+1,j-1:j+1,k-1:k+1))
                  end if
               end do
            end do
         end do
         drho(sc%cfg%imin_:sc%cfg%imax_,sc%cfg%jmin_:sc%cfg%jmax_,sc%cfg%kmin_:sc%cfg%kmax_)=Fdrho(sc%cfg%imin_:sc%cfg%imax_,sc%cfg%jmin_:sc%cfg%jmax_,sc%cfg%kmin_:sc%cfg%kmax_)
         ! Sync drho
         call sc%cfg%sync(drho)
      end do

   end subroutine filter_drho


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_getsize,param_read
      implicit none


      ! Read in inputs
      call param_read('Z jet',Z_jet)
      call param_read('D jet',D_jet)
      call param_read('U jet',U_jet)
      call param_read('Z coflow',Z_cof)
      call param_read('D coflow',D_cof)
      call param_read('U coflow',U_cof)
      call param_read('Filtering levels',nfilter)
      call param_read('Density limiter',rho_limiter)
      n_Y=param_getsize('Ensight output species')


      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use geometry,        only: latbc
         use hypre_str_class, only: pcg_pfmg2,gmres_pfmg
         use lowmach_class,   only: dirichlet,clipped_neumann,slip
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Define BCs
         call fs%add_bcond(name='jet'    ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=jet)
         call fs%add_bcond(name='coflow' ,type=dirichlet      ,face='x',dir=-1,canCorrect=.false.,locator=coflow)
         call fs%add_bcond(name='outflow',type=clipped_neumann,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         if (trim(latbc).eq.'slip') then
            call fs%add_bcond(name='yp',type=slip,face='y',dir=+1,canCorrect=.false.,locator=yp_locator)
            call fs%add_bcond(name='ym',type=slip,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
            call fs%add_bcond(name='zp',type=slip,face='z',dir=+1,canCorrect=.false.,locator=zp_locator)
            call fs%add_bcond(name='zm',type=slip,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)
         end if
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=18
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         vs=hypre_str(cfg=cfg,name='Velocity',method=gmres_pfmg,nst=7)
         call param_read('Implicit iteration',vs%maxit)
         call param_read('Implicit tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_velocity_solver


      ! Create a scalar solver
      create_scalar_solver: block
         use vdscalar_class, only: dirichlet,neumann,quick
         real(WP) :: diffusivity
         ! Create scalar solver
         sc=vdscalar(cfg=cfg,scheme=quick,name='MixFrac')
         ! Define BCs
         call sc%add_bcond(name='jet'   ,type=dirichlet,locator=jetsc   )
         call sc%add_bcond(name='coflow',type=dirichlet,locator=coflowsc)
         ! Outflow on the right
         call sc%add_bcond(name='outflow',type=neumann,locator=xp_locator,dir='+x')
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=13)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
      end block create_scalar_solver


      ! Create a combustion model
      create_combustion: block
         use flameletLib_class, only: sfm
         use tabulation,        only: tabulate_flamelet
         character(len=str_medium) :: chfname
         logical :: mkchmtbl
         ! Create the chemtable
         call param_read('Chemtable file name',chfname)
         call param_read('Create chemtable',mkchmtbl)
         if ((cfg%amRoot).and.(mkchmtbl)) call tabulate_flamelet(model=sfm,chfname=chfname)
         ! Construct the flamelet object
         flm=flamelet(cfg=cfg,flmModel=sfm,tablefile=trim(chfname),name='Steady flamelet model')
         call flm%print()
      end block create_combustion


      ! Create a SGS model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
         ! call param_read('SGS Schmidt number',SchmidtSGS)
      end block create_sgs


      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(resU(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
         ! Turbulence
         allocate(gradU(1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR(1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         ! Combustion
         allocate(drho      (sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))    ; drho=0.0_WP
         allocate(ZgradMagSq(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))    ; ZgradMagSq=0.0_WP
         allocate(T         (sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))    ; T=0.0_WP
         allocate(Y         (sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_,n_Y)); Y=0.0_WP
         allocate(Y_name(n_Y))
      end block allocate_work_arrays


      ! Initialize time tracker
      initialize_timetracker: block
         use messager, only: die
         time=timetracker(amRoot=fs%cfg%amRoot,name='dlra')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub iterations',time%itmax)
         call param_read('Time stepping',time_stepping)
         if ((trim(time_stepping).ne.'explicit').and.(trim(time_stepping).ne.'implicit')) call die('Time stepping must be either explicit or implicit.')
         time%dt=time%dtmax
      end block initialize_timetracker


      ! Initialize our mixture fraction field
      initialize_scalar: block
         use vdscalar_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         ! Zero initial field
         sc%SC=0.0_WP
         ! Apply BCs
         call sc%get_bcond('jet',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc%SC(i,j,k)=2.0_WP*Z_jet-sc%SC(i+1,j,k)
         end do
         call sc%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            sc%SC(i,j,k)=2.0_WP*Z_cof-sc%SC(i+1,j,k)
         end do
      end block initialize_scalar


      ! Initialize combustion
      initialize_combustion: block
         ! Number of cells
         ncells=cfg%nxo_*cfg%nyo_*cfg%nzo_
         ! Lookup viscosity
         call flm%chmtbl%lookup('viscosity',fs%visc,sc%SC,flm%Zvar,flm%chi,ncells)
         ! Lookup diffusivity
         call flm%chmtbl%lookup('diffusivity',sc%diff,sc%SC,flm%Zvar,flm%chi,ncells)
         ! Mixture fraction gradient
         call sc%get_gradient(itpr_x=fs%itpr_x,itpr_y=fs%itpr_y,itpr_z=fs%itpr_z,SCgradMagSq=ZgradMagSq)
         ! Mixture fraction variance
         call flm%get_Zvar(delta=sgs%delta,ZgradMagSq=ZgradMagSq,Z=sc%SC)
         ! Lookup density
         call flm%chmtbl%lookup('density',sc%rho,sc%SC,flm%Zvar,flm%chi,ncells)
         ! Find the min and max for rho
         call flm%chmtbl%lookup_max('density',rho_max)
         call flm%chmtbl%lookup_min('density',rho_min)
         ! Post-processed thermochemical quantities
         call flm%chmtbl%lookup('temperature',T,sc%SC,flm%Zvar,flm%chi,ncells)
         call param_read('Ensight output species',Y_name)
         do iY=1,n_Y
            Y_name(iY)=trim('Y_'//Y_name(iY))
            call flm%chmtbl%lookup(Y_name(iY),Y(:,:,:,iY),sc%SC,flm%Zvar,flm%chi,ncells)
         end do
      end block initialize_combustion


      ! Debug
      debug: block
         integer :: ftest,i,j,nZmean,nZvar
         real(WP) :: dZmean,dZvar
         real(WP), dimension(:,:), allocatable :: Zmean,Zvar,chi,Tout
         if (flm%cfg%amRoot) then
            ! dZmean=0.01_WP
            ! dZvar =0.0025_WP
            ! nZmean=int(1.0_WP/dZmean+1)
            ! nZvar =int(1.0_WP/dZvar +1)
            nZmean=flm%chmtbl%n1
            nZvar =flm%chmtbl%n2
            print*,'nZmean = ',nZmean,', nZvar = ',nZvar
            allocate(Zmean(nZmean,nZvar),Zvar(nZmean,nZvar),Tout(nZmean,nZvar),chi(nZmean,nZvar))
            ! do j=1,nZvar
            !    do i=1,nZmean
            !       Zmean(i,j)=0.0_WP+real(i-1,WP)*dZmean
            !       Zvar(i,j) =0.0_WP+real(j-1,WP)*dZvar
            !    end do
            ! end do
            do j=1,nZvar
               Zmean(:,j)=flm%chmtbl%x1
            end do
            do i=1,nZmean
               Zvar(i,:)=flm%chmtbl%x2
            end do
            chi=0.0_WP
            call flm%chmtbl%lookup('temperature',Tout,Zmean,Zvar,chi,nZmean*nZvar)
            open(newunit=ftest,file='T_lookedup.dat',status='replace',form='formatted',position='rewind')
            do j=1,nZvar
               do i=1,nZmean
                  write(ftest,'(3es12.5)') Zmean(i,j),Zvar(i,j),Tout(i,j)
               end do
            end do
            close(ftest)
         end if
      end block debug


      ! Get the inlet data
      inlet_data: block
         use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX
         use parallel, only: MPI_REAL_WP
         use vdscalar_class, only: bcond
         use, intrinsic :: iso_fortran_env, only: output_unit
         real(WP) :: myRe_jet,Re_jet,visc_jet
         integer  :: i,j,k,ierr
         type(bcond), pointer :: mybc
         visc_jet=0.0_WP
         rho_jet=0.0_WP
         rho_cof=0.0_WP
         myRe_jet=-1.0_WP
         call sc%get_bcond('jet',mybc)
         if (mybc%itr%n_.gt.0) then
            i=mybc%itr%map(1,1); j=mybc%itr%map(2,1); k=mybc%itr%map(3,1)
            ! Note: sc%SC(i,j,k) is 2, which is greater than 1 and thus, the following two statements give us the density and viscosity corresponding to Z=1
            rho_jet=sc%rho(i,j,k)
            visc_jet=fs%visc(i,j,k)
            myRe_jet=rho_jet*U_jet*D_jet/visc_jet
         end if
         call sc%get_bcond('coflow',mybc)
         if (mybc%itr%n_.gt.0) then
            i=mybc%itr%map(1,1); j=mybc%itr%map(2,1); k=mybc%itr%map(3,1)
            ! Note: As long as sc%SC(i,j,k) is between 0 and 1 (Z_cof<=0.5), the following statement is correct (accurate for Z_cof, but probabely not very good for Z_cof>0)
            rho_cof=0.5_WP*(sc%rho(i,j,k)+sc%rho(i+1,j,k))
         end if
         call MPI_ALLREDUCE(myRe_jet,Re_jet,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         if (cfg%amRoot) write(output_unit,'("Jet Reynolds = ",es12.5)') Re_jet
      end block inlet_data


      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         integer :: n,i,j,k
         type(bcond), pointer :: mybc
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Set density from scalar
         fs%rho=sc%rho
         ! Form momentum
         call fs%rho_multiply()
         ! Apply BCs
         call fs%apply_bcond(time%t,time%dt)
         call fs%get_bcond('jet',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=U_jet
            fs%rhoU(i,j,k)=rho_jet*U_jet
         end do
         call fs%get_bcond('coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%U(i,j,k)=U_cof
            fs%rhoU(i,j,k)=rho_cof*U_cof
         end do
         ! Get cell-centered velocities and continuity residual
         call fs%interp_vel(Ui,Vi,Wi)
         resSC=0.0_WP
         call fs%get_div(drhodt=resSC)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity


      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='dlra')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure'   ,fs%P)
         call ens_out%add_vector('velocity'   ,Ui,Vi,Wi)
         call ens_out%add_scalar('divergence' ,fs%div)
         call ens_out%add_scalar('density'    ,sc%rho)
         call ens_out%add_scalar('mixfrac'    ,sc%SC)
         call ens_out%add_scalar('temperature',T)
         call ens_out%add_scalar('ZgradMagSq',ZgradMagSq)
         call ens_out%add_scalar('chi',flm%chi)
         call ens_out%add_scalar('Zvar',flm%Zvar)
         do iY=1,n_Y
            call ens_out%add_scalar(Y_name(iY),Y(:,:,:,iY))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
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
         call mfile%add_column(sc%SCmax,'Zmax')
         call mfile%add_column(sc%SCmin,'Zmin')
         call mfile%add_column(sc%rhomax,'RHOmax')
         call mfile%add_column(sc%rhomin,'RHOmin')
         call mfile%add_column(int_RP,'Int(RP)')
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
         ! Create conservation monitor
         consfile=monitor(fs%cfg%amRoot,'conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(sc%SCint,'SC integral')
         call consfile%add_column(sc%rhoint,'RHO integral')
         call consfile%add_column(sc%rhoSCint,'rhoSC integral')
         call consfile%write()
      end block create_monitor


   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      integer :: i,j,k


      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold =sc%SC

         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Get transport properties and chemtable inputs
         call flm%chmtbl%lookup('viscosity',fs%visc,sc%SC,flm%Zvar,flm%chi,ncells)
         call flm%chmtbl%lookup('diffusivity',sc%diff,sc%SC,flm%Zvar,flm%chi,ncells)
         call sc%get_gradient(itpr_x=fs%itpr_x,itpr_y=fs%itpr_y,itpr_z=fs%itpr_z,SCgradMagSq=ZgradMagSq)
         call flm%get_Zvar(delta=sgs%delta,ZgradMagSq=ZgradMagSq,Z=sc%SC)
         ! Not sure which of the following is good
         call flm%get_chi(mueff=fs%visc,rho=sc%rho,ZgradMagSq=ZgradMagSq)
         ! call flm%get_chi(mueff=sc%diff,rho=sc%rho,ZgradMagSq=ZgradMagSq)
         
         ! Turbulence modeling
         sgs_modeling: block
            use sgsmodel_class, only: constant_smag
            resU=fs%rho
            ! call fs%get_gradu(gradU)
            call fs%get_strainrate(SR)
            ! call sgs%get_visc(type=vreman,dt=time%dtold,rho=resU,gradu=gradU)
            call sgs%get_visc(type=constant_smag,dt=time%dtold,rho=resU,SR=SR)
            call sgs%get_diff(type=constant_smag,diff_mol=sc%diff,rho=resU,SR=SR)
            fs%visc=fs%visc+sgs%visc
            sc%diff=sc%diff+sgs%diff
         end block sgs_modeling

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            if (cfg%amRoot) write(output_unit,'(" >  it/ itmax = ",i0,"/",i0)') time%it,time%itmax

            ! ============= SCALAR SOLVER ============= !
            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)

            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Assemble explicit residual
            resSC=time%dt*resSC-(2.0_WP*sc%rho*sc%SC-(sc%rho+sc%rhoold)*sc%SCold)

            ! Get the residual
            if (time_stepping.eq.'implicit') then
               ! Form implicit residual
               call sc%solve_implicit(time%dt,resSC,fs%rhoU,fs%rhoV,fs%rhoW)
            else
               ! Divide by density
               resSC=resSC/sc%rho
            end if

            ! Apply these residuals
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC

            ! Apply BCs on the resulting field
            call sc%apply_bcond(time%t,time%dt)
            dirichlet_scalar: block
               use vdscalar_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n
               call sc%get_bcond('jet',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k)=2.0_WP*Z_jet-sc%SC(i+1,j,k)
               end do
               call sc%get_bcond('coflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  sc%SC(i,j,k)=2.0_WP*Z_cof-sc%SC(i+1,j,k)
               end do
            end block dirichlet_scalar
            ! ========================================= !

            ! ============ UPDATE DENSITY ============ !
            ! Lookup density
            call flm%chmtbl%lookup('density',sc%rho,sc%SC,flm%Zvar,flm%chi,ncells)
            ! Smooth density
            if (nfilter.gt.0) then
               ! Compute density changes
               drho=sc%rho-sc%rhoold
               ! Filter density changes
               call filter_drho()
               ! Recompute new density
               sc%rho=sc%rhoold+drho
               ! Bound density
               if (rho_limiter) then
                  sc%rho=max(min(sc%rho,rho_max),rho_min)
                  call sc%cfg%sync(sc%rho)
               end if
            end if
            ! ========================================= !

            ! ============ VELOCITY SOLVER ============ !
            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)

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

            ! Get the residual
            if (time_stepping.eq.'implicit') then
               ! Form implicit residuals
               call fs%solve_implicit(time%dtmid,resU,resV,resW)
            else
               ! Divide by density
               do k=cfg%kmin_,cfg%kmax_+1
                  do j=cfg%jmin_,cfg%jmax_+1
                     do i=cfg%imin_,cfg%imax_+1
                        resU(i,j,k)=resU(i,j,k)/sum(fs%itpr_x(:,i,j,k)*fs%rho(i-1:i,j,k))
                        resV(i,j,k)=resV(i,j,k)/sum(fs%itpr_y(:,i,j,k)*fs%rho(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)/sum(fs%itpr_z(:,i,j,k)*fs%rho(i,j,k-1:k))
                     end do
                  end do
               end do
               call fs%cfg%sync(resU)
               call fs%cfg%sync(resV)
               call fs%cfg%sync(resW)
            end if

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Update momentum
            call fs%rho_multiply()

            ! Apply BCs
            call fs%apply_bcond(time%tmid,time%dtmid)
            dirichlet_velocity: block
               use lowmach_class, only: bcond
               type(bcond), pointer :: mybc
               integer :: n
               call fs%get_bcond('jet',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=U_jet
                  fs%rhoU(i,j,k)=rho_jet*U_jet
               end do
               call fs%get_bcond('coflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%U(i,j,k)=U_cof
                  fs%rhoU(i,j,k)=rho_cof*U_cof
               end do
            end block dirichlet_velocity

            ! Solve Poisson equation
            call sc%get_drhodt(dt=time%dt,drhodt=resSC)
            call fs%correct_mfr(drhodt=resSC)
            call fs%get_div(drhodt=resSC)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            call cfg%integrate(A=fs%psolv%rhs,integral=int_RP)
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide()
            ! ========================================= !

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resSC)
         call fs%get_div(drhodt=resSC)

         ! Output to ensight
         if (ens_evt%occurs()) then
            ! Combustion post-process
            call flm%chmtbl%lookup('temperature',T,sc%SC,flm%Zvar,flm%chi,ncells)
            do iY=1,n_Y
               call flm%chmtbl%lookup(Y_name(iY),Y(:,:,:,iY),sc%SC,flm%Zvar,flm%chi,ncells)
            end do
            ! wrtie
            call ens_out%write_data(time%t)
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%rho_multiply()
         call sc%get_int()
         call mfile%write()
         call cflfile%write()
         call consfile%write()

      end do

      
   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none


      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi,gradU,ZgradMagSq,drho,T,Y,Y_name)


   end subroutine simulation_final


end module simulation
