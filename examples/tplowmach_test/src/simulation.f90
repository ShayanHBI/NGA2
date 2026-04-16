!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg,nx,Lx
   use hypre_str_class,   only: hypre_str
   use vfs_class,         only: vfs,VFlo
   use tplowmach_class,   only: tplowmach
   use tpvdscalar_class,  only: tpvdscalar,Lphase,Gphase
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private

   !> Solvers and time tracker
   type(hypre_str),   public :: ps,ss,vs
   type(vfs),         public :: vf
   type(tplowmach),   public :: fs
   type(tpvdscalar),  public :: sc
   type(timetracker), public :: time

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor files
   type(monitor) :: mfile,cflfile,hitfile,cvgfile,scfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: rho,div_src_rho,rho_lOld,rho_gOld,div_src_l,div_src_g
   real(WP), dimension(:,:,:,:), allocatable :: SR
   real(WP), dimension(:,:,:,:,:), allocatable :: gradU

   !> Fluid parameters
   real(WP) :: visc_l,visc_g           !< Dynamic viscosities
   real(WP) :: rho_l0,rho_l1           !< Liquid density EOS limits (Z=0 and Z=1)
   real(WP) :: rho_g0,rho_g1           !< Gas density EOS limits

   !> HIT forcing parameters and monitoring
   real(WP) :: meanU,meanV,meanW
   real(WP) :: Urms0,TKE0,EPS0,Re_max
   real(WP) :: TKE,URMS,EPS
   real(WP) :: tauinf,G,Gdtau,Gdtaui,dx
   real(WP) :: Re_L,Re_lambda,eta,ell
   real(WP) :: dx_eta,ell_Lx,Re_ratio,eps_ratio,tke_ratio,nondtime

   !> Problem definition
   real(WP) :: center(3),radius
   integer  :: iZl!,iZg

   ! Debug
   real(WP) :: prhs_int,mass,mass_old,dmass,div_mean


contains


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x- boundary for scalar fields
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator


   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y- boundary for scalar fields
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc


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


   !> Function that localizes z- boundary for scalar fields
   function zm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin-1) isIn=.true.
   end function zm_locator_sc


   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator


   function levelset_sphere(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=radius-sqrt(sum((xyz-center)**2))
   end function levelset_sphere


   !> Linear EOS: density varies linearly with a mixture fraction Z in [0,1]
   real(WP) function linear_rho(Z,rho0,rho1)
      real(WP), intent(in) :: Z,rho0,rho1
      real(WP) :: Zclip
      Zclip=min(max(Z,0.0_WP),1.0_WP)
      linear_rho=rho0+(rho1-rho0)*Zclip
   end function linear_rho


   !> Get phasic densities from the scalar field via EOS
   subroutine get_rho()
      use vfs_class, only: VFlo,VFhi
      integer :: i,j,k
      rho_lOld=fs%rho_l
      rho_gOld=fs%rho_g
      do k=cfg%kmino_,cfg%kmaxo_
         do j=cfg%jmino_,cfg%jmaxo_
            do i=cfg%imino_,cfg%imaxo_
               ! Liquid density
               if (vf%VF(i,j,k).ge.VFlo) then
                  fs%rho_l(i,j,k)=linear_rho(sc%SC(i,j,k,iZl),rho_l0,rho_l1)
                  sc%Prho(i,j,k,Lphase)=fs%rho_l(i,j,k)
               end if
               ! Gas density
               if (vf%VF(i,j,k).le.VFhi) then
                  ! fs%rho_g(i,j,k)=linear_rho(sc%SC(i,j,k,iZg),rho_g0,rho_g1)
                  fs%rho_g(i,j,k)=rho_g1
                  sc%Prho(i,j,k,Gphase)=fs%rho_g(i,j,k)
               end if
               rho(i,j,k)=vf%VF(i,j,k)*fs%rho_l(i,j,k)+(1.0_WP-vf%VF(i,j,k))*fs%rho_g(i,j,k)
            end do
         end do
      end do
   end subroutine get_rho


   !> Get the material derivative of a field
   ! function get_DDt(A,Aold) result(DADt)
   !    real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: A,Aold
   !    real(WP), dimension(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_) :: DADt
   !    integer :: i,j,k
   !    do k=cfg%kmin_,cfg%kmax_
   !       do j=cfg%jmin_,cfg%jmax_
   !          do i=cfg%imin_,cfg%imax_
   !             DADt(i,j,k)=(A(i,j,k)-Aold(i,j,k))/time%dt+Ui(i,j,k)*(A(i+1,j,k)-A(i-1,j,k))/(cfg%xm(i+1)-cfg%xm(i-1))+&
   !             &                                          Vi(i,j,k)*(A(i,j+1,k)-A(i,j-1,k))/(cfg%ym(j+1)-cfg%ym(j-1))+&
   !             &                                          Wi(i,j,k)*(A(i,j,k+1)-A(i,j,k-1))/(cfg%zm(k+1)-cfg%zm(k-1))
   !          end do
   !       end do
   !    end do
   ! end function get_DDt


   !> Get the velocity divergence induced by density variations
   function get_div_src_rho(myRHO,myRHOold) result(divsrc)
      real(WP), dimension(cfg%imino_:,cfg%jmino_:,cfg%kmino_:), intent(in) :: myRHO,myRHOold
      real(WP), dimension(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_) :: divsrc
      integer :: i,j,k
      divsrc=0.0_WP
      do k=cfg%kmin_,cfg%kmax_
         do j=cfg%jmin_,cfg%jmax_
            do i=cfg%imin_,cfg%imax_
               if (myRHO(i,j,k).gt.0.0_WP) then
                  divsrc(i,j,k)=divsrc(i,j,k)                                                      +&
                  &             (myRHO(i,j,k)-myRHOold(i,j,k))/time%dt                             +&
                  &             Ui(i,j,k)*(myRHO(i+1,j,k)-myRHO(i-1,j,k))/(cfg%xm(i+1)-cfg%xm(i-1))+&
                  &             Vi(i,j,k)*(myRHO(i,j+1,k)-myRHO(i,j-1,k))/(cfg%ym(j+1)-cfg%ym(j-1))+&
                  &             Wi(i,j,k)*(myRHO(i,j,k+1)-myRHO(i,j,k-1))/(cfg%zm(k+1)-cfg%zm(k-1))
                  divsrc(i,j,k)=-divsrc(i,j,k)/myRHO(i,j,k)
               end if
            end do
         end do
      end do
      ! Probably no need for syncing
      call cfg%sync(divsrc)
   end function get_div_src_rho


   !> Compute turbulence stats
   ! subroutine compute_stats()
   !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
   !    use parallel, only: MPI_REAL_WP
   !    real(WP) :: myTKE,myEPS
   !    integer :: i,j,k,ierr
   !    ! Interpolate to cell centers
   !    call fs%interp_vel(Ui,Vi,Wi)
   !    ! Compute mean velocities
   !    call fs%cfg%integrate(A=Ui,integral=meanU); meanU=meanU/fs%cfg%vol_total
   !    call fs%cfg%integrate(A=Vi,integral=meanV); meanV=meanV/fs%cfg%vol_total
   !    call fs%cfg%integrate(A=Wi,integral=meanW); meanW=meanW/fs%cfg%vol_total
   !    ! Compute strainrate
   !    call fs%get_strainrate(SR=SR)
   !    ! Compute TKE and dissipation (using mixture density)
   !    myTKE=0.0_WP; myEPS=0.0_WP
   !    do k=fs%cfg%kmin_,fs%cfg%kmax_
   !       do j=fs%cfg%jmin_,fs%cfg%jmax_
   !          do i=fs%cfg%imin_,fs%cfg%imax_
   !             myTKE=myTKE+0.5_WP*rho(i,j,k)*((Ui(i,j,k)-meanU)**2+(Vi(i,j,k)-meanV)**2+(Wi(i,j,k)-meanW)**2)*fs%cfg%vol(i,j,k)
   !             myEPS=myEPS+2.0_WP*fs%visc(i,j,k)*fs%cfg%vol(i,j,k)*(SR(1,i,j,k)**2+SR(2,i,j,k)**2+SR(3,i,j,k)**2+ &
   !             &    2.0_WP*(SR(4,i,j,k)**2+SR(5,i,j,k)**2+SR(6,i,j,k)**2))
   !          end do
   !       end do
   !    end do
   !    call MPI_ALLREDUCE(myTKE,TKE,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE=TKE/fs%cfg%vol_total
   !    call MPI_ALLREDUCE(myEPS,EPS,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPS=EPS/fs%cfg%vol_total
   !    ! Compute length/velocity scales (using mean kinematic viscosity estimate)
   !    Urms=sqrt(max(2.0_WP/3.0_WP*TKE,0.0_WP))
   !    Re_L=TKE**2.0_WP/max(EPS,1.0e-12_WP)/(visc_l/rho_l1)
   !    Re_lambda=sqrt(max(20.0_WP*Re_L/3.0_WP,0.0_WP))
   !    eta=((visc_l/rho_l1)**3.0_WP/max(EPS,1.0e-12_WP))**0.25_WP
   !    ell=(0.6667_WP*TKE)**1.5_WP/max(EPS,1.0e-12_WP)
   !    nondtime =time%t/tauinf
   !    dx_eta   =dx/max(eta,1.0e-12_WP)
   !    eps_ratio=EPS/max(EPS0,1.0e-12_WP)
   !    tke_ratio=TKE/max(TKE0,1.0e-12_WP)
   !    ell_Lx   =ell/Lx
   !    Re_ratio =Re_lambda/max(Re_max,1.0e-12_WP)
   ! end subroutine compute_stats


   !> Initialization of problem solver
   subroutine simulation_init
      use param,      only: param_read
      use mathtools,  only: Pi
      implicit none

      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1))
         allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(SR    (1:6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(gradU (1:3,1:3,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_src_rho(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); div_src_rho=0.0_WP
         allocate(rho_lOld(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(rho_gOld(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(div_src_l(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); div_src_l=0.0_WP
         allocate(div_src_g(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); div_src_g=0.0_WP
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker


      ! Initialize VOF solver with a droplet at the center
      create_and_initialize_vof: block
         use vfs_class, only: neumann
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,flux_storage
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3)   :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,name='VOF')
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Read droplet geometry
         call param_read('Droplet center',center,default=[0.5_WP*Lx,0.5_WP*Lx,0.5_WP*Lx])
         call param_read('Droplet radius',radius)
         ! Initialize VF field
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_sphere,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Correct mask outside physical domain
         if (cfg%iproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imin_-1
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%iproc.eq.cfg%npx) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imax_+1,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmin_-1
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.cfg%npy) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmax_+1,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.1) then
            do k=cfg%kmino_,cfg%kmin_-1
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.cfg%npz) then
            do k=cfg%kmax_+1,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         ! Apply boundary conditions
         call vf%apply_bcond(time%t,time%dt)
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()

      end block create_and_initialize_vof


      ! Create flow solver
      create_flow_solver: block
         use tplowmach_class, only: clipped_neumann
         use hypre_str_class, only: pcg_smg,pcg_pfmg2
         real(WP) :: r,m
         ! Create the flow solver
         fs=tplowmach(cfg=cfg,name='Two-phase LowMach NS')
         ! Read density and viscosity ratios
         call param_read('Density ratio',r)    !< rho_l / rho_g
         call param_read('Viscosity ratio',m)  !< mu_l / mu_g
         ! Densities
         call param_read('Liquid density',rho_l1)
         rho_l0=0.5_WP*rho_l1
         rho_g1=rho_l1/r
         rho_g0=rho_g1
         ! Viscosities
         call param_read('Liquid viscosity',visc_l)
         visc_g=visc_l/m
         fs%visc_l=visc_l; fs%visc_g=visc_g
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Boundary conditions
         call fs%add_bcond(name='xm',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)
         call fs%add_bcond(name='ym',type=clipped_neumann,face='y',dir=-1,canCorrect=.true.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=clipped_neumann,face='y',dir=+1,canCorrect=.true.,locator=yp_locator)
         call fs%add_bcond(name='zm',type=clipped_neumann,face='z',dir=-1,canCorrect=.true.,locator=zm_locator)
         call fs%add_bcond(name='zp',type=clipped_neumann,face='z',dir=+1,canCorrect=.true.,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_smg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         call param_read('Max coarsening levels',ps%maxlevel)
         ! ps=fft3d(cfg=cfg,name='Pressure',nst=7)
         ! Implicit velocity solver
         ! vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_pfmg2,nst=7)
         call param_read('Velocity iteration',vs%maxit)
         call param_read('Velocity tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
      end block create_flow_solver


      ! Create a one-sided scalar solver
      create_scalar: block
         use tpvdscalar_class, only: neumann
         use hypre_str_class,  only: pcg_pfmg2
         integer :: i,j,k
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=1,name='Scalar')
         iZl=1
         sc%SCname=['Zl']
         sc%phase(iZl)=Lphase
         ! Phasic VOF
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Diffusivity
         sc%diff(:,:,:,iZl)=fs%visc_l
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call sc%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call sc%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call sc%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call sc%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize the linear solver
         ! ss=hypre_str(cfg=cfg,name='Scalar',method=DDADIk,nst=7)
         ss=hypre_str(cfg=cfg,name='Scalar',method=pcg_pfmg2,nst=7)
         call param_read('Scalar iteration',ss%maxit)
         call param_read('Scalar tolerance',ss%rcvg)
         call sc%setup(implicit_solver=ss)
         ! Initialize scalar
         do k=sc%cfg%kmin_,sc%cfg%kmax_
            if ((sc%cfg%zm(k).gt.-0.25_WP).and.(sc%cfg%zm(k).lt.0.25_WP)) then
               do j=sc%cfg%jmin_,sc%cfg%jmax_
                  if ((sc%cfg%ym(j).gt.-0.25_WP).and.(sc%cfg%ym(j).lt.0.25_WP)) then
                     do i=sc%cfg%imin_,sc%cfg%imax_
                        ! if ((sc%cfg%xm(i).gt.-0.25_WP).and.(sc%cfg%xm(i).lt.0.25_WP)) sc%SC(i,j,k,iZl)=exp(-sc%cfg%xm(i)**2)
                        if ((sc%cfg%xm(i).gt.-0.25_WP).and.(sc%cfg%xm(i).lt.0.25_WP)) sc%SC(i,j,k,iZl)=1.0_WP
                     end do
                  end if
               end do
            end if
         end do
         call sc%cfg%sync(sc%SC(:,:,:,iZl))
         ! Compute face apertures
         call sc%get_face_apt()
         ! Compute density
         call get_rho()
         mass_old=mass
         call cfg%integrate(rho,mass)
         dmass=0.0_WP
      end block create_scalar


      ! Prepare initial velocity field: Gaussian random HIT
      initialize_velocity: block
         use random,    only: random_normal
         use param,     only: param_exists
         use messager,  only: die
         use tplowmach_class, only: static_contact,harmonic_visc
         use, intrinsic :: iso_fortran_env, only: output_unit
         integer :: i,j,k
         real(WP) :: taueta,nu
         nu=visc_l/rho_l1   !< Reference kinematic viscosity
         call param_read('Forcing constant',G)
         dx=Lx/real(nx,WP)
         if (param_exists('Steady-state TKE')) then
            call param_read('Steady-state TKE',TKE0)
            EPS0=5.0_WP*(0.6667_WP*TKE0)**1.5_WP/Lx
         else
            EPS0=nu**3*(Pi*cfg%nx/(1.5_WP*Lx))**4
            TKE0=1.5_WP*(0.2_WP*Lx*EPS0)**(0.6667_WP)
         end if
         Re_max=sqrt(15.0_WP*sqrt(0.6667_WP*TKE0)*0.2_WP*Lx/nu)
         tauinf=2.0_WP*TKE0/(3.0_WP*EPS0)
         taueta=sqrt(nu/EPS0)
         Gdtau =G/tauinf
         Gdtaui=1.0_WP/Gdtau
         if (Gdtaui.lt.time%dt) call die('[tplowmach_test] Forcing time constant < dt')
         Urms0=sqrt(0.6667_WP*TKE0)
         if (fs%cfg%amRoot) then
            write(output_unit,'("Expected turbulence properties:")')
            write(output_unit,'("Re_lambda = ",es12.5)') Re_max
            write(output_unit,'("tau_eddy  = ",es12.5)') tauinf
            write(output_unit,'("Urms      = ",es12.5)') Urms0
            write(output_unit,'("tau_eta   = ",es12.5)') taueta
         end if
         ! Random Gaussian velocity field
         do k=fs%cfg%kmin_,fs%cfg%kmax_
            do j=fs%cfg%jmin_,fs%cfg%jmax_
               do i=fs%cfg%imin_,fs%cfg%imax_
                  fs%U(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%V(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
                  fs%W(i,j,k)=random_normal(m=0.0_WP,sd=Urms0)
               end do
            end do
         end do
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         call fs%cfg%sync(fs%W)
         ! Remove mean
         call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
         call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
         fs%U=fs%U-meanU
         fs%V=fs%V-meanV
         fs%W=fs%W-meanW
         ! Projection
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         call fs%get_olddensity(vf=vf)
         fs%rho_U=fs%rho_Uold
         fs%rho_V=fs%rho_Vold
         fs%rho_W=fs%rho_Wold
         call fs%interp_vel(Ui,Vi,Wi)
         rho_lOld=fs%rho_l
         rho_gOld=fs%rho_g
         div_src_l=get_div_src_rho(fs%rho_l,rho_lOld)
         div_src_g=get_div_src_rho(fs%rho_g,rho_gOld)
         div_src_rho=vf%VF*div_src_l+(1.0_WP-vf%VF)*div_src_g
         call fs%update_laplacian()
         call fs%correct_mfr(src=div_src_rho)
         call fs%get_div(src=div_src_rho)
         call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
         ! call cfg%integrate(fs%div,div_mean)
         ! fs%div=fs%div-div_mean
         fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
         fs%psolv%sol=0.0_WP
         call fs%psolv%solve()
         call fs%shift_p(fs%psolv%sol)
         call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
         call cfg%integrate(fs%psolv%rhs,prhs_int)
         fs%P=fs%P+fs%psolv%sol
         fs%U=fs%U-time%dt*resU/fs%rho_U
         fs%V=fs%V-time%dt*resV/fs%rho_V
         fs%W=fs%W-time%dt*resW/fs%rho_W
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(src=div_src_rho)
      end block initialize_velocity


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: isc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='tplowmach_test')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('rho_l',fs%rho_l)
         call ens_out%add_scalar('rho_g',fs%rho_g)
         call ens_out%add_scalar('rho',rho)
         call ens_out%add_scalar('div_src_rho',div_src_rho)
         call ens_out%add_surface('plic',smesh)
         do isc=1,sc%nscalar
            call ens_out%add_scalar(trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create monitors
      create_monitor: block
         integer :: isc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         ! call compute_stats()
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
         call mfile%add_column(prhs_int,'prhs_int')
         call mfile%add_column(mass,'mass')
         call mfile%add_column(dmass,'dmass')
         call mfile%write()
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
         ! Create scalar monitor
         scfile=monitor(sc%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do isc=1,sc%nscalar
            call scfile%add_column(sc%SCmin(isc),trim(sc%SCname(isc))//' min')
            call scfile%add_column(sc%SCmax(isc),trim(sc%SCname(isc))//' max')
            call scfile%add_column(sc%SCint(isc),trim(sc%SCname(isc))//' int')
         end do
         call scfile%write()
         ! hitfile=monitor(fs%cfg%amRoot,'hit')
         ! call hitfile%add_column(time%n,'Timestep number')
         ! call hitfile%add_column(time%t,'Time')
         ! call hitfile%add_column(Re_L,'Re_L')
         ! call hitfile%add_column(Re_lambda,'Re_lambda')
         ! call hitfile%add_column(eta,'eta')
         ! call hitfile%add_column(TKE,'TKE')
         ! call hitfile%add_column(URMS,'Urms')
         ! call hitfile%add_column(EPS,'EPS')
         ! call hitfile%add_column(ell,'L')
         ! call hitfile%write()
         ! cvgfile=monitor(fs%cfg%amRoot,'convergence')
         ! call cvgfile%add_column(time%n,'Timestep number')
         ! call cvgfile%add_column(time%t,'Time')
         ! call cvgfile%add_column(nondtime,'Time/t_int')
         ! call cvgfile%add_column(Re_ratio,'Re_ratio')
         ! call cvgfile%add_column(eps_ratio,'EPS_ratio')
         ! call cvgfile%add_column(tke_ratio,'TKE_ratio')
         ! call cvgfile%add_column(dx_eta,'dx/eta')
         ! call cvgfile%add_column(ell_Lx,'ell/Lx')
         ! call cvgfile%write()
      end block create_monitor


   end subroutine simulation_init


   !> Time integrate our problem
   subroutine simulation_run
      use tplowmach_class, only: static_contact,harmonic_visc
      implicit none
      integer :: i,j,k
      integer :: isc,p

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         ! call time%adjust_dt()
         time%dtold=time%dt
         if (time%cfl.gt.0.0_WP) time%dt=time%dtold*time%cflmax/time%cfl
         call time%increment()

         ! Remember old velocity and density
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         rho_lOld=fs%rho_l
         rho_gOld=fs%rho_g

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! ================== VOF ================== !
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         call vf%apply_bcond(time%t,time%dt)

         ! ================== SCALAR ================== !
         advance_scalar: block
            ! Remember old SC
            sc%SCold =sc%SC
            sc%PVFold=sc%PVF

            ! Update phasic VOF from new interface
            sc%PVF(:,:,:,Lphase)=vf%VF
            sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
            call sc%get_face_apt()

            ! Correct for the emptied out and new interfacial cells
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               do k=cfg%kmino_,cfg%kmaxo_
                  do j=cfg%jmino_,cfg%jmaxo_
                     do i=cfg%imino_,cfg%imaxo_
                        if (sc%PVF(i,j,k,p).eq.0.0_WP) sc%SC(i,j,k,isc)=0.0_WP
                     end do
                  end do
               end do
               ! call lg%pure_zero_interfacial_extp(p,sc%SC(:,:,:,isc))
            end do

            ! Explicit calculation of d(PVOFrhoSC)/dt from scalar advection
            call sc%get_drhoSCdt_adv(drhoSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)

            ! Advance scalar advection
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               where (sc%mask.eq.0.and.sc%PVF(:,:,:,p).gt.0.0_WP) sc%SC(:,:,:,isc)=(sc%PVFold(:,:,:,p)*sc%SCold(:,:,:,isc)+time%dt*resSC(:,:,:,isc)/sc%Prho(:,:,:,p))/sc%PVF(:,:,:,p)
               where (sc%PVF(:,:,:,p).eq.0.0_WP) sc%SC(:,:,:,isc)=0.0_WP
            end do

            ! Explicit calculation of dVOFSC/dt from scalar diffusion
            call sc%get_drhoSCdt_dff(drhoSCdt=resSC)
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               where (sc%mask.eq.0.and.sc%PVF(:,:,:,p).gt.0.0_WP) resSC(:,:,:,isc)=time%dt*resSC(:,:,:,isc)/sc%PVF(:,:,:,p)
               where (sc%PVF(:,:,:,p).eq.0.0_WP) resSC(:,:,:,isc)=0.0_WP
            end do

            ! Form implicit diffusive residual
            call sc%solve_implicit_dff(time%dt,resSC)

            ! Advance scalar diffusion
            sc%SC=sc%SC+resSC

            ! Apply boundary conditions
            call sc%apply_bcond(time%t,time%dt)

         end block advance_scalar

         ! ================== UPDATE PROPERTIES ================== !

         ! Backup rhoSC
         ! do isc=1,sc%nscalar
         !    p=sc%phase(isc)
         !    resSC(:,:,:,isc)=sc%Prho(:,:,:,p)*sc%SC(:,:,:,isc)
         ! end do
         ! Update density
         call get_rho()
         mass_old=mass
         call cfg%integrate(rho,mass)
         dmass=mass-mass_old
         div_src_l=get_div_src_rho(fs%rho_l,rho_lOld)
         div_src_g=get_div_src_rho(fs%rho_g,rho_gOld)
         div_src_rho=vf%VF*div_src_l+(1.0_WP-vf%VF)*div_src_g
         ! Rescale scalar for conservation
         ! do isc=1,sc%nscalar
         !    p=sc%phase(isc)
         !    sc%SC(:,:,:,isc)=resSC(:,:,:,isc)/sc%Prho(:,:,:,p)
         ! end do
         ! UPDATE THE VISCOSITY
         ! UPDATE THE DIFFUSIVITY

         ! ================== VELOCITY ================== !

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)

            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Add momentum mass fluxes
            call fs%addsrc_gravity(resU,resV,resW)

            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

            ! Add linear forcing (Bassenne et al. 2016) to maintain TKE
            ! linear_forcing: block
            !    use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
            !    use parallel, only: MPI_REAL_WP
            !    real(WP) :: myTKE,myEPSp,EPSp,A
            !    integer :: ii,jj,kk,ierr
            !    call fs%cfg%integrate(A=fs%U,integral=meanU); meanU=meanU/fs%cfg%vol_total
            !    call fs%cfg%integrate(A=fs%V,integral=meanV); meanV=meanV/fs%cfg%vol_total
            !    call fs%cfg%integrate(A=fs%W,integral=meanW); meanW=meanW/fs%cfg%vol_total
            !    call fs%interp_vel(Ui,Vi,Wi)
            !    call fs%get_gradu(gradu=gradU)
            !    myTKE=0.0_WP; myEPSp=0.0_WP
            !    do kk=fs%cfg%kmin_,fs%cfg%kmax_
            !       do jj=fs%cfg%jmin_,fs%cfg%jmax_
            !          do ii=fs%cfg%imin_,fs%cfg%imax_
            !             myTKE =myTKE +0.5_WP*rho(ii,jj,kk)*((Ui(ii,jj,kk)-meanU)**2+(Vi(ii,jj,kk)-meanV)**2+(Wi(ii,jj,kk)-meanW)**2)*fs%cfg%vol(ii,jj,kk)
            !             myEPSp=myEPSp+fs%cfg%vol(ii,jj,kk)*fs%visc(ii,jj,kk)*(gradU(1,1,ii,jj,kk)**2+gradU(1,2,ii,jj,kk)**2+gradU(1,3,ii,jj,kk)**2+ &
            !             &                                                        gradU(2,1,ii,jj,kk)**2+gradU(2,2,ii,jj,kk)**2+gradU(2,3,ii,jj,kk)**2+ &
            !             &                                                        gradU(3,1,ii,jj,kk)**2+gradU(3,2,ii,jj,kk)**2+gradU(3,3,ii,jj,kk)**2)
            !          end do
            !       end do
            !    end do
            !    call MPI_ALLREDUCE(myTKE ,TKE ,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); TKE =TKE /fs%cfg%vol_total
            !    call MPI_ALLREDUCE(myEPSp,EPSp,1,MPI_REAL_WP,MPI_SUM,fs%cfg%comm,ierr); EPSp=EPSp/fs%cfg%vol_total
            !    A=(EPSp-Gdtau*(TKE-TKE0))/(2.0_WP*max(TKE,1.0e-12_WP))
            !    ! Force is applied at the momentum level
            !    resU=resU+time%dt*(fs%rho_U*fs%U-fs%rho_U*meanU)*A
            !    resV=resV+time%dt*(fs%rho_V*fs%V-fs%rho_V*meanV)*A
            !    resW=resW+time%dt*(fs%rho_W*fs%W-fs%rho_W*meanW)*A
            ! end block linear_forcing

            ! Form implicit residuals
            call fs%solve_implicit(time%dt,resU,resV,resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply boundary conditions
            call fs%apply_bcond(time%t,time%dt)

            ! Solve pressure Poisson
            call fs%update_laplacian()
            call fs%correct_mfr(src=div_src_rho)
            call fs%get_div(src=div_src_rho)
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            ! call cfg%integrate(fs%div,div_mean)
            ! fs%div=fs%div-div_mean
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct velocity (velocity-level, divide by face density)
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            call cfg%integrate(fs%psolv%rhs,prhs_int)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W

            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div(src=div_src_rho)

         ! Ensight output
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Monitoring
         ! call compute_stats()
         call fs%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         ! call hitfile%write()
         ! call cvgfile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      deallocate(resSC,resU,resV,resW,Ui,Vi,Wi,rho,div_src_rho,rho_lOld,rho_gOld,SR,gradU)
   end subroutine simulation_final


end module simulation
