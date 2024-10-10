!> Evaporation class:
module evap_class
   use precision,         only: WP
   use string,            only: str_medium
   use config_class,      only: config
   use timetracker_class, only: timetracker
   use vfs_class,         only: vfs
   ! Debug
   use monitor_class,     only: monitor
   use ensight_class,     only: ensight
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: evap
   
   ! Index shift
   integer, dimension(1:3,1:3) :: ind_shift=reshape([1,0,0,0,1,0,0,0,1], shape(ind_shift))

   type :: arr_ptr_4d
      real(WP), pointer :: arr(:,:,:,:)
   end type arr_ptr_4d

   !> Evaporation object definition
   type :: evap
      
      ! Debug
      type(monitor) :: mfile
      type(ensight) :: ens_out
      real(WP), dimension(:,:,:), allocatable :: mfluxL_errField,mfluxG_errField

      ! This is our config
      class(config), pointer :: cfg                                    !< This is the config the object is build for
      
      ! This is the name of the object
      character(len=str_medium) :: name='UNNAMED_EVAP'                 !< Object name (default=UNNAMED_EVAP)

      ! This is the corresponding VOF solver
      class(vfs), pointer :: vf

      ! Liquid and gas densities
      real(WP) :: rho_l,rho_g

      ! Evaporation mass flux data
      real(WP), dimension(:,:,:), allocatable :: mdotdp                !< Evaporation mass flux (m dot double prime of dimension M/(T*L^2))
      real(WP), dimension(:,:,:), allocatable :: mflux                 !< Evaporation mass flux scaled by the surface density
      real(WP), dimension(:,:,:), allocatable :: mfluxL,mfluxL_old     !< Liquid side shifted evaporation mass flux scaled by the surface density
      real(WP), dimension(:,:,:), allocatable :: mfluxG,mfluxG_old     !< Gas side shifted evaporation mass flux scaled by the surface density
      real(WP), dimension(:,:,:), allocatable :: evp_div               !< Evaporatin source term (div(U) = evp_div)

      ! Pseudo time over which the mflux is being shifted
      type(timetracker), public :: pseudo_time

      ! Phase-change and interfacial velocity
      real(WP), dimension(:,:,:,:), allocatable :: normal              !< Interface normal vector
      real(WP), dimension(:,:,:,:), allocatable :: vel_pc              !< Phase-change velocity
      real(WP), dimension(:,:,:),   allocatable :: U_itf,V_itf,W_itf   !< Interfacial velocity components
      real(WP), dimension(:,:,:,:), allocatable :: pseudo_vel          !< Pseudo velocity for shifting mflux

      ! Metrics (point to flow solver metrics)
      type(arr_ptr_4d), dimension(:), allocatable :: itp               !< Cell to face interpolation coefficients
      type(arr_ptr_4d), dimension(:), allocatable :: div               !< Divergence operator (cell-centered)

      ! Number of cells in each direction
      integer, dimension(3) :: nCell
      
      ! Monitoring quantities
      real(WP) :: mflux_int,mflux_tol                                  !< Integral and tolerence of the scaled evap mass flux
      real(WP) :: mfluxL_int,mfluxL_err,mfluxL_int_err                 !< Liquid side scaled evap mass flux maximum, integral, and error
      real(WP) :: mfluxG_int,mfluxG_err,mfluxG_int_err                 !< Gas side scaled evap mass flux maximum, integral, and error
      real(WP) :: Upcmax,Vpcmax,Wpcmax                                 !< Maximum of phase-change velocity components
      
   contains

      procedure :: initialize                                          !< Class initializer
      procedure :: get_normal                                          !< Get the interface normal vector
      procedure :: get_vel_pc                                          !< Get the phase-change velocity
      procedure :: get_max_vel_pc                                      !< Calculate maximum field values
      procedure :: get_pseudo_vel                                      !< Get the face-centered normilized gradient of VOF
      procedure :: get_dmfluxdt                                        !< Get the time derivative of the evaporation mass fllux
      procedure :: shift_mflux                                         !< Shift the evaporation mass flux
      procedure :: get_div                                             !< Get the evaporation source term
      procedure :: get_cfl                                             !< Get the CFL

   end type evap

   
contains
   
   
   !> Default constructor for evap class
   subroutine initialize(this,cfg,vf,itp_x,itp_y,itp_z,div_x,div_y,div_z,name)
      implicit none
      class(evap), intent(inout) :: this
      class(config), target, intent(in) :: cfg
      class(vfs), target, intent(in) :: vf
      real(WP), target, dimension(-1: 0,cfg%imino_+1:cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: itp_x
      real(WP), target, dimension(-1: 0,cfg%imino_  :cfg%imaxo_,cfg%jmino_+1:cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: itp_y
      real(WP), target, dimension(-1: 0,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_+1:cfg%kmaxo_), intent(in) :: itp_z
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_x
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_y
      real(WP), target, dimension( 0:+1,cfg%imino_  :cfg%imaxo_,cfg%jmino_  :cfg%jmaxo_,cfg%kmino_  :cfg%kmaxo_), intent(in) :: div_z
      character(len=*), optional :: name
      integer :: i,j,k
      
      ! Set the name for the solver
      if (present(name)) this%name=trim(adjustl(name))
      
      ! Point to pgrid object
      this%cfg=>cfg

      ! Point to VOF solver
      this%vf=>vf
      
      ! Allocate variables
      allocate(this%mdotdp    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mdotdp    =0.0_WP
      allocate(this%mflux     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mflux     =0.0_WP
      allocate(this%mfluxL    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mfluxL    =0.0_WP
      allocate(this%mfluxG    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mfluxG    =0.0_WP
      allocate(this%mfluxL_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mfluxL_old=0.0_WP
      allocate(this%mfluxG_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%mfluxG_old=0.0_WP
      allocate(this%evp_div   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%evp_div   =0.0_WP
      allocate(this%normal    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3)); this%normal    =0.0_WP
      allocate(this%vel_pc    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3)); this%vel_pc    =0.0_WP
      allocate(this%pseudo_vel(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:3)); this%pseudo_vel=0.0_WP
      allocate(this%U_itf     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%U_itf     =0.0_WP
      allocate(this%V_itf     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%V_itf     =0.0_WP
      allocate(this%W_itf     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_));     this%W_itf     =0.0_WP
      allocate(this%itp(3))
      allocate(this%div(3))
      ! Debug
      allocate(this%mfluxL_errField(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); this%mfluxL_errField=0.0_WP
      allocate(this%mfluxG_errField(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); this%mfluxG_errField=0.0_WP

      ! Set metrics
      this%itp(1)%arr=>itp_x
      this%itp(2)%arr=>itp_y
      this%itp(3)%arr=>itp_z
      this%div(1)%arr=>div_x
      this%div(2)%arr=>div_y
      this%div(3)%arr=>div_z

      ! Number of cells in each direction
      this%nCell(1)=this%cfg%nx
      this%nCell(2)=this%cfg%ny
      this%nCell(3)=this%cfg%nz

      ! Create a pseudo time
      this%pseudo_time=timetracker(amRoot=this%cfg%amRoot,name='Pseudo',print_info=.false.)

      ! Debug
      ! this%mfile=monitor(cfg%amRoot,'mflux')
      ! call this%mfile%add_column(this%pseudo_time%n,'Pseudo time step')
      ! call this%mfile%add_column(this%pseudo_time%t,'Pseudo time')
      ! call this%mfile%add_column(this%mfluxL_err,'Liquid error')
      ! call this%mfile%add_column(this%mfluxG_err,'Gas error')
      ! this%ens_out=ensight(cfg=cfg,name='mflux')
      ! call this%ens_out%add_scalar('mflux',this%mflux)
      ! call this%ens_out%add_scalar('mfluxL',this%mfluxL)
      ! call this%ens_out%add_scalar('mfluxG',this%mfluxG)
      ! call this%ens_out%add_scalar('mfluxL_err',this%mfluxL_errField)
      ! call this%ens_out%add_scalar('mfluxG_err',this%mfluxG_errField)

   end subroutine initialize
   

   !> Calculate the interfacial normal vector
   subroutine get_normal(this)
      use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
      implicit none
      class(evap), intent(inout) :: this
      real(WP), dimension(3) :: n1,n2
      integer :: i,j,k,dir
      ! Loop over cells
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               n1=calculateNormal(this%vf%interface_polygon(1,i,j,k))
               if (getNumberOfVertices(this%vf%interface_polygon(2,i,j,k)).gt.0) then
                  n2=calculateNormal(this%vf%interface_polygon(2,i,j,k))
                  n1=0.5_WP*(n1+n2)
               end if
               do dir=1,3
                  this%normal(i,j,k,dir)=n1(dir)
               end do
            end do
         end do
      end do
   end subroutine get_normal

   
   !> Calculate the face-centered phase-change velocity
   subroutine get_vel_pc(this)
      use vfs_class, only: vfs
      implicit none
      class(evap), intent(inout) :: this
      real(WP) :: VFm,VFp
      logical  :: is_interfacial_m,is_interfacial_p
      integer  :: i,j,k,dir
      integer  :: im,jm,km
      integer  :: ip,jp,kp
      
      ! Initialize with zeros
      this%vel_pc=0.0_WP
      
      ! Loop over directions
      do dir=1,3
         ! Skip if needed
         if (this%nCell(dir).eq.1) cycle
         ! Loop over cell faces
         do k=this%cfg%kmin_,this%cfg%kmax_+ind_shift(3,dir)
            do j=this%cfg%jmin_,this%cfg%jmax_+ind_shift(2,dir)
               do i=this%cfg%imin_,this%cfg%imax_+ind_shift(1,dir)
                  ! Prepare indices for the adjacent cells
                  im=i-ind_shift(1,dir); ip=i
                  jm=j-ind_shift(2,dir); jp=j
                  km=k-ind_shift(3,dir); kp=k
                  ! Get the corresponding VOF values
                  VFm=this%vf%VF(im,jm,km)
                  VFp=this%vf%VF(ip,jp,kp)
                  ! Check if the adjacent cells are interfacial
                  is_interfacial_m=VFm.gt.0.0_WP.and.VFm.lt.1.0_WP
                  is_interfacial_p=VFp.gt.0.0_WP.and.VFp.lt.1.0_WP
                  if (is_interfacial_m) then
                     if (is_interfacial_p) then
                        ! Both cells are interfacial, linear interpolation of the phase-change velocity
                        this%vel_pc(i,j,k,dir)=(this%itp(dir)%arr(-1,i,j,k)*this%mdotdp(im,jm,km)*this%normal(im,jm,km,dir) &
                        &                      +this%itp(dir)%arr( 0,i,j,k)*this%mdotdp(ip,jp,kp)*this%normal(ip,jp,kp,dir))/this%rho_l
                     else
                        ! The plus cell is not interfacial, use the minus cell's phase-change velocity
                        this%vel_pc(i,j,k,dir)=this%mdotdp(im,jm,km)*this%normal(im,jm,km,dir)/this%rho_l
                     end if
                  else if (is_interfacial_p) then
                     ! The minus cell is not interfacial, use the plus cell's phase-change velocity
                     this%vel_pc(i,j,k,dir)=this%mdotdp(ip,jp,kp)*this%normal(ip,jp,kp,dir)/this%rho_l
                  end if
               end do
            end do
         end do
         ! Sync the phase-change velocity component
         call this%cfg%sync(this%vel_pc(:,:,:,dir))
      end do

   end subroutine get_vel_pc


   !> Calculate the pseudo velocity used to shift the evaporation mass flux
   subroutine get_pseudo_vel(this)
      implicit none
      class(evap), intent(inout) :: this
      integer :: i,j,k,dir
      integer :: im,jm,km
      integer :: ip,jp,kp
      real(WP):: VFm,VFp
      logical :: is_interfacial_m,is_interfacial_p
      
      ! Initialize with zeros
      this%pseudo_vel=0.0_WP

      ! Loop over directions
      do dir=1,3
         ! Skip if needed
         if (this%nCell(dir).eq.1) cycle
         ! Loop over cell faces
         do k=this%cfg%kmin_,this%cfg%kmax_+ind_shift(3,dir)
            do j=this%cfg%jmin_,this%cfg%jmax_+ind_shift(2,dir)
               do i=this%cfg%imin_,this%cfg%imax_+ind_shift(1,dir)
                  ! Prepare indices for the adjacent cells
                  im=i-ind_shift(1,dir); ip=i
                  jm=j-ind_shift(2,dir); jp=j
                  km=k-ind_shift(3,dir); kp=k
                  ! Get the corresponding VOF values
                  VFm=this%vf%VF(im,jm,km)
                  VFp=this%vf%VF(ip,jp,kp)
                  ! Check if the adjacent cells are interfacial
                  is_interfacial_m=VFm.gt.0.0_WP.and.VFm.lt.1.0_WP
                  is_interfacial_p=VFp.gt.0.0_WP.and.VFp.lt.1.0_WP
                  if (is_interfacial_m) then
                     if (is_interfacial_p) then
                        ! Both cells are interfacial, linear interpolation of the normal vector
                        this%pseudo_vel(i,j,k,dir)=-this%itp(dir)%arr(-1,i,j,k)*this%normal(im,jm,km,dir) &
                        &                          -this%itp(dir)%arr( 0,i,j,k)*this%normal(ip,jp,kp,dir)
                     else
                        ! The plus cell is not interfacial, use the minus cell's normal vector
                        this%pseudo_vel(i,j,k,dir)=-this%normal(im,jm,km,dir)
                     end if
                  else if (is_interfacial_p) then
                     ! The minus cell is not interfacial, use the plus cell's normal vector
                     this%pseudo_vel(i,j,k,dir)=-this%normal(ip,jp,kp,dir)
                  end if
               end do
            end do
         end do
      end do

   end subroutine get_pseudo_vel
   

   !> Calculate the explicit mflux time derivative
   subroutine get_dmfluxdt(this,vel,mflux_old,dmfluxdt)
      implicit none
      class(evap), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:,1:), intent(in)  :: vel        !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_,1:3)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(in)  :: mflux_old  !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:),    intent(out) :: dmfluxdt   !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      integer :: i,j,k,dir
      integer :: im,jm,km
      integer :: ip,jp,kp
      real(WP) :: F
      ! Zero out dmflux/dt array
      dmfluxdt=0.0_WP
      ! Fluxes of mflux
      do dir=1,3
         ! Loop over cell faces
         do k=this%cfg%kmin_,this%cfg%kmax_+ind_shift(3,dir)
            do j=this%cfg%jmin_,this%cfg%jmax_+ind_shift(2,dir)
               do i=this%cfg%imin_,this%cfg%imax_+ind_shift(1,dir)
                  ! Prepare indices for the adjacent cells
                  im=i-ind_shift(1,dir); ip=i
                  jm=j-ind_shift(2,dir); jp=j
                  km=k-ind_shift(3,dir); kp=k
                  ! Compute the face flux
                  F=-0.5_WP*(vel(i,j,k,dir)+abs(vel(i,j,k,dir)))*mflux_old(im,jm,km) &
                  & -0.5_WP*(vel(i,j,k,dir)-abs(vel(i,j,k,dir)))*mflux_old(ip,jp,kp)
                  ! Add it to the time derivative of mflux
                  dmfluxdt(im,jm,km)=dmfluxdt(im,jm,km)+this%div(dir)%arr(1,im,jm,km)*F
                  dmfluxdt(ip,jp,kp)=dmfluxdt(ip,jp,kp)+this%div(dir)%arr(0,ip,jp,kp)*F
               end do
            end do
         end do
      end do
      ! Sync residual
      call this%cfg%sync(dmfluxdt)
   end subroutine get_dmfluxdt


   !> Shift mflux away from the interface
   subroutine shift_mflux(this)
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_MAX
      use parallel,  only: MPI_REAL_WP
      use vfs_class, only: vfs
      implicit none
      class(evap), intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: resmfluxL,resmfluxG
      integer  :: ierr
      real(WP) :: mflux_max,my_mflux_max,mflux_err,my_mflux_err

      ! Allocate memory for mflux residuals
      allocate(resmfluxL(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(resmfluxG(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Get the normalized gradient of VOF
      call this%get_pseudo_vel()

      ! Get the CFL based on the gradient of the VOF
      call this%get_cfl(this%pseudo_time%dt,this%pseudo_vel(:,:,:,1),this%pseudo_vel(:,:,:,2),this%pseudo_vel(:,:,:,3),this%pseudo_time%cfl)
      
      ! Reset the pseudo time
      call this%pseudo_time%reset()
      
      ! Adjust the pseudo time step
      call this%pseudo_time%adjust_dt()
      
      ! Initialize the evaporation mass fluxes on the liquid and gas sides
      this%mfluxL=this%mflux; this%mfluxG=this%mflux

      ! Move the evaporation mass flux away from the interface
      do while (.not.this%pseudo_time%done())
         
         ! Remember old mflux
         this%mfluxL_old=this%mfluxL
         this%mfluxG_old=this%mfluxG
         
         ! Increment pseudo time
         call this%pseudo_time%increment()
         
         ! Assemble explicit residual
         call this%get_dmfluxdt(vel= this%pseudo_vel,mflux_old=this%mfluxL_old,dmfluxdt=resmfluxL)
         call this%get_dmfluxdt(vel=-this%pseudo_vel,mflux_old=this%mfluxG_old,dmfluxdt=resmfluxG)
         
         ! Apply these residuals
         this%mfluxL=this%mfluxL_old+this%pseudo_time%dt*resmfluxL
         this%mfluxG=this%mfluxG_old+this%pseudo_time%dt*resmfluxG

         ! Get the maximum of mflux
         my_mflux_max=maxval(this%mflux)
         call MPI_ALLREDUCE(my_mflux_max,mflux_max,1,MPI_REAL_WP,MPI_Max,this%cfg%comm,ierr)

         ! Calculate the error on the liquid side
         my_mflux_err=maxval(abs((this%mfluxL-this%mfluxL_old)/mflux_max))
         call MPI_ALLREDUCE(my_mflux_err,this%mfluxL_err,1,MPI_REAL_WP,MPI_Max,this%cfg%comm,ierr)
         
         ! Calculate the error on the gas side
         my_mflux_err=maxval(abs((this%mfluxG-this%mfluxG_old)/mflux_max))
         call MPI_ALLREDUCE(my_mflux_err,this%mfluxG_err,1,MPI_REAL_WP,MPI_Max,this%cfg%comm,ierr)

         ! Debug
         ! this%mfluxL_errField=abs((this%mfluxL-this%mfluxL_old)/mflux_max)
         ! this%mfluxG_errField=abs((this%mfluxG-this%mfluxG_old)/mflux_max)
         ! call this%mfile%write()
         ! call this%ens_out%write_data(this%pseudo_time%t)

         ! Check convergence
         mflux_err=max(this%mfluxL_err,this%mfluxG_err)
         if (mflux_err.lt.this%mflux_tol) exit

      end do

      ! Integral of mflux
      call this%cfg%integrate(this%mflux,this%mflux_int)
      call this%cfg%integrate(this%mfluxL,this%mfluxL_int)
      call this%cfg%integrate(this%mfluxG,this%mfluxG_int)
      this%mfluxL_int_err=abs(this%mfluxL_int-this%mflux_int)
      this%mfluxG_int_err=abs(this%mfluxG_int-this%mflux_int)

      ! Deallocate mflux residuals
      deallocate(resmfluxL,resmfluxG)

   end subroutine shift_mflux


   !> Calculate the divergence induced by phase change
   subroutine get_div(this)
      implicit none
      class(evap), intent(inout) :: this
      this%evp_div=(this%mfluxG/this%rho_g-this%mfluxL/this%rho_l)
   end subroutine get_div


   !> Calculate the CFL
   subroutine get_cfl(this,dt,U,V,W,cfl)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
      use parallel, only: MPI_REAL_WP
      implicit none
      class(evap), intent(inout) :: this
      real(WP), intent(in)  :: dt
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: U     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: V     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), intent(in) :: W     !< Needs to be (imino_:imaxo_,jmino_:jmaxo_,kmino_:kmaxo_)
      real(WP), intent(out) :: cfl
      integer :: i,j,k,ierr
      real(WP) :: my_CFL
      
      ! Set the CFL to zero
      my_CFL=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_CFL=max(my_CFL,abs(U(i,j,k))*this%cfg%dxmi(i))
               my_CFL=max(my_CFL,abs(V(i,j,k))*this%cfg%dymi(j))
               my_CFL=max(my_CFL,abs(W(i,j,k))*this%cfg%dzmi(k))
            end do
         end do
      end do
      my_CFL=my_CFL*dt
      
      ! Get the parallel max
      call MPI_ALLREDUCE(my_CFL,cfl,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      
   end subroutine get_cfl


   !> Calculate the maximum of the phase-change velocity
   subroutine get_max_vel_pc(this)
      use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(evap), intent(inout) :: this
      real(WP) :: my_Upcmax,my_Vpcmax,my_Wpcmax
      integer  :: i,j,k,ierr
      ! Set all to zero
      my_Upcmax=0.0_WP; my_Vpcmax=0.0_WP; my_Wpcmax=0.0_WP
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               my_Upcmax=max(my_Upcmax,abs(this%vel_pc(i,j,k,1)))
               my_Vpcmax=max(my_Vpcmax,abs(this%vel_pc(i,j,k,2)))
               my_Wpcmax=max(my_Wpcmax,abs(this%vel_pc(i,j,k,3)))
            end do
         end do
      end do
      ! Get the parallel max
      call MPI_ALLREDUCE(my_Upcmax,this%Upcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Vpcmax,this%Vpcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_Wpcmax,this%Wpcmax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
   end subroutine get_max_vel_pc


end module evap_class
