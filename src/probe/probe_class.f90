module probe_class
   use string,       only: str_medium
   use precision,    only: WP
   use config_class, only: config
   implicit none
   private

   ! Expose type/constructor/methods
   public :: probe

   ! List of known available interpolation methods for probe
   integer, parameter, public :: nearest_cell=1

   !> Probe object definition
   type :: probe
      character(len=str_medium) :: name='UNNAMED_PROBE' !< Probe name (default=UNNAMED_PROBE)
      class(config), pointer    :: cfg                  !< The config that the probe is build for
      integer, pointer          :: method               !< Interpolation methods for probe
      integer                   :: rank                 !< Processor rank that contains the probe
      integer                   :: funit                !< File unit
      real(WP)                  :: x,y,z                !< The physical coordinates of the probe
      integer                   :: i,j,k                !< The local indices of the nearest cell to the probe
   contains
      procedure :: fopen  !< Open the probe data file
      procedure :: fwrite !< Write the probe data to file
   end type probe

   !> Declare probe constructor
   interface probe
      procedure constructor
   end interface probe

contains

   !> Default constructor for probe
   function constructor(cfg,x,y,z,method,name) result(self)
      use mpi_f08,  only: MPI_MIN,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      use messager, only: die
      implicit none
      type(probe)                       :: self
      class(config), target, intent(in) :: cfg
      real(WP), intent(in)              :: x,y,z
      integer, target, intent(in)       :: method
      character(len=*), optional        :: name
      real(WP)                          :: d,my_d_min,d_min
      integer                           :: i,j,k,ierr,my_rank
      ! Set the name for the probe
      if (present(name)) self%name=trim(adjustl(name))
      ! Point to pgrid object
      self%cfg=>cfg
      ! Interpolation method
      self%method=>method
      select case(self%method)
       case (nearest_cell)
         ! Find the nearest cell to the physical probe location and the corresponding processor ID
         self%rank=-1
         my_d_min=huge(my_d_min)
         do k=self%cfg%kmin_,self%cfg%kmax_
            do j=self%cfg%jmin_,self%cfg%jmax_
               do i=self%cfg%imin_,self%cfg%imax_
                  d=norm2([self%cfg%xm(i),self%cfg%ym(j),self%cfg%zm(k)]-[self%x,self%y,self%z])
                  if (d.lt.my_d_min) then
                     self%i=i; self%j=j; self%k=k;
                     my_d_min=d
                  end if
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_d_min,d_min,1,MPI_REAL_WP,MPI_MIN,self%cfg%comm,ierr)
         if (my_d_min.eq.d_min) then
            self%rank=self%cfg%rank
         end if
       case default
         call die('[probe constructor] Unknown unknown interpolation method')
      end select
   end function constructor

   !> Open the probe data file
   subroutine fopen(this,fdir,fname,reset,fext)
      implicit none
      class(probe), intent(inout)  :: this
      character(len=*), intent(in) :: fdir,fname,fext
      logical, intent(in)          :: reset
      logical                      :: fexists
      integer                      :: dir_status
      character(len=str_medium)    :: fdirname
      if (this%rank.eq.this%cfg%rank) then
         dir_status=system("test -d "//trim(fdir))
         if (dir_status.ne.0) then
            dir_status=system("mkdir -p "//trim(fdir))
         end if
         fdirname=trim(adjustl(trim(fdir)//trim(fname)//'.'//fext))
         inquire(file=trim(fdirname),exist=fexists)
         if ((.not.reset).and.(fexists)) then
            open(newunit=this%funit,file=trim(fdirname),status='old',action='write',position='append')
         else
            open(newunit=this%funit,file=trim(fdirname),status='replace',action='write',position='rewind')
         end if
         print*,"file unit ", this%funit
      end if
   end subroutine fopen

   !> Write the probe data
   subroutine fwrite(this,arr)
      implicit none
      class(probe), intent(inout)        :: this
      real(WP), dimension(:), intent(in) :: arr
      integer                            :: arr_ind
      if (this%rank.eq.this%cfg%rank) then
         do arr_ind=1,size(arr)
            write(this%funit,'(es12.5)',advance='no') arr(arr_ind)
            write(this%funit,'(a4)',advance='no') '    '
         end do
         write(this%funit,*)
         flush(this%funit)
      end if
   end subroutine fwrite
end module probe_class
