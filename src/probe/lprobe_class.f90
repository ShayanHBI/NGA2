module lprobe_class
   use string,       only: str_medium
   use precision,    only: WP
   use config_class, only: config
   use probe_class,  only: probe
   implicit none
   private

   public :: lprobe

   type :: lprobe
      character(len=str_medium)              :: name='UNNAMED_LINE_PROBE' !< Line probe name (default=UNNAMED_LINE_PROBE)
      class(config), pointer                 :: cfg                       !< The config that the line probe is build for
      integer                                :: np                        !< Number of points on the lprobe
      type(probe), dimension(:), allocatable :: probes                    !< Array of probes on the line
      integer, pointer                       :: method                    !< Interpolation methods for probe
   contains
      procedure :: fopen !< Open the lprobe data file
   end type lprobe

   !> Declare lprobe constructor
   interface lprobe
      procedure constructor
   end interface lprobe

contains

   !> Default constructor for lprobe
   function constructor(cfg,p1,p2,np,method,name) result(self)
      use mpi_f08 , only: MPI_MIN,MPI_ALLREDUCE
      use parallel, only: MPI_REAL_WP
      implicit none
      type(lprobe)                        :: self
      class(config), target, intent(in)   :: cfg
      real(WP), dimension(3), intent (in) :: p1, p2
      integer, intent(in)                 :: np
      integer, target, intent(in)         :: method
      character(len=*), optional          :: name
      integer                             :: i
      real(WP)                            :: x,y,z,L,dL,dx,dy,dz
      character(len=5)                    :: istr
      ! Set the name for the lprobe
      if (present(name)) self%name=trim(adjustl(name))
      ! Point to pgrid object
      self%cfg=>cfg
      ! Number of points on the lprobe
      self%np=np
      ! Interplation method
      self%method=>method
      ! Allocate memory for probes
      allocate(self%probes(np))
      ! Construct np probes on the line
      L=norm2(p2-p1)
      dL=L/real(np-1,WP)
      dx=dL*(p2(1)-p1(1))/L
      dy=dL*(p2(2)-p1(2))/L
      dz=dL*(p2(3)-p1(3))/L
      do i=1,self%np
         x=p1(1)+real(i-1,WP)*dx
         y=p1(2)+real(i-1,WP)*dy
         z=p1(3)+real(i-1,WP)*dz
         write(istr,"(I5)") i
         self%probes(i)=probe(self%cfg,x,y,z,self%method,trim(adjustl(self%name//'_p'//istr)))
      end do
   end function constructor

   !> Open the lprobe data file
   subroutine fopen(this,fdir,fname,fext,funit)
      implicit none
      class(lprobe), intent(inout) :: this
      character(len=*), intent(in) :: fdir,fname,fext
      integer, intent(out)         :: funit
      integer                      :: dir_status
      character(len=str_medium)    :: fdirname
      fdirname=trim(adjustl(trim(fdir)//trim(fname)//'.'//fext))
      dir_status = system("test -d "//trim(fdir))
      if (dir_status.ne.0) then
         dir_status = system("mkdir -p "//trim(fdir))
      end if
      open(newunit=funit,file=fdirname,status='replace',action='write',position='rewind')
   end subroutine fopen
end module lprobe_class
