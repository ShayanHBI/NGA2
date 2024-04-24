!> Chemistry table class:
module chemtable_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   use messager,     only: die
   implicit none
   private
   
   ! Expose type/constructor/methods
   public :: chemtable

   !> Chemistry table object definition
   type :: chemtable

      ! This is our config
      class(config), pointer :: cfg                         !< This is the config the solver is build for

      ! This is the name of the flamelet object
      character(len=str_medium) :: name='UNNAMED_chemtable' !< Name (default=UNNAMED_chemtable)

      ! Chemistry table file name
      character(len=str_medium) :: filename                 !< Name of datafile to read/write

      ! Diffusion model name
      integer :: comb_model

      ! Mapping variables of the chemtable
      integer :: n1,n2,n3
      real(WP) :: x1min,x1max,x2min,x2max,x3min,x3max
      real(WP), dimension(:), allocatable :: x1,x2,x3

      ! Number of variables tabulated
      integer :: nvar
      character(len=str_medium), dimension(:), allocatable :: chem_name

      ! Arrays of mapped variables
      real(WP), dimension(:,:,:,:), allocatable :: table

      ! Index of density
      integer  :: index_rho

   contains
      procedure :: lookup
      procedure :: lookup_rho
      procedure :: lookup_max
      procedure :: lookup_min
      procedure :: get_var_ind
   end type chemtable

   !> Declare chemtable constructor
   interface chemtable
      procedure constructor
   end interface chemtable


contains


   !> Default constructor for chemtable object
   function constructor(cfg,fdata) result(self)
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      type(chemtable) :: self
      class(config), target, intent(in) :: cfg
      character(len=*), intent(in) :: fdata
      integer :: ierr,n,i,j,k
      type(MPI_File) :: ifile
      type(MPI_Status) :: status
      integer, dimension(4) :: dims

      ! Point to config object
      self%cfg=>cfg

      ! Open the chemtable file
      self%filename=trim(adjustl(fdata))
      call MPI_FILE_OPEN(self%cfg%comm,trim(self%filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[chemtable constructor] Problem encountered while reading chemistry table: '//trim(self%filename))

      ! Read the headers
      call MPI_FILE_READ_ALL(ifile,self%n1  ,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,self%n2  ,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,self%n3  ,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,self%nvar,1,MPI_INTEGER,status,ierr)

      ! Allocate the corresponding arrays
      allocate(self%x1(self%n1),self%x2(self%n2),self%x3(self%n3))
      allocate(self%chem_name(self%nvar))
      allocate(self%table(self%n1,self%n2,self%n3,self%nvar))

      ! Read the mapping variables
      call MPI_FILE_READ_ALL(ifile,self%x1,self%n1,MPI_REAL_WP,status,ierr)
      call MPI_FILE_READ_ALL(ifile,self%x2,self%n2,MPI_REAL_WP,status,ierr)
      call MPI_FILE_READ_ALL(ifile,self%x3,self%n3,MPI_REAL_WP,status,ierr)

      ! Read the combustion model used
      call MPI_FILE_READ_ALL(ifile,self%comb_model,1,MPI_INTEGER,status,ierr)

      ! Read the names of the mapped variables
      do n=1,self%nvar
         call MPI_FILE_READ_ALL(ifile,self%chem_name(n),str_medium,MPI_CHARACTER,status,ierr)
      end do

      ! Read the mapped variables
      do n=1,self%nvar
         call MPI_FILE_READ_ALL(ifile,self%table(:,:,:,n),self%n1*self%n2*self%n3,MPI_REAL_WP,status,ierr)
         ! Store index for density
         if (trim(self%chem_name(n)).eq.'density') self%index_rho=n
      end do

      ! Get some properties of the mapping
      self%x1min=minval(self%x1)
      self%x1max=maxval(self%x1)
      self%x2min=minval(self%x2)
      self%x2max=maxval(self%x2)
      self%x3min=minval(self%x3)
      self%x3max=maxval(self%x3)

      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
   end function constructor


   !< Look in the chemtable for the variable named 'tag'
   subroutine lookup(this,tag,R,A1,A2,A3,n)
      implicit none
      class(chemtable),       intent(inout) :: this
      character(len=*),       intent(in)    :: tag
      real(WP), dimension(n), intent(in)    :: A1,A2,A3
      real(WP), dimension(n), intent(out)   :: R
      integer, intent(in)                   :: n
      integer  :: var_ind,i,j
      integer  :: i1,i2,i3
      real(WP) :: w11,w12,w21,w22,w31,w32
      ! If density call another routine
      if (trim(tag).eq.'density') then
         call this%lookup_rho(R,A1,A2,A3,n)
         return
      end if
      ! Get the index of the variable
      var_ind=this%get_var_ind(tag)
      ! Trilinear interpolation
      do i=1,n
         ! First direction
         if (A1(i).lt.this%x1min) then
            i1=1
            w11=1.0_WP
         else if (A1(i).ge.this%x1max) then
            i1=this%n1-1
            w11=0.0_WP
         else
            loop1: do j=1,this%n1-1
               if (A1(i).lt.this%x1(j+1)) then
                  i1=j
                  exit loop1
               end if
            end do loop1
            w11=(this%x1(i1+1)-A1(i))/(this%x1(i1+1)-this%x1(i1))
         end if
         w12=1.0_WP-w11
         ! Second direction
         if (A2(i).lt.this%x2min) then
            i2=1
            w21=1.0_WP
         else if (A2(i).ge.this%x2max) then
            i2=this%n2-1
            w21=0.0_WP
         else
            loop2: do j=1,this%n2-1
               if (A2(i).lt.this%x2(j+1)) then
                  i2=j
                  exit loop2
               end if
            end do loop2
            w21=(this%x2(i2+1)-A2(i))/(this%x2(i2+1)-this%x2(i2))
         end if
         w22=1.0_WP-w21
         ! Third direction
         if (A3(i).lt.this%x3min) then
            i3=1
            w31=1.0_WP
         else if (A3(i).ge.this%x3max) then
            i3=this%n3-1
            w31=0.0_WP
         else
            loop3: do j=1,this%n3-1
               if (A3(i).lt.this%x3(j+1)) then
                  i3=j
                  exit loop3
               end if
            end do loop3
            w31=(this%x3(i3+1)-A3(i))/(this%x3(i3+1)-this%x3(i3))
         end if
         w32=1.0_WP-w31
         ! Interpolation
         R(i)=&
            w31*(w21*(w11*this%table(i1,i2,i3  ,var_ind)+w12*this%table(i1+1,i2,i3  ,var_ind))+w22*(w11*this%table(i1,i2+1,i3  ,var_ind)+w12*this%table(i1+1,i2+1,i3  ,var_ind)))+&
            w32*(w21*(w11*this%table(i1,i2,i3+1,var_ind)+w12*this%table(i1+1,i2,i3+1,var_ind))+w22*(w11*this%table(i1,i2+1,i3+1,var_ind)+w12*this%table(i1+1,i2+1,i3+1,var_ind)))
      end do
   end subroutine lookup


   !< Look in the chemtable for the density with different interpolation than for other variables
   subroutine lookup_rho(this,R,A1,A2,A3,n)
      implicit none
      class(chemtable),       intent(in)  :: this
      real(WP), dimension(n), intent(in)  :: A1,A2,A3
      real(WP), dimension(n), intent(out) :: R
      integer, intent(in)                 :: n
      integer  :: i,j
      integer  :: i1,i2,i3
      real(WP) :: w11,w12,w21,w22,w31,w32
      ! Trilinear interpolation
      do i=1,n
         ! First direction
         if (A1(i).lt.this%x1min) then
            i1=1
            w11=1.0_WP
         else if (A1(i).ge.this%x1max) then
            i1=this%n1-1
            w11=0.0_WP
         else
            loop1:do j=1,this%n1-1
               if (A1(i).lt.this%x1(j+1)) then
                  i1=j
                  exit loop1
               end if
            end do loop1
            w11=(this%x1(i1+1)-A1(i))/(this%x1(i1+1)-this%x1(i1))
         end if
         w12=1.0_WP-w11
         ! Second direction
         if (A2(i).lt.this%x2min) then
            i2=1
            w21=1.0_WP
         else if (A2(i).ge.this%x2max) then
            i2=this%n2-1
            w21=0.0_WP
         else
            loop2:do j=1,this%n2-1
               if (A2(i).lt.this%x2(j+1)) then
                  i2=j
                  exit loop2
               end if
            end do loop2
            w21=(this%x2(i2+1)-A2(i))/(this%x2(i2+1)-this%x2(i2))
         end if
         w22=1.0_WP-w21
         ! Third direction
         if (A3(i).lt.this%x3min) then
            i3=1
            w31=1.0_WP
         else if (A3(i).ge.this%x3max) then
            i3=this%n3-1
            w31=0.0_WP
         else
            loop3:do j=1,this%n3-1
               if (A3(i).lt.this%x3(j+1)) then
                  i3=j
                  exit loop3
               end if
            end do loop3
            w31=(this%x3(i3+1)-A3(i))/(this%x3(i3+1)-this%x3(i3))
         end if
         w32=1.0_WP-w31
         ! Interpolation of 1/rho
         R(i)=&
            w31*(w21*(w11/this%table(i1,i2,i3  ,this%index_rho)+w12/this%table(i1+1,i2,i3  ,this%index_rho))+w22*(w11/this%table(i1,i2+1,i3  ,this%index_rho)+w12/this%table(i1+1,i2+1,i3  ,this%index_rho)))+&
            w32*(w21*(w11/this%table(i1,i2,i3+1,this%index_rho)+w12/this%table(i1+1,i2,i3+1,this%index_rho))+w22*(w11/this%table(i1,i2+1,i3+1,this%index_rho)+w12/this%table(i1+1,i2+1,i3+1,this%index_rho)))
         R(i)=1.0_WP/R(i)
      end do
   end subroutine lookup_rho


   !< Find the maximum of a variable in the chemtable
   subroutine lookup_max(this,tag,R)
      implicit none
      class(chemtable), intent(inout) :: this
      character(len=*), intent(in) :: tag
      real(WP), intent(out) :: R
      integer :: var_ind
      ! Get the index of the variable
      var_ind=this%get_var_ind(tag)
      ! Return the max
      R=maxval(this%table(:,:,:,var_ind))
   end subroutine lookup_max


   !< Find the maximum of a variable in the chemtable
   subroutine lookup_min(this,tag,R)
      implicit none
      class(chemtable), intent(inout) :: this
      character(len=*), intent(in) :: tag
      real(WP), intent(out) :: R
      integer :: var_ind
      ! Get the index of the variable
      var_ind=this%get_var_ind(tag)
      ! Return the min
      R=minval(this%table(:,:,:,var_ind))
   end subroutine lookup_min


   !< Get the index of a variable
   function get_var_ind(this,tag) result(var_ind)
      implicit none
      class(chemtable), intent(inout) :: this
      character(len=*), intent(in) :: tag
      integer :: var_ind
      ! Get the index of the variable
      do var_ind=1,this%nvar
         if (trim(this%chem_name(var_ind)).eq.trim(tag)) exit
      end do
      if (var_ind.gt.this%nvar) then
         call die('[chemtable get_var_ind] unknown variable : '//trim(tag))
      end if
   end function


end module chemtable_class
