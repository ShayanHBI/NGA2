!> Flamelet library class:
module flameletLib_class
   use precision, only: WP
   use string,    only: str_medium,str_long
   use messager,  only: die
   implicit none
   private

   !< Expose type/constructor/methods
   public :: flameletLib,modify_varname

   ! List of known available flamelet models
   integer, parameter, public :: sfm=1                       !< Steady Flamelet Model (SFM)

   !> Flamelet library object definition
   type :: flameletLib

      ! Combustiopn model
      integer :: combModel

      ! Number of points in the flamelet
      integer :: nPoints

      ! Coordinate in mixture fraction space
      real(WP), dimension(:), pointer :: Z

      ! List of names/variables to get from FlameMaster files
      integer :: nvar_in
      character(len=str_medium), dimension(:), allocatable :: input_name
      real(WP), dimension(:,:), pointer :: input_data
      logical, dimension(:), allocatable :: found

      ! List of Flamelet files
      integer :: nfiles
      character(len=str_medium), dimension(:), allocatable :: files

   contains
      procedure :: readfile
      procedure :: cleanup
   end type flameletLib

   !> Declare flameletLib constructor
   interface flameletLib
      procedure constructor
   end interface flameletLib

contains


   !> Default constructor for flameletLib object
   function constructor(model) result(self)
      use param, only: param_getsize,param_read
      implicit none
      type(flameletLib) :: self
      integer :: model,var

      ! Read the list of files
      self%nfiles=param_getsize('List of Flamelets')
      allocate(self%files(self%nfiles))
      call param_read('List of Flamelets',self%files)

      ! Get number of additionnal variables to store in the table
      self%nvar_in=param_getsize('FlameMaster variables')

      self%combModel=model
      select case(self%combModel)
       case (sfm)
         allocate(self%input_name(self%nvar_in+1))
         call param_read('FlameMaster variables',self%input_name(1:self%nvar_in))
         self%nvar_in=self%nvar_in+1
         self%input_name(self%nvar_in)='chi'
       case default
         call die('[flameletlib constructor] Unknown combustion model')
      end select

      ! Make the mass fraction names compatible with FlameMaster files
      do var=1,self%nvar_in
         call modify_varname(self%input_name(var),'Y_','massfraction-')
         call modify_varname(self%input_name(var),'Y-','massfraction-')
         call modify_varname(self%input_name(var),'massfraction_','massfraction-')
      end do

      ! Allocate array to specify wether the variables have been found
      allocate(self%found(self%nvar_in))
   end function constructor

   subroutine readfile(this,ifile)
      implicit none
      class(flameletLib), intent(inout) :: this
      integer, intent(in) :: ifile
      integer :: iunit,ierr,var,n,nlines,index1
      character(len=str_long) :: buffer
      character(len=str_long) :: line
      character(len=str_medium) :: varname
      real(WP), dimension(:), pointer :: tmp

      ! Open the file
      open(newunit=iunit,file=trim(this%files(ifile)),form='formatted',status='old',iostat=ierr)
      if (ierr.ne.0) then
         call die("[flameletLib readfile] Error opening the file : " // trim(this%files(ifile)))
      end if

      nlines = 0
      this%found = .false.
      ierr = 0
      buffer = ''
      do while(index(buffer,'body').eq.0)

         ! First get some parameters
         read(iunit,'(a)',iostat=ierr) buffer

         ! Get nPoints and allocate arrays
         if (index(buffer,'gridPoints').ne.0) then
            read(buffer(index(buffer,'=')+1:),*) this%nPoints
            nlines=ceiling(this%nPoints/5.0_WP)

            allocate(this%input_data(this%nPoints,this%nvar_in))
            allocate(tmp(this%nPoints))
            allocate(this%Z(this%nPoints))
         end if

      end do

      ! Test
      if (nlines.eq.0) call die("[flameletLib readfile] missing gridPoints in flamemet file")

      ! Preset diffusivity to 1
      loop0: do var=1,this%nvar_in
         if (trim(this%input_name(var)).eq.'diffusivity') exit loop0
      end do loop0
      if (var.le.this%nvar_in) this%input_data(:,var) = 1.0_WP

      loop1:do while (ierr.eq.0)

         ! Read name of variable
         read(iunit,'(a)',iostat=ierr) buffer
         if (trim(buffer).eq.'trailer') exit loop1

         ! Read name of variable
         read(buffer,'(a)',iostat=ierr) varname
         index1=index(varname,' ')
         if (index1.ne.0) varname(index1:)=''

         ! Read the array
         line = ''
         do n=1,nlines
            read(iunit,'(a)',iostat=ierr) buffer
            line=trim(line)//adjustl(trim(buffer))
         end do

         ! Is it the coordinate Z?
         if (trim(varname).eq.'Z') then
            read(line,*) this%Z
         end if

         ! Is it part of the diffusivity?
         if (trim(varname).eq.'lambda') then
            read(line,*) tmp
            loop4: do var=1,this%nvar_in
               if (trim(this%input_name(var)).eq.'diffusivity') exit loop4
            end do loop4
            if (var.le.this%nvar_in) then
               this%input_data(:,var)=this%input_data(:,var)*tmp
               this%found(var)=.true.
            end if
         end if
         if (trim(varname).eq.'cp') then
            read(line,*) tmp
            loop5:do var=1,this%nvar_in
               if (trim(this%input_name(var)).eq.'diffusivity') exit loop5
            end do loop5
            if (var.le.this%nvar_in) then
               this%input_data(:,var)=this%input_data(:,var)/tmp
               this%found(var)=.true.
            end if
         end if

         ! Do we want that variable?
         loop3:do var=1,this%nvar_in
            if (trim(this%input_name(var)).eq.varname) then
               read(line,*) this%input_data(:,var)
               this%found(var) = .true.
               exit loop3
            end if
         end do loop3

      end do loop1

      ! Do we have all the variables?
      do var=1,this%nvar_in
         if (.not.this%found(var)) then
            print*,"Variable ",trim(this%input_name(var))," not found in flamelet file"
            stop
         end if
      end do

      ! Force 0 at Z=1 for chi
      if (this%combModel.eq.sfm) then
         this%input_data(this%nPoints,this%nvar_in)=0.0_WP
      end if

      ! Deallocate
      deallocate(tmp)
      nullify(tmp)
      close(iunit)
   end subroutine readfile

   subroutine cleanup(this)
      implicit none
      class(flameletLib), intent(inout) :: this

      deallocate(this%input_data)
      deallocate(this%Z)
      nullify(this%input_data)
      nullify(this%Z)
   end subroutine cleanup

   subroutine modify_varname(input_str,search_str,replace_str)
      character(len=*), intent(inout) :: input_str
      character(len=*), intent(in)    :: search_str,replace_str
      integer :: pos_sub,len_input,len_search
      len_input =len_trim(input_str)
      len_search=len_trim(search_str)
      pos_sub=index(input_str,search_str)
      if (pos_sub.gt.0) then
         input_str=trim(replace_str//input_str(pos_sub+len_search:len_input))
      end if
   end subroutine
end module flameletLib_class
