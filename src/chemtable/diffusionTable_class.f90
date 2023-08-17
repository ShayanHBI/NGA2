module diffusionTable_class
   use precision,         only: WP
   use string,            only: str_medium,str_long
   use flameletlib_class, only: flameletlib,sfm
   use messager,          only: die
   implicit none
   private

   ! Expose type/constructor/methods
   public :: diffusionTable

   !> Flamelet table object definition
   type :: diffusionTable

      ! Flamelet library
      type(flameletlib), pointer :: flmlib

      ! Stoichiometric mixture fraction
      real(WP) :: Zst

      ! Mesh stretching for the third direction
      character(len=str_long) :: Z3_scale

      ! Table parameters
      integer :: nZMean, nZVar, n3
      real(WP), dimension(:), allocatable :: ZMean, ZVar, Z3

      ! Beta pdf
      real(WP), dimension(:), pointer :: pdf

      ! Variables to be wrtitten into the table
      integer :: nvar_out
      character(len=str_medium), dimension(:), allocatable :: output_name
      real(WP), dimension(:,:,:,:), allocatable :: output_data
      integer,  dimension(:,:,:),   allocatable :: mask

      ! Variable after convolution
      real(WP), dimension(:,:,:,:), pointer :: postconv

      ! Name conversion
      ! character(len=str_medium), dimension(:), allocatable :: conversion

      ! Binary file name of the table
      character(len=str_medium) :: filename

   contains
      procedure :: create_Zmean
      procedure :: create_Zvar
      procedure :: create_Z3
      procedure :: create_beta_pdf
      procedure :: convolute
      procedure :: convert_names
      procedure :: setup
      procedure :: stats
      procedure :: write => diffusionTable_write
   end type diffusionTable

   !> Declare flameletlib constructor
   interface diffusionTable
      procedure constructor
   end interface diffusionTable

contains


   !> Default constructor for diffusionTable object
   function constructor(flmlib) result(self)
      use param, only: param_read
      implicit none
      type(diffusionTable) :: self
      class(flameletlib), target, intent(in) :: flmlib
  
      ! Point to the flmlib object
      self%flmlib=>flmlib

      ! Read the dimension of the final table
      call param_read('Number of points for mean Z',self%nZMean)
      call param_read('Number of points for variance of Z',self%nZVar)
      call param_read('Number of points for third direction',self%n3)

      call param_read('Stoichiometric mixture fraction',self%Zst)

      ! call param_read('Name conversion',self%conversion)

      ! Default name as in FlameMaster
      self%nvar_out=self%flmlib%nvar_in
      allocate(self%output_name(self%nvar_out))
      self%output_name=self%flmlib%input_name(1:self%flmlib%nvar_in)

      call param_read('Table filename',self%filename)

      call param_read('Scale for third direction',self%Z3_scale)

      ! Dimensions of the final table
      allocate(self%ZMean(self%nZMean))
      allocate(self%ZVar(self%nZVar))
      allocate(self%Z3(self%n3))

      ! Create the first two directions of the table
      call self%create_Zmean
      call self%create_Zvar

      ! Allocate arrays
      allocate(self%postconv(self%flmlib%nvar_in,self%nZMean,self%nZVar,self%flmlib%nfiles))
   end function constructor

   subroutine create_Zmean(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer  :: zm,zcut
      real(WP) :: m11,m12,m21,m22,r1,r2,delta
      real(WP) :: a,b,c
      real(WP) :: dz

      ! Create such mesh only if Zst sufficiently small
      if (this%Zst.gt.0.30_WP) then
         dz=1.0_WP/real(this%nZMean-1,WP)
         do zm=1,this%nZMean
            this%ZMean(zm)=real(zm-1,WP)*dz
         end do
         return
      end if

      ! Linear mesh for [0,Zst]
      zcut=this%nZMean/3+1
      dz=this%Zst/real(zcut-1,WP)
      do zm=1,zcut
         this%ZMean(zm)=real(zm-1,WP)*dz
      end do

      ! Mesh with linear growth to reach Z=1
      m11=real(this%nZMean**2-zcut**2,WP)
      m12=real(this%nZMean-zcut,WP)
      m21=real(2*zcut+1,WP)
      m22=real(1,WP)
      r1=1.0_WP-this%Zst
      r2=dz
      delta=m11*m22-m12*m21
      a=(+m22*r1-m12*r2)/delta
      b=(-m21*r1+m11*r2)/delta
      c=this%Zst-a*zcut**2-b*zcut

      do zm=zcut+1,this%nZMean
         this%ZMean(zm)=a*real(zm,WP)**2+b*real(zm,WP)+c
      end do
   end subroutine create_Zmean

   subroutine create_Zvar(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer :: izv

      do izv=1,this%nZVar
         this%ZVar(izv)=0.25_WP*(real(izv-1,WP)/real(this%nZVar-1,WP))**2
      end do
   end subroutine create_Zvar

   subroutine create_Z3(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer  :: i3
      real(WP) :: min3,max3

      ! Find min and max
      max3=maxval(this%postconv(this%flmlib%nvar_in,:,:,:))
      min3=minval(this%postconv(this%flmlib%nvar_in,:,:,:))

      ! Linear or log progression
      select case(trim(this%Z3_scale))
       case ('lin')
         do i3=1,this%n3
            this%Z3(i3)=min3+real(i3-1,WP)*(max3-min3)/real(this%n3-1,WP)
         end do
       case ('log')
         do i3=1,this%n3
            this%Z3(i3)=min3*(max3/min3)**(real(i3-1,WP)/real(this%n3-1,WP))
         end do
       case default
         call die("[create_Z3] Unknown Scale for third direction")
      end select
   end subroutine create_Z3

   subroutine create_beta_pdf(this,zm,zv)
      use mathtools, only: gammaln
      implicit none
      class(diffusionTable), intent(inout) :: this
      real(WP), intent(in) :: zm,zv
      real(WP) :: a,b,factor,tmp,dz
      integer  :: index1,n

      this%pdf=0.0_WP

      ! Zero mean : delta at Z=0
      if (zm.le.1.0E-10_WP) then
         this%pdf(1)=1.0_WP
         return
      end if

      ! Max mean : delta at Z=1
      if (zm.ge.1.0_WP-1.0E-10_WP) then
         this%pdf(this%flmlib%nPoints)=1.0_WP
         return
      end if

      ! Zero variance : delta at Z=zm
      if (zv.le.1.0E-10_WP) then
         index1=1
         do while (this%flmlib%Z(index1).lt.zm)
            index1=index1+1
         end do
         this%pdf(index1-1)=(this%flmlib%Z(index1)-zm)/(this%flmlib%Z(index1)-this%flmlib%Z(index1-1))
         this%pdf(index1)=(zm-this%flmlib%Z(index1-1))/(this%flmlib%Z(index1)-this%flmlib%Z(index1-1))
         return
      end if

      ! Impossible cases => two delta at 0 and 1
      if (zv.ge.zm*(1.0_WP-zm)) then
         this%pdf(1)=1.0_WP-zm
         this%pdf(this%flmlib%nPoints)=zm
         return
      end if

      a=zm*(zm*(1.0_WP-zm)/zv-1.0_WP)
      b=a/zm-a
      factor=gammaln(a+b)-gammaln(a)-gammaln(b)

      ! Left BC : explicit integration
      dz=0.5_WP*(this%flmlib%Z(2)-this%flmlib%Z(1))
      tmp=a*log(dz)+factor
      this%pdf(1)=exp(tmp)/a
      ! Right BC : explicit integration
      dz=0.5_WP*(this%flmlib%Z(this%flmlib%nPoints)-this%flmlib%Z(this%flmlib%nPoints-1))
      tmp=b*log(dz)+factor
      this%pdf(this%flmlib%nPoints)=exp(tmp)/b
      ! Other Points
      do n=2,this%flmlib%nPoints-1
         dz=0.5_WP*(this%flmlib%Z(n+1)-this%flmlib%Z(n-1))
         tmp=(a-1.0_WP)*log(this%flmlib%Z(n))+(b-1.0_WP)*log(1.0_WP-this%flmlib%Z(n))
         tmp=tmp+factor
         this%pdf(n)=exp(tmp)*dz
      end do

      ! Normalize the pdf
      this%pdf=this%pdf/sum(this%pdf)
   end subroutine create_beta_pdf

   subroutine convolute(this,ifile)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer, intent(in) :: ifile
      integer :: izm, izv, k, var
      real(WP) :: meanVal

      ! Prepare the convolution
      allocate(this%pdf(this%flmlib%nPoints))
      
      ! Convolution
      do izv=1,this%nZVar
         do izm=1,this%nZMean
            call this%create_beta_pdf(this%ZMean(izm),this%ZVar(izv))

            do var=1,this%flmlib%nvar_in
               meanVal=0.0_WP

               if (trim(this%flmlib%input_name(var)).eq.'density') then
                  do k=1,this%flmlib%nPoints
                     meanVal=meanVal+this%pdf(k)/this%flmlib%input_data(k,var)
                  end do
                  this%postconv(var,izm,izv,ifile)=1.0_WP/meanVal
               else
                  do k=1,this%flmlib%nPoints
                     meanVal=meanVal+this%pdf(k)*this%flmlib%input_data(k,var)
                  end do
                  this%postconv(var,izm,izv,ifile)=meanVal
               end if
            end do
         end do
      end do

      ! Finish the convolution
      deallocate(this%pdf)
      nullify(this%pdf)
   end subroutine convolute

   subroutine convert_names(this)
      use flameletLib_class, only:modify_varname
      implicit none
      class(diffusionTable), intent(inout) :: this
      character(len=str_medium) :: varname
      integer :: i,n,var

      ! ! Get the number of name conversions
      !    call parser_getsize('Name conversion',n)
      !    if (mod(n,3).ne.0) stop "diffusion_table_convert_names: Problem in the definition of conversion names"
         
      !    Allocate array and read
      !    allocate(conversion(n))
      !    call parser_read('Name conversion',conversion)

      ! ! Convert the names
      ! n = n / 3
      ! do i=1,n
      !    varname = trim(this%conversion((i-1)*3+3))
      !    loop1:do var=1,this%flmlib%nvar_in
      !       if (trim(this%flmlib%input_name(var)).eq.trim(varname)) exit loop1
      !    end do loop1
      !    if (var.eq.this%flmlib%nvar_in+1) then
      !       call die("[diffusionTable convert_names] Unknown variable name : " // varname)
      !    end if
      !    this%output_name(var) = trim(this%conversion((i-1)*3+1))
      ! end do

      ! Modify the mass fraction names
      do var=1,this%nvar_out
         call modify_varname(this%output_name(var),'massfraction-','Y_')
      end do
   end subroutine convert_names

   subroutine setup(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer  :: izm,izv,i3,var
      integer  :: file, file_up, file_down
      real(WP) :: err, err_up, err_down
      real(WP) :: alpha_up,alpha_down
      real(WP), dimension(:), pointer :: tmp

      ! Convert to a chi=0 flamelet
      if (this%flmlib%combModel.eq.sfm) then
         file=minloc(maxval(maxval(this%postconv(this%flmlib%nvar_in,:,:,:),dim=1),dim=1),dim=1)
         print*,''
         print*,'Flamelet #',file,'used as chi=0 flamelet'
         this%postconv(this%flmlib%nvar_in,:,:,file) = 0.0_WP
      end if

      ! Allocate final table
      allocate(this%output_data(this%nZMean,this%nZVar,this%n3,this%nvar_out))
      allocate(this%mask(this%nZMean,this%nZVar,this%n3))
      this%mask=0

      ! Create mesh in thrid direction
      call this%create_Z3

      ! Loop over the three mapping directions
      do i3=1,this%n3
         do izv=1,this%nZVar
            do izm=1,this%nZMean

               tmp=>this%postconv(this%flmlib%nvar_in,izm,izv,:)

               ! Find the two files right above and right below
               err_up  =+huge(1.0_WP)
               err_down=-huge(1.0_WP)

               file_up  =0
               file_down=0

               do file=1,this%flmlib%nfiles
                  err=tmp(file)-this%Z3(i3)

                  if ((err.ge.0.0_WP).and.(err.le.err_up)) then
                     file_up=file
                     err_up=err
                  end if
                  if ((err.le.0.0_WP).and.(err.ge.err_down)) then
                     file_down=file
                     err_down=err
                  end if
               end do

               ! Interpolate
               if (file_up.eq.0.or.file_down.eq.0) then
                  if (file_up.eq.0) then
                     alpha_up  =0.0_WP
                     alpha_down=1.0_WP
                     file_up=1
                  end if
                  if (file_down.eq.0) then
                     alpha_up  =1.0_WP
                     alpha_down=0.0_WP
                     file_down=1
                  end if
                  ! Mask it
                  this%mask(izm,izv,i3)=1
               else
                  if (file_up.eq.file_down) then
                     alpha_up  =1.0_WP
                     alpha_down=0.0_WP
                  else
                     alpha_up  =(this%Z3(i3)-tmp(file_down))/(tmp(file_up)-tmp(file_down))
                     alpha_down=(tmp(file_up)-this%Z3(i3))  /(tmp(file_up)-tmp(file_down))
                  end if
               end if

               do var=1,this%nvar_out
                  if (trim(this%flmlib%input_name(var)).eq.'density') then
                     this%output_data(izm,izv,i3,var)= 1.0_WP/(&
                        alpha_up  /this%postconv(var,izm,izv,file_up)+&
                        alpha_down/this%postconv(var,izm,izv,file_down))
                  else
                     this%output_data(izm,izv,i3,var)=&
                        alpha_up  *this%postconv(var,izm,izv,file_up)+&
                        alpha_down*this%postconv(var,izm,izv,file_down)
                  end if
               end do

            end do
         end do
      end do
   end subroutine setup

   subroutine stats(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer :: var

      print*,''

      ! Min and Max of the coordinates
      print*, '** Coordinates of the table **'
      write(*,10) 'Variable       ','Min       ','Max       '
      write(*,10) '---------------','---------------','---------------'
      write(*,11) 'ZMEAN       ', minval(this%ZMean), maxval(this%ZMean)
      write(*,11) 'ZVAR        ', minval(this%ZVar), maxval(this%ZVar)
      write(*,11) this%flmlib%input_name(this%flmlib%nvar_in), minval(this%Z3), maxval(this%Z3)
      print*,''
      
      ! Min and Max of all the mapped quantities
      print*, '** Mapped quantities **'
      write(*,10) 'Variable       ','Min       ','Max       '
      write(*,10) '---------------','---------------','---------------'
      do var=1,this%nvar_out
         write(*,11) this%output_name(var),minval(this%output_data(:,:,:,var)),maxval(this%output_data(:,:,:,var))
      end do
      print*,''

10    format (A15,'  ',A12,'  ',A12)
11    format (A15,'  ',ES12.4,'  ',ES12.4)
   end subroutine stats

   subroutine diffusionTable_write(this)
      implicit none
      class(diffusionTable), intent(inout) :: this
      integer :: ierr,var,iunit

      ! Open the data file
      open(newunit=iunit,file=trim(this%filename),form='unformatted',status='replace',access='stream',iostat=ierr)
      ! Write sizes
      write(iunit) this%nZMean
      write(iunit) this%nZVar
      write(iunit) this%n3
      write(iunit) this%nvar_out
      ! Write the axis coordinates
      write(iunit) this%ZMean
      write(iunit) this%ZVar
      write(iunit) this%Z3
      ! Masks
      ! write(iunit) this%mask
      ! Write additional stuff
      write(iunit) this%flmlib%combModel
      ! Write variable names
      do var=1,this%nvar_out
         write(iunit) this%output_name(var)
      end do
      ! Write data field
      do var=1,this%nvar_out
         write(iunit) this%output_data(:,:,:,var)
      end do
      ! Close the data file
      close(iunit)
   end subroutine diffusionTable_write
end module diffusionTable_class
