!> Artificial neural network class:
module ann_class
   use precision,    only: WP
   use string,       only: str_medium
   use config_class, only: config
   use mathtools,    only: ReLU
   implicit none
   private

   ! Expose type/constructor/methods
   public :: ann


   !> Hidden layer for the ANN
   type :: hidden_layer
      real(WP), dimension(:,:), allocatable :: weight                    !< Weight matrix
      real(WP), dimension(:),   allocatable :: bias                      !< Bias vector
   end type hidden_layer


   !> ANN object definition
   type :: ann

      ! Pointer to the config
      type(config), pointer :: cfg

      ! The name of the ann
      character(len=str_medium) :: name='UNNAMED_ANN'

      ! Data file name
      character(len=str_medium) :: filename 

      ! Input transformers
      logical :: has_inp_trn
      real(WP), dimension(:), allocatable :: inp_shift,inp_scale

      ! Input sub-array indices
      logical :: has_sub_ind
      integer, dimension(:), allocatable :: inp_sub_ind

      ! Encoder
      logical :: has_encoder
      real(WP), dimension(:,:), allocatable :: encoder_weight

      ! Hidden layers
      integer :: n_hid_lay
      type(hidden_layer), dimension(:), allocatable :: hid_lays

      ! Output layer
      real(WP), dimension(:,:), allocatable :: out_weight
      real(WP), dimension(:),   allocatable :: out_bias

      ! Output transformers
      logical :: has_out_trn
      real(WP), dimension(:), allocatable :: out_shift,out_scale

      ! Temporary array
      real(WP), dimension(:), allocatable :: tmp_arr                     !< Temporary array with the greatest size needed
      integer,  dimension(:), allocatable :: tmp_sz                      !< Stores the size of tmp_arr for all the layers if we had multiple tmp_arr

      ! Forward pass and encode
      procedure(forward_pass_interface), pointer :: forward_pass=>NULL() !< Take input to forward pass and get the output
      procedure(encode_interface),       pointer :: encode      =>NULL() !< Map the input to the latent variables

   contains
      generic            :: read_vector=>read_ivector,read_rvector       !< Parallel read a vector
      procedure, private :: read_ivector                                 !< Parallel read an integer vector
      procedure, private :: read_rvector                                 !< Parallel read a real vector
      generic            :: read_matrix=>read_imatrix,read_rmatrix       !< Parallel read a matrix
      procedure, private :: read_imatrix                                 !< Parallel read an integer matrix
      procedure, private :: read_rmatrix                                 !< Parallel read a real matrix
      procedure          :: print=>ann_print                             !< Print information
   end type ann


   !> Declare ANN constructor
   interface ann
      procedure constructor
   end interface ann


   !> Interface for encode and forward pass subroutines
   interface
      subroutine forward_pass_interface(this,input,output)
         use precision, only: WP
         import ann
         class(ann), intent(inout)           :: this
         real(WP), dimension(:), intent(in)  :: input
         real(WP), dimension(:), intent(out) :: output
      end subroutine forward_pass_interface
      subroutine encode_interface(this,input,output)
         use precision, only: WP
         import ann
         class(ann), intent(inout) :: this
         real(WP), dimension(:), intent(in)  :: input
         real(WP), dimension(:), intent(out) :: output
      end subroutine encode_interface
   end interface


contains


   !> Default constructor for ANN
   function constructor(cfg,fdata,name) result(self)
      use messager, only: die
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      type(ann)                              :: self
      class(config), target, intent(in)      :: cfg
      character(len=*), intent(in)           :: fdata
      character(len=*), intent(in), optional :: name
      integer                                :: ihid
      integer                                :: ierr
      type(MPI_File)                         :: ifile
      type(MPI_Status)                       :: status

      ! Point to pgrid object
      self%cfg=>cfg

      ! Set the name for the ann object
      if (present(name)) self%name=trim(adjustl(name))

      ! Set the file name
      self%filename=trim(adjustl(fdata))

      ! Open the file
      call MPI_FILE_OPEN(self%cfg%comm,trim(self%filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
      if (ierr.ne.0) call die('[ann constructor] Problem encountered while reading file: '//trim(self%filename))

      ! Read the input transformers
      call MPI_FILE_READ_ALL(ifile,self%has_inp_trn,1,MPI_LOGICAL,status,ierr)
      if (self%has_inp_trn) then
         ! Shift vector
         call self%read_vector(vector=self%inp_shift,ifile=ifile,status=status,ierr=ierr)
         ! Scale vector
         call self%read_vector(vector=self%inp_scale,ifile=ifile,status=status,ierr=ierr)
      end if

      ! Read the input sub-array indices
      call MPI_FILE_READ_ALL(ifile,self%has_sub_ind,1,MPI_LOGICAL,status,ierr)
      if (self%has_sub_ind) then
         call self%read_vector(vector=self%inp_sub_ind,ifile=ifile,status=status,ierr=ierr)
      end if

      ! Read the encoder
      call MPI_FILE_READ_ALL(ifile,self%has_encoder,1,MPI_LOGICAL,status,ierr)
      if (self%has_encoder) then
            ! Read the weight matrix
            call self%read_matrix(matrix=self%encoder_weight,ifile=ifile,status=status,ierr=ierr)
            ! Point to the correspoding encode subsoutine
            if (self%has_inp_trn) then
               self%encode=>encode_T
            else
               self%encode=>encode_F
            end if
      else
         ! Point to the correspoding encode subsoutine
         self%encode=>NULL()
      end if

      ! Read the hidden layers
      call MPI_FILE_READ_ALL(ifile,self%n_hid_lay,1,MPI_INTEGER,status,ierr)
      allocate(self%hid_lays(self%n_hid_lay),self%tmp_sz(self%n_hid_lay+1))
      do ihid=1,self%n_hid_lay
         ! Read the bias vector
         call self%read_vector(vector=self%hid_lays(ihid)%bias,ifile=ifile,status=status,ierr=ierr)
         ! Read the weight matrix
         call self%read_matrix(matrix=self%hid_lays(ihid)%weight,ifile=ifile,status=status,ierr=ierr)
         ! Update the temporary array length
         self%tmp_sz(ihid)=size(self%hid_lays(ihid)%weight,dim=1)
      end do

      ! Read the output layer
      call self%read_vector(vector=self%out_bias,ifile=ifile,status=status,ierr=ierr)
      call self%read_matrix(matrix=self%out_weight,ifile=ifile,status=status,ierr=ierr)
      
      ! Update the temporary array length and allocate
      self%tmp_sz(self%n_hid_lay+1)=size(self%out_weight,dim=1)
      allocate(self%tmp_arr(maxval(self%tmp_sz))); self%tmp_arr=0.0_WP

      ! Read the output transformers
      call MPI_FILE_READ_ALL(ifile,self%has_out_trn,1,MPI_LOGICAL,status,ierr)
      if (self%has_out_trn) then
         ! Shift vector
         call self%read_vector(vector=self%out_shift,ifile=ifile,status=status,ierr=ierr)
         ! Scale vector
         call self%read_vector(vector=self%out_scale,ifile=ifile,status=status,ierr=ierr)
         ! Point to the corresponding forward pass subroutine
         self%forward_pass=>forward_pass_T
      else
         ! Point to the corresponding forward pass subroutine
         self%forward_pass=>forward_pass_F
      end if

      ! Close the file
      call MPI_FILE_CLOSE(ifile,ierr)
   end function constructor


   !> Parallel read an integer vector
   subroutine read_ivector(this,vector,ifile,status,ierr)
      use parallel, only: info_mpiio
      use mpi_f08
      implicit none
      class(ann), intent(inout)                         :: this
      integer, dimension(:), allocatable, intent(out)   :: vector
      type(MPI_File),                     intent(inout) :: ifile
      type(MPI_Status),                   intent(inout) :: status
      integer,                            intent(inout) :: ierr
      integer :: n
      call MPI_FILE_READ_ALL(ifile,n,1,MPI_INTEGER,status,ierr)
      allocate(vector(n))
      call MPI_FILE_READ_ALL(ifile,vector,n,MPI_INTEGER,status,ierr)
   end subroutine read_ivector


   !> Parallel read a real vector
   subroutine read_rvector(this,vector,ifile,status,ierr)
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      class(ann), intent(inout)                          :: this
      real(WP), dimension(:), allocatable, intent(out)   :: vector
      type(MPI_File),                      intent(inout) :: ifile
      type(MPI_Status),                    intent(inout) :: status
      integer,                             intent(inout) :: ierr
      integer :: n
      call MPI_FILE_READ_ALL(ifile,n,1,MPI_INTEGER,status,ierr)
      allocate(vector(n))
      call MPI_FILE_READ_ALL(ifile,vector,n,MPI_REAL_WP,status,ierr)
   end subroutine read_rvector


   !> Parallel read an integer matrix
   subroutine read_imatrix(this,matrix,ifile,status,ierr)
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      class(ann), intent(inout)                           :: this
      integer, dimension(:,:), allocatable, intent(out)   :: matrix
      type(MPI_File),                       intent(inout) :: ifile
      type(MPI_Status),                     intent(inout) :: status
      integer,                              intent(inout) :: ierr
      integer :: n_row,n_col
      call MPI_FILE_READ_ALL(ifile,n_row,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,n_col,1,MPI_INTEGER,status,ierr)
      allocate(matrix(n_row,n_col))
      call MPI_FILE_READ_ALL(ifile,matrix,n_row*n_col,MPI_INTEGER,status,ierr)
   end subroutine read_imatrix


   !> Parallel read a real matrix
   subroutine read_rmatrix(this,matrix,ifile,status,ierr)
      use parallel, only: info_mpiio,MPI_REAL_WP
      use mpi_f08
      implicit none
      class(ann), intent(inout)                            :: this
      real(WP), dimension(:,:), allocatable, intent(out)   :: matrix
      type(MPI_File),                        intent(inout) :: ifile
      type(MPI_Status),                      intent(inout) :: status
      integer,                               intent(inout) :: ierr
      integer :: n_row,n_col
      call MPI_FILE_READ_ALL(ifile,n_row,1,MPI_INTEGER,status,ierr)
      call MPI_FILE_READ_ALL(ifile,n_col,1,MPI_INTEGER,status,ierr)
      allocate(matrix(n_row,n_col))
      call MPI_FILE_READ_ALL(ifile,matrix,n_row*n_col,MPI_REAL_WP,status,ierr)
   end subroutine read_rmatrix


   !> Encoder without input transform
   subroutine encode_F(this,input,output)
      implicit none
      class(ann), intent(inout)           :: this
      real(WP), dimension(:), intent(in)  :: input
      real(WP), dimension(:), intent(out) :: output
      output=matmul(input,this%encoder_weight)
   end subroutine encode_F


   !> Encoder without input transform
   subroutine encode_T(this,input,output)
      implicit none
      class(ann), intent(inout)           :: this
      real(WP), dimension(:), intent(in)  :: input
      real(WP), dimension(:), intent(out) :: output
      output=matmul((input-this%inp_shift)/this%inp_scale,this%encoder_weight)
   end subroutine encode_T


   !> Forward pass without output transform
   subroutine forward_pass_F(this,input,output)
      implicit none
      class(ann), intent(inout)           :: this
      real(WP), dimension(:), intent(in)  :: input
      real(WP), dimension(:), intent(out) :: output
      integer :: ihid
      this%tmp_arr(1:this%tmp_sz(1))=input
      ! Hidden layer calculations
      do ihid=1,this%n_hid_lay
         this%tmp_arr(1:this%tmp_sz(ihid+1))=ReLU(matmul(this%tmp_arr(1:this%tmp_sz(ihid)),this%hid_lays(ihid)%weight)+this%hid_lays(ihid)%bias)
      end do
      ! Output layer calculations
      output=matmul(this%tmp_arr(1:this%tmp_sz(this%n_hid_lay+1)),this%out_weight)+this%out_bias
   end subroutine forward_pass_F


   !> Forward pass with output transform
   subroutine forward_pass_T(this,input,output)
      implicit none
      class(ann), intent(inout)           :: this
      real(WP), dimension(:), intent(in)  :: input
      real(WP), dimension(:), intent(out) :: output
      integer :: ihid
      this%tmp_arr(1:this%tmp_sz(1))=input
      ! Hidden layer calculations
      do ihid=1,this%n_hid_lay
         this%tmp_arr(1:this%tmp_sz(ihid+1))=ReLU(matmul(this%tmp_arr(1:this%tmp_sz(ihid)),this%hid_lays(ihid)%weight)+this%hid_lays(ihid)%bias)
      end do
      ! Output layer calculations
      output=(matmul(this%tmp_arr(1:this%tmp_sz(this%n_hid_lay+1)),this%out_weight)+this%out_bias)*this%out_scale+this%out_shift
   end subroutine forward_pass_T


   !> Print out info for ann object
   subroutine ann_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(ann), intent(in) :: this
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("ANN object [",a,"] was read from file [",a,"].")') trim(this%name),trim(this%filename)
      end if
   end subroutine ann_print


end module ann_class