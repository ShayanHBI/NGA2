!> Flamelet class:
module flamelet_class
   use precision,       only: WP
   use string,          only: str_medium
   use config_class,    only: config
   use chemtable_class, only: chemtable
   use messager,        only: die
   implicit none
   private

   ! Expose type/constructor/methods
   public :: flamelet

   ! List of known available flamelet models
   integer, parameter, public :: sfm=1                       !< Steady Flamelet Model (SFM)
   integer, parameter, public :: bs =1                       !< Burke–Schumann

   !> flamelet object definition
   type :: flamelet

      ! This is our config
      class(config), pointer :: cfg                          !< This is the config the solver is build for

      ! chemtable
      type(chemtable) :: chmtbl

      ! This is the name of the flamelet object
      character(len=str_medium) :: name='UNNAMED_flamelet'   !< Name (default=UNNAMED_flamelet)

      ! flamelet model
      integer :: flmModel

      ! Field variables
      real(WP), dimension(:,:,:), allocatable :: Zvar !< Favre-spatially-filtered variance of the mixture fraction
      real(WP), dimension(:,:,:), allocatable :: chi  !< Favre-spatially-filtered scalar dissipation rate
   contains
      procedure :: compute_Zvar
      procedure :: compute_chi
      procedure :: print=>flamelet_print
   end type flamelet


   !> Declare flamelet constructor
   interface flamelet
      procedure constructor
   end interface flamelet

contains


   !> Default constructor for flamelet object
   function constructor(cfg,flmModel,tablefile,name) result(self)
      implicit none
      type(flamelet) :: self
      class(config), target, intent(in) :: cfg
      integer, intent(in) :: flmModel
      character(len=*), intent(in) :: tablefile
      character(len=*), optional :: name
      integer :: i,j,k

      ! Set the name for the solver
      if (present(name)) self%name=trim(adjustl(name))

      ! Point to config object
      self%cfg=>cfg

      ! Build the chmistry table
      self%chmtbl=chemtable(cfg=cfg,fdata=tablefile)

      ! Set the flamelet model
      self%flmModel=flmModel
      if (self%chmtbl%comb_model.ne.self%flmModel) call die('[flamelet constructor] Chemtable model is incompatible with the flamelet model')
      select case (self%flmModel)
       case (sfm)
       case default
         call die('[flamelet constructor] Unknown flamelet model')
      end select

      ! Allocate arrays
      allocate(self%Zvar(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%Zvar=0.0_WP
      allocate(self%chi (self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_)); self%chi =0.0_WP
   end function constructor


   subroutine compute_Zvar(this,delta,ZgradMagSq,Z)
      class(flamelet), intent(inout) :: this
      real(WP), dimension(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_), intent(in) :: delta
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in) :: ZgradMagSq,Z
      real(WP), parameter :: Cz=0.15_WP**2/1.0_WP
      
      this%Zvar=Cz*delta**2*ZgradMagSq
      ! Clip the computed Variance
      this%Zvar=max(0.0_WP,min(Z*(1.0_WP-Z),this%Zvar))
   end subroutine compute_Zvar


   subroutine compute_chi(this,mueff,rho,ZgradMagSq)
      class(flamelet), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in) :: mueff,rho,ZgradMagSq
      real(WP), parameter :: Cchi=2.0_WP

      this%chi=Cchi*mueff/rho*ZgradMagSq
   end subroutine compute_chi


   !> Print out info for flamelet model
   subroutine flamelet_print(this)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      class(flamelet), intent(in) :: this
      
      ! Output
      if (this%cfg%amRoot) then
         write(output_unit,'("Flamelet model [",a,"] for config [",a,"]")') trim(this%name),trim(this%cfg%name)
      end if
      
   end subroutine flamelet_print

   
end module flamelet_class
