!> Flamelet class:
module flamelet_class
   use precision,       only: WP
   use string,          only: str_medium
   use config_class,    only: config
   use chemtable_class, only: chemtable
   use messager,        only: die
   implicit none
   private

   real(WP), parameter :: Cs=0.15_WP**2/1.0_WP

   ! Expose type/constructor/methods
   public :: flamelet

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
      real(WP), dimension(:,:,:), allocatable :: Zvar        !< Favre-spatially-filtered variance of the mixture fraction
      real(WP), dimension(:,:,:), allocatable :: chi         !< Favre-spatially-filtered scalar dissipation rate

   contains
      procedure :: get_Zvar
      procedure :: get_chi
      procedure :: print=>flamelet_print
   end type flamelet


   !> Declare flamelet constructor
   interface flamelet
      procedure constructor
   end interface flamelet


contains


   !> Default constructor for flamelet object
   function constructor(cfg,flmModel,tablefile,name) result(self)
      use flameletLib_class, only: sfm
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


   subroutine get_Zvar(this,delta,ZgradMagSq,Z)
      use sgsmodel_class, only: Cs_diff
      class(flamelet), intent(inout) :: this
      real(WP), dimension(this%cfg%imin_:this%cfg%imax_,this%cfg%jmin_:this%cfg%jmax_,this%cfg%kmin_:this%cfg%kmax_), intent(in) :: delta
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in) :: ZgradMagSq,Z
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%Zvar(i,j,k)=Cs_diff*delta(i,j,k)**2*ZgradMagSq(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(this%Zvar)
      ! Clip the computed Variance
      this%Zvar=max(0.0_WP,min(Z*(1.0_WP-Z),this%Zvar))
   end subroutine get_Zvar


   subroutine get_chi(this,mueff,rho,ZgradMagSq)
      class(flamelet), intent(inout) :: this
      real(WP), dimension(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_), intent(in) :: mueff,rho,ZgradMagSq
      real(WP), parameter :: Cchi=2.0_WP
      integer :: i,j,k
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%chi(i,j,k)=Cchi*mueff(i,j,k)/rho(i,j,k)*ZgradMagSq(i,j,k)
            end do
         end do
      end do
      ! Sync it
      call this%cfg%sync(this%chi)
   end subroutine get_chi


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
