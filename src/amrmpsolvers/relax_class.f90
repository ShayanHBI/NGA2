!> Abstract base class for two-phase compressible relaxation models.
module relax_class
   use precision, only: WP
   implicit none
   private

   public :: relax

   !> Abstract relaxation model type
   type, abstract :: relax

   contains
      procedure(relax_apply_iface), deferred :: apply
   end type relax

   !> Interface for the deferred apply procedure
   abstract interface
      subroutine relax_apply_iface(this,VF,Q,Pjump)
         import :: relax,WP
         class(relax), intent(inout) :: this
         real(WP), intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
         real(WP), intent(in)    :: Pjump
      end subroutine relax_apply_iface
   end interface

end module relax_class
