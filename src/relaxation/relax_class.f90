!> Abstract base class for two-phase compressible relaxation models.
module relax_class
   use precision, only: WP
   implicit none
   private

   public :: relax,dbg_i,dbg_j

   !> Debug: current cell index, set by the caller right before invoking relaxation,
   !> so prints deep inside relax/relax_pTg implementations can be filtered to one cell
   integer :: dbg_i=-1,dbg_j=-1

   !> Abstract relaxation model type
   type, abstract :: relax

   contains
      procedure(relax_iface), deferred :: relax_p
      procedure(relax_iface), deferred :: relax_pT
      procedure(relax_iface), deferred :: relax_pTg
   end type relax

   !> Shared interface for liquid-gas relaxation
   abstract interface
      subroutine relax_iface(this,VF,Q,Pjump)
         import :: relax,WP
         class(relax), intent(inout)      :: this
         real(WP),               intent(inout) :: VF
         real(WP), dimension(:), intent(inout) :: Q
         real(WP),               intent(in)    :: Pjump
      end subroutine relax_iface
   end interface

end module relax_class
