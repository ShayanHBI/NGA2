!> Abstract base class for two-phase compressible relaxation models.
module relax_class
   use precision, only: WP
   implicit none
   private

   public :: relax

   !> Abstract relaxation model type
   type, abstract :: relax

   contains
      procedure(relax_iface), deferred :: relax_p
      procedure(relax_iface), deferred :: relax_pT
      procedure(relax_iface), deferred :: relax_pTg
   end type relax

   !> Shared interface for all three relaxation steps.
   !> relax_pT calls relax_p internally; relax_pTg calls relax_pT internally.
   !> The caller selects the level: relax_p for P only, relax_pT for P+T, relax_pTg for P+T+g.
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
