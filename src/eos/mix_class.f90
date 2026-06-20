!> Abstract mixture thermodynamic model.
module mix_class
   use precision, only: WP
   use eos_class, only: eos
   implicit none
   private

   public :: mix

   type, abstract :: mix
      integer :: ns=0
   contains
      procedure(mix_rho_e_iface), deferred :: get_p_from_rho_e
      procedure(mix_p_rho_iface), deferred :: get_T_from_p_rho
      procedure(mix_p_rho_iface), deferred :: get_c_from_p_rho
      procedure(mix_p_rho_iface), deferred :: get_e_from_p_rho
      procedure(mix_p_T_iface  ), deferred :: get_e_from_p_T
      procedure(mix_rho_T_iface), deferred :: get_p_from_rho_T
      procedure(mix_p_T_iface  ), deferred :: get_rho_from_p_T
      procedure(mix_p_T_iface  ), deferred :: get_h_from_p_T
      procedure(mix_p_T_iface  ), deferred :: get_s_from_p_T
      procedure(mix_rho_e_iface), deferred :: get_gruneisen_from_rho_e
      procedure(mix_p_rho_iface), deferred :: get_rhoe_from_p_rho
      procedure(mix_p_T_iface  ), deferred :: get_rhoe_from_p_T
      procedure(mix_p_T_iface  ), deferred :: get_g_from_p_T
      procedure(mix_rho_T_iface), deferred :: get_drhodT_const_p_from_rho_T
      procedure(mix_rho_T_iface), deferred :: get_drhodp_const_T_from_rho_T
   end type mix

   abstract interface
      real(WP) function mix_rho_e_iface(this,rho,e,y)
         import :: WP,mix
         class(mix), intent(in) :: this
         real(WP), intent(in) :: rho,e
         real(WP), dimension(:), intent(in) :: y
      end function mix_rho_e_iface

      real(WP) function mix_p_rho_iface(this,p,rho,y)
         import :: WP,mix
         class(mix), intent(in) :: this
         real(WP), intent(in) :: p,rho
         real(WP), dimension(:), intent(in) :: y
      end function mix_p_rho_iface

      real(WP) function mix_p_T_iface(this,p,T,y)
         import :: WP,mix
         class(mix), intent(in) :: this
         real(WP), intent(in) :: p,T
         real(WP), dimension(:), intent(in) :: y
      end function mix_p_T_iface

      real(WP) function mix_rho_T_iface(this,rho,T,y)
         import :: WP,mix
         class(mix), intent(in) :: this
         real(WP), intent(in) :: rho,T
         real(WP), dimension(:), intent(in) :: y
      end function mix_rho_T_iface
   end interface

end module mix_class
