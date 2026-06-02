!> Abstract pure-substance equation-of-state base class.
module eos_class
   use precision, only: WP
   use string,    only: str_medium
   implicit none
   private

   public :: eos

   type, abstract :: eos
      character(len=str_medium) :: name = 'UNNAMED_EOS'
   contains
      procedure(eos_rho_e_iface), deferred :: get_p_from_rho_e
      procedure(eos_p_rho_iface), deferred :: get_T_from_p_rho
      procedure(eos_p_rho_iface), deferred :: get_c_from_p_rho
      procedure(eos_p_rho_iface), deferred :: get_e_from_p_rho
      procedure(eos_p_T_iface  ), deferred :: get_e_from_p_T
      procedure(eos_rho_T_iface), deferred :: get_p_from_rho_T
      procedure(eos_p_T_iface  ), deferred :: get_rho_from_p_T
      procedure(eos_p_T_iface  ), deferred :: get_h_from_p_T
      procedure(eos_p_T_iface  ), deferred :: get_s_from_p_T
      procedure(eos_p_T_iface  ), deferred :: get_g_from_p_T
      procedure(eos_rho_e_iface), deferred :: get_gruneisen_from_rho_e
      procedure(eos_p_rho_iface), deferred :: get_rhoe_from_p_rho
      procedure(eos_p_T_iface  ), deferred :: get_rhoe_from_p_T
   end type eos

   abstract interface
      real(WP) function eos_rho_e_iface(this,rho,e)
         import :: WP,eos
         class(eos), intent(in) :: this
         real(WP), intent(in) :: rho,e
      end function eos_rho_e_iface

      real(WP) function eos_p_rho_iface(this,p,rho)
         import :: WP,eos
         class(eos), intent(in) :: this
         real(WP), intent(in) :: p,rho
      end function eos_p_rho_iface

      real(WP) function eos_p_T_iface(this,p,T)
         import :: WP,eos
         class(eos), intent(in) :: this
         real(WP), intent(in) :: p,T
      end function eos_p_T_iface

      real(WP) function eos_rho_T_iface(this,rho,T)
         import :: WP,eos
         class(eos), intent(in) :: this
         real(WP), intent(in) :: rho,T
      end function eos_rho_T_iface
   end interface

end module eos_class
