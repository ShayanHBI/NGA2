!> Calorically perfect ideal-gas EOS class.
module ig_class
   use precision, only: WP
   use eos_class, only: eos
   implicit none
   private

   public :: ig

   type, extends(eos) :: ig
      real(WP) :: gamma=0.0_WP
      real(WP) :: cv   =0.0_WP
      real(WP) :: cp   =0.0_WP   !< cp = gamma*cv, stored at init
      real(WP) :: R    =0.0_WP   !< R  = (gamma-1)*cv, stored at init
      real(WP) :: q    =0.0_WP
      real(WP) :: qp   =0.0_WP
   contains
      procedure, private :: ig_initialize
      generic   :: initialize              =>ig_initialize
      procedure :: get_p_from_rho_e        =>ig_get_p_from_rho_e
      procedure :: get_T_from_p_rho        =>ig_get_T_from_p_rho
      procedure :: get_c_from_p_rho        =>ig_get_c_from_p_rho
      procedure :: get_e_from_p_rho        =>ig_get_e_from_p_rho
      procedure :: get_e_from_p_T          =>ig_get_e_from_p_T
      procedure :: get_p_from_rho_T        =>ig_get_p_from_rho_T
      procedure :: get_rho_from_p_T        =>ig_get_rho_from_p_T
      procedure :: get_h_from_p_T          =>ig_get_h_from_p_T
      procedure :: get_s_from_p_T          =>ig_get_s_from_p_T
      procedure :: get_g_from_p_T          =>ig_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e=>ig_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     =>ig_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       =>ig_get_rhoe_from_p_T
   end type ig

contains

   subroutine ig_initialize(this,gamma,cv,q,qp)
      class(ig), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv,q,qp
      this%gamma=gamma
      this%cv   =cv
      this%cp   =gamma*cv
      this%R    =(gamma-1.0_WP)*cv
      this%q    =q
      this%qp   =qp
   end subroutine ig_initialize

   real(WP) function ig_get_p_from_rho_e(this,rho,e) result(p)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)
   end function ig_get_p_from_rho_e

   real(WP) function ig_get_T_from_p_rho(this,p,rho) result(T)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=p/(this%R*rho)
   end function ig_get_T_from_p_rho

   real(WP) function ig_get_c_from_p_rho(this,p,rho) result(c)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,rho
      c=sqrt(max(0.0_WP,this%gamma*p/rho))
   end function ig_get_c_from_p_rho

   real(WP) function ig_get_e_from_p_rho(this,p,rho) result(e)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,rho
      e=p/((this%gamma-1.0_WP)*rho)+this%q
   end function ig_get_e_from_p_rho

   real(WP) function ig_get_e_from_p_T(this,p,T) result(e)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      e=this%cv*T+this%q
   end function ig_get_e_from_p_T

   real(WP) function ig_get_p_from_rho_T(this,rho,T) result(p)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: rho,T
      p=this%R*rho*T
   end function ig_get_p_from_rho_T

   real(WP) function ig_get_rho_from_p_T(this,p,T) result(rho)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=p/(this%R*T)
   end function ig_get_rho_from_p_T

   real(WP) function ig_get_h_from_p_T(this,p,T) result(h)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      h=this%cp*T+this%q
   end function ig_get_h_from_p_T

   real(WP) function ig_get_s_from_p_T(this,p,T) result(s)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      s=this%cp*log(T)-this%R*log(p)+this%qp
   end function ig_get_s_from_p_T

   real(WP) function ig_get_g_from_p_T(this,p,T) result(g)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      g=(this%cp-this%qp)*T-T*(this%cp*log(T)-this%R*log(p))+this%q
   end function ig_get_g_from_p_T

   real(WP) function ig_get_gruneisen_from_rho_e(this,rho,e) result(gruneisen)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: rho,e
      gruneisen=this%gamma-1.0_WP
   end function ig_get_gruneisen_from_rho_e

   real(WP) function ig_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,rho
      rhoe=p/(this%gamma-1.0_WP)+rho*this%q
   end function ig_get_rhoe_from_p_rho

   real(WP) function ig_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      rhoe=p*(this%cv*T+this%q)/(this%R*T)
   end function ig_get_rhoe_from_p_T

end module ig_class
