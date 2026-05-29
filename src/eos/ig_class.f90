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
      real(WP) :: q    =0.0_WP
      real(WP) :: qp   =0.0_WP
   contains
      procedure :: initialize => ig_initialize

      procedure :: get_p_from_rho_e         => ig_get_p_from_rho_e
      procedure :: get_T_from_p_rho         => ig_get_T_from_p_rho
      procedure :: get_c_from_p_rho         => ig_get_c_from_p_rho
      procedure :: get_e_from_p_rho         => ig_get_e_from_p_rho
      procedure :: get_e_from_p_T           => ig_get_e_from_p_T
      procedure :: get_p_from_rho_T         => ig_get_p_from_rho_T
      procedure :: get_rho_from_p_T         => ig_get_rho_from_p_T
      procedure :: get_h_from_p_T           => ig_get_h_from_p_T
      procedure :: get_s_from_p_T           => ig_get_s_from_p_T
      procedure :: get_gruneisen_from_rho_e => ig_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho      => ig_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T        => ig_get_rhoe_from_p_T
      procedure :: get_g_from_p_T           => ig_get_g_from_p_T

      ! Caloric parameter accessors — inherited by all subclasses.
      procedure :: get_gamma => ig_get_gamma
      procedure :: get_q     => ig_get_q
      procedure :: get_qp    => ig_get_qp
      procedure :: get_cv    => ig_get_cv
      procedure :: get_cp    => ig_get_cp
   end type ig

contains

   subroutine ig_initialize(this,gamma,cv,q,qp,pinf,b)
      class(ig), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv
      real(WP), intent(in), optional :: q,qp,pinf,b
      this%gamma=gamma
      this%cv   =cv
      if (present(q))  this%q =q
      if (present(qp)) this%qp=qp
   end subroutine ig_initialize

   real(WP) function ig_get_p_from_rho_e(this,rho,e) result(p)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)
   end function ig_get_p_from_rho_e

   real(WP) function ig_get_T_from_p_rho(this,p,rho) result(T)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=p/((this%gamma-1.0_WP)*this%cv*rho)
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
      p=(this%gamma-1.0_WP)*this%cv*rho*T
   end function ig_get_p_from_rho_T

   real(WP) function ig_get_rho_from_p_T(this,p,T) result(rho)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=p/((this%gamma-1.0_WP)*this%cv*T)
   end function ig_get_rho_from_p_T

   real(WP) function ig_get_h_from_p_T(this,p,T) result(h)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      h=this%gamma*this%cv*T+this%q
   end function ig_get_h_from_p_T

   real(WP) function ig_get_s_from_p_T(this,p,T) result(s)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      s=this%cv*(this%gamma*log(T)-(this%gamma-1.0_WP)*log(p))+this%qp
   end function ig_get_s_from_p_T

   real(WP) function ig_get_g_from_p_T(this,p,T) result(g)
      class(ig), intent(in) :: this
      real(WP), intent(in) :: p,T
      g=(this%gamma*this%cv-this%qp)*T-this%cv*T*(this%gamma*log(T)-(this%gamma-1.0_WP)*log(p))
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
      real(WP) :: rho
      rho=p/((this%gamma-1.0_WP)*this%cv*T)
      rhoe=rho*(this%cv*T+this%q)
   end function ig_get_rhoe_from_p_T

   real(WP) function ig_get_gamma(this) result(gamma)
      class(ig), intent(in) :: this
      gamma=this%gamma
   end function ig_get_gamma

   real(WP) function ig_get_cv(this) result(cv)
      class(ig), intent(in) :: this
      cv=this%cv
   end function ig_get_cv

   real(WP) function ig_get_q(this) result(q)
      class(ig), intent(in) :: this
      q=this%q
   end function ig_get_q

   real(WP) function ig_get_qp(this) result(qp)
      class(ig), intent(in) :: this
      qp=this%qp
   end function ig_get_qp

   real(WP) function ig_get_cp(this) result(cp)
      class(ig), intent(in) :: this
      cp=this%gamma*this%cv
   end function ig_get_cp

end module ig_class
