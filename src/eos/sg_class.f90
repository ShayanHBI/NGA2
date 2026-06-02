!> Stiffened-gas EOS class.
module sg_class
   use precision, only: WP
   use ig_class,  only: ig
   implicit none
   private

   public :: sg

   type, extends(ig) :: sg
      real(WP) :: pinf=0.0_WP
   contains
      procedure, private :: sg_initialize
      generic   :: initialize         =>sg_initialize
      procedure :: get_p_from_rho_e   =>sg_get_p_from_rho_e
      procedure :: get_T_from_p_rho   =>sg_get_T_from_p_rho
      procedure :: get_c_from_p_rho   =>sg_get_c_from_p_rho
      procedure :: get_e_from_p_rho   =>sg_get_e_from_p_rho
      procedure :: get_e_from_p_T     =>sg_get_e_from_p_T
      procedure :: get_p_from_rho_T   =>sg_get_p_from_rho_T
      procedure :: get_rho_from_p_T   =>sg_get_rho_from_p_T
      procedure :: get_s_from_p_T     =>sg_get_s_from_p_T
      procedure :: get_rhoe_from_p_rho=>sg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T  =>sg_get_rhoe_from_p_T
      procedure :: get_g_from_p_T     =>sg_get_g_from_p_T
   end type sg

contains

   subroutine sg_initialize(this,pinf,gamma,cv,q,qp)
      class(sg), intent(inout) :: this
      real(WP), intent(in) :: pinf,gamma,cv,q,qp
      call this%ig%initialize(gamma=gamma,cv=cv,q=q,qp=qp)
      this%pinf=pinf
   end subroutine sg_initialize

   real(WP) function sg_get_p_from_rho_e(this,rho,e) result(p)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)-this%gamma*this%pinf
   end function sg_get_p_from_rho_e

   real(WP) function sg_get_T_from_p_rho(this,p,rho) result(T)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=(p+this%pinf)/(this%R*rho)
   end function sg_get_T_from_p_rho

   real(WP) function sg_get_c_from_p_rho(this,p,rho) result(c)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      c=sqrt(max(0.0_WP,this%gamma*(p+this%pinf)/rho))
   end function sg_get_c_from_p_rho

   real(WP) function sg_get_e_from_p_rho(this,p,rho) result(e)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      e=(p+this%gamma*this%pinf)/((this%gamma-1.0_WP)*rho)+this%q
   end function sg_get_e_from_p_rho

   real(WP) function sg_get_e_from_p_T(this,p,T) result(e)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      e=this%cv*T*(p+this%gamma*this%pinf)/(p+this%pinf)+this%q
   end function sg_get_e_from_p_T

   real(WP) function sg_get_p_from_rho_T(this,rho,T) result(p)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: rho,T
      p=this%R*rho*T-this%pinf
   end function sg_get_p_from_rho_T

   real(WP) function sg_get_rho_from_p_T(this,p,T) result(rho)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=(p+this%pinf)/(this%R*T)
   end function sg_get_rho_from_p_T

   real(WP) function sg_get_s_from_p_T(this,p,T) result(s)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      s=this%cp*log(T)-this%R*log(p+this%pinf)+this%qp
   end function sg_get_s_from_p_T

   real(WP) function sg_get_g_from_p_T(this,p,T) result(g)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      g=(this%cp-this%qp)*T-T*(this%cp*log(T)-this%R*log(p+this%pinf))+this%q
   end function sg_get_g_from_p_T

   real(WP) function sg_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      rhoe=(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function sg_get_rhoe_from_p_rho

   real(WP) function sg_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rhoe=((p+this%gamma*this%pinf)*this%cv*T+this%q*(p+this%pinf))/(this%R*T)
   end function sg_get_rhoe_from_p_T

end module sg_class
