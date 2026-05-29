!> Stiffened-gas EOS class.
module sg_class
   use precision, only: WP
   use ig_class, only: ig
   implicit none
   private

   public :: sg

   type, extends(ig) :: sg
      real(WP) :: pinf = 0.0_WP
   contains
      procedure :: initialize => sg_initialize

      procedure :: get_p_from_rho_e    => sg_get_p_from_rho_e
      procedure :: get_T_from_p_rho    => sg_get_T_from_p_rho
      procedure :: get_c_from_p_rho    => sg_get_c_from_p_rho
      procedure :: get_e_from_p_rho    => sg_get_e_from_p_rho
      procedure :: get_e_from_p_T      => sg_get_e_from_p_T
      procedure :: get_p_from_rho_T    => sg_get_p_from_rho_T
      procedure :: get_rho_from_p_T    => sg_get_rho_from_p_T
      procedure :: get_s_from_p_T      => sg_get_s_from_p_T
      procedure :: get_rhoe_from_p_rho => sg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T   => sg_get_rhoe_from_p_T
      procedure :: get_g_from_p_T      => sg_get_g_from_p_T

      procedure :: get_pinf => sg_get_pinf
   end type sg

contains

   subroutine sg_initialize(this,gamma,cv,q,qp,pinf,b)
      class(sg), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv
      real(WP), intent(in), optional :: q,qp,pinf,b
      call this%ig%initialize(gamma,cv,q=q,qp=qp)
      this%pinf = 0.0_WP; if (present(pinf)) this%pinf = pinf
   end subroutine sg_initialize

   real(WP) function sg_get_p_from_rho_e(this,rho,e) result(p)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)-this%gamma*this%pinf
   end function sg_get_p_from_rho_e

   real(WP) function sg_get_T_from_p_rho(this,p,rho) result(T)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=(p+this%pinf)/((this%gamma-1.0_WP)*this%cv*rho)
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
      p=(this%gamma-1.0_WP)*this%cv*rho*T-this%pinf
   end function sg_get_p_from_rho_T

   real(WP) function sg_get_rho_from_p_T(this,p,T) result(rho)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=(p+this%pinf)/((this%gamma-1.0_WP)*this%cv*T)
   end function sg_get_rho_from_p_T

   real(WP) function sg_get_s_from_p_T(this,p,T) result(s)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP) :: pshift,Tsafe
      pshift=max(p+this%pinf,tiny(1.0_WP))
      Tsafe =max(T,tiny(1.0_WP))
      s=this%cv*(this%gamma*log(Tsafe)-(this%gamma-1.0_WP)*log(pshift))+this%qp
   end function sg_get_s_from_p_T

   real(WP) function sg_get_g_from_p_T(this,p,T) result(g)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP) :: pshift,Tsafe
      pshift=max(p+this%pinf,tiny(1.0_WP))
      Tsafe =max(T,tiny(1.0_WP))
      g=this%gamma*this%cv*T+this%q-T*(this%cv*(this%gamma*log(Tsafe)-(this%gamma-1.0_WP)*log(pshift))+this%qp)
   end function sg_get_g_from_p_T

   !> rho*e = (p + gamma*pinf)/(gamma-1) + rho*q
   real(WP) function sg_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      rhoe=(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function sg_get_rhoe_from_p_rho

   !> rho*e = rho*(cv*T*(p+gamma*pinf)/(p+pinf) + q),  rho = (p+pinf)/((gamma-1)*cv*T)
   real(WP) function sg_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(sg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP) :: rho
      rho=(p+this%pinf)/((this%gamma-1.0_WP)*this%cv*T)
      rhoe=rho*(this%cv*T*(p+this%gamma*this%pinf)/(p+this%pinf)+this%q)
   end function sg_get_rhoe_from_p_T

   real(WP) function sg_get_pinf(this) result(pinf)
      class(sg), intent(in) :: this
      pinf=this%pinf
   end function sg_get_pinf

end module sg_class
