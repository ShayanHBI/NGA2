!> Noble-Abel stiffened-gas EOS class.
module nasg_class
   use precision, only: WP
   use sg_class, only: sg
   implicit none
   private

   public :: nasg

   type, extends(sg) :: nasg
      real(WP) :: b = 0.0_WP
   contains
      procedure :: initialize => nasg_initialize

      procedure :: get_p_from_rho_e         => nasg_get_p_from_rho_e
      procedure :: get_T_from_p_rho         => nasg_get_T_from_p_rho
      procedure :: get_c_from_p_rho         => nasg_get_c_from_p_rho
      procedure :: get_e_from_p_rho         => nasg_get_e_from_p_rho
      procedure :: get_p_from_rho_T         => nasg_get_p_from_rho_T
      procedure :: get_rho_from_p_T         => nasg_get_rho_from_p_T
      procedure :: get_h_from_p_T           => nasg_get_h_from_p_T
      procedure :: get_gruneisen_from_rho_e => nasg_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho      => nasg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T        => nasg_get_rhoe_from_p_T
      procedure :: get_g_from_p_T           => nasg_get_g_from_p_T

      procedure :: get_b => nasg_get_b
   end type nasg

contains

   subroutine nasg_initialize(this,gamma,cv,q,qp,pinf,b)
      class(nasg), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv
      real(WP), intent(in), optional :: q,qp,pinf,b
      call this%sg%initialize(gamma,cv,q=q,qp=qp,pinf=pinf)
      this%b = 0.0_WP; if (present(b)) this%b = b
   end subroutine nasg_initialize

   real(WP) function nasg_get_p_from_rho_e(this,rho,e) result(p)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)/(1.0_WP-this%b*rho)-this%gamma*this%pinf
   end function nasg_get_p_from_rho_e

   real(WP) function nasg_get_T_from_p_rho(this,p,rho) result(T)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=(p+this%pinf)*(1.0_WP-this%b*rho)/((this%gamma-1.0_WP)*this%cv*rho)
   end function nasg_get_T_from_p_rho

   real(WP) function nasg_get_c_from_p_rho(this,p,rho) result(c)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      c=sqrt(max(0.0_WP,this%gamma*(p+this%pinf)/(rho*(1.0_WP-this%b*rho))))
   end function nasg_get_c_from_p_rho

   real(WP) function nasg_get_e_from_p_rho(this,p,rho) result(e)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      e=(1.0_WP-this%b*rho)*(p+this%gamma*this%pinf)/((this%gamma-1.0_WP)*rho)+this%q
   end function nasg_get_e_from_p_rho

   real(WP) function nasg_get_p_from_rho_T(this,rho,T) result(p)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,T
      p=(this%gamma-1.0_WP)*this%cv*rho*T/(1.0_WP-this%b*rho)-this%pinf
   end function nasg_get_p_from_rho_T

   real(WP) function nasg_get_rho_from_p_T(this,p,T) result(rho)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=(p+this%pinf)/((this%gamma-1.0_WP)*this%cv*T+this%b*(p+this%pinf))
   end function nasg_get_rho_from_p_T

   real(WP) function nasg_get_h_from_p_T(this,p,T) result(h)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      h=this%gamma*this%cv*T+this%b*p+this%q
   end function nasg_get_h_from_p_T

   real(WP) function nasg_get_gruneisen_from_rho_e(this,rho,e) result(gruneisen)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      gruneisen=(this%gamma-1.0_WP)/(1.0_WP-this%b*rho)
   end function nasg_get_gruneisen_from_rho_e

   !> rho*e = (1 - b*rho)*(p + gamma*pinf)/(gamma-1) + rho*q
   real(WP) function nasg_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      rhoe=(1.0_WP-this%b*rho)*(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function nasg_get_rhoe_from_p_rho

   !> rho*e = rho*(cv*T*(p+gamma*pinf)/(p+pinf) + q),  rho = (p+pinf)/((gamma-1)*cv*T + b*(p+pinf))
   real(WP) function nasg_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP) :: rho
      rho=(p+this%pinf)/((this%gamma-1.0_WP)*this%cv*T+this%b*(p+this%pinf))
      rhoe=rho*(this%cv*T*(p+this%gamma*this%pinf)/(p+this%pinf)+this%q)
   end function nasg_get_rhoe_from_p_T

   !> g = h - T*s = gamma*cv*T + b*p + q - T*(cv*(gamma*ln T - (gamma-1)*ln(p+pinf)) + qp)
   real(WP) function nasg_get_g_from_p_T(this,p,T) result(g)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP) :: pshift,Tsafe
      pshift=max(p+this%pinf,tiny(1.0_WP))
      Tsafe =max(T,tiny(1.0_WP))
      g=this%gamma*this%cv*T+this%b*p+this%q-T*(this%cv*(this%gamma*log(Tsafe)-(this%gamma-1.0_WP)*log(pshift))+this%qp)
   end function nasg_get_g_from_p_T

   real(WP) function nasg_get_b(this) result(b)
      class(nasg), intent(in) :: this
      b=this%b
   end function nasg_get_b

end module nasg_class
