!> Noble-Abel stiffened-gas EOS class.
module nasg_class
   use precision, only: WP
   use sg_class,  only: sg
   implicit none
   private

   public :: nasg

   type, extends(sg) :: nasg
      real(WP) :: b=0.0_WP
   contains
      procedure, private :: nasg_initialize
      generic   :: initialize                   =>nasg_initialize
      procedure :: get_p_from_rho_e             =>nasg_get_p_from_rho_e
      procedure :: get_T_from_p_rho             =>nasg_get_T_from_p_rho
      procedure :: get_T_from_p_v               =>nasg_get_T_from_p_v
      procedure :: get_c_from_p_rho             =>nasg_get_c_from_p_rho
      procedure :: get_e_from_p_rho             =>nasg_get_e_from_p_rho
      procedure :: get_p_from_rho_T             =>nasg_get_p_from_rho_T
      procedure :: get_rho_from_p_T             =>nasg_get_rho_from_p_T
      procedure :: get_h_from_p_T               =>nasg_get_h_from_p_T
      procedure :: get_gruneisen_from_rho_e     =>nasg_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho          =>nasg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T            =>nasg_get_rhoe_from_p_T
      procedure :: get_g_from_p_T               =>nasg_get_g_from_p_T
      procedure :: get_drhodT_const_p_from_rho_T=>nasg_get_drhodT_const_p_from_rho_T
      procedure :: get_drhodp_const_T_from_rho_T=>nasg_get_drhodp_const_T_from_rho_T
   end type nasg

contains

   subroutine nasg_initialize(this,pinf,b,gamma,cv,q,qp)
      class(nasg), intent(inout) :: this
      real(WP), intent(in) :: pinf,b,gamma,cv,q,qp
      call this%sg%initialize(pinf=pinf,gamma=gamma,cv=cv,q=q,qp=qp)
      this%b=b
   end subroutine nasg_initialize

   real(WP) function nasg_get_p_from_rho_e(this,rho,e) result(p)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      p=(this%gamma-1.0_WP)*rho*(e-this%q)/(1.0_WP-this%b*rho)-this%gamma*this%pinf
   end function nasg_get_p_from_rho_e

   real(WP) function nasg_get_T_from_p_rho(this,p,rho) result(T)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      T=(p+this%pinf)*(1.0_WP-this%b*rho)/(this%R*rho)
   end function nasg_get_T_from_p_rho

   real(WP) function nasg_get_T_from_p_v(this,p,v) result(T)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,v
      T=(p+this%pinf)*(v-this%b)/this%R
   end function nasg_get_T_from_p_v

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
      p=this%R*rho*T/(1.0_WP-this%b*rho)-this%pinf
   end function nasg_get_p_from_rho_T

   real(WP) function nasg_get_rho_from_p_T(this,p,T) result(rho)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rho=(p+this%pinf)/(this%R*T+this%b*(p+this%pinf))
   end function nasg_get_rho_from_p_T

   real(WP) function nasg_get_h_from_p_T(this,p,T) result(h)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      h=this%cp*T+this%b*p+this%q
   end function nasg_get_h_from_p_T

   real(WP) function nasg_get_gruneisen_from_rho_e(this,rho,e) result(gruneisen)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,e
      gruneisen=(this%gamma-1.0_WP)/(1.0_WP-this%b*rho)
   end function nasg_get_gruneisen_from_rho_e

   real(WP) function nasg_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,rho
      rhoe=(1.0_WP-this%b*rho)*(p+this%gamma*this%pinf)/(this%gamma-1.0_WP)+rho*this%q
   end function nasg_get_rhoe_from_p_rho

   real(WP) function nasg_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      rhoe=((p+this%gamma*this%pinf)*this%cv*T+this%q*(p+this%pinf))/(this%R*T+this%b*(p+this%pinf))
   end function nasg_get_rhoe_from_p_T

   real(WP) function nasg_get_g_from_p_T(this,p,T) result(g)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: p,T
      g=this%cp*T+this%b*p+this%q-T*(this%cp*log(T)-this%R*log(p+this%pinf)+this%qp)
   end function nasg_get_g_from_p_T

   real(WP) function nasg_get_drhodT_const_p_from_rho_T(this,rho,T) result(drhodT)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,T
      drhodT=rho*(this%b*rho-1.0_WP)/T
   end function nasg_get_drhodT_const_p_from_rho_T

   real(WP) function nasg_get_drhodp_const_T_from_rho_T(this,rho,T) result(drhodp)
      class(nasg), intent(in) :: this
      real(WP), intent(in) :: rho,T
      drhodp=(1.0_WP-this%b*rho)**2/(this%R*T)
   end function nasg_get_drhodp_const_T_from_rho_T

end module nasg_class
