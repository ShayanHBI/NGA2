!> Extended Noble-Abel Stiffened-Gas (ENASG) EOS class.
!> Chiapolino et al. (2018), Fluids 3(3):48.
!> Extends NASG with temperature-dependent attractive pressure:
!>   p_inf(T)=p_inf1*T+p_inf0
!> and volume-dependent covolume:
!>   b(v)=b1*v+b0
!> Reduces to NASG when p_inf1=0, b1=0.
module enasg_class
   use precision,  only: WP
   use nasg_class, only: nasg
   implicit none
   private

   public :: enasg

   type, extends(nasg) :: enasg
      !> Temperature coefficient of the attractive pressure [Pa/K]
      real(WP) :: pinf1=0.0_WP
      !> Linear coefficient of the covolume [-]
      real(WP) :: b1=0.0_WP
      !> Derived constants (computed once at init, never change)
      !> pp_inf0=gamma*pinf0*(1-b1)/(gamma-b1)   [Pa]
      real(WP) :: pp_inf0=0.0_WP
   contains
      procedure, private :: enasg_initialize
      generic   :: initialize                   =>enasg_initialize
      procedure :: get_p_from_rho_e             =>enasg_get_p_from_rho_e
      procedure :: get_T_from_p_rho             =>enasg_get_T_from_p_rho
      procedure :: get_T_from_p_v               =>enasg_get_T_from_p_v
      procedure :: get_c_from_p_rho             =>enasg_get_c_from_p_rho
      procedure :: get_e_from_p_rho             =>enasg_get_e_from_p_rho
      procedure :: get_e_from_p_T               =>enasg_get_e_from_p_T
      procedure :: get_p_from_rho_T             =>enasg_get_p_from_rho_T
      procedure :: get_rho_from_p_T             =>enasg_get_rho_from_p_T
      procedure :: get_h_from_p_T               =>enasg_get_h_from_p_T
      procedure :: get_s_from_p_T               =>enasg_get_s_from_p_T
      procedure :: get_g_from_p_T               =>enasg_get_g_from_p_T
      procedure :: get_gruneisen_from_rho_e     =>enasg_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho          =>enasg_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T            =>enasg_get_rhoe_from_p_T
      procedure :: get_pinf_T                   =>enasg_get_pinf_T
      procedure :: get_ppinf_T                  =>enasg_get_ppinf_T
      procedure :: get_drhodT_const_p_from_rho_T=>enasg_get_drhodT_const_p_from_rho_T
      procedure :: get_drhodp_const_T_from_rho_T=>enasg_get_drhodp_const_T_from_rho_T
   end type enasg

contains

   !> Initialize with all ENASG parameters.
   !> pinf0 maps to the parent nasg%pinf (constant part of p_inf).
   !> b0    maps to the parent nasg%b    (constant part of covolume).
   subroutine enasg_initialize(this,gamma,cv,pinf0,pinf1,b0,b1,q,qp)
      class(enasg), intent(inout) :: this
      real(WP), intent(in) :: gamma,cv
      real(WP), intent(in) :: pinf0   !< p_inf0 [Pa]
      real(WP), intent(in) :: pinf1   !< p_inf1 [Pa/K]
      real(WP), intent(in) :: b0      !< b0 [m^3/kg]
      real(WP), intent(in) :: b1      !< b1 [-]
      real(WP), intent(in) :: q,qp
      ! Parent stores pinf0->pinf, b0->b, and computes cp, R
      call this%nasg%initialize(pinf=pinf0,b=b0,gamma=gamma,cv=cv,q=q,qp=qp)
      this%pinf1  =pinf1
      this%b1     =b1
      ! Derived constant: tilde_p'_inf=gamma*pinf0*(1-b1)/(gamma-b1)
      if (abs(gamma-b1).gt.tiny(1.0_WP)) then
         this%pp_inf0=gamma*pinf0*(1.0_WP-b1)/(gamma-b1)
      else
         this%pp_inf0=0.0_WP
      end if
   end subroutine enasg_initialize

   ! ---------------------------------------------------------------------------
   ! ENASG helper functions
   ! ---------------------------------------------------------------------------

   !> p_inf(T)=pinf1*T+pinf0
   real(WP) function enasg_get_pinf_T(this,T) result(pinf_T)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: T
      pinf_T=this%pinf1*T+this%pinf
   end function enasg_get_pinf_T

   !> p'_inf(T)=gamma*pinf1*T+tilde_p'_inf
   real(WP) function enasg_get_ppinf_T(this,T) result(ppinf_T)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: T
      ppinf_T=this%gamma*this%pinf1*T+this%pp_inf0
   end function enasg_get_ppinf_T

   !> T(p, v): invert the ENASG thermal EOS at fixed (p, v). Inverting
   !>   v=(gamma-1)*cv*T/((1-b1)*(p+p'_inf(T)))+b0/(1-b1)
   !> for T at fixed v gives, with w=(1-b1)*v-b0:
   !>   T=(p+tilde_p'_inf)*w / ((gamma-1)*cv-gamma*pinf1*w)
   real(WP) function enasg_get_T_from_p_v(this,p,v) result(T)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,v
      real(WP) :: w,denom
      w     =(1.0_WP-this%b1)*v-this%b
      denom =(this%gamma-1.0_WP)*this%cv-this%gamma*this%pinf1*w
      if (abs(denom).gt.tiny(1.0_WP)) then
         T=(p+this%pp_inf0)*w/denom
      else
         T=0.0_WP
      end if
   end function enasg_get_T_from_p_v

   ! ---------------------------------------------------------------------------
   ! EOS interface implementations
   ! ---------------------------------------------------------------------------

   !> p(rho, e): thermal EOS inverted from the caloric form.
   !>   p=(gamma-1)*cv*T/(v-b(v))-p'_inf(T)
   !> Because p'_inf depends on T(e,v), we solve for p from
   !>   T(e,v) using the caloric identity then evaluate the thermal EOS.
   !> The caloric identity, obtained by eliminating p between
   !>   e=(p+gamma*p_inf(T))*cv*T/(p+p'_inf(T))+q  and
   !>   v-b(v)=(gamma-1)*cv*T/(p+p'_inf(T)),
   !> is e=cv*T+gamma*pinf0*(v-b(v))/(gamma-b1)+q, with v=1/rho and
   !> b(v)=b1*v+b0, giving p via the thermal EOS above.
   real(WP) function enasg_get_p_from_rho_e(this,rho,e) result(p)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: rho,e
      real(WP) :: v,bv,vmbv,T,ppinf_T
      v     =1.0_WP/rho
      bv    =this%b1*v+this%b
      vmbv  =v-bv
      T=(e-this%q-this%gamma*this%pinf*vmbv/(this%gamma-this%b1))/this%cv
      ppinf_T=this%gamma*this%pinf1*T+this%pp_inf0
      p=(this%gamma-1.0_WP)*this%cv*T/vmbv-ppinf_T
   end function enasg_get_p_from_rho_e

   !> T(p, rho): invert ENASG thermal EOS at fixed (p, rho=1/v).
   real(WP) function enasg_get_T_from_p_rho(this,p,rho) result(T)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,rho
      T=this%get_T_from_p_v(p, 1.0_WP/rho)
   end function enasg_get_T_from_p_rho

   !> c(p, rho): frozen sound speed (eq. 42 in Chiapolino 2018).
   !>   c^2=v^2*(gamma-b1)*(gamma-1)*cv*(p+pp_inf0)
   !>        /[(v-bv)*cv*((gamma-1)*cv-gamma*pinf1*(v-bv))]
   real(WP) function enasg_get_c_from_p_rho(this,p,rho) result(c)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,rho
      real(WP) :: v,bv,vmbv,num,den,c2
      v   =1.0_WP/rho
      bv  =this%b1*v+this%b
      vmbv=v-bv
      num =v**2*(this%gamma-this%b1)*(this%gamma-1.0_WP)*this%cv*(p+this%pp_inf0)
      den =vmbv*this%cv*((this%gamma-1.0_WP)*this%cv-this%gamma*this%pinf1*vmbv)
      if (abs(den).gt.tiny(1.0_WP)) then
         c2=num/den
      else
         c2=0.0_WP
      end if
      c=sqrt(max(0.0_WP,c2))
   end function enasg_get_c_from_p_rho

   !> e(p, rho): specific internal energy at fixed (p, rho).
   real(WP) function enasg_get_e_from_p_rho(this,p,rho) result(e)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,rho
      real(WP) :: T
      T=this%get_T_from_p_rho(p,rho)
      e=this%get_e_from_p_T(p,T)
   end function enasg_get_e_from_p_rho

   !> e(p, T): ENASG caloric EOS (eq. 22 in Chiapolino 2018).
   !>   e=(p+gamma*p_inf(T))/(p+p'_inf(T))*cv*T+q
   real(WP) function enasg_get_e_from_p_T(this,p,T) result(e)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: pinf_T,ppinf_T,denom
      pinf_T =this%get_pinf_T(T)
      ppinf_T=this%get_ppinf_T(T)
      denom  =p+ppinf_T
      if (abs(denom).gt.tiny(1.0_WP)) then
         e=(p+this%gamma*pinf_T)*this%cv*T/denom+this%q
      else
         e=this%q
      end if
   end function enasg_get_e_from_p_T

   !> p(rho, T): ENASG thermal EOS (eq. 18 in Chiapolino 2018).
   !>   p=(gamma-1)*cv*T/(v-b(v))-p'_inf(T)
   real(WP) function enasg_get_p_from_rho_T(this,rho,T) result(p)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: rho,T
      real(WP) :: v,bv,ppinf_T
      v      =1.0_WP/rho
      bv     =this%b1*v+this%b
      ppinf_T=this%get_ppinf_T(T)
      p=(this%gamma-1.0_WP)*this%cv*T/(v-bv)-ppinf_T
   end function enasg_get_p_from_rho_T

   !> rho(p, T): ENASG specific volume inverted (eq. 18 in Chiapolino 2018).
   !>   v=(gamma-1)*cv*T/((1-b1)*(p+p'_inf(T)))+b0/(1-b1)
   real(WP) function enasg_get_rho_from_p_T(this,p,T) result(rho)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: ppinf_T,v
      ppinf_T=this%get_ppinf_T(T)
      v=(this%gamma-1.0_WP)*this%cv*T/((1.0_WP-this%b1)*(p+ppinf_T)) &
       +this%b/(1.0_WP-this%b1)
      rho=1.0_WP/v
   end function enasg_get_rho_from_p_T

   !> h(p, T): specific enthalpy (eq. 47 in Chiapolino 2018).
   real(WP) function enasg_get_h_from_p_T(this,p,T) result(h)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: pinf_T,ppinf_T,one_mb1
      pinf_T =this%get_pinf_T(T)
      ppinf_T=this%get_ppinf_T(T)
      one_mb1=1.0_WP-this%b1
      if (abs(ppinf_T+p).gt.tiny(1.0_WP)) then
         h=this%cv*T*(this%gamma*(p+pinf_T)-p*this%b1-this%gamma*this%b1*pinf_T) &
          /(one_mb1*(p+ppinf_T)) &
          +p*this%b/(one_mb1)+this%q
      else
         h=this%q
      end if
   end function enasg_get_h_from_p_T

   !> s(p, T): specific entropy (eq. 38 in Chiapolino 2018).
   real(WP) function enasg_get_s_from_p_T(this,p,T) result(s)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: ppinf_T,one_mb1,exp_T,exp_p
      ppinf_T=this%get_ppinf_T(T)
      one_mb1=1.0_WP-this%b1
      exp_T  =(this%gamma-this%b1)/one_mb1
      exp_p  =(this%gamma-1.0_WP) /one_mb1
      s=this%cv*( exp_T*log(max(T,tiny(1.0_WP)))-exp_p*log(max(p+ppinf_T,tiny(1.0_WP))) ) &
       -this%gamma*this%pinf1*(this%gamma-1.0_WP)*this%cv*T/(one_mb1*(p+ppinf_T))              &
       +this%qp
   end function enasg_get_s_from_p_T

   !> g(p, T): Gibbs free energy (eq. 50 in Chiapolino 2018).
   real(WP) function enasg_get_g_from_p_T(this,p,T) result(g)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: pinf_T,ppinf_T,one_mb1,exp_T,exp_p,bracket
      pinf_T =this%get_pinf_T(T)
      ppinf_T=this%get_ppinf_T(T)
      one_mb1=1.0_WP-this%b1
      exp_T  =(this%gamma-this%b1)/one_mb1
      exp_p  =(this%gamma-1.0_WP) /one_mb1
      ! Bracket: cv/(1-b1)*(gamma*(p+pinf)-p*b1-gamma*b1*pinf)/(p+ppinf)-qp
      if (abs(ppinf_T+p).gt.tiny(1.0_WP)) then
         bracket=this%cv/one_mb1*(this%gamma*(p+pinf_T)-p*this%b1-this%gamma*this%b1*pinf_T) &
                  /(p+ppinf_T)-this%qp
      else
         bracket=-this%qp
      end if
      g=bracket*T &
       -this%cv*T*( exp_T*log(max(T,tiny(1.0_WP)))-exp_p*log(max(p+ppinf_T,tiny(1.0_WP))) ) &
       +p*this%b/one_mb1+this%q &
       +this%gamma*this%pinf1*(this%gamma-1.0_WP)*this%cv*T**2 &
         /(one_mb1*(p+ppinf_T))
   end function enasg_get_g_from_p_T

   !> Gruneisen coefficient: Gamma=(gamma-1)/(1-b1-b0*rho)
   !> Note: for ENASG with b(v)=b1*v+b0, the effective covolume fraction is
   !>   rho*(b1*v+b0)=b1+b0*rho, so 1-rho*b(v)=1-b1-b0*rho.
   real(WP) function enasg_get_gruneisen_from_rho_e(this,rho,e) result(gruneisen)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: rho,e
      real(WP) :: one_mb
      one_mb=1.0_WP-this%b1-this%b*rho
      if (abs(one_mb).gt.tiny(1.0_WP)) then
         gruneisen=(this%gamma-1.0_WP)/one_mb
      else
         gruneisen=0.0_WP
      end if
   end function enasg_get_gruneisen_from_rho_e

   !> rho*e(p, rho): volumetric internal energy.
   real(WP) function enasg_get_rhoe_from_p_rho(this,p,rho) result(rhoe)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,rho
      rhoe=rho*this%get_e_from_p_rho(p,rho)
   end function enasg_get_rhoe_from_p_rho

   !> rho*e(p, T): volumetric internal energy at fixed (p,T).
   real(WP) function enasg_get_rhoe_from_p_T(this,p,T) result(rhoe)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: p,T
      real(WP) :: rho
      rho =this%get_rho_from_p_T(p,T)
      rhoe=rho*this%get_e_from_p_T(p,T)
   end function enasg_get_rhoe_from_p_T

   real(WP) function enasg_get_drhodT_const_p_from_rho_T(this,rho,T) result(drhodT)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: rho,T
      real(WP) :: one_mb1,one_mb,A,drhodp
      one_mb1=1.0_WP-this%b1
      one_mb =one_mb1-this%b*rho
      A      =this%R*T*rho/one_mb
      drhodp =one_mb**2/(one_mb1*this%R*T)
      drhodT =-drhodp*(A-this%gamma*this%pinf1*T)/T
   end function enasg_get_drhodT_const_p_from_rho_T

   real(WP) function enasg_get_drhodp_const_T_from_rho_T(this,rho,T) result(drhodp)
      class(enasg), intent(in) :: this
      real(WP),     intent(in) :: rho,T
      real(WP) :: one_mb1,one_mb
      one_mb1=1.0_WP-this%b1
      one_mb =one_mb1-this%b*rho
      drhodp =one_mb**2/(one_mb1*this%R*T)
   end function enasg_get_drhodp_const_T_from_rho_T

end module enasg_class
