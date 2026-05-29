!> Arbitrary-species ideal-gas mixture class.
module igmix_class
   use precision, only: WP
   use ig_class, only: ig
   use mix_class, only: mix
   implicit none
   private

   public :: igmix

   type :: ig_ptr
      class(ig), pointer :: ptr => null()
   end type ig_ptr

   type, extends(mix) :: igmix
      type(ig_ptr), allocatable :: species(:)
   contains
      procedure :: initialize  => igmix_initialize
      procedure :: set_species => igmix_set_species

      procedure :: get_p_from_rho_e         => igmix_get_p_from_rho_e
      procedure :: get_T_from_p_rho         => igmix_get_T_from_p_rho
      procedure :: get_c_from_p_rho         => igmix_get_c_from_p_rho
      procedure :: get_e_from_p_rho         => igmix_get_e_from_p_rho
      procedure :: get_e_from_p_T           => igmix_get_e_from_p_T
      procedure :: get_p_from_rho_T         => igmix_get_p_from_rho_T
      procedure :: get_rho_from_p_T         => igmix_get_rho_from_p_T
      procedure :: get_h_from_p_T           => igmix_get_h_from_p_T
      procedure :: get_s_from_p_T           => igmix_get_s_from_p_T
      procedure :: get_gruneisen_from_rho_e => igmix_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho      => igmix_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T        => igmix_get_rhoe_from_p_T
      procedure :: get_g_from_p_T           => igmix_get_g_from_p_T
   end type igmix

contains

   subroutine igmix_initialize(this,ns)
      class(igmix), intent(inout) :: this
      integer, intent(in) :: ns
      if (allocated(this%species)) deallocate(this%species)
      this%ns=ns
      allocate(this%species(ns))
   end subroutine igmix_initialize

   subroutine igmix_set_species(this,eos_models)
      class(igmix), intent(inout) :: this
      class(ig), target, intent(in) :: eos_models(:)
      integer :: is
      if (.not.allocated(this%species)) error stop '[igmix set_species] mix not initialized'
      if (size(eos_models).ne.this%ns) error stop '[igmix set_species] eos_models size mismatch'
      do is=1,this%ns
         this%species(is)%ptr=>eos_models(is)
      end do
   end subroutine igmix_set_species

   subroutine check_ready(this,y)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      integer :: n
      if (.not.allocated(this%species)) error stop '[igmix] mix not initialized'
      if (size(y).ne.this%ns) error stop '[igmix] y has wrong size'
      do n=1,this%ns
         if (.not.associated(this%species(n)%ptr)) error stop '[igmix] unset species pointer'
      end do
   end subroutine check_ready

   subroutine get_normalized_y(this,y,yn)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(this%ns), intent(out) :: yn
      real(WP) :: ysum
      call check_ready(this,y)
      ysum=sum(y)
      if (ysum.le.tiny(1.0_WP)) error stop '[igmix] non-positive species mass-fraction sum'
      yn=y/ysum
   end subroutine get_normalized_y

   real(WP) function species_R(this,n) result(R)
      class(igmix), intent(in) :: this
      integer, intent(in) :: n
      R=this%species(n)%ptr%get_cp()-this%species(n)%ptr%get_cv()
   end function species_R

   subroutine get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      real(WP), intent(out) :: cv,R,q,cp,gamma
      real(WP), dimension(this%ns) :: yn
      integer :: n
      call get_normalized_y(this,y,yn)
      cv=0.0_WP
      R =0.0_WP
      q =0.0_WP
      do n=1,this%ns
         cv=cv+yn(n)*this%species(n)%ptr%get_cv()
         R =R +yn(n)*species_R(this,n)
         q =q +yn(n)*this%species(n)%ptr%get_q()
      end do
      cp=cv+R
      gamma=cp/cv
   end subroutine get_mix_coeffs

   subroutine get_mole_fractions(this,y,x)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(this%ns), intent(out) :: x
      real(WP), dimension(this%ns) :: yn
      real(WP) :: denom,Ri
      integer :: n
      call get_normalized_y(this,y,yn)
      denom=0.0_WP
      do n=1,this%ns
         Ri=species_R(this,n)
         if (Ri.le.0.0_WP) error stop '[igmix] non-positive species gas constant'
         x(n)=yn(n)*Ri
         denom=denom+x(n)
      end do
      x=x/denom
   end subroutine get_mole_fractions

   real(WP) function igmix_get_p_from_rho_e(this,rho,e,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma,T
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      T=(e-q)/cv
      p=rho*R*T
   end function igmix_get_p_from_rho_e

   real(WP) function igmix_get_T_from_p_rho(this,p,rho,y) result(T)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      T=p/(rho*R)
   end function igmix_get_T_from_p_rho

   real(WP) function igmix_get_c_from_p_rho(this,p,rho,y) result(c)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      c=sqrt(max(0.0_WP,gamma*p/rho))
   end function igmix_get_c_from_p_rho

   real(WP) function igmix_get_e_from_p_rho(this,p,rho,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma,T
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      T=p/(rho*R)
      e=cv*T+q
   end function igmix_get_e_from_p_rho

   real(WP) function igmix_get_e_from_p_T(this,p,T,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      e=cv*T+q
   end function igmix_get_e_from_p_T

   real(WP) function igmix_get_p_from_rho_T(this,rho,T,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      p=rho*R*T
   end function igmix_get_p_from_rho_T

   real(WP) function igmix_get_rho_from_p_T(this,p,T,y) result(rho)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      rho=p/(R*T)
   end function igmix_get_rho_from_p_T

   real(WP) function igmix_get_h_from_p_T(this,p,T,y) result(h)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      h=cp*T+q
   end function igmix_get_h_from_p_T

   real(WP) function igmix_get_s_from_p_T(this,p,T,y) result(s)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(this%ns) :: yn,x
      real(WP) :: pi
      integer :: n
      call get_normalized_y(this,y,yn)
      call get_mole_fractions(this,y,x)
      s=0.0_WP
      do n=1,this%ns
         pi=max(x(n)*p,tiny(1.0_WP))
         s=s+yn(n)*this%species(n)%ptr%get_s_from_p_T(pi,T)
      end do
   end function igmix_get_s_from_p_T

   real(WP) function igmix_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      gruneisen=gamma-1.0_WP
   end function igmix_get_gruneisen_from_rho_e

   !> rho*e = cv*p/R + rho*q
   real(WP) function igmix_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      rhoe=cv*p/R+rho*q
   end function igmix_get_rhoe_from_p_rho

   !> rho*e = rho*(cv*T + q),  rho = p/(R*T)
   real(WP) function igmix_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q,cp,gamma,rho
      call get_mix_coeffs(this,y,cv,R,q,cp,gamma)
      rho=p/(R*T)
      rhoe=rho*(cv*T+q)
   end function igmix_get_rhoe_from_p_T

   !> g = h - T*s  (entropy sum over partial pressures prevents a simpler closed form)
   real(WP) function igmix_get_g_from_p_T(this,p,T,y) result(g)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=this%get_h_from_p_T(p,T,y)-T*this%get_s_from_p_T(p,T,y)
   end function igmix_get_g_from_p_T

end module igmix_class
