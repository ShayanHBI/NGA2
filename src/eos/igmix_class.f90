!> Arbitrary-species ideal-gas mixture class.
module igmix_class
   use precision, only: WP
   use messager,  only: die
   use ig_class,  only: ig
   use mix_class, only: mix
   implicit none
   private

   public :: igmix

   type :: ig_ptr
      class(ig), pointer :: eos=>null()
   end type ig_ptr

   type, extends(mix) :: igmix
      type(ig_ptr), allocatable :: species(:)
   contains
      procedure :: initialize              =>igmix_initialize
      procedure :: set_species             =>igmix_set_species
      procedure :: get_species_cv          =>igmix_get_species_cv
      procedure :: get_species_cp          =>igmix_get_species_cp
      procedure :: get_species_gamma       =>igmix_get_species_gamma
      procedure :: get_species_q           =>igmix_get_species_q
      procedure :: get_species_qp          =>igmix_get_species_qp
      procedure :: get_p_from_rho_e        =>igmix_get_p_from_rho_e
      procedure :: get_T_from_p_rho        =>igmix_get_T_from_p_rho
      procedure :: get_c_from_p_rho        =>igmix_get_c_from_p_rho
      procedure :: get_e_from_p_rho        =>igmix_get_e_from_p_rho
      procedure :: get_e_from_p_T          =>igmix_get_e_from_p_T
      procedure :: get_p_from_rho_T        =>igmix_get_p_from_rho_T
      procedure :: get_rho_from_p_T        =>igmix_get_rho_from_p_T
      procedure :: get_h_from_p_T          =>igmix_get_h_from_p_T
      procedure :: get_s_from_p_T          =>igmix_get_s_from_p_T
      procedure :: get_gruneisen_from_rho_e=>igmix_get_gruneisen_from_rho_e
      procedure :: get_rhoe_from_p_rho     =>igmix_get_rhoe_from_p_rho
      procedure :: get_rhoe_from_p_T       =>igmix_get_rhoe_from_p_T
      procedure :: get_g_from_p_T          =>igmix_get_g_from_p_T
      procedure :: get_mix_coeffs
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
      use eos_class, only: eos
      class(igmix), intent(inout) :: this
      class(eos), target, intent(in) :: eos_models(:)
      integer :: is
      if (.not.allocated(this%species)) call die('[igmix set_species] Undefined species')
      if (size(eos_models).ne.this%ns) call die('[igmix set_species] Incompatible number of species and eos models')
      do is=1,this%ns
         select type (m=>eos_models(is))
         class is (ig)
            this%species(is)%eos=>m
         class default
            call die('[igmix set_species] Species EOS must be ig')
         end select
      end do
   end subroutine igmix_set_species

   subroutine check_ready(this,y)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      integer :: n
      if (.not.allocated(this%species)) call die('[igmix] mix not initialized')
      if (size(y).ne.this%ns) call die('[igmix] y has wrong size')
      do n=1,this%ns
         if (.not.associated(this%species(n)%eos)) call die('[igmix] unset species pointer')
      end do
   end subroutine check_ready

   subroutine get_normalized_y(this,y,yn)
      class(igmix), intent(in) :: this
      real(WP), dimension(this%ns), intent(in) :: y
      real(WP), dimension(this%ns), intent(out) :: yn
      real(WP) :: ysum
      call check_ready(this,y)
      ysum=sum(y)
      if (ysum.le.tiny(1.0_WP)) call die('[igmix] non-positive species mass-fraction sum')
      yn=y/ysum
   end subroutine get_normalized_y

   real(WP) function species_R(this,n) result(R)
      class(igmix), intent(in) :: this
      integer, intent(in) :: n
      R=this%species(n)%eos%R
   end function species_R

   subroutine get_mix_coeffs(this,y,cv,R,q,qp,cp,gamma)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      real(WP), optional, intent(out) :: cv,R,q,qp,cp,gamma
      real(WP) :: cv_,cp_,q_,qp_
      integer :: n
      cv_=0.0_WP; cp_=0.0_WP; q_=0.0_WP; qp_=0.0_WP
      do n=1,this%ns
         cv_=cv_+y(n)*this%species(n)%eos%cv
         cp_=cp_+y(n)*this%species(n)%eos%cp
         q_ =q_ +y(n)*this%species(n)%eos%q
         qp_=qp_+y(n)*this%species(n)%eos%qp
      end do
      if (present(cv))    cv    =cv_
      if (present(cp))    cp    =cp_
      if (present(q))     q     =q_
      if (present(qp))    qp    =qp_
      if (present(R))     R     =cp_-cv_
      if (present(gamma)) gamma =cp_/cv_
   end subroutine get_mix_coeffs

   subroutine get_mole_fractions(this,y,x)
      class(igmix), intent(in) :: this
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(this%ns), intent(out) :: x
      real(WP) :: denom,R
      integer :: n
      denom=0.0_WP
      do n=1,this%ns
         R=species_R(this,n)
         if (R.le.0.0_WP) call die('[igmix] non-positive species gas constant')
         x(n)=y(n)*R
         denom=denom+x(n)
      end do
      x=x/denom
   end subroutine get_mole_fractions

   real(WP) function igmix_get_p_from_rho_e(this,rho,e,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: q,gamma
      call this%get_mix_coeffs(y=y,q=q,gamma=gamma)
      p=(gamma-1.0_WP)*rho*(e-q)
   end function igmix_get_p_from_rho_e

   real(WP) function igmix_get_T_from_p_rho(this,p,rho,y) result(T)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: R
      call this%get_mix_coeffs(y=y,R=R)
      T=p/(rho*R)
   end function igmix_get_T_from_p_rho

   real(WP) function igmix_get_c_from_p_rho(this,p,rho,y) result(c)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gamma
      call this%get_mix_coeffs(y=y,gamma=gamma)
      c=sqrt(max(0.0_WP,gamma*p/rho))
   end function igmix_get_c_from_p_rho

   real(WP) function igmix_get_e_from_p_rho(this,p,rho,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: q,gamma
      call this%get_mix_coeffs(y=y,q=q,gamma=gamma)
      e=p/((gamma-1.0_WP)*rho)+q
   end function igmix_get_e_from_p_rho

   real(WP) function igmix_get_e_from_p_T(this,p,T,y) result(e)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,q
      call this%get_mix_coeffs(y=y,cv=cv,q=q)
      e=cv*T+q
   end function igmix_get_e_from_p_T

   real(WP) function igmix_get_p_from_rho_T(this,rho,T,y) result(p)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: R
      call this%get_mix_coeffs(y=y,R=R)
      p=rho*R*T
   end function igmix_get_p_from_rho_T

   real(WP) function igmix_get_rho_from_p_T(this,p,T,y) result(rho)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: R
      call this%get_mix_coeffs(y=y,R=R)
      rho=p/(R*T)
   end function igmix_get_rho_from_p_T

   real(WP) function igmix_get_h_from_p_T(this,p,T,y) result(h)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cp,q
      call this%get_mix_coeffs(y=y,cp=cp,q=q)
      h=cp*T+q
   end function igmix_get_h_from_p_T

   real(WP) function igmix_get_s_from_p_T(this,p,T,y) result(s)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP), dimension(this%ns) :: x
      real(WP) :: pp
      integer :: n
      call get_mole_fractions(this,y,x)
      s=0.0_WP
      do n=1,this%ns
         pp=max(x(n)*p,tiny(1.0_WP))
         s=s+y(n)*this%species(n)%eos%get_s_from_p_T(pp,T)
      end do
   end function igmix_get_s_from_p_T

   real(WP) function igmix_get_g_from_p_T(this,p,T,y) result(g)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      g=this%get_h_from_p_T(p,T,y)-T*this%get_s_from_p_T(p,T,y)
   end function igmix_get_g_from_p_T

   real(WP) function igmix_get_gruneisen_from_rho_e(this,rho,e,y) result(gruneisen)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: rho,e
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: gamma
      call this%get_mix_coeffs(y=y,gamma=gamma)
      gruneisen=gamma-1.0_WP
   end function igmix_get_gruneisen_from_rho_e

   real(WP) function igmix_get_rhoe_from_p_rho(this,p,rho,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,rho
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: q,gamma
      call this%get_mix_coeffs(y=y,gamma=gamma,q=q)
      rhoe=p/(gamma-1.0_WP)+rho*q
   end function igmix_get_rhoe_from_p_rho

   real(WP) function igmix_get_rhoe_from_p_T(this,p,T,y) result(rhoe)
      class(igmix), intent(in) :: this
      real(WP), intent(in) :: p,T
      real(WP), dimension(:), intent(in) :: y
      real(WP) :: cv,R,q
      call this%get_mix_coeffs(y=y,cv=cv,R=R,q=q)
      rhoe=p*(cv*T+q)/(R*T)
   end function igmix_get_rhoe_from_p_T

   real(WP) function igmix_get_species_cv(this,is) result(cv)
      class(igmix), intent(in) :: this
      integer, intent(in) :: is
      cv=this%species(is)%eos%cv
   end function igmix_get_species_cv

   real(WP) function igmix_get_species_cp(this,is) result(cp)
      class(igmix), intent(in) :: this
      integer, intent(in) :: is
      cp=this%species(is)%eos%cp
   end function igmix_get_species_cp

   real(WP) function igmix_get_species_gamma(this,is) result(gamma)
      class(igmix), intent(in) :: this
      integer, intent(in) :: is
      gamma=this%species(is)%eos%gamma
   end function igmix_get_species_gamma

   real(WP) function igmix_get_species_q(this,is) result(q)
      class(igmix), intent(in) :: this
      integer, intent(in) :: is
      q=this%species(is)%eos%q
   end function igmix_get_species_q

   real(WP) function igmix_get_species_qp(this,is) result(qp)
      class(igmix), intent(in) :: this
      integer, intent(in) :: is
      qp=this%species(is)%eos%qp
   end function igmix_get_species_qp

end module igmix_class
