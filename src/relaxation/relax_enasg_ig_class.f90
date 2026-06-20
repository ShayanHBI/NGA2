!> ENASG liquid and ideal gas relaxation model
module relax_enasg_ig_class
   use precision,           only: WP
   use relax_nasg_ig_class, only: relax_nasg_ig
   use sg_class,            only: sg
   use enasg_class,         only: enasg
   use igmix_class,         only: igmix
   implicit none
   private

   public :: relax_enasg_ig

   type, extends(relax_nasg_ig) :: relax_enasg_ig
      !> Typed pointer for direct ENASG parameter access
      type(enasg), pointer :: liq_enasg=>null()
   contains
      procedure :: initialize   =>relax_enasg_ig_initialize
      procedure :: relax_p      =>relax_enasg_ig_relax_p
      procedure :: relax_pT     =>relax_enasg_ig_relax_pT
      procedure :: get_T_lvg    =>relax_enasg_ig_get_T_lvg
      procedure :: get_coeffs_lv=>relax_enasg_ig_get_coeffs_lv
      procedure :: pTsat        =>relax_enasg_ig_pTsat
      procedure :: dpTsatdT     =>relax_enasg_ig_dpTsatdT
      procedure :: dpTsatdp_lv  =>relax_enasg_ig_dpTsatdp_lv
      procedure :: dpTsatdlnp   =>relax_enasg_ig_dpTsatdlnp
      procedure :: get_pvsat    =>relax_enasg_ig_get_pvsat
   end type relax_enasg_ig

contains

   ! ===========================================================================
   ! Initialization
   ! ===========================================================================

   !> Initialize: call grandparent (relax_nasg_ig), which sets liq, gas,
   !> indV/A, and the AS-ES saturation-curve coefficients (used only as an
   !> initial guess by activate_chem; the exact ENASG saturation curve is
   !> implemented below via pTsat/dpTsatdT/dpTsatdp_lv/dpTsatdlnp/get_pvsat).
   !> Then set the typed ENASG pointer for direct parameter access (pinf1,
   !> b1, pp_inf0).
   subroutine relax_enasg_ig_initialize(this,liq,gas,indV,indA)
      class(relax_enasg_ig), intent(inout) :: this
      class(sg),    target,  intent(in)    :: liq
      class(igmix), target,  intent(in)    :: gas
      integer,               intent(in)    :: indV,indA
      call this%relax_nasg_ig%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      ! Set typed ENASG pointer
      select type (liq)
      type is (enasg)
         this%liq_enasg=>liq
      end select
   end subroutine relax_enasg_ig_initialize

   ! ===========================================================================
   ! Mechanical (p) relaxation
   ! ===========================================================================

   !> Mechanical relaxation for ENASG liquid+ideal gas mixture.
   !>
   !> From the gas energy-jump equation with e_g(p,v_g) linear in p,
   !> alpha*(p*) is the explicit rational function:
   !>
   !>   alpha*(p*)=[p*/(gammaG-1)+pIbar*alpha0-Ag]/[p*/(gammaG-1)+pIbar]
   !>
   !> where Ag=pG0*(1-alpha0)/(gammaG-1) and pIbar=(pL0+p*)/2 (pIbar depends
   !> on p*, so the solve is implicit but scalar).
   !>
   !> The liquid energy-jump equation with the ENASG e_l(p*,v_l*(p*)) gives
   !> the residual:
   !>
   !>   f(p*)=(alpha*rho_l)^0*e_l*(p*)-(alpha*rho_l*e_l)^0+pIbar*(alpha*(p*)-alpha0)=0
   !>
   !> Solved by damped Newton in ln(p*), with the Jacobian from finite
   !> differences.
   subroutine relax_enasg_ig_relax_p(this,VF,Q,Pjump)
      class(relax_enasg_ig),   intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: y
      real(WP) :: arhoL0,arhoG0,arhoeL0,arhoeG0
      real(WP) :: alpha0,pL0,pG0,cvG,cpG,qG,gammaG,Yv
      real(WP) :: Ag,pIbar,alphastar,vLstar,TLstar,eLstar
      real(WP) :: fp,fpp,dlnp,lnp,p,pold,alpha_try,vL_try,TL_try,eL_try,f_try
      integer  :: it
      real(WP), parameter :: fd_eps=1.0e-7_WP,tol=1.0e-8_WP,p_eps=1.0e-10_WP
      integer,  parameter :: itmax=40
      ! Gas composition
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Pre-relaxation invariants
      arhoL0 =Q(1)
      arhoG0 =Q(2)
      arhoeL0=Q(3)
      arhoeG0=Q(4)
      alpha0 =VF
      pL0=this%liq%get_p_from_rho_e(rho=arhoL0/alpha0,e=arhoeL0/arhoL0)
      pG0=this%gas%get_p_from_rho_e(rho=arhoG0/(1.0_WP-alpha0),e=arhoeG0/arhoG0,y=y)
      ! Handle limit cases
      if (pL0.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED! (ENASG)",pL0,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(8)=Yv*Q(2); deallocate(y); return
      end if
      if (pG0.le.0.0_WP) then
         print*,"*** GAS CLIPPED! (ENASG)",pG0,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP; deallocate(y); return
      end if
      ! Ag=pG0*(1-alpha0)/(gammaG-1), constant through the Newton iteration
      Ag=pG0*(1.0_WP-alpha0)/(gammaG-1.0_WP)
      ! Initial guess: arithmetic mean pressure
      p=0.5_WP*(pL0+pG0)
      p=max(p,p_eps)
      ! Newton iterations in ln(p)
      do it=1,itmax
         ! Interface pressure (depends on p*)
         pIbar=0.5_WP*(pL0+p)
         ! alpha*(p) from gas energy-jump eq.
         alphastar=(p/(gammaG-1.0_WP)+pIbar*alpha0-Ag) &
                  /(p/(gammaG-1.0_WP)+pIbar)
         alphastar=max(0.0_WP,min(1.0_WP,alphastar))
         ! Liquid specific volume and energy at (p*,alpha*)
         vLstar=alphastar/arhoL0
         TLstar=this%liq%get_T_from_p_rho(p=p,rho=1.0_WP/vLstar)
         eLstar=this%liq%get_e_from_p_T(p=p,T=TLstar)
         ! Residual f(p*)=arhoL0*eL*-arhoeL0+pIbar*(alpha*-alpha0)
         fp=arhoL0*eLstar-arhoeL0+pIbar*(alphastar-alpha0)
         if (abs(fp).lt.tol*max(abs(arhoeL0),1.0_WP)) exit
         ! Derivative df/dp via finite differences
         pold=p
         p=pold*(1.0_WP+fd_eps)
         pIbar=0.5_WP*(pL0+p)
         alpha_try=(p/(gammaG-1.0_WP)+pIbar*alpha0-Ag) &
                  /(p/(gammaG-1.0_WP)+pIbar)
         alpha_try=max(0.0_WP,min(1.0_WP,alpha_try))
         vL_try=alpha_try/arhoL0
         TL_try=this%liq%get_T_from_p_rho(p=p,rho=1.0_WP/vL_try)
         eL_try=this%liq%get_e_from_p_T(p=p,T=TL_try)
         f_try=arhoL0*eL_try-arhoeL0+pIbar*(alpha_try-alpha0)
         ! df/d(lnp)=p*df/dp
         fpp=pold*(f_try-fp)/(pold*fd_eps)
         p=pold
         if (abs(fpp).lt.tiny(1.0_WP)) exit
         ! Newton step in ln(p), capped at 80% change
         dlnp=-fp/fpp
         if (abs(dlnp).gt.0.8_WP) dlnp=sign(0.8_WP,dlnp)
         lnp=log(p)+dlnp
         p=exp(lnp)
         p=max(p,p_eps)
         if (abs(dlnp).lt.tol) exit
      end do
      ! Final alpha* and state
      pIbar=0.5_WP*(pL0+p)
      alphastar=(p/(gammaG-1.0_WP)+pIbar*alpha0-Ag) &
               /(p/(gammaG-1.0_WP)+pIbar)
      alphastar=max(0.0_WP,min(1.0_WP,alphastar))
      if (alphastar.le.0.0_WP.or.alphastar.ge.1.0_WP) then
         deallocate(y); return
      end if
      ! Update conserved quantities (masses unchanged)
      VF=alphastar
      Q(3)=VF      *this%liq%get_rhoe_from_p_rho(p=p,rho=arhoL0/VF)
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=p,rho=arhoG0/(1.0_WP-VF),y=y)
      deallocate(y)
   end subroutine relax_enasg_ig_relax_p

   ! ===========================================================================
   ! Mechanical-thermal (pT) relaxation
   ! ===========================================================================

   !> pT relaxation for ENASG liquid+ideal gas mixture.
   !>
   !> Step 1: mechanical relaxation (relax_p).
   !>
   !> Step 2: the saturation constraint (sum of volume fractions=1) with the
   !> ENASG v_l(p,T) and ideal-gas v_g(p,T) gives a quadratic in T**(p**):
   !>
   !>   aT*T**^2+bT*T**+cT=0
   !>   aT=(arhoG)^0*(gammaG-1)*cvG*gamma_l*pinf1
   !>   bT=(arhoG)^0*(gammaG-1)*cvG*(p**+pp_inf0)+(arhoL)^0*hatCl*p**
   !>      -gamma_l*pinf1*p**(1-(arhoL)^0*hatbl)
   !>   cT=-p**(p**+pp_inf0)(1-(arhoL)^0*hatbl)
   !>
   !> with hatCl=(gamma_l-1)*cv_l/(1-b1) and hatbl=b0/(1-b1). The positive
   !> root gives T**(p**) (linear fallback in the NASG limit pinf1=0).
   !>
   !> Energy conservation
   !>   f(p**)=(arhoL)^0*e_l(p**,T**)+(arhoG)^0*e_g(T**)-rhoe0=0
   !> is then solved by damped Newton in ln(p**), same FD-Jacobian pattern
   !> as relax_p.
   subroutine relax_enasg_ig_relax_pT(this,VF,Q,Pjump)
      class(relax_enasg_ig),   intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: y
      real(WP) :: arhoL0,arhoG0,rhoe0
      real(WP) :: cvG,cpG,qG,gammaG,Yv
      real(WP) :: hatCl,hatbl,pp_inf0,gam_pinf1,one_mb1
      real(WP) :: Tstar,VFeq
      real(WP) :: eLstar,eGstar
      real(WP) :: fp,fpp,p,pold,dlnp,lnp
      real(WP) :: Ts_p,eL_p,eG_p,f_p
      integer  :: it
      real(WP), parameter :: fd_eps=1.0e-7_WP,tol=1.0e-8_WP,p_eps=1.0e-10_WP
      integer,  parameter :: itmax=40
      ! ---- Step 1: mechanical relaxation ----
      call this%relax_p(VF,Q,Pjump)
      ! ---- Step 2: thermal relaxation ----
      ! Gas composition (unchanged by relax_p)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Invariants (partial densities and total energy from pre-relaxation)
      arhoL0=Q(1)
      arhoG0=Q(2)
      rhoe0 =Q(3)+Q(4)
      ! ENASG shorthands
      one_mb1  =1.0_WP-this%liq_enasg%b1
      hatCl    =(this%liq%gamma-1.0_WP)*this%liq%cv/one_mb1
      hatbl    =this%liq_enasg%b/one_mb1
      pp_inf0  =this%liq_enasg%pp_inf0
      gam_pinf1=this%liq%gamma*this%liq_enasg%pinf1
      ! Initial guess: pressure from p-relaxation output
      if (VF.gt.0.5_WP) then
         p=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1))
      else
         p=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
      end if
      p=max(p,p_eps)
      ! Newton iterations in ln(p**)
      do it=1,itmax
         ! T**(p**) from the saturation-constraint quadratic
         Tstar=get_Tstar(p)
         if (Tstar.le.0.0_WP) then
            deallocate(y); return
         end if
         ! Phase energies at (p**,T**)
         eLstar=this%liq%get_e_from_p_T(p=p,T=Tstar)
         eGstar=cvG*Tstar+qG
         ! Residual: total energy conservation
         fp=arhoL0*eLstar+arhoG0*eGstar-rhoe0
         if (abs(fp).lt.tol*max(abs(rhoe0),1.0_WP)) exit
         ! Derivative via finite difference in ln(p)
         pold=p
         p=pold*(1.0_WP+fd_eps)
         Ts_p=get_Tstar(p)
         if (Ts_p.le.0.0_WP) then
            p=pold; exit
         end if
         eL_p=this%liq%get_e_from_p_T(p=p,T=Ts_p)
         eG_p=cvG*Ts_p+qG
         f_p=arhoL0*eL_p+arhoG0*eG_p-rhoe0
         fpp=pold*(f_p-fp)/(pold*fd_eps)
         p=pold
         if (abs(fpp).lt.tiny(1.0_WP)) exit
         dlnp=-fp/fpp
         if (abs(dlnp).gt.0.8_WP) dlnp=sign(0.8_WP,dlnp)
         lnp=log(p)+dlnp
         p=exp(lnp)
         p=max(p,p_eps)
         if (abs(dlnp).lt.tol) exit
      end do
      ! Final T** and alpha**
      Tstar=get_Tstar(p)
      if (Tstar.le.0.0_WP) then
         deallocate(y); return
      end if
      if (p.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y); return
      end if
      ! alpha**=(arhoL)^0*v_l(p**,T**)=arhoL0/rho_l(p**,T**)
      VFeq=arhoL0/this%liq%get_rho_from_p_T(p=p,T=Tstar)
      VFeq=max(0.0_WP,min(1.0_WP,VFeq))
      ! Update conserved quantities (masses unchanged)
      VF=VFeq
      Q(3)=VF      *this%liq%get_rhoe_from_p_T(p=p,T=Tstar)
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_T(p=p,T=Tstar,y=y)
      deallocate(y)
   contains
      !> T**(p) from the saturation-constraint quadratic aT*T^2+bT*T+cT=0
      !> (linear fallback for aT~0, NASG limit pinf1=0). For aT<0 (pinf1<0),
      !> both roots are positive but only the smaller one lies inside the
      !> EOS's physical branch T<Tcap(p) (one can show Q(Tcap)=arhoL0*hatCl
      !> *p_*Tcap>0, so Tcap always falls strictly between the two roots);
      !> for aT>0, exactly one root is positive and Tcap=huge admits it.
      real(WP) function get_Tstar(p_) result(T_)
         real(WP), intent(in) :: p_
         real(WP) :: aT_,bT_,cT_,disc_,one_mb_alpha_b,Tcap_,Tlo_,Thi_
         one_mb_alpha_b=1.0_WP-arhoL0*hatbl
         aT_=arhoG0*(gammaG-1.0_WP)*cvG*gam_pinf1
         bT_=arhoG0*(gammaG-1.0_WP)*cvG*(p_+pp_inf0) &
            +arhoL0*hatCl*p_                          &
            -gam_pinf1*p_*one_mb_alpha_b
         cT_=-p_*(p_+pp_inf0)*one_mb_alpha_b
         if (abs(aT_).lt.tiny(1.0_WP)) then
            if (abs(bT_).gt.tiny(1.0_WP)) then
               T_=-cT_/bT_
            else
               T_=0.0_WP
            end if
         else
            disc_=bT_**2-4.0_WP*aT_*cT_
            if (disc_.lt.0.0_WP) then
               T_=0.0_WP
            else
               Tlo_=(-bT_-sqrt(disc_))/(2.0_WP*aT_)
               Thi_=(-bT_+sqrt(disc_))/(2.0_WP*aT_)
               if (Tlo_.gt.Thi_) call swap(Tlo_,Thi_)
               Tcap_=relax_enasg_ig_Tcap(this,p_)
               if (Tlo_.gt.0.0_WP.and.Tlo_.lt.Tcap_) then
                  T_=Tlo_
               else if (Thi_.gt.0.0_WP.and.Thi_.lt.Tcap_) then
                  T_=Thi_
               else
                  T_=0.0_WP
               end if
            end if
         end if
      end function get_Tstar
      pure subroutine swap(a,b)
         real(WP), intent(inout) :: a,b
         real(WP) :: tmp
         tmp=a; a=b; b=tmp
      end subroutine swap
   end subroutine relax_enasg_ig_relax_pT

   ! ===========================================================================
   ! Chemical (pTg) relaxation helpers
   ! ===========================================================================

   !> Equilibrium temperature from the saturation constraint (volume conservation)
   !> at fixed (p, Yv) for the general LVG case.
   !>
   !> The ENASG v_l(p,T) is nonlinear in T through p'_inf(T), so the
   !> saturation constraint yields a quadratic in T (eq. ptg_Tquad in notes).
   !> Coefficients:
   !>   a_T*T^2+b_T*T+c_T=0
   !> with the physically admissible positive root returned.
   !> In the NASG limit (pinf1=0, b1=0) the quadratic degenerates to linear.
   real(WP) function relax_enasg_ig_get_T_lvg(this,p_,Yv_,rho0,rhoA0) result(T)
      class(relax_enasg_ig), intent(in) :: this
      real(WP),              intent(in) :: p_,Yv_,rho0,rhoA0
      real(WP) :: rhoLm,rhoGm,cvG_,gammaG_,one_mb1,hatCl,hatbl
      real(WP) :: aT,bT,cT,disc,T1,T2,Tlo,Thi,Tcap
      real(WP) :: pp_inf0,gam_pinf1
      ! Liquid and gas partial densities from conservation
      rhoLm =rho0*(1.0_WP-Yv_)-rhoA0    ! (alpha_l*rho_l)^0
      rhoGm =rhoA0                      ! (alpha_g*rho_g)^0
      if (rhoLm.le.0.0_WP) then
         T=0.0_WP; return
      end if
      ! Gas mixture parameters at Yv_
      cvG_  =Yv_*this%gas%get_species_cv(this%indV)+(1.0_WP-Yv_)*this%gas%get_species_cv(this%indA)
      gammaG_= (Yv_*this%gas%get_species_cp(this%indV)+(1.0_WP-Yv_)*this%gas%get_species_cp(this%indA))/cvG_
      ! ENASG shorthands
      one_mb1 =1.0_WP-this%liq_enasg%b1
      hatCl   =(this%liq%gamma-1.0_WP)*this%liq%cv/one_mb1
      hatbl   =this%liq_enasg%b/one_mb1
      pp_inf0 =this%liq_enasg%pp_inf0
      gam_pinf1= this%liq%gamma*this%liq_enasg%pinf1
      ! Quadratic coefficients (eq. ptg_Tquad in derivation notes)
      ! a_T=rhoGm*(gammaG-1)*cvG*gamma_l*pinf1/p
      aT=rhoGm*(gammaG_-1.0_WP)*cvG_*gam_pinf1/p_
      ! b_T=rhoLm*hatCl+rhoGm*(gammaG-1)*cvG*(p+pp_inf0)/p
      !      -gamma_l*pinf1*(1-rhoLm*hatbl)
      bT=rhoLm*hatCl &
        +rhoGm*(gammaG_-1.0_WP)*cvG_*(p_+pp_inf0)/p_ &
        -gam_pinf1*((1.0_WP-Yv_)-rhoLm*hatbl)
      ! c_T=-(p+pp_inf0)*((1-Yv)-rhoLm*hatbl)
      cT=-(p_+pp_inf0)*((1.0_WP-Yv_)-rhoLm*hatbl)
      ! Solve
      if (abs(aT).lt.tiny(1.0_WP)) then
         ! Linear case (NASG limit or pinf1=0)
         if (abs(bT).gt.tiny(1.0_WP)) then
            T=-cT/bT
         else
            T=0.0_WP
         end if
      else
         disc=bT**2-4.0_WP*aT*cT
         if (disc.lt.0.0_WP) then
            T=0.0_WP
         else
            ! For ENASG (pinf1<0, aT<0), both roots are positive; only the
            ! one inside (0,Tcap) lies in the EOS's physical branch
            ! p+p'_inf(T)>0 (same selection as get_Tstar in relax_pT).
            T1=(-bT-sqrt(disc))/(2.0_WP*aT)
            T2=(-bT+sqrt(disc))/(2.0_WP*aT)
            Tlo=min(T1,T2); Thi=max(T1,T2)
            Tcap=relax_enasg_ig_Tcap(this,p_)
            if (Tlo.gt.0.0_WP.and.Tlo.lt.Tcap) then
               T=Tlo
            else if (Thi.gt.0.0_WP.and.Thi.lt.Tcap) then
               T=Thi
            else
               T=0.0_WP
            end if
         end if
      end if
   end function relax_enasg_ig_get_T_lvg

   !> LV energy-conservation residual at fixed (p,T,rho0,rhoe0):
   !>   (rho0-rho_v)*rho_l*e_l+(rho_l-rho0)*rho_v*e_v-rhoe0*(rho_l-rho_v)
   !> rho_l,e_l are the exact ENASG liquid EOS; rho_v=p/(Rv*T), e_v=cvV*T+qV
   !> for the ideal-gas vapor.
   real(WP) function relax_enasg_ig_Phi_lv(this,p_,T_,rho0,rhoe0,cvV_,qV_,Rv) result(Phi)
      class(relax_enasg_ig), intent(in) :: this
      real(WP),              intent(in) :: p_,T_,rho0,rhoe0,cvV_,qV_,Rv
      real(WP) :: rho_l,e_l,rho_v,e_v
      rho_l=this%liq%get_rho_from_p_T(p=p_,T=T_)
      e_l  =this%liq%get_e_from_p_T(p=p_,T=T_)
      rho_v=p_/(Rv*T_)
      e_v  =cvV_*T_+qV_
      Phi=(rho0-rho_v)*rho_l*e_l+(rho_l-rho0)*rho_v*e_v-rhoe0*(rho_l-rho_v)
   end function relax_enasg_ig_Phi_lv

   !> Equilibrium temperature T(p) for the pure liquid-vapor case (Yv=1, no
   !> inert gas) at fixed (rho0,rhoe0): the root of relax_enasg_ig_Phi_lv,
   !> found by Newton's method with a finite-difference derivative. T is
   !> clamped to Tcap(p) (see relax_enasg_ig_Tcap) so that the exact ENASG
   !> EOS calls inside Phi_lv stay in their physical branch. The initial
   !> guess assumes the (typically tiny) liquid mass fraction does not
   !> affect the energy balance: e_v(T0)=rhoe0/rho0.
   real(WP) function relax_enasg_ig_Teq_lv(this,p_,rho0,rhoe0,cvV_,qV_,Rv) result(T)
      class(relax_enasg_ig), intent(in) :: this
      real(WP),              intent(in) :: p_,rho0,rhoe0,cvV_,qV_,Rv
      real(WP), parameter :: tol=1.0e-12_WP,fd_eps=1.0e-7_WP
      integer,  parameter :: itmax=50
      real(WP) :: Tcap_,Phi,dPhidT,dT,h,scale_
      integer  :: it
      Tcap_=0.999_WP*relax_enasg_ig_Tcap(this,p_)
      T=min(max((rhoe0/rho0-qV_)/cvV_,1.0_WP),Tcap_)
      scale_=abs(rho0*rhoe0)
      do it=1,itmax
         Phi=relax_enasg_ig_Phi_lv(this,p_,T,rho0,rhoe0,cvV_,qV_,Rv)
         if (abs(Phi).lt.tol*scale_) exit
         h=fd_eps*T
         dPhidT=(relax_enasg_ig_Phi_lv(this,p_,T+h,rho0,rhoe0,cvV_,qV_,Rv)-Phi)/h
         if (abs(dPhidT).lt.tiny(1.0_WP)) exit
         dT=-Phi/dPhidT
         T=min(max(T+dT,1.0_WP),Tcap_)
      end do
   end function relax_enasg_ig_Teq_lv

   !> Coefficients for the equilibrium temperature T(p) in the pure
   !> liquid-vapor case (Yv=1, no inert gas), used by solve_lv via
   !> T_eq=get_T_lv(ap,bp,dp) and dT/dp=get_dTdp_lv(...).
   !>
   !> The exact LV energy-conservation constraint (relax_enasg_ig_Phi_lv=0)
   !> is cubic in T for ENASG (rho_l,e_l from the exact ENASG EOS have a
   !> p'_inf(T)=gamma*pinf1*T+pp_inf0 term linear in T in their denominator),
   !> so it does not reduce to a quadratic. Rather than approximate, T_eq(p)
   !> is found by relax_enasg_ig_Teq_lv, and (ap,bp,dp) are constructed so
   !> that get_T_lv/get_dTdp_lv (which assume ap*T^2+bp*T+dp=0) recover T_eq
   !> and dT_eq/dp exactly: the quadratic (T-T_eq)*(T-Tref)=0 has roots T_eq
   !> and Tref=1K, and get_T_lv picks whichever root is closer to the host T
   !> (always T_eq, since T_eq is O(300-500K) while Tref=1K is not).
   subroutine relax_enasg_ig_get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG, &
                                            ap,bp,dp,dapdp,dbpdp,ddpdp)
      class(relax_enasg_ig), intent(in)  :: this
      real(WP),              intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP),              intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      real(WP), parameter :: Tref=1.0_WP,fd_eps=1.0e-7_WP
      real(WP) :: Rv,Teq,Teq_p,h
      Rv=(GammaG-1.0_WP)*cvG
      Teq  =relax_enasg_ig_Teq_lv(this,p_eq,rho0,rhoe0,cvG,qG,Rv)
      h    =fd_eps*p_eq
      Teq_p=(relax_enasg_ig_Teq_lv(this,p_eq+h,rho0,rhoe0,cvG,qG,Rv)-Teq)/h
      ap=1.0_WP
      bp=-(Teq+Tref)
      dp=Teq*Tref
      dapdp=0.0_WP
      dbpdp=-Teq_p
      ddpdp= Teq_p*Tref
   end subroutine relax_enasg_ig_get_coeffs_lv

   ! ===========================================================================
   ! ENASG p-T saturation curve (Gibbs-free-energy matching + Maxwell relations)
   ! ===========================================================================

   !> Largest temperature at which the ENASG liquid EOS remains in its
   !> physical branch at pressure pl_, i.e. pl_+p'_inf(T)>0. p'_inf(T)=
   !> gamma*pinf1*T+pp_inf0 is decreasing in T only when pinf1<0, in which
   !> case it crosses zero at Tcap; for pinf1>=0 (NASG/SG limit) the
   !> denominator never vanishes and Tcap=huge. A small safety margin keeps
   !> get_g_from_p_T/get_h_from_p_T/get_rho_from_p_T finite when get_Tsat's
   !> bracket search probes T far outside the saturation region.
   real(WP) function relax_enasg_ig_Tcap(this,pl_) result(Tcap)
      class(relax_enasg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_
      real(WP), parameter :: safety=0.999_WP
      real(WP) :: gam_pinf1
      gam_pinf1=this%liq%gamma*this%liq_enasg%pinf1
      if (gam_pinf1.lt.0.0_WP) then
         Tcap=safety*(pl_+this%liq_enasg%pp_inf0)/(-gam_pinf1)
      else
         Tcap=huge(1.0_WP)
      end if
   end function relax_enasg_ig_Tcap

   !> p-T saturation residual: [g_l(pl,T)-g_v(pv,T)]/(T*RV), with g_l the
   !> exact ENASG liquid Gibbs free energy and g_v the ideal-gas vapor Gibbs
   !> free energy. Reduces exactly to the inherited AS+(BS+ES*pl)/T+CS*ln(T)
   !> +DS*ln(pl+pinf)-ln(pv) form when liq is SG/NASG. T is clamped to Tcap
   !> so that get_Tsat's wide bracket search never evaluates g_l outside the
   !> EOS's physical branch.
   real(WP) function relax_enasg_ig_pTsat(this,pl_,pv_,T_) result(Fsat)
      class(relax_enasg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,pv_,T_
      real(WP) :: cpV,cvV,RV,qV,qpV,gv,Teff
      cpV=this%gas%get_species_cp(this%indV)
      cvV=this%gas%get_species_cv(this%indV)
      qV =this%gas%get_species_q(this%indV)
      qpV=this%gas%get_species_qp(this%indV)
      RV =cpV-cvV
      Teff=min(T_,relax_enasg_ig_Tcap(this,pl_))
      gv =(cpV-qpV)*Teff-cpV*Teff*log(Teff)+RV*Teff*log(pv_)+qV
      Fsat=(this%liq%get_g_from_p_T(p=pl_,T=Teff)-gv)/(Teff*RV)
   end function relax_enasg_ig_pTsat

   !> d(pTsat)/dT at fixed (pl,pv): (h_v(T)-h_l(pl,T))/(T^2*RV). The
   !> ideal-gas vapor enthalpy depends only on T, so this is exactly
   !> pv-independent. T is clamped to Tcap, see relax_enasg_ig_pTsat.
   real(WP) function relax_enasg_ig_dpTsatdT(this,pl_,T_) result(dFdT)
      class(relax_enasg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,T_
      real(WP) :: cpV,cvV,RV,qV,hv,Teff
      cpV=this%gas%get_species_cp(this%indV)
      cvV=this%gas%get_species_cv(this%indV)
      qV =this%gas%get_species_q(this%indV)
      RV =cpV-cvV
      Teff=min(T_,relax_enasg_ig_Tcap(this,pl_))
      hv =cpV*Teff+qV
      dFdT=(hv-this%liq%get_h_from_p_T(p=pl_,T=Teff))/(Teff**2*RV)
   end function relax_enasg_ig_dpTsatdT

   !> d(pTsat)/dp for the LV solver (pv_=pl_). T is clamped to Tcap, see
   !> relax_enasg_ig_pTsat.
   real(WP) function relax_enasg_ig_dpTsatdp_lv(this,p_,T_,dTdp_) result(dFdp)
      class(relax_enasg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdp_
      real(WP) :: cpV,cvV,RV,vl,Teff
      cpV=this%gas%get_species_cp(this%indV)
      cvV=this%gas%get_species_cv(this%indV)
      RV =cpV-cvV
      Teff=min(T_,relax_enasg_ig_Tcap(this,p_))
      vl =1.0_WP/this%liq%get_rho_from_p_T(p=p_,T=Teff)
      dFdp=this%dpTsatdT(p_,T_)*dTdp_+vl/(Teff*RV)-1.0_WP/p_
   end function relax_enasg_ig_dpTsatdp_lv

   !> d(pTsat)/dlnp at fixed Yv, with T=T(p,Yv). T is clamped to Tcap, see
   !> relax_enasg_ig_pTsat.
   real(WP) function relax_enasg_ig_dpTsatdlnp(this,p_,T_,dTdlnp_) result(dFdlnp)
      class(relax_enasg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdlnp_
      real(WP) :: cpV,cvV,RV,vl,Teff
      cpV=this%gas%get_species_cp(this%indV)
      cvV=this%gas%get_species_cv(this%indV)
      RV =cpV-cvV
      Teff=min(T_,relax_enasg_ig_Tcap(this,p_))
      vl =1.0_WP/this%liq%get_rho_from_p_T(p=p_,T=Teff)
      dFdlnp=this%dpTsatdT(p_,T_)*dTdlnp_+p_*vl/(Teff*RV)-1.0_WP
   end function relax_enasg_ig_dpTsatdlnp

   !> Saturated vapor pressure at given T and p_l: invert pTsat(pl,pv,T)=0
   !> for pv by evaluating pTsat at pv=1 and exponentiating.
   real(WP) function relax_enasg_ig_get_pvsat(this,pl_,T_) result(pvsat)
      class(relax_enasg_ig), intent(inout) :: this
      real(WP), intent(in)  :: pl_,T_
      pvsat=exp(this%pTsat(pl_,1.0_WP,T_))
   end function relax_enasg_ig_get_pvsat

end module relax_enasg_ig_class
