!> SG liquid and ideal gas relaxation model
module relax_sg_ig_class
   use precision,   only: WP
   use relax_class, only: relax
   use sg_class,    only: sg
   use igmix_class, only: igmix
   implicit none
   private

   public :: relax_sg_ig,Mv,Ma

   !> Molar mass of vapor and air [kg/mol]
   real(WP), parameter :: Mv=0.0180153_WP
   real(WP), parameter :: Ma=0.02897_WP

   type, extends(relax) :: relax_sg_ig
      !> Typed EOS pointers (liq is class(sg) so nasg IS-A sg works for child class)
      class(sg),    pointer :: liq => null()
      class(igmix), pointer :: gas => null()
      !> Species indices for vapor and air in the gas mixture
      integer :: indV=0,indA=0
      !> Saturation curve coefficients
      real(WP) :: AS=0.0_WP,BS=0.0_WP,CS=0.0_WP,DS=0.0_WP,ES=0.0_WP
      !> Convergence tolerances
      real(WP) :: p_tol     =1.0e-5_WP
      real(WP) :: Yv_tol    =1.0e-5_WP
      real(WP) :: Yv_tol_abs=1.0e-8_WP
      real(WP) :: Tsat_tol  =1.0e-5_WP
      real(WP) :: rho_tol   =1.0e-5_WP
      real(WP) :: rhoe_tol  =1.0e-5_WP
      real(WP) :: F1_tol    =1.0e-5_WP
      real(WP) :: F2_tol    =1.0e-5_WP
      !> Iteration limits
      integer :: Tsat_itmax=40
      integer :: NR_itmax  =40
   contains
      procedure :: initialize   =>relax_sg_ig_initialize
      procedure :: relax_p      =>relax_sg_ig_relax_p
      procedure :: relax_pT     =>relax_sg_ig_relax_pT
      procedure :: relax_pTg    =>relax_sg_ig_relax_pTg
      procedure :: get_T_lvg    =>relax_sg_ig_get_T_lvg
      procedure :: get_coeffs_lv=>relax_sg_ig_get_coeffs_lv
      procedure :: pTsat        =>relax_sg_ig_pTsat
      procedure :: dpTsatdT     =>relax_sg_ig_dpTsatdT
      procedure :: dpTsatdp_lv  =>relax_sg_ig_dpTsatdp_lv
      procedure :: dpTsatdlnp   =>relax_sg_ig_dpTsatdlnp
      procedure :: get_Tsat     =>relax_sg_ig_get_Tsat
      procedure :: get_pvsat    =>relax_sg_ig_get_pvsat
      procedure :: get_xv       =>relax_sg_ig_get_xv
   end type relax_sg_ig

contains

   !> Initialize: store EOS pointers and compute saturation curve coefficients
   subroutine relax_sg_ig_initialize(this,liq,gas,indV,indA)
      class(relax_sg_ig), intent(inout)    :: this
      class(sg),    target, intent(in)     :: liq
      class(igmix), target, intent(in)     :: gas
      integer, intent(in) :: indV,indA
      real(WP) :: cpV,cvV,RV
      this%liq=>liq
      this%gas=>gas
      this%indV=indV
      this%indA=indA
      cpV=gas%get_species_cp(indV)
      cvV=gas%get_species_cv(indV)
      RV=cpV-cvV
      this%AS=(liq%cp-cpV+gas%get_species_qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%get_species_q(indV))            /RV
      this%CS=(cpV-liq%cp)                                /RV
      this%DS=(liq%cp-liq%cv)                             /RV
      ! ES stays 0.0_WP (SG has no co-volume b); NASG overrides this after calling parent init
   end subroutine relax_sg_ig_initialize

   !> Mechanical relaxation
   subroutine relax_sg_ig_relax_p(this,VF,Q,Pjump)
      class(relax_sg_ig),      intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: y
      real(WP) :: PL,PG,ZL,ZG,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Get phasic pressures
      PL=this%liq%get_p_from_rho_e(Q(1)/(       VF),Q(3)/Q(1))
      PG=this%gas%get_p_from_rho_e(Q(2)/(1.0_WP-VF),Q(4)/Q(2),y)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(8)=Yv*Q(2)
         deallocate(y)
         return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP
         deallocate(y)
         return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(PL,Q(1)/(       VF))**2
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(PG,Q(2)/(1.0_WP-VF),y)**2
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup quadratic problem
      coeffL=(this%liq%gamma-1.0_WP)*Pint+2.0_WP*this%liq%gamma*this%liq%pinf
      coeffG=(gammaG        -1.0_WP)*Pint
      a=1.0_WP+gammaG*VF+this%liq%gamma*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+gammaG)*VF*PL-(1.0_WP+this%liq%gamma)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Get equilibrium volume fraction
      VFeq=VF*((this%liq%gamma-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+this%liq%gamma)*Peq+coeffL)
      ! Adjust conserved quantities
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      deallocate(y)
   end subroutine relax_sg_ig_relax_p

   !> Mechanical and thermal relaxation
   subroutine relax_sg_ig_relax_pT(this,VF,Q,Pjump)
      class(relax_sg_ig), intent(inout)      :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:), allocatable :: y
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! ================ First step for mechanical relaxation ================
      call this%relax_p(VF,Q,Pjump)
      ! ================= Second step for thermal relaxation =================
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Setup quadratic problem
      a=Q(1)*this%liq%cv+Q(2)*cvG
      b=this%liq%q*this%liq%cv*(this%liq%gamma-1.0_WP)*Q(1)**2+qG*cvG*(gammaG-1.0_WP)*Q(2)**2+&
      &  Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf+Q(2)*cvG*this%liq%pinf                +&
      &  Q(1)*Q(2)*(this%liq%q*cvG*(gammaG-1.0_WP)+qG*this%liq%cv*(this%liq%gamma-1.0_WP))   -&
      &  sum(Q(3:4))*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*cvG*(gammaG-1.0_WP))
      d=cvG*(gammaG-1.0_WP)*this%liq%pinf*(qG*Q(2)**2+this%liq%q*Q(1)*Q(2)-sum(Q(3:4))*Q(2))
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y)
         return
      end if
      ! Get equilibrium volume fraction
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq/(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*cvG*(gammaG-1.0_WP)*(Peq+this%liq%pinf))
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then
         VFeq=0.0_WP
         Peq=max(Peq,-this%liq%pinf)
      end if
      if (VFeq.gt.1.0_WP) then
         VFeq=1.0_WP
         Peq=max(Peq,0.0_WP)
      end if
      ! Adjust conserved quantities
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(Peq,Q(1)/max(VFeq,tiny(1.0_WP)))
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(Peq,Q(2)/max(1.0_WP-VFeq,tiny(1.0_WP)),y)
      VF=VFeq
      deallocate(y)
   end subroutine relax_sg_ig_relax_pT

   !> Mechanical, thermal, and chemical relaxation
   subroutine relax_sg_ig_relax_pTg(this,VF,Q,Pjump)
      class(relax_sg_ig),      intent(inout) :: this
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,Qin,y
      real(WP) :: VF0,VFin,p,T,Yv
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: RHOL,RHOG
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP), parameter :: p_eps=1.0e-10_WP,VFmin=1.0e-5_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Yv_dry=1.0e-5_WP,pv_dry=1.0_WP,Yv_pure=0.999_WP
      real(WP), parameter :: fd_eps=1.0e-7_WP,F_line_search_tol=0.3_WP
      logical :: chem_relax,nucleated
      allocate(Qin(size(Q)))
      Qin=Q
      VFin=VF
      nucleated=.false.
      ! Nucleation: Conservatively move a little mass and energy so the chemical relaxation starts from a non-stiff initial condition.
      ! This handles both cavitation and condensation.
      nucleation: block
         real(WP), parameter :: VF_nuc=1.0e-7_WP
         real(WP) :: rhoL_nuc,pL_nuc,TL_nuc,pv_sat,rhoV_nuc,eV_nuc
         real(WP) :: rhoG_nuc,pG_nuc,TG_nuc,Yv_nuc,xv_nuc,pv_nuc,Tsat_nuc
         real(WP) :: rhoL_new,eL_new,drho,de
         real(WP) :: y_nuc(this%gas%ns)
         logical  :: conv_nuc
         integer  :: Tsat_it_nuc
         ! Almost pure liquid: check if liquid is metastable and needs vapor nucleation
         if (VF.ge.1.0_WP-VF_nuc) then
            ! Get current liquid state
            rhoL_nuc=Q(1)/VF
            pL_nuc=this%liq%get_p_from_rho_e(rhoL_nuc,Q(3)/Q(1))
            TL_nuc=this%liq%get_T_from_p_rho(pL_nuc,rhoL_nuc)
            ! Check if inside EOS validity range
            if (pL_nuc.le.-this%liq%pinf.or.TL_nuc.le.0.0_WP) return
            ! Saturation vapor pressure at current liquid state
            pv_sat=this%get_pvsat(pL_nuc,TL_nuc)
            ! Check if metastable
            if (pv_sat.le.pL_nuc) return
            ! Nucleate a tiny vapor phase
            y_nuc           =0.0_WP
            y_nuc(this%indV)=1.0_WP
            ! Approximate vapor state
            rhoV_nuc=this%gas%get_rho_from_p_T(pv_sat,TL_nuc,y_nuc)
            eV_nuc=this%gas%get_e_from_p_T(pv_sat,TL_nuc,y_nuc)
            ! Transfer mass and energy from liquid to vapor
            drho=VF_nuc*rhoV_nuc
            de=drho*eV_nuc
            Q(1)=Q(1)-drho; Q(2)=Q(2)+drho
            Q(3)=Q(3)-de;   Q(4)=Q(4)+de
            Q(8)=Q(8)+drho
            VF=1.0_WP-VF_nuc
            nucleated=.true.
         ! Almost pure gas: check if vapor is metastable and needs liquid nucleation.
         else if (VF.le.VF_nuc) then
            ! Get vapor mass fraction
            Yv_nuc=Q(8)/Q(2)
            Yv_nuc=max(Yvmin,min(Yvmax,Yv_nuc))
            ! Return if not enough vapor to condense
            if (Yv_nuc.le.Yv_dry) return
            ! Get gas composition and state
            y_nuc(this%indV)=Yv_nuc
            y_nuc(this%indA)=1.0_WP-Yv_nuc
            rhoG_nuc=Q(2)/max(1.0_WP-VF,tiny(1.0_WP))
            pG_nuc=this%gas%get_p_from_rho_e(rhoG_nuc,Q(4)/Q(2),y_nuc)
            TG_nuc=this%gas%get_T_from_p_rho(pG_nuc,rhoG_nuc,y_nuc)
            ! Return if not physical
            if (pG_nuc.le.0.0_WP.or.TG_nuc.le.0.0_WP) return
            ! Get vapor partial pressure
            xv_nuc=this%get_xv(Yv_nuc)
            pv_nuc=xv_nuc*pG_nuc
            if (.not.check_pv(pv_nuc)) return
            ! Calculate saturation temperature
            call this%get_Tsat(pG_nuc,pv_nuc,TG_nuc,Tsat_nuc,conv_nuc,Tsat_it_nuc)
            if (.not.conv_nuc) return
            ! Check if metastable
            if (TG_nuc.ge.Tsat_nuc) return
            ! Nucleate a tiny liquid phase
            rhoL_new=this%liq%get_rho_from_p_T(pG_nuc,TG_nuc)
            eL_new=this%liq%get_e_from_p_T(pG_nuc,TG_nuc)
            ! Transfer mass and energy from liquid to vapor
            drho=VF_nuc*rhoL_new
            drho=min(drho,0.5_WP*Q(8),0.5_WP*Q(2))
            if (drho.le.0.0_WP) return
            de=drho*eL_new
            Q(1)=Q(1)+drho; Q(2)=Q(2)-drho
            Q(3)=Q(3)+de;   Q(4)=Q(4)-de
            Q(8)=Q(8)-drho
            VF=drho/rhoL_new
            nucleated=.true.
         end if
      end block nucleation
      ! ================ First and second steps for mechanical and thermal relaxation ================
      call this%relax_pT(VF,Q,Pjump)
      ! ================= Third step for chemical relaxation =================
      ! Set gas EOS coefficients from vapor mass fraction (Q(8), Q(2) unchanged by relax_pT)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Recover p and T from the dominant phase for better numerical stability
      if (VF.gt.0.5_WP) then
         p=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1))
         T=this%liq%get_T_from_p_rho(p=p,rho=Q(1)/VF)
      else
         p=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
         T=this%gas%get_T_from_p_rho(p=p,rho=Q(2)/(1.0_WP-VF),y=y)
      end if
      ! Store input state to the chemical relaxation algorithm
      VF0=VF
      allocate(Q0(size(Q)))
      Q0=Q
      rho0 =sum(Q0(1:2))
      rhoe0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      ! Check whether a stable pure phase state exists
      pure_phase_bounds: block
         real(WP), parameter :: rhoA_pure=1.0e-12_WP
         real(WP) :: rhoL_pure,eL_pure,pL_pure,TL_pure,Tsat_pure
         real(WP) :: rhoV_pure,eV_pure,pV_pure,TV_pure
         real(WP) :: rhoV0,Yv_gas,pG_pure,TG_pure,xv_gas,pv_gas
         integer  :: Tsat_it_pure
         logical  :: conv_pure
         if (rho0.le.0.0_WP.or.rhoe0.le.0.0_WP) exit pure_phase_bounds
         ! Try pure liquid: rho_l = rho0, e_l = rhoe0/rho0.
         if (rhoA0/rho0.le.rhoA_pure) then
            rhoL_pure=rho0
            eL_pure=rhoe0/rho0
            pL_pure=this%liq%get_p_from_rho_e(rhoL_pure,eL_pure)
            TL_pure=this%liq%get_T_from_p_rho(pL_pure,rhoL_pure)
            if (pL_pure.gt.p_eps.and.TL_pure.gt.0.0_WP) then
               call this%get_Tsat(pL_pure,pL_pure,TL_pure,Tsat_pure,conv_pure,Tsat_it_pure)
               if (conv_pure.and.TL_pure.le.Tsat_pure*(1.0_WP+this%Tsat_tol)) then
                  VF=1.0_WP
                  Q(1)=rho0;  Q(2)=0.0_WP
                  Q(3)=rhoe0; Q(4)=0.0_WP
                  Q(8)=0.0_WP
                  call dealloc()
                  return
               end if
            end if
            ! Try pure vapor: rho_v = rho0, e_v = rhoe0/rho0.
            rhoV_pure=rho0
            eV_pure=rhoe0/rho0
            y=0.0_WP
            y(this%indV)=1.0_WP
            pV_pure=this%gas%get_p_from_rho_e(rhoV_pure,eV_pure,y)
            TV_pure=this%gas%get_T_from_p_rho(pV_pure,rhoV_pure,y)
            if (pV_pure.gt.p_eps.and.TV_pure.gt.0.0_WP) then
               call this%get_Tsat(pV_pure,pV_pure,TV_pure,Tsat_pure,conv_pure,Tsat_it_pure)
               if (conv_pure.and.TV_pure.ge.Tsat_pure*(1.0_WP-this%Tsat_tol)) then
                  VF=0.0_WP
                  Q(1)=0.0_WP;   Q(2)=rho0
                  Q(3)=0.0_WP;   Q(4)=rhoe0
                  Q(8)=rho0
                  call dealloc()
                  return
               end if
            end if
         else
            ! With non-condensable gas
            rhoV0=rho0-rhoA0
            if (rhoV0.le.0.0_WP) exit pure_phase_bounds
            Yv_gas=rhoV0/rho0
            Yv_gas=max(Yvmin,min(Yvmax,Yv_gas))
            y=0.0_WP
            y(this%indV)=Yv_gas
            y(this%indA)=1.0_WP-Yv_gas
            pG_pure=this%gas%get_p_from_rho_e(rho0,rhoe0/rho0,y)
            TG_pure=this%gas%get_T_from_p_rho(pG_pure,rho0,y)
            if (pG_pure.gt.p_eps.and.TG_pure.gt.0.0_WP) then
               xv_gas=this%get_xv(Yv_gas)
               pv_gas=xv_gas*pG_pure
               if (check_pv(pv_gas)) then
                  call this%get_Tsat(pG_pure,pv_gas,TG_pure,Tsat_pure,conv_pure,Tsat_it_pure)
                  if (conv_pure.and.TG_pure.ge.Tsat_pure*(1.0_WP-this%Tsat_tol)) then
                     VF=0.0_WP
                     Q(1)=0.0_WP;   Q(2)=rho0
                     Q(3)=0.0_WP;   Q(4)=rhoe0
                     Q(8)=rhoV0
                     call dealloc()
                     return
                  end if
               end if
            end if
         end if
      end block pure_phase_bounds
      ! Check if chemical relaxation should be activated
      chem_relax=activate_chem(p,T,Yv)
      if (.not.chem_relax) then
         call restore()
         call dealloc()
         return
      end if
      ! Solve chemical equilibrium
      if (Yv.gt.Yv_pure) then
         Yv=Yvmax
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
         call solve_lv(p,T,chem_relax)
      else
         call solve_lvg(p,T,Yv,chem_relax)
      end if
      ! Skip if not converged
      if (.not.chem_relax) then
         if (nucleated) then
            VF=VFin
            Q=Qin
         else
            call restore()
         end if
         call dealloc()
         return
      end if
      ! Adjust vapor mass fraction and update gas EOS parameters
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Adjust densities
      RHOL=this%liq%get_rho_from_p_T(p=p,T=T)
      RHOG=this%gas%get_rho_from_p_T(p=p,T=T,y=y)
      ! Adjust VF
      VF=(rho0-RHOG)/(RHOL-RHOG)
      ! Clean up solution
      if (VF.lt.0.0_WP) then
         VF=0.0_WP
         p=max(p,-this%liq%pinf)
      end if
      if (VF.gt.1.0_WP) then
         VF=1.0_WP
         p=max(p,0.0_WP)
      end if
      ! Adjust conserved quantities (velocities remain unchanged)
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*this%liq%get_e_from_p_T(p=p,T=T)
      Q(4)=Q(2)*this%gas%get_e_from_p_T(p=p,T=T,y=y)
      Q(8)=Q(2)*Yv
      ! Evaluate conservation
      if (.not.check_cons()) then
         call restore()
         call dealloc()
         return
      end if
      ! Release memory
      call dealloc()
   contains
      !> Reset the output to the initial values fed into chemical relaxation algorithm
      subroutine restore()
         VF=VF0
         Q=Q0
      end subroutine restore
      !> Release memory allocated for Q0 and y
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
      end subroutine dealloc
      !> Check sanity of a vapor partial pressure
      logical function check_pv(pv_)
         real(WP), intent(in) :: pv_
         check_pv=(pv_.gt.p_eps)
      end function check_pv
      !> Check mass and energy conservation
      logical function check_cons()
         real(WP) :: re,ee
         re=(sum(Q(1:2))-rho0)/rho0
         ee=(sum(Q(3:4))-rhoe0)/rhoe0
         check_cons=(abs(re).le.this%rho_tol).and.(abs(ee).le.this%rhoe_tol)
      end function check_cons
      !> Decide if chemical relaxation needs to be activated and adjust Yv if needed
      logical function activate_chem(p_,T_,Yv_)
         real(WP), intent(in)    :: p_,T_
         real(WP), intent(inout) :: Yv_
         real(WP) :: xv,pv_,Fsat
         activate_chem=.false.
         ! Get vapor mole fraction and partial pressure
         xv=this%get_xv(Yv_)
         pv_=xv*p_
         if (Yv_.le.Yv_dry) then
            ! Dry/nearly-dry air edge case: pv_ is zero or so tiny that solving Tsat(pv_) is log-singular/ill-conditioned.
            ! Seed Yv from saturation at the current thermally-relaxed state and then continue with the ordinary LVG NR.
            pv_=exp(this%AS+(this%BS+this%ES*p_)/T_)*T_**this%CS*(p_+this%liq%pinf)**this%DS
            if (.not.check_pv(pv_)) then
               return
            end if
            if (pv_.ge.p_) then
               ! pv_sat >= p: flash regime. Seed with a small Yv so that get_T stays well-conditioned.
               Yv_=0.01_WP
            else
               xv=pv_/p_
               Yv_=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
            end if
            Yv_=max(Yvmin,min(Yvmax,Yv_))
         else if ((pv_.le.pv_dry).or.(.not.check_pv(pv_))) then
            ! Yv_ is already known and is NOT dry (e.g. Yv_=1, no inert gas at all),
            ! yet pv_=xv*p_ comes out tiny/negative -- this happens when the *total*
            ! pressure p_ itself is small or negative (a liquid-dominated cell in
            ! tension/cavitation), not because the vapor content is small. The
            ! log-based saturation residual pTsat(p_,pv_,T_) is ill-conditioned or
            ! undefined for pv_<=0, but a physical saturation pressure pv_sat(T_) is
            ! always positive while p_<=0 here, so the mixture cannot possibly be at
            ! equilibrium: flag relaxation as needed without overwriting Yv_, which
            ! is already correct.
         else
            ! Direct saturation residual. No Tsat solve needed.
            Fsat=this%pTsat(p_,pv_,T_)
            if (abs(Fsat).lt.this%F1_tol) return
         end if
         activate_chem=.true.
      end function activate_chem
      !> Equilibrium T for LV (Both quadratic roots for T are computed)
      real(WP) function get_T_lv(ap,bp,dp)
         real(WP), intent(in) :: ap,bp,dp
         real(WP) :: T1,T2
         T1=(-bp-sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
         T2=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
         if (abs(T1-T).lt.abs(T2-T)) then
            get_T_lv=T1
         else
            get_T_lv=T2
         end if
      end function get_T_lv
      !> dT/dp from the quadratic coefficients for LV (Both quadratic roots for T are considered)
      real(WP) function get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
         real(WP), intent(in) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: T1,T2,sgn
         T1=(-bp-sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
         T2=(-bp+sqrt(bp**2-4.0_WP*ap*dp))/(2.0_WP*ap)
         if (abs(T1-T).lt.abs(T2-T)) then
            sgn=-1.0_WP
         else
            sgn= 1.0_WP
         end if
         get_dTdp_lv=(ap*(-dbpdp+sgn*(bp*dbpdp-2.0_WP*(dapdp*dp+ap*ddpdp))/sqrt(bp**2-4.0_WP*ap*dp))-dapdp*(-bp+sgn*sqrt(bp**2-4.0_WP*ap*dp)))/(2.0_WP*ap**2)
      end function get_dTdp_lv
      !> Energy residual for LVG
      real(WP) function rhoe_res_lvg(p_,T_,Yv_)
         real(WP), intent(in) :: p_,T_,Yv_
         rhoe_res_lvg=(rho0*(1.0_WP-Yv_)-rhoA0)*this%liq%get_e_from_p_T(p_,T_)+&
         &             rhoA0*this%gas%get_e_from_p_T(p_,T_,[Yv_,1.0_WP-Yv_])  -&
         &             rhoe0*(1.0_WP-Yv_)
      end function rhoe_res_lvg
      !> Energy residual for LV
      real(WP) function rhoe_res_lv(p_,T_)
         real(WP), intent(in) :: p_,T_
         real(WP) :: rho_l,rho_g
         rho_l=this%liq%get_rho_from_p_T(p=p_,T=T_)
         rho_g=this%gas%get_rho_from_p_T(p=p_,T=T_,y=[1.0_WP,0.0_WP])
         rhoe_res_lv=(rho0-rho_g)/(rho_l-rho_g)*this%liq%get_e_from_p_T(p=p_,T=T_)+                  &
         &           (rho_l-rho0)/(rho_l-rho_g)*this%gas%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP,0.0_WP])-&
         &            rhoe0
      end function rhoe_res_lv
      !> d(ln(xv))/dYv
      real(WP) function dlnxvdYv(Yv_)
         real(WP), intent(in) :: Yv_
         real(WP) :: Ys,den
         Ys=max(Yv_,Yvmin+fd_eps)
         den=Ys*Ma+(1.0_WP-Ys)*Mv
         dlnxvdYv=1.0_WP/Ys-(Ma-Mv)/den
      end function dlnxvdYv
      !> Pure liquid and vapor chemical relaxation (Solves for ln(p))
      subroutine solve_lv(p_eq,T_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq
         logical,  intent(out)   :: conv
         real(WP) :: pOld,lnpOld,p_try,T_try
         real(WP) :: ap,bp,dp,dapdp,dbpdp,ddpdp
         real(WP) :: dTdp,dTdlnp,dF1dlnp
         real(WP) :: F1,F1_try,dlnp_nr,p_err,alpha
         integer  :: it
         logical  :: accepted
         ! Reset small pressure to saturated vapor pressure
         if (p_eq.le.p_eps) p_eq=this%get_pvsat(p_eq,T_eq)
         ! Iteratively solve for the equilibrium log-pressure
         conv=.false.
         p_err=10.0_WP*this%p_tol
         do it=1,this%NR_itmax
            ! Get the coefficients
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
            ! Get temperature and dT/dlnp
            T_eq=get_T_lv(ap,bp,dp)
            if (T_eq.le.0.0_WP) return
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp)
            dTdlnp=p_eq*dTdp
            ! Get residual and derivative with respect to lnp
            F1=this%pTsat(p_eq,p_eq,T_eq)
            dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp)
            if (abs(dF1dlnp).lt.1.0e-30_WP) exit
            ! Newton-Raphson update in lnp
            dlnp_nr=-F1/dF1dlnp
            ! Damped Newton-Raphson update
            pOld=p_eq
            lnpOld=log(pOld)
            alpha=1.0_WP
            accepted=.false.
            do while (alpha.gt.1.0e-8_WP)
               p_try=exp(lnpOld+alpha*dlnp_nr)
               if (p_try.le.p_eps) then
                  alpha=0.5_WP*alpha
                  cycle
               end if
               call this%get_coeffs_lv(p_try,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
               T_try=get_T_lv(ap,bp,dp)
               if (T_try.le.0.0_WP) then
                  alpha=0.5_WP*alpha
                  cycle
               end if
               F1_try=this%pTsat(p_try,p_try,T_try)
               if (abs(F1_try).lt.abs(F1)) then
                  p_eq=p_try
                  T_eq=T_try
                  accepted=.true.
                  exit
               end if
               alpha=0.5_WP*alpha
            end do
            if (.not.accepted) then
               exit
            end if
            ! Refresh temperature and residual from the accepted pressure
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp)
            if (T_eq.le.0.0_WP) return
            ! Evaluate the error
            p_err=abs(log(p_eq/pOld))
            F1=this%pTsat(p_eq,p_eq,T_eq)
            if ((p_err.lt.this%p_tol).and.(abs(F1).lt.this%F1_tol)) then
               conv=.true.
               exit
            end if
         end do
         ! Check convergence
         if (.not.conv) then
            return
         end if
         ! Update equilibrium temperature
         call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
         T_eq=get_T_lv(ap,bp,dp)
      end subroutine solve_lv
      !> Pure liquid and vapor-gas mixture chemical relaxation (Solves for ln(p) and Yv)
      subroutine solve_lvg(p_eq,T_eq,Yv_eq,conv)
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: conv
         real(WP) :: xv,pv
         real(WP) :: F1,F2,F2p,F2Y,dF1dlnp,dF1dYv,dF2dlnp,dF2dYv,detJ
         real(WP) :: lnp_pert,p_pert,Yv_pert,T_pert,xv_pert,pv_pert,dTdlnp,dTdYv,dlnp_nr,dYv_nr
         real(WP) :: pOld,YvOld,lnpOld,p_err,Yv_err
         real(WP) :: alpha,res0,res_try
         real(WP) :: p_try,Yv_try,T_try,xv_try,pv_try,F1_try,F2_try
         real(WP) :: Yv_max_phys,Yv_hi
         integer  :: it,lsit
         logical  :: accepted
         ! Calculate the absolute physical ceiling for Yv based on available liquid
         Yv_max_phys=1.0_WP-(rhoA0/rho0)
         Yv_hi=min(Yvmax,Yv_max_phys)
         ! Iteratively solve for the equilibrium ln(p) and Yv
         conv=.false.
         p_err=10.0_WP*this%p_tol
         Yv_err=10.0_WP*this%Yv_tol
         do it=1,this%NR_itmax
            ! Get temperature
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0)
            if (T_eq.le.0.0_WP) return
            ! Get vapor partial pressure
            xv=this%get_xv(Yv_eq)
            pv=xv*p_eq
            if (.not.check_pv(pv)) return
            ! Get residuals at current state
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            res0=sqrt(F1**2+F2**2)
            ! Perturbation in lnp
            lnp_pert=log(p_eq)+fd_eps
            p_pert=exp(lnp_pert)
            pv_pert=xv*p_pert
            if (.not.check_pv(pv_pert)) return
            ! Get the corresponding T
            T_pert=this%get_T_lvg(p_pert,Yv_eq,rho0,rhoA0)
            if (T_pert.le.0.0_WP) return
            ! Get residuals and derivatives with respect to lnp
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)/rhoe0
            dTdlnp=(T_pert-T_eq)/fd_eps
            dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp)
            dF2dlnp=(F2p-F2)/fd_eps
            ! Perturbation in Yv
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yv_hi-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            if ((Yv_pert.le.Yvmin+fd_eps).or.(Yv_pert.ge.Yv_hi-fd_eps)) return
            xv_pert=this%get_xv(Yv_pert)
            pv_pert=xv_pert*p_eq
            if (.not.check_pv(pv_pert)) return
            ! Get the corresponding T
            T_pert=this%get_T_lvg(p_eq,Yv_pert,rho0,rhoA0)
            if (T_pert.le.0.0_WP) return
            ! Get residuals and derivatives with respect to Yv
            F2Y=rhoe_res_lvg(p_eq,T_pert,Yv_pert)/rhoe0
            dTdYv=(T_pert-T_eq)/(Yv_pert-Yv_eq)
            dF1dYv=this%dpTsatdT(p_eq,T_eq)*dTdYv-dlnxvdYv(Yv_eq)
            dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
            ! Solve 2x2 system: J*[dlnp; dYv]=-[F1; F2]
            detJ=dF1dlnp*dF2dYv-dF1dYv*dF2dlnp
            if (abs(detJ).lt.1.0e-30_WP) exit
            dlnp_nr=-( dF2dYv *F1-dF1dYv *F2)/detJ
            dYv_nr =-(-dF2dlnp*F1+dF1dlnp*F2)/detJ
            ! Direction preserving step limiter
            step_limit: block
               real(WP) :: ms,lnp_up,lnp_dn
               ms=1.0_WP
               ! 1. Prevent pressure from changing by more than 50% in a single step
               lnp_up=log(1.5_WP)
               lnp_dn=log(0.5_WP)
               if (dlnp_nr.gt.lnp_up) ms=min(ms,lnp_up/dlnp_nr)
               if (dlnp_nr.lt.lnp_dn) ms=min(ms,lnp_dn/dlnp_nr)
               ! 2. Prevent Yv from crossing physical boundaries
               if (dYv_nr.gt.0.0_WP) then
                  if (Yv_eq+dYv_nr.ge.Yv_hi-fd_eps) ms=min(ms,0.9_WP*(Yv_hi-fd_eps-Yv_eq)/dYv_nr)
               else if (dYv_nr.lt.0.0_WP) then
                  if (Yv_eq+dYv_nr.le.Yvmin+fd_eps) ms=min(ms,0.9_WP*(Yv_eq-Yvmin-fd_eps)/abs(dYv_nr))
               end if
               ! Scale the Newton-Raphson step uniformly to keep pointing directly at the root
               ms=max(0.0_WP,min(1.0_WP,ms))
               dlnp_nr=dlnp_nr*ms
               dYv_nr=dYv_nr*ms
            end block step_limit
            ! Damped Newton-Raphson update
            pOld=p_eq
            YvOld=Yv_eq
            lnpOld=log(pOld)
            alpha=1.0_WP
            lsit=0
            if ((abs(F1).lt.F_line_search_tol).and.(abs(F2).lt.F_line_search_tol)) then
               p_eq=exp(lnpOld+dlnp_nr)
               Yv_eq=YvOld+dYv_nr
            else
               accepted=.false.
               do while (alpha.gt.1.0e-8_WP)
                  lsit=lsit+1
                  p_try=exp(lnpOld+alpha*dlnp_nr)
                  Yv_try=YvOld+alpha*dYv_nr
                  ! Keep the trial state inside the physical/log-safe domain
                  if (p_try.le.p_eps) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  if ((Yv_try.le.Yvmin+fd_eps).or.(Yv_try.ge.Yv_hi-fd_eps)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  if ((rho0*(1.0_WP-Yv_try)-rhoA0).le.0.0_WP) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  T_try=this%get_T_lvg(p_try,Yv_try,rho0,rhoA0)
                  if (T_try.le.0.0_WP) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  xv_try=this%get_xv(Yv_try)
                  pv_try=xv_try*p_try
                  if (.not.check_pv(pv_try)) then
                     alpha=0.5_WP*alpha
                     cycle
                  end if
                  F1_try=this%pTsat(p_try,pv_try,T_try)
                  F2_try=rhoe_res_lvg(p_try,T_try,Yv_try)/rhoe0
                  res_try=sqrt(F1_try**2+F2_try**2)
                  if (res_try.lt.res0) then
                     p_eq=p_try
                     Yv_eq=Yv_try
                     T_eq=T_try
                     accepted=.true.
                     exit
                  end if
                  alpha=0.5_WP*alpha
               end do
               if (.not.accepted) then
                  exit
               end if
            end if
            ! Evaluate errors
            p_err=abs(log(p_eq/pOld))
            Yv_err=abs(Yv_eq-YvOld)
            ! Refresh xv and pv from the accepted solution
            xv=this%get_xv(Yv_eq)
            pv=xv*p_eq
            if (.not.check_pv(pv)) return
            ! Get temperature
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0)
            if (T_eq.le.0.0_WP) return
            ! Evaluate errors
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            if ((p_err.lt.this%p_tol).and.Yv_err.lt.this%Yv_tol_abs+this%Yv_tol*max(abs(YvOld),abs(Yv_eq)).and.(abs(F1).lt.this%F1_tol).and.(abs(F2).lt.this%F2_tol)) then
               conv=.true.
               exit
            end if
         end do
         ! Check convergence
         if (.not.conv) then
            return
         end if
      end subroutine solve_lvg
   end subroutine relax_sg_ig_relax_pTg

   !> p-T saturation curve
   real(WP) function relax_sg_ig_pTsat(this,pl_,pv_,T_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,pv_,T_
      relax_sg_ig_pTsat=this%AS+(this%BS+this%ES*pl_)/T_+this%CS*log(T_)+this%DS*log(pl_+this%liq%pinf)-log(pv_)
   end function relax_sg_ig_pTsat

   !> d(pTsat)/dT
   real(WP) function relax_sg_ig_dpTsatdT(this,pl_,T_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: pl_,T_
      relax_sg_ig_dpTsatdT=-(this%BS+this%ES*pl_)/T_**2+this%CS/T_
   end function relax_sg_ig_dpTsatdT

   !> d(pTsat)/dp for the LV solver (pv_ = pl_)
   real(WP) function relax_sg_ig_dpTsatdp_lv(this,p_,T_,dTdp_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdp_
      relax_sg_ig_dpTsatdp_lv=this%dpTsatdT(p_,T_)*dTdp_+this%ES/T_+this%DS/(p_+this%liq%pinf)-1.0_WP/p_
   end function relax_sg_ig_dpTsatdp_lv

   !> d(pTsat)/dlnp at fixed Yv, with T=T(p,Yv)
   real(WP) function relax_sg_ig_dpTsatdlnp(this,p_,T_,dTdlnp_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdlnp_
      relax_sg_ig_dpTsatdlnp=this%dpTsatdT(p_,T_)*dTdlnp_+this%ES*p_/T_+this%DS*p_/(p_+this%liq%pinf)-1.0_WP
   end function relax_sg_ig_dpTsatdlnp

   !> Safeguarded Newton-Raphson for saturation temperature at fixed (pl_, pv_)
   subroutine relax_sg_ig_get_Tsat(this,pl_,pv_,Tguess,Tsat,conv,Tsat_it)
      class(relax_sg_ig), intent(inout) :: this
      real(WP), intent(in)  :: pl_,pv_,Tguess
      real(WP), intent(out) :: Tsat
      logical,  intent(out) :: conv
      integer,  intent(out) :: Tsat_it
      real(WP) :: Tlo,Thi,Told,Tnew,Flo,Fhi,Fold,Fnew,dFold
      integer  :: it,expand_it
      conv=.false.
      Tsat_it=0
      ! Broad physical bracket for the current water EOS.  Newton is accepted only when it remains inside this bracket;
      ! otherwise the update falls back to bisection.
      Tlo=250.0_WP
      Thi=900.0_WP
      Flo=this%pTsat(pl_,pv_,Tlo)
      Fhi=this%pTsat(pl_,pv_,Thi)
      ! Expand the bracket if needed.  In the normal liquid-vapor range pTsat is monotone in T, so these one-sided expansions are enough.
      expand_it=0
      do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
         if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
            Tlo=max(1.0_WP,0.8_WP*Tlo)
            Flo=this%pTsat(pl_,pv_,Tlo)
         else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
            Thi=1.2_WP*Thi
            Fhi=this%pTsat(pl_,pv_,Thi)
         else
            exit
         end if
         expand_it=expand_it+1
      end do
      if (Flo*Fhi.gt.0.0_WP) then
         return
      end if
      ! Use the caller's guess only after clamping it to the safe bracket.
      Tsat=max(Tlo,min(Thi,Tguess))
      do it=1,this%Tsat_itmax
         Told=Tsat
         Fold=this%pTsat(pl_,pv_,Told)
         dFold=this%dpTsatdT(pl_,Told)
         if (abs(dFold).gt.tiny(1.0_WP)) then
            Tnew=Told-Fold/dFold
         else
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         ! Safeguard: reject Newton steps that leave the bracket or are NaN.
         if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) then
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         Fnew=this%pTsat(pl_,pv_,Tnew)
         if (Fnew.ne.Fnew) then
            Tnew=0.5_WP*(Tlo+Thi)
            Fnew=this%pTsat(pl_,pv_,Tnew)
         end if
         ! Update the bracket around the root.
         if (Flo*Fnew.le.0.0_WP) then
            Thi=Tnew
            Fhi=Fnew
         else
            Tlo=Tnew
            Flo=Fnew
         end if
         Tsat_it=it
         Tsat=Tnew
         ! Evaluate error
         if ((abs((Tnew-Told)/max(abs(Told),tiny(1.0_WP))).lt.this%Tsat_tol).or.(abs(Fnew).lt.this%F1_tol)) then
            conv=.true.
            return
         end if
      end do
   end subroutine relax_sg_ig_get_Tsat

   !> Get saturated vapor pressure at given T and p_l
   real(WP) function relax_sg_ig_get_pvsat(this,pl_,T_)
      class(relax_sg_ig), intent(inout) :: this
      real(WP), intent(in)  :: pl_,T_
      relax_sg_ig_get_pvsat=exp(this%AS+(this%BS+this%ES*pl_)/T_)*T_**this%CS*(pl_+this%liq%pinf)**this%DS
   end function relax_sg_ig_get_pvsat

   !> Vapor mole fraction from mass fraction
   real(WP) function relax_sg_ig_get_xv(this,Yv_)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: Yv_
      relax_sg_ig_get_xv=Yv_*Ma/(Yv_*Ma+(1.0_WP-Yv_)*Mv)
   end function relax_sg_ig_get_xv

   !> Equilibrium T from energy conservation
   real(WP) function relax_sg_ig_get_T_lvg(this,p_,Yv_,rho0,rhoA0)
      class(relax_sg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,Yv_,rho0,rhoA0
      relax_sg_ig_get_T_lvg=(1.0_WP-Yv_)/                                                                                   &
      &         ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)+                          &
      &           rhoA0*((this%gas%get_species_gamma(this%indV)-1.0_WP)*this%gas%get_species_cv(this%indV)*Yv_+             &
      &                  (this%gas%get_species_gamma(this%indA)-1.0_WP)*this%gas%get_species_cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function relax_sg_ig_get_T_lvg

   !> Update the quadratic coefficients of equilibrium temperature equation
   subroutine relax_sg_ig_get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      class(relax_sg_ig), intent(in)  :: this
      real(WP), intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      ! Coefficients (rho0 = sum of phase masses; rhoe0 = sum of phase energies — both invariant during steps 1-2)
      ap=rho0*this%liq%cv*cvG*((GammaG-1.0_WP)*(p_eq+this%liq%gamma*this%liq%pinf)-(this%liq%gamma-1.0_WP)*p_eq)
      bp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv*p_eq-(GammaG-1.0_WP)*cvG*(p_eq+this%liq%pinf))              +&
      &  rho0*((GammaG-1.0_WP)*cvG*this%liq%q*(p_eq+this%liq%pinf)-(this%liq%gamma-1.0_WP)*this%liq%cv*qG*p_eq) +&
      &  cvG*p_eq*(p_eq+this%liq%pinf)-this%liq%cv*p_eq*(p_eq+this%liq%gamma*this%liq%pinf)
      dp=(qG-this%liq%q)*(p_eq+this%liq%pinf)*p_eq
      ! Pressure derivative of the coefficients
      dapdp=rho0*this%liq%cv*cvG*(GammaG-this%liq%gamma)
      dbpdp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(GammaG-1.0_WP)*cvG)                                     +&
      &     rho0*((GammaG-1.0_WP)*cvG*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qG)                        +&
      &     cvG*(2.0_WP*p_eq+this%liq%pinf)-this%liq%cv*(2.0_WP*p_eq+this%liq%gamma*this%liq%pinf)
      ddpdp=(qG-this%liq%q)*(2.0_WP*p_eq+this%liq%pinf)
   end subroutine relax_sg_ig_get_coeffs_lv


end module relax_sg_ig_class