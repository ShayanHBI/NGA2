!> NASG liquid and ideal gas relaxation model
module relax_nasg_ig_class
   use precision,         only: WP
   use relax_sg_ig_class, only: relax_sg_ig,Mv,Ma
   use sg_class,          only: sg
   use nasg_class,        only: nasg
   use igmix_class,       only: igmix
   implicit none
   private

   public :: relax_nasg_ig

   type, extends(relax_sg_ig) :: relax_nasg_ig
      !> Typed pointer for direct b access (avoids select type in relax_p, relax_pT, helpers)
      type(nasg), pointer :: liq_nasg=>null()
      !> Relaxation parameters
      real(WP) :: mu,nu
   contains
      procedure :: initialize   =>relax_nasg_ig_initialize
      procedure :: relax_p      =>relax_nasg_ig_relax_p
      procedure :: relax_pT     =>relax_nasg_ig_relax_pT
      procedure :: get_T_lvg    =>relax_nasg_ig_get_T_lvg
      procedure :: get_coeffs_lv=>relax_nasg_ig_get_coeffs_lv
      procedure :: get_p_eq     =>relax_nasg_ig_get_p_eq
   end type relax_nasg_ig

contains

   !> Initialize: call parent, then add NASG-specific fields.
   !> Takes class(sg) to match parent interface; select type extracts the nasg-specific fields.
   subroutine relax_nasg_ig_initialize(this,liq,gas,indV,indA,p_cav,VF_nuc)
      class(relax_nasg_ig), intent(inout) :: this
      class(sg),    target, intent(in)    :: liq
      class(igmix), target, intent(in)    :: gas
      integer, intent(in) :: indV,indA
      real(WP), intent(in), optional :: p_cav,VF_nuc
      real(WP) :: cpV,cvV,RV
      ! Parent sets this%liq, this%gas, this%indV/A, this%AS-DS, this%ES=0, this%p_cav/VF_nuc
      call this%relax_sg_ig%initialize(liq=liq,gas=gas,indV=indV,indA=indA,p_cav=p_cav,VF_nuc=VF_nuc)
      ! Set typed pointer and override ES (select type needed to get type(nasg) from class(sg))
      select type (liq)
      type is (nasg)
         this%liq_nasg=>liq
         cpV=gas%get_species_cp(indV)
         cvV=gas%get_species_cv(indV)
         RV=cpV-cvV
         this%ES=liq%b/RV
      end select
   end subroutine relax_nasg_ig_initialize

   !> Mechanical relaxation
   !> Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
   subroutine relax_nasg_ig_relax_p(this,dt,VF,Q,Pjump)
      class(relax_nasg_ig),    intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,y
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,PL,PG,IL,IG,ZL,ZG,Pint
      real(WP) :: xiL,xiG,xiLinv,xiGinv
      real(WP) :: VFeq,VF0,Peq
      real(WP) :: Kp
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      if ((VF.eq.0.0_WP).or.(VF.eq.1.0_WP)) return
      ! Store the input state
      allocate(Q0(size(Q)))
      VF0=VF
      Q0=Q
      ! Set gas EOS coefficients from vapor mass fraction
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Get phasic thermodynamic quantities
      RHOL=Q(1)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF)
      IL=Q(3)/Q(1)
      IG=Q(4)/Q(2)
      PL=this%liq%get_p_from_rho_e(rho=RHOL,e=IL)
      PG=this%gas%get_p_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=PL,rho=RHOL)
      CG=this%gas%get_c_from_p_rho(p=PG,rho=RHOG,y=y)
      ! Handle limit cases - should mass/energy be transfered or lost? - this should probably never happen...
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         Q(8)=Yv*Q(2)
         call dealloc()
         return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(8)=0.0_WP
         call dealloc()
         return
      end if
      ! Get phasic impedances
      ZL=Q(1)/(       VF)*CL
      ZG=Q(2)/(1.0_WP-VF)*CG
      ! Calculate model interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Setup the ODE coefficients
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL)
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      xiL=         VF/(GL*(Pint-PL)+RHOL*CL**2)
      xiG=(1.0_WP-VF)/(GG*(Pint-PG)+RHOG*CG**2)
      xiLinv=1.0_WP/xiL
      xiGinv=1.0_WP/xiG
      Kp=this%mu*(xiLinv+xiGinv)
      ! Get equilibrium volume fraction
      VFeq=VF-(PG-PL)/(xiLinv+xiGinv)*(1.0_WP-exp(-Kp*dt))
      if ((VFeq.lt.0.0_WP).or.(VFeq.gt.1.0_WP)) then
         call restore()
         call dealloc()
         return
      end if
      ! Get equilibrium pressure
      Peq=this%get_p_eq(VFeq,Q0,qG,gammaG)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         call restore()
         call dealloc()
         return
      end if
      ! Adjust densities
      RHOL=Q0(1)/(       VFeq)
      RHOG=Q0(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL)
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      call dealloc()
   contains
      !> Reset the output to the initial values
      subroutine restore()
         VF=VF0
         Q=Q0
      end subroutine restore
      !> Release memory allocated for Q0 and y
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
      end subroutine dealloc
   end subroutine relax_nasg_ig_relax_p

   !> Mechanical and thermal relaxation
   !> Pelanti 2022: https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
   subroutine relax_nasg_ig_relax_pT(this,dt,VF,Q,Pjump)
      class(relax_nasg_ig),    intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      real(WP), dimension(:),  allocatable   :: Q0,y
      real(WP) :: RHOL,RHOG,CL,CG,GL,GG,TL,TG,IL,IG
      real(WP) :: PHIL,PHIG,zetaL,zetaG,Z,D,COF
      real(WP) :: xiTL,xiTG,xiTLinv,xiTGinv
      real(WP) :: VFeq,Peq!,Teq
      real(WP) :: KT
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      if ((VF.eq.0.0_WP).or.(VF.eq.1.0_WP)) return
      ! Store initial Q before relax_p so get_p_eq sees the same Q0 as the original monolithic code
      allocate(Q0(size(Q)))
      Q0=Q
      ! ================ First step for mechanical relaxation ================
      call this%relax_p(dt,VF,Q,Pjump)
      ! ================= Second step for thermal relaxation =================
      ! Set gas EOS coefficients (Q(8) and Q(2) unchanged by relax_p)
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      call this%gas%get_mix_coeffs(y=y,cv=cvG,cp=cpG,q=qG,gamma=gammaG)
      ! Get densities
      RHOL=Q(1)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF)
      ! Recover p from the dominant phase for better numerical stability
      if (VF.gt.0.5_WP) then
         Peq=this%liq%get_p_from_rho_e(rho=RHOL,e=Q(3)/Q(1))
      else
         Peq=this%gas%get_p_from_rho_e(rho=RHOG,e=Q(4)/Q(2),y=y)
      end if
      ! Return if unphysical pressure
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(Q0,y)
         return
      end if
      ! Update phasic thermodynamic quantities
      IL=Q(3)/Q(1)
      IG=Q(4)/Q(2)
      TL=this%liq%get_T_from_p_rho(p=Peq,rho=RHOL)
      TG=this%gas%get_T_from_p_rho(p=Peq,rho=RHOG,y=y)
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL)
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=Peq,rho=RHOL)
      CG=this%gas%get_c_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! Setup ODE coefficients
      Z=(1.0_WP-VF)*GL+VF*GG
      D=VF*RHOG*CG**2+(1.0_WP-VF)*RHOL*CL**2
      PHIL=this%liq%get_drhodT_const_p_from_rho_T(rho=RHOL,T=TL)
      PHIG=this%gas%get_drhodT_const_p_from_rho_T(rho=RHOG,T=TG,y=y)
      zetaL=this%liq%get_drhodp_const_T_from_rho_T(rho=RHOL,T=TL)
      zetaG=this%gas%get_drhodp_const_T_from_rho_T(rho=RHOG,T=TG,y=y)
      COF=GL*RHOG*CG**2-GG*RHOL*CL**2
      xiTL=-PHIL*D/(RHOL/(       VF)*Z+zetaL*COF)
      xiTG=-PHIG*D/(RHOG/(1.0_WP-VF)*Z-zetaG*COF)
      xiTLinv=1.0_WP/xiTL
      xiTGinv=1.0_WP/xiTG
      KT=this%nu*(xiTLinv+xiTGinv)
      ! Get equilibrium VF, T ,and p
      VFeq=VF+Z/D*(TG-TL)/(xiTLinv+xiTGinv)*(1.0_WP-exp(-KT*dt))
      ! Teq=(xiTL*TL+xiTG*TG)/(xiTL+xiTG)
      Peq=this%get_p_eq(VFeq,Q0,qG,gammaG)
      ! Check if pressure is sound
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y)
         deallocate(Q0)
         return
      end if
      ! Clean up solution
      if (VFeq.lt.0.0_WP) then
         VFeq=0.0_WP
         VF=VFeq
         Peq=max(Peq,-this%liq%pinf)
         Q(2)=sum(Q(1:2)); Q(8)=Q(8)+Q(1); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         deallocate(Q0,y)
         return
      end if
      if (VFeq.gt.1.0_WP) then
         VFeq=1.0_WP
         VF=VFeq
         Peq=max(Peq,0.0_WP)
         Q(1)=sum(Q(1:2)); Q(8)=0.0_WP; Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         deallocate(Q0,y)
         return
      end if
      ! Adjust densities
      RHOL=Q(1)/(       VFeq)
      RHOG=Q(2)/(1.0_WP-VFeq)
      ! Adjust conserved quantities (Masses and velocities remain unchanged)
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=RHOL)
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=RHOG,y=y)
      ! Release memory
      deallocate(Q0,y)
   end subroutine relax_nasg_ig_relax_pT

   !> Energy-conserving equilibrium pressure at given VF_ (Used for p and pT relaxation steps only)
   real(WP) function relax_nasg_ig_get_p_eq(this,VF_,Q0_,qG_,gammaG_)
      class(relax_nasg_ig), intent(in) :: this
      real(WP), intent(in) :: VF_
      real(WP), dimension(:), intent(in) :: Q0_
      real(WP), intent(in) :: qG_,gammaG_
      real(WP) :: one_brho
      one_brho=1.0_WP-this%liq_nasg%b*Q0_(1)/VF_
      relax_nasg_ig_get_p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-one_brho*VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_)/&
      &                      (one_brho*VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function relax_nasg_ig_get_p_eq

   !> Equilibrium T from energy conservation
   real(WP) function relax_nasg_ig_get_T_lvg(this,p_,Yv_,rho0,rhoA0)
      class(relax_nasg_ig), intent(in) :: this
      real(WP), intent(in) :: p_,Yv_,rho0,rhoA0
      relax_nasg_ig_get_T_lvg=(1.0_WP-Yv_-this%liq_nasg%b*(rho0*(1.0_WP-Yv_)-rhoA0))/&
      & ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)                           +&
          rhoA0*((this%gas%get_species_gamma(this%indV)-1.0_WP)*this%gas%get_species_cv(this%indV)*Yv_              +&
      &          (this%gas%get_species_gamma(this%indA)-1.0_WP)*this%gas%get_species_cv(this%indA)*(1.0_WP-Yv_))/p_)
   end function relax_nasg_ig_get_T_lvg

   !> Update the quadratic coefficients of equilibrium temperature equation
   subroutine relax_nasg_ig_get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,ap,bp,dp,dapdp,dbpdp,ddpdp)
      class(relax_nasg_ig), intent(in)  :: this
      real(WP), intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG
      real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      real(WP) :: cvV_,gammaV_,qV_
      cvV_   =this%gas%get_species_cv(this%indV)
      gammaV_=this%gas%get_species_gamma(this%indV)
      qV_    =this%gas%get_species_q(this%indV)
      ! Coefficients
      ap=rho0*this%liq%cv*cvV_*((gammaV_-this%liq%gamma)*p_eq+this%liq%gamma*(gammaV_-1.0_WP)*this%liq%pinf)
      bp=(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq**2                                                    +&
      &  (this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)                                                          -&
      &   this%liq%gamma*this%liq%cv)+rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
      &   rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_eq                                    +&
      &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(rho0*this%liq%q-rhoe0)
      dp=p_eq*(p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
      ! Pressure derivative of the coefficients
      dapdp=rho0*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
      dbpdp=2.0_WP*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq                                             +&
      &     this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%gamma*this%liq%cv)                            +&
      &     rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)                          +&
      &     rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)
      ddpdp=(2.0_WP*p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
   end subroutine relax_nasg_ig_get_coeffs_lv

end module relax_nasg_ig_class
