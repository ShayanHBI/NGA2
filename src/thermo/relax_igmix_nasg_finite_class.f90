!> NASG-liquid + ideal-gas-mixture relaxation (extends SG-IG mixture with co-volume b corrections)
!> Inherits pTg_relax unchanged from parent; virtual dispatch into the overridden helpers
!> (relax_p, pT_relax, get_T_lvg, get_coeffs_lv, get_p_eq) does the right thing.
!> References Pelanti 2022 https://doi.org/10.1016/j.ijmultiphaseflow.2022.104097
module relax_igmix_nasg_finite_class
   use precision,               only: WP
   use relax_igmix_sg_class,    only: relax_igmix_sg,sg_apply=>apply,Prelax,PTgrelax,Mv,Ma,dbg_cell,dbg_i,dbg_j
   use stiffened_gas_class,     only: stiffened_gas
   use nasg_class,              only: nasg
   use igmix_class,             only: igmix
   use thermorelax_class,       only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE,RELAX_NUC_FAILED
   implicit none
   private

   public :: relax_igmix_nasg_finite

   type, extends(relax_igmix_sg) :: relax_igmix_nasg_finite
      !> Typed pointer for direct b access (avoids select type in hot path)
      type(nasg), pointer :: liq_nasg=>null()
      !> matm-style p/pT relaxation parameters (mirror relax_ig_nasg)
      real(WP) :: RHOGmin=1.0e-2_WP   !< Skip mechanical relax when gas density falls below this
      real(WP) :: phist  =1.0_WP      !< Temporal weighting on equilibrium pressure
      real(WP) :: phi0   =0.0_WP      !< Temporal weighting on interface pressure
      !> Relaxation parameters
      real(WP) :: lambda=huge(1.0_WP) !< Mechanical rate [1/(Pa s)] (Pelanti mu); huge = instantaneous
      real(WP) :: theta =huge(1.0_WP) !< Thermal rate [W/(m^3 K)] (Pelanti vartheta); huge = instantaneous
      real(WP) :: nu    =huge(1.0_WP) !< Chemical rate [kg s/m^5] (Pelanti nu); huge = instantaneous
   contains
      procedure :: initialize
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: get_chem_transfer
      procedure :: get_T_lvg
      procedure :: get_coeffs_lv
      procedure :: get_p_eq
   end type relax_igmix_nasg_finite

contains

   !> Initialize: call parent (sets liq, gas, indV/A, AS-DS, ES=0), then set NASG-specific fields
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_nasg_finite),       intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      real(WP) :: cpV,cvV,RV
      ! Parent: sets this%liq=>liq, this%gas=>gas, indV/A, saturation coeffs AS-DS, ES=0
      call this%relax_igmix_sg%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      ! Set typed pointer and override ES (select type to extract type(nasg) from class(stiffened_gas))
      select type (liq)
      type is (nasg)
         this%liq_nasg=>liq
         cpV=gas%cp(indV); cvV=gas%cv(indV); RV=cpV-cvV
         this%ES=liq%b/RV
      end select
   end subroutine initialize

   !> Finite-rate Prelax leaves the phases out of pressure equilibrium; the parent's fixed-VF
   !> energy split (stage 2) would undo that, so this path bypasses it. Everything else is the parent's apply.
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      use messager,     only: die
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),                       intent(in)    :: dt
      real(WP),                       intent(inout) :: VF
      real(WP), dimension(:),         intent(inout) :: Q
      real(WP),                       intent(in)    :: Pjump
      integer,  optional,             intent(out)   :: ierr
      integer :: ier
      real(WP) :: VF0                   ! debug
      real(WP), dimension(size(Q)) :: Q0 ! debug
      ! Each step is posed on the previous step's equilibrium manifold: pT needs instantaneous p, pTg also instantaneous T
      if (this%model.ne.Prelax.and.this%lambda.lt.huge(1.0_WP)) call die('[relax_igmix_nasg_finite apply] pT, pTg and pThybrid require lambda=huge')
      if (this%model.eq.PTgrelax.and.this%theta.lt.huge(1.0_WP)) call die('[relax_igmix_nasg_finite apply] pTg requires theta=huge')
      if (this%model.ne.Prelax.or.this%lambda.ge.huge(1.0_WP)) then
         ! Parent's apply on this itself, so p_relax/get_p_eq still dispatch to this class
         VF0=VF; Q0=Q ! debug
         call sg_apply(this,dt,VF,Q,Pjump,ierr)
         ! debug: report the entry state of a call that leaves NaN or multiplies the internal energy by >10
         if (any(Q.ne.Q).or.VF.ne.VF.or.Q(3)+Q(4).gt.10.0_WP*max(Q0(3)+Q0(4),tiny(1.0_WP))) then ! debug
            print'(a,es25.17,a,es25.17)','BADEXIT entry VF=',VF0,' Pjump=',Pjump ! debug
            print'(a,*(es25.17))','BADEXIT entry Q=',Q0 ! debug
            print'(a,es25.17,a,*(es25.17))','BADEXIT exit  VF=',VF,' Q=',Q ! debug
         end if ! debug
         return
      end if
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         if (present(ierr)) ierr=RELAX_DEGENERATE
         return
      end if
      if (dbg_cell) call this%debug_dump('ENTRY apply (finite Prelax)',VF,Q,Pjump)
      call this%p_relax(dt,VF,Q,Pjump,ier)
      ! Ledger as the parent: 1 proposal succeeded, 5 non-positive phase mass, 2 proposal failed (cell left as is)
      if (ier.eq.RELAX_OK) then
         this%acc(1)=this%acc(1)+1.0_WP
      else if (Q(1).le.0.0_WP.or.Q(2).le.0.0_WP) then
         this%acc(5)=this%acc(5)+1.0_WP
      else
         this%acc(2)=this%acc(2)+1.0_WP
      end if
      if (dbg_cell) call this%debug_dump('EXIT apply (finite Prelax)',VF,Q,Pjump,ier)
      if (present(ierr)) ierr=ier
   end subroutine apply

   !> Mechanical relaxation
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      integer,  optional,      intent(out)   :: ierr
      real(WP) :: PG,PL,ZG,ZL,Pint,bL,rhoL,rhoG,cL,cG,GamL,GamG
      real(WP) :: xiLinv,xiGinv,Kdt,frac,Peq,VFeq,dW
      logical  :: inst
      real(WP) :: cvG,cpG,qG,gammaG,gL,pinfL,qL,Yv
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Instantaneous (lambda=huge): closed-form mechanical equilibrium instead of the exponential step
      if (this%lambda.ge.huge(1.0_WP)) then
         call p_exact(this,VF,Q,Pjump,ierr)
         return
      end if
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) then; if (present(ierr)) ierr=RELAX_VACUUM_GAS; return; end if
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      ! Yv=Q(iVQ)/Q(2)
      ! Clamp Yv to [0,1] as get_primitive does: near-total VF collapse can leave Q(iVQ)>Q(2) ! debug
      Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Liquid parameters (b via the typed nasg pointer; EOS calls dispatch to nasg)
      bL=this%liq_nasg%b; gL=this%liq%gamma; pinfL=this%liq%pinf; qL=this%liq%q
      rhoL=Q(1)/(VF); rhoG=Q(2)/(1.0_WP-VF)
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=y)
      ! Handle limit cases
      if (PL.le.-pinfL.or.1.0_WP-bL*rhoL.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (PG.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      cL=this%liq%get_c_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
      cG=this%gas%get_c_from_p_rho(p=PG,rho=rhoG,y=y)
      ! Phasic acoustic impedances (rho*c)
      ZL=rhoL*cL
      ZG=rhoG*cG
      ! Calculate model interface pressure
      ! Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      Pint=(ZG*PL+ZL*(PG+Pjump))/(ZG+ZL)
      ! Setup the ODE coefficients: inverse capacities 1/xi_k (Pelanti 2022 eq. 51)
      GamL=this%liq%get_gruneisen_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      GamG=this%gas%get_gruneisen_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=y)
      xiLinv=(GamL*(Pint-PL)+rhoL*cL**2)/(       VF)
      xiGinv=(GamG*(Pint-Pjump-PG)+rhoG*cG**2)/(1.0_WP-VF)
      if (xiLinv+xiGinv.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Fraction of the gap crossed in dt, 1-exp(-Kp*dt) with Kp=lambda*(1/xi_l+1/xi_g) (eq. 55).
      ! Instantaneous when lambda=huge or lambda*Kdt>=50 (exp(-50)<eps); compared in logs so lambda*Kdt never overflows
      Kdt=(xiLinv+xiGinv)*dt
      inst=.true.; frac=1.0_WP
      if (this%lambda.lt.huge(1.0_WP)) then
         if (this%lambda.le.0.0_WP.or.Kdt.le.0.0_WP) then
            inst=.false.; frac=0.0_WP
         else if (log(this%lambda)+log(Kdt).lt.log(50.0_WP)) then
            inst=.false.; frac=1.0_WP-exp(-this%lambda*Kdt)
         end if
      end if
      ! Instantaneous by rate (lambda*Kdt>=50): closed form as well
      if (inst) then
         call p_exact(this,VF,Q,Pjump,ierr)
         return
      end if
      ! Get equilibrium volume fraction
      VFeq=VF-(PG+Pjump-PL)/(xiLinv+xiGinv)*frac
      if (VFeq.le.0.0_WP.or.VFeq.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Clamp per-call phase-volume change to a factor VFratmax (bounds consistent for VFratmax>=1)
      VFeq=max(1.0_WP-this%VFratmax*(1.0_WP-VF),VF/this%VFratmax,min(VFeq,1.0_WP-(1.0_WP-VF)/this%VFratmax,this%VFratmax*VF))
      if (inst) then
         ! Get equilibrium pressure at the final VF: liquid at Peq, gas at Peq-Pjump (eq. 59)
         ! Peq=this%get_p_eq(VFeq,Q,qG,gammaG,Pjump)
         Peq=this%get_p_eq(VFeq,[Q(1:3),Q(4)-Pjump*(VFeq-VF)],qG,gammaG,Pjump)
         ! Soundness check (pinf_g=0)
         if (Peq.le.-pinfL) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
         if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
         ! Adjust densities
         rhoL=Q(1)/(       VFeq)
         rhoG=Q(2)/(1.0_WP-VFeq)
         ! Adjust VF and conserved quantities (Masses and velocities remain unchanged)
         VF=VFeq
         Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,      rho=rhoL,y=[1.0_WP])
         Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq-Pjump,rho=rhoG,y=y)
      else
         ! Finite rate: phasic energies change only by the p_I*dV work at the frozen interface pressure
         ! (eqs. 1e-1f); the phases are left out of pressure equilibrium
         dW=Pint*(VFeq-VF)
         ! Gas works at Pint-Pjump; Pjump*dV goes to the surface energy
         if (Q(3)-dW.le.0.0_WP.or.Q(4)+dW-Pjump*(VFeq-VF).le.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
         Q(3)=Q(3)-dW
         Q(4)=Q(4)+dW-Pjump*(VFeq-VF)
         VF=VFeq
      end if
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_relax

   !> Mechanical + thermal relaxation
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),                intent(in)    :: dt
      real(WP),                intent(inout) :: VF
      real(WP), dimension(:),  intent(inout) :: Q
      real(WP),                intent(in)    :: Pjump
      integer,  optional,      intent(out)   :: ierr
      real(WP) :: Peq,VFeq,bL,cvG,cpG,qG,gammaG,Yv
      real(WP) :: RHOL,RHOG,PL,PG,IL,IG,TL,TG,GL,GG,CL,CG
      real(WP) :: Z,D,COF,PHIL,PHIG,zetaL,zetaG,xiTLinv,xiTGinv,Kdt,frac
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Mechanical relaxation first (its own ierr is not propagated; pT owns the final verdict)
      call this%p_relax(dt,VF,Q,Pjump)
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Instantaneous (theta=huge): closed-form p,T equilibrium instead of the exponential step
      if (this%theta.ge.huge(1.0_WP)) then
         call pT_exact(this,VF,Q,Pjump,ierr)
         return
      end if
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      ! Clamp Yv to [0,1]
      Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP))
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Liquid shorthands (b via typed nasg pointer)
      bL=this%liq_nasg%b
      RHOL=Q(1)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF)
      if (1.0_WP-bL*RHOL.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if   ! liquid past co-volume packing limit
      ! Post-mechanical pressure (eq. 59): liquid at PL, gas at PG=PL-Pjump
      PL=this%get_p_eq(VF,Q,qG,gammaG,Pjump)
      PG=PL-Pjump
      ! Update phasic thermodynamic quantities (each phase at its own pressure)
      IL=Q(3)/Q(1)
      IG=Q(4)/Q(2)
      TL=this%liq%get_T_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
      TG=this%gas%get_T_from_p_rho(p=PG,rho=RHOG,y=y)
      if (TL.le.0.0_WP.or.TG.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      GL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=IL,y=[1.0_WP])
      GG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=IG,y=y)
      CL=this%liq%get_c_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
      CG=this%gas%get_c_from_p_rho(p=PG,rho=RHOG,y=y)
      ! Setup ODE coefficients
      Z=(1.0_WP-VF)*GL+VF*GG
      D=VF*RHOG*CG**2+(1.0_WP-VF)*RHOL*CL**2
      ! phi=(drho/dT)_p and zeta=(drho/dp)_T (eqs. 16, 34); liquid branches match nasg get_rho_from_p_T
      if (bL*RHOL.gt.this%liq_nasg%brhomax) then
         PHIL =-RHOL/TL
         zetaL= RHOL/(PL+this%liq%pinf)
      else
         PHIL =-RHOL*(1.0_WP-bL*RHOL)/TL
         zetaL= RHOL*(1.0_WP-bL*RHOL)/(PL+this%liq%pinf)
      end if
      PHIG =-RHOG/TG
      zetaG= RHOG/PG
      COF=GL*RHOG*CG**2-GG*RHOL*CL**2
      ! Inverse thermal capacities 1/xiT_k (eq. 62)
      xiTLinv=-(RHOL/(       VF)*Z+zetaL*COF)/(PHIL*D)
      xiTGinv=-(RHOG/(1.0_WP-VF)*Z-zetaG*COF)/(PHIG*D)
      if (xiTLinv+xiTGinv.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Fraction of the gap crossed in dt, 1-exp(-KT*dt) with KT=theta*(1/xiT_l+1/xiT_g) (eq. 66).
      ! Instantaneous when theta=huge or theta*Kdt>=50 (exp(-50)<eps); compared in logs so theta*Kdt never overflows
      Kdt=(xiTLinv+xiTGinv)*dt
      frac=1.0_WP
      if (this%theta.lt.huge(1.0_WP)) then
         if (this%theta.le.0.0_WP.or.Kdt.le.0.0_WP) then
            frac=0.0_WP
         else if (log(this%theta)+log(Kdt).lt.log(50.0_WP)) then
            frac=1.0_WP-exp(-this%theta*Kdt)
         end if
      end if
      ! Instantaneous by rate (theta*Kdt>=50): closed form as well
      if (frac.ge.1.0_WP) then
         call pT_exact(this,VF,Q,Pjump,ierr)
         return
      end if
      ! Get equilibrium VF (eq. 67) and p (eq. 70)
      VFeq=VF+Z/D*(TG-TL)/(xiTLinv+xiTGinv)*frac
      if (VFeq.le.0.0_WP.or.VFeq.ge.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Peq=this%get_p_eq(VFeq,Q,qG,gammaG,Pjump)
      Peq=this%get_p_eq(VFeq,[Q(1:3),Q(4)-Pjump*(VFeq-VF)],qG,gammaG,Pjump)
      ! Soundness check (pinf_g=0): liquid p=Peq>-pinf, gas p=Peq-Pjump>0
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP)   then; if (present(ierr)) ierr=RELAX_BAD_GAS;    return; end if
      ! Adjust densities
      RHOL=Q(1)/(       VFeq)
      RHOG=Q(2)/(1.0_WP-VFeq)
      ! Adjust VF and conserved quantities
      VF=VFeq
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,      rho=RHOL,y=[1.0_WP])
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq-Pjump,rho=RHOG,y=y)
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_relax

   !> Instantaneous mechanical relaxation: matm-style quadratic of relax_igmix_nasg%p_relax (PL=Peq, PG=Peq-Pjump),
   !> with surface energy conserved: the liquid works at Pw, the gas receives Pw-Pjump, so Q(3)+Q(4)+Pjump*VF is kept
   subroutine p_exact(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),                       intent(inout) :: VF
      real(WP), dimension(:),         intent(inout) :: Q
      real(WP),                       intent(in)    :: Pjump
      integer,  optional,             intent(out)   :: ierr
      real(WP) :: PG,PL,ZG,ZL,Pint,cJ,bL,rhoL,rhoG
      real(WP) :: a,b,d,n1,n0,d1,d0,Peq,VFeq,cw,cg
      real(WP) :: cvG,cpG,qG,gammaG,gL,pinfL,qL,Yv
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Skip near-pure-liquid cells (gas density too low)
      if (Q(2)/(1.0_WP-VF).lt.this%RHOGmin) then; if (present(ierr)) ierr=RELAX_VACUUM_GAS; return; end if
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP))
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Liquid parameters (b via the typed nasg pointer; EOS calls dispatch to nasg)
      bL=this%liq_nasg%b; gL=this%liq%gamma; pinfL=this%liq%pinf; qL=this%liq%q
      rhoL=Q(1)/(VF); rhoG=Q(2)/(1.0_WP-VF)
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=rhoL,e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=rhoG,e=Q(4)/Q(2),y=y)
      ! No cavitation model: leave sub-vacuum / over-packed cells untouched (avoid NaN impedances)
      if (PL.le.-pinfL.or.1.0_WP-bL*rhoL.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (PG.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Phasic acoustic impedances (rho*c)
      ZL=rhoL*this%liq%get_c_from_p_rho(p=PL,rho=rhoL,y=[1.0_WP])
      ZG=rhoG*this%gas%get_c_from_p_rho(p=PG,rho=rhoG,y=y)
      cJ=ZL/(ZG+ZL)
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Work pressure Pw=phist*Peq+cw for the liquid; the gas receives Pw-Pjump, i.e. cg=cw-Pjump in the gas equation
      cw=this%phi0*Pint-this%phist*cJ*Pjump
      cg=cw-Pjump
      ! Quadratic (pinf_g=0; liquid co-volume b in n1,n0; mixture gas gamma/q)
      n1=VF*this%phist+Q(1)*bL/(gL-1.0_WP)
      n0=VF*cw+Q(3)-Q(1)*qL+Q(1)*bL*gL/(gL-1.0_WP)*pinfL
      d1=this%phist+1.0_WP/(gL-1.0_WP)
      d0=cw+gL/(gL-1.0_WP)*pinfL
      a=d1*(1.0_WP/(gammaG-1.0_WP)+this%phist*VF)+n1*(-1.0_WP/(gammaG-1.0_WP)-this%phist)
      b=d1*(-Pjump/(gammaG-1.0_WP)+Q(2)*qG-Q(4)+VF*cg)+n1*(Pjump/(gammaG-1.0_WP)-cg)+d0*(1.0_WP/(gammaG-1.0_WP)+this%phist*VF)+n0*(-1.0_WP/(gammaG-1.0_WP)-this%phist)
      d=d0*(-Pjump/(gammaG-1.0_WP)+Q(2)*qG-Q(4)+VF*cg)+n0*(Pjump/(gammaG-1.0_WP)-cg)
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g=0)
      if (Peq.le.-pinfL) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Equilibrium VF (strict bounds)
      VFeq=(n1*Peq+n0)/(d1*Peq+d0)
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Clamp per-call phase-volume change to a factor VFratmax (bounds consistent for VFratmax>=1)
      VFeq=max(1.0_WP-this%VFratmax*(1.0_WP-VF),VF/this%VFratmax, &
      &        min(VFeq,1.0_WP-(1.0_WP-VF)/this%VFratmax,this%VFratmax*VF))
      ! Update Q with p*dV work at the relaxed interface pressure; Pjump*dV goes to the surface energy
      Q(3)=Q(3)-(this%phist*Peq+cw)*(VFeq-VF)
      Q(4)=Q(4)+(this%phist*Peq+cg)*(VFeq-VF)
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_exact

   !> Instantaneous thermal relaxation on a pressure-relaxed state: matm-style quadratic of relax_igmix_nasg%pT_relax
   !> (PL=Peq, PG=Peq-Pjump, TL=TG), with surface energy conserved: the energy target is E'=Eth+Pjump*(VF-b*Q(1))
   subroutine pT_exact(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),                       intent(inout) :: VF
      real(WP), dimension(:),         intent(inout) :: Q
      real(WP),                       intent(in)    :: Pjump
      integer,  optional,             intent(out)   :: ierr
      real(WP) :: a,b,d,Peq,VFeq,Eth,C,S
      real(WP) :: cv1,cv2,g1,g2,pinf,R1,R2,bL,ombm,cpG,qG,Yv
      real(WP) :: y(this%gas%ns)
      integer  :: iVQ
      ! Skip if any conserved quantity is non-positive
      if (any(Q(1:4).le.0.0_WP)) then; if (present(ierr)) ierr=RELAX_DEGENERATE; return; end if
      ! Frozen gas composition -> mixture-effective ideal-gas parameters
      iVQ=7+this%liq%ns+this%indV-1
      Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP))
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cv2=sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG=sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      g2 =cpG/cv2
      ! Liquid shorthands (b via typed nasg pointer)
      cv1=this%liq%cv; g1=this%liq%gamma; pinf=this%liq%pinf
      R1=cv1*(g1-1.0_WP); R2=cv2*(g2-1.0_WP)
      bL=this%liq_nasg%b
      ombm=1.0_WP-bL*Q(1)
      if (ombm.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if   ! liquid past co-volume packing limit
      ! Energy target (formation energies removed), including the surface energy released as VF changes
      Eth=Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*qG+Pjump*(VF-bL*Q(1))
      ! Quadratic for liquid equilibrium pressure Peq (gas pressure Peq-Pjump), TL=TG, pinf_g=0
      C=Q(1)*cv1+Q(2)*cv2
      S=Q(1)*cv1*g1*pinf+Q(2)*cv2*pinf+Q(1)*R1*Pjump
      a=ombm*C
      b=ombm*(S-Pjump*C)-Eth*(Q(1)*R1+Q(2)*R2)
      d=-ombm*Pjump*S+Pjump*Eth*Q(1)*R1-Eth*Q(2)*R2*pinf
      ! Equilibrium pressure
      if (b**2-4.0_WP*a*d.lt.0.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Soundness check (pinf_g=0): liquid p=Peq>-pinf, gas p=Peq-Pjump>0
      if (Peq.le.-this%liq%pinf) then; if (present(ierr)) ierr=RELAX_BAD_LIQUID; return; end if
      if (Peq-Pjump.le.0.0_WP) then; if (present(ierr)) ierr=RELAX_BAD_GAS; return; end if
      ! Equilibrium VF (strict bounds); co-volume floor b*Q(1) appears naturally
      VFeq=bL*Q(1)+ombm*Q(1)*R1*(Peq-Pjump)/(Q(1)*R1*(Peq-Pjump)+Q(2)*R2*(Peq+pinf))
      if (VFeq.lt.0.0_WP.or.VFeq.gt.1.0_WP) then; if (present(ierr)) ierr=RELAX_FAILED; return; end if
      ! Update Q with the new equilibrium state (co-volume in liquid rhoe, formation energies re-added)
      Q(3)=(VFeq-bL*Q(1))*(Peq+g1*pinf)/(g1-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VFeq)*(Peq-Pjump  )/(g2-1.0_WP)+Q(2)*qG
      VF=VFeq
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_exact

   !> Mechanical + thermal + chemical (phase change) relaxation.
   !> Nucleates a tiny opposite phase in metastable pure cells; calls pT_relax;
   !> runs pure-phase admissibility tests; falls back to LV (pure water) or LVG (with non-condensable) Newton solves.
   !> Finite-rate chemical step when nu<huge (Pelanti 2022 eqs. 71-83).
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      implicit none
      class(relax_igmix_nasg_finite), intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(:), allocatable   :: Qin,Q0,y
      real(WP) :: VFin,VF0,p,T,Yv,Yvin,Yv0
      real(WP) :: rho0
      real(WP) :: rhoL,rhoG
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP), parameter :: p_eps=1.0e-10_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Y_small=0.001_WP
      real(WP), parameter :: VF_nuc=0.01_WP
      logical :: chem_relax,near_pure,nucleated,cavitated,condensed
      integer :: ier
      if (dbg_cell) then
         print*,'--------------------------------------------------'
         print*,'inside pTg_relax, VF=',VF
      end if
      allocate(y(this%gas%ns))
      ! Default to success
      ier=RELAX_OK
      ! Check if pure phase is already stable
      pure_phase_stability: block
         real(WP) :: TL,pL,eL
         real(WP) :: TG,pG,eG
         real(WP) :: pV,Tsat,pVsat
         real(WP), parameter :: pure_VFlo=0.05_WP,pure_VFhi=0.95_WP
         logical  :: conv
         integer  :: it
         ! Instantaneous-equilibrium shortcut: skipped at finite nu
         if (this%nu.lt.huge(1.0_WP)) exit pure_phase_stability
         if (VF.lt.pure_VFlo) then
            if (dbg_cell) print*,'attempting pure gas'
            ! Transfer all liquid mass and energy to vapor
            rhoG=Q(1)+Q(2)
            eG=(Q(3)+Q(4)+Pjump*VF)/rhoG
            Yv=(Q(1)+Q(8))/rhoG
            y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
            pG=this%gas%get_p_from_rho_e(rho=rhoG,e=eG,y=y)
            TG=this%gas%get_T_from_rho_e(rho=rhoG,e=eG,y=y)
            pV=this%get_xv(Yv)*pG
            call this%get_Tsat(pG,pV,TG,Tsat,conv,it)
            if (dbg_cell) then
               print*,'rhoG=',rhoG,'eG=',eG,'TG=',TG,', pG=',pG,', pV=',pV,', Yv=',Yv
               print*,'Tsat=',Tsat,', conv=',conv,', it=',it
            end if
            if (conv.and.TG.ge.Tsat) then
               Q(4)=Q(3)+Q(4)+Pjump*VF
               VF=0.0_WP
               Q(2)=Q(1)+Q(2)
               Q(8)=Q(1)+Q(8)
               Q(1)=0.0_WP
               Q(3)=0.0_WP
               if (present(ierr)) ierr=RELAX_OK
               if (dbg_cell) print*,'pure gas is stable'
               return
            else
               if (dbg_cell) print*,'pure gas not stable: conv=',conv,' TG=',TG,' Tsat=',Tsat
            end if
         else if (VF.gt.pure_VFhi) then
            ! Near-pure liquid. This path folds the WHOLE gas phase into the liquid, so it is only
            ! valid when that gas is essentially all condensable: a cell carrying real air keeps a
            ! gas phase however little room it has, and falls through to the general solve, which
            ! retains the air (see boundary_all_liquid). The old gate was Q(8).eq.Q(2) -- an exact
            ! floating-point equality on two independently advected variables, which almost never
            ! holds, and which sat on the OUTER branch so a near-pure-liquid cell with any air
            ! never even reached this VF test.
            if (Q(2).gt.0.0_WP.and.Q(8).lt.(1.0_WP-Y_small)*Q(2)) then
               if (dbg_cell) print*,'near-pure liquid but gas carries inert air -- falling through to general chem_relax'
            else
               if (dbg_cell) print*,'attempting pure liquid'
               ! Transfer all vapor mass and energy to liquid
               rhoL=Q(1)+Q(2)
               eL=(Q(3)+Q(4)-Pjump*(1.0_WP-VF))/rhoL
               pL=this%liq%get_p_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
               TL=this%liq%get_T_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
               ! Liquid compressed above its saturation pressure is stable -- no vapor survives.
               ! Testing pL against get_pvsat(pL,TL) puts p on both sides of the comparison, and
               ! pvsat grows exponentially in p, so that test stops being monotone in pressure and
               ! calls strongly compressed liquid superheated (see get_psat).
               call this%get_psat(TL,pVsat,conv)
               if (conv.and.pL.ge.pVsat) then
                  Q(3)=Q(3)+Q(4)-Pjump*(1.0_WP-VF)
                  VF=1.0_WP
                  Q(1)=Q(1)+Q(2)
                  Q(2)=0.0_WP
                  Q(8)=0.0_WP
                  Q(4)=0.0_WP
                  if (present(ierr)) ierr=RELAX_OK
                  if (dbg_cell) print*,'pure liquid is stable: pL=',pL,' psat=',pVsat
                  return
               else
                  if (dbg_cell) print*,'pure liquid not stable: pL=',pL,' psat=',pVsat,' conv=',conv
               end if
            end if
         else
            if (dbg_cell) print*,'falling through to general chemical relaxation'
         end if
      end block pure_phase_stability
      ! Store input
      VFin=VF
      allocate(Qin(size(Q))); Qin=Q
      if (Qin(2).gt.0.0_WP) then
         Yvin=Qin(7+this%liq%ns+this%indV-1)/Qin(2)
         Yvin=max(Yvmin,min(Yvmax,Yvin))
      else
         Yvin=0.0_WP
      end if
      ! Near pure cell flag
      near_pure=(VF.lt.VF_nuc).or.(VF.gt.1.0_WP-VF_nuc)
      ! Nucleation flags
      nucleated=.false.
      cavitated=.false.
      condensed=.false.
      ! Nucleation: seed a tiny opposite phase in metastable pure-ish cells so pTg starts well-conditioned.
      ! Gated by do_nucleate -- off defers cavitation/condensation of pure cells (first phase-change pass).
      ! No nucleation at nu=0: the seed is itself a phase change
      if (this%do_nucleate.and.near_pure.and.this%nu.gt.0.0_WP) then
         nucleation: block
            real(WP), parameter :: dT_nuc_max=0.05_WP ! max self-heating (relative to TL_nuc/TG_nuc) a nucleation seed may cause
            real(WP) :: rhoL_nuc,pL_nuc,TL_nuc,pv_sat,rhoV_nuc,eV_nuc
            real(WP) :: rhoG_nuc,pG_nuc,TG_nuc,Yv_nuc,xv_nuc,pv_nuc,Tsat_nuc
            real(WP) :: rhoL_new,eL_new,eG_nuc,cvG_nuc,drho,de
            real(WP) :: y_nuc(this%gas%ns),Yv_hi_nuc
            logical  :: conv_nuc
            integer  :: Tsat_it_nuc
            ! Vapor mass fraction in gas phase
            if (Q(2).gt.0.0_WP) then
               Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
               Yv=max(Yvmin,min(Yvmax,Yv))
            else
               Yv=0.0_WP
            end if
            if (VF.gt.1.0_WP-VF_nuc) then
               if (dbg_cell) print*,'Near pure liquid, attempting cavitation nucleation'
               ! Near-pure liquid: check cavitation
               rhoL_nuc=Q(1)/VF
               pL_nuc=this%liq%get_p_from_rho_e(rho=rhoL_nuc,e=Q(3)/Q(1),y=[1.0_WP])
               TL_nuc=this%liq%get_T_from_p_rho(p=pL_nuc,rho=rhoL_nuc,y=[1.0_WP])
               if (dbg_cell) then
                  print*,'rhoL_nuc=',rhoL_nuc,', pL_nuc=',pL_nuc,', TL_nuc=',TL_nuc
               end if
               ! if (pL_nuc.le.-this%liq%pinf.or.TL_nuc.le.0.0_WP) exit nucleation ! debug
               if (pL_nuc.le.-this%liq%pinf.or.TL_nuc.le.0.0_WP) then
                  if (dbg_cell) print*,'cavitation rejected: liquid state unphysical (pL_nuc<=-pinf or TL_nuc<=0)'
                  exit nucleation
               end if
               pv_sat=this%get_pvsat(pL_nuc,TL_nuc)
               if (dbg_cell) then
                  print*,'pv_sat=',pv_sat,', p_cav=',this%p_cav
               end if
               ! if (pv_sat.le.pL_nuc) exit nucleation  ! stable ! debug
               if (pv_sat.le.pL_nuc) then
                  if (dbg_cell) print*,'cavitation rejected: liquid already stable (pv_sat<=pL_nuc, not superheated)'
                  exit nucleation
               end if
               ! if (pL_nuc.gt.this%p_cav) exit nucleation  ! metastable but not deep enough to nucleate ! debug
               if (pL_nuc.gt.this%p_cav) then
                  if (dbg_cell) print*,'cavitation rejected: metastable but pL_nuc>p_cav, not deep enough yet'
                  exit nucleation
               end if
               ! Skip if too little liquid
               if (Q(1)/sum(Q(1:2)).lt.Y_small) then
                  if (dbg_cell) print*,'too little liquid, rejecting cavitation.'
                  exit nucleation
               end if
               ! Superheated: estimate nucleated vapor state
               y_nuc=0.0_WP; y_nuc(this%indV)=1.0_WP
               rhoV_nuc=this%gas%get_rho_from_p_T(p=pv_sat,T=TL_nuc,y=y_nuc)
               eV_nuc  =this%gas%get_e_from_p_T  (p=pv_sat,T=TL_nuc,y=y_nuc)
               drho=VF_nuc*rhoV_nuc
               if (dbg_cell) then
                  print*,'rhoV_nuc=',rhoV_nuc,', eV_nuc=',eV_nuc,', drho=',drho
               end if
               drho=min(drho,0.5_WP*Q(1))
               if (dbg_cell) print*,'capped drho=',drho
               ! Cap further so the seed's own latent-heat release can't overheat a sparse liquid phase
               if ((Q(3)/Q(1)).gt.eV_nuc) then
                  drho=min(drho,dT_nuc_max*TL_nuc*this%liq%cv*Q(1)/(Q(3)/Q(1)-eV_nuc))
                  if (dbg_cell) print*,'further capped drho=',drho
               end if
               ! Cap to the room the stretched liquid frees at pv_sat (isothermal NASG: (v-b)(p+pinf) fixed)
               drho=min(drho,rhoV_nuc*(VF-this%liq_nasg%b*Q(1))*(pv_sat-pL_nuc)/(pv_sat+this%liq%pinf))
               if (dbg_cell) print*,'tension-room capped drho=',drho,' room=',(VF-this%liq_nasg%b*Q(1))*(pv_sat-pL_nuc)/(pv_sat+this%liq%pinf) ! debug
               ! if (drho.le.0.0_WP) exit nucleation ! debug
               if (drho.le.0.0_WP) then
                  if (dbg_cell) print*,'cavitation rejected: capped nucleation mass drho<=0, nothing to seed'
                  exit nucleation
               end if
               de=drho*eV_nuc
               Q(1)=Q(1)-drho; Q(2)=Q(2)+drho
               Q(3)=Q(3)-de;   Q(4)=Q(4)+de
               Q(7+this%liq%ns+this%indV-1)=Q(7+this%liq%ns+this%indV-1)+drho
               ! VF=VF-VF_nuc
               Q(4)=Q(4)-Pjump*(Q(1)/rhoL_nuc-VF)
               VF=Q(1)/rhoL_nuc
               nucleated=.true.; cavitated=.true.
               if (dbg_cell) print*,'Cavitating: drho=',drho,', de=',de
            else
               if (dbg_cell) print*,'Near pure gas, attempting condensation nucleation'
               ! Near-pure gas: check condensation
               ! if (Q(2).le.0.0_WP) exit nucleation ! debug
               if (Q(2).le.0.0_WP) then
                  if (dbg_cell) print*,'condensation rejected: no gas mass present (Q(2)<=0)'
                  exit nucleation
               end if
               Yv_nuc=Q(7+this%liq%ns+this%indV-1)/Q(2); Yv_nuc=max(Yvmin,min(Yvmax,Yv_nuc))
               ! Skip if too little vapor
               Yv_hi_nuc=min(Yvmax,1.0_WP-(1.0_WP-Yv_nuc)*Q(2)/(Q(1)+Q(2)))
               if (dbg_cell) print*,'Yv_nuc=',Yv_nuc,', Yv_hi_nuc=',Yv_hi_nuc
               ! if ((Yv_nuc.lt.Y_small).or.(Yv_hi_nuc.lt.Y_small)) exit nucleation ! debug
               if ((Yv_nuc.lt.Y_small).or.(Yv_hi_nuc.lt.Y_small)) then
                  if (dbg_cell) print*,'condensation rejected: too little vapor headroom (Yv_nuc or Yv_hi_nuc < Y_small)'
                  exit nucleation
               end if
               y_nuc(this%indV)=Yv_nuc; y_nuc(this%indA)=1.0_WP-Yv_nuc
               rhoG_nuc=Q(2)/max(1.0_WP-VF,tiny(1.0_WP))
               pG_nuc=this%gas%get_p_from_rho_e(rho=rhoG_nuc,e=Q(4)/Q(2),y=y_nuc)
               TG_nuc=this%gas%get_T_from_p_rho(p=pG_nuc,rho=rhoG_nuc,y=y_nuc)
               if (dbg_cell) print*,'rhoG_nuc=',rhoG_nuc,', pG_nuc=',pG_nuc,', TG_nuc=',TG_nuc
               ! if (pG_nuc.le.0.0_WP.or.TG_nuc.le.0.0_WP) exit nucleation ! debug
               if (pG_nuc.le.0.0_WP.or.TG_nuc.le.0.0_WP) then
                  if (dbg_cell) print*,'condensation rejected: gas state unphysical (pG_nuc<=0 or TG_nuc<=0)'
                  exit nucleation
               end if
               xv_nuc=this%get_xv(Yv_nuc); pv_nuc=xv_nuc*pG_nuc
               if (dbg_cell) print*,'xv_nuc=',xv_nuc,', pv_nuc=',pv_nuc
               ! if (.not.check_pv(pv_nuc)) exit nucleation ! debug
               if (.not.check_pv(pv_nuc)) then
                  if (dbg_cell) print*,'condensation rejected: vapor partial pressure pv_nuc too small (check_pv failed)'
                  exit nucleation
               end if
               call this%get_Tsat(pG_nuc,pv_nuc,TG_nuc,Tsat_nuc,conv_nuc,Tsat_it_nuc)
               if (dbg_cell) print*,'saturation calculation: Tsat_nuc=',Tsat_nuc,', conv_nuc=',conv_nuc,', Tsat_it_nuc=',Tsat_it_nuc
               ! if (.not.conv_nuc) exit nucleation ! debug
               if (.not.conv_nuc) then
                  if (dbg_cell) print*,'condensation rejected: get_Tsat saturation solve did not converge'
                  exit nucleation
               end if
               ! if (TG_nuc.ge.Tsat_nuc) exit nucleation  ! stable pure vapor/gas ! debug
               if (TG_nuc.ge.Tsat_nuc) then
                  if (dbg_cell) print*,'condensation rejected: gas already stable (TG_nuc>=Tsat_nuc, not supersaturated)'
                  exit nucleation
               end if
               ! if (TG_nuc.gt.Tsat_nuc-this%Tctol) exit nucleation  ! metastable but not deep enough to nucleate ! debug
               if (TG_nuc.gt.Tsat_nuc-this%Tctol) then
                  if (dbg_cell) print*,'condensation rejected: metastable but within Tctol of Tsat, not deep enough yet'
                  exit nucleation
               end if
               ! Supersaturated: nucleate tiny liquid
               rhoL_new=this%liq%get_rho_from_p_T(p=pG_nuc,T=TG_nuc,y=[1.0_WP])
               eL_new  =this%liq%get_e_from_p_T  (p=pG_nuc,T=TG_nuc,y=[1.0_WP])
               eG_nuc  =this%gas%get_e_from_p_T  (p=pG_nuc,T=TG_nuc,y=y_nuc)
               cvG_nuc =sum(y_nuc(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
               drho=VF_nuc*rhoL_new
               if (dbg_cell) print*,'drho=',drho
               drho=min(drho,0.5_WP*Q(7+this%liq%ns+this%indV-1),0.5_WP*Q(2))
               if (dbg_cell) then
                  print*,'drho capped=',drho,', rhoL_new=',rhoL_new,', eL_new=',eL_new,', eG_nuc=',eG_nuc
               end if
               ! Cap further so the seed's own latent-heat release can't overheat a sparse vapor phase
               if (eG_nuc.gt.eL_new) then
                  drho=min(drho,dT_nuc_max*TG_nuc*cvG_nuc*Q(2)/(eG_nuc-eL_new))
                  if (dbg_cell) print*,'capping further, drho=',drho
                  ! Latent heat may warm the cell by at most half its subcooling, so the seed stays on the condensing side
                  drho=min(drho,0.5_WP*(Tsat_nuc-TG_nuc)*cvG_nuc*Q(2)/(eG_nuc-eL_new))
                  if (dbg_cell) print*,'subcooling cap, drho=',drho ! debug
               end if
               if (drho.le.0.0_WP) then
                  if (dbg_cell) print*,'drho negative; drho=',drho
                  exit nucleation
               end if
               de=drho*eL_new
               Q(1)=Q(1)+drho; Q(2)=Q(2)-drho
               Q(3)=Q(3)+de;   Q(4)=Q(4)-de
               Q(7+this%liq%ns+this%indV-1)=Q(7+this%liq%ns+this%indV-1)-drho
               ! VF=drho/rhoL_new
               Q(4)=Q(4)-Pjump*drho/rhoL_new
               VF=VF+drho/rhoL_new
               nucleated=.true.; condensed=.true.
               if (dbg_cell) print*,'Condensing: drho=',drho,', de=',de
            end if
         end block nucleation
         ! Proceed to next step only if the cell is either interfacial or successfully nucleated from a near pure state
         if ((.not.nucleated).and.(VF.lt.VFlo.or.VF.gt.VFhi)) then
            if(dbg_cell) print*,'Pure cell nucleation failed, returning.'
            if (present(ierr)) ierr=RELAX_NUC_FAILED
            return
         end if
      end if
      call this%pT_relax(dt,VF,Q,Pjump,ier)
      if (ier.ne.RELAX_OK) then
         if(dbg_cell) print*,'pT_relax failed ier=',ier
         if (nucleated) then
            call undo_nucleation()
            if (present(ierr)) ierr=RELAX_NUC_FAILED
         else
            if (present(ierr)) ierr=ier
         end if
         return
      end if
      ! A phase volume below round-off is left to the solver's pure-cell snap, as in the parent's apply
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         if (dbg_cell) print*,'pT_relax left VF outside [VFlo,VFhi], no chemical step: VF=',VF ! debug
         if (present(ierr)) ierr=RELAX_OK
         return
      end if
      ! Step 3: chemical (phase change)
      if (Q(2).gt.0.0_WP) then
         Yv=max(Yvmin,min(Yvmax,Q(7+this%liq%ns+this%indV-1)/Q(2)))
      else
         Yv=0.0_WP
      end if
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Recover p, T from dominant phase (p is always the liquid-side pressure convention
      ! used throughout activate_chem/solve_lv/solve_lvg below: p_l=p_g+Pjump post pT_relax)
      if (VF.gt.0.5_WP) then
         p=this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
         T=this%liq%get_T_from_p_rho(p=p,rho=Q(1)/VF,y=[1.0_WP])
      else
         p=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
         T=this%gas%get_T_from_p_rho(p=p,rho=Q(2)/(1.0_WP-VF),y=y)
         p=p+Pjump
      end if
      if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' pTg post-pT_relax VF=',VF,' p=',p,' TL=TG=',T,' Yv=',Yv ! debug
      if (dbg_cell) print*,'   RHOL=',Q(1)/VF,' RHOG=',merge(Q(2)/(1.0_WP-VF),-1.0_WP,VF.lt.1.0_WP),' Q=',Q ! debug
      ! Conserve totals
      VF0=VF
      allocate(Q0(size(Q))); Q0=Q
      rho0 =sum(Q0(1:2))
      if (Q0(2).gt.0.0_WP) then
         Yv0=Q0(7+this%liq%ns+this%indV-1)/Q0(2)
         Yv0=max(Yvmin,min(Yvmax,Yv0))
      else
         Yv0=0.0_WP
      end if
      if (nucleated) then
         ! Liquid wants to cavitate but there is not enough of it
         ! if (cavitated.and.(Qin(1)/rho0.lt.Y_small)) then
         if (cavitated.and.((Q0(1)/rho0.lt.Y_small).or.(Qin(1)/rho0.lt.Y_small))) then
            chem_relax=.false.
            ier=RELAX_NUC_FAILED
            if (dbg_cell) print*,'Liquid wants to cavitate but there is not enough of it'
            ! Vapor wants to condense but there is not enough of it
            ! else if (condensed.and.(Yvin.lt.Y_small)) then
         else if (condensed.and.((Yv0.lt.Y_small).or.(Yvin.lt.Y_small))) then
            chem_relax=.false.
            ier=RELAX_NUC_FAILED
            if (dbg_cell) print*,'Vapor wants to condense but there is not enough of it'
         else
            chem_relax=activate_chem(p,T,Yv)
         end if
      else if ((VF.lt.VFlo).or.(VF.gt.VFhi)) then
         chem_relax=.false.
         ! QUESTION: Do we need to set a flag or doing anything more here? Is there any ambigiouity here?
         ! ANSWER: No. This branch is only reached when .not.nucleated, so ier is left at
         ! pT_relax's success value (RELAX_OK) -- exactly the right signal below ("no
         ! chemical work needed, accept pT_relax's result"). No ambiguity: RELAX_NUC_FAILED
         ! is only ever set by the nucleated guards above, so it can't leak in here.
      else
         chem_relax=activate_chem(p,T,Yv)
      end if
      ! if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' pTg activate_chem chem_relax=',chem_relax,' p=',p,' TL=TG=',T,' Yv=',Yv,&
      ! &                    ' nucleated=',nucleated,' cavitated=',cavitated,' condensed=',condensed,' RHOL=',Q(1)/VF ! debug
      if (.not.chem_relax) then
         if (ier.eq.RELAX_OK) then
            ! Already relaxed
            call restore_pT_state()
            if (present(ierr)) ierr=RELAX_OK
         ! Nucleation happened but chemical relaxation cannot happen
         else if (nucleated) then
            call undo_nucleation()
            if (present(ierr)) ierr=RELAX_NUC_FAILED
         ! Nucleation did not happen and chemical relaxation cannot happen; this is probably never reached
         else
            call restore_pT_state()
            if (present(ierr)) ierr=ier
         end if
         ! Return here
         call dealloc(); return
      end if
      ! Chemical step: Pelanti's exponential solution for any nu (nu=huge gives frac=1), no iterative solve
      finite_chemical: block
         logical :: done
         call finite_chem(done)
         if (dbg_cell) print*,'finite_chem done=',done,' VF=',VF
         if (.not.done) then
            if (nucleated) then
               call undo_nucleation()
               if (present(ierr)) ierr=RELAX_NUC_FAILED
            else
               call restore_pT_state()
               if (present(ierr)) ierr=RELAX_FAILED
            end if
            call dealloc(); return
         end if
         ! Confirm we had enough of the reducing phase to begin with
         if (.not.nucleated) then
            if ((VF.lt.VF0.and.Q0(1)/rho0.lt.Y_small).or.(VF.gt.VF0.and.Yv0.lt.Y_small)) then
               call restore_pT_state(); call dealloc()
               if (present(ierr)) ierr=RELAX_FAILED
               return
            end if
         end if
         call dealloc()
         if (present(ierr)) ierr=ier
      end block finite_chemical
   contains
      !> Finite-rate chemical relaxation (Pelanti 2022 Sec. 7.2.3, eqs. 71-83), gas = vapor+air mixture, surface energy conserved
      subroutine finite_chem(done)
         logical, intent(out) :: done
         real(WP) :: dm,VFeq,Ysat,Peq
         integer :: iVQ
         logical :: ok
         done=.false.
         iVQ=7+this%liq%ns+this%indV-1
         ! Transfer over dt at ** (eqs. 71-78, or the exact dry-gas step)
         call this%get_chem_transfer(dt,VF,Q,Pjump,dm,VFeq,Ysat,ok)
         if (dbg_cell) print*,'[finite_chem] ok=',ok,' dm=',dm,' (liquid gained; >0 condenses) VF=',VF,' VFeq=',VFeq,' Ysat=',Ysat ! debug
         if (dbg_cell) print*,'[finite_chem] vapor mass Q(iVQ)=',Q(iVQ),' gas mass Q(2)=',Q(2),' dm/Q(iVQ)=',dm/max(Q(iVQ),tiny(1.0_WP)) ! debug
         ! A non-finite transfer would pass every bounds test below, so it is rejected here
         if (.not.ok.or.dm.ne.dm.or.VFeq.ne.VFeq.or.abs(dm).gt.huge(1.0_WP)) return
         if (dm.eq.0.0_WP) then
            done=.true.
            return
         end if
         ! A phase that runs out (mass or volume) transfers completely; no eq. 83 clamp
         ! Liquid that no longer fits its new volume (b*rhoL at the packing knee) has run out too
         if (dm.lt.0.0_WP.and.(Q(1)+dm.le.0.0_WP.or.VFeq.le.0.0_WP.or.this%liq_nasg%b*(Q(1)+dm).ge.this%liq_nasg%brhomax*VFeq)) then
            ! Complete evaporation: all liquid joins the gas as vapor, its surface energy Pjump*VF goes to the gas
            if (dbg_cell) print*,'[finite_chem] COMPLETE EVAPORATION: Q(1)+dm=',Q(1)+dm,' VFeq=',VFeq ! debug
            Q(2)=Q(2)+Q(1)
            Q(iVQ)=Q(iVQ)+Q(1)
            Q(4)=Q(3)+Q(4)+Pjump*VF
            Q(1)=0.0_WP
            Q(3)=0.0_WP
            VF=0.0_WP
            done=.true.
            return
         end if
         if (dm.gt.0.0_WP.and.(Q(iVQ)-dm.le.0.0_WP.or.VFeq.ge.1.0_WP)) then
            if (dbg_cell) print*,'[finite_chem] vapor runs out: Q(iVQ)-dm=',Q(iVQ)-dm,' VFeq=',VFeq,' (dm/Q(iVQ)=',dm/Q(iVQ),')' ! debug
            if (Q(iVQ).ge.(1.0_WP-Y_small)*Q(2)) then
               ! Complete condensation of an (essentially) air-free gas: pure liquid, as the pure-liquid snap
               if (dbg_cell) print*,'[finite_chem] COMPLETE CONDENSATION (air-free gas): all gas folded into the liquid, VF=1' ! debug
               Q(1)=Q(1)+Q(2)
               Q(3)=Q(3)+Q(4)-Pjump*(1.0_WP-VF)
               Q(2)=0.0_WP
               Q(4)=0.0_WP
               Q(iVQ)=0.0_WP
               VF=1.0_WP
               done=.true.
               return
            end if
            ! Air remains: all vapor condenses, VF from eq. 78
            VFeq=VF+(VFeq-VF)*Q(iVQ)/dm
            dm=Q(iVQ)
            if (dbg_cell) print*,'[finite_chem] air remains: all vapor condenses, dm=',dm,' VFeq=',VFeq ! debug
         end if
         if (VFeq.le.0.0_WP.or.VFeq.ge.1.0_WP.or.Q(1)+dm.le.0.0_WP.or.Q(iVQ)-dm.lt.0.0_WP.or. &
         &   this%liq_nasg%b*(Q(1)+dm).ge.this%liq_nasg%brhomax*VFeq) then
            if (dbg_cell) print*,'[finite_chem] REJECTED: VFeq=',VFeq,' Q(1)+dm=',Q(1)+dm,' Q(iVQ)-dm=',Q(iVQ)-dm ! debug
            return
         end if
         if (dbg_cell) print*,'[finite_chem] PARTIAL transfer: dm=',dm,' VF ',VF,' -> ',VFeq ! debug
         ! Update masses, and the surface energy Pjump*dV_l
         Q(1)=Q(1)+dm
         Q(2)=Q(2)-dm
         Q(iVQ)=Q(iVQ)-dm
         Q(4)=Q(4)-Pjump*(VFeq-VF)
         VF=VFeq
         ! Pressure from the mixture energy relation at the new masses and composition (eq. 82)
         Yv=max(Yvmin,min(Yvmax,Q(iVQ)/Q(2)))
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
         cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
         qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
         gammaG=cpG/cvG
         Peq=this%get_p_eq(VF,Q,qG,gammaG,Pjump)
         if (dbg_cell) print*,'[finite_chem] after transfer: Peq=',Peq,' RHOL=',Q(1)/VF,' RHOG=',Q(2)/(1.0_WP-VF),' Yv=',Yv ! debug
         if (Peq.le.-this%liq%pinf.or.Peq-Pjump.le.0.0_WP) return
         Q(3)=(       VF)*this%liq%get_rhoe_from_p_rho(p=Peq,      rho=Q(1)/(       VF),y=[1.0_WP])
         Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=Peq-Pjump,rho=Q(2)/(1.0_WP-VF),y=y)
         if (dbg_cell) print*,'[finite_chem] after transfer: TL=',this%liq%get_T_from_p_rho(p=Peq,rho=Q(1)/VF,y=[1.0_WP]),' TG=',this%gas%get_T_from_p_rho(p=Peq-Pjump,rho=Q(2)/(1.0_WP-VF),y=y) ! debug
         done=.true.
      end subroutine finite_chem
      subroutine restore_pT_state()
         VF=VF0; Q=Q0
      end subroutine restore_pT_state
      subroutine undo_nucleation()
         VF=VFin; Q=Qin
      end subroutine undo_nucleation
      subroutine dealloc()
         if (allocated(Q0)) deallocate(Q0)
         if (allocated(y))  deallocate(y)
         if (allocated(Qin))deallocate(Qin)
      end subroutine dealloc
      logical function check_pv(pv_)
         real(WP), intent(in) :: pv_
         check_pv=(pv_.gt.p_eps)
      end function check_pv
      logical function activate_chem(p_,T_,Yv_)
         real(WP), intent(in)    :: p_,T_
         real(WP), intent(in)    :: Yv_
         real(WP) :: xv,pv_,Fsat
         if (dbg_cell) print*,'--------------------------------------------------'
         if (dbg_cell) print*,'inside activate_chem'
         ! Default to false
         activate_chem=.false.
         ! Get vapor mole fraction and partial pressure (gas-side pressure is p_-Pjump)
         xv=this%get_xv(Yv_); pv_=xv*(p_-Pjump)
         if (pv_.lt.this%pv_min) then
            ! Dry gas: F_C=+inf, liquid always evaporates (finite_chem integrates the singular term)
            if (dbg_cell) print*,'pv too small (dry gas), heading to chemical relaxation.'
         else
            ! Direct saturation residual
            Fsat=this%pTsat(p_,pv_,T_)
            if (dbg_cell) print*,'pv normal, Fsat=',Fsat
            if (abs(Fsat).lt.this%F1_tol) then
               ier=RELAX_OK
               if (dbg_cell) print*,'Already saturated, no phase change needed, RELAX_OK.'
               return
            else
               if (dbg_cell) print*,'Everything is fine and heading to chemical relaxation.'
            end if
         end if
         ! Chemical relaxation should fire
         activate_chem=.true.
      end function activate_chem
   end subroutine pTg_relax

   !> Predicted chemical transfer over dt from the state (VF,Q), which is not modified: liquid mass gained dm
   !> (Pelanti eqs. 71-77, or the exact dry-gas step), linearized VFeq (eq. 78) and saturation Ysat at the cell p,T.
   subroutine get_chem_transfer(this,dt,VF,Q,Pjump,dm,VFeq,Ysat,ok)
      implicit none
      class(relax_igmix_nasg_finite), intent(in)  :: this
      real(WP),                       intent(in)  :: dt
      real(WP),                       intent(in)  :: VF
      real(WP), dimension(:),         intent(in)  :: Q
      real(WP),                       intent(in)  :: Pjump
      real(WP),                       intent(out) :: dm,VFeq,Ysat
      logical,                        intent(out) :: ok
      real(WP) :: RHOL,RHOG,EL,EG,PL,PG,TL,TG,CL,CG,GamL,GamG,HL,HG,chiL,chiG
      real(WP) :: PHIL,PHIG,zetaL,zetaG,CEp,CET,DTP,SVF,ST,Sp
      real(WP) :: YvS,xv,pv,RV,RA,Rm,sL,sV,g1,g2,KCC,xi1inv,xi2inv,dg,Kdt,frac,lnpvs,xs
      real(WP) :: mAir,ms,u0,u1,tgt
      real(WP), dimension(3,3) :: M,Mc
      real(WP), dimension(3)   :: rhs
      real(WP), dimension(this%gas%ns) :: yS
      integer :: iVQ
      logical :: dry
      dm=0.0_WP
      VFeq=VF
      Ysat=0.0_WP
      ok=.false.
      iVQ=7+this%liq%ns+this%indV-1
      if (this%nu.le.0.0_WP) then; ok=.true.; return; end if
      ! Phasic state at ** (each phase at its own p and T); EL, EG are internal energies per unit volume
      RHOL=Q(1)/(       VF); EL=Q(3)/(       VF)
      RHOG=Q(2)/(1.0_WP-VF); EG=Q(4)/(1.0_WP-VF)
      YvS=max(0.0_WP,min(1.0_WP,Q(iVQ)/Q(2)))
      yS=0.0_WP; yS(this%indV)=YvS; yS(this%indA)=1.0_WP-YvS
      PL=this%liq%get_p_from_rho_e(rho=RHOL,e=EL/RHOL,y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=RHOG,e=EG/RHOG,y=yS)
      TL=this%liq%get_T_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
      TG=this%gas%get_T_from_p_rho(p=PG,rho=RHOG,y=yS)
      if (PL.le.-this%liq%pinf.or.PG.le.0.0_WP.or.TL.le.0.0_WP.or.TG.le.0.0_WP) return
      ! Vapor partial pressure p_v=x_v*p_g; dry gas (no vapor) is integrated separately below
      xv=this%get_xv(YvS); pv=xv*PG
      dry=(Q(iVQ).le.0.0_WP.or.pv.lt.this%pv_min)
      if (dbg_cell) print*,'[get_chem_transfer] state: VF=',VF,' RHOL=',RHOL,' RHOG=',RHOG,' PL=',PL,' PG=',PG,' TL=',TL,' TG=',TG ! debug
      if (dbg_cell) print*,'[get_chem_transfer] Yv=',YvS,' xv=',xv,' pv=',pv,' dry=',dry,' nu=',this%nu,' dt=',dt ! debug
      ! Phasic coefficients: Gamma, c, chi=c^2-Gamma*h, phi=(drho/dT)_p, zeta=(drho/dp)_T (eqs. 16, 34)
      GamL=this%liq%get_gruneisen_from_rho_e(rho=RHOL,e=EL/RHOL,y=[1.0_WP])
      GamG=this%gas%get_gruneisen_from_rho_e(rho=RHOG,e=EG/RHOG,y=yS)
      CL=this%liq%get_c_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
      CG=this%gas%get_c_from_p_rho(p=PG,rho=RHOG,y=yS)
      HL=(EL+PL)/RHOL
      HG=(EG+PG)/RHOG
      chiL=CL**2-GamL*HL
      chiG=CG**2-GamG*HG
      if (this%liq_nasg%b*RHOL.gt.this%liq_nasg%brhomax) then
         PHIL =-RHOL/TL
         zetaL= RHOL/(PL+this%liq%pinf)
      else
         PHIL =-RHOL*(1.0_WP-this%liq_nasg%b*RHOL)/TL
         zetaL= RHOL*(1.0_WP-this%liq_nasg%b*RHOL)/(PL+this%liq%pinf)
      end if
      PHIG =-RHOG/TG
      zetaG= RHOG/PG
      ! C_Ep=sum alpha_k (dE_k/dT)_p, C_ET=sum alpha_k (dE_k/dp)_T (eq. 101)
      CEp=VF*(-chiL*PHIL/GamL)+(1.0_WP-VF)*(-chiG*PHIG/GamG)
      CET=VF*(1.0_WP-chiL*zetaL)/GamL+(1.0_WP-VF)*(1.0_WP-chiG*zetaG)/GamG
      ! dq/dw of eq. 97 for w=[alpha_1,T,p]; surface energy adds Pjump to the alpha_1 entry of the energy row
      M(1,:)=[ RHOL       ,         VF *PHIL,          VF*zetaL]
      M(2,:)=[-RHOG       , (1.0_WP-VF)*PHIG, (1.0_WP-VF)*zetaG]
      M(3,:)=[ EL-EG+Pjump, CEp             , CET              ]
      ! Right-hand side of eq. 100: [1,-1,0] for pure vapor; transferred mass enters the mixture as vapor
      RV=this%gas%cp(this%indV)-this%gas%cv(this%indV)
      RA=this%gas%cp(this%indA)-this%gas%cv(this%indA)
      Rm=YvS*RV+(1.0_WP-YvS)*RA
      rhs(1)=1.0_WP
      rhs(2)=-RV/Rm
      rhs(3)=(1.0_WP-YvS)*(-(RV-RA)*(EG/RHOG)/Rm+(this%gas%cv(this%indV)-this%gas%cv(this%indA))*TG+this%gas%q(this%indV)-this%gas%q(this%indA))
      ! S_alpha, S_T, S_p by Cramer's rule (eqs. 23)
      DTP=det3(M)
      if (abs(DTP).le.tiny(1.0_WP)) return
      Mc=M; Mc(:,1)=rhs; SVF=det3(Mc)/DTP
      Mc=M; Mc(:,2)=rhs; ST =det3(Mc)/DTP
      Mc=M; Mc(:,3)=rhs; Sp =det3(Mc)/DTP
      if (dbg_cell) print*,'[get_chem_transfer] S_alpha=',SVF,' S_T=',ST,' S_p=',Sp,' (d alpha_l, dT, dp per unit dm)' ! debug
      if (dry) then
         ! Dry gas: K_CC->inf and eq. 77 gives dm->0, yet the ODE evaporates. Keep only the singular mixing
         ! term F_C=Rv*T*ln(m_s/m_v) and integrate it exactly: E1(ln(m_s/m_v))=E1(ln(m_s/m_v**))+nu*Rv*T*dt/m_s
         g1=this%liq%get_g_from_p_T(p=PL,T=TL,y=[1.0_WP])
         g2=this%gas%cp(this%indV)*TG+this%gas%q(this%indV)-TG*(this%gas%cp(this%indV)*log(TG)+this%gas%qp(this%indV))
         ! m_s: vapor mass at which F_C=0, from p_vs=exp((g_l-g_v0)/(Rv*T)) and the dilute x_v=(Ma/Mv)*m_v/m_a
         mAir=Q(2)-max(Q(iVQ),0.0_WP)
         ms=(Mv/Ma)*mAir*exp((g1-g2)/(RV*TG)-log(PG))
         u0=max(Q(iVQ),0.0_WP)/ms
         if (u0.ge.1.0_WP) then
            ok=.true.
            return
         end if
         ! u1 from E1(ln(1/u1))=E1(ln(1/u0))+tau, tau=nu*Rv*T*dt/m_s; u1=1 when tau exceeds E1(1e-300)~690
         if (this%nu.ge.huge(1.0_WP).or.log(this%nu)+log(RV*TG*dt/ms).gt.log(690.0_WP)) then
            u1=1.0_WP
         else
            tgt=this%nu*RV*TG*dt/ms
            if (u0.gt.0.0_WP) tgt=tgt+exp(log_E1(-log(u0)))
            if (tgt.gt.690.0_WP) then
               u1=1.0_WP
            else
               u1=exp(-inv_log_E1(log(tgt)))
            end if
         end if
         dm=-(max(u1,u0)-u0)*ms
         lnpvs=(g1-g2)/(RV*TG)
         if (dbg_cell) print*,'[get_chem_transfer] DRY branch: g1=',g1,' g2(p0)=',g2,' m_s=',ms,' u0=',u0,' u1=',u1,' dm=',dm ! debug
      else
         ! Phasic entropies and chemical potentials; vapor at its partial pressure
         sL=this%liq%get_s_from_p_T(p=PL,T=TL,y=[1.0_WP])
         sV=this%gas%cp(this%indV)*log(TG)-RV*log(pv)+this%gas%qp(this%indV)
         g1=this%liq%get_g_from_p_T(p=PL,T=TL,y=[1.0_WP])
         g2=this%gas%cp(this%indV)*TG+this%gas%q(this%indV)-TG*sV
         ! Inverse chemical capacities (eq. 72, corrected); 1/rho_2 -> Rv*T/p_g and s_2 -> s_v(p_v,T), plus composition term K_CC
         KCC=RV*TG*(1.0_WP-xv)/Q(iVQ)
         xi1inv= Sp/RHOL-sL*ST
         xi2inv=-(Sp*RV*TG/PG-sV*ST)+KCC
         if (xi1inv+xi2inv.le.0.0_WP) return
         ! Fraction 1-exp(-Kg*dt) with Kg=nu*(1/xi_1+1/xi_2) (eq. 76); instantaneous when nu=huge or nu*Kdt>=50
         Kdt=(xi1inv+xi2inv)*dt
         frac=1.0_WP
         if (this%nu.lt.huge(1.0_WP).and.log(this%nu)+log(Kdt).lt.log(50.0_WP)) frac=1.0_WP-exp(-this%nu*Kdt)
         ! Partial density (eq. 77), Delta g=g_2-g_1
         dg=g2-g1
         dm=dg/(xi1inv+xi2inv)*frac
         lnpvs=(g1-g2)/(RV*TG)+log(pv)
         if (dbg_cell) print*,'[get_chem_transfer] g_liq=',g1,' g_vap=',g2,' dg=g_vap-g_liq=',dg,' (>0 condenses)' ! debug
         if (dbg_cell) print*,'[get_chem_transfer] sL=',sL,' sV=',sV,' KCC=',KCC,' 1/xi1=',xi1inv,' 1/xi2=',xi2inv,' nu*Kdt=',this%nu*Kdt,' frac=',frac ! debug
         if (dbg_cell) print*,'[get_chem_transfer] dm=dg/(1/xi1+1/xi2)*frac=',dm,' psat(T)~',exp(lnpvs),' vs pv=',pv ! debug
      end if
      ! Volume fraction (eq. 78)
      VFeq=VF+SVF*dm
      ! Saturation vapor mass fraction at p_g: x_vs=p_vs/p_g (1 if the liquid boils at this p,T)
      xs=exp(lnpvs-log(PG))
      Ysat=1.0_WP
      if (xs.lt.1.0_WP) Ysat=xs*Mv/(xs*Mv+(1.0_WP-xs)*Ma)
      if (dbg_cell) print*,'[get_chem_transfer] VFeq=VF+S_alpha*dm=',VFeq,' Ysat=',Ysat ! debug
      ok=.true.
   end subroutine get_chem_transfer

   !> 3x3 determinant
   pure real(WP) function det3(A)
      implicit none
      real(WP), dimension(3,3), intent(in) :: A
      det3=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
   end function det3

   !> ln E1(x), E1(x)=int_x^inf exp(-t)/t dt: series for x<=1, continued fraction (modified Lentz) for x>1
   pure real(WP) function log_E1(x) result(lE)
      implicit none
      real(WP), intent(in) :: x
      real(WP), parameter :: euler=0.5772156649015329_WP
      real(WP) :: s,term,b,c,d,h,an,del
      integer :: k
      if (x.le.1.0_WP) then
         ! E1=-euler-ln(x)-sum_k (-x)^k/(k*k!)
         s=0.0_WP
         term=1.0_WP
         do k=1,100
            term=-term*x/real(k,WP)
            s=s+term/real(k,WP)
            if (abs(term).lt.1.0e-2_WP*epsilon(1.0_WP)) exit
         end do
         lE=log(-euler-log(x)-s)
      else
         b=x+1.0_WP
         c=huge(1.0_WP)
         d=1.0_WP/b
         h=d
         do k=1,200
            an=-real(k,WP)**2
            b=b+2.0_WP
            d=1.0_WP/(an*d+b)
            c=b+an/c
            del=c*d
            h=h*del
            if (abs(del-1.0_WP).lt.epsilon(1.0_WP)) exit
         end do
         lE=log(h)-x
      end if
   end function log_E1

   !> Root x>0 of ln E1(x)=lt (E1 decreases from +inf to 0): Newton in ln x, bisection-safeguarded
   real(WP) function inv_log_E1(lt) result(x)
      implicit none
      real(WP), intent(in) :: lt
      real(WP), parameter :: euler=0.5772156649015329_WP
      real(WP) :: ylo,yhi,y,g,dg,ynew,lE
      integer :: it
      ylo=log(1.0e-300_WP)
      yhi=log(745.0_WP)
      ! Initial guess from the asymptotes E1~-euler-ln(x) (x->0) and E1~exp(-x)/x (x->inf)
      if (lt.gt.0.0_WP) then
         y=log(max(exp(-euler-exp(lt)),1.0e-300_WP))
      else
         y=log(max(-lt-log(max(-lt,1.0_WP)),1.0e-3_WP))
      end if
      do it=1,100
         x=exp(y)
         lE=log_E1(x)
         g=lE-lt
         if (g.gt.0.0_WP) then
            ylo=y
         else
            yhi=y
         end if
         ! d(ln E1)/d(ln x)=-exp(-x)/E1
         dg=-exp(-x-lE)
         ynew=y-g/dg
         if (ynew.le.ylo.or.ynew.ge.yhi) ynew=0.5_WP*(ylo+yhi)
         if (abs(ynew-y).lt.1.0e-14_WP*max(1.0_WP,abs(y))) then
            y=ynew
            exit
         end if
         y=ynew
      end do
      x=exp(y)
   end function inv_log_E1

   !> NASG-form energy-conserving equilibrium pressure at given (frozen) VF, liquid pressure PL
   !> with PL-PG=Pjump_ (co-volume clamped consistently with the nasg accessors)
   real(WP) function get_p_eq(this,VF_,Q0_,qG_,gammaG_,Pjump_) result(p_eq)
      implicit none
      class(relax_igmix_nasg_finite), intent(in) :: this
      real(WP),                intent(in) :: VF_
      real(WP), dimension(:),  intent(in) :: Q0_
      real(WP),                intent(in) :: qG_,gammaG_,Pjump_
      real(WP) :: one_brho
      one_brho=max(1.0_WP-this%liq_nasg%b*Q0_(1)/VF_,1.0_WP-this%liq_nasg%brhomax)
      p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-one_brho*VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_+&
      &     (1.0_WP-VF_)*Pjump_/(gammaG_-1.0_WP))/&
      &    (one_brho*VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function get_p_eq

   !> Equilibrium T from energy conservation (NASG form: includes co-volume b correction)
   real(WP) function get_T_lvg(this,p_,Yv_,rho0,rhoA0,Pjump_) result(T)
      implicit none
      class(relax_igmix_nasg_finite), intent(in) :: this
      real(WP),                intent(in) :: p_,Yv_,rho0,rhoA0,Pjump_
      ! debug: clamp covolume term like nasg_class does
      ! T=(1.0_WP-Yv_-this%liq_nasg%b*(rho0*(1.0_WP-Yv_)-rhoA0))/                                                       &
      ! & ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)                          +  &
      ! &   rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_                                     +  &
      ! &          (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/(p_-Pjump_))
      T=(1.0_WP-Yv_-min(this%liq_nasg%b*(rho0*(1.0_WP-Yv_)-rhoA0),this%liq_nasg%brhomax))/                             &
      & ((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)                          +  &
      &   rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_                                     +  &
      &          (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/(p_-Pjump_))
   end function get_T_lvg

   !> Quadratic coefficients for the equilibrium-T equation (NASG form: includes co-volume terms)
   subroutine get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
      implicit none
      class(relax_igmix_nasg_finite), intent(in)  :: this
      real(WP),                intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG,Pjump
      real(WP),                intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      real(WP) :: cvV_,gammaV_,qV_,ombm_lv
      cvV_   =this%gas%cv   (this%indV)
      gammaV_=this%gas%gamma(this%indV)
      qV_    =this%gas%q    (this%indV)
      ! debug: clamp covolume term like nasg_class does
      ombm_lv=max(1.0_WP-rho0*this%liq_nasg%b,1.0_WP-this%liq_nasg%brhomax)
      ap=rho0*this%liq%cv*cvV_*((gammaV_-this%liq%gamma)*p_eq+this%liq%gamma*(gammaV_-1.0_WP)*this%liq%pinf       +&
      &  (this%liq%gamma-1.0_WP)*Pjump)
      ! bp=(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq**2                                                    +&
      ! &  (this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)                                                          -&
      ! &   this%liq%gamma*this%liq%cv)+rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
      ! &   rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_eq                                    +&
      ! &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(rho0*this%liq%q-rhoe0)
      ! dp=p_eq*(p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
      ! dapdp=rho0*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
      ! dbpdp=2.0_WP*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%cv)*p_eq                                             +&
      ! &     this%liq%pinf*(cvV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%gamma*this%liq%cv)                            +&
      ! &     rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)                          +&
      ! &     rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)
      ! ddpdp=(2.0_WP*p_eq+this%liq%pinf)*(qV_*(1.0_WP-rho0*this%liq_nasg%b)-this%liq%q+this%liq_nasg%b*rhoe0)
      bp=(cvV_*ombm_lv-this%liq%cv)*p_eq**2                                                                          +&
      &  (this%liq%pinf*(cvV_*ombm_lv                                                                                -&
      &   this%liq%gamma*this%liq%cv)+rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)+&
      &   rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_))*p_eq                                    +&
      &   (gammaV_-1.0_WP)*cvV_*this%liq%pinf*(rho0*this%liq%q-rhoe0)                                                -&
      &   Pjump*(cvV_*ombm_lv-this%liq%cv)*p_eq-Pjump*this%liq%pinf*(cvV_*ombm_lv-this%liq%gamma*this%liq%cv)        +&
      &   Pjump*this%liq%cv*(this%liq%gamma-1.0_WP)*(rho0*qV_-rhoe0)
      dp=(p_eq-Pjump)*(p_eq+this%liq%pinf)*(qV_*ombm_lv-this%liq%q+this%liq_nasg%b*rhoe0)
      dapdp=rho0*this%liq%cv*cvV_*(gammaV_-this%liq%gamma)
      dbpdp=2.0_WP*(cvV_*ombm_lv-this%liq%cv)*p_eq                                                                   +&
      &     this%liq%pinf*(cvV_*ombm_lv-this%liq%gamma*this%liq%cv)                                                  +&
      &     rho0*((gammaV_-1.0_WP)*cvV_*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qV_)                          +&
      &     rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(gammaV_-1.0_WP)*cvV_)                                        -&
      &     Pjump*(cvV_*ombm_lv-this%liq%cv)
      ddpdp=(qV_*ombm_lv-this%liq%q+this%liq_nasg%b*rhoe0)*(2.0_WP*p_eq+this%liq%pinf-Pjump)
   end subroutine get_coeffs_lv

end module relax_igmix_nasg_finite_class
