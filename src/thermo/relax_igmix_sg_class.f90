!> SG-liquid + ideal-gas-mixture relaxation: mechanical (p), thermal (pT),
!> chemical/phase-change (pTg) with non-condensable gas
!> Liquid: class(stiffened_gas) (accepts SG or NASG via inheritance)
!> Gas: class(igmix) (multi-species mixture; ns=1 reduces to pure vapor)
!> Assumes vapor is the transported gas species at Q(7+liq%ns+indV-1); air (or other carrier) is the implicit ns-th species
module relax_igmix_sg_class
   use precision,             only: WP
   use messager,              only: die
   use thermorelax_class,     only: thermorelax,RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_VACUUM_GAS,RELAX_DEGENERATE,RELAX_NUC_FAILED,RELAX_VACUUM_VAPOR
   use stiffened_gas_class,   only: stiffened_gas
   use igmix_class,           only: igmix
   implicit none
   private

   public :: relax_igmix_sg
   public :: Prelax,PTrelax,PTgrelax,PThybrid
   public :: Mv,Ma
   ! debug: (dbg_i,dbg_j,dbg_k) identifies the cell CURRENTLY being traced; caller sets it to
   ! the active cell's own indices and dbg_cell=.true. right before apply(), so every print
   ! inside this module (which has no i,j,k of its own) labels itself correctly even when the
   ! caller is sweeping a whole box of cells, not just one. The box itself is
   ! [dbg_ilo,dbg_ihi]x[dbg_jlo,dbg_jhi]x[dbg_klo,dbg_khi] (inclusive); an empty range (lo>hi,
   ! the default) matches no cell.
   logical, public :: dbg_cell=.false.
   integer, public :: dbg_i=-100000000 ! 60
   integer, public :: dbg_j=-100000000 ! 1032
   integer, public :: dbg_k=-100000000 ! 0
   integer, public :: dbg_ilo=100000000, dbg_ihi=-100000000
   integer, public :: dbg_jlo=100000000, dbg_jhi=-100000000
   integer, public :: dbg_klo=100000000, dbg_khi=-100000000

   !> Molar masses of vapor and air [kg/mol]
   real(WP), parameter :: Mv=0.0180153_WP
   real(WP), parameter :: Ma=0.02897_WP

   !> Model enum
   integer, parameter :: Prelax  =1   !< Mechanical only
   integer, parameter :: PTrelax =2   !< Mechanical + thermal
   integer, parameter :: PTgrelax=3   !< Mechanical + thermal + chemical/phase change
   integer, parameter :: PThybrid=4   !< Mechanical below Tratmax temperature contrast, mechanical+thermal above

   type, extends(thermorelax) :: relax_igmix_sg
      class(stiffened_gas), pointer :: liq => null()
      class(igmix),         pointer :: gas => null()
      !> Species indices for vapor and air in the gas mixture
      integer  :: indV =0
      integer  :: indA =0
      !> Saturation curve coefficients (computed in initialize from material properties)
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
      integer  :: Tsat_itmax=40
      integer  :: NR_itmax  =40
      !> Phase-change controls (defaults preserve the original dimensional behaviour; override per case)
      real(WP) :: pv_min     =1.0e-8_WP    !< Vapor partial-pressure floor for the dry-edge reseed in activate_chem
      logical  :: do_nucleate=.true.       !< Seed opposite phase in near-pure metastable cells (cavitation/condensation)
      real(WP) :: p_cav=huge(1.0_WP)       !< Cavitation delay: nucleate only when pL < p_cav (default: nucleate at p_sat)
      real(WP) :: Tctol=0.0_WP             !< Condensation temperature tolerance: nucleate only when TG < Tsat-Tctol (default: nucleate at Tsat)
      !> Dispatch
      integer  :: model=Prelax
      real(WP) :: Tratmax=10.0_WP          !< PThybrid: temperature contrast max(TG/TL,TL/TG) above which pT_relax is used
      real(WP) :: VFratmax=10.0_WP         !< Max per-call phase-volume change factor in p_relax (partial relax beyond)
      !> Case-owned stability guards (defaults are no-ops; a case opts in by setting these)
      real(WP) :: Pmin_liq=-1.0e30_WP   !< Liquid pressure floor (e.g. max sustainable tension); off by default
      real(WP) :: Tmin_liq=-1.0_WP      !< Liquid temperature floor; off by default
      real(WP) :: Pmin_gas=-1.0e30_WP   !< Gas pressure floor (e.g. ~saturation pressure); off by default
      real(WP) :: Tmin_gas=-1.0_WP      !< Gas temperature floor; off by default
      real(WP) :: diss_P=1.0e30_WP      !< Dissolution: absorb gas where phasic gas pressure exceeds this (1e30=off)
      real(WP) :: diss_RHO=1.0e30_WP    !< Dissolution: absorb gas where phasic gas density exceeds this (1e30=off)
      real(WP) :: vol=1.0_WP            !< Cell volume for ledger units (set by the case; apply() runs on the finest level only)
      !> Ledger: rank-local cumulative accumulators (1-2 dissolution n/dm; 3 quadratic proposal
      !> succeeded; 4 proposal failed, completed by the swap alone; 5-6 floor n/dE; 7 untouched
      !> cells, non-positive phase mass). The model only counts; reduction across ranks is the
      !> monitoring code's business (the case reduces acc itself)
      real(WP), dimension(7) :: acc=0.0_WP
   contains
      procedure :: initialize
      procedure :: debug_dump
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: get_T_lvg
      procedure :: get_coeffs_lv
      procedure :: get_p_eq
      procedure :: pTsat
      procedure :: dpTsatdT
      procedure :: dpTsatdp_lv
      procedure :: dpTsatdlnp
      procedure :: get_Tsat
      procedure :: get_psat
      procedure :: get_pvsat
      procedure :: get_xv
   end type relax_igmix_sg

contains

   !> Store EOS pointers, vapor/air species indices, and compute saturation-curve coefficients
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_sg),         intent(inout) :: this
      class(stiffened_gas),  target, intent(in)    :: liq
      class(igmix),          target, intent(in)    :: gas
      integer,                       intent(in)    :: indV,indA
      real(WP) :: cpV,cvV,RV
      this%liq=>liq
      this%gas=>gas
      this%indV=indV
      this%indA=indA
      ! Vapor properties from the gas mixture (flat-array access)
      cvV=gas%cv(indV)
      cpV=gas%cp(indV)
      RV =cpV-cvV
      ! Saturation curve coefficients (integrated Clapeyron form)
      this%AS=(liq%cp-cpV+gas%qp(indV)-liq%qp)/RV
      this%BS=(liq%q -gas%q (indV))            /RV
      this%CS=(cpV-liq%cp)                     /RV
      this%DS=(liq%cp-liq%cv)                  /RV
      ! ES stays 0 for SG; NASG override sets ES=liq%b/RV after calling parent initialize
   end subroutine initialize

   subroutine debug_dump(this,label,VF,Q,Pjump,ier)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      character(len=*),       intent(in)    :: label
      real(WP),                intent(in)   :: VF
      real(WP), dimension(:),  intent(in)   :: Q
      real(WP),                intent(in)   :: Pjump
      integer,  optional,      intent(in)   :: ier
      real(WP) :: RHOL,PL,TL,RHOG,PG,TG,Yv
      real(WP), dimension(this%gas%ns) :: y
      integer :: iVQ
      RHOL=0.0_WP; PL=0.0_WP; TL=0.0_WP
      if (VF.gt.0.0_WP) then
         RHOL=Q(1)/VF
         PL=this%liq%get_p_from_rho_e(rho=RHOL,e=Q(3)/Q(1),y=[1.0_WP])
         TL=this%liq%get_T_from_rho_e(rho=RHOL,e=Q(3)/Q(1),y=[1.0_WP])
      end if
      RHOG=0.0_WP; PG=0.0_WP; TG=0.0_WP; Yv=0.0_WP
      if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
         iVQ=7+this%liq%ns+this%indV-1
         ! Yv=Q(iVQ)/Q(2)
         Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
         y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         RHOG=Q(2)/(1.0_WP-VF)
         PG=this%gas%get_p_from_rho_e(rho=RHOG,e=Q(4)/Q(2),y=y)
         TG=this%gas%get_T_from_rho_e(rho=RHOG,e=Q(4)/Q(2),y=y)
      end if
      print*,trim(label)
      print*,'   i=',dbg_i,' j=',dbg_j
      if (present(ier)) print*,'   ier=',ier
      print*,'   VF=',VF,' Pjump=',Pjump
      print*,'   Q=',Q
      print*,'   RHOL=',RHOL,' PL=',PL,' TL=',TL
      print*,'   RHOG=',RHOG,' PG=',PG,' TG=',TG,' Yv=',Yv
   end subroutine debug_dump

   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: Yv,cvG,cpG,qG,gammaG,PG,PL,Ptar,Eold,TL,TG,RHOG
      integer  :: iVQ,ier
      logical :: run_pTg,interfacial
      interfacial=((VF.ge.VFlo).and.(VF.le.VFhi))
      if (dbg_cell) call this%debug_dump('ENTRY apply',VF,Q,Pjump)

      if (this%model.ne.PTgrelax) then
         ! Only mixture cells; the solver's pure-cell snap owns the rest
         if (.not.interfacial) then
            if (present(ierr)) ierr=RELAX_DEGENERATE
            return
         end if
      end if
      iVQ=7+this%liq%ns+this%indV-1
      ! Absorb (seed culling, high-pressure extreme): a gas packet above diss_P is
      ! supercritical and mixes into the liquid; conserves cell totals exactly, the cell
      ! becomes pure liquid, and the caller's pure-cell snap completes the PLIC reset
      if ((this%diss_P.lt.1.0e30_WP).and.(Q(2).gt.0.0_WP).and.(Q(4).gt.0.0_WP)) then
         ! Yv=Q(iVQ)/Q(2)
         Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
         if (PG.gt.this%diss_P) then
            this%acc(1)=this%acc(1)+1.0_WP; this%acc(2)=this%acc(2)+Q(2)*this%vol
            Q(1)=Q(1)+Q(2); Q(3)=Q(3)+Q(4)
            Q(2)=0.0_WP; Q(4)=0.0_WP; Q(iVQ)=0.0_WP
            VF=1.0_WP
            if (dbg_cell) call this%debug_dump('ABSORB',VF,Q,Pjump)
            if (present(ierr)) ierr=RELAX_OK
            return
         end if
      end if
      if ((this%model.eq.PTgrelax).and.(this%diss_RHO.lt.1.0e30_WP).and.(Q(2).gt.0.0_WP).and.(Q(4).gt.0.0_WP)) then
         RHOG=Q(2)/(1.0_WP-VF)
         if (RHOG.gt.this%diss_RHO) then
            this%acc(1)=this%acc(1)+1.0_WP; this%acc(2)=this%acc(2)+Q(2)*this%vol
            Q(1)=Q(1)+Q(2); Q(3)=Q(3)+Q(4)
            Q(2)=0.0_WP; Q(4)=0.0_WP; Q(iVQ)=0.0_WP
            VF=1.0_WP
            if (dbg_cell) call this%debug_dump('ABSORB-RHOG',VF,Q,Pjump)
            if (present(ierr)) ierr=RELAX_OK
            return
         end if
      end if
      ! Stage 1 — propose VF: stock quadratic relaxation (Prelax/PTrelax/PThybrid dispatch).
      ! Its return code only feeds the ledger; the energy split is set unconditionally below.
      ier=RELAX_OK
      select case (this%model)
      case (Prelax);   call this%p_relax  (dt,VF,Q,Pjump,ier)
      case (PTrelax);  call this%pT_relax (dt,VF,Q,Pjump,ier)
      case (PTgrelax)
         if (any(Q(1:4).lt.0.0_WP)) then
            if (dbg_cell) print*, 'Running p instead of pTg; Q=',Q
            call this%p_relax(dt,VF,Q,Pjump,ier)
         else
            run_pTg=.true.
            if (Q(1).gt.0.0_WP.and.Q(3).gt.0.0_WP) then
               if (this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP]).le.0.0_WP) run_pTg=.false.
            end if
            if (Q(2).gt.0.0_WP.and.Q(4).gt.0.0_WP) then
               ! Yv=Q(iVQ)/Q(2)
               Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
               y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
               if (this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y).le.0.0_WP) run_pTg=.false.
            end if
            if (run_pTg) then
               if (dbg_cell) print*, 'Running pTg'
               call this%pTg_relax(dt,VF,Q,Pjump,ier)
               if (dbg_cell) call this%debug_dump('POST-pTg_relax',VF,Q,Pjump,ier)
               if (ier.ne.RELAX_OK.and.interfacial) then
                  if (dbg_cell) print*, 'pTg relaxation was not ok'
                  if (ier.eq.RELAX_NUC_FAILED) then
                     if (dbg_cell) print*, 'Nucleation failure'
                     if (all(Q(1:4).gt.0.0_WP)) then
                        if (dbg_cell) print*, 'all Q(1:4) are positive'
                        ! Yv=Q(iVQ)/Q(2)
                        Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
                        y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
                        TL=this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
                        TG=this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
                        if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' NUC_FAILED-fallback TL=',TL,' TG=',TG ! debug
                        if (dbg_cell) print*,'   ratio=',merge(max(TG/TL,TL/TG),-1.0_WP,TL.gt.0.0_WP.and.TG.gt.0.0_WP),' Tratmax=',this%Tratmax ! debug
                        if (dbg_cell) print*,'   choice=',merge('pT_relax','p_relax ',TL.gt.0.0_WP.and.TG.gt.0.0_WP.and.max(TG/TL,TL/TG).gt.this%Tratmax) ! debug
                        ! &      ' choice=',merge('pT_relax','p_relax ',TL.gt.0.0_WP.and.TG.gt.0.0_WP) ! debug
                        ! if (TL.gt.0.0_WP.and.TG.gt.0.0_WP.and.max(TG/TL,TL/TG).gt.this%Tratmax) then
                        if (TL.gt.0.0_WP.and.TG.gt.0.0_WP) then
                           if (dbg_cell) print*, 'Temperatures positive; running pT'
                           call this%pT_relax(dt,VF,Q,Pjump,ier)
                           if (dbg_cell) call this%debug_dump('POST-pT_relax',VF,Q,Pjump,ier)
                        else
                           if (dbg_cell) print*, 'Negative T; running p.', '   TL=',TL,', TG=',TG
                           call this%p_relax(dt,VF,Q,Pjump,ier)
                           if (dbg_cell) call this%debug_dump('POST-p_relax',VF,Q,Pjump,ier)
                        end if
                     else
                        if (dbg_cell) print*, 'Some Q(1:4) are NOT positive, Q=',Q
                        call this%p_relax(dt,VF,Q,Pjump,ier)
                        if (dbg_cell) call this%debug_dump('POST-p_relax(negQ)',VF,Q,Pjump,ier)
                     end if
                  end if
                  ! Non-NUC_FAILED codes (pT_relax's own failure, or chem declined/failed after
                  ! pT_relax already succeeded): current VF/Q is the accepted answer already --
                  ! either pT_relax's own p_relax-consistent result, or the preserved
                  ! calling p_relax again here would be redundant.
                  if (dbg_cell) call this%debug_dump('POST-fallback to pT_relax',VF,Q,Pjump,ier)
               end if
            else
               if (dbg_cell) print*, 'Cant run pTg'
               ! if (interfacial) call this%p_relax(dt,VF,Q,Pjump,ier)
               if (interfacial) then
                  if (dbg_cell) print*, 'Running p'
                  call this%p_relax(dt,VF,Q,Pjump,ier)
                  if (dbg_cell) call this%debug_dump('POST-p_relax',VF,Q,Pjump,ier)
               end if
            end if
         end if
      case (PThybrid)
         if (any(Q(1:4).le.0.0_WP)) then
            call this%p_relax(dt,VF,Q,Pjump,ier)
         else
            ! Yv=Q(iVQ)/Q(2)
            Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
            y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
            TL=this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP])
            TG=this%gas%get_T_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
            if (TL.gt.0.0_WP.and.TG.gt.0.0_WP.and.max(TG/TL,TL/TG).gt.this%Tratmax) then
               call this%pT_relax(dt,VF,Q,Pjump,ier)
            else
               call this%p_relax(dt,VF,Q,Pjump,ier)
            end if
         end if
      case default; call die('[relax_igmix_sg apply] unknown model')
      end select
      ! Ledger the proposal outcome; a non-positive phase mass is the only untouched exit
      if (ier.eq.RELAX_OK) then
         this%acc(3)=this%acc(3)+1.0_WP
      else if (Q(1).le.0.0_WP.or.Q(2).le.0.0_WP) then
         this%acc(7)=this%acc(7)+1.0_WP
         ! if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' EXIT stuck (acc7) VF=',VF,' Q=',Q,&
         ! &                    ' RHOL=',Q(1)/VF,' PL=',this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP]),&
         ! &                    ' TL=',this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP]) ! debug
         if (present(ierr)) ierr=ier
         return
      else
         this%acc(4)=this%acc(4)+1.0_WP
      end if
      ! The proposal may exit at the VF bounds; the solver's pure-cell snap owns those
      if (VF.lt.VFlo.or.VF.gt.VFhi) then
         ! if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' EXIT at VF bounds VF=',VF,' Q=',Q,&
         ! &                    ' RHOL=',Q(1)/VF,' PL=',this%liq%get_p_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP]),&
         ! &                    ' TL=',this%liq%get_T_from_rho_e(rho=Q(1)/VF,e=Q(3)/Q(1),y=[1.0_WP]) ! debug
         if (present(ierr)) ierr=RELAX_OK
         return
      end if
      ! Gas composition
      ! Yv=Q(iVQ)/Q(2)
      Yv=max(0.0_WP,min(Q(iVQ)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave Q(iVQ)>Q(2))
      y=0.0_WP; y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Stage 2 — set energies (unconditional fixed-VF swap): whatever VF stage 1 produced
      ! (converged, VFratmax-clamped, or unchanged on failure), impose the unique energy
      ! split with PL-PG=Pjump at frozen VF and masses. Conserves phasic masses, total
      ! energy, momentum; a no-op to roundoff where the quadratic fully converged; completes
      ! the equilibration where it was clamped or failed (also the vacuum-runaway cutoff).
      ! Co-volume factor VF*(1-b*rhoL) clamped consistently with the nasg accessors.
      Eold=Q(3)+Q(4)
      PL=this%get_p_eq(VF,Q,qG,gammaG,Pjump)
      ! Stage 3 — floor: if the shared pressure sits below any user limit (PL,TL,PG,TG),
      ! raise it to the binding limit (all four rise monotonically with the shared pressure);
      ! the energy added is Asum*(Ptar-PL). Ledgered.
      Ptar=-1.0e30_WP
      if (this%Pmin_liq.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_liq)
      if (this%Pmin_gas.gt.-1.0e29_WP) Ptar=max(Ptar,this%Pmin_gas+Pjump)
      if (this%Tmin_liq.gt.0.0_WP) Ptar=max(Ptar,this%liq%get_p_from_rho_T(rho=Q(1)/VF,T=this%Tmin_liq,y=[1.0_WP]))
      if (this%Tmin_gas.gt.0.0_WP) Ptar=max(Ptar,this%gas%get_p_from_rho_T(rho=Q(2)/(1.0_WP-VF),T=this%Tmin_gas,y=y)+Pjump)
      if (PL.lt.Ptar) then
         this%acc(5)=this%acc(5)+1.0_WP
         this%acc(6)=this%acc(6)+(VF*this%liq%get_rhoe_from_p_rho(p=Ptar,rho=Q(1)/VF,y=[1.0_WP])+(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=Ptar-Pjump,rho=Q(2)/(1.0_WP-VF),y=y)-Eold)*this%vol
         PL=Ptar
      end if
      ! Write the split
      Q(3)=VF*this%liq%get_rhoe_from_p_rho(p=PL,rho=Q(1)/VF,y=[1.0_WP])
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=PL-Pjump,rho=Q(2)/(1.0_WP-VF),y=y)
      ! Final verdict: Prelax/PTrelax/PThybrid are always left in a consistent state by the
      ! swap above; PTgrelax's own convergence verdict (ier) is propagated unchanged.
      if (present(ierr)) then
         if (this%model.eq.PTgrelax) then
            ierr=ier
         else
            ierr=RELAX_OK
         end if
      end if
      if (dbg_cell) call this%debug_dump('EXIT-final',VF,Q,Pjump,ier)
   end subroutine apply

   !> Mechanical relaxation (Pelanti quadratic). Has clipping for unphysical phasic pressures.
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(:), allocatable   :: y
      real(WP) :: PL,PG,ZL,ZG,Pint
      real(WP) :: a,b,d,coeffL,coeffG,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Vapor mass fraction (from Q layout)
      if (Q(2).gt.0.0_WP) then
         ! Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
         Yv=max(0.0_WP,min(Q(7+this%liq%ns+this%indV-1)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave species mass>Q(2))
      else
         Yv=0.0_WP
      end if
      ! Gas composition vector (vapor + carrier-by-closure)
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      ! Inline mixture coefficients
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Phasic pressures
      PL=this%liq%get_p_from_rho_e(rho=Q(1)/(       VF),e=Q(3)/Q(1),y=[1.0_WP])
      PG=this%gas%get_p_from_rho_e(rho=Q(2)/(1.0_WP-VF),e=Q(4)/Q(2),y=y)
      ! Hard clipping for unphysical phasic pressure (cavitation, collapse)
      if (PL.le.-this%liq%pinf) then
         print*,"*** LIQUID CLIPPED!",PL,VF,Q
         VF=0.0_WP
         ! All flashed liquid mass becomes vapor (matches pTg_relax's pure-gas branch
         ! pattern Q(8)=Q(1)+Q(8)), not the old vapor fraction spread over the new total
         ! mass -- old code zeroed Q(1) before using it below, so it couldn't do that right
         ! Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         ! Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         ! Q(7+this%liq%ns+this%indV-1)=Yv*Q(2)
         Q(7+this%liq%ns+this%indV-1)=Q(1)+Q(7+this%liq%ns+this%indV-1)
         Q(2)=sum(Q(1:2)); Q(1)=0.0_WP
         Q(4)=sum(Q(3:4)); Q(3)=0.0_WP
         deallocate(y); if (present(ierr)) ierr=RELAX_BAD_LIQUID; return
      end if
      if (PG.le.0.0_WP) then
         print*,"*** GAS CLIPPED!",PG,VF,Q
         VF=1.0_WP
         Q(1)=sum(Q(1:2)); Q(2)=0.0_WP
         Q(3)=sum(Q(3:4)); Q(4)=0.0_WP
         Q(7+this%liq%ns+this%indV-1)=0.0_WP
         deallocate(y); if (present(ierr)) ierr=RELAX_BAD_GAS; return
      end if
      ! Phasic acoustic impedances (rho*c)
      ZL=Q(1)/(       VF)*this%liq%get_c_from_p_rho(p=PL,rho=Q(1)/(       VF),y=[1.0_WP])
      ZG=Q(2)/(1.0_WP-VF)*this%gas%get_c_from_p_rho(p=PG,rho=Q(2)/(1.0_WP-VF),y=y)
      ! Interface pressure
      Pint=(ZG*PL+ZL*PG)/(ZG+ZL)
      ! Quadratic for Peq
      coeffL=(this%liq%gamma-1.0_WP)*Pint+2.0_WP*this%liq%gamma*this%liq%pinf
      coeffG=(gammaG       -1.0_WP)*Pint
      a=1.0_WP+gammaG*VF+this%liq%gamma*(1.0_WP-VF)
      b=coeffL*(1.0_WP-VF)+coeffG*VF-(1.0_WP+gammaG)*VF*PL-(1.0_WP+this%liq%gamma)*(1.0_WP-VF)*PG
      d=-(coeffG*VF*PL+coeffL*(1.0_WP-VF)*PG)
      Peq =(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      VFeq=VF*((this%liq%gamma-1.0_WP)*Peq+2.0_WP*PL+coeffL)/((1.0_WP+this%liq%gamma)*Peq+coeffL)
      ! Clamp per-call phase-volume change to a factor VFratmax (bounds consistent for VFratmax>=1)
      VFeq=max(1.0_WP-this%VFratmax*(1.0_WP-VF),VF/this%VFratmax, &
      &        min(VFeq,1.0_WP-(1.0_WP-VF)/this%VFratmax,this%VFratmax*VF))
      ! Update conservatives
      Q(3)=Q(3)-0.5_WP*(Pint+Peq)*(VFeq-VF)
      Q(4)=Q(4)+0.5_WP*(Pint+Peq)*(VFeq-VF)
      VF=VFeq
      deallocate(y)
      if (present(ierr)) ierr=RELAX_OK
   end subroutine p_relax

   !> Mechanical + thermal relaxation. Calls p_relax first, then enforces TL=TG.
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(:), allocatable :: y
      real(WP) :: a,b,d,Peq,VFeq
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP) :: Yv
      ! Step 1: mechanical
      call this%p_relax(dt,VF,Q,Pjump)
      ! Step 2: thermal
      if (Q(2).gt.0.0_WP) then
         ! Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
         Yv=max(0.0_WP,min(Q(7+this%liq%ns+this%indV-1)/Q(2),1.0_WP)) ! debug: clamp as get_primitive does (near-total VF collapse can leave species mass>Q(2))
      else
         Yv=0.0_WP
      end if
      allocate(y(this%gas%ns))
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      ! Quadratic for the equilibrium pressure under PL=PG, TL=TG (SG liquid + IG mixture gas)
      a=Q(1)*this%liq%cv+Q(2)*cvG
      b=this%liq%q*this%liq%cv*(this%liq%gamma-1.0_WP)*Q(1)**2+qG*cvG*(gammaG-1.0_WP)*Q(2)**2+&
      &  Q(1)*this%liq%cv*this%liq%gamma*this%liq%pinf+Q(2)*cvG*this%liq%pinf                +&
      &  Q(1)*Q(2)*(this%liq%q*cvG*(gammaG-1.0_WP)+qG*this%liq%cv*(this%liq%gamma-1.0_WP))   -&
      &  sum(Q(3:4))*(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)+Q(2)*cvG*(gammaG-1.0_WP))
      d=cvG*(gammaG-1.0_WP)*this%liq%pinf*(qG*Q(2)**2+this%liq%q*Q(1)*Q(2)-sum(Q(3:4))*Q(2))
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      if (Peq.le.max(0.0_WP,-this%liq%pinf)) then
         deallocate(y); if (present(ierr)) ierr=RELAX_BAD_LIQUID; return
      end if
      VFeq=Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq &
      &   /(Q(1)*this%liq%cv*(this%liq%gamma-1.0_WP)*Peq+Q(2)*cvG*(gammaG-1.0_WP)*(Peq+this%liq%pinf))
      ! Clamp
      if (VFeq.lt.0.0_WP) then; VFeq=0.0_WP; Peq=max(Peq,-this%liq%pinf); end if
      if (VFeq.gt.1.0_WP) then; VFeq=1.0_WP; Peq=max(Peq,0.0_WP);          end if
      ! Update conservatives
      Q(3)=(       VFeq)*this%liq%get_rhoe_from_p_rho(p=Peq,rho=Q(1)/max(VFeq,tiny(1.0_WP)),y=[1.0_WP])
      Q(4)=(1.0_WP-VFeq)*this%gas%get_rhoe_from_p_rho(p=Peq,rho=Q(2)/max(1.0_WP-VFeq,tiny(1.0_WP)),y=y)
      VF=VFeq
      deallocate(y)
      if (present(ierr)) ierr=RELAX_OK
   end subroutine pT_relax

   !> Mechanical + thermal + chemical (phase change) relaxation.
   !> Nucleates a tiny opposite phase in metastable pure cells; calls pT_relax;
   !> runs pure-phase admissibility tests; falls back to LV (pure water) or LVG (with non-condensable) Newton solves.
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      use amrvof_class, only: VFlo,VFhi
      implicit none
      class(relax_igmix_sg),  intent(inout) :: this
      real(WP),               intent(in)    :: dt
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP),               intent(in)    :: Pjump
      integer,  optional,     intent(out)   :: ierr
      real(WP), dimension(:), allocatable   :: Qin,Q0,y
      real(WP) :: VFin,VF0,p,T,Yv,Yvin,Yv0,xv_eq
      real(WP) :: rho0,rhoe0,rhoA0
      real(WP) :: rhoL,rhoG
      real(WP) :: cvG,cpG,qG,gammaG
      real(WP), parameter :: p_eps=1.0e-10_WP,Yvmin=0.0_WP,Yvmax=1.0_WP
      real(WP), parameter :: Y_small=0.001_WP,Yv_pure=0.999_WP,Y_seed=0.01_WP
      real(WP), parameter :: fd_eps=1.0e-7_WP,F_line_search_tol=0.3_WP
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
      ! Old 0.5-split version, replaced by a narrow-band version below (misfired on pooled ~0.5 VF cells):
      ! pure_phase_stability: block
      !    real(WP) :: TL,pL,eL
      !    real(WP) :: TG,pG,eG
      !    real(WP) :: pV,Tsat,pVsat
      !    logical  :: conv
      !    integer  :: it
      !    if (VF.lt.0.5_WP) then
      !       if (dbg_cell) print*,'attempting pure gas'
      !       rhoG=Q(1)+Q(2)
      !       eG=(Q(3)+Q(4))/rhoG
      !       Yv=(Q(1)+Q(8))/rhoG
      !       y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      !       pG=this%gas%get_p_from_rho_e(rho=rhoG,e=eG,y=y)
      !       TG=this%gas%get_T_from_rho_e(rho=rhoG,e=eG,y=y)
      !       pV=this%get_xv(Yv)*pG
      !       call this%get_Tsat(pG,pV,TG,Tsat,conv,it)
      !       if (conv.and.TG.ge.Tsat) then
      !          VF=0.0_WP
      !          Q(2)=Q(1)+Q(2)
      !          Q(8)=Q(1)+Q(8)
      !          Q(1)=0.0_WP
      !          Q(4)=Q(3)+Q(4)
      !          Q(3)=0.0_WP
      !          if (present(ierr)) ierr=RELAX_OK
      !          return
      !       end if
      !    else
      !       if (Q(8).eq.Q(2)) then
      !          rhoL=Q(1)+Q(2)
      !          eL=(Q(3)+Q(4))/rhoL
      !          pL=this%liq%get_p_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
      !          TL=this%liq%get_T_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
      !          pVsat=this%get_pvsat(pL,TL)
      !          if (pL.ge.pVsat) then
      !             VF=1.0_WP
      !             Q(1)=Q(1)+Q(2)
      !             Q(2)=0.0_WP
      !             Q(8)=0.0_WP
      !             Q(3)=Q(3)+Q(4)
      !             Q(4)=0.0_WP
      !             if (present(ierr)) ierr=RELAX_OK
      !             return
      !          end if
      !       end if
      !    end if
      ! end block pure_phase_stability
      pure_phase_stability: block
         real(WP) :: TL,pL,eL
         real(WP) :: TG,pG,eG
         real(WP) :: pV,Tsat,pVsat
         real(WP), parameter :: pure_VFlo=0.05_WP,pure_VFhi=0.95_WP
         logical  :: conv
         integer  :: it
         if (VF.lt.pure_VFlo) then
            if (dbg_cell) print*,'attempting pure gas'
            ! Transfer all liquid mass and energy to vapor
            rhoG=Q(1)+Q(2)
            eG=(Q(3)+Q(4))/rhoG
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
               VF=0.0_WP
               Q(2)=Q(1)+Q(2)
               Q(8)=Q(1)+Q(8)
               Q(1)=0.0_WP
               Q(4)=Q(3)+Q(4)
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
               eL=(Q(3)+Q(4))/rhoL
               pL=this%liq%get_p_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
               TL=this%liq%get_T_from_rho_e(rho=rhoL,e=eL,y=[1.0_WP])
               ! Liquid compressed above its saturation pressure is stable -- no vapor survives.
               ! Testing pL against get_pvsat(pL,TL) puts p on both sides of the comparison, and
               ! pvsat grows exponentially in p, so that test stops being monotone in pressure and
               ! calls strongly compressed liquid superheated (see get_psat).
               call this%get_psat(TL,pVsat,conv)
               if (conv.and.pL.ge.pVsat) then
                  VF=1.0_WP
                  Q(1)=Q(1)+Q(2)
                  Q(2)=0.0_WP
                  Q(8)=0.0_WP
                  Q(3)=Q(3)+Q(4)
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
      if (this%do_nucleate.and.near_pure) then
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
      ! Step 3: chemical (phase change)
      if (Q(2).gt.0.0_WP) then
         ! Yv=Q(7+this%liq%ns+this%indV-1)/Q(2)
         Yv=max(Yvmin,min(Yvmax,Q(7+this%liq%ns+this%indV-1)/Q(2))) ! debug: clamp (near-total VF collapse can leave species mass>Q(2))
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
      rhoe0=sum(Q0(3:4))
      rhoA0=(1.0_WP-Yv)*Q0(2)
      if (Q0(2).gt.0.0_WP) then
         Yv0=Q0(7+this%liq%ns+this%indV-1)/Q0(2)
         Yv0=max(Yvmin,min(Yvmax,Yv0))
      else
         Yv0=0.0_WP
      end if
      ! Better Newton starting guess for liquid-dominated interfacial cells: Raoult/Dalton Yv from pVsat/p.
      ! Pjump deliberately omitted here (unlike elsewhere) -- its curvature noise is the same order as
      ! the pL-pVsat margin in this branch and checkers the Yv field; -Pjump tested, confirmed worse.
      ! if ((.not.nucleated).and.(VF.gt.0.5_WP).and.(Q(2).gt.0.0_WP).and.(Q(2).ne.Q(7+this%liq%ns+this%indV-1))) then
      if ((.not.nucleated).and.(VF.gt.0.5_WP).and.(Q(2).gt.0.0_WP).and.(Yv.lt.Y_small)) then
         xv_eq=this%get_pvsat(p,T)/p
         Yv=xv_eq*Mv/(xv_eq*Mv+(1.0_WP-xv_eq)*Ma)
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
      ! Solve chemical equilibrium for p, T, Yv (without modifying Q yet)
      if (Yv.gt.Yv_pure) then
         Yv=Yvmax
         y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
         cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
         cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
         qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
         gammaG=cpG/cvG
         call solve_lv(p,T,chem_relax)
      else
         call solve_lvg(p,T,Yv,chem_relax)
      end if
      ! Skip if not converged
      if (.not.chem_relax) then
         if (nucleated) then
            call undo_nucleation()
            if (present(ierr)) ierr=RELAX_NUC_FAILED
         else
            call restore_pT_state()
            if (present(ierr)) ierr=RELAX_FAILED
         end if
         call dealloc(); return
      end if
      ! Update densities and VF
      y(this%indV)=Yv; y(this%indA)=1.0_WP-Yv
      cvG   =sum(y(1:this%gas%ns)*this%gas%cv(1:this%gas%ns))
      cpG   =sum(y(1:this%gas%ns)*this%gas%cp(1:this%gas%ns))
      qG    =sum(y(1:this%gas%ns)*this%gas%q (1:this%gas%ns))
      gammaG=cpG/cvG
      RHOL=this%liq%get_rho_from_p_T(p=p,T=T,y=[1.0_WP])
      RHOG=this%gas%get_rho_from_p_T(p=p-Pjump,T=T,y=y)
      VF=(rho0-RHOG)/(RHOL-RHOG)
      if (dbg_cell) print*,'i=',dbg_i,' j=',dbg_j,' pTg converged: RHOL=',RHOL,', RHOG=',RHOG
      if (dbg_cell) print*,'   VF=',VF,', p=',p,', T=',T
      if (VF.lt.0.0_WP) then; VF=0.0_WP; p=max(p,-this%liq%pinf); end if
      if (VF.gt.1.0_WP) then; VF=1.0_WP; p=max(p,0.0_WP);         end if
      ! Confirm we had enough of the reducing phase to begin with
      if (.not.nucleated) then
         if (VF.lt.VF0) then
            ! Liquid turned into vapor while initial liquid content was smaller than the threshold
            if (Q0(1)/rho0.lt.Y_small) then
               call restore_pT_state(); call dealloc()
               if (present(ierr)) ierr=RELAX_FAILED
               return
            end if
         else if (VF.gt.VF0) then
            ! Vapor turned into liquid while initial vapor content was smaller than the threshold
            if (Yv0.lt.Y_small) then
               call restore_pT_state(); call dealloc()
               if (present(ierr)) ierr=RELAX_FAILED
               return
            end if
         end if
      end if
      ! Update Q with the converged equilibrium state
      Q(1)=(       VF)*RHOL
      Q(2)=(1.0_WP-VF)*RHOG
      Q(3)=Q(1)*this%liq%get_e_from_p_T(p=p,T=T,y=[1.0_WP])
      Q(4)=Q(2)*this%gas%get_e_from_p_T(p=p-Pjump,T=T,y=y)
      Q(7+this%liq%ns+this%indV-1)=Q(2)*Yv
      ! Evaluate conservation
      if (.not.check_cons()) then
         if (nucleated) then
            call undo_nucleation()
            if (present(ierr)) ierr=RELAX_NUC_FAILED
         else
            call restore_pT_state()
            if (present(ierr)) ierr=RELAX_FAILED
         end if
         call dealloc(); return
      end if
      call dealloc()
      if (present(ierr)) ierr=ier
   contains
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
      logical function check_cons()
         real(WP) :: re,ee
         re=(sum(Q(1:2))-rho0)/rho0
         ee=(sum(Q(3:4))-rhoe0)/rhoe0
         check_cons=(abs(re).le.this%rho_tol).and.(abs(ee).le.this%rhoe_tol)
      end function check_cons
      logical function activate_chem(p_,T_,Yv_)
         real(WP), intent(in)    :: p_,T_
         real(WP), intent(inout) :: Yv_
         real(WP) :: xv,pv_,Fsat
         if (dbg_cell) print*,'--------------------------------------------------'
         if (dbg_cell) print*,'inside activate_chem'
         ! Default to false
         activate_chem=.false.
         ! Get vapor mole fraction and partial pressure (gas-side pressure is p_-Pjump)
         xv=this%get_xv(Yv_); pv_=xv*(p_-Pjump)
         if (pv_.lt.this%pv_min) then
            if (dbg_cell) print*,'pv too small'
            pv_=exp(this%AS+(this%BS+this%ES*p_)/T_)*T_**this%CS*(p_+this%liq%pinf)**this%DS
            if (pv_.lt.this%pv_min) then
               ier=RELAX_VACUUM_VAPOR
               if (dbg_cell) print*,'RELAX_VACUUM_VAPOR. failure'
               return
            else if (pv_.ge.(p_-Pjump)) then
               if (dbg_cell) print*,'Too much vapor pressure, seeding Yv=',Y_seed
               Yv_=Y_seed
            else
               xv=pv_/(p_-Pjump); Yv_=xv*Mv/(xv*Mv+(1.0_WP-xv)*Ma)
               if (Yv_.lt.Y_small) then
                  if (dbg_cell) print*,'Too little Yv, seeding Yv=',Y_seed
                  Yv_=Y_seed
               end if
            end if
            Yv_=max(Yvmin,min(Yvmax,Yv_))
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
      real(WP) function rhoe_res_lvg(p_,T_,Yv_)
         real(WP), intent(in) :: p_,T_,Yv_
         real(WP), dimension(this%gas%ns) :: y_
         ! debug: composition must be indexed by this%indV/this%indA -- species order is not
         ! guaranteed to be [vapor,air], and this%gas%ns may exceed 2. Old hardcoded literal:
         ! rhoe_res_lvg=(rho0*(1.0_WP-Yv_)-rhoA0)*this%liq%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP])+&
         ! &             rhoA0*this%gas%get_e_from_p_T(p=p_,T=T_,y=[Yv_,1.0_WP-Yv_])           -&
         ! &             rhoe0*(1.0_WP-Yv_)
         y_(this%indV)=Yv_; y_(this%indA)=1.0_WP-Yv_
         rhoe_res_lvg=(rho0*(1.0_WP-Yv_)-rhoA0)*this%liq%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP])+&
         &             rhoA0*this%gas%get_e_from_p_T(p=p_,T=T_,y=y_)                          -&
         &             rhoe0*(1.0_WP-Yv_)
      end function rhoe_res_lvg
      real(WP) function dlnxvdYv(Yv_)
         real(WP), intent(in) :: Yv_
         real(WP) :: Ys,den
         Ys=max(Yv_,Yvmin+fd_eps)
         den=Ys*Ma+(1.0_WP-Ys)*Mv
         dlnxvdYv=1.0_WP/Ys-(Ma-Mv)/den
      end function dlnxvdYv
      !> LV solve: damped Newton in ln(p)
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
         conv=.false.; p_err=10.0_WP*this%p_tol
         do it=1,this%NR_itmax
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp); if (T_eq.le.0.0_WP) return
            dTdp=get_dTdp_lv(ap,bp,dp,dapdp,dbpdp,ddpdp); dTdlnp=p_eq*dTdp
            F1=this%pTsat(p_eq,p_eq-Pjump,T_eq); dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp,Pjump)
            if (abs(dF1dlnp).lt.1.0e-30_WP) exit
            dlnp_nr=-F1/dF1dlnp
            dlnp_nr=max(log(0.5_WP),min(log(1.5_WP),dlnp_nr))
            pOld=p_eq; lnpOld=log(pOld); alpha=1.0_WP; accepted=.false.
            do while (alpha.gt.1.0e-8_WP)
               p_try=exp(lnpOld+alpha*dlnp_nr)
               if (p_try.le.p_eps) then; alpha=0.5_WP*alpha; cycle; end if
               call this%get_coeffs_lv(p_try,rho0,rhoe0,cvG,gammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
               T_try=get_T_lv(ap,bp,dp)
               if (T_try.le.0.0_WP) then; alpha=0.5_WP*alpha; cycle; end if
               F1_try=this%pTsat(p_try,p_try-Pjump,T_try)
               if (abs(F1_try).lt.abs(F1)) then
                  p_eq=p_try; T_eq=T_try; accepted=.true.; exit
               end if
               alpha=0.5_WP*alpha
            end do
            if (.not.accepted) exit
            p_err=abs(log(p_eq/pOld))
            call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
            T_eq=get_T_lv(ap,bp,dp); if (T_eq.le.0.0_WP) return
            F1=this%pTsat(p_eq,p_eq-Pjump,T_eq)
            if ((p_err.lt.this%p_tol).and.(abs(F1).lt.this%F1_tol)) then
               conv=.true.; exit
            end if
         end do
         if (.not.conv) return
         call this%get_coeffs_lv(p_eq,rho0,rhoe0,cvG,gammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
         T_eq=get_T_lv(ap,bp,dp)
      end subroutine solve_lv
      !> LVG solve: 2x2 damped Newton in (ln(p), Yv) with step limiter
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
         integer  :: it,lsit,hi_lock
         logical  :: accepted,bnd_ok
         Yv_max_phys=1.0_WP-(rhoA0/rho0)
         Yv_hi=min(Yvmax,Yv_max_phys)
         if (dbg_cell) then
            print*,'--------------------------------------------------'
            print*,'inside solve_lvg, p_eq=',p_eq,' T_eq=',T_eq,' Yv_eq=',Yv_eq
            print*,'inside solve_lvg, rhoA0=',rhoA0,', rho0=',rho0,', Yv_hi=',Yv_hi
         end if
         conv=.false.
         p_err=10.0_WP*this%p_tol; Yv_err=10.0_WP*this%Yv_tol
         hi_lock=0
         do it=1,this%NR_itmax
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0,Pjump)
            if (T_eq.le.0.0_WP) then
               if (dbg_cell) print*,'EXIT-A T_eq<=0 it=',it,' p_eq=',p_eq,' Yv_eq=',Yv_eq ! debug
               return
            end if
            xv=this%get_xv(Yv_eq); pv=xv*(p_eq-Pjump)
            if (.not.check_pv(pv)) then
               if (dbg_cell) print*,'EXIT-B check_pv(pv) it=',it,' pv=',pv ! debug
               return
            end if
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            res0=sqrt(F1**2+F2**2)
            if (dbg_cell) print*,'ITER-START it=',it
            if (dbg_cell) print*,'   p_eq=',p_eq,' T_eq=',T_eq,' Yv_eq=',Yv_eq
            if (dbg_cell) print*,'   xv=',xv,' pv=',pv,' F1=',F1,' F2=',F2,' res0=',res0
            if (dbg_cell) print*,'   rho0=',rho0,' rhoA0=',rhoA0,' rhoe0=',rhoe0
            ! d/dlnp via forward difference
            lnp_pert=log(p_eq)+fd_eps; p_pert=exp(lnp_pert); pv_pert=xv*(p_pert-Pjump)
            if (.not.check_pv(pv_pert)) then
               if (dbg_cell) print*,'EXIT-C check_pv(pv_pert,lnp) it=',it,' pv_pert=',pv_pert ! debug
               return
            end if
            T_pert=this%get_T_lvg(p_pert,Yv_eq,rho0,rhoA0,Pjump)
            if (T_pert.le.0.0_WP) then
               if (dbg_cell) print*,'EXIT-D T_pert<=0(lnp) it=',it,' p_pert=',p_pert ! debug
               return
            end if
            F2p=rhoe_res_lvg(p_pert,T_pert,Yv_eq)/rhoe0
            dTdlnp=(T_pert-T_eq)/fd_eps
            dF1dlnp=this%dpTsatdlnp(p_eq,T_eq,dTdlnp,Pjump)
            dF2dlnp=(F2p-F2)/fd_eps
            ! d/dYv via forward (or backward) difference, staying inside domain
            Yv_pert=Yv_eq+fd_eps
            if (Yv_pert.gt.Yv_hi-fd_eps) Yv_pert=Yv_eq-fd_eps
            if (Yv_pert.lt.Yvmin+fd_eps) Yv_pert=Yv_eq+fd_eps
            if ((Yv_pert.le.Yvmin+fd_eps).or.(Yv_pert.ge.Yv_hi-fd_eps)) then
               if (dbg_cell) print*,'EXIT-E Yv_pert out of [Yvmin,Yv_hi] it=',it
               if (dbg_cell) print*,'   Yv_eq=',Yv_eq,' Yv_pert=',Yv_pert,' Yvmin=',Yvmin,' Yv_hi=',Yv_hi
               return
            end if
            xv_pert=this%get_xv(Yv_pert); pv_pert=xv_pert*(p_eq-Pjump)
            if (.not.check_pv(pv_pert)) then
               if (dbg_cell) print*,'EXIT-F check_pv(pv_pert,Yv) it=',it,' pv_pert=',pv_pert ! debug
               return
            end if
            T_pert=this%get_T_lvg(p_eq,Yv_pert,rho0,rhoA0,Pjump)
            if (T_pert.le.0.0_WP) then
               if (dbg_cell) print*,'EXIT-G T_pert<=0(Yv) it=',it,' Yv_pert=',Yv_pert ! debug
               return
            end if
            F2Y=rhoe_res_lvg(p_eq,T_pert,Yv_pert)/rhoe0
            dTdYv=(T_pert-T_eq)/(Yv_pert-Yv_eq)
            dF1dYv=this%dpTsatdT(p_eq,T_eq)*dTdYv-dlnxvdYv(Yv_eq)
            dF2dYv=(F2Y-F2)/(Yv_pert-Yv_eq)
            ! 2x2 Newton
            detJ=dF1dlnp*dF2dYv-dF1dYv*dF2dlnp
            if (dbg_cell) print*,'JACOBIAN it=',it
            if (dbg_cell) print*,'   dF1dlnp=',dF1dlnp,' dF1dYv=',dF1dYv
            if (dbg_cell) print*,'   dF2dlnp=',dF2dlnp,' dF2dYv=',dF2dYv,' detJ=',detJ
            if (abs(detJ).lt.1.0e-30_WP) then
               if (dbg_cell) print*,'EXIT-H detJ singular it=',it,' detJ=',detJ ! debug
               exit
            end if
            dlnp_nr=-( dF2dYv *F1-dF1dYv  *F2)/detJ
            dYv_nr =-(-dF2dlnp*F1+dF1dlnp*F2)/detJ
            if (dbg_cell) print*,'NR-STEP-RAW it=',it,' dlnp_nr=',dlnp_nr,' dYv_nr=',dYv_nr ! debug
            ! Direction-preserving step limiter
            step_limit: block
               real(WP) :: ms,lnp_up,lnp_dn
               logical  :: hit_hi
               ms=1.0_WP; hit_hi=.false.
               lnp_up=log(1.5_WP); lnp_dn=log(0.5_WP)
               if (dlnp_nr.gt.lnp_up) ms=min(ms,lnp_up/dlnp_nr)
               if (dlnp_nr.lt.lnp_dn) ms=min(ms,lnp_dn/dlnp_nr)
               if (dYv_nr.gt.0.0_WP) then
                  if (Yv_eq+dYv_nr.ge.Yv_hi-fd_eps) then
                     ms=min(ms,0.9_WP*(Yv_hi-fd_eps-Yv_eq)/dYv_nr); hit_hi=.true.
                  end if
               else if (dYv_nr.lt.0.0_WP) then
                  if (Yv_eq+dYv_nr.le.Yvmin+fd_eps) ms=min(ms,0.9_WP*(Yv_eq-Yvmin-fd_eps)/abs(dYv_nr))
               end if
               ms=max(0.0_WP,min(1.0_WP,ms))
               dlnp_nr=dlnp_nr*ms; dYv_nr=dYv_nr*ms
               if (hit_hi) then; hi_lock=hi_lock+1; else; hi_lock=0; end if
               if (dbg_cell) print*,'NR-STEP-LIMITED it=',it,' ms=',ms,' hi_lock=',hi_lock
               if (dbg_cell) print*,'   dlnp_nr=',dlnp_nr,' dYv_nr=',dYv_nr
            end block step_limit
            ! Persistently pinned against the evaporation ceiling (Yv_hi) for several iterations
            ! straight: the raw step keeps demanding more vapor than the cell's liquid can ever
            ! supply, so no interior two-phase equilibrium exists. The true equilibrium is then the
            ! boundary itself -- all liquid evaporated (VF=0) -- which has a direct, non-iterative
            ! solution (fixed mass rho0, fixed energy rhoe0, fixed composition Yv_hi, single gas
            ! phase, no saturation condition left to satisfy since no liquid remains).
            if (hi_lock.ge.3) then
               call boundary_all_vapor(Yv_hi,p_eq,T_eq,Yv_eq,bnd_ok)
               if (bnd_ok) then
                  if (dbg_cell) print*,'EXIT-K pinned at Yv_hi for',hi_lock,' iterations -- resolving as complete evaporation' ! debug
                  conv=.true.
                  return
               end if
               ! The cell cannot supply the latent heat, so complete evaporation is not the answer.
               ! The saturation residual is positive on BOTH sides of the two-phase region, so a
               ! Newton pinned against the vapor ceiling here is really sitting on compressed
               ! liquid: resolve at the opposite boundary instead.
               call boundary_all_liquid(p_eq,T_eq,Yv_eq,bnd_ok)
               if (bnd_ok) then
                  if (dbg_cell) print*,'EXIT-K pinned at Yv_hi but all-vapor infeasible -- resolving as complete condensation' ! debug
                  conv=.true.
                  return
               end if
               if (dbg_cell) print*,'EXIT-K pinned at Yv_hi, neither boundary feasible -- rejecting' ! debug
               conv=.false.
               return
            end if
            ! Damped update with line search
            pOld=p_eq; YvOld=Yv_eq; lnpOld=log(pOld); alpha=1.0_WP; lsit=0
            if ((abs(F1).lt.F_line_search_tol).and.(abs(F2).lt.F_line_search_tol)) then
               p_eq=exp(lnpOld+dlnp_nr); Yv_eq=YvOld+dYv_nr
               if (dbg_cell) print*,'FULL-STEP(no-line-search) it=',it
               if (dbg_cell) print*,'   p_eq=',p_eq,' Yv_eq=',Yv_eq
            else
               accepted=.false.
               do while (alpha.gt.1.0e-8_WP)
                  lsit=lsit+1
                  p_try=exp(lnpOld+alpha*dlnp_nr); Yv_try=YvOld+alpha*dYv_nr
                  if (p_try.le.p_eps) then; alpha=0.5_WP*alpha; cycle; end if
                  if ((Yv_try.le.Yvmin+fd_eps).or.(Yv_try.ge.Yv_hi-fd_eps)) then
                     alpha=0.5_WP*alpha; cycle
                  end if
                  if ((rho0*(1.0_WP-Yv_try)-rhoA0).le.0.0_WP) then
                     alpha=0.5_WP*alpha; cycle
                  end if
                  T_try=this%get_T_lvg(p_try,Yv_try,rho0,rhoA0,Pjump)
                  if (T_try.le.0.0_WP) then; alpha=0.5_WP*alpha; cycle; end if
                  xv_try=this%get_xv(Yv_try); pv_try=xv_try*(p_try-Pjump)
                  if (.not.check_pv(pv_try)) then; alpha=0.5_WP*alpha; cycle; end if
                  F1_try=this%pTsat(p_try,pv_try,T_try)
                  F2_try=rhoe_res_lvg(p_try,T_try,Yv_try)/rhoe0
                  res_try=sqrt(F1_try**2+F2_try**2)
                  if (dbg_cell) print*,'LINESEARCH it=',it,' lsit=',lsit,' alpha=',alpha
                  if (dbg_cell) print*,'   p_try=',p_try,' Yv_try=',Yv_try,' T_try=',T_try
                  if (dbg_cell) print*,'   F1_try=',F1_try,' F2_try=',F2_try,' res_try=',res_try,' res0=',res0
                  if (res_try.lt.res0) then
                     p_eq=p_try; Yv_eq=Yv_try; T_eq=T_try; accepted=.true.; exit
                  end if
                  alpha=0.5_WP*alpha
               end do
               if (.not.accepted) then
                  if (dbg_cell) print*,'EXIT-I line search not accepted it=',it
                  if (dbg_cell) print*,'   lsit=',lsit,' alpha=',alpha,' res0=',res0
                  exit
               end if
            end if
            p_err=abs(log(p_eq/pOld))
            Yv_err=abs(Yv_eq-YvOld)
            xv=this%get_xv(Yv_eq); pv=xv*(p_eq-Pjump)
            if (.not.check_pv(pv)) then
               if (dbg_cell) print*,'EXIT-J check_pv(pv,post-step) it=',it,' pv=',pv ! debug
               return
            end if
            T_eq=this%get_T_lvg(p_eq,Yv_eq,rho0,rhoA0,Pjump)
            if (T_eq.le.0.0_WP) then
               if (dbg_cell) print*,'EXIT-K T_eq<=0(post-step) it=',it,' p_eq=',p_eq,' Yv_eq=',Yv_eq ! debug
               return
            end if
            F1=this%pTsat(p_eq,pv,T_eq)
            F2=rhoe_res_lvg(p_eq,T_eq,Yv_eq)/rhoe0
            if (dbg_cell) print*,'ITER-END it=',it,' p_eq=',p_eq,' T_eq=',T_eq,' Yv_eq=',Yv_eq
            if (dbg_cell) print*,'   p_err=',p_err,' p_tol=',this%p_tol,' Yv_err=',Yv_err
            if (dbg_cell) print*,'   Yv_tol_eff=',this%Yv_tol_abs+this%Yv_tol*max(abs(YvOld),abs(Yv_eq))
            if (dbg_cell) print*,'   F1=',F1,' F1_tol=',this%F1_tol,' F2=',F2,' F2_tol=',this%F2_tol
            if ((p_err.lt.this%p_tol).and.Yv_err.lt.this%Yv_tol_abs+this%Yv_tol*max(abs(YvOld),abs(Yv_eq)).and.(abs(F1).lt.this%F1_tol).and.(abs(F2).lt.this%F2_tol)) then
               conv=.true.; exit
            end if
         end do
         if (dbg_cell.and..not.conv) print*,'EXIT-L exhausted NR_itmax=',this%NR_itmax,' final it=',it-1
         if (dbg_cell.and..not.conv) print*,'   p_eq=',p_eq,' T_eq=',T_eq,' Yv_eq=',Yv_eq,' p_err=',p_err,' Yv_err=',Yv_err
      end subroutine solve_lvg
      !> Direct (non-iterative) equilibrium for the evaporation ceiling: with Yv pinned at
      !> Yv_hi_, all liquid is gone (VF=0), so the cell is single-phase gas at fixed mass rho0,
      !> fixed energy rhoe0, and fixed composition Yv_hi_ -- there is no liquid left for a
      !> saturation condition to hold against, so p,T follow directly from the gas EOS, the same
      !> (rho,e)->p->T route apply() itself uses for a single-phase gas cell.
      subroutine boundary_all_vapor(Yv_hi_,p_eq,T_eq,Yv_eq,feasible)
         real(WP), intent(in)    :: Yv_hi_
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: feasible
         real(WP) :: p_gas,T_gas
         y(this%indV)=Yv_hi_; y(this%indA)=1.0_WP-Yv_hi_
         p_gas=this%gas%get_p_from_rho_e(rho=rho0,e=rhoe0/rho0,y=y)
         T_gas=this%gas%get_T_from_p_rho(p=p_gas,rho=rho0,y=y)
         ! Boiling off the last of the liquid costs the vapor's formation energy q. A cell whose
         ! internal energy sits below it cannot BE vapor: e-q is negative, so T and p come out
         ! negative and this boundary is not an equilibrium at all -- it must be refused, not
         ! reported as converged.
         feasible=(T_gas.gt.0.0_WP.and.p_gas.gt.0.0_WP)
         if (.not.feasible) return
         Yv_eq=Yv_hi_; T_eq=T_gas; p_eq=p_gas+Pjump
      end subroutine boundary_all_vapor
      !> Mirror of boundary_all_vapor at the condensation floor: all vapor condensed (VF=1), so the
      !> cell is single-phase liquid at fixed mass rho0 and fixed energy rhoe0, with no vapor left
      !> for a saturation condition to hold against. This is the equilibrium for a cell denser than
      !> saturated liquid -- compressed, subcooled liquid, where no two-phase split exists at any
      !> temperature and the vapor boundary is energetically out of reach.
      subroutine boundary_all_liquid(p_eq,T_eq,Yv_eq,feasible)
         real(WP), intent(inout) :: p_eq,T_eq,Yv_eq
         logical,  intent(out)   :: feasible
         real(WP) :: p_l,T_l,r1,r2,r1p,r2p,r1T,r2T,J11,J12,J21,J22,det,dp,dT,escl
         real(WP), parameter :: fd=1.0e-7_WP,rtol=1.0e-12_WP
         integer :: itb
         feasible=.false.
         if (rho0-rhoA0.le.0.0_WP) return
         ! Start from the air-free answer: with no air the cell is single-phase liquid and (p,T)
         ! follow directly from the liquid EOS at fixed mass and energy. Exact when rhoA0=0.
         T_l=this%liq%get_T_from_rho_e(rho=rho0,e=rhoe0/rho0,y=[1.0_WP])
         p_l=this%liq%get_p_from_rho_e(rho=rho0,e=rhoe0/rho0,y=[1.0_WP])
         if (T_l.le.0.0_WP.or.p_l+this%liq%pinf.le.0.0_WP) return
         if (rhoA0.le.0.0_WP) then
            feasible=.true.; Yv_eq=Yvmin; T_eq=T_l; p_eq=p_l; return
         end if
         ! Air is non-condensable: it cannot enter the liquid, so it stays behind as a pure-air gas
         ! phase and this boundary is NOT single-phase. Solve the volume closure and the energy
         ! balance for (p,T) with the liquid holding rho0-rhoA0 and the gas holding rhoA0. Folding
         ! the air into the liquid instead is exact only as rhoA0->0 (3% error in p at 1% air mass).
         escl=max(abs(rhoe0),1.0_WP)
         do itb=1,50
            call resid_all_liquid(p_l,T_l,r1,r2)
            if (max(abs(r1),abs(r2)/escl).lt.rtol) then; feasible=.true.; exit; end if
            call resid_all_liquid(p_l+fd*(abs(p_l)+this%liq%pinf),T_l,r1p,r2p)
            call resid_all_liquid(p_l,T_l*(1.0_WP+fd),r1T,r2T)
            J11=(r1p-r1)/(fd*(abs(p_l)+this%liq%pinf)); J21=(r2p-r2)/(fd*(abs(p_l)+this%liq%pinf))
            J12=(r1T-r1)/(fd*T_l);                      J22=(r2T-r2)/(fd*T_l)
            det=J11*J22-J12*J21
            if (abs(det).le.tiny(1.0_WP)) return
            dp=-( J22*r1-J12*r2)/det
            dT=-(-J21*r1+J11*r2)/det
            ! keep the step inside the liquid's pressure floor and positive temperature
            if (p_l+dp+this%liq%pinf.le.0.0_WP) dp=-0.5_WP*(p_l+this%liq%pinf)
            if (T_l+dT.le.0.0_WP)               dT=-0.5_WP*T_l
            p_l=p_l+dp; T_l=T_l+dT
         end do
         if (.not.feasible) return
         feasible=(T_l.gt.0.0_WP.and.p_l+this%liq%pinf.gt.0.0_WP)
         if (.not.feasible) return
         Yv_eq=Yvmin; T_eq=T_l; p_eq=p_l
      end subroutine boundary_all_liquid
      !> Volume-closure and energy residuals for boundary_all_liquid: liquid mass rho0-rhoA0 and
      !> air mass rhoA0 must together fill the cell and carry its energy at (p_,T_).
      subroutine resid_all_liquid(p_,T_,r1,r2)
         real(WP), intent(in)  :: p_,T_
         real(WP), intent(out) :: r1,r2
         real(WP) :: mL,rl,rg
         r1=0.0_WP; r2=0.0_WP
         mL=rho0-rhoA0
         rl=this%liq%get_rho_from_p_T(p=p_,T=T_,y=[1.0_WP])
         y(this%indV)=0.0_WP; y(this%indA)=1.0_WP
         rg=this%gas%get_rho_from_p_T(p=p_-Pjump,T=T_,y=y)
         if (rl.le.0.0_WP.or.rg.le.0.0_WP) then; r1=huge(1.0_WP); return; end if
         r1=mL/rl+rhoA0/rg-1.0_WP
         r2=mL*this%liq%get_e_from_p_T(p=p_,T=T_,y=[1.0_WP]) &
         & +rhoA0*this%gas%get_e_from_p_T(p=p_-Pjump,T=T_,y=y)-rhoe0
      end subroutine resid_all_liquid
   end subroutine pTg_relax

   !> p-T saturation residual (general form: ES=0 for SG, ES=b/RV for NASG)
   real(WP) function pTsat(this,pl_,pv_,T_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: pl_,pv_,T_
      pTsat=this%AS+(this%BS+this%ES*pl_)/T_+this%CS*log(T_)+this%DS*log(pl_+this%liq%pinf)-log(pv_)
   end function pTsat

   real(WP) function dpTsatdT(this,pl_,T_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: pl_,T_
      dpTsatdT=-(this%BS+this%ES*pl_)/T_**2+this%CS/T_
   end function dpTsatdT

   real(WP) function dpTsatdp_lv(this,p_,T_,dTdp_,Pjump_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdp_,Pjump_
      dpTsatdp_lv=this%dpTsatdT(p_,T_)*dTdp_+this%ES/T_+this%DS/(p_+this%liq%pinf)-1.0_WP/(p_-Pjump_)
   end function dpTsatdp_lv

   real(WP) function dpTsatdlnp(this,p_,T_,dTdlnp_,Pjump_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,T_,dTdlnp_,Pjump_
      dpTsatdlnp=this%dpTsatdT(p_,T_)*dTdlnp_+this%ES*p_/T_+this%DS*p_/(p_+this%liq%pinf)-p_/(p_-Pjump_)
   end function dpTsatdlnp

   !> Safeguarded Newton on the saturation curve for Tsat(pl, pv)
   subroutine get_Tsat(this,pl_,pv_,Tguess,Tsat,conv,Tsat_it)
      implicit none
      class(relax_igmix_sg), intent(inout) :: this
      real(WP), intent(in)  :: pl_,pv_,Tguess
      real(WP), intent(out) :: Tsat
      logical,  intent(out) :: conv
      integer,  intent(out) :: Tsat_it
      real(WP) :: Tlo,Thi,Told,Tnew,Flo,Fhi,Fold,Fnew,dFold,Tpeak
      integer  :: it,expand_it
      ! Tsat is intent(out), so it must be defined on the early-return paths below too -- leaving it
      ! undefined there returns whatever was on the stack (a denormal, in practice)
      conv=.false.; Tsat_it=0; Tsat=Tguess
      if (Tguess.le.0.0_WP) return
      ! pTsat is NOT monotonic: it rises from -infinity as T->0, peaks where dpTsatdT=0, then falls.
      ! Only below that peak does the fitted vapor pressure increase with T, so that branch is the
      ! physical one and carries the root we want. Bracketing off the two ends of a fixed window
      ! straddles the peak and reports "no root" even when two exist.
      Tpeak=-1.0_WP
      if (this%CS.ne.0.0_WP) Tpeak=(this%BS+this%ES*pl_)/this%CS
      if (Tpeak.gt.0.0_WP) then
         ! search the rising branch only, starting from the guess (or the peak, if the guess is past it)
         Thi=min(Tguess,Tpeak); Fhi=this%pTsat(pl_,pv_,Thi)
         if (Fhi.lt.0.0_WP.and.Thi.lt.Tpeak) then
            Thi=Tpeak; Fhi=this%pTsat(pl_,pv_,Thi)
         end if
         ! the peak is the largest pTsat can ever be: still negative there means no saturation state
         if (Fhi.lt.0.0_WP) return
         Tlo=Thi; Flo=Fhi; expand_it=0
         do while ((Flo.ge.0.0_WP).and.(expand_it.lt.60))
            Tlo=0.8_WP*Tlo;             Flo=this%pTsat(pl_,pv_,Tlo)
            expand_it=expand_it+1
         end do
         if (Flo.ge.0.0_WP) return
      else
         ! degenerate fit (CS=0): pTsat is monotonic, so the original two-sided expansion applies
         Tlo=0.3_WP*Tguess; Thi=3.0_WP*Tguess
         Flo=this%pTsat(pl_,pv_,Tlo)
         Fhi=this%pTsat(pl_,pv_,Thi)
         expand_it=0
         do while ((Flo*Fhi.gt.0.0_WP).and.(expand_it.lt.20))
            if ((Flo.gt.0.0_WP).and.(Fhi.gt.0.0_WP)) then
               Tlo=0.8_WP*Tlo;             Flo=this%pTsat(pl_,pv_,Tlo)
            else if ((Flo.lt.0.0_WP).and.(Fhi.lt.0.0_WP)) then
               Thi=1.2_WP*Thi;             Fhi=this%pTsat(pl_,pv_,Thi)
            else
               exit
            end if
            expand_it=expand_it+1
         end do
         if (Flo*Fhi.gt.0.0_WP) return
      end if
      Tsat=max(Tlo,min(Thi,Tguess))
      do it=1,this%Tsat_itmax
         Told=Tsat
         Fold=this%pTsat(pl_,pv_,Told); dFold=this%dpTsatdT(pl_,Told)
         if (abs(dFold).gt.tiny(1.0_WP)) then
            Tnew=Told-Fold/dFold
         else
            Tnew=0.5_WP*(Tlo+Thi)
         end if
         if ((Tnew.ne.Tnew).or.(Tnew.le.Tlo).or.(Tnew.ge.Thi)) Tnew=0.5_WP*(Tlo+Thi)
         Fnew=this%pTsat(pl_,pv_,Tnew)
         if (Fnew.ne.Fnew) then
            Tnew=0.5_WP*(Tlo+Thi); Fnew=this%pTsat(pl_,pv_,Tnew)
         end if
         if (Flo*Fnew.le.0.0_WP) then
            Thi=Tnew; Fhi=Fnew
         else
            Tlo=Tnew; Flo=Fnew
         end if
         Tsat_it=it; Tsat=Tnew
         if ((abs((Tnew-Told)/max(abs(Told),tiny(1.0_WP))).lt.this%Tsat_tol).or.(abs(Fnew).lt.this%F1_tol)) then
            conv=.true.; return
         end if
      end do
   end subroutine get_Tsat

   !> Physical saturation pressure at T_: the LOWER fixed point of pvsat(p,T)=p.
   !>
   !> The residual ln(pvsat(p,T))-ln(p) falls from +infinity as p->0, passes through a single
   !> minimum, then climbs back to +infinity, so it has TWO roots. Only the lower one is physical --
   !> the upper one is the covolume/Poynting terms (ES*p/T, DS*ln(p+pinf)) taking over far outside
   !> the range the Clausius-Clapeyron fit was calibrated on. Any stability test that compares pL
   !> against get_pvsat(pL,T) instead flips sign at that spurious root and reports strongly
   !> compressed liquid as boiling. The minimum is available in closed form (the residual's
   !> derivative is a quadratic in p), so the lower root can be bracketed on the monotone side.
   subroutine get_psat(this,T_,psat,conv)
      implicit none
      class(relax_igmix_sg), intent(inout) :: this
      real(WP), intent(in)  :: T_
      real(WP), intent(out) :: psat
      logical,  intent(out) :: conv
      real(WP) :: a,b,c,disc,pmin,plo,phi,pmid
      integer  :: it
      conv=.false.; psat=0.0_WP
      if (T_.le.0.0_WP) return
      ! d/dp of the residual is ES/T + DS/(p+pinf) - 1/p; clearing denominators gives a quadratic
      a=this%ES/T_
      b=this%ES*this%liq%pinf/T_+this%DS-1.0_WP
      c=-this%liq%pinf
      if (abs(a).gt.tiny(1.0_WP)) then
         disc=b*b-4.0_WP*a*c
         if (disc.lt.0.0_WP) return
         pmin=0.5_WP*(-b+sqrt(disc))/a
      else
         if (abs(b).le.tiny(1.0_WP)) return
         pmin=-c/b
      end if
      if (pmin.le.0.0_WP) return
      ! the minimum is the lowest the residual ever gets: still positive means no saturation state
      if (psat_res(pmin).ge.0.0_WP) return
      ! walk down from the minimum until the residual turns positive, then bisect
      plo=pmin
      do it=1,200
         plo=0.5_WP*plo
         if (psat_res(plo).ge.0.0_WP) exit
      end do
      if (psat_res(plo).lt.0.0_WP) return
      phi=pmin
      do it=1,200
         pmid=0.5_WP*(plo+phi)
         if (psat_res(pmid).ge.0.0_WP) then; plo=pmid; else; phi=pmid; end if
         if (phi-plo.le.1.0e-14_WP*max(phi,1.0_WP)) exit
      end do
      psat=0.5_WP*(plo+phi); conv=.true.
   contains
      real(WP) function psat_res(p_)
         real(WP), intent(in) :: p_
         psat_res=log(this%get_pvsat(p_,T_))-log(p_)
      end function psat_res
   end subroutine get_psat

   real(WP) function get_pvsat(this,pl_,T_)
      implicit none
      class(relax_igmix_sg), intent(inout) :: this
      real(WP), intent(in)  :: pl_,T_
      get_pvsat=exp(this%AS+(this%BS+this%ES*pl_)/T_)*T_**this%CS*(pl_+this%liq%pinf)**this%DS
   end function get_pvsat

   real(WP) function get_xv(this,Yv_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: Yv_
      get_xv=Yv_*Ma/(Yv_*Ma+(1.0_WP-Yv_)*Mv)
   end function get_xv

   !> Equilibrium T from energy conservation in liquid-vapor-gas mixture (SG form: no co-volume)
   real(WP) function get_T_lvg(this,p_,Yv_,rho0,rhoA0,Pjump_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP), intent(in) :: p_,Yv_,rho0,rhoA0,Pjump_
      get_T_lvg=(1.0_WP-Yv_)/((rho0*(1.0_WP-Yv_)-rhoA0)*(this%liq%gamma-1.0_WP)*this%liq%cv/(p_+this%liq%pinf)+&
      &           rhoA0*((this%gas%gamma(this%indV)-1.0_WP)*this%gas%cv(this%indV)*Yv_+ &
      &                  (this%gas%gamma(this%indA)-1.0_WP)*this%gas%cv(this%indA)*(1.0_WP-Yv_))/(p_-Pjump_))
   end function get_T_lvg

   !> Quadratic coefficients for the equilibrium-T equation (SG form: PinfG=0)
   subroutine get_coeffs_lv(this,p_eq,rho0,rhoe0,cvG,GammaG,qG,Pjump,ap,bp,dp,dapdp,dbpdp,ddpdp)
      implicit none
      class(relax_igmix_sg), intent(in)  :: this
      real(WP), intent(in)  :: p_eq,rho0,rhoe0,cvG,GammaG,qG,Pjump
      real(WP), intent(out) :: ap,bp,dp,dapdp,dbpdp,ddpdp
      ap=rho0*this%liq%cv*cvG*((GammaG-1.0_WP)*(p_eq+this%liq%gamma*this%liq%pinf)-(this%liq%gamma-1.0_WP)*p_eq   +&
      &  (this%liq%gamma-1.0_WP)*Pjump)
      bp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv*p_eq-(GammaG-1.0_WP)*cvG*(p_eq+this%liq%pinf))              +&
      &  rho0*((GammaG-1.0_WP)*cvG*this%liq%q*(p_eq+this%liq%pinf)-(this%liq%gamma-1.0_WP)*this%liq%cv*qG*p_eq) +&
      &  cvG*p_eq*(p_eq+this%liq%pinf)-this%liq%cv*p_eq*(p_eq+this%liq%gamma*this%liq%pinf)                     -&
      &  Pjump*(cvG-this%liq%cv)*p_eq-Pjump*this%liq%pinf*(cvG-this%liq%gamma*this%liq%cv)                      +&
      &  Pjump*this%liq%cv*(this%liq%gamma-1.0_WP)*(rho0*qG-rhoe0)
      dp=(qG-this%liq%q)*(p_eq+this%liq%pinf)*(p_eq-Pjump)
      dapdp=rho0*this%liq%cv*cvG*(GammaG-this%liq%gamma)
      dbpdp=rhoe0*((this%liq%gamma-1.0_WP)*this%liq%cv-(GammaG-1.0_WP)*cvG)                                     +&
      &     rho0*((GammaG-1.0_WP)*cvG*this%liq%q-(this%liq%gamma-1.0_WP)*this%liq%cv*qG)                        +&
      &     cvG*(2.0_WP*p_eq+this%liq%pinf)-this%liq%cv*(2.0_WP*p_eq+this%liq%gamma*this%liq%pinf)              -&
      &     Pjump*(cvG-this%liq%cv)
      ddpdp=(qG-this%liq%q)*(2.0_WP*p_eq+this%liq%pinf-Pjump)
   end subroutine get_coeffs_lv

   !> SG-form energy-conserving equilibrium pressure at given (frozen) VF, liquid pressure PL
   !> with PL-PG=Pjump_ (used by NASG override, and by apply()'s fixed-VF energy-split stage)
   real(WP) function get_p_eq(this,VF_,Q0_,qG_,gammaG_,Pjump_)
      implicit none
      class(relax_igmix_sg), intent(in) :: this
      real(WP),               intent(in) :: VF_
      real(WP), dimension(:), intent(in) :: Q0_
      real(WP),               intent(in) :: qG_,gammaG_,Pjump_
      get_p_eq=(sum(Q0_(3:4))-Q0_(1)*this%liq%q-VF_*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)-Q0_(2)*qG_+&
      &         (1.0_WP-VF_)*Pjump_/(gammaG_-1.0_WP))/&
      &        (VF_/(this%liq%gamma-1.0_WP)+(1.0_WP-VF_)/(gammaG_-1.0_WP))
   end function get_p_eq

end module relax_igmix_sg_class
