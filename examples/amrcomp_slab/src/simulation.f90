!> AMR compressible slab test case
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use amrgrid_class,     only: amrgrid
   use amrmpcomp_class,   only: amrmpcomp
   use amrviz_class,      only: amrviz
   use amrdata_class,     only: amrdata
   use timetracker_class, only: timetracker
   use event_class,       only: event
   use monitor_class,     only: monitor
   use amrio_class,       only: amrio
   implicit none
   private
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach,PMIX,RHOMIX
   
   !> Visualization
   type(event) :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time
   
   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,relaxfile

   !> Relaxation diagnostics (module-level for monitor binding)
   integer :: rlx_n_mixture,rlx_n_pass1_ok,rlx_n_need_cluster
   integer :: rlx_n_cluster_ok,rlx_n_skip_no_normal,rlx_n_skip_no_gas,rlx_n_skip_cluster_fail
   integer :: rlx_stage
   
   !> Stiffened gas EOS parameters (liquid and gas)
   real(WP) :: GammaL,PinfL,CvL
   real(WP) :: GammaG,PinfG,CvG

   !> Flow parameters
   real(WP) :: rhoG1,pG1,u1           !< Pre-shock gas state
   real(WP) :: rhoG2,pG2,u2           !< Post-shock gas state
   real(WP) :: rhoL1,pL1              !< Initial liquid state
   real(WP) :: M2,Xs                  !< Post-shock Mach and shock location
   real(WP) :: Ms                     !< Shock Mach number
   real(WP) :: density_ratio          !< rhoL1/rhoG1
   real(WP) :: ML                     !< Liquid Mach number
   real(WP) :: Reynolds,visc_ratio    !< Viscosity 
   real(WP) :: Prandtl ,diff_ratio    !< Heat diffusivity
   
   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Slab geometry
   real(WP) :: slab_left,slab_right

   !> Sponge parameters
   real(WP) :: Y_spg=3.0_WP              !< Sponge start distance from y-center
   real(WP) :: L_spg=1.0_WP              !< Sponge width

   !> Tagging parameters
   real(WP) :: vorticity_tag=huge(1.0_WP)
   real(WP) :: rho_ratio_tag=huge(1.0_WP)

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for slab (between slab_left and slab_right in x)
   function slab_levelset(xyz,t) result(G)
      use mathtools, only: twoPi
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,ymid,sl,sr,sin_pert,amp
      integer  :: per
      amp=0.2_WP
      per=4
      ymid=0.5_WP*(amr%yhi+amr%ylo)
      sin_pert=1.0_WP+amp*sin(twoPi*real(per,WP)*(xyz(2)-ymid)/(amr%yhi-amr%ylo))
      sl=slab_left *sin_pert
      sr=slab_right*sin_pert
      G=min(xyz(1)-sl,sr-xyz(1))
   end function slab_levelset

   !> Liquid EOS: P=f(RHO,I) - Stiffened gas
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> Liquid EOS: T=f(RHO,P)
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> Liquid EOS: C=f(RHO,P)
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(max(0.0_WP,GammaL*(P+PinfL)/RHO))
   end function get_CL
   !> Liquid EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IL=(P+GammaL*PinfL)/(RHO*(GammaL-1.0_WP))
   end function get_IL

   !> Gas EOS: P=f(RHO,I) - Ideal gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> Gas EOS: T=f(RHO,P)
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> Gas EOS: C=f(RHO,P)
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(max(0.0_WP,GammaG*(P+PinfG)/RHO))
   end function get_CG
   !> Gas EOS: I=f(RHO,P) (used for initialization)
   pure real(WP) function get_IG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_IG=(P+GammaG*PinfG)/(RHO*(GammaG-1.0_WP))
   end function get_IG

   !> Implicit mechanical relaxation for stiffened gas EOS pair
   !> Solves quadratic for equilibrium pressure Peq where PL=PG=Peq,
   !> then computes adjustments to VF and internal energies via p*dV work exchange.
   !> Conserves: phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7)
   ! subroutine P_relax_implicit(VF,Q,Peq,dVF,dQ)
   !    use amrmpcomp_class, only: VFlo,VFhi
   !    implicit none
   !    real(WP),               intent(in) :: VF
   !    real(WP), dimension(:), intent(in) :: Q
   !    real(WP), intent(inout) :: Peq
   !    real(WP), intent(out), optional :: dVF
   !    real(WP), dimension(:), intent(out), optional :: dQ
   !    real(WP) :: invG1G,invG1L,d0,d1,facG,facL,a,b,d
   !    ! Skip if any conserved quantity is non-positive (EOS undefined)
   !    if (any(Q(1:4).le.0.0_WP)) return
   !    ! Precompute EOS constants
   !    invG1L=1.0_WP/(GammaL-1.0_WP); d0=GammaL*PinfL*invG1L; d1=1.0_WP+invG1L
   !    ! Switch behavior depending on how the subroutine is called
   !    if (.not.present(dVF)) then
   !       ! First mode: return Peq by solving a*Peq^2 + b*Peq + d = 0
   !       invG1G=1.0_WP/(GammaG-1.0_WP); facG=GammaG*PinfG*invG1G; facL=invG1G+VF
   !       a=d1*facL-VF*(invG1G+1.0_WP)
   !       b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
   !       d=d0*(facG-Q(4))-Q(3)*facG
   !       if (b**2-4.0_WP*a*d.lt.0.0_WP) return
   !       Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
   !    else
   !       ! Second mode: use provided Peq to compute dVF and dQ
   !       dVF=(VF*Peq+Q(3))/(d1*Peq+d0)-VF
   !       dQ=0.0_WP; dQ(3)=-Peq*dVF; dQ(4)=+Peq*dVF
   !    end if
   ! end subroutine P_relax_implicit

   !> Implicit mechanical relaxation for stiffened gas EOS pair
   !> Solves quadratic for equilibrium pressure Peq where PL=PG=Peq,
   !> then adjusts VF and internal energies via p*dV work exchange.
   !> Conserves: phasic masses Q(1:2), total internal energy Q(3)+Q(4), momentum Q(5:7)
   !> If success is present, it is set to .true./.false. and the routine
   !> returns early on failure WITHOUT modifying VF or Q.
   subroutine P_relax_implicit(VF,Q,success)
      use amrmpcomp_class, only: VFlo,VFhi
      implicit none
      real(WP),               intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      logical, intent(out), optional :: success
      real(WP) :: a,b,disc,Peq,VFeq
      real(WP) :: invG1G,invG1L,d0,d1,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-2_WP
      if (present(success)) success=.false.
      ! Skip if any conserved quantity is non-positive (EOS undefined)
      if (any(Q(1:4).le.0.0_WP)) return
      ! Skip near-pure-liquid cells (gas density too low)
      if (VF.gt.VFhi) return
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Precompute EOS constants
      invG1G=1.0_WP/(GammaG-1.0_WP)
      invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=GammaL*PinfL*invG1L
      d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G
      facL=invG1G+VF
      ! Quadratic coefficients: a*Peq^2 + b*Peq + disc_term = 0
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      disc=d0*(facG-Q(4))-Q(3)*facG
      ! Solve for equilibrium pressure (positive root)
      if (b**2-4.0_WP*a*disc.lt.0.0_WP) return
      Peq=(-b+sqrt(b**2-4.0_WP*a*disc))/(2.0_WP*a)
      ! Bail if pressure is unphysical
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Equilibrium volume fraction from liquid energy constraint
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Update internal energies via p*dV work exchange
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
      if (present(success)) success=.true.
   end subroutine P_relax_implicit

   !> Cluster-based pressure relaxation (generic 3D)
   !> Pass 1: Attempt single-cell relaxation; mark failures
   !> Pass 2: For failures, find a pure-gas neighbor along the interface normal,
   !>         form a 2-cell cluster, relax the aggregate, distribute back
   subroutine apply_cluster_relax()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use amrmpcomp_class,  only: VFlo,VFhi
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_INTEGER
      use parallel,         only: MPI_REAL_WP
      implicit none
      integer :: lvl,i,j,k,i1,j1,k1,nf,i2,j2,k2
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pSD,pQ,pCL,pCG,pPLIC
      real(WP) :: dx,dy,dz
      ! Flags
      logical, dimension(:,:,:), allocatable :: needs_clustering, processed
      integer :: ilo,ihi,jlo,jhi,klo,khi
      ! Normal computation
      real(WP) :: nx,ny,nz,nmag
      ! Neighbor search
      integer :: in,jn,kn
      integer, dimension(3,6) :: offsets
      real(WP) :: alignment,best_alignment
      logical :: found_gas,ok
      ! Cluster work variables
      real(WP) :: VF_c,Q_c(7),VF_orig(2),Q_orig(7,2)
      real(WP) :: VF_c_old,dVF_cluster,dQ3,dQ4
      real(WP) :: w1,w2,wsum,dVF1,dVF2
      real(WP), parameter :: wmin=1.0e-16_WP
      ! Temporaries for try
      real(WP) :: VF_tmp,Q_tmp(7)
      ! ===== Diagnostics =====
      integer :: ierr,n_print_cluster,n_print_skip
      integer, parameter :: max_print=5
      real(WP) :: rhoG

      ! 6 face-connected neighbor offsets
      offsets(:,1)=[+1, 0, 0]
      offsets(:,2)=[-1, 0, 0]
      offsets(:,3)=[ 0,+1, 0]
      offsets(:,4)=[ 0,-1, 0]
      offsets(:,5)=[ 0, 0,+1]
      offsets(:,6)=[ 0, 0,-1]

      ! Operate on finest level only
      lvl=fs%amr%maxlvl
      dx=fs%amr%dx(lvl); dy=fs%amr%dy(lvl); dz=fs%amr%dz(lvl)

      ! Initialize counters
      rlx_n_mixture=0; rlx_n_pass1_ok=0; rlx_n_need_cluster=0
      rlx_n_cluster_ok=0; rlx_n_skip_no_normal=0; rlx_n_skip_no_gas=0; rlx_n_skip_cluster_fail=0
      n_print_cluster=0; n_print_skip=0

      call fs%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
         pSD  =>fs%SD%dataptr(mfi)
         pQ   =>fs%Q%mf(lvl)%dataptr(mfi)
         pCL  =>fs%CL%dataptr(mfi)
         pCG  =>fs%CG%dataptr(mfi)
         pPLIC=>fs%PLIC%dataptr(mfi)

         ! Grown box (includes ghost cells)
         bx=mfi%growntilebox(fs%nover)
         ilo=bx%lo(1); ihi=bx%hi(1)
         jlo=bx%lo(2); jhi=bx%hi(2)
         klo=bx%lo(3); khi=bx%hi(3)

         ! Allocate per-box flag arrays
         allocate(needs_clustering(ilo:ihi,jlo:jhi,klo:khi)); needs_clustering=.false.
         allocate(processed       (ilo:ihi,jlo:jhi,klo:khi)); processed=.false.

         ! ===============================================================
         ! Pass 1: Try single-cell relaxation on temporaries, mark failures
         ! ===============================================================
         do k=klo,khi; do j=jlo,jhi; do i=ilo,ihi
            ! Only relax interfacial cells
            if (pVF(i,j,k,1).lt.VFlo.or.pVF(i,j,k,1).gt.VFhi) cycle
            rlx_n_mixture=rlx_n_mixture+1
            ! Try relaxation on copies
            VF_tmp=pVF(i,j,k,1)
            Q_tmp =pQ(i,j,k,:)
            call P_relax_implicit(VF_tmp,Q_tmp,ok)
            if (ok) then
               ! Apply result directly
               pVF(i,j,k,1)=VF_tmp
               pQ(i,j,k,:) =Q_tmp
               call cleanup_cell(i,j,k,pVF,pQ,pCL,pCG,dx,dy,dz)
               processed(i,j,k)=.true.
               rlx_n_pass1_ok=rlx_n_pass1_ok+1
            else
               needs_clustering(i,j,k)=.true.
               rlx_n_need_cluster=rlx_n_need_cluster+1
            end if
         end do; end do; end do

         ! Sync needs_clustering? Use fill boundary if sure we are away from the coarse-fine boundary and we are at the max level

         ! ===============================================================
         ! Pass 2: Cluster failures with neighbor along normal
         ! ===============================================================
         do k1=klo,khi; do j1=jlo,jhi; do i1=ilo,ihi
            if (.not.needs_clustering(i1,j1,k1)) cycle
            if (processed(i1,j1,k1)) cycle

            ! ============================================================
            ! Get interface normal from PLIC reconstruction
            ! PLIC stores (nx, ny, nz, d); normal points from liquid to gas
            ! ============================================================
            nx=pPLIC(i1,j1,k1,1)
            ny=pPLIC(i1,j1,k1,2)
            nz=pPLIC(i1,j1,k1,3)
            nmag=sqrt(nx**2+ny**2+nz**2)
            if (nmag.gt.0.0_WP) then
               nx=nx/nmag; ny=ny/nmag; nz=nz/nmag
            else
               ! Degenerate: no normal available, skip
               rlx_n_skip_no_normal=rlx_n_skip_no_normal+1
               if (n_print_skip.lt.max_print) then
                  n_print_skip=n_print_skip+1
                  rhoG=0.0_WP; if (pVF(i1,j1,k1,1).lt.1.0_WP) rhoG=pQ(i1,j1,k1,2)/(1.0_WP-pVF(i1,j1,k1,1))
                  print '(A,I6,A,3I5,A,ES12.5,A,ES12.5,A)', &
                     '  [SKIP:no_norm] n=',time%n,' ijk=',i1,j1,k1,' VF=',pVF(i1,j1,k1,1),' rhoG=',rhoG,' PLIC=(0,0,0)'
               end if
               cycle
            end if

            ! ============================================================
            ! Search 6 face-connected neighbors for the best gas cell
            ! "Best" = pure gas (VF < VFlo) AND most aligned with normal
            ! ============================================================
            found_gas=.false.
            best_alignment=-huge(1.0_WP)
            i2=i1; j2=j1; k2=k1   ! will hold the chosen gas neighbor
            do nf=1,6
               in=i1+offsets(1,nf)
               jn=j1+offsets(2,nf)
               kn=k1+offsets(3,nf)
               ! Bounds check (Remove after figuring out how to sync needs_clustering)
               if (in.lt.lbound(pVF,1).or.in.gt.ubound(pVF,1)) cycle
               if (jn.lt.lbound(pVF,2).or.jn.gt.ubound(pVF,2)) cycle
               if (kn.lt.lbound(pVF,3).or.kn.gt.ubound(pVF,3)) cycle
               ! Must be pure gas
               ! if (pVF(in,jn,kn,1).ge.VFlo) cycle
               ! Skip if liquid
               if (pVF(in,jn,kn,1).gt.VFhi) cycle
               ! Compute alignment with normal
               alignment=abs(real(offsets(1,nf),WP)*nx+real(offsets(2,nf),WP)*ny+real(offsets(3,nf),WP)*nz)
               if (alignment.gt.best_alignment) then
                  best_alignment=alignment
                  i2=in; j2=jn; k2=kn
                  found_gas=.true.
               end if
            end do

            ! If no pure-gas neighbor found, skip (leave unrelaxed)
            if (.not.found_gas) then
               rlx_n_skip_no_gas=rlx_n_skip_no_gas+1
               if (n_print_skip.lt.max_print) then
                  n_print_skip=n_print_skip+1
                  rhoG=0.0_WP; if (pVF(i1,j1,k1,1).lt.1.0_WP) rhoG=pQ(i1,j1,k1,2)/(1.0_WP-pVF(i1,j1,k1,1))
                  print '(A,I6,A,3I5,A,ES12.5,A,ES12.5,A,3F8.4,A)', &
                     '  [SKIP:no_gas] n=',time%n,' ijk=',i1,j1,k1,' VF=',pVF(i1,j1,k1,1),' rhoG=',rhoG, &
                     ' n=(',nx,ny,nz,')'
               end if
               cycle
            end if

            ! ============================================================
            ! Form 2-cell cluster and relax the aggregate
            ! Both cells have equal volume -> cluster average = simple mean
            ! ============================================================
            VF_orig(1)=pVF(i1,j1,k1,1)
            VF_orig(2)=pVF(i2,j2,k2,1)
            Q_orig(:,1)=pQ(i1,j1,k1,:)
            Q_orig(:,2)=pQ(i2,j2,k2,:)

            VF_c=0.5_WP*(VF_orig(1)+VF_orig(2))
            Q_c =0.5_WP*(Q_orig(:,1)+Q_orig(:,2))
            VF_c_old=VF_c

            ! Relax the cluster aggregate
            call P_relax_implicit(VF_c,Q_c,ok)
            if (.not.ok) then
               ! Cluster relaxation also failed, skip
               rlx_n_skip_cluster_fail=rlx_n_skip_cluster_fail+1
               if (n_print_skip.lt.max_print) then
                  n_print_skip=n_print_skip+1
                  rhoG=0.0_WP; if (VF_c.lt.1.0_WP) rhoG=Q_c(2)/(1.0_WP-VF_c)
                  print '(A,I6,A,3I5,A,ES12.5,A,3I5,A,ES12.5,A,ES12.5)', &
                     '  [SKIP:clst_fail] n=',time%n,' ijk=',i1,j1,k1,' VF=',VF_orig(1), &
                     ' gas=',i2,j2,k2,' VF_c=',VF_c,' rhoG_c=',rhoG
               end if
               cycle
            end if

            ! Cluster-level deltas
            dVF_cluster=VF_c-VF_c_old
            dQ3=Q_c(3)-0.5_WP*(Q_orig(3,1)+Q_orig(3,2))
            dQ4=Q_c(4)-0.5_WP*(Q_orig(4,1)+Q_orig(4,2))

            ! Print successful clustering info
            rlx_n_cluster_ok=rlx_n_cluster_ok+1
            if (n_print_cluster.lt.max_print) then
               n_print_cluster=n_print_cluster+1
               rhoG=0.0_WP; if (VF_orig(1).lt.1.0_WP) rhoG=Q_orig(2,1)/(1.0_WP-VF_orig(1))
               print '(A,I6,A,3I5,A,ES12.5,A,ES12.5,A,3I5,A,ES12.5,A,ES12.5)', &
                  '  [CLUSTER] n=',time%n,' ijk=',i1,j1,k1,' VF=',VF_orig(1),' rhoG=',rhoG, &
                  ' gas=',i2,j2,k2,' dVF=',dVF_cluster,' dQ3=',dQ3

            end if

            ! ============================================================
            ! Distribute VF change
            ! ============================================================
            w1=max(amr%cell_vol(amr%maxlvl)*pSD(i1,j1,k1,1)*min(VF_orig(1),1.0_WP-VF_orig(1)),wmin)
            w2=max(amr%cell_vol(amr%maxlvl)*pSD(i2,j2,k2,1)*min(VF_orig(2),1.0_WP-VF_orig(2)),wmin)
            wsum=w1+w2

            dVF1=dVF_cluster*(w1/wsum)
            dVF2=dVF_cluster*(w2/wsum)

            ! Clip and redistribute remainder
            if (VF_orig(1)+dVF1.gt.1.0_WP) then
               dVF1=1.0_WP-VF_orig(1)
               dVF2=dVF_cluster-dVF1
            else if (VF_orig(1)+dVF1.lt.0.0_WP) then
               dVF1=-VF_orig(1)
               dVF2=dVF_cluster-dVF1
            end if
            ! if (VF_orig(2)+dVF2.gt.1.0_WP) then
            !    dVF2=1.0_WP-VF_orig(2)
            ! else if (VF_orig(2)+dVF2.lt.0.0_WP) then
            !    dVF2=-VF_orig(2)
            ! end if

            ! Apply VF changes
            pVF(i1,j1,k1,1)=VF_orig(1)+dVF1
            pVF(i2,j2,k2,1)=VF_orig(2)+dVF2

            ! ============================================================
            ! Distribute Q changes: gradient-preserving (additive delta)
            ! Only Q(3) and Q(4) change from relaxation
            ! ============================================================
            pQ(i1,j1,k1,3)=Q_orig(3,1)+dQ3
            pQ(i1,j1,k1,4)=Q_orig(4,1)+dQ4
            pQ(i2,j2,k2,3)=Q_orig(3,2)+dQ3
            pQ(i2,j2,k2,4)=Q_orig(4,2)+dQ4

            ! Cleanup both cells
            ! call cleanup_cell(i1,j1,k1,pVF,pQ,pCL,pCG,dx,dy,dz)
            ! call cleanup_cell(i2,j2,k2,pVF,pQ,pCL,pCG,dx,dy,dz)
            processed(i1,j1,k1)=.true.
            processed(i2,j2,k2)=.true.
            ! How do I make other processors aware of this?

         end do; end do; end do

         deallocate(needs_clustering,processed)
      end do
      call fs%amr%mfiter_destroy(mfi)

      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_mixture,          1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_pass1_ok,         1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_need_cluster,     1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_cluster_ok,       1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_skip_no_normal,   1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_skip_no_gas,      1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,rlx_n_skip_cluster_fail,1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      if (rlx_n_need_cluster.gt.0) call relaxfile%write()

   contains

      !> Cleanup a cell after relaxation: handle pure-liquid/pure-gas transitions
      subroutine cleanup_cell(ic,jc,kc,pVF,pQ,pCL,pCG,dx,dy,dz)
         implicit none
         integer, intent(in) :: ic,jc,kc
         real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pCL,pCG
         real(WP), intent(in) :: dx,dy,dz
         real(WP) :: xc,yc,zc
         xc=fs%amr%xlo+(real(ic,WP)+0.5_WP)*dx
         yc=fs%amr%ylo+(real(jc,WP)+0.5_WP)*dy
         zc=fs%amr%zlo+(real(kc,WP)+0.5_WP)*dz
         if (pVF(ic,jc,kc,1).lt.VFlo) then
            ! Pure gas
            pVF(ic,jc,kc,1)=0.0_WP
            pCL(ic,jc,kc,1:3)=[xc,yc,zc]
            pCG(ic,jc,kc,1:3)=[xc,yc,zc]
            pQ(ic,jc,kc,1)=0.0_WP
            pQ(ic,jc,kc,3)=0.0_WP
         else if (pVF(ic,jc,kc,1).gt.VFhi) then
            ! Pure liquid
            pVF(ic,jc,kc,1)=1.0_WP
            pCL(ic,jc,kc,1:3)=[xc,yc,zc]
            pCG(ic,jc,kc,1:3)=[xc,yc,zc]
            pQ(ic,jc,kc,2)=0.0_WP
            pQ(ic,jc,kc,4)=0.0_WP
         end if
      end subroutine cleanup_cell

   end subroutine apply_cluster_relax

   !> Compute viscosity: Sutherland for gas, VF-weighted blend with liquid
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pQ,pVisc,pBeta,pDiff,pRHOL,pRHOG
      real(WP) :: mu_g,mu_l,k_g,k_l,y_cc,y_dist,blend,nu_spg,mu_spg
      real(WP), parameter :: Tmax_visc=10.0_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cdiff=0.1_WP
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Get maximum allowable kinematic viscosity in the sponge at finest level
      nu_spg=max_cfl*amr%min_meshsize(amr%clvl())**2/(4.0_WP*time%dt)
      ! Loop over levels
      do lvl=0,amr%clvl()
         ! Loop over domain
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pQ=>fs%Q%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Gas viscosity from Sutherland
               mu_g=(1.0_WP+Suth_T)*min(pTG(i,j,k,1),Tmax_visc)**Suth_n/(Reynolds*(min(pTG(i,j,k,1),Tmax_visc)+Suth_T))
               ! Liquid viscosity from ratio
               mu_l=visc_ratio*Reynolds**(-1.0_WP)
               ! Mixture viscosity
               !pVisc(i,j,k,1)=pVF(i,j,k,1)*mu_l+(1.0_WP-pVF(i,j,k,1))*mu_g ! Arithmetic averaging
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps)) ! Harmonic averaging
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Gas heat diffusivity: k=Cv*Gamma*mu/Pr
               k_g=GammaG*CvG*mu_g/Prandtl
               ! Liquid heat diffusivity from ratio
               k_l=diff_ratio*GammaG*CvG/(Reynolds*Prandtl)
               ! Mixture diffusivity
               !pDiff(i,j,k,1)=pVF(i,j,k,1)*k_l+(1.0_WP-pVF(i,j,k,1))*k_g ! Arithmetic averaging
               pDiff(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(k_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(k_g,myeps)) ! Harmonic averaging
               ! Apply sponge layer viscosity (distance from domain y-center)
               y_cc=amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl)
               y_dist=abs(y_cc-0.5_WP*(amr%ylo+amr%yhi))
               if (y_dist.gt.Y_spg) then
                  blend=min((y_dist-Y_spg)/L_spg,1.0_WP)**2
                  mu_spg=nu_spg/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),myeps)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),myeps))
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*mu_spg)
                  pDiff(i,j,k,1)=max(pDiff(i,j,k,1),Cdiff*blend*mu_spg)
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> Compute mixture quantities: PMIX = VF*PL + (1-VF)*PG, RHOMIX = VF*RHOL + (1-VF)*RHOG
   subroutine compute_mix()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pPL,pPG,pRHOL,pRHOG,pPMIX,pRHOMIX
      real(WP) :: vf
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF    =>fs%VF%mf(lvl)%dataptr(mfi)
            pPL    =>fs%PL%mf(lvl)%dataptr(mfi)
            pPG    =>fs%PG%mf(lvl)%dataptr(mfi)
            pRHOL  =>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG  =>fs%RHOG%mf(lvl)%dataptr(mfi)
            pPMIX  =>PMIX%mf(lvl)%dataptr(mfi)
            pRHOMIX=>RHOMIX%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               vf=pVF(i,j,k,1)
               pPMIX(i,j,k,1)  =vf*pPL(i,j,k,1)+(1.0_WP-vf)*pPG(i,j,k,1)
               pRHOMIX(i,j,k,1)=vf*pRHOL(i,j,k,1)+(1.0_WP-vf)*pRHOG(i,j,k,1)
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine compute_mix

   !> User init callback - set Q and VF/barycenters for a slab at rest with a shock
   subroutine shockslab_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom, only: initialize_volume_moments
      use amrmpcomp_class, only: VFlo
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,IEL,x_cc,rhoG,pG,uG,H
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get internal energy of liquid
      IEL=get_IL(rhoL1,pL1)
      ! Use passed ba/dm since grid is being constructed
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         ! Get pointers to data
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         ! Loop over grown tilebox
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute VF and barycenters from levelset
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=slab_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Compute local gas state from shock profile
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
            rhoG=rhoG1+(rhoG2-rhoG1)*H
            pG  =pG1  +(pG2  -pG1  )*H
            uG  =u1   +(u2   -u1   )*H
            ! Set conserved variables: Q=(VF*rhoL, (1-VF)*rhoG, VF*rhoL*IL, (1-VF)*rhoG*IG, rho_mix*U, 0, 0)
            pQ(i,j,k,1)=(       myVF)*rhoL1
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*get_IG(rhoG,pG)
            pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*uG
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockslab_init

   !> Apply inflow BC at low-x (face=1)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,pQ)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      integer :: i,j,k
      select case (face)
       case (1)  ! X-LOW: Dirichlet inflow with post-shock (gas only, no liquid)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            pQ(i,j,k,1)=0.0_WP                  ! No liquid
            pQ(i,j,k,2)=rhoG2                   ! Gas density
            pQ(i,j,k,3)=0.0_WP                  ! No liquid energy
            pQ(i,j,k,4)=rhoG2*get_IG(rhoG2,pG2) ! Gas internal energy
            pQ(i,j,k,5)=rhoG2*u2                ! X-momentum
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
         end do; end do; end do
      end select
   end subroutine shock_dirichlet

   !> Tagger based on rho ratio (1D)
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP) :: irho_cc,irho_xp,irho_xm,irho_yp,irho_ym,irho_zp,irho_zm
      real(WP) :: vort_x,vort_y,vort_z,vort_mag
      real(WP) :: rho_max,rho_min,rho_nb,rho_ratio,y_dist
      integer :: i,j,k,ii,jj,kk
      ! Get mesh size
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz
      ! Recast tags
      tags=tags_ptr
      ! Compute tags
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         ! Get pointers to data
         tagarr=>tags%dataPtr(mfi)
         pQ=>solver%Q%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Compute vorticity and tag based on it
            irho_cc=1.0_WP/max(sum(pQ(i  ,j,  k,  1:2)),solver%rho_floor)
            irho_xp=1.0_WP/max(sum(pQ(i+1,j,  k,  1:2)),solver%rho_floor)
            irho_xm=1.0_WP/max(sum(pQ(i-1,j,  k,  1:2)),solver%rho_floor)
            irho_yp=1.0_WP/max(sum(pQ(i,  j+1,k,  1:2)),solver%rho_floor)
            irho_ym=1.0_WP/max(sum(pQ(i,  j-1,k,  1:2)),solver%rho_floor)
            irho_zp=1.0_WP/max(sum(pQ(i,  j,  k+1,1:2)),solver%rho_floor)
            irho_zm=1.0_WP/max(sum(pQ(i,  j,  k-1,1:2)),solver%rho_floor)
            vort_x=(pQ(i,j+1,k,7)*irho_yp-pQ(i,j-1,k,7)*irho_ym)*0.5_WP*dyi-(pQ(i,j,k+1,6)*irho_zp-pQ(i,j,k-1,6)*irho_zm)*0.5_WP*dzi
            vort_y=(pQ(i,j,k+1,5)*irho_zp-pQ(i,j,k-1,5)*irho_zm)*0.5_WP*dzi-(pQ(i+1,j,k,7)*irho_xp-pQ(i-1,j,k,7)*irho_xm)*0.5_WP*dxi
            vort_z=(pQ(i+1,j,k,6)*irho_xp-pQ(i-1,j,k,6)*irho_xm)*0.5_WP*dxi-(pQ(i,j+1,k,5)*irho_yp-pQ(i,j-1,k,5)*irho_ym)*0.5_WP*dyi
            vort_mag=sqrt(vort_x**2+vort_y**2+vort_z**2)
            if (vort_mag.gt.vorticity_tag) tagarr(i,j,k,1)=SETtag
            ! Compute density ratio in 3x3x3 stencil and tag based on it
            rho_max=solver%rho_floor; rho_min=huge(1.0_WP)
            do kk=-1,1; do jj=-1,1; do ii=-1,1
               rho_nb=sum(pQ(i+ii,j+jj,k+kk,1:2))
               rho_max=max(rho_max,rho_nb)
               rho_min=min(rho_min,max(rho_nb,solver%rho_floor))
            end do; end do; end do
            rho_ratio=rho_max/rho_min
            ! Only tag outside sponge for density ratio (prevent over-refining in sponge)
            y_dist=abs(solver%amr%ylo+(real(j,WP)+0.5_WP)*dy-0.5_WP*(solver%amr%ylo+solver%amr%yhi))
            if (rho_ratio.gt.rho_ratio_tag.and.(y_dist.lt.Y_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: A,B,C
         ! Gas EoS parameters (ideal gas = stiffened gas with Pinf=0)
         call param_read('GammaG',GammaG)
         PinfG=0.0_WP
         ! Liquid EoS: gamma only, PinfL is computed below
         call param_read('GammaL',GammaL)
         ! Shock parameters (gas phase, uses GammaG)
         call param_read('Gas Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rhoG2=1, Deltau=1, T2=1
         rhoG2=1.0_WP
         pG2=1.0_WP/(GammaG*M2**2)
         ! Quadratic for rhoG1: A*rhoG1^2 - B*rhoG1 + C = 0
         A=2.0_WP*GammaG*pG2+(GammaG-1.0_WP)
         B=4.0_WP*GammaG*pG2+(GammaG+1.0_WP)
         C=2.0_WP*GammaG*pG2
         rhoG1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)  ! smaller root for compression
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rhoG1)
         u2=u1-1.0_WP
         pG1=pG2-rhoG1/(1.0_WP-rhoG1)
         if (pG1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         ! Shock Mach number
         Ms=u1/sqrt(GammaG*pG1/rhoG1)
         ! Shift to lab frame: pre-shock stationary
         u2=1.0_WP
         u1=0.0_WP
         ! CvG from T2=1
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         ! Liquid state from density ratio and liquid Mach number
         call param_read('Density ratio',density_ratio)
         call param_read('Liquid Mach number',ML)
         rhoL1=density_ratio
         PinfL=rhoL1/(GammaL*ML**2)-pG1
         pL1=pG1                                                        ! Default: pressure equilibrium
         call param_read('Liquid pressure offset',pL1,default=0.0_WP)    ! Optional: shift to negative
         pL1=pG1+pL1                                                     ! Apply offset (negative = tension)
         CvL=(pL1+PinfL)/(rhoL1*(GammaL-1.0_WP)*get_TG(rhoG1,pG1))      ! Force thermal equilibrium
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Viscosity ratio',visc_ratio)
         call param_read('Diffusivity ratio',diff_ratio)
         call param_read('Sutherland exponent',Suth_n)
         call param_read('Sutherland temperature',Suth_T)
         ! Log
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5)') rhoG2,pG2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," ML=",es12.5)') rhoL1,pL1,ML; call log(message)
         write(message,'("[Liquid] GammaL=",es12.5," PinfL=",es12.5," CvL=",es12.5)') GammaL,PinfL,CvL; call log(message)
         write(message,'("[Gas]    GammaG=",es12.5," PinfG=",es12.5," CvG=",es12.5)') GammaG,PinfG,CvG; call log(message)
         write(message,'("[Visc]   Re=",es12.5," mu*=",es12.5," Suth_n=",es12.5," Suth_T=",es12.5)') Reynolds,visc_ratio,Suth_n,Suth_T; call log(message)
      end block init_eos_and_flow
      
      ! Read slab geometry
      read_slab: block
         call param_read('Slab left',slab_left)
         call param_read('Slab right',slab_right)
      end block read_slab
      
      ! Initialize AMR grid
      create_amrgrid: block
         amr%name='amrcomp_slab'
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=-05.0_WP;    amr%xhi=+15.0_WP
         amr%ylo=-2.5_WP;     amr%yhi=+2.5_WP
         amr%zlo=-0.15625_WP; amr%zhi=+0.15625_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         integer :: restart_step
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap
         use amrmpcomp_class,  only: BC_GAS
         ! Create flow solver
         call fs%initialize(amr=amr,name='slab')
         ! Provide thermodynamic model (6 EOS pointers)
         fs%getPL=>get_PL; fs%getCL=>get_CL; fs%getTL=>get_TL
         fs%getPG=>get_PG; fs%getCG=>get_CG; fs%getTG=>get_TG
         ! Provide pressure relaxation model
         fs%relax=>P_relax_implicit
         ! Set initial conditions
         fs%user_mpcomp_init=>shockslab_init
         ! Set BCs
         if (.not.amr%xper) then
            fs%lo_bc(1)=BC_GAS
            fs%Q%lo_bc(1,:)=amrex_bc_ext_dir
            fs%Q%hi_bc(1,:)=amrex_bc_foextrap
            fs%user_mpcomp_bc=>shock_dirichlet
         end if
      end block create_solver
      
      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: amrex_interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=7,ng=0,interp=amrex_interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=amrex_interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=amrex_interp_none); call Mach%register()
         call PMIX%initialize(amr,name='PMIX',ncomp=1,ng=0,interp=amrex_interp_none); call PMIX%register()
         call RHOMIX%initialize(amr,name='RHOMIX',ncomp=1,ng=0,interp=amrex_interp_none); call RHOMIX%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_mpcomp_tagging=>my_tagger
         call param_read('Tagging rho ratio',rho_ratio_tag)
         call param_read('Tagging vorticity',vorticity_tag,default=huge(1.0_WP))
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC
            call fs%build_plic(time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         call compute_mix()
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint
      
      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,'slab',use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%U,1,'U')
         call viz%add_scalar(fs%V,1,'V')
         call viz%add_scalar(fs%W,1,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(PMIX,1,'PMIX')
         call viz%add_scalar(RHOMIX,1,'RHOMIX')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz
      
      ! Create monitors
      create_monitors: block
         ! Get solver info and cfl
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%VFmin,'VFmin')
         call mfile%add_column(fs%VFmax,'VFmax')
         call mfile%add_column(fs%VFint,'VFint')
         call mfile%add_column(fs%dPmax,'dPmax')
         call mfile%write()
         ! Create relaxation monitor
         relaxfile=monitor(amRoot=amr%amRoot,name='relaxation')
         call relaxfile%add_column(time%n,'Timestep number')
         call relaxfile%add_column(time%t,'Time')
         call relaxfile%add_column(rlx_stage,'Stage')
         call relaxfile%add_column(rlx_n_mixture,'Interfacial cells')
         call relaxfile%add_column(rlx_n_pass1_ok,'Pass 1 OK')
         call relaxfile%add_column(rlx_n_need_cluster,'Need cluster')
         call relaxfile%add_column(rlx_n_cluster_ok,'Cluster OK')
         call relaxfile%add_column(rlx_n_skip_no_normal,'Skip no normal')
         call relaxfile%add_column(rlx_n_skip_no_gas,'Skip no gas')
         call relaxfile%add_column(rlx_n_skip_cluster_fail,'Skip cluster fail')
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLc_y,'CFLc_y')
         call cflfile%add_column(fs%CFLc_z,'CFLc_z')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLa_y,'CFLa_y')
         call cflfile%add_column(fs%CFLa_z,'CFLa_z')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLv_y,'CFLv_y')
         call cflfile%add_column(fs%CFLv_z,'CFLv_z')
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
         ! Create grid monitor
         gridfile=monitor(amRoot=amr%amRoot,name='grid')
         call gridfile%add_column(time%n,'Timestep')
         call gridfile%add_column(time%t,'Time')
         call gridfile%add_column(amr%nlevels,'Nlvl')
         call gridfile%add_column(amr%nboxes,'Nbox')
         call gridfile%add_column(amr%ncells,'Ncell')
         call gridfile%add_column(amr%compression,'Compression')
         call gridfile%add_column(amr%maxRSS,'Maximum RSS')
         call gridfile%add_column(amr%minRSS,'Minimum RSS')
         call gridfile%add_column(amr%avgRSS,'Average RSS')
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         ! Full routine times (max across ranks = wall-clock cost)
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         ! Compute loop times (max = slowest rank, min = fastest rank)
         call tfile%add_column(fs%wtmax_prim,'prim_max')
         call tfile%add_column(fs%wtmin_prim,'prim_min')
         call tfile%add_column(fs%wtmax_sl,'sl_max')
         call tfile%add_column(fs%wtmin_sl,'sl_min')
         call tfile%add_column(fs%wtmax_fv,'fv_max')
         call tfile%add_column(fs%wtmin_fv,'fv_min')
         call tfile%add_column(fs%wtmax_div,'div_max')
         call tfile%add_column(fs%wtmin_div,'div_min')
         call tfile%add_column(fs%wtmax_plicnet,'plicnet_max')
         call tfile%add_column(fs%wtmin_plicnet,'plicnet_min')
         call tfile%add_column(fs%wtmax_polygon,'polygon_max')
         call tfile%add_column(fs%wtmin_polygon,'polygon_min')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%add_column(fs%nmixed_min,'mixed_min')
         call tfile%write()
      end block create_monitors

   end subroutine simulation_init
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      
      ! Perform time integration
      do while (.not.time%done())
         
         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Remember old state
         call fs%store_old()
         
         ! ===== RK2 Stage 1: dQdt = f(t, Q) =====
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=0.5_WP*time%dt,time=time%t)
         
         ! ===== RK2 Stage 2: Q* = Qold + dt/2*dQdt, dQdt* = f(t+dt/2, Q*) =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call check_Q('RK1   ')
         rlx_stage=1
         call apply_cluster_relax()
         call check_Q('RELAX1')
         call fs%get_dQdt(Q=fs%Q,dQdt=dQdt,dt=time%dt,time=time%t+0.5_WP*time%dt)

         ! ===== RK2 Final: Q = Qold + dt*dQdt* =====
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call check_Q('RK2   ')
         rlx_stage=2
         call apply_cluster_relax()
         call check_Q('RELAX2')

         ! Rebuild PLIC
         call fs%build_plic(time%t)

         ! Recompute primitive variables
         call fs%get_primitive(fs%Q)

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(fs%U,fs%V,fs%W)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         call compute_mix()

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/slab_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         
      end do
      
   end subroutine simulation_run
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      ! Finalize time
      call time%finalize()
      ! Finalize grid
      call amr%finalize()
      call regrid_evt%finalize()
      ! Finalize solver
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      call PMIX%finalize()
      call RHOMIX%finalize()
      ! Finalize visualization
      call viz%finalize()
      call viz_evt%finalize()
      ! Finalize checkpoint
      call save_evt%finalize()
      call io%finalize()
      ! Finalize monitoring
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
      call relaxfile%finalize()
   end subroutine simulation_final
   
   !> Diagnostic: scan Q/primitives for extreme values
   subroutine check_Q(label)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use amrmpcomp_class,  only: VFlo,VFhi
      use ieee_arithmetic,  only: ieee_is_nan
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM,MPI_MAX,MPI_MIN,MPI_INTEGER
      use parallel,         only: MPI_REAL_WP
      implicit none
      character(len=*), intent(in) :: label
      integer :: lvl,i,j,k,nbad,nnan,nnan_print,nbad_print,ierr,flvl
      integer :: nbadF,nnanF
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF
      real(WP) :: IL,IG,PG,PL
      real(WP) :: Q1min,Q2min,Q3min,Q4min,VFmin,VFmax,dPmax,PGmax,PLmin
      real(WP) :: Q1minF,Q2minF,Q3minF,Q4minF,VFminF,VFmaxF,dPmaxF,PGmaxF,PLminF
      logical  :: is_bad,has_nan,is_finest
      integer, parameter :: max_nan_print=10,max_bad_print=5
      flvl=fs%amr%clvl()
      ! Initialize coarse-level counters
      nbad=0; nnan=0; nnan_print=0; nbad_print=0
      Q1min=huge(1.0_WP); Q2min=huge(1.0_WP); Q3min=huge(1.0_WP); Q4min=huge(1.0_WP)
      VFmin=huge(1.0_WP); VFmax=-huge(1.0_WP)
      dPmax=0.0_WP; PGmax=-huge(1.0_WP); PLmin=huge(1.0_WP)
      ! Initialize finest-level counters
      nbadF=0; nnanF=0
      Q1minF=huge(1.0_WP); Q2minF=huge(1.0_WP); Q3minF=huge(1.0_WP); Q4minF=huge(1.0_WP)
      VFminF=huge(1.0_WP); VFmaxF=-huge(1.0_WP)
      dPmaxF=0.0_WP; PGmaxF=-huge(1.0_WP); PLminF=huge(1.0_WP)
      ! Loop over ALL levels
      do lvl=0,flvl
         is_finest=(lvl.eq.flvl)
         call fs%amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ =>fs%Q%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! NaN check
               has_nan=ieee_is_nan(pQ(i,j,k,1)).or.ieee_is_nan(pQ(i,j,k,2)).or. &
               &       ieee_is_nan(pQ(i,j,k,3)).or.ieee_is_nan(pQ(i,j,k,4)).or. &
               &       ieee_is_nan(pQ(i,j,k,5)).or.ieee_is_nan(pQ(i,j,k,6)).or. &
               &       ieee_is_nan(pQ(i,j,k,7)).or.ieee_is_nan(pVF(i,j,k,1))
               if (has_nan) then
                  if (is_finest) then; nnanF=nnanF+1; else; nnan=nnan+1; end if
                  if (nnan_print.lt.max_nan_print) then
                     nnan_print=nnan_print+1
                     print '(A,A,A,I6,A,I2,A,3I5)', '  [',label,'] NaN n=',time%n,' l=',lvl,' ijk=',i,j,k
                     print '(4(A,ES20.13))', '    Q1=',pQ(i,j,k,1),' Q2=',pQ(i,j,k,2),' Q3=',pQ(i,j,k,3),' Q4=',pQ(i,j,k,4)
                     print '(4(A,ES20.13))', '    Q5=',pQ(i,j,k,5),' Q6=',pQ(i,j,k,6),' Q7=',pQ(i,j,k,7),' VF=',pVF(i,j,k,1)
                  end if
                  cycle
               end if
               ! Track Q1/Q3 min (liquid present)
               if (pVF(i,j,k,1).ge.VFlo) then
                  if (is_finest) then
                     if (pQ(i,j,k,1).lt.Q1minF) Q1minF=pQ(i,j,k,1)
                     if (pQ(i,j,k,3).lt.Q3minF) Q3minF=pQ(i,j,k,3)
                  else
                     if (pQ(i,j,k,1).lt.Q1min) Q1min=pQ(i,j,k,1)
                     if (pQ(i,j,k,3).lt.Q3min) Q3min=pQ(i,j,k,3)
                  end if
               end if
               ! Track Q2/Q4 min (gas present)
               if (pVF(i,j,k,1).le.VFhi) then
                  if (is_finest) then
                     if (pQ(i,j,k,2).lt.Q2minF) Q2minF=pQ(i,j,k,2)
                     if (pQ(i,j,k,4).lt.Q4minF) Q4minF=pQ(i,j,k,4)
                  else
                     if (pQ(i,j,k,2).lt.Q2min) Q2min=pQ(i,j,k,2)
                     if (pQ(i,j,k,4).lt.Q4min) Q4min=pQ(i,j,k,4)
                  end if
               end if
               ! Track VF extrema
               if (is_finest) then
                  if (pVF(i,j,k,1).lt.VFminF) VFminF=pVF(i,j,k,1)
                  if (pVF(i,j,k,1).gt.VFmaxF) VFmaxF=pVF(i,j,k,1)
               else
                  if (pVF(i,j,k,1).lt.VFmin) VFmin=pVF(i,j,k,1)
                  if (pVF(i,j,k,1).gt.VFmax) VFmax=pVF(i,j,k,1)
               end if
               ! Compute phasic pressures
               PG=0.0_WP; PL=0.0_WP
               if (pVF(i,j,k,1).le.VFhi.and.pQ(i,j,k,2).gt.0.0_WP) then
                  IG=pQ(i,j,k,4)/pQ(i,j,k,2)
                  PG=(GammaG-1.0_WP)*pQ(i,j,k,2)/(1.0_WP-pVF(i,j,k,1))*IG-GammaG*PinfG
               end if
               if (pVF(i,j,k,1).ge.VFlo.and.pQ(i,j,k,1).gt.0.0_WP) then
                  IL=pQ(i,j,k,3)/pQ(i,j,k,1)
                  PL=(GammaL-1.0_WP)*pQ(i,j,k,1)/pVF(i,j,k,1)*IL-GammaL*PinfL
               end if
               if (is_finest) then
                  if (PG.gt.PGmaxF) PGmaxF=PG
                  if (PL.lt.PLminF.and.pVF(i,j,k,1).ge.VFlo) PLminF=PL
                  if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
                     if (abs(PL-PG).gt.dPmaxF) dPmaxF=abs(PL-PG)
                  end if
               else
                  if (PG.gt.PGmax) PGmax=PG
                  if (PL.lt.PLmin.and.pVF(i,j,k,1).ge.VFlo) PLmin=PL
                  if (pVF(i,j,k,1).ge.VFlo.and.pVF(i,j,k,1).le.VFhi) then
                     if (abs(PL-PG).gt.dPmax) dPmax=abs(PL-PG)
                  end if
               end if
               ! Bad cell check
               is_bad=pQ(i,j,k,1).lt.-1.0e-10_WP.or.pQ(i,j,k,2).lt.-1.0e-10_WP.or. &
               &      pQ(i,j,k,3).lt.-1.0e-10_WP.or.pQ(i,j,k,4).lt.-1.0e-10_WP
               if (is_bad) then
                  if (is_finest) then; nbadF=nbadF+1; else; nbad=nbad+1; end if
                  if (nbad_print.lt.max_bad_print) then
                     nbad_print=nbad_print+1
                     print '(A,A,A,I6,A,I2,A,3I5)', '  [',label,'] BAD n=',time%n,' l=',lvl,' ijk=',i,j,k
                     print '(4(A,ES20.13))', '    Q1=',pQ(i,j,k,1),' Q2=',pQ(i,j,k,2),' Q3=',pQ(i,j,k,3),' Q4=',pQ(i,j,k,4)
                     print '(4(A,ES20.13))', '    Q5=',pQ(i,j,k,5),' Q6=',pQ(i,j,k,6),' Q7=',pQ(i,j,k,7),' VF=',pVF(i,j,k,1)
                  end if
               end if
            end do; end do; end do
         end do
         call fs%amr%mfiter_destroy(mfi)
      end do
      ! Reductions for coarse levels
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q1min ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q2min ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q3min ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q4min ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFmin ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFmax ,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,PLmin ,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,PGmax ,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dPmax ,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nbad  ,1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nnan  ,1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      ! Reductions for finest level
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q1minF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q2minF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q3minF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,Q4minF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFminF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,VFmaxF,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,PLminF,1,MPI_REAL_WP,MPI_MIN,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,PGmaxF,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,dPmaxF,1,MPI_REAL_WP,MPI_MAX,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nbadF ,1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,nnanF ,1,MPI_INTEGER,MPI_SUM,fs%amr%comm,ierr)
      ! Summary (rank 0 only): F=finest, C=coarse
      if (fs%amr%amRoot) then
         print '(A,A,A,I6,A,4(A,ES12.5),2(A,ES20.13),3(A,ES12.5),A,I6,A,I6)', &
            '[',label,'] n=',time%n,' F', &
            ' Q1m=',Q1minF,' Q2m=',Q2minF,' Q3m=',Q3minF,' Q4m=',Q4minF, &
            ' VFm=',VFminF,' VFM=',VFmaxF, &
            ' PLm=',PLminF,' PGM=',PGmaxF,' dPM=',dPmaxF, &
            ' bad:',nbadF,' nan:',nnanF
         print '(A,A,A,I6,A,4(A,ES12.5),2(A,ES20.13),3(A,ES12.5),A,I6,A,I6)', &
            '[',label,'] n=',time%n,' C', &
            ' Q1m=',Q1min,' Q2m=',Q2min,' Q3m=',Q3min,' Q4m=',Q4min, &
            ' VFm=',VFmin,' VFM=',VFmax, &
            ' PLm=',PLmin,' PGM=',PGmax,' dPM=',dPmax, &
            ' bad:',nbad,' nan:',nnan
      end if
   end subroutine check_Q

end module simulation
