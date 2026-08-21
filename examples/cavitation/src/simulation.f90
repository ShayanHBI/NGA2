!> Cavitation/condensation in tensioned liquid
module simulation
   use precision,              only: WP
   use string,                 only: str_medium
   use amrgrid_class,          only: amrgrid
   use amrmpcomp_class,        only: amrmpcomp
   use amrviz_class,           only: amrviz
   use amrdata_class,          only: amrdata
   use timetracker_class,      only: timetracker
   use event_class,            only: event
   use monitor_class,          only: monitor
   use amrio_class,            only: amrio
   use stiffened_gas_class,    only: stiffened_gas
   use nasg_class,             only: nasg
   use igmix_class,            only: igmix
   use relax_igmix_sg_class,   only: relax_igmix_sg,Prelax,PTrelax,PTgrelax
   use relax_igmix_nasg_class, only: relax_igmix_nasg
   implicit none
   private

   public :: simulation_init,simulation_run,simulation_final

   !> AMR grid
   type(amrgrid), target :: amr

   !> Timetracker and compressible multiphase solver
   type(timetracker) :: time
   type(amrmpcomp), target :: fs
   type(amrdata) :: dQdt,Umag,Mach

   !> Visualization
   type(event)  :: viz_evt
   type(amrviz) :: viz

   ! Regrid parameters
   type(event) :: regrid_evt

   ! Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical  :: restarted
   real(WP) :: restart_time
   integer  :: restart_step

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> EOS parameters (SG/NASG liquid and ideal gas)
   real(WP) :: GammaL,PinfL,qL,qpL,CvL,CpL,bL
   real(WP) :: GammaA,qA,qpA,CvA,CpA
   real(WP) :: GammaV,qV,qpV,CvV,CpV

   !> EOS and relaxation
   class(stiffened_gas), allocatable,    target, save :: eosL         !< Liquid pure-substance EOS (SG or NASG)
   type(igmix),                          target, save :: mixG         !< Gas mixture (vapor + air)
   class(relax_igmix_nasg), allocatable, target, save :: relax_model  !< Relaxation model (SG or NASG liquid)
   character(len=str_medium), save :: liquid_eos_type,relaxation_type
   character(len=str_medium), save :: case_name

   !> Boundary condition mode: 'dirichlet' (cavitation only, sustained expansion)
   !> or 'wall' (cavitation+condensation, reflective walls bounce the wave back)
   character(len=str_medium), save :: bc_type

   !> Case parameters
   real(WP) :: T0,p0              !< Uniform initial liquid temperature/pressure
   real(WP) :: rhoL0,eL0          !< Liquid density/energy at (p0,T0), used by the Dirichlet BC
   real(WP) :: U0,r0              !< Gaussian pulse radial profile: amplitude [1/s] and radius [m]
   real(WP) :: p_cav              !< Cavitation onset pressure threshold
   real(WP) :: Tctol              !< Condensation temperature tolerance
   real(WP) :: VFratmax           !< Max per-relax-call VF change factor (relax_igmix_sg default: 10.0)
   real(WP) :: muG,muL            !< Dynamic viscosities
   real(WP) :: PrL,PrG            !< Prandtl numbers

   !> Domain dimensions
   real(WP) :: Lx,Ly              !< Domain lengths [m]

   !> Tagging parameters
   real(WP) :: vorticity_tag=huge(1.0_WP)
   real(WP) :: rho_ratio_tag=huge(1.0_WP)
   real(WP) :: divergence_tag=huge(1.0_WP)

   !> Time stepping
   real(WP) :: dt_init

contains

   !> Levelset function for a domain that is liquid everywhere
   function levelset_liquid(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP
   end function levelset_liquid

   function get_radial_velocity(r) result(Ur)
      real(WP), intent(in) :: r
      real(WP) :: Ur
      Ur=U0*r*exp(-(r/r0)**2)
   end function get_radial_velocity

   !> debug: print Q, RHOL/RHOG, PL/PG, TL/TG for one target cell at a given pipeline stage,
   !> gated to a single timestep so we can see exactly where an unphysical value first appears
   subroutine debug_probe(label)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      character(len=*), intent(in) :: label
      integer :: lvl
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ
      real(WP) :: VFc,RHOL,RHOG,PL,PG,TL,TG,eL,eG,Yv
      real(WP), dimension(2) :: y
      ! debug: nucleation onset, domain-center cell (i,j)=(128,128) at maxlvl=3 (Base nx=32 ->
      ! 256 finest cells over Lx=0.03m -> center sits on the face between i,j=127 and 128),
      ! n=6-16 spans first VF<1 in the monitor (n=7) through several relaxation steps, to check
      ! whether PL/PG and TL/TG converge to equality (full pTg) or stay split (fallback path)
      integer, parameter :: ti=128,tj=128,tk=0
      if (time%n.lt.6.or.time%n.gt.16) return
      lvl=amr%maxlvl
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF=>fs%VF%mf(lvl)%dataptr(mfi)
         pQ =>fs%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         if (ti.ge.bx%lo(1).and.ti.le.bx%hi(1).and.tj.ge.bx%lo(2).and.tj.le.bx%hi(2).and.tk.ge.bx%lo(3).and.tk.le.bx%hi(3)) then
            VFc=pVF(ti,tj,tk,1)
            RHOL=-1.0_WP; RHOG=-1.0_WP; PL=0.0_WP; PG=0.0_WP; TL=0.0_WP; TG=0.0_WP
            if (VFc.gt.0.0_WP.and.pQ(ti,tj,tk,1).gt.0.0_WP) then
               RHOL=pQ(ti,tj,tk,1)/VFc
               eL=pQ(ti,tj,tk,3)/pQ(ti,tj,tk,1)
               PL=fs%liq%get_p_from_rho_e(rho=RHOL,e=eL,y=[1.0_WP])
               TL=fs%liq%get_T_from_p_rho(p=PL,rho=RHOL,y=[1.0_WP])
            end if
            if (VFc.lt.1.0_WP.and.pQ(ti,tj,tk,2).gt.0.0_WP) then
               RHOG=pQ(ti,tj,tk,2)/(1.0_WP-VFc)
               eG=pQ(ti,tj,tk,4)/pQ(ti,tj,tk,2)
               Yv=pQ(ti,tj,tk,8)/pQ(ti,tj,tk,2)
               y=[Yv,1.0_WP-Yv]
               PG=fs%gas%get_p_from_rho_e(rho=RHOG,e=eG,y=y)
               TG=fs%gas%get_T_from_p_rho(p=PG,rho=RHOG,y=y)
            end if
            print*,'PROBE[',trim(label),'] n=',time%n,' VF=',VFc,' Q=',pQ(ti,tj,tk,1:8),&
            &      ' RHOL=',RHOL,' RHOG=',RHOG,' PL=',PL,' PG=',PG,' TL=',TL,' TG=',TG ! debug
         end if
      end do
      call amr%mfiter_destroy(mfi)
   end subroutine debug_probe

   !> debug: scan for the interfacial cell with the smallest RHOG this step (the small-cell
   !> noise culprit), then dump its own VF and its 4 face-neighbor VF plus whether any neighbor
   !> qualifies as a merge_Q gas reservoir (VF<merge_VFhi) -- tests whether the noise correlates
   !> with merge_Q's claimant search coming up empty (gpd=0 -> gclaim=.false. -> cell never rescued)
   subroutine debug_scan_worst_gas()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      integer :: lvl,i,j,k,wi,wj,wk
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ
      real(WP) :: VFc,RHOG,wRHOG
      logical :: hasRes
      if (time%n.lt.7.or.time%n.gt.30) return
      lvl=amr%maxlvl
      wRHOG=huge(1.0_WP); wi=-huge(0); wj=-huge(0); wk=0
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF=>fs%VF%mf(lvl)%dataptr(mfi)
         pQ =>fs%Q%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            VFc=pVF(i,j,k,1)
            if (VFc.gt.0.0_WP.and.VFc.lt.1.0_WP.and.pQ(i,j,k,2).gt.0.0_WP) then
               RHOG=pQ(i,j,k,2)/(1.0_WP-VFc)
               if (RHOG.lt.wRHOG) then; wRHOG=RHOG; wi=i; wj=j; wk=k; end if
            end if
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      if (wRHOG.eq.huge(1.0_WP)) return
      ! Re-locate the worst cell to safely read its face neighbors (needs ghost access within
      ! the same FAB/tile)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF=>fs%VF%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         if (wi.ge.bx%lo(1).and.wi.le.bx%hi(1).and.wj.ge.bx%lo(2).and.wj.le.bx%hi(2).and.wk.ge.bx%lo(3).and.wk.le.bx%hi(3)) then
            hasRes=(pVF(wi-1,wj,wk,1).lt.fs%merge_VFhi).or.(pVF(wi+1,wj,wk,1).lt.fs%merge_VFhi).or.&
            &      (pVF(wi,wj-1,wk,1).lt.fs%merge_VFhi).or.(pVF(wi,wj+1,wk,1).lt.fs%merge_VFhi)
            print*,'WORSTGAS n=',time%n,' (i,j)=',wi,wj,' RHOG=',wRHOG,' VF=',pVF(wi,wj,wk,1),&
            &      ' nbrVF(xm,xp,ym,yp)=',pVF(wi-1,wj,wk,1),pVF(wi+1,wj,wk,1),pVF(wi,wj-1,wk,1),pVF(wi,wj+1,wk,1),&
            &      ' merge_VFhi=',fs%merge_VFhi,' hasReservoirNeighbor=',hasRes ! debug
         end if
      end do
      call amr%mfiter_destroy(mfi)
   end subroutine debug_scan_worst_gas

   !> Compute viscosity: constant gas and liquid, VF-weighted blend
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc,pBeta,pDiffL,pDiffG
      real(WP) :: mu_g,mu_l
      real(WP), parameter :: myeps=1.0e-15_WP
      ! Loop over levels
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Constant gas and liquid viscosity
               mu_g=muG; mu_l=muL
               ! Mixture viscosity (harmonic averaging)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps))
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Phasic heat diffusivities (zero when Pr=0, i.e. no heat/species diffusion)
               if (PrL.gt.0.0_WP) then; pDiffL(i,j,k,1)=mu_l*CpL/PrL; else; pDiffL(i,j,k,1)=0.0_WP; end if
               if (PrG.gt.0.0_WP) then; pDiffG(i,j,k,1)=mu_g*CpV/PrG; else; pDiffG(i,j,k,1)=0.0_WP; end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - set Q and VF/barycenters for a uniform liquid domain
   subroutine cavitation_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
      use amrmpcomp_class,  only: VFlo
      use param,            only: param_read
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz,myVF,rhoG_local,eG_local,rhoL_local,eL_local,Yv0,y(1:2)
      real(WP) :: x_cc,y_cc,r,Ur,rho_mix
      integer :: i,j,k
      integer, parameter :: nref=3
      ! Get mesh size
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      ! Get initial vapor mass fraction
      call param_read('Initial vapor mass fraction',Yv0)
      y=[Yv0,1.0_WP-Yv0]
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
            ! Compute VF and barycenters from levelset (liquid everywhere)
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=levelset_liquid,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            ! Store volume fraction
            pVF(i,j,k,1)=myVF
            ! Store barycenters
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            ! Get liquid density and internal energy from p and T
            eL_local  =eosL%get_e_from_p_T  (p=p0,T=T0,y=[1.0_WP])
            rhoL_local=eosL%get_rho_from_p_T(p=p0,T=T0,y=[1.0_WP])
            ! Get gas density and internal energy from p and T
            eG_local  =mixG%get_e_from_p_T(p=p0,T=T0,y=y)
            rhoG_local=mixG%get_rho_from_p_T(p=p0,T=T0,y=y)
            ! Set conserved variables
            pQ(i,j,k,1)=(       myVF)*rhoL_local
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG_local
            pQ(i,j,k,3)=pQ(i,j,k,1)*eL_local
            pQ(i,j,k,4)=pQ(i,j,k,2)*eG_local
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            y_cc=solver%amr%ylo+(real(j,WP)+0.5_WP)*dy
            r=sqrt(x_cc**2+y_cc**2)
            Ur=get_radial_velocity(r)
            rho_mix=pQ(i,j,k,1)+pQ(i,j,k,2)
            pQ(i,j,k,5)=rho_mix*Ur*x_cc/r
            pQ(i,j,k,6)=rho_mix*Ur*y_cc/r
            pQ(i,j,k,7)=0.0_WP
            pQ(i,j,k,8)=pQ(i,j,k,2)*Yv0
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine cavitation_init

   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrmpcomp_class,  only: VFhi
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF
      real(WP) :: dx,dy,dz,dxi,dyi,dzi
      real(WP) :: irho_cc,irho_xp,irho_xm,irho_yp,irho_ym,irho_zp,irho_zm
      real(WP) :: vort_x,vort_y,vort_z,vort_mag
      real(WP) :: div_mag
      real(WP) :: rho_max,rho_min,rho_nb,rho_ratio
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
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         ! Loop over tile
         bx=mfi%tilebox()
         if (time.eq.0.0_WP) then
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Get local inverse densities
               irho_xp=1.0_WP/max(sum(pQ(i+1,j,  k,  1:2)),solver%rho_floor)
               irho_xm=1.0_WP/max(sum(pQ(i-1,j,  k,  1:2)),solver%rho_floor)
               irho_yp=1.0_WP/max(sum(pQ(i,  j+1,k,  1:2)),solver%rho_floor)
               irho_ym=1.0_WP/max(sum(pQ(i,  j-1,k,  1:2)),solver%rho_floor)
               irho_zp=1.0_WP/max(sum(pQ(i,  j,  k+1,1:2)),solver%rho_floor)
               irho_zm=1.0_WP/max(sum(pQ(i,  j,  k-1,1:2)),solver%rho_floor)
               ! Compute divergence and tag based on it
               div_mag=(pQ(i+1,j,k,5)*irho_xp-pQ(i-1,j,k,5)*irho_xm)*0.5_WP*dxi &
               &      +(pQ(i,j+1,k,6)*irho_yp-pQ(i,j-1,k,6)*irho_ym)*0.5_WP*dyi &
               &      +(pQ(i,j,k+1,7)*irho_zp-pQ(i,j,k-1,7)*irho_zm)*0.5_WP*dzi
               if (abs(div_mag).gt.divergence_tag) tagarr(i,j,k,1)=SETtag
            end do; end do; end do
         end if
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Get local inverse densities
            irho_cc=1.0_WP/max(sum(pQ(i  ,j,  k,  1:2)),solver%rho_floor)
            irho_xp=1.0_WP/max(sum(pQ(i+1,j,  k,  1:2)),solver%rho_floor)
            irho_xm=1.0_WP/max(sum(pQ(i-1,j,  k,  1:2)),solver%rho_floor)
            irho_yp=1.0_WP/max(sum(pQ(i,  j+1,k,  1:2)),solver%rho_floor)
            irho_ym=1.0_WP/max(sum(pQ(i,  j-1,k,  1:2)),solver%rho_floor)
            irho_zp=1.0_WP/max(sum(pQ(i,  j,  k+1,1:2)),solver%rho_floor)
            irho_zm=1.0_WP/max(sum(pQ(i,  j,  k-1,1:2)),solver%rho_floor)
            ! Compute vorticity and tag based on it
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
            if (rho_ratio.gt.rho_ratio_tag) tagarr(i,j,k,1)=SETtag
            ! Gas and interfacial cells
            if (pVF(i,j,k,1).lt.VFhi) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> User BC callback - radial outward velocity Dirichlet condition on every boundary
   !> Only used when Boundary condition = dirichlet (sustained expansion, cavitation only)
   subroutine radial_dirichlet_bc(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: dx,dy,dz,x,y,z,r,U_r,ur,vr
      integer :: i,j,k
      dx=amr%dx(lvl); dy=amr%dy(lvl); dz=amr%dz(lvl)
      do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
         ! Position: staggered in the component's own direction, cell-centered otherwise
         if (comp.eq.'U') then
            x=amr%xlo+real(i,WP)*dx
         else
            x=amr%xlo+(real(i,WP)+0.5_WP)*dx
         end if
         if (comp.eq.'V') then
            y=amr%ylo+real(j,WP)*dy
         else
            y=amr%ylo+(real(j,WP)+0.5_WP)*dy
         end if
         r=sqrt(x**2+y**2)
         U_r=get_radial_velocity(r)
         ur=U_r*x/r; vr=U_r*y/r
         select case (comp)
         case ('U')
            p(i,j,k,1)=ur
         case ('V')
            p(i,j,k,1)=vr
         case ('W')
            p(i,j,k,1)=0.0_WP
         case ('Q')
            p(i,j,k,1)=rhoL0; p(i,j,k,2)=0.0_WP
            p(i,j,k,3)=rhoL0*eL0; p(i,j,k,4)=0.0_WP
            p(i,j,k,5)=rhoL0*ur; p(i,j,k,6)=rhoL0*vr; p(i,j,k,7)=0.0_WP
            p(i,j,k,8)=0.0_WP
         end select
      end do; end do; end do
   end subroutine radial_dirichlet_bc

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log,die
         use string,   only: str_long
         character(len=str_long) :: message
         ! Get liquid EOS parameters
         call param_read('Liquid specific heat capacity at constant volume',CvL)
         call param_read('Liquid specific heat capacity ratio',GammaL)
         call param_read('Liquid reference energy shift',qL)
         call param_read('Liquid reference entropy shift',qpL)
         call param_read('Liquid stiffening pressure',PinfL)
         call param_read('Liquid co-volume',bL,default=0.0_WP)
         ! Get vapor EOS parameters
         call param_read('Vapor specific heat capacity at constant volume',CvV)
         call param_read('Vapor specific heat capacity ratio',GammaV)
         call param_read('Vapor reference energy shift',qV)
         call param_read('Vapor reference entropy shift',qpV)
         ! Get air EOS parameters
         call param_read('Air specific heat capacity at constant volume',CvA)
         call param_read('Air specific heat capacity ratio',GammaA)
         call param_read('Air reference energy shift',qA)
         call param_read('Air reference entropy shift',qpA)
         ! Derive Cp from Gamma and Cv
         CpL=GammaL*CvL
         CpA=GammaA*CvA
         CpV=GammaV*CvV
         call param_read('Liquid viscosity',muL)
         call param_read('Gas viscosity',muG)
         ! Prandtl numbers
         call param_read('Liquid Prandtl number',PrL)
         call param_read('Gas Prandtl number',PrG)
         ! Uniform liquid state and boundary velocity
         call param_read('Liquid temperature',T0)
         call param_read('Liquid pressure',p0)
         call param_read('Pulse amplitude',U0)
         call param_read('Pulse radius',r0)
         ! Cavitation onset: nucleate when liquid pressure drops below p_cav
         call param_read('Cavitation pressure threshold',p_cav,default=huge(1.0_WP))
         ! Condensation onset: nucleate only when vapor is subcooled by more than Tctol below Tsat
         call param_read('Condensation temperature tolerance',Tctol,default=0.0_WP)
         call param_read('Max VF change factor',VFratmax,default=10.0_WP)
         ! Domain dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         ! Select boundary condition
         call param_read('Boundary condition',bc_type)
         select case (trim(bc_type))
         case ('dirichlet','wall')
         case default
            call die('Boundary condition has to be either dirichlet or wall')
         end select
         ! Select liquid eos
         call param_read('Liquid EOS type',liquid_eos_type)
         select case (trim(liquid_eos_type))
         ! case ('SG')
         !    allocate(stiffened_gas  :: eosL)
         !    allocate(relax_igmix_sg :: relax_model)
         case ('NASG')
            allocate(nasg              :: eosL)
            allocate(relax_igmix_nasg  :: relax_model)
         case default
            call die('[simulation] Unknown Liquid EOS type: '//trim(liquid_eos_type))
         end select
         ! Select relaxation model
         call param_read('Relaxation type',relaxation_type)
         select case(relaxation_type)
         case('p','pT','pTg')
         case default
            call die('Relaxation type has to be either p, pT, or pTg')
         end select
         ! Build case name
         case_name='cavitation_'//trim(bc_type)//'_'//trim(liquid_eos_type)//'_'//trim(relaxation_type)
         ! Initialize EOS objects
         select type (eosL)
         type is (stiffened_gas)
            call eosL%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=qL,qp=qpL,name='water')
         type is (nasg)
            call eosL%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=qL,qp=qpL,name='water')
            eosL%brhomax=1e10_WP
         end select
         ! Liquid density and energy at the uniform initial state (used by the Dirichlet BC)
         rhoL0=eosL%get_rho_from_p_T(p=p0,T=T0,y=[1.0_WP])
         eL0  =eosL%get_e_from_p_T  (p=p0,T=T0,y=[1.0_WP])
         ! Gas mixture: species 1=vapor (transported), species 2=air (carrier)
         call mixG%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,qA],qp=[qpV,qpA], &
         &                     species_names=['vapor','air  '],name='gas')
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name=trim(case_name)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('AMReX blocking factor',amr%nbloc,default=8)
         amr%xlo=-0.5_WP*Lx; amr%xhi=+0.5_WP*Lx
         amr%ylo=-0.5_WP*Ly; amr%yhi=+0.5_WP*Ly
         amr%zlo=-0.5_WP*Ly; amr%zhi=+0.5_WP*Ly
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         call param_read('Max level',amr%maxlvl)
         ! Enable quasi-2D
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zper=.true.
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves here
      handle_restart: block
         ! Initialize IO object
         call io%initialize(amr=amr,nfiles=1)
         ! Check if restarting
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         if (restarted) restart_dir='restart/'//trim(case_name)//'/'//trim(adjustl(restart_dir))
         ! If restarting, read header
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Initial dt',dt_init,default=time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=dt_init
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use param,            only: param_read
         use amrex_amr_module, only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrmpcomp_class,  only: BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         ! Use piecewise-linear face interpolation -- FaceDivFree requires ratio==2 in all dirs
         ! but this case is quasi-2D with ref_ratio_z=1
         fs%interp_vel=interp_face_lin
         ! Assign materials and create flow solver
         fs%liq=>eosL; fs%gas=>mixG
         call fs%initialize(amr=amr,name=trim(case_name))
         ! Get surface tension
         call param_read('Surface tension',fs%sigma)
         ! Provide relaxation model
         fs%relax=>relax_model
         ! Set initial conditions via cavitation callback
         fs%user_init=>cavitation_init
         select case (trim(bc_type))
         case ('dirichlet')
            ! Radial outward velocity is continuously re-imposed at the boundary
            ! -> sustained expansion, vapor keeps growing (cavitation only)
            if (.not.amr%xper) then
               fs%Q%lo_bc(1,:)=amrex_bc_ext_dir; fs%Q%hi_bc(1,:)=amrex_bc_ext_dir
               fs%U%lo_bc(1,1)=amrex_bc_ext_dir; fs%U%hi_bc(1,1)=amrex_bc_ext_dir
               fs%V%lo_bc(1,1)=amrex_bc_ext_dir; fs%V%hi_bc(1,1)=amrex_bc_ext_dir
               fs%W%lo_bc(1,1)=amrex_bc_ext_dir; fs%W%hi_bc(1,1)=amrex_bc_ext_dir
            end if
            if (.not.amr%yper) then
               fs%Q%lo_bc(2,:)=amrex_bc_ext_dir; fs%Q%hi_bc(2,:)=amrex_bc_ext_dir
               fs%U%lo_bc(2,1)=amrex_bc_ext_dir; fs%U%hi_bc(2,1)=amrex_bc_ext_dir
               fs%V%lo_bc(2,1)=amrex_bc_ext_dir; fs%V%hi_bc(2,1)=amrex_bc_ext_dir
               fs%W%lo_bc(2,1)=amrex_bc_ext_dir; fs%W%hi_bc(2,1)=amrex_bc_ext_dir
            end if
            if (.not.amr%zper) then
               fs%Q%lo_bc(3,:)=amrex_bc_ext_dir; fs%Q%hi_bc(3,:)=amrex_bc_ext_dir
               fs%U%lo_bc(3,1)=amrex_bc_ext_dir; fs%U%hi_bc(3,1)=amrex_bc_ext_dir
               fs%V%lo_bc(3,1)=amrex_bc_ext_dir; fs%V%hi_bc(3,1)=amrex_bc_ext_dir
               fs%W%lo_bc(3,1)=amrex_bc_ext_dir; fs%W%hi_bc(3,1)=amrex_bc_ext_dir
            end if
            fs%user_bc=>radial_dirichlet_bc
         case ('wall')
            ! Solid free-slip walls on every non-periodic boundary. Unlike a Dirichlet
            ! inflow/outflow BC, a reflective wall bounces the outgoing rarefaction wave
            ! (launched once by the initial expansion pulse) back inward as a compression
            ! wave, which recompresses -- and via the pTg relaxation model, condenses --
            ! the vapor generated at the center
            if (.not.amr%xper) then
               fs%lo_bc(1)=BC_REFLECT; fs%hi_bc(1)=BC_REFLECT
               fs%Q%lo_bc(1,:)=amrex_bc_foextrap;    fs%Q%hi_bc(1,:)=amrex_bc_foextrap
               fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd; fs%Q%hi_bc(1,5)=amrex_bc_reflect_odd
               fs%U%lo_bc(1,:)=amrex_bc_reflect_odd; fs%U%hi_bc(1,:)=amrex_bc_reflect_odd
               fs%V%lo_bc(1,:)=amrex_bc_foextrap;    fs%V%hi_bc(1,:)=amrex_bc_foextrap
               fs%W%lo_bc(1,:)=amrex_bc_foextrap;    fs%W%hi_bc(1,:)=amrex_bc_foextrap
            end if
            if (.not.amr%yper) then
               fs%lo_bc(2)=BC_REFLECT; fs%hi_bc(2)=BC_REFLECT
               fs%Q%lo_bc(2,:)=amrex_bc_foextrap;    fs%Q%hi_bc(2,:)=amrex_bc_foextrap
               fs%Q%lo_bc(2,6)=amrex_bc_reflect_odd; fs%Q%hi_bc(2,6)=amrex_bc_reflect_odd
               fs%V%lo_bc(2,:)=amrex_bc_reflect_odd; fs%V%hi_bc(2,:)=amrex_bc_reflect_odd
               fs%U%lo_bc(2,:)=amrex_bc_foextrap;    fs%U%hi_bc(2,:)=amrex_bc_foextrap
               fs%W%lo_bc(2,:)=amrex_bc_foextrap;    fs%W%hi_bc(2,:)=amrex_bc_foextrap
            end if
            if (.not.amr%zper) then
               fs%lo_bc(3)=BC_REFLECT; fs%hi_bc(3)=BC_REFLECT
               fs%Q%lo_bc(3,:)=amrex_bc_foextrap;    fs%Q%hi_bc(3,:)=amrex_bc_foextrap
               fs%Q%lo_bc(3,7)=amrex_bc_reflect_odd; fs%Q%hi_bc(3,7)=amrex_bc_reflect_odd
               fs%W%lo_bc(3,:)=amrex_bc_reflect_odd; fs%W%hi_bc(3,:)=amrex_bc_reflect_odd
               fs%U%lo_bc(3,:)=amrex_bc_foextrap;    fs%U%hi_bc(3,:)=amrex_bc_foextrap
               fs%V%lo_bc(3,:)=amrex_bc_foextrap;    fs%V%hi_bc(3,:)=amrex_bc_foextrap
            end if
         end select
         ! Relaxation model: initialize then set cavitation/condensation thresholds and dispatch model
         call relax_model%initialize(liq=eosL,gas=mixG,indV=1,indA=2)
         relax_model%p_cav=p_cav
         relax_model%Tctol=Tctol
         select case (trim(relaxation_type))
         case ('p');   relax_model%model=Prelax
         case ('pT');  relax_model%model=PTrelax
         case ('pTg'); relax_model%model=PTgrelax
         end select
         relax_model%RHOGmin=0.01_WP
         relax_model%VFratmax=VFratmax
         relax_model%vol=amr%cell_vol(amr%maxlvl)
         fs%merge_sick=100.0_WP
         ! fs%rho_floor=1.0e-3_WP
         fs%Pmin_liq=-0.9_WP*eosL%pinf
         fs%Tmin_liq=300.0_WP
         fs%Pmin_gas=1.0e-4_WP
         fs%Tmin_gas=0.1_WP
         relax_model%Pmin_liq=fs%Pmin_liq; relax_model%Tmin_liq=fs%Tmin_liq
         relax_model%Pmin_gas=fs%Pmin_gas; relax_model%Tmin_gas=fs%Tmin_gas
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Initialize regridding
      init_regridding: block
         ! KnapSack load balancing
         amr%lb_strat=1
         ! Create regridding event
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         ! Set case-specific tagging
         fs%user_tagging=>my_tagger
         call param_read('Tagging vorticity',vorticity_tag)
         call param_read('Tagging rho ratio',rho_ratio_tag)
         call param_read('Tagging divergence',divergence_tag)
         ! Build the grid
         if (restarted) then
            ! Restore grid hierarchy from checkpoint
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            ! Restore solver state
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
            ! Force a regrid immediately after restoring checkpoint
            call amr%regrid(baselvl=0,time=time%t)
         else
            ! Fresh start
            call amr%init_from_scratch(time=time%t)
            ! Build PLIC and reset moments
            call fs%build_plic(time%t)
            call fs%reset_moments()
            call fs%Q%fill(time=time%t)
            call fs%get_primitive(Q=fs%Q)
            call fs%build_subVF()
            call fs%get_face_velocity()
            call fs%average_down_velocity()
            call fs%fill_velocity(time=time%t)
         end if
         ! Compute viscosities
         call get_viscosities()
         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)
         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         ! Create checkpoint save event
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period', save_evt%tper,default=-1.0_WP)
         ! Let solver self-register for checkpointing
         call fs%register_checkpoint(io)
         ! Add dt to checkpoint save
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         ! Create visualization object
         call viz%initialize(amr,trim(case_name),use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%visc,1,'visc')
         call viz%add_scalar(fs%beta,1,'beta')
         call viz%add_scalar(fs%C,1,'C')
         call viz%add_scalar(fs%Yg,1,'Yv')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         call viz%add_scalar(fs%IL,1,'IL')
         call viz%add_scalar(fs%IG,1,'IG')
         call viz%add_surfmesh(fs%smesh,'plic')
         ! Create visualization output event
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         ! Write initial state
         if (viz_evt%occurs().or.restarted) call viz%write(time=time%t)
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
         call mfile%add_column(fs%Ygmin(1),'Yvmin')
         call mfile%add_column(fs%Ygmax(1),'Yvmax')
         call mfile%add_column(fs%Qint(8),'Vapor mass')
         call mfile%write()
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
         call consfile%add_column(fs%Qint(8),'Vapor Mass')
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
         if (time%n.le.20) then
            time%dtold=time%dt
            time%dt=dt_init
         else
            call time%adjust_dt()
         end if
         call time%increment()

         ! Remember old state
         call fs%store_old()
         call debug_probe('00-start-of-step-Qold')

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call debug_probe('01-post-lincomb-raw-advection')
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call debug_probe('02-post-avgdown-fill')
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         call debug_probe('03-post-build_plic(merge_Q+clean_Q)')
         ! Relax and clean up
         call fs%apply_relax(dt=0.5_WP*time%dt,time=time%tmid)
         call debug_probe('04-post-apply_relax')
         call fs%clean_Q()
         call debug_probe('05-post-clean_Q-2nd-call')
         call fs%get_primitive(Q=fs%Q)
         call debug_probe('06-post-get_primitive')
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         call debug_probe('07-post-build_subVF')
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         call debug_probe('08-post-add_phasic_pressure')
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call debug_probe('09-post-final-avgdown-fill')
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         call debug_probe('10-END-of-stage1')

         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         ! Rebuild PLIC
         call fs%build_plic(time=time%t)
         ! Relax and clean up
         call fs%apply_relax(dt=time%dt,time=time%t)
         call fs%clean_Q()
         call fs%get_primitive(Q=fs%Q)
         ! Rebuild sub-cell VF
         call fs%build_subVF()
         ! Compute face velocities and ensure C/F consistency
         call fs%get_face_velocity(); call fs%average_down_velocity()
         ! Add pressure term
         call fs%add_phasic_pressure(scale=time%dt)
         ! Average down and fill ghosts
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         ! Get primitive variables
         call fs%get_primitive(Q=fs%Q)
         call debug_scan_worst_gas()
         ! ======================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities
         call get_viscosities()

         ! Add SGS models
         call fs%add_viscartif(dt=time%dt,Cvisc=1e-2_WP)
         call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         ! if (viz_evt%occurs()) call viz%write(time=time%t)
         call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/'//trim(case_name)//'/'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
         call fs%get_info()
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()

      end do

      ! Save the final checkpoint
      save_final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/'//trim(case_name)//'/'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block save_final_checkpoint

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
   end subroutine simulation_final

end module simulation
