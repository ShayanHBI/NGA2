!> AMR compressible venturi cavitation case
!> Case B (cavitation-induced choked flow) of Liu, Wang & Liu (2026), Int. J. Multiphase
!> Flow 195:105524, geometry per Long et al. (2017), Int. J. Multiphase Flow 89:290-298.
!>
!> True 3D round duct: the converging-diverging wall is an internal immersed-boundary
!> solid (polygon-defined profile in the (x,r) half-plane, r=sqrt(y^2+z^2), revolved
!> about the x-axis), rasterized into a fluid volume fraction field VFib (1=fluid,
!> 0=solid) and enforced by momentum direct-forcing every RK substep -- the same
!> technique as examples/amrsimplex, adapted to amrmpcomp's conserved-variable Q array
!> (only the momentum components Q(5:7) are forced; Q(1:4) and Q(8) -- phasic mass and
!> energy -- are left alone so density is never multiplied by a zero VFib).
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
   use nasg_class,             only: nasg
   use igmix_class,            only: igmix
   use relax_igmix_nasg_class, only: relax_igmix_nasg
   use polygon_class,          only: polygon
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
   type(event) :: viz_evt
   type(amrviz) :: viz

   !> Regrid parameters
   type(event) :: regrid_evt

   !> Restart parameters
   type(amrio) :: io
   type(event) :: save_evt
   character(len=str_medium) :: restart_dir
   logical :: restarted
   real(WP) :: restart_time
   integer  :: restart_step

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,rescfile
   !> Relaxation-model census (relax_model%acc reduced across ranks for the rescue monitor)
   real(WP) :: diss_n=0.0_WP,diss_m=0.0_WP
   real(WP) :: quad_n=0.0_WP,swap_n=0.0_WP,flr_n=0.0_WP,flr_e=0.0_WP,stuck_n=0.0_WP

   !> Materials
   type(nasg),  target :: water
   type(igmix), target :: gas

   !> Relaxation model
   type(relax_igmix_nasg), target :: relax_model
   character(len=str_medium) :: relaxation_type,case_name
   real(WP) :: p_cav_cof                           !< Cavitation pressure threshold coefficient
   real(WP) :: Tctol                               !< Condensation temperature tolerance

   !> Venturi geometry (Case B, following Long et al. 2017 actual dimensions)
   real(WP) :: Dth                                 !< Throat diameter [m], read from input
   real(WP) :: Rth,Rin,Rout                        !< Throat/inlet/outlet "radius" (half-height)
   real(WP) :: Lstraight,Lconv,Ldiv,Lthroat        !< Axial segment lengths
   real(WP) :: x0,x1,x2,x3,x4,x5                   !< Axial breakpoints
   real(WP), parameter :: conv_half_angle=22.5_WP  !< Convergent half-angle [deg]
   real(WP), parameter :: div_half_angle =6.0_WP   !< Divergent half-angle [deg]

   !> Immersed-boundary wall: polygon-defined flow-passage profile + fluid VF mask
   type(polygon) :: poly
   type(amrdata), target :: VFib

   !> Flow parameters
   real(WP) :: T0                                  !< Reference (isothermal) temperature [K]
   real(WP) :: p_in,p_out,p_r                       !< Inlet/outlet reservoir pressures and ratio
   real(WP) :: rhoL_in,eL_in, rhoL_out,eL_out       !< Reservoir liquid states (for the BC)
   real(WP) :: muL,muG                              !< Dynamic viscosities

   !> Wall BC type ('slip' or 'noslip')
   character(len=str_medium) :: wall_bc_type

   !> Tagging parameters
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)
   real(WP) :: VF_tag=huge(1.0_WP)

   !> Throat mass flow rate diagnostic (module-level: monitor%add_column stores a pointer
   !> to this, so it must be a persistent variable, not a subroutine-local one)
   real(WP) :: mdot_throat=0.0_WP

   !> Pre/post-projection-correction divergence diagnostics (module-level: fs%divmax itself
   !> gets overwritten by every get_div() call, so the monitor's pointer to it would only
   !> ever show whichever get_div() ran last -- these snapshot each value right after its
   !> own get_div() call so the two columns stay distinct)
   real(WP) :: divmax_precorr=0.0_WP
   real(WP) :: divmax_postcorr=0.0_WP

   !> Relative pressure-solver residual (res/bnorm) -- module-level for the same reason
   real(WP) :: pressure_residual_rel=0.0_WP

   !> Time stepping
   real(WP) :: dt_init

contains

   !> Levelset for the venturi flow passage: positive outside (solid), negative inside (fluid)
   !> per this polygon_class implementation's winding-number sign convention (verified
   !> standalone: interior points return negative distance for this vertex ordering).
   function venturi_wall_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=-poly%get_distance([xyz(1),sqrt(xyz(2)**2+xyz(3)**2)])
   end function venturi_wall_levelset

   !> Initialize VFib (1=fluid,0=solid) from the venturi wall levelset
   subroutine init_VFib(data,lvl,time,ba,dm)
      use mms_geom, only: initialize_volume_moments
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
      class(amrdata), intent(inout) :: data
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVFib
      real(WP), dimension(3) :: BL,BG
      real(WP) :: dx,dy,dz
      integer :: i,j,k
      real(WP), parameter :: VFlo=1.0e-12_WP
      integer, parameter :: nref=3
      dx=data%amr%dx(lvl); dy=data%amr%dy(lvl); dz=data%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         bx=mfi%growntilebox(data%ng)
         pVFib=>data%mf(lvl)%dataptr(mfi)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            call initialize_volume_moments(lo=[data%amr%xlo+real(i  ,WP)*dx,data%amr%ylo+real(j  ,WP)*dy,data%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[data%amr%xlo+real(i+1,WP)*dx,data%amr%ylo+real(j+1,WP)*dy,data%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=venturi_wall_levelset,time=time,level=nref,VFlo=VFlo,VF=pVFib(i,j,k,1),BL=BL,BG=BG)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine init_VFib

   !> Post-regrid dispatcher for automatic VFib refilling
   subroutine vfib_postregrid(ctx,lbase,time)
      use iso_c_binding, only: c_ptr,c_f_pointer
      type(c_ptr), intent(in) :: ctx
      integer, intent(in) :: lbase
      real(WP), intent(in) :: time
      type(amrdata), pointer :: this
      call c_f_pointer(ctx,this)
      call this%fill(time=time,lbase=lbase)
   end subroutine vfib_postregrid

   !> Zero the momentum components (Q5:7) and face velocities in/near solid cells, scaled
   !> by the local fluid fraction VFib -- same technique as amrsimplex's apply_ib_forcing,
   !> restricted to momentum so phasic mass/energy (Q1:4,Q8) are never multiplied by zero.
   subroutine apply_ib_forcing()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pU,pV,pW,pVFib
      integer :: i,j,k,lvl
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pQ   =>fs%Q%mf(lvl)%dataptr(mfi)
            pU   =>fs%U%mf(lvl)%dataptr(mfi)
            pV   =>fs%V%mf(lvl)%dataptr(mfi)
            pW   =>fs%W%mf(lvl)%dataptr(mfi)
            pVFib=>VFib%mf(lvl)%dataptr(mfi)
            ! Force cell-centered momentum
            bx=mfi%tilebox()
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pQ(i,j,k,5:7)=pVFib(i,j,k,1)*pQ(i,j,k,5:7)
            end do; end do; end do
            ! Force staggered face velocities
            bx=mfi%nodaltilebox(1)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pU(i,j,k,1)=0.5_WP*sum(pVFib(i-1:i,j,k,1))*pU(i,j,k,1)
            end do; end do; end do
            bx=mfi%nodaltilebox(2)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pV(i,j,k,1)=0.5_WP*sum(pVFib(i,j-1:j,k,1))*pV(i,j,k,1)
            end do; end do; end do
            bx=mfi%nodaltilebox(3)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pW(i,j,k,1)=0.5_WP*sum(pVFib(i,j,k-1:k,1))*pW(i,j,k,1)
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine apply_ib_forcing

   !> Compute viscosity: constant molecular values (ambient, near-isothermal), VF-blended
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pVisc,pBeta,pDiffL,pDiffG
      real(WP), parameter :: myeps=1.0e-15_WP
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pVF   =>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc =>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta =>fs%beta%mf(lvl)%dataptr(mfi)
            pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Harmonic averaging of the two molecular viscosities
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(muL,myeps)+(1.0_WP-pVF(i,j,k,1))/max(muG,myeps))
               pBeta(i,j,k,1)=0.0_WP
               ! Near-isothermal case: heat diffusion is not the controlling physics here
               pDiffL(i,j,k,1)=0.0_WP
               pDiffG(i,j,k,1)=0.0_WP
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback: quiescent liquid everywhere, pressure linearly ramped from the
   !> inlet to outlet reservoir value along x (a simple, non-equilibrium starting guess)
   subroutine venturi_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(amrex_boxarray), intent(in) :: ba
      type(amrex_distromap), intent(in) :: dm
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVF,pCL,pCG
      real(WP) :: dx,dy,dz,x_cc,p_loc,rho_loc,e_loc
      integer :: i,j,k
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      call amrex_mfiter_build(mfi,ba,dm,tiling=.false.)
      do while (mfi%next())
         pQ =>solver%Q%mf(lvl)%dataptr(mfi)
         pVF=>solver%VF%mf(lvl)%dataptr(mfi)
         if (lvl.eq.solver%amr%maxlvl) then
            pCL=>solver%CL%dataptr(mfi)
            pCG=>solver%CG%dataptr(mfi)
         end if
         bx=mfi%growntilebox(solver%nover)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Pure liquid everywhere (fs%VF is the liquid/vapor VOF, unrelated to VFib);
            ! cavitation nucleates entirely through the pTg relaxation model as pressure
            ! drops toward p_sat
            pVF(i,j,k,1)=1.0_WP
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,:)=pCL(i,j,k,:)
            end if
            ! Linear pressure ramp from inlet to outlet reservoir value
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            p_loc=p_in+(p_out-p_in)*min(max((x_cc-x0)/(x5-x0),0.0_WP),1.0_WP)
            rho_loc=water%get_rho_from_p_T(p=p_loc,T=T0,y=[1.0_WP])
            e_loc  =water%get_e_from_p_rho(p=p_loc,rho=rho_loc,y=[1.0_WP])
            pQ(i,j,k,1)=rho_loc
            pQ(i,j,k,2)=0.0_WP
            pQ(i,j,k,3)=rho_loc*e_loc
            pQ(i,j,k,4)=0.0_WP
            pQ(i,j,k,5)=0.0_WP
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
            pQ(i,j,k,8)=0.0_WP
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine venturi_init

   !> BC: fixed-pressure liquid reservoir at x-lo (inlet) and x-hi (outlet), with
   !> zero-gradient (self-adjusting) streamwise velocity -- matches the paper's own
   !> OpenFOAM Table 6 (velocity zeroGradient, pressure fixedValue at both ends)
   subroutine venturi_dirichlet_bc(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: rho_loc,e_loc
      integer :: i,j,k
      select case (face)
       case (1); rho_loc=rhoL_in;  e_loc=eL_in   ! X-LOW: inlet reservoir
       case (2); rho_loc=rhoL_out; e_loc=eL_out  ! X-HIGH: outlet reservoir
       case default; return
      end select
      select case (comp)
       case ('Q')  ! Cell-centered: pure liquid at the reservoir state, zero transverse momentum
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            p(i,j,k,1)=rho_loc; p(i,j,k,2)=0.0_WP
            p(i,j,k,3)=rho_loc*e_loc; p(i,j,k,4)=0.0_WP
            p(i,j,k,6)=0.0_WP; p(i,j,k,7)=0.0_WP; p(i,j,k,8)=0.0_WP
         end do; end do; end do
       case ('V','W')  ! Staggered V,W = 0 (no transverse inflow)
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            p(i,j,k,1)=0.0_WP
         end do; end do; end do
      end select
   end subroutine venturi_dirichlet_bc

   !> Tagger: SGS Reynolds number, density/pressure errors, Ducros sensor (ported from
   !> amrcomp_impact's my_tagger; the SGS Reynolds criterion is adapted to use the local
   !> dimensional viscosity/density fields here instead of a single non-dimensional
   !> Reynolds-number input), plus IB-wall proximity and cavitating-cell criteria specific
   !> to this case
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrtag,           only: lap_error,grd_error
      use amrmpcomp_class,  only: VFlo,VFhi
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pPL,pVF,pVFib,pUVW,pC,pVisc
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,Re
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      integer :: i,j,k
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Peps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      dx=solver%amr%dx(lvl); dxi=1.0_WP/dx; dxi2=1.0_WP/dx**2
      dy=solver%amr%dy(lvl); dyi=1.0_WP/dy; dyi2=1.0_WP/dy**2
      dz=solver%amr%dz(lvl); dzi=1.0_WP/dz; dzi2=1.0_WP/dz**2
      delta=solver%amr%min_meshsize(lvl); delta2=delta**2
      tags=tags_ptr
      call solver%amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         tagarr=>tags%dataPtr(mfi)
         pQ  =>solver%Q%mf(lvl)%dataptr(mfi)
         pPL =>solver%PL%mf(lvl)%dataptr(mfi)
         pVF =>solver%VF%mf(lvl)%dataptr(mfi)
         pVFib=>VFib%mf(lvl)%dataptr(mfi)
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         pC  =>solver%C%mf(lvl)%dataptr(mfi)
         pVisc=>solver%visc%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            ! Tag cut cells at the immersed wall
            if (pVFib(i,j,k,1).gt.0.0_WP.and.pVFib(i,j,k,1).lt.1.0_WP) tagarr(i,j,k,1)=SETtag

            ! Tag cavitating (mixed liquid/vapor) cells
            if (pVF(i,j,k,1).gt.VFlo.and.pVF(i,j,k,1).lt.VFhi) then
               if (VF_tag.lt.huge(1.0_WP)) tagarr(i,j,k,1)=SETtag
            end if

            ! Mixture density laplacian error
            rho_cc=sum(pQ(i  ,j,  k,  1:2))
            rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
            rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
            rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
            if (lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps).gt.Rho_tag) tagarr(i,j,k,1)=SETtag

            ! Liquid pressure gradient
            if (pVF(i,j,k,1).gt.0.0_WP) then
               if (grd_error(pPL(i,j,k,1),pPL(i-1,j,k,1),pPL(i+1,j,k,1),pPL(i,j-1,k,1),pPL(i,j+1,k,1),pPL(i,j,k-1,1),pPL(i,j,k+1,1),Peps).gt.P_tag) tagarr(i,j,k,1)=SETtag
            end if

            ! SGS cell Reynolds number: Re=u_sgs*delta*rho/mu using the local dimensional
            ! density and viscosity (amrcomp_impact instead used u_sgs*delta*Reynolds,
            ! since its viscosity was set via a single non-dimensional Reynolds-number input)
            lapU=(pUVW(i+1,j,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i-1,j,k,1))*dxi2+(pUVW(i,j+1,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j-1,k,1))*dyi2+(pUVW(i,j,k+1,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j,k-1,1))*dzi2
            lapV=(pUVW(i+1,j,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i-1,j,k,2))*dxi2+(pUVW(i,j+1,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j-1,k,2))*dyi2+(pUVW(i,j,k+1,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j,k-1,2))*dzi2
            lapW=(pUVW(i+1,j,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i-1,j,k,3))*dxi2+(pUVW(i,j+1,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j-1,k,3))*dyi2+(pUVW(i,j,k+1,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j,k-1,3))*dzi2
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            Re=rho_cc*u_sgs*delta/max(pVisc(i,j,k,1),tiny(1.0_WP))
            if (Re.gt.Re_tag) tagarr(i,j,k,1)=SETtag

            ! Ducros compression switch
            divu =0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))+0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))+0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            vortx=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))-0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            vorty=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))-0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            vortz=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))-0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            vort=sqrt(vortx**2+vorty**2+vortz**2)
            Deps=(Cduc*pC(i,j,k,1)/delta)**2
            Ducros=divu**2/max(divu**2+vort**2+Deps,tiny(1.0_WP))
            if (divu.lt.0.0_WP.and.Ducros.gt.Ducros_tag) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> Mass flow rate through the throat centerplane (x=0.5*(x2+x3)), weighted by VFib so
   !> only the open flow passage contributes -- for comparison with the paper's Fig. 20
   subroutine get_mdot(mdot)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      use mpi_f08,          only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
      use parallel,         only: MPI_REAL_WP
      real(WP), intent(out) :: mdot
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pVFib
      real(WP) :: dy,dz,x_throat,x_cc,dx
      integer :: i,j,k,lvl
      x_throat=0.5_WP*(x2+x3)
      mdot=0.0_WP
      lvl=amr%clvl()
      dx=amr%dx(lvl); dy=amr%dy(lvl); dz=amr%dz(lvl)
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pQ   =>fs%Q%mf(lvl)%dataptr(mfi)
         pVFib=>VFib%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            x_cc=amr%xlo+(real(i,WP)+0.5_WP)*dx
            if (abs(x_cc-x_throat).le.0.5_WP*dx) then
               ! Streamwise momentum (=rho*u) integrated across the throat, times dy*dz
               mdot=mdot+pVFib(i,j,k,1)*pQ(i,j,k,5)*dy*dz
            end if
         end do; end do; end do
      end do
      call amr%mfiter_destroy(mfi)
      call MPI_ALLREDUCE(MPI_IN_PLACE,mdot,1,MPI_REAL_WP,MPI_SUM,amr%comm)
   end subroutine get_mdot

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Determine case name from relaxation type (read early: needed by init_eos_and_flow
      ! below for the cavitation/condensation threshold parameters)
      call param_read('Relaxation type',relaxation_type)
      case_name='venturi_'//trim(relaxation_type)

      ! Set up venturi geometry (Case B dimensions, following Long et al. 2017)
      compute_geometry: block
         real(WP), parameter :: Pi=3.14159265358979_WP
         call param_read('Throat diameter',Dth)
         Rth=0.5_WP*Dth; Rin=2.5_WP*Dth; Rout=Rin
         Lstraight=2.0_WP*Dth; Lthroat=Dth
         Lconv=(Rin-Rth)/tan(conv_half_angle*Pi/180.0_WP)
         Ldiv =(Rout-Rth)/tan(div_half_angle *Pi/180.0_WP)
         x0=0.0_WP; x1=x0+Lstraight; x2=x1+Lconv; x3=x2+Lthroat; x4=x3+Ldiv; x5=x4+Lstraight
         ! Build the flow-passage polygon (see venturi_wall_levelset for the sign convention)
         call poly%initialize(nvert=8,name='venturi')
         poly%vert(:,1)=[x0,0.0_WP]; poly%vert(:,2)=[x0,Rin]; poly%vert(:,3)=[x1,Rin]
         poly%vert(:,4)=[x2,Rth];    poly%vert(:,5)=[x3,Rth]
         poly%vert(:,6)=[x4,Rout];   poly%vert(:,7)=[x5,Rout]; poly%vert(:,8)=[x5,0.0_WP]
      end block compute_geometry

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name=trim(case_name)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=x0; amr%xhi=x5
         amr%ylo=-Rin; amr%yhi=+Rin
         amr%zlo=-Rin; amr%zhi=+Rin  ! true 3D: z spans the same extent as y (round duct, r=sqrt(y^2+z^2))
         amr%xper=.false.; amr%yper=.false.; amr%zper=.false.
         call param_read('Max level',amr%maxlvl)
         call amr%initialize()
      end block create_amrgrid

      ! Read EoS and flow parameters
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: GammaL,PinfL,bL,CvL,qL,qpL
         real(WP) :: GammaA,CvA
         real(WP) :: GammaV,CvV,qV,qpV
         call param_read('Liquid specific heat capacity ratio',GammaL)
         call param_read('Liquid specific heat capacity at constant volume',CvL)
         call param_read('Liquid stiffening pressure',PinfL)
         call param_read('Liquid co-volume',bL)
         call param_read('Liquid reference energy shift',qL)
         call param_read('Liquid reference entropy shift',qpL)
         call param_read('Vapor specific heat capacity ratio',GammaV)
         call param_read('Vapor specific heat capacity at constant volume',CvV)
         call param_read('Vapor reference energy shift',qV)
         call param_read('Vapor reference entropy shift',qpV)
         call param_read('Air specific heat capacity ratio',GammaA)
         call param_read('Air specific heat capacity at constant volume',CvA)
         call water%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=qL,qp=qpL,name='water')
         call gas%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,0.0_WP],qp=[qpV,0.0_WP], &
         &                    species_names=['vapor','air  '],name='gas')
         if (trim(relaxation_type).eq.'pTg') then
            call param_read('Cavitation pressure threshold coefficient',p_cav_cof)
            call param_read('Condensation temperature tolerance',Tctol)
            relax_model%p_cav=-p_cav_cof*PinfL
            relax_model%Tctol=Tctol
         end if
         ! Reservoir states
         call param_read('Liquid temperature',T0)
         call param_read('Inlet pressure',p_in)
         call param_read('Pressure ratio',p_r)
         p_out=p_r*p_in
         rhoL_in =water%get_rho_from_p_T(p=p_in, T=T0,y=[1.0_WP]); eL_in =water%get_e_from_p_rho(p=p_in, rho=rhoL_in, y=[1.0_WP])
         rhoL_out=water%get_rho_from_p_T(p=p_out,T=T0,y=[1.0_WP]); eL_out=water%get_e_from_p_rho(p=p_out,rho=rhoL_out,y=[1.0_WP])
         ! Viscosities
         call param_read('Liquid viscosity',muL)
         call param_read('Gas viscosity',muG)
         ! Log
         write(message,'("[Geometry] Dth=",es12.5," Rth=",es12.5," Rin=",es12.5," Lx=",es12.5)') Dth,Rth,Rin,x5-x0; call log(message)
         write(message,'("[Reservoirs] p_in=",es12.5," p_out=",es12.5," (p_r=",es12.5,")")') p_in,p_out,p_r; call log(message)
         call water%print(); call gas%print()
      end block init_eos_and_flow

      ! Handle restart/saves here
      handle_restart: block
         call io%initialize(amr=amr,nfiles=1)
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         if (restarted) restart_dir='restart/'//trim(case_name)//'/'//trim(adjustl(restart_dir))
         if (restarted) call io%read_header(dirname=trim(restart_dir),time=restart_time,step=restart_step)
      end block handle_restart

      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=amr%amRoot)
         call param_read('Max time',time%tmax)
         call param_read('Max dt',time%dtmax)
         call param_read('Initial dt',dt_init,default=time%dtmax)
         call param_read('Max CFL',time%cflmax)
         time%dt=time%dtmax
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module,     only: amrex_bc_ext_dir,amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrmpcomp_class,      only: BC_REFLECT,BC_LIQ
         use amrdata_class,        only: interp_face_lin
         use messager,             only: die
         use relax_igmix_sg_class, only: Prelax,PTrelax,PTgrelax,PThybrid
         use amrmg_class,          only: amrmg_outer_pcg_mlmg
         fs%liq=>water; fs%gas=>gas
         call param_read('Use projection',fs%use_projection,default=.true.)
         call fs%initialize(amr=amr,name=trim(case_name))
         if (fs%use_projection) then
            ! fs%psolver%outer_solver=amrmg_outer_pcg_mlmg
            ! fs%psolver%max_iter=300
            fs%psolver%tol_rel=1.0e-5_WP
            fs%psolver%verbose=2
         end if
         call param_read('Surface tension',fs%sigma,default=0.07_WP)
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Relaxation model: species 1=vapor (transported), 2=air (carrier)
         call relax_model%initialize(liq=water,gas=gas,indV=1,indA=2)
         select case (trim(relaxation_type))
         case ('p');        relax_model%model=Prelax
         case ('pT');       relax_model%model=PTrelax
         case ('pTg');      relax_model%model=PTgrelax
         case ('pThybrid'); relax_model%model=PThybrid
         case default; call die('[simulation_init] Relaxation type must be p, pT, pThybrid, or pTg')
         end select
         relax_model%RHOGmin=0.0_WP
         relax_model%vol=amr%cell_vol(amr%maxlvl)
         fs%merge_sick=100.0_WP
         relax_model%diss_P=1.0e30_WP
         fs%Pmin_liq=-0.98_WP*water%pinf
         fs%Tmin_liq=250.0_WP
         fs%Pmin_gas=1.0e-4_WP
         fs%Tmin_gas=250.0_WP
         relax_model%Pmin_liq=fs%Pmin_liq; relax_model%Tmin_liq=fs%Tmin_liq
         relax_model%Pmin_gas=fs%Pmin_gas; relax_model%Tmin_gas=fs%Tmin_gas
         fs%relax=>relax_model
         fs%user_init=>venturi_init

         ! X-LOW (inlet) and X-HIGH (outlet): both are liquid reservoirs
         fs%lo_bc(1)=BC_LIQ; fs%hi_bc(1)=BC_LIQ
         fs%Q%lo_bc(1,1:4)=amrex_bc_ext_dir;    fs%Q%hi_bc(1,1:4)=amrex_bc_ext_dir
         fs%Q%lo_bc(1,5)  =amrex_bc_foextrap;   fs%Q%hi_bc(1,5)  =amrex_bc_foextrap
         fs%Q%lo_bc(1,6:8)=amrex_bc_ext_dir;    fs%Q%hi_bc(1,6:8)=amrex_bc_ext_dir
         fs%U%lo_bc(1,1)  =amrex_bc_foextrap;   fs%U%hi_bc(1,1)  =amrex_bc_foextrap
         fs%V%lo_bc(1,1)  =amrex_bc_ext_dir;    fs%V%hi_bc(1,1)  =amrex_bc_ext_dir
         fs%W%lo_bc(1,1)  =amrex_bc_ext_dir;    fs%W%hi_bc(1,1)  =amrex_bc_ext_dir
         fs%user_bc=>venturi_dirichlet_bc

         ! Y-LO/Y-HI and Z-LO/Z-HI: reflective wall. The domain's y and z extent both equal
         ! Rin=Rout exactly, so these domain-edge BCs ARE the physical wall along the
         ! inlet/outlet straight sections (where the profile radius equals Rin); along the
         ! converging/throat/diverging sections the domain-edge cells there are already
         ! masked solid by VFib, so these BCs are redundant-but-harmless there (IB forcing
         ! zeroes their momentum anyway)
         call param_read('Wall BC',wall_bc_type,default='noslip')

         fs%lo_bc(2)=BC_REFLECT; fs%hi_bc(2)=BC_REFLECT
         fs%Q%lo_bc(2,:)=amrex_bc_foextrap;    fs%Q%hi_bc(2,:)=amrex_bc_foextrap
         fs%Q%lo_bc(2,6)=amrex_bc_reflect_odd; fs%Q%hi_bc(2,6)=amrex_bc_reflect_odd
         fs%V%lo_bc(2,:)=amrex_bc_reflect_odd; fs%V%hi_bc(2,:)=amrex_bc_reflect_odd
         select case (trim(wall_bc_type))
         case ('noslip')
            fs%Q%lo_bc(2,5)=amrex_bc_reflect_odd; fs%Q%hi_bc(2,5)=amrex_bc_reflect_odd
            fs%Q%lo_bc(2,7)=amrex_bc_reflect_odd; fs%Q%hi_bc(2,7)=amrex_bc_reflect_odd
            fs%U%lo_bc(2,:)=amrex_bc_reflect_odd; fs%U%hi_bc(2,:)=amrex_bc_reflect_odd
            fs%W%lo_bc(2,:)=amrex_bc_reflect_odd; fs%W%hi_bc(2,:)=amrex_bc_reflect_odd
         case ('slip')
            fs%Q%lo_bc(2,5)=amrex_bc_foextrap;    fs%Q%hi_bc(2,5)=amrex_bc_foextrap
            fs%Q%lo_bc(2,7)=amrex_bc_foextrap;    fs%Q%hi_bc(2,7)=amrex_bc_foextrap
            fs%U%lo_bc(2,:)=amrex_bc_foextrap;    fs%U%hi_bc(2,:)=amrex_bc_foextrap
            fs%W%lo_bc(2,:)=amrex_bc_foextrap;    fs%W%hi_bc(2,:)=amrex_bc_foextrap
         case default
            call die('[simulation_init] Wall BC must be slip or noslip')
         end select

         fs%lo_bc(3)=BC_REFLECT; fs%hi_bc(3)=BC_REFLECT
         fs%Q%lo_bc(3,:)=amrex_bc_foextrap;    fs%Q%hi_bc(3,:)=amrex_bc_foextrap
         fs%Q%lo_bc(3,7)=amrex_bc_reflect_odd; fs%Q%hi_bc(3,7)=amrex_bc_reflect_odd
         fs%W%lo_bc(3,:)=amrex_bc_reflect_odd; fs%W%hi_bc(3,:)=amrex_bc_reflect_odd
         select case (trim(wall_bc_type))
         case ('noslip')
            fs%Q%lo_bc(3,5)=amrex_bc_reflect_odd; fs%Q%hi_bc(3,5)=amrex_bc_reflect_odd
            fs%Q%lo_bc(3,6)=amrex_bc_reflect_odd; fs%Q%hi_bc(3,6)=amrex_bc_reflect_odd
            fs%U%lo_bc(3,:)=amrex_bc_reflect_odd; fs%U%hi_bc(3,:)=amrex_bc_reflect_odd
            fs%V%lo_bc(3,:)=amrex_bc_reflect_odd; fs%V%hi_bc(3,:)=amrex_bc_reflect_odd
         case ('slip')
            fs%Q%lo_bc(3,5)=amrex_bc_foextrap;    fs%Q%hi_bc(3,5)=amrex_bc_foextrap
            fs%Q%lo_bc(3,6)=amrex_bc_foextrap;    fs%Q%hi_bc(3,6)=amrex_bc_foextrap
            fs%U%lo_bc(3,:)=amrex_bc_foextrap;    fs%U%hi_bc(3,:)=amrex_bc_foextrap
            fs%V%lo_bc(3,:)=amrex_bc_foextrap;    fs%V%hi_bc(3,:)=amrex_bc_foextrap
         case default
            call die('[simulation_init] Wall BC must be slip or noslip')
         end select
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1    ,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1    ,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

      ! Create the immersed-boundary fluid-fraction field for the venturi wall
      create_VFib: block
         use amrdata_class,    only: interp_const
         use amrex_amr_module, only: amrex_bc_foextrap
         use iso_c_binding,    only: c_loc
         call VFib%initialize(amr,name='VFib',ncomp=1,ng=fs%nover,interp=interp_const); call VFib%register()
         call amr%add_postregrid(vfib_postregrid,c_loc(VFib))
         VFib%user_init=>init_VFib
         VFib%lo_bc(1,1)=amrex_bc_foextrap; VFib%hi_bc(1,1)=amrex_bc_foextrap
         VFib%lo_bc(2,1)=amrex_bc_foextrap; VFib%hi_bc(2,1)=amrex_bc_foextrap
      end block create_VFib

      ! Initialize regridding
      init_regridding: block
         amr%lb_strat=1
         regrid_evt=event(time=time,name='Regrid')
         call param_read('Regrid nsteps',regrid_evt%nper)
         fs%user_tagging=>my_tagger
         call param_read('Tag Reynolds value',Re_tag)
         call param_read('Tag density error' ,Rho_tag)
         call param_read('Tag pressure error',P_tag)
         call param_read('Tag Ducros value',Ducros_tag)
         call param_read('Tag cavitating cells',VF_tag,default=1.0_WP)
         if (restarted) then
            call amr%init_from_checkpoint(dirname=trim(restart_dir),time=time%t)
            call fs%restore_checkpoint(io=io,dirname=trim(restart_dir),time=time%t)
         else
            call amr%init_from_scratch(time=time%t)
            call fs%build_plic(time%t)
            call fs%build_subVF()
            call fs%get_primitive(Q=fs%Q)
            call fs%get_face_velocity()
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         end if
         call apply_ib_forcing()
         call get_viscosities()
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
      end block init_regridding

      ! Initialize checkpoint save event
      init_checkpoint: block
         save_evt=event(time=time,name='Checkpoint')
         call param_read('Checkpoint period',save_evt%tper,default=-1.0_WP)
         call fs%register_checkpoint(io)
         call io%add_scalar(name='dt',value=time%dt)
      end block init_checkpoint

      ! Initialize visualization
      create_viz: block
         call viz%initialize(amr,trim(case_name),use_hdf5=.false.)
         call viz%add_scalar(fs%VF,1,'VF')
         call viz%add_scalar(VFib,1,'VFib')
         call viz%add_scalar(fs%RHOL,1,'RHOL')
         call viz%add_scalar(fs%RHOG,1,'RHOG')
         call viz%add_scalar(fs%PL,1,'PL')
         call viz%add_scalar(fs%PG,1,'PG')
         call viz%add_scalar(fs%TL,1,'TL')
         call viz%add_scalar(fs%TG,1,'TG')
         call viz%add_scalar(fs%UVW,1,'U')
         call viz%add_scalar(fs%UVW,2,'V')
         call viz%add_scalar(fs%UVW,3,'W')
         call viz%add_scalar(Umag,1,'Umag')
         call viz%add_scalar(Mach,1,'Mach')
         call viz%add_scalar(fs%Yg,1,'Yv')
         call viz%add_surfmesh(fs%smesh,'plic')
         viz_evt=event(time=time,name='Visualization output')
         call param_read('Output period',viz_evt%tper)
         if (viz_evt%occurs()) call viz%write(time=time%t)
      end block create_viz

      ! Create monitors
      create_monitors: block
         call fs%get_info()
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call get_mdot(mdot_throat)
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(time%cfl,'CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%RHOLmin,'rhoLmin')
         call mfile%add_column(fs%RHOLmax,'rhoLmax')
         call mfile%add_column(fs%PLmin,'PLmin')
         call mfile%add_column(fs%PLmax,'PLmax')
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
         call mfile%add_column(fs%dPmax,'dPmax')
         call mfile%add_column(mdot_throat,'mdot')
         call mfile%add_column(fs%psolver%res,'P res')
         call mfile%add_column(pressure_residual_rel,'P res rel')
         call mfile%add_column(fs%psolver%niter,'Pressure iterations')
         call mfile%add_column(divmax_precorr,'Divergence')
         call mfile%add_column(divmax_postcorr,'Divergence (post)')
         call mfile%write()
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
         call cflfile%add_column(fs%CFLst ,'CFLst' )
         call cflfile%write()
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%VFint,'VFint')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
         call consfile%add_column(fs%Qint(6),'V Momentum')
         call consfile%add_column(fs%Qint(7),'W Momentum')
         call consfile%add_column(fs%rhoKint,'Kinetic energy')
         call consfile%write()
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
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
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
         rescfile=monitor(amRoot=amr%amRoot,name='rescue')
         call rescfile%add_column(time%n,'Timestep')
         call rescfile%add_column(time%t,'Time')
         call rescfile%add_column(fs%resc_nl,'LiqResc n')
         call rescfile%add_column(fs%resc_ml,'LiqResc dm')
         call rescfile%add_column(fs%resc_el,'LiqResc dE')
         call rescfile%add_column(fs%resc_ng,'GasResc n')
         call rescfile%add_column(fs%resc_mg,'GasResc dm')
         call rescfile%add_column(fs%resc_eg,'GasResc dE')
         call rescfile%add_column(diss_n,'Diss n')
         call rescfile%add_column(diss_m,'Diss dm')
         call rescfile%add_column(quad_n,'Quad n')
         call rescfile%add_column(swap_n,'Swap n')
         call rescfile%add_column(flr_n,'Floor n')
         call rescfile%add_column(flr_e,'Floor dE')
         call rescfile%add_column(stuck_n,'Stuck n')
         call rescfile%add_column(fs%pool_n,'Pool n')
         call rescfile%write()
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      do while (.not.time%done())

         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         if (time%n.lt.50) then
            time%dtold=time%dt
            time%dt=dt_init
         else
            call time%adjust_dt()
         end if
         call time%increment()

         call fs%store_old()

         ! ======================= RK2 Stage 1 =======================
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%build_plic(time=time%t)
         call fs%apply_relax(dt=0.5_WP*time%dt,time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         call apply_ib_forcing()
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         if (fs%use_projection) then
            call fs%get_div(); divmax_precorr=fs%divmax
            call fs%div%mult(val=1.0_WP/(0.5_WP*time%dt))
            ! TEMP DEBUG: independent raw max(abs(div)) check right before handing rhs to psolver
            debug_rhscheck1: block
               use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
               use parallel, only: MPI_REAL_WP
               integer :: lvl2,i2,j2,k2,ierr2
               type(amrex_mfiter) :: mfi2
               type(amrex_box) :: bx2
               real(WP), dimension(:,:,:,:), contiguous, pointer :: pDiv2
               real(WP) :: mydivmax
               mydivmax=0.0_WP
               do lvl2=0,amr%clvl()
                  call amr%mfiter_build(lvl2,mfi2)
                  do while (mfi2%next())
                     pDiv2=>fs%div%mf(lvl2)%dataptr(mfi2)
                     bx2=mfi2%tilebox()
                     do k2=bx2%lo(3),bx2%hi(3); do j2=bx2%lo(2),bx2%hi(2); do i2=bx2%lo(1),bx2%hi(1)
                        mydivmax=max(mydivmax,abs(pDiv2(i2,j2,k2,1)))
                     end do; end do; end do
                  end do
                  call amr%mfiter_destroy(mfi2)
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,mydivmax,1,MPI_REAL_WP,MPI_MAX,amr%comm,ierr2)
               if (amr%amRoot) print '("[DEBUGRHS1] n=",i0," raw maxabs(fs%div before solve)=",es14.6)', time%n, mydivmax
               flush(6)
            end block debug_rhscheck1
            call fs%prepare_psolver(dt=0.5_WP*time%dt)
            call fs%psolver%solve(rhs=fs%div)
            pressure_residual_rel=fs%psolver%res/max(fs%psolver%bnorm,epsilon(1.0_WP))
            call fs%add_pressure_correction(scale=0.5_WP*time%dt)
            call apply_ib_forcing()
            call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
            call fs%get_div(); divmax_postcorr=fs%divmax
         end if
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2 =======================
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%build_plic(time=time%t)
         call fs%apply_relax(dt=time%dt,time=time%t)
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=time%dt)
         call fs%add_surface_tension(scale=time%dt)
         call apply_ib_forcing()
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         if (fs%use_projection) then
            call fs%get_div(); divmax_precorr=fs%divmax
            call fs%div%mult(val=1.0_WP/time%dt)
            ! TEMP DEBUG: independent raw max(abs(div)) check right before handing rhs to psolver
            debug_rhscheck2: block
               use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_mfiter_build,amrex_mfiter_destroy
               use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_IN_PLACE
               use parallel, only: MPI_REAL_WP
               integer :: lvl2,i2,j2,k2,ierr2
               type(amrex_mfiter) :: mfi2
               type(amrex_box) :: bx2
               real(WP), dimension(:,:,:,:), contiguous, pointer :: pDiv2
               real(WP) :: mydivmax
               mydivmax=0.0_WP
               do lvl2=0,amr%clvl()
                  call amr%mfiter_build(lvl2,mfi2)
                  do while (mfi2%next())
                     pDiv2=>fs%div%mf(lvl2)%dataptr(mfi2)
                     bx2=mfi2%tilebox()
                     do k2=bx2%lo(3),bx2%hi(3); do j2=bx2%lo(2),bx2%hi(2); do i2=bx2%lo(1),bx2%hi(1)
                        mydivmax=max(mydivmax,abs(pDiv2(i2,j2,k2,1)))
                     end do; end do; end do
                  end do
                  call amr%mfiter_destroy(mfi2)
               end do
               call MPI_ALLREDUCE(MPI_IN_PLACE,mydivmax,1,MPI_REAL_WP,MPI_MAX,amr%comm,ierr2)
               if (amr%amRoot) print '("[DEBUGRHS2] n=",i0," raw maxabs(fs%div before solve)=",es14.6)', time%n, mydivmax
               flush(6)
            end block debug_rhscheck2
            call fs%prepare_psolver(dt=time%dt)
            call fs%psolver%solve(rhs=fs%div)
            pressure_residual_rel=fs%psolver%res/max(fs%psolver%bnorm,epsilon(1.0_WP))
            call fs%add_pressure_correction(scale=time%dt)
            call apply_ib_forcing()
            call fs%Q%average_down(); call fs%Q%fill(time=time%t)
            call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
            call fs%get_div(); divmax_postcorr=fs%divmax
         end if
         call fs%get_primitive(Q=fs%Q)
         ! =============================================================

         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call apply_ib_forcing()
            call gridfile%write()
         end if

         call get_viscosities()
         call fs%add_viscartif(dt=time%dt,Cvisc=1.0e-2_WP)
         call fs%add_vreman(dt=time%dt)

         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         if (viz_evt%occurs()) call viz%write(time=time%t)

         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/'//trim(case_name)//'/'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         call fs%get_info()
         relax_census: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_IN_PLACE,MPI_SUM
            use parallel, only: MPI_REAL_WP
            real(WP), dimension(7) :: tmp
            integer :: ierr
            tmp=relax_model%acc
            call MPI_ALLREDUCE(MPI_IN_PLACE,tmp,7,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
            diss_n=tmp(1); diss_m=tmp(2); quad_n=tmp(3); swap_n=tmp(4)
            flr_n=tmp(5); flr_e=tmp(6); stuck_n=tmp(7)
         end block relax_census
         call get_mdot(mdot_throat)
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         call rescfile%write()

      end do

      save_final_checkpoint: block
         use string, only: rtoa
         call io%write(dirname='restart/'//trim(case_name)//'/'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
      end block save_final_checkpoint

   end subroutine simulation_run

   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      call time%finalize()
      call amr%finalize()
      call regrid_evt%finalize()
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      call VFib%finalize()
      call water%finalize()
      call gas%finalize()
      call viz%finalize()
      call viz_evt%finalize()
      call save_evt%finalize()
      call io%finalize()
      call mfile%finalize()
      call consfile%finalize()
      call cflfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
      call rescfile%finalize()
   end subroutine simulation_final

end module simulation
