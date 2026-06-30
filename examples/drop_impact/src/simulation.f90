module simulation
   use precision,           only: WP
   use string,              only: str_medium
   use amrgrid_class,       only: amrgrid
   use amrmpcomp_class,     only: amrmpcomp
   use amrviz_class,        only: amrviz
   use amrdata_class,       only: amrdata
   use timetracker_class,   only: timetracker
   use event_class,         only: event
   use monitor_class,       only: monitor
   use amrio_class,         only: amrio
   use nasg_class,          only: nasg
   use ig_class,            only: ig
   use igmix_class,         only: igmix
   use relax_class,         only: relax
   use relax_nasg_ig_class, only: relax_nasg_ig
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
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile

   !> EOS and relaxation
   type(nasg), target, save :: water        !< Liquid (NASG water)
   type(ig), allocatable, target, save :: eosG(:)   !< Gas-phase species EOS array (2: inert-vapor slot + air)
   type(igmix), target, save :: mixG         !< Gas mixture
   type(relax_nasg_ig), target, save :: relax_model  !< Relaxation model
   character(len=str_medium), save   :: relaxation_type
   character(len=str_medium), save   :: case_name

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
   real(WP) :: Weber                  !< Weber number

   !> EOS scalars needed in get_viscosities
   real(WP) :: GammaG,CvG,CpG        !< Gas EOS scalars
   real(WP) :: GammaL,PinfL,bL,CvL,qpL !< Liquid EOS scalars

   !> Drop disturbances
   real(WP) :: dist_amp=0.005_WP      !< Disturbance amplitude
   real(WP) :: dist_num=64.0_WP       !< Disturbance wavenumber

   !> Sutherland viscosity parameters: mu_g = (1+Suth_T)*T^Suth_n / (Re*(T+Suth_T))
   real(WP) :: Suth_n=1.5_WP          !< Sutherland exponent (1.0 for constant)
   real(WP) :: Suth_T=0.4042_WP       !< Sutherland temperature (0.0 for constant)

   !> Drop initial location
   real(WP) :: x_drop

   !> Wall BC type ('slip' or 'noslip')
   character(len=str_medium) :: wall_bc_type

   !> Sponge parameters
   real(WP) :: R_spg=3.0_WP
   real(WP) :: L_spg=1.0_WP

   !> Tagging parameters
   real(WP) :: Re_tag=huge(1.0_WP)
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)

contains

   !> Smooth Heaviside function
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock

   !> Levelset function for drop (centered at x=x_drop)
   function sphere_levelset(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G,r,theta,r_perturbed
      theta=atan2(xyz(2),xyz(1)-x_drop)
      r_perturbed=0.5_WP*(1.0_WP+dist_amp*cos(real(dist_num,WP)*theta))
      if (amr%nz.eq.1) then
         r=sqrt((xyz(1)-x_drop)**2+xyz(2)**2)
      else
         r=sqrt((xyz(1)-x_drop)**2+xyz(2)**2+xyz(3)**2)
      end if
      G=r_perturbed-r
   end function sphere_levelset

   !> Relaxation step wrapper for p
   subroutine relax_p(VF,Q,Pjump)
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      call relax_model%relax_p(VF=VF,Q=Q,Pjump=Pjump)
   end subroutine relax_p

   !> Relaxation step wrapper for p and T
   subroutine relax_pT(VF,Q,Pjump)
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      call relax_model%relax_pT(VF=VF,Q=Q,Pjump=Pjump)
   end subroutine relax_pT

   !> Relaxation step wrapper for p, T, and g
   subroutine relax_pTg(VF,Q,Pjump)
      implicit none
      real(WP), intent(inout) :: VF
      real(WP), dimension(:), intent(inout) :: Q
      real(WP), intent(in) :: Pjump
      call relax_model%relax_pTg(VF=VF,Q=Q,Pjump=Pjump)
   end subroutine relax_pTg

   !> Compute viscosity: Sutherland for gas, sponge layer in outer radial region
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pYg,pVisc,pBeta,pCond,pDiff,pRHOL,pRHOG
      real(WP) :: r_cyl,blend,nu_spg,mu_spg,mu_g,mu_l
      real(WP), parameter :: Tmax_visc=10.0_WP
      real(WP), parameter :: myeps=1.0e-15_WP
      real(WP), parameter :: max_cfl=0.5_WP
      real(WP), parameter :: Cdiff=0.1_WP
      nu_spg=max_cfl*amr%min_meshsize(amr%clvl())**2/(4.0_WP*time%dt)
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            pTG  =>fs%TG%mf(lvl)%dataptr(mfi)
            pVF  =>fs%VF%mf(lvl)%dataptr(mfi)
            pYg  =>fs%Yg%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pCond=>fs%cond%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            pRHOL=>fs%RHOL%mf(lvl)%dataptr(mfi)
            pRHOG=>fs%RHOG%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! Gas viscosity from Sutherland
               mu_g=(1.0_WP+Suth_T)*min(pTG(i,j,k,1),Tmax_visc)**Suth_n/(Reynolds*(min(pTG(i,j,k,1),Tmax_visc)+Suth_T))
               ! Liquid viscosity from ratio
               mu_l=visc_ratio*Reynolds**(-1.0_WP)
               ! Mixture viscosity (harmonic)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps))
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Phasic heat conductivities: k = cp*mu/Pr
               pCond(i,j,k,1)=diff_ratio*CpG/(Reynolds*Prandtl)   ! liquid (using Cp_G as reference, cf. develop)
               pCond(i,j,k,2)=CpG*mu_g/Prandtl                    ! gas (Sutherland)
               ! Vapor mass diffusivity: zero (no real phase-change physics in this case)
               pDiff(i,j,k,1)=0.0_WP
               ! Sponge layer in outer radial region (y-z plane)
               r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2+(amr%zlo+(real(k,WP)+0.5_WP)*amr%dz(lvl))**2)
               if (amr%nz.eq.1) r_cyl=sqrt((amr%ylo+(real(j,WP)+0.5_WP)*amr%dy(lvl))**2)
               if (r_cyl.gt.R_spg) then
                  blend=min((r_cyl-R_spg)/L_spg,1.0_WP)**2
                  mu_spg=nu_spg/(pVF(i,j,k,1)/max(pRHOL(i,j,k,1),myeps)+(1.0_WP-pVF(i,j,k,1))/max(pRHOG(i,j,k,1),myeps))
                  pVisc(i,j,k,1)=max(pVisc(i,j,k,1),blend*mu_spg)
                  pCond(i,j,k,1)=max(pCond(i,j,k,1),Cdiff*blend*mu_spg)
                  pCond(i,j,k,2)=max(pCond(i,j,k,2),Cdiff*blend*mu_spg)
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback - set Q and VF/barycenters for a drop moving into a static shock
   subroutine shockdrop_init(solver,lvl,time,ba,dm)
      use amrex_amr_module, only: amrex_boxarray,amrex_distromap,amrex_mfiter,amrex_box
      use amrex_amr_module, only: amrex_mfiter_build,amrex_mfiter_destroy
      use mms_geom,         only: initialize_volume_moments
      use amrmpcomp_class,  only: VFlo
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
      dx=solver%amr%dx(lvl); dy=solver%amr%dy(lvl); dz=solver%amr%dz(lvl)
      IEL=water%get_e_from_p_rho(p=pL1,rho=rhoL1)
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
            call initialize_volume_moments(lo=[solver%amr%xlo+real(i  ,WP)*dx,solver%amr%ylo+real(j  ,WP)*dy,solver%amr%zlo+real(k  ,WP)*dz], &
            &                              hi=[solver%amr%xlo+real(i+1,WP)*dx,solver%amr%ylo+real(j+1,WP)*dy,solver%amr%zlo+real(k+1,WP)*dz], &
            &                              levelset=sphere_levelset,time=time,level=nref,VFlo=VFlo,VF=myVF,BL=BL,BG=BG)
            pVF(i,j,k,1)=myVF
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=BL
               pCG(i,j,k,:)=BG
            end if
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            H=Hshock(x=Xs-x_cc,delta=0.5_WP*dx)
            rhoG=rhoG1+(rhoG2-rhoG1)*H
            pG  =pG1  +(pG2  -pG1  )*H
            uG  =u1   +(u2   -u1   )*H
            pQ(i,j,k,1)=(       myVF)*rhoL1
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG
            pQ(i,j,k,3)=pQ(i,j,k,1)*IEL
            pQ(i,j,k,4)=pQ(i,j,k,2)*eosG(2)%get_e_from_p_rho(p=pG,rho=rhoG)
            pQ(i,j,k,5)=(pQ(i,j,k,1)+pQ(i,j,k,2))*uG
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
            pQ(i,j,k,8)=0.0_WP   ! inert vapor slot: Yv=0
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine shockdrop_init

   !> BC routine: x-HIGH Dirichlet inflow with pre-shock state (disabled — using BC_REFLECT instead)
   subroutine shock_dirichlet(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      integer :: i,j,k
      select case (face)
       case (2)  ! X-HIGH: Dirichlet inflow with pre-shock state (gas only, no liquid)
         select case (comp)
          case ('U')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=u1
            end do; end do; end do
          case ('V','W')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
            end do; end do; end do
          case ('Q')
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               p(i,j,k,1)=0.0_WP
               p(i,j,k,2)=rhoG1
               p(i,j,k,3)=0.0_WP
               p(i,j,k,4)=rhoG1*eosG(2)%get_e_from_p_rho(p=pG1,rho=rhoG1)
               p(i,j,k,5)=rhoG1*u1
               p(i,j,k,6)=0.0_WP
               p(i,j,k,7)=0.0_WP
               p(i,j,k,8)=0.0_WP   ! inert vapor slot
            end do; end do; end do
         end select
      end select
   end subroutine shock_dirichlet

   !> Tagger based on SGS Reynolds number, density/pressure errors, Ducros sensor
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrtag,           only: lap_error,grd_error
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      type(c_ptr), intent(in) :: tags_ptr
      type(amrex_tagboxarray) :: tags
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      character(kind=c_char), dimension(:,:,:,:), contiguous, pointer :: tagarr
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pQ,pPL,pVF,pUVW,pC
      real(WP) :: dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,delta,delta2
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: lapU,lapV,lapW,u_sgs,Re
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      real(WP) :: r_cyl
      logical  :: in_zone
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
         pUVW=>solver%UVW%mf(lvl)%dataptr(mfi)
         pC  =>solver%C%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
            r_cyl=sqrt((solver%amr%ylo+(real(j,WP)+0.5_WP)*dy)**2+(solver%amr%zlo+(real(k,WP)+0.5_WP)*dz)**2)
            in_zone=(r_cyl.lt.R_spg+L_spg.or.lvl.lt.solver%amr%maxlvl-1)
            ! Mixture density laplacian error
            rho_cc=sum(pQ(i  ,j,  k,  1:2))
            rho_xp=sum(pQ(i+1,j,  k,  1:2)); rho_xm=sum(pQ(i-1,j,  k,  1:2))
            rho_yp=sum(pQ(i,  j+1,k,  1:2)); rho_ym=sum(pQ(i,  j-1,k,  1:2))
            rho_zp=sum(pQ(i,  j,  k+1,1:2)); rho_zm=sum(pQ(i,  j,  k-1,1:2))
            if (lap_error(rho_cc,rho_xm,rho_xp,rho_ym,rho_yp,rho_zm,rho_zp,Reps).gt.Rho_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
            ! Liquid pressure gradient
            if (pVF(i,j,k,1).gt.0.0_WP) then
               if (grd_error(pPL(i,j,k,1),pPL(i-1,j,k,1),pPL(i+1,j,k,1),pPL(i,j-1,k,1),pPL(i,j+1,k,1),pPL(i,j,k-1,1),pPL(i,j,k+1,1),Peps).gt.P_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
            end if
            ! SGS cell Reynolds number
            lapU=(pUVW(i+1,j,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i-1,j,k,1))*dxi2+(pUVW(i,j+1,k,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j-1,k,1))*dyi2+(pUVW(i,j,k+1,1)-2.0_WP*pUVW(i,j,k,1)+pUVW(i,j,k-1,1))*dzi2
            lapV=(pUVW(i+1,j,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i-1,j,k,2))*dxi2+(pUVW(i,j+1,k,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j-1,k,2))*dyi2+(pUVW(i,j,k+1,2)-2.0_WP*pUVW(i,j,k,2)+pUVW(i,j,k-1,2))*dzi2
            lapW=(pUVW(i+1,j,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i-1,j,k,3))*dxi2+(pUVW(i,j+1,k,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j-1,k,3))*dyi2+(pUVW(i,j,k+1,3)-2.0_WP*pUVW(i,j,k,3)+pUVW(i,j,k-1,3))*dzi2
            u_sgs=0.2_WP*sqrt(lapU**2+lapV**2+lapW**2)*delta2
            Re=Reynolds*u_sgs*delta
            if (Re.gt.Re_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
            ! Ducros compression switch
            divu =0.5_WP*dxi*(pUVW(i+1,j,k,1)-pUVW(i-1,j,k,1))+0.5_WP*dyi*(pUVW(i,j+1,k,2)-pUVW(i,j-1,k,2))+0.5_WP*dzi*(pUVW(i,j,k+1,3)-pUVW(i,j,k-1,3))
            vortx=0.5_WP*dyi*(pUVW(i,j+1,k,3)-pUVW(i,j-1,k,3))-0.5_WP*dzi*(pUVW(i,j,k+1,2)-pUVW(i,j,k-1,2))
            vorty=0.5_WP*dzi*(pUVW(i,j,k+1,1)-pUVW(i,j,k-1,1))-0.5_WP*dxi*(pUVW(i+1,j,k,3)-pUVW(i-1,j,k,3))
            vortz=0.5_WP*dxi*(pUVW(i+1,j,k,2)-pUVW(i-1,j,k,2))-0.5_WP*dyi*(pUVW(i,j+1,k,1)-pUVW(i,j-1,k,1))
            vort=sqrt(vortx**2+vorty**2+vortz**2)
            Deps=(Cduc*pC(i,j,k,1)/delta)**2
            Ducros=divu**2/max(divu**2+vort**2+Deps,tiny(1.0_WP))
            if (divu.lt.0.0_WP.and.Ducros.gt.Ducros_tag.and.in_zone) tagarr(i,j,k,1)=SETtag
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
         real(WP) :: T_G
         ! Shock and EOS parameters
         call param_read('GammaG',GammaG)
         call param_read('GammaL',GammaL)
         call param_read('Liquid pinf',PinfL)
         call param_read('Liquid covolume',bL)
         call param_read('Liquid cv',CvL)
         call param_read('Liquid qp',qpL)
         ! Shock parameters (gas phase)
         call param_read('Gas Mach number',M2)
         call param_read('Shock location',Xs)
         ! Post-shock normalization: rhoG2=1, Deltau=1, T2=1
         rhoG2=1.0_WP
         pG2=1.0_WP/(GammaG*M2**2)
         ! Quadratic for rhoG1 from Rankine-Hugoniot (Deltau=1)
         A=2.0_WP*GammaG*pG2+(GammaG-1.0_WP)
         B=4.0_WP*GammaG*pG2+(GammaG+1.0_WP)
         C=2.0_WP*GammaG*pG2
         rhoG1=(B-sqrt(B**2-4.0_WP*A*C))/(2.0_WP*A)
         ! Shock-fixed frame velocities and pressure
         u1=1.0_WP/(1.0_WP-rhoG1)
         u2=u1-1.0_WP
         pG1=pG2-rhoG1/(1.0_WP-rhoG1)
         if (pG1.le.0.0_WP) call die('[simulation_init] Cannot achieve requested Mach number - negative pre-shock pressure')
         Ms=u1/sqrt(GammaG*pG1/rhoG1)
         ! Shift to lab frame then wall frame
         u2=1.0_WP; u1=0.0_WP
         u2=u2-1.0_WP; u1=u1-1.0_WP
         ! Drop initial location
         call param_read('Drop location',x_drop)
         ! CvG from T2=1 normalization
         CvG=pG2/(rhoG2*(GammaG-1.0_WP))
         CpG=GammaG*CvG
         ! Surface tension
         call param_read('Weber number',Weber)
         ! Pre-shock gas temperature
         T_G=pG1/(rhoG1*(GammaG-1.0_WP)*CvG)
         ! Liquid pressure (Laplace jump)
         pL1=pG1+4.0_WP/Weber
         if (amr%nz.eq.1) pL1=pG1+2.0_WP/Weber
         ! Liquid density from NASG EOS at (T_L=T_G, p=pL1)
         rhoL1=(pL1+PinfL)/((GammaL-1.0_WP)*CvL*T_G+bL*(pL1+PinfL))
         density_ratio=rhoL1/rhoG1
         ML=1.0_WP/sqrt(GammaL*(pL1+PinfL)/(rhoL1*(1.0_WP-bL*rhoL1)))
         ! Build EOS objects
         call water%initialize(pinf=PinfL,b=bL,gamma=GammaL,cv=CvL,q=0.0_WP,qp=qpL)
         allocate(eosG(2))
         call eosG(1)%initialize(gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP)  ! inert vapor slot (mirrors air)
         call eosG(2)%initialize(gamma=GammaG,cv=CvG,q=0.0_WP,qp=0.0_WP)  ! air (carrier)
         call mixG%initialize(ns=2)
         call mixG%set_species(eosG)
         ! Relaxation model
         call relax_model%initialize(liq=water,gas=mixG,indV=1,indA=2)
         ! Viscous parameters
         call param_read('Reynolds number',Reynolds)
         call param_read('Prandtl number',Prandtl)
         call param_read('Viscosity ratio',visc_ratio)
         call param_read('Diffusivity ratio',diff_ratio)
         call param_read('Sutherland exponent',Suth_n)
         call param_read('Sutherland temperature',Suth_T)
         ! Relaxation type
         call param_read('Relaxation type',relaxation_type,default='p')
         select case (trim(relaxation_type))
         case ('p','pT','pTg')
            case_name='drop_impact_relax_'//trim(relaxation_type)
         case default
            call die('[simulation_init] Relaxation type must be p, pT, or pTg')
         end select
         ! Log
         write(message,'("[Post-shock Mach] M2=",es12.5)') M2; call log(message)
         write(message,'("[Shock Mach]      Ms=",es12.5)') Ms; call log(message)
         write(message,'("[Pre-shock]  rhoG1=",es12.5," pG1=",es12.5)') rhoG1,pG1; call log(message)
         write(message,'("[Post-shock] rhoG2=",es12.5," pG2=",es12.5)') rhoG2,pG2; call log(message)
         write(message,'("[Liquid] rhoL1=",es12.5," pL1=",es12.5," ML=",es12.5)') rhoL1,pL1,ML; call log(message)
         write(message,'("[Temp]   TL=",es12.5," TG=",es12.5)') water%get_T_from_p_rho(p=pL1,rho=rhoL1),T_G; call log(message)
         write(message,'("[Visc]   Re=",es12.5," mu*=",es12.5," Suth_n=",es12.5," Suth_T=",es12.5)') Reynolds,visc_ratio,Suth_n,Suth_T; call log(message)
         write(message,'("[Surface tension] We=",es12.5)') Weber; call log(message)
      end block init_eos_and_flow

      ! Initialize AMR grid
      create_amrgrid: block
         amr%name=trim(case_name)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         amr%xlo=  0.0_WP; call param_read('Domain length',amr%xhi)
         amr%ylo=-10.0_WP; amr%yhi=+10.0_WP
         amr%zlo=-10.0_WP; amr%zhi=+10.0_WP
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call param_read('Max level',amr%maxlvl)
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*(amr%yhi-amr%ylo)/real(amr%ny*2**amr%maxlvl,WP)
         end if
         call amr%initialize()
      end block create_amrgrid

      ! Handle restart/saves
      handle_restart: block
         call io%initialize(amr=amr,nfiles=1)
         call param_read('Restart from',restart_dir,default='')
         restarted=(len_trim(restart_dir).gt.0)
         if (restarted) restart_dir='restart/'//trim(case_name)//'_'//trim(adjustl(restart_dir))
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
            time%n=restart_step
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module, only: amrex_bc_foextrap,amrex_bc_reflect_odd
         use amrmpcomp_class,  only: BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         use messager,         only: die
         call fs%set_thermo(water,mixG)
         call fs%initialize(amr=amr,name=trim(case_name))
         fs%sigma=1.0_WP/Weber
         if (amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Wire relaxation step
         select case (trim(relaxation_type))
         case ('p')
            fs%relax=>relax_p
         case ('pT')
            fs%relax=>relax_pT
         case ('pTg')
            fs%relax=>relax_pTg
         end select
         fs%user_init=>shockdrop_init

         ! Neumann (zero-gradient) at x+
         fs%hi_bc(1)=BC_REFLECT
         fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%hi_bc(1,:)=amrex_bc_foextrap
         fs%V%hi_bc(1,:)=amrex_bc_foextrap
         fs%W%hi_bc(1,:)=amrex_bc_foextrap

         ! Wall BC at x- (90-degree contact)
         fs%lo_bc(1)=BC_REFLECT
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap
         fs%Q%lo_bc(1,5)=amrex_bc_reflect_odd
         fs%U%lo_bc(1,:)=amrex_bc_reflect_odd
         call param_read('Wall BC',wall_bc_type,default='noslip')
         select case (trim(wall_bc_type))
         case ('noslip')
            fs%Q%lo_bc(1,6:7)=amrex_bc_reflect_odd
            fs%V%lo_bc(1,:)  =amrex_bc_reflect_odd
            fs%W%lo_bc(1,:)  =amrex_bc_reflect_odd
         case ('slip')
            fs%Q%lo_bc(1,6:7)=amrex_bc_foextrap
            fs%V%lo_bc(1,:)  =amrex_bc_foextrap
            fs%W%lo_bc(1,:)  =amrex_bc_foextrap
         case default
            call die('[simulation_init] Unknown Wall BC type: must be slip or noslip')
         end select
      end block create_solver

      ! Initialize workspaces
      create_workspace: block
         use amrdata_class, only: interp_none
         call dQdt%initialize(amr,name='dQdt',ncomp=fs%nQ,ng=0,interp=interp_none); call dQdt%register()
         call Umag%initialize(amr,name='Umag',ncomp=1    ,ng=0,interp=interp_none); call Umag%register()
         call Mach%initialize(amr,name='Mach',ncomp=1    ,ng=0,interp=interp_none); call Mach%register()
      end block create_workspace

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
         call mfile%add_column(fs%TLmin,'TLmin')
         call mfile%add_column(fs%TLmax,'TLmax')
         call mfile%add_column(fs%RHOGmin,'rhoGmin')
         call mfile%add_column(fs%RHOGmax,'rhoGmax')
         call mfile%add_column(fs%PGmin,'PGmin')
         call mfile%add_column(fs%PGmax,'PGmax')
         call mfile%add_column(fs%TGmin,'TGmin')
         call mfile%add_column(fs%TGmax,'TGmax')
         call mfile%add_column(fs%dPmax,'dPmax')
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
         call consfile%add_column(time%n,'Timestep number')
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
      end block create_monitors

   end subroutine simulation_init

   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old state
         call fs%store_old()

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%build_plic(time=time%t)
         call fs%apply_relax(time=time%tmid)
         call fs%clean_Q()
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%build_plic(time=time%t)
         call fs%apply_relax(time=time%t)
         call fs%clean_Q()
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=time%dt)
         call fs%add_surface_tension(scale=time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         call fs%get_primitive(Q=fs%Q)
         ! ======================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
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
               call io%write(dirname='restart/'//trim(case_name)//'_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

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
      call time%finalize()
      call amr%finalize()
      call regrid_evt%finalize()
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      if (allocated(eosG)) deallocate(eosG)
      call viz%finalize()
      call viz_evt%finalize()
      call save_evt%finalize()
      call io%finalize()
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
   end subroutine simulation_final

end module simulation
