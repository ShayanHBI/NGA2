!> Cavitation of an energized liquid droplet (1D version)
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
   use sg_class,            only: sg
   use ig_class,            only: ig
   use igmix_class,         only: igmix
   use relax_class,         only: relax
   use relax_sg_ig_class,   only: relax_sg_ig
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
   class(ig),    allocatable, target, save  :: eosL         !< Liquid pure-substance EOS
   type(ig),     allocatable, target, save  :: eosG(:)      !< Gas-phase species EOS array
   type(igmix),  target, save               :: mixG         !< Gas mixture
   class(relax), allocatable, target, save  :: relax_model  !< Relaxation model
   character(len=str_medium), save          :: liquid_eos_type,relaxation_type
   character(len=str_medium), save          :: case_name

   !> Molar mass of the gas
   real(WP), parameter :: Mv=0.0180153_WP,Ma=0.02897_WP

   !> Case parameters
   real(WP) :: T0,p0              !< Uniform initial liquid temperature/pressure
   real(WP) :: rhoL0,eL0          !< Liquid density/energy at (p0,T0)
   real(WP) :: U_exit,U_core,U_slp
   real(WP) :: p_cav,VF_nuc       !< Cavitation onset pressure threshold and nucleation seed VF
   real(WP) :: muG,muL            !< Dynamic viscosities
   real(WP) :: PrL,PrG,ScV        !< Prandtle and Schmidt numbers

   !> Domain dimensions
   real(WP) :: Lx,Ly              !< Domain lengths [m]

   !> Tagging parameters
   real(WP) :: vorticity_tag=huge(1.0_WP)
   real(WP) :: rho_ratio_tag=huge(1.0_WP)

   !> Vapor mass fraction integral within a fixed radius of the droplet center
   real(WP), parameter :: Yv_int_radius=0.01_WP
   real(WP) :: Yv_int_r=0.0_WP

   !> Time stepping
   real(WP) :: dt_init

contains

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

   !> Levelset function for a domain that is liquid everywhere
   function levelset_liquid(xyz,t) result(G)
      real(WP), dimension(3), intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      G=1.0_WP
   end function levelset_liquid


   function get_radial_velocity(x) result(Ux)
      real(WP), intent(in) :: x
      real(WP) :: Ux
      ! Ux=U_exit*x/sqrt(x**2+Rc**2)
      ! Ux=U_exit
      Ux=sign(U_slp*abs(x)+U_core,x)
   end function get_radial_velocity

   !> Compute viscosity: constant gas and liquid, VF-weighted blend
   !> Contains commented-out Sutherland law for variable gas viscosity (dimensional form)
   subroutine get_viscosities()
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      integer :: lvl,i,j,k
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pTG,pVF,pYg,pVisc,pBeta,pCond,pDiff
      real(WP) :: mu_g,mu_l,cp_g
      real(WP), parameter :: myeps=1.0e-15_WP
      !> Sutherland's law parameters (dimensional, SI units)
      !> mu_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K
      ! real(WP), parameter :: mu_ref=1.716e-5_WP   !< Reference viscosity [Pa·s]
      ! real(WP), parameter :: T_ref=273.15_WP      !< Reference temperature [K]
      ! real(WP), parameter :: S_suth=110.4_WP      !< Sutherland constant [K]
      ! real(WP) :: T_gas                           !< Local gas temperature [K]
      ! Loop over levels
      do lvl=0,amr%clvl()
         call amr%mfiter_build(lvl,mfi)
         do while (mfi%next())
            ! Get pointers to data
            pTG=>fs%TG%mf(lvl)%dataptr(mfi)
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pYg=>fs%Yg%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pCond=>fs%cond%mf(lvl)%dataptr(mfi)
            pDiff=>fs%diff%mf(lvl)%dataptr(mfi)
            ! Get tilebox with overlap
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               ! --- Gas viscosity ---
               ! Constant gas viscosity (current default)
               mu_g=muG
               ! Sutherland's law for variable gas viscosity (dimensional):
               ! T_gas=pTG(i,j,k,1)  ! Gas temperature from EOS [K]
               ! mu_g=mu_ref*(T_gas/T_ref)**1.5_WP*(T_ref+S_suth)/(T_gas+S_suth)
               ! --- Liquid viscosity ---
               mu_l=muL
               ! Mixture viscosity (harmonic averaging)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(mu_l,myeps)+(1.0_WP-pVF(i,j,k,1))/max(mu_g,myeps))
               ! Zero bulk viscosity
               pBeta(i,j,k,1)=0.0_WP
               ! Thermal conductivity (zero when Pr=0 i.e. no heat conduction)
               cp_g=pYg(i,j,k,1)*CpV+(1.0_WP-pYg(i,j,k,1))*CpA
               if (PrL.gt.0.0_WP) then
                  pCond(i,j,k,1)=mu_l*CpL/PrL
               else
                  pCond(i,j,k,1)=0.0_WP
               end if
               if (PrG.gt.0.0_WP) then
                  pCond(i,j,k,2)=mu_g*cp_g/PrG
               else
                  pCond(i,j,k,2)=0.0_WP
               end if
               ! Vapor mass diffusivity (zero when Sc=0 i.e. no mass diffusion)
               if (ScV.gt.0.0_WP) then
                  pDiff(i,j,k,1)=mu_g/ScV
               else
                  pDiff(i,j,k,1)=0.0_WP
               end if
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback – set Q and VF/barycenters for a uniform liquid domain
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
      real(WP) :: x_cc,Ux,rho_mix
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
            eL_local  =eosL%get_e_from_p_T(p=p0,T=T0)
            rhoL_local=eosL%get_rho_from_p_T(p=p0,T=T0)
            ! Get gas density and internal energy from p and T
            eG_local  =mixG%get_e_from_p_T(p=p0,T=T0,y=y)
            rhoG_local=mixG%get_rho_from_p_T(p=p0,T=T0,y=y)
            ! Set conserved variables
            pQ(i,j,k,1)=(       myVF)*rhoL_local
            pQ(i,j,k,2)=(1.0_WP-myVF)*rhoG_local
            pQ(i,j,k,3)=pQ(i,j,k,1)*eL_local
            pQ(i,j,k,4)=pQ(i,j,k,2)*eG_local
            x_cc=solver%amr%xlo+(real(i,WP)+0.5_WP)*dx
            Ux=get_radial_velocity(x_cc)
            rho_mix=pQ(i,j,k,1)+pQ(i,j,k,2)
            pQ(i,j,k,5)=rho_mix*Ux
            pQ(i,j,k,6)=0.0_WP
            pQ(i,j,k,7)=0.0_WP
            pQ(i,j,k,8)=pQ(i,j,k,2)*Yv0
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine cavitation_init

   !> Tagger based on vorticity and density ratio (from amrcomp_drop)
   subroutine my_tagger(solver,lvl,time,tags_ptr)
      use iso_c_binding,    only: c_ptr,c_char
      use amrex_amr_module, only: amrex_mfiter,amrex_box,amrex_tagboxarray
      use amrgrid_class,    only: SETtag
      use amrmpcomp_class,  only: VFlo
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
            ! Always keep any liquid-containing cell at the finest level, so that
            ! apply_relax (which only processes lvl=maxlvl) sees every droplet cell
            ! every step -- bulk interior cells in tension must never be stranded
            ! on a coarser level where cavitation nucleation can never trigger
            if (pVF(i,j,k,1).gt.VFlo) tagarr(i,j,k,1)=SETtag
         end do; end do; end do
      end do
      call solver%amr%mfiter_destroy(mfi)
   end subroutine my_tagger

   !> User BC callback – outward expansion velocity Dirichlet condition on every boundary
   subroutine radial_dirichlet_bc(solver,lvl,time,face,bx,comp,p)
      use amrex_amr_module, only: amrex_box
      class(amrmpcomp), intent(inout) :: solver
      integer, intent(in) :: lvl
      real(WP), intent(in) :: time
      integer, intent(in) :: face
      type(amrex_box), intent(in) :: bx
      character(len=1), intent(in) :: comp
      real(WP), dimension(:,:,:,:), contiguous, pointer :: p
      real(WP) :: dx,x,Ux
      integer :: i,j,k
      dx=amr%dx(lvl)
      do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
         ! Position: staggered in the component's own direction, cell-centered otherwise
         if (comp.eq.'U') then
            x=amr%xlo+real(i,WP)*dx
         else
            x=amr%xlo+(real(i,WP)+0.5_WP)*dx
         end if
         Ux=get_radial_velocity(x)
         select case (comp)
         case ('U')
            p(i,j,k,1)=Ux
         case ('V')
            p(i,j,k,1)=0.0_WP
         case ('W')
            p(i,j,k,1)=0.0_WP
         case ('Q')
            p(i,j,k,1)=rhoL0; p(i,j,k,2)=0.0_WP
            p(i,j,k,3)=rhoL0*eL0; p(i,j,k,4)=0.0_WP
            p(i,j,k,5)=rhoL0*Ux; p(i,j,k,6)=0.0_WP; p(i,j,k,7)=0.0_WP
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
         ! Prandtl and Schmidt numbers
         call param_read('Liquid Prandtl number',PrL)
         call param_read('Gas Prandtl number',PrG)
         call param_read('Vapor Schmidt number',ScV)
         ! Uniform liquid state and boundary velocity
         call param_read('Liquid temperature',T0)
         call param_read('Liquid pressure',p0)
         call param_read('Exit velocity',U_exit)
         call param_read('Core velocity',U_core)
         ! Delayed cavitation onset: nucleate at p_cav (deeper than p_sat) with a VF_nuc-sized seed
         call param_read('Cavitation pressure threshold',p_cav,default=huge(1.0_WP))
         call param_read('Cavitation seed VF',VF_nuc,default=1.0e-7_WP)
         ! Domain dimensions
         call param_read('Lx',Lx)
         call param_read('Ly',Ly)
         U_slp=(U_exit-U_core)/Lx
         ! Select liquid EOS type and allocate eosL and relax_model
         call param_read('Liquid EOS type',liquid_eos_type)
         call param_read('Relaxation type',relaxation_type)
         select case(relaxation_type)
         case('p','pT','pTg')
            case_name='cavitation1D_'//trim(liquid_eos_type)//'_relax_'//trim(relaxation_type)
         case default
            call die('Relaxation type has to be either p, pT, or pTg')
         end select
         select case (trim(liquid_eos_type))
         case ('SG')
            allocate(sg :: eosL)
            allocate(relax_sg_ig :: relax_model)
         case ('NASG')
            allocate(nasg :: eosL)
            allocate(relax_nasg_ig :: relax_model)
         case default
            call die('[simulation] Unknown Liquid EOS type: '//trim(liquid_eos_type))
         end select
         ! Initialize EOS objects
         select type (eosL)
         type is (sg)
            call eosL%initialize(pinf=PinfL,gamma=GammaL,cv=CvL,q=qL,qp=qpL)
         type is (nasg)
            call eosL%initialize(pinf=PinfL,b=bL,gamma=GammaL,cv=CvL,q=qL,qp=qpL)
         end select
         ! Liquid density and energy at the uniform initial state
         rhoL0=eosL%get_rho_from_p_T(p=p0,T=T0)
         eL0  =eosL%get_e_from_p_T  (p=p0,T=T0)
         allocate(eosG(2))
         call eosG(1)%initialize(gamma=GammaV,cv=CvV,q=qV,qp=qpV)   ! species 1 = vapor (transported)
         call eosG(2)%initialize(gamma=GammaA,cv=CvA,q=qA,qp=qpA)   ! species 2 = air (carrier)
         call mixG%initialize(ns=2)
         call mixG%set_species(eosG)
         select type (relax_model)
         type is (relax_sg_ig)
            select type (eosL)
            type is (sg)
               call relax_model%initialize(liq=eosL,gas=mixG,indV=1,indA=2,p_cav=p_cav,VF_nuc=VF_nuc)
            end select
         type is (relax_nasg_ig)
            select type (eosL)
            type is (nasg)
               call relax_model%initialize(liq=eosL,gas=mixG,indV=1,indA=2,p_cav=p_cav,VF_nuc=VF_nuc)
            end select
         end select
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
         ! Enable quasi-1D: collapse y and z to a single periodic cell matching the finest x spacing
         if (amr%ny.eq.1) then
            amr%ylo=-0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP)
            amr%yhi=+0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP)
            amr%yper=.true.
         end if
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP)
            amr%zhi=+0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP)
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
         if (restarted) restart_dir='restart/'//trim(case_name)//'_'//trim(adjustl(restart_dir))
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
         use amrex_amr_module, only: amrex_bc_ext_dir, amrex_bc_reflect_even, amrex_bc_reflect_odd
         use amrmpcomp_class,  only: BC_REFLECT
         use amrdata_class,    only: interp_face_lin
         ! Use piecewise-linear face interpolation — FaceDivFree requires ratio==2 in all dirs
         ! but this case is quasi-2D with ref_ratio_z=1
         fs%interp_vel=interp_face_lin
         ! Assign EOS objects and create flow solver
         call fs%set_thermo(eosL,mixG)
         call fs%initialize(amr=amr,name=trim(case_name))
         ! Provide relaxation step
         select case(relaxation_type)
         case('p')
            fs%relax=>relax_p
         case('pT')
            fs%relax=>relax_pT
         case('pTg')
            fs%relax=>relax_pTg
         end select
         ! Set initial conditions via cavitation callback
         fs%user_init=>cavitation_init
         ! Set BCs: radial outward velocity Dirichlet on every non-periodic boundary
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
            ! Q(8) = (1-VF)*rhoG*Yv is restored from the checkpoint; just recompute primitives
            call fs%Q%fill(time=time%t)
            call fs%get_primitive(fs%Q)
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
         call Umag%get_magnitude(fs%UVW,fs%UVW,fs%UVW,compX=1,compY=2,compZ=3)
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
         ! Compute Yv integral within fixed radius
         ! call get_Yv_int_r()
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
         call mfile%add_column(fs%Yvmin,'Yvmin')
         call mfile%add_column(fs%Yvmax,'Yvmax')
         call mfile%add_column(fs%Qint(8),'Vapor mass')
         ! call mfile%add_column(Yv_int_r,'Yv_int_r01')
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
            time%dt=dt_init
         else
            call time%adjust_dt()
         end if
         call time%increment()

         ! Remember old state
         call fs%Qold%copy(src=fs%Q)
         call fs%store_old()

         ! ===== RK2 Stage 1: advective flux at t =====
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%t)
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=0.5_WP*time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%build_plic(time%t)
         call fs%apply_relax(time=time%t+0.5_WP*time%dt)
         call fs%clean_Q()
         call fs%get_primitive(fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity()
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t+0.5_WP*time%dt)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t+0.5_WP*time%dt)
         call fs%get_primitive(fs%Q)

         ! ===== RK2 Stage 2: advective flux at midpoint =====
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t+0.5_WP*time%dt)
         call fs%Q%copy(src=fs%Qold); call fs%Q%saxpy(a=time%dt,src=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%build_plic(time%t)
         call fs%apply_relax(time=time%t)
         call fs%clean_Q()
         call fs%get_primitive(fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity()
         call fs%add_phasic_pressure(scale=time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         call fs%get_primitive(fs%Q)

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
         call Umag%get_magnitude(fs%UVW,fs%UVW,fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)
         ! if (fs%amr%amRoot) print*, 'time: ', time%t, ' dt: ', time%dt,' occurs: ',viz_evt%occurs()
         ! Visualization output
         if (viz_evt%occurs()) then
            call viz%write(time%t)
         end if
         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/'//trim(case_name)//'_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
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
         call io%write(dirname='restart/'//trim(case_name)//'_'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
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
