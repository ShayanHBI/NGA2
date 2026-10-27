!> AMR compressible water cavitation tube (Pelanti 2022, Sec. 8.2.1), starting from pure liquid
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
   ! use relax_igmix_nasg_finite_class, only: relax_igmix_nasg_finite
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
   integer  :: restart_step

   !> Simulation monitoring
   type(monitor) :: mfile,consfile,cflfile,gridfile,tfile,rescfile
   !> Relaxation-model census (relax_model%acc reduced across ranks for the rescue monitor)
   real(WP) :: diss_n=0.0_WP,diss_m=0.0_WP
   real(WP) :: quad_n=0.0_WP,swap_n=0.0_WP,flr_n=0.0_WP,flr_e=0.0_WP,stuck_n=0.0_WP
   real(WP) :: sliverL_n=0.0_WP,sliverLf_n=0.0_WP,chemskip_n=0.0_WP

   !> Materials
   type(nasg),  target :: water
   type(igmix), target :: gas

   !> Relaxation model
   ! type(relax_igmix_nasg_finite), target :: relax_model
   type(relax_igmix_nasg), target :: relax_model
   character(len=str_medium) :: relaxation_type,case_name
   integer :: nsteps_diag                            !< diag run: bounded number of steps (-1 = unlimited)

   !> Flow parameters
   real(WP) :: p0,T0                  !< Initial liquid pressure and temperature
   real(WP) :: rhoL0,eL0              !< Initial liquid density and specific internal energy
   real(WP) :: u0,x0                  !< Expansion velocity magnitude and location of the velocity jump

   !> Transport properties (constant) and subgrid models
   real(WP) :: muL,muG                !< Dynamic viscosities
   real(WP) :: kL,kG                  !< Thermal conductivities
   real(WP) :: Cartif                 !< Artificial bulk viscosity coefficient (0: off)
   real(WP) :: Cvisc                  !< Artificial shear-to-bulk viscosity ratio
   logical  :: use_vreman             !< Vreman SGS model

   !> Tagging parameter
   real(WP) :: Rho_tag=huge(1.0_WP)
   real(WP) :: P_tag=huge(1.0_WP)
   real(WP) :: Ducros_tag=huge(1.0_WP)

   !> debug: line probe along x (j=k=0)
   integer  :: probe_ilo,probe_ihi,probe_stride ! debug
   real(WP) :: probe_t0 ! debug

contains


   !> debug: one CSV row per cell i in [probe_ilo,probe_ihi] (j=k=0) at a pipeline stage, every probe_stride steps from probe_t0
   subroutine probe_line(label)
      use amrex_amr_module, only: amrex_mfiter,amrex_box
      implicit none
      character(len=*), intent(in) :: label
      integer :: lvl,i
      type(amrex_mfiter) :: mfi
      type(amrex_box) :: bx
      real(WP), dimension(:,:,:,:), contiguous, pointer :: pVF,pQ,pU,pSub,pPLs,pPGs
      real(WP) :: VFc,RHOL,RHOG,PL,PG,TL,TG,Yv
      integer, save :: iunit
      logical, save :: opened=.false.
      character(len=32) :: fname
      if (probe_ihi.lt.probe_ilo.or.time%t.lt.probe_t0) return
      if (mod(time%n,probe_stride).ne.0) return
      if (.not.opened) then
         opened=.true.
         write(fname,'("probe_r",i0,".csv")') amr%rank
         open(newunit=iunit,file=trim(fname),status='replace',action='write')
         write(iunit,'(a)') 'label,n,t,i,VF,Q1,Q2,Q3,Q4,Q5,Q8,Uface_lo,Uface_hi,Ucell,RHOL,RHOG,PL,PG,TL,TG,subVF_lo,subVF_hi,PLs,PGs'
      end if
      lvl=amr%maxlvl
      call amr%mfiter_build(lvl,mfi)
      do while (mfi%next())
         pVF=>fs%VF%mf(lvl)%dataptr(mfi); pQ=>fs%Q%mf(lvl)%dataptr(mfi); pU=>fs%U%mf(lvl)%dataptr(mfi)
         pSub=>fs%subVF%dataptr(mfi); pPLs=>fs%PL%mf(lvl)%dataptr(mfi); pPGs=>fs%PG%mf(lvl)%dataptr(mfi)
         bx=mfi%tilebox()
         if (bx%lo(2).gt.0.or.bx%hi(2).lt.0.or.bx%lo(3).gt.0.or.bx%hi(3).lt.0) cycle
         do i=max(probe_ilo,bx%lo(1)),min(probe_ihi,bx%hi(1))
            VFc=pVF(i,0,0,1)
            RHOL=0.0_WP; RHOG=0.0_WP; PL=0.0_WP; PG=0.0_WP; TL=0.0_WP; TG=0.0_WP
            if (VFc.gt.0.0_WP.and.pQ(i,0,0,1).gt.0.0_WP) then
               RHOL=pQ(i,0,0,1)/VFc
               PL=fs%liq%get_p_from_rho_e(rho=RHOL,e=pQ(i,0,0,3)/pQ(i,0,0,1),y=[1.0_WP])
               TL=fs%liq%get_T_from_rho_e(rho=RHOL,e=pQ(i,0,0,3)/pQ(i,0,0,1),y=[1.0_WP])
            end if
            if (VFc.lt.1.0_WP.and.pQ(i,0,0,2).gt.0.0_WP) then
               RHOG=pQ(i,0,0,2)/(1.0_WP-VFc)
               Yv=max(0.0_WP,min(pQ(i,0,0,8)/pQ(i,0,0,2),1.0_WP))
               PG=fs%gas%get_p_from_rho_e(rho=RHOG,e=pQ(i,0,0,4)/pQ(i,0,0,2),y=[Yv,1.0_WP-Yv])
               TG=fs%gas%get_T_from_rho_e(rho=RHOG,e=pQ(i,0,0,4)/pQ(i,0,0,2),y=[Yv,1.0_WP-Yv])
            end if
            write(iunit,'(a,",",i0,",",es16.8,",",i0,20(",",es16.8))') trim(label),time%n,time%t,i,VFc,pQ(i,0,0,1),pQ(i,0,0,2), &
            &  pQ(i,0,0,3),pQ(i,0,0,4),pQ(i,0,0,5),pQ(i,0,0,8),pU(i,0,0,1),pU(i+1,0,0,1),pQ(i,0,0,5)/max(pQ(i,0,0,1)+pQ(i,0,0,2),tiny(1.0_WP)), &
            &  RHOL,RHOG,PL,PG,TL,TG,pSub(i,0,0,1),pSub(i,0,0,2),pPLs(i,0,0,1),pPGs(i,0,0,1)
         end do
      end do
      call amr%mfiter_destroy(mfi)
      flush(iunit)
   end subroutine probe_line

   !> Compute viscosities and heat diffusivities: constant phasic values, harmonic VF blend for viscosity
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
            pVF=>fs%VF%mf(lvl)%dataptr(mfi)
            pVisc=>fs%visc%mf(lvl)%dataptr(mfi)
            pBeta=>fs%beta%mf(lvl)%dataptr(mfi)
            pDiffL=>fs%diffL%mf(lvl)%dataptr(mfi)
            pDiffG=>fs%diffG%mf(lvl)%dataptr(mfi)
            bx=mfi%growntilebox(fs%nover)
            do k=bx%lo(3),bx%hi(3); do j=bx%lo(2),bx%hi(2); do i=bx%lo(1),bx%hi(1)
               pVisc(i,j,k,1)=1.0_WP/(pVF(i,j,k,1)/max(muL,myeps)+(1.0_WP-pVF(i,j,k,1))/max(muG,myeps))
               pBeta(i,j,k,1)=0.0_WP
               pDiffL(i,j,k,1)=kL
               pDiffG(i,j,k,1)=kG
            end do; end do; end do
         end do
         call amr%mfiter_destroy(mfi)
      end do
   end subroutine get_viscosities

   !> User init callback: pure liquid at (p0,T0), velocity -u0 left of x0 and +u0 right of it; the cell holding x0 gets the volume-weighted momentum
   subroutine cavtube_init(solver,lvl,time,ba,dm)
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
      real(WP) :: dx,dy,dz,fL
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
            pVF(i,j,k,1)=1.0_WP
            if (lvl.eq.solver%amr%maxlvl) then
               pCL(i,j,k,:)=[solver%amr%xlo+(real(i,WP)+0.5_WP)*dx,solver%amr%ylo+(real(j,WP)+0.5_WP)*dy,solver%amr%zlo+(real(k,WP)+0.5_WP)*dz]
               pCG(i,j,k,:)=pCL(i,j,k,:)
            end if
            ! Fraction of the cell left of x0
            fL=min(max((x0-solver%amr%xlo-real(i,WP)*dx)/dx,0.0_WP),1.0_WP)
            pQ(i,j,k,:)=0.0_WP
            pQ(i,j,k,1)=rhoL0
            pQ(i,j,k,3)=rhoL0*eL0
            pQ(i,j,k,5)=rhoL0*u0*(1.0_WP-2.0_WP*fL)
         end do; end do; end do
      end do
      call amrex_mfiter_destroy(mfi)
   end subroutine cavtube_init

   !> Tagger based on density/pressure errors and Ducros sensor
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
      real(WP) :: dxi,dyi,dzi,delta
      real(WP) :: rho_cc,rho_xp,rho_xm,rho_yp,rho_ym,rho_zp,rho_zm
      real(WP) :: divu,vortx,vorty,vortz,vort,Ducros,Deps
      integer :: i,j,k
      real(WP), parameter :: Reps=1.0e-2_WP
      real(WP), parameter :: Peps=1.0e-2_WP
      real(WP), parameter :: Cduc=0.05_WP
      dxi=1.0_WP/solver%amr%dx(lvl); dyi=1.0_WP/solver%amr%dy(lvl); dzi=1.0_WP/solver%amr%dz(lvl)
      delta=solver%amr%min_meshsize(lvl)
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

   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none

      ! Determine case name from relaxation type
      call param_read('Relaxation type',relaxation_type)
      case_name='cavtube_relax_'//trim(relaxation_type)

      ! Initialize AMR grid: tube along x; a single-cell y or z direction is one finest cell thick
      create_amrgrid: block
         real(WP) :: Lx
         amr%name=trim(case_name)
         call param_read('Base nx',amr%nx)
         call param_read('Base ny',amr%ny)
         call param_read('Base nz',amr%nz)
         call param_read('Max level',amr%maxlvl)
         call param_read('Domain length',Lx)
         amr%xlo=0.0_WP; amr%xhi=Lx
         ! Square base cells in y and z
         amr%ylo=-0.5_WP*Lx*real(amr%ny,WP)/real(amr%nx,WP); amr%yhi=-amr%ylo
         amr%zlo=-0.5_WP*Lx*real(amr%nz,WP)/real(amr%nx,WP); amr%zhi=-amr%zlo
         if (amr%ny.eq.1) then
            amr%ylo=-0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP); amr%yhi=-amr%ylo
         end if
         if (amr%nz.eq.1) then
            amr%zlo=-0.5_WP*Lx/real(amr%nx*2**amr%maxlvl,WP); amr%zhi=-amr%zlo
         end if
         amr%xper=.false.; amr%yper=.true.; amr%zper=.true.
         call amr%initialize()
      end block create_amrgrid

      ! Read EoS and flow parameters (SI units)
      init_eos_and_flow: block
         use messager, only: log
         use string,   only: str_long
         character(len=str_long) :: message
         real(WP) :: GammaL,PinfL,bL,CvL,qL,qpL
         real(WP) :: GammaV,CvV,qV,qpV
         real(WP) :: GammaA,CvA
         ! Liquid water (NASG)
         call param_read('Liquid gamma',GammaL)
         call param_read('Liquid pinf',PinfL)
         call param_read('Liquid covolume',bL)
         call param_read('Liquid cv',CvL)
         call param_read('Liquid q',qL)
         call param_read('Liquid qp',qpL)
         ! Vapor (ideal gas, transported species 1)
         call param_read('Vapor gamma',GammaV)
         call param_read('Vapor cv',CvV)
         call param_read('Vapor q',qV)
         call param_read('Vapor qp',qpV)
         ! Air (ideal gas, carrier species 2)
         call param_read('Air gamma',GammaA)
         call param_read('Air cv',CvA)
         call gas%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,0.0_WP],qp=[qpV,0.0_WP],species_names=['vapor','air  '],name='gas')
         call water%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=qL,qp=qpL,name='water')
         ! Initial liquid state and expansion
         call param_read('Initial pressure',p0)
         call param_read('Initial temperature',T0)
         call param_read('Expansion velocity',u0)
         call param_read('Expansion location',x0)
         rhoL0=water%get_rho_from_p_T(p=p0,T=T0,y=[1.0_WP])
         eL0=water%get_e_from_p_rho(p=p0,rho=rhoL0,y=[1.0_WP])
         ! Transport properties and subgrid models
         call param_read('Liquid viscosity',muL,default=0.0_WP)
         call param_read('Gas viscosity',muG,default=0.0_WP)
         call param_read('Liquid thermal conductivity',kL,default=0.0_WP)
         call param_read('Gas thermal conductivity',kG,default=0.0_WP)
         call param_read('Artificial bulk viscosity coefficient',Cartif,default=0.0_WP)
         call param_read('Artificial viscosity coefficient',Cvisc,default=0.0_WP)
         call param_read('Use Vreman',use_vreman,default=.false.)
         ! Log
         write(message,'("[Liquid] p0=",es12.5," T0=",es12.5," rhoL0=",es12.5," cL0=",es12.5)') p0,T0,rhoL0,water%get_c_from_p_rho(p=p0,rho=rhoL0,y=[1.0_WP]); call log(message)
         write(message,'("[Expansion] u0=",es12.5," x0=",es12.5)') u0,x0; call log(message)
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
         time%dt=time%dtmax
         call param_read('Max CFL',time%cflmax)
         if (restarted) then
            call io%get_scalar('dt',time%dt)
            time%t=restart_time
            time%n=restart_step
            call param_read('Max steps',nsteps_diag,default=-1)
            if (nsteps_diag.gt.0) time%nmax=restart_step+nsteps_diag   ! diag run: bounded step count
         end if
      end block initialize_timetracker

      ! Initialize compressible multiphase solver
      create_solver: block
         use amrex_amr_module,     only: amrex_bc_foextrap
         use amrmpcomp_class,      only: BC_LIQ
         use amrdata_class,        only: interp_face_lin
         use messager,             only: log,die
         use string,               only: str_long
         use relax_igmix_sg_class, only: Prelax,PTrelax,PTgrelax,PThybrid,dbg_ilo,dbg_ihi,dbg_jlo,dbg_jhi,dbg_klo,dbg_khi
         character(len=str_long) :: message
         real(WP) :: sigma,psat0,Tsat0
         logical :: conv
         integer :: it
         ! Assign materials and create flow solver
         fs%liq=>water; fs%gas=>gas; call fs%initialize(amr=amr,name=trim(case_name))
         ! Surface tension coefficient
         call param_read('Surface tension',sigma,default=0.0_WP)
         fs%sigma=sigma
         ! Use face-linear interp if a direction is a single cell (divfree requires ratio=2 in all dirs)
         if (amr%ny.eq.1.or.amr%nz.eq.1) fs%interp_vel=interp_face_lin
         ! Initialize relaxation model: species 1=vapor (transported), 2=air (carrier)
         call relax_model%initialize(liq=water,gas=gas,indV=1,indA=2)
         select case (trim(relaxation_type))
         case ('p');        relax_model%model=Prelax
         case ('pT');       relax_model%model=PTrelax
         case ('pTg');      relax_model%model=PTgrelax
         case ('pThybrid'); relax_model%model=PThybrid
         case default; call die('[simulation_init] Relaxation type must be p, pT, pThybrid, or pTg')
         end select
         ! Finite relaxation rates (Pelanti 2022: lambda=mu, theta=vartheta, nu); omitted = instantaneous
         ! call param_read('Mechanical relaxation rate',relax_model%lambda,default=huge(1.0_WP))
         ! call param_read('Thermal relaxation rate',relax_model%theta,default=huge(1.0_WP))
         ! call param_read('Chemical relaxation rate',relax_model%nu,default=huge(1.0_WP))
         ! write(message,'("[Relaxation rates] lambda=",es12.5," theta=",es12.5," nu=",es12.5)') relax_model%lambda,relax_model%theta,relax_model%nu
         ! call log(message)
         ! Nucleation thresholds in SI: liquid pressure [Pa] below which cavitation seeds, condensation tolerance [K]
         if (trim(relaxation_type).eq.'pTg') then
            call param_read('Cavitation pressure threshold',relax_model%p_cav)
            call param_read('Condensation temperature tolerance',relax_model%Tctol)
         end if
         ! Saturation reference of the EOS pair: psat at T0 and Tsat at p0 (pure vapor)
         call relax_model%get_psat(T0,psat0,conv)
         call relax_model%get_Tsat(p0,p0,T0,Tsat0,conv,it)
         write(message,'("[Saturation] psat(T0)=",es12.5," Tsat(p0)=",es12.5," p_cav=",es12.5)') psat0,Tsat0,relax_model%p_cav
         call log(message)
         relax_model%RHOGmin=0.0_WP
         ! call param_read('Relax skip gas density',relax_model%rhoG_skip,default=0.0_WP)
         ! call param_read('Relax sliver gas fraction',relax_model%sliver_aG,default=0.0_WP)
         ! call param_read('Relax liquid thermal tolerance',relax_model%dTtol_L,default=0.0_WP)
         ! call param_read('Relax vapor fraction tolerance',relax_model%dYtol,default=0.0_WP)
         ! Pass A0: a sliver whose gas, and a neighbour whose gas, sit below this fraction of the sliver's PL are mixed isentropically (0: off)
         call param_read('Relax vacuum pressure ratio',fs%vacmix_pratio,default=0.0_WP)
         ! fs%vacmix_aG=relax_model%sliver_aG
         relax_model%vol=amr%cell_vol(amr%maxlvl)
         fs%merge_sick=100.0_WP
         ! Floors in SI
         fs%Pmin_liq=-0.98_WP*water%pinf
         fs%Tmin_liq=250.0_WP
         fs%Pmin_gas=1.0e-4_WP
         fs%Tmin_gas=250.0_WP
         relax_model%Pmin_liq=fs%Pmin_liq; relax_model%Tmin_liq=fs%Tmin_liq
         relax_model%Pmin_gas=fs%Pmin_gas; relax_model%Tmin_gas=fs%Tmin_gas
         ! Set relaxation model for the flow solver
         fs%relax=>relax_model
         fs%cluster_on=.false.
         call param_read('Mechanical clustering',fs%cluster_mech_on,default=.false.) ! Pass A1
         call param_read('Thermal clustering',fs%cluster_therm_on,default=.false.)   ! Pass A2
         call param_read('Chemical clustering',fs%cluster_chem_on,default=.false.)  ! Pass A3
         call param_read('Dissolve stranded gas',fs%dissolve_on,default=.false.)
         fs%relax_diag_on=.false.
         call param_read('Diag stride',fs%relax_diag_stride,default=20)
         ! debug trace box: inclusive index range the relax/cluster prints trace (empty by default)
         call param_read('Debug ilo',dbg_ilo,default=100000000)
         call param_read('Debug ihi',dbg_ihi,default=-100000000)
         call param_read('Debug jlo',dbg_jlo,default=100000000)
         call param_read('Debug jhi',dbg_jhi,default=-100000000)
         call param_read('Debug klo',dbg_klo,default=100000000)
         call param_read('Debug khi',dbg_khi,default=-100000000)
         ! debug: line probe (off unless Probe ihi>=Probe ilo)
         call param_read('Probe ilo',probe_ilo,default=1) ! debug
         call param_read('Probe ihi',probe_ihi,default=0) ! debug
         call param_read('Probe start time',probe_t0,default=0.0_WP) ! debug
         call param_read('Probe stride',probe_stride,default=1) ! debug
         ! Set initial conditions
         fs%user_init=>cavtube_init
         ! Transmissive (foextrap) liquid ends in x
         fs%lo_bc(1)=BC_LIQ; fs%hi_bc(1)=BC_LIQ
         fs%Q%lo_bc(1,:)=amrex_bc_foextrap; fs%Q%hi_bc(1,:)=amrex_bc_foextrap
         fs%U%lo_bc(1,:)=amrex_bc_foextrap; fs%U%hi_bc(1,:)=amrex_bc_foextrap
         fs%V%lo_bc(1,:)=amrex_bc_foextrap; fs%V%hi_bc(1,:)=amrex_bc_foextrap
         fs%W%lo_bc(1,:)=amrex_bc_foextrap; fs%W%hi_bc(1,:)=amrex_bc_foextrap
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
         call param_read('Tag density error' ,Rho_tag,default=huge(1.0_WP))
         call param_read('Tag pressure error',P_tag,default=huge(1.0_WP))
         call param_read('Tag Ducros value',Ducros_tag,default=huge(1.0_WP))
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
         if (Cartif.gt.0.0_WP) call fs%add_viscartif(dt=time%dt,Cartif=Cartif,Cvisc=Cvisc)
         if (use_vreman) call fs%add_vreman(dt=time%dt)
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
         ! Create simulation monitor
         mfile=monitor(amRoot=amr%amRoot,name='simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
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
         ! Create CFL monitor
         cflfile=monitor(amRoot=amr%amRoot,name='cfl')
         call cflfile%add_column(time%n,'Timestep')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(time%dt,'dt')
         call cflfile%add_column(fs%CFLc_x,'CFLc_x')
         call cflfile%add_column(fs%CFLa_x,'CFLa_x')
         call cflfile%add_column(fs%CFLv_x,'CFLv_x')
         call cflfile%add_column(fs%CFLst ,'CFLst' )
         call cflfile%write()
         ! Create conservation monitor
         consfile=monitor(amRoot=amr%amRoot,name='conservation')
         call consfile%add_column(time%n,'Timestep number')
         call consfile%add_column(time%t,'Time')
         call consfile%add_column(fs%VFint,'VFint')
         call consfile%add_column(fs%Qint(1),'Liquid Mass')
         call consfile%add_column(fs%Qint(2),'Gas Mass')
         call consfile%add_column(fs%Qint(3),'Liquid IntEnergy')
         call consfile%add_column(fs%Qint(4),'Gas IntEnergy')
         call consfile%add_column(fs%Qint(5),'U Momentum')
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
         call gridfile%write()
         ! Create timing monitor
         tfile=monitor(amRoot=amr%amRoot,name='timing')
         call tfile%add_column(time%n,'Timestep')
         call tfile%add_column(time%t,'Time')
         call tfile%add_column(fs%wtmax_dQdt,'dQdt_max')
         call tfile%add_column(fs%wtmax_plic,'plic_max')
         call tfile%add_column(fs%wtmax_relax,'relax_max')
         call tfile%add_column(fs%wtmax_visc,'visc_max')
         call tfile%add_column(fs%nmixed_max,'mixed_max')
         call tfile%write()
         ! Create rescue-census monitor (cumulative counters/amounts per mechanism)
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
         call rescfile%add_column(fs%vacpair_n,'A0 pairs')
         call rescfile%add_column(sliverL_n,'L sliver n')
         call rescfile%add_column(sliverLf_n,'L sliver fail')
         call rescfile%add_column(chemskip_n,'Chem skip n')
         call rescfile%write()
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
         fs%relax_diag_step=time%n

         ! Remember old state
         call fs%store_old()
         call probe_line('A0_start') ! debug

         ! ======================= RK2 Stage 1: Q*=Q[n]+dt/2*dQdt(t,Q[n]) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=0.5_WP*time%dt,time=time%tmid)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=0.5_WP*time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call probe_line('S1a_transport') ! debug
         call fs%build_plic(time=time%t)
         call probe_line('S1b_plic') ! debug
         ! call fs%apply_relax(dt=0.5_WP*time%dt,time=time%tmid)
         call fs%apply_relax_splitremix(dt=0.5_WP*time%dt,time=time%tmid)
         call probe_line('S1c_relax') ! debug
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=0.5_WP*time%dt)
         call fs%add_surface_tension(scale=0.5_WP*time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%tmid)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%tmid)
         call fs%get_primitive(Q=fs%Q)
         call probe_line('S1d_end') ! debug
         ! ======================= RK2 Stage 2: Q[n+1]=Q[n]+dt*dQdt(t,Q*) =======================
         call fs%get_dQdt(dQdt=dQdt,dt=time%dt,time=time%t)
         call fs%Q%lincomb(a=1.0_WP,src1=fs%Qold,b=time%dt,src2=dQdt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call probe_line('S2a_transport') ! debug
         call fs%build_plic(time=time%t)
         call probe_line('S2b_plic') ! debug
         ! call fs%apply_relax(dt=time%dt,time=time%t)
         call fs%apply_relax_splitremix(dt=time%dt,time=time%t)
         call probe_line('S2c_relax') ! debug
         call fs%get_primitive(Q=fs%Q)
         call fs%build_subVF()
         call fs%get_face_velocity(); call fs%average_down_velocity()
         call fs%add_phasic_pressure(scale=time%dt)
         call fs%add_surface_tension(scale=time%dt)
         call fs%Q%average_down(); call fs%Q%fill(time=time%t)
         call fs%average_down_velocity(); call fs%fill_velocity(time=time%t)
         call fs%get_primitive(Q=fs%Q)
         call probe_line('S2d_end') ! debug
         ! ======================================================================================

         ! Regrid if event triggers
         if (regrid_evt%occurs()) then
            call amr%regrid(baselvl=0,time=time%t)
            call gridfile%write()
         end if

         ! Compute viscosities and subgrid models
         call get_viscosities()
         if (Cartif.gt.0.0_WP) call fs%add_viscartif(dt=time%dt,Cartif=Cartif,Cvisc=Cvisc)
         if (use_vreman) call fs%add_vreman(dt=time%dt)

         ! Compute Umag and Mach number
         call Umag%get_magnitude(srcX=fs%UVW,srcY=fs%UVW,srcZ=fs%UVW,compX=1,compY=2,compZ=3)
         call Mach%copy(src=Umag); call Mach%divide(src=fs%C)

         ! Visualization output
         if (viz_evt%occurs()) call viz%write(time=time%t)

         ! Checkpoint save
         if (save_evt%occurs()) then
            save_checkpoint: block
               use string, only: rtoa
               call io%write(dirname='restart/'//trim(case_name)//'/'//trim(adjustl(rtoa(time%t))),time=time%t,step=time%n)
            end block save_checkpoint
         end if

         ! Perform and output monitoring
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
            ! sliverL_n=relax_model%nsliver_L
            ! call MPI_ALLREDUCE(MPI_IN_PLACE,sliverL_n,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
            ! sliverLf_n=relax_model%nsliver_Lf
            ! call MPI_ALLREDUCE(MPI_IN_PLACE,sliverLf_n,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
            ! chemskip_n=relax_model%nchem_skip
            ! call MPI_ALLREDUCE(MPI_IN_PLACE,chemskip_n,1,MPI_REAL_WP,MPI_SUM,amr%comm,ierr)
         end block relax_census
         call mfile%write()
         call consfile%write()
         call cflfile%write()
         call tfile%write()
         call rescfile%write()

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
      call time%finalize()
      call amr%finalize()
      call regrid_evt%finalize()
      call fs%finalize()
      call dQdt%finalize()
      call Umag%finalize()
      call Mach%finalize()
      call water%finalize()
      call gas%finalize()
      call viz%finalize()
      call viz_evt%finalize()
      call save_evt%finalize()
      call io%finalize()
      call mfile%finalize()
      call cflfile%finalize()
      call consfile%finalize()
      call gridfile%finalize()
      call tfile%finalize()
      call rescfile%finalize()
   end subroutine simulation_final

end module simulation
