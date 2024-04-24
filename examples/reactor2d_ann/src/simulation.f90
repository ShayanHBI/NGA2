!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,           only: WP
   use string,              only: str_medium
   use geometry,            only: cfg,Lx
   use hypre_str_class,     only: hypre_str
   use lowmach_class,       only: lowmach
   use multivdscalar_class, only: multivdscalar
   use ann_class,           only: ann
   use timetracker_class,   only: timetracker
   use ensight_class,       only: ensight
   use event_class,         only: event
   use monitor_class,       only: monitor

   implicit none
   private

   !> Single phase low Mach flow solver, scalar solver, and corresponding time tracker
   type(hypre_str),     public :: ps
   type(lowmach),       public :: fs
   type(multivdscalar), public :: sc
   type(timetracker),   public :: time

   !> Artificial neural networks
   type(ann) :: aen    !< Auto-encoder
   type(ann) :: csn    !< Source terms
   type(ann) :: trn    !< Transport properties
   integer   :: nY_sub !< The number of input species to the auto-encoder

   !> Chemistry
   character(len=str_medium), dimension(:), allocatable :: spec_name

   !> Ensight postprocessing
   type(ensight) :: ens_out
   type(event)   :: ens_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,statfile

   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC,SCtmp,SC_src    !< Scalar solver arrays
   logical , dimension(:,:,:,:), allocatable :: bqflag                !< Flag for bquick scheme
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW,resRHO !< Residuals
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi              !< Cell-centred velocity components
   real(WP), dimension(:,:,:),   allocatable :: SC_init,SC_init_,U_init,V_init !< Initial condition for scalar and velocity fields
   real(WP), dimension(:,:,:),   allocatable :: T                     !< Temperature
   real(WP), dimension(:,:,:),   allocatable :: tmpfield              !< Temporary field for statistics
   real(WP), dimension(:),       allocatable :: Y_sub                 !< Mass fraction of species that are used by networks
   real(WP), dimension(:),       allocatable :: Y_init                !< Initial mass fractions
   real(WP), dimension(:),       allocatable :: TYS                   !< Temperature, mass fractions, and source terms of sub-species
   real(WP), dimension(4)                    :: trnprop               !< Transport properties: Temperature, logarithm of density, viscosity, and scalar diffusivity

   !> Latent variables
   integer :: nlatent

   !> Post process for species mass fractions
   integer :: n_Y                                                     !< Number of output species
   character(len=str_medium), dimension(:), allocatable :: Y_name     !< Names of output species
   real(WP), dimension(:,:,:,:), allocatable            :: Y          !< Output species mass fractions (must exist in Y_sub)
   integer, dimension(:), allocatable                   :: iY_in_sub  !< The indices of post-process species in the sub-species list

   !> Scalar and species indices
   integer :: isc,iY

   !> Statistics
   real(WP) :: rhomean,Tmean,Trms,Tmax

   !> Simulation sub-routines
   public :: simulation_init,simulation_run,simulation_final


contains


   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes the x+ boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator


   !> Function that localizes y- boundary
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc


   !> Function that localizes y+ boundary
   function yp_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator_sc


   !> Function that localizes jet at -x
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the right domain boundary
   function xp_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid),intent(in) :: pg
      integer,intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator_sc


   !> Initialize a double delta scalar field
   subroutine doubledelta_SCinit(Lbu,imin,imax,jmin,jmax,kmin,kmax)
      use precision
      use param,   only: param_read
      use random,  only: random_normal,random_uniform
      use, intrinsic :: iso_c_binding
      implicit none
      ! Buffer region lenght
      real(WP),intent(in) :: Lbu
      ! Bounds
      integer :: imin,imax,jmin,jmax,kmin,kmax
      ! Complex and real buffer
      complex(WP), dimension(:,:,:), pointer :: Cbuf
      real(WP),    dimension(:,:,:), pointer :: Rbuf
      ! Spectrum computation
      real(WP) :: spec_amp,eps,amp_disc,energy_spec
      complex(WP), dimension(:,:,:), pointer :: ak,bk
      ! Fourier coefficients
      integer(KIND=8) :: plan_r2c,plan_c2r
      ! Other
      integer     :: i,j,k,nk,nx,ny,nz
      complex(WP) :: ii=(0.0_WP,1.0_WP)
      real(WP)    :: rand,pi,ke,dk,kc,ks,ksk0ratio,kcksratio,kx,ky,kz,kk,f_phi,kk2
      include 'fftw3.f03'

      ! Create pi
      pi=acos(-1.0_WP)

      ! Step size for wave number
      dk=2.0_WP*pi/(Lx-2.0_WP*Lbu)

      ! Number of cells iniside the initialization region
      nx=imax-imin+1
      ny=jmax-jmin+1
      nz=kmax-kmin+1
      nk=nx/2+1

      ! Initialize in similar manner to Eswaran and Pope 1988
      call param_read('ks/ko',ksk0ratio)
      ks=ksk0ratio*dk
      call param_read('kc/ks',kcksratio)
      kc=kcksratio*ks

      ! Allocate Cbuf and Rbuf
      allocate(Cbuf(nk,ny,nz))
      allocate(Rbuf(nx,ny,nz))

      ! Compute the Fourier coefficients
      do k=1,nz
         do j=1,ny
            do i=1,nk

               ! Wavenumbers
               kx=real(i-1,WP)*dk
               ky=real(j-1,WP)*dk
               if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
               kz=real(k-1,WP)*dk
               if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
               kk =sqrt(kx**2+ky**2+kz**2)
               kk2=sqrt(kx**2+ky**2)

               ! Compute the Fourier coefficients
               if ((ks-dk/2.0_WP.le.kk).and.(kk.le.ks+dk/2.0_WP)) then
                  f_phi=1.0_WP
               else
                  f_phi=0.0_WP
               end if
               call random_number(rand)
               if (kk.lt.1e-10) then
                  Cbuf(i,j,k)=0.0_WP
               else
                  Cbuf(i,j,k)=sqrt(f_phi/(4.0_WP*pi*kk**2))*exp(ii*2.0_WP*pi*rand)
               end if

            end do
         end do
      end do

      ! Oddball and setting up plans based on geometry
      do j=nk+1,ny
         Cbuf(1,j,1)=conjg(Cbuf(1,ny+2-j,1))
      end do
      call dfftw_plan_dft_c2r_2d(plan_c2r,nx,ny,Cbuf,Rbuf,FFTW_ESTIMATE)
      call dfftw_plan_dft_r2c_2d(plan_r2c,nx,ny,Rbuf,Cbuf,FFTW_ESTIMATE)

      ! Inverse Fourier transform
      call dfftw_execute(plan_c2r)

      ! Force 'double-delta' pdf on scalar field
      do k=1,nz
         do j=1,ny
            do i=1,nx
               if (Rbuf(i,j,k).le.0.0_WP) then
                  Rbuf(i,j,k)=0.0_WP
               else
                  Rbuf(i,j,k)=1.0_WP
               end if
            end do
         end do
      end do

      ! Fourier Transform and filter to smooth
      call dfftw_execute(plan_r2c)

      do k=1,nz
         do j=1,ny
            do i=1,nk

               ! Wavenumbers
               kx=real(i-1,WP)*dk
               ky=real(j-1,WP)*dk
               if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
               kz=real(k-1,WP)*dk
               if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
               kk =sqrt(kx**2+ky**2+kz**2)
               kk2=sqrt(kx**2+ky**2)

               ! Filter to remove high wavenumber components
               if (kk.le.kc) then
                  Cbuf(i,j,k)=Cbuf(i,j,k)*1.0_WP
               else
                  Cbuf(i,j,k)=Cbuf(i,j,k)*(kc/kk)**2
               end if

            end do
         end do
      end do

      ! Oddball
      do j=nk+1,ny
         Cbuf(1,j,1)=conjg(Cbuf(1,ny+2-j,1))
      end do

      ! Fourier Transform back to real
      call dfftw_execute(plan_c2r)

      ! Set zero value for buffer region
      SC_init=0.0_WP
      ! Set the internal scalar field
      SC_init(imin:imax,jmin:jmax,kmin:kmax)=Rbuf/real(nx*ny*nz,WP)

      ! Destroy the plans
      call dfftw_destroy_plan(plan_c2r)
      call dfftw_destroy_plan(plan_r2c)

      ! Clean up
      deallocate(Cbuf)
      deallocate(Rbuf)
   end subroutine doubledelta_SCinit


   !> Initialize PP spectrum for velocity
   subroutine PP_spectrum(Lbu,Lfd,Ut,le)
      use precision
      use param,    only: param_read
      use random,   only: random_normal,random_uniform
      use, intrinsic :: iso_c_binding
      implicit none
      ! Buffer and faded region lenght
      real(WP), intent(in) :: Lbu,Lfd
      ! Turbulent velocity
      real(WP) :: Ut
      ! Spectrum type
      real(WP) :: le
      ! Spectrum computation
      real(WP) :: psr,ps1,ps2,ke,dk,kc,kk,kx,ky,kz,kk2
      real(WP) :: spec_amp,eps,amp_disc,energy_spec
      complex(WP), dimension(:,:,:), pointer :: ak,bk
      ! Cutoff wave number
      integer  :: nk
      ! Complex buffer
      complex(WP), dimension(:,:,:), pointer :: Cbuf
      ! Real buffer
      real(WP), dimension(:,:,:), pointer :: Rbuf
      ! Other
      integer :: i,j,k
      integer :: imin,imax,jmin,jmax,kmin,kmax,nx,ny,nz
      complex(WP) :: ii=(0.0_WP,1.0_WP)
      real(WP) :: rand,pi
      ! Fourier coefficients
      integer(KIND=8) :: plan_r2c,plan_c2r
      complex(WP), dimension(:,:,:), pointer :: Uk,Vk
      include 'fftw3.f03'

      ! Create pi
      pi=acos(-1.0_WP)

      ! Find bounds of the region to be initialized
      call get_borders(Lbu,imin,imax,jmin,jmax,kmin,kmax)

      ! Number of cells iniside the initialization region
      nx=imax-imin+1
      ny=jmax-jmin+1
      nz=kmax-kmin+1
      nk=nx/2+1

      ! Spectrum computation
      ke=2.0_WP*pi/le
      dk=2.0_WP*pi/(Lx-2.0_WP*Lbu)
      kc=real(nx/2,WP)*dk

      eps=ke/1000000.0_WP
      spec_amp=(32.0_WP/3.0_WP)*sqrt(2.0_WP/pi)*Ut**2/ke
      amp_disc=sqrt(dk)**3

      ! Compute spectrum
      allocate(ak(nk,ny,nz),bk(nk,ny,nz))
      do k=1,nz
         do j=1,ny
            do i=1,nk
               ! Random numbers
               call random_number(rand)
               psr=2.0_WP*pi*(rand-0.5_WP)
               call random_number(rand)
               ps1=2.0_WP*pi*(rand-0.5_WP)
               call random_number(rand)
               ps2=2.0_WP*pi*(rand-0.5_WP)
               ! Wavenumbers
               kx=real(i-1,WP)*dk
               ky=real(j-1,WP)*dk
               if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
               kz=real(k-1,WP)*dk
               if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
               kk=sqrt(kx**2+ky**2+kz**2)
               ! Spectrums
               energy_spec=spec_amp*(kk/ke)**4*exp(-2.0_WP*(kk/ke)**2)
               ! Coeff
               ak(i,j,k)=0.0_WP
               bk(i,j,k)=0.0_WP
               if ((kk.gt.eps).and.(kk.le.kc)) then
                  ak(i,j,k)=dk*sqrt(energy_spec/(1.0_WP*pi*kk**1))*exp(ii*ps1)
               end if
            end do
         end do
      end do

      ! Compute 3D field
      allocate(Uk(nk,ny,nz))
      allocate(Vk(nk,ny,nz))
      Uk=(0.0_WP,0.0_WP)
      Vk=(0.0_WP,0.0_WP)

      ! Compute the Fourier coefficients
      do k=1,nz
         do j=1,ny
            do i=1,nk
               ! Wavenumbers
               kx=real(i-1,WP)*dk
               ky=real(j-1,WP)*dk
               if (j.gt.nk) ky=-real(nx+1-j,WP)*dk
               kz=real(k-1,WP)*dk
               if (k.gt.nk) kz=-real(nx+1-k,WP)*dk
               kk =sqrt(kx**2+ky**2+kz**2)
               kk2=sqrt(kx**2+ky**2)

               if ((kk.gt.eps).and.(kk.le.kc)) then
                  if (kk2.lt.eps) then
                     Uk(i,j,k)=(ak(i,j,k)+bk(i,j,k))/sqrt(2.0_WP)
                     Vk(i,j,k)=(bk(i,j,k)-ak(i,j,k))/sqrt(2.0_WP)
                  else
                     Uk(i,j,k)=(ak(i,j,k)*kk*ky+bk(i,j,k)*kx*kz)/(kk*kk2)
                     Vk(i,j,k)=(bk(i,j,k)*ky*kz-ak(i,j,k)*kk*kx)/(kk*kk2)
                  end if
               end if
            end do
         end do
      end do

      ! Oddball
      do j=nk+1,ny
         Uk(1,j,1)=conjg(Uk(1,ny+2-j,1))
         Vk(1,j,1)=conjg(Vk(1,ny+2-j,1))
      end do

      ! Inverse Fourier transform
      allocate(Cbuf(nk,ny,nz))
      allocate(Rbuf(nx,ny,nz))
      call dfftw_plan_dft_c2r_2d(plan_c2r,nx,ny,Cbuf,Rbuf,FFTW_ESTIMATE)
      call dfftw_plan_dft_r2c_2d(plan_r2c,nx,ny,Rbuf,Cbuf,FFTW_ESTIMATE)

      ! Set zero in the buffer region
      U_init=0.0_WP
      V_init=0.0_WP

      ! Execute the plans
      Cbuf=Uk
      call dfftw_execute(plan_c2r)
      U_init(imin:imax,jmin:jmax,kmin:kmax)=Rbuf
      Cbuf=Vk
      call dfftw_execute(plan_c2r)
      V_init(imin:imax,jmin:jmax,kmin:kmax)=Rbuf

      ! Clean up
      deallocate(Uk)
      deallocate(Vk)
      deallocate(ak)
      deallocate(bk)

      ! Fade to zero in the buffer region
      ! call fade_borders(U_init,Lbu,Lfd,imin,imax,jmin,jmax,kmin,kmax)
      ! call fade_borders(V_init,Lbu,Lfd,imin,imax,jmin,jmax,kmin,kmax)
   end subroutine PP_spectrum


   !> Find the border indices of the initialization region
   subroutine get_borders(Lbu,imin,imax,jmin,jmax,kmin,kmax)
      implicit none
      real(WP),intent(in) :: Lbu
      integer,intent(out) :: imin,imax,jmin,jmax,kmin,kmax
      integer :: i,j,k

      ! Find the x bounds
      imin=cfg%imin
      do i=cfg%imin,cfg%imax
         if (cfg%xm(i).gt.cfg%x(cfg%imin)+Lbu) then
            imin=i
            exit
         end if
      end do
      imax=cfg%imax
      do i=cfg%imax,cfg%imin,-1
         if (cfg%xm(i).lt.cfg%x(cfg%imax+1)-Lbu) then
            imax=i
            exit
         end if
      end do

      ! Find the y bounds
      jmin=cfg%jmin
      do j=cfg%jmin,cfg%jmax
         if (cfg%ym(j).gt.cfg%y(cfg%jmin)+Lbu) then
            jmin=j
            exit
         end if
      end do
      jmax=cfg%jmax
      do j=cfg%jmax,cfg%jmin,-1
         if (cfg%ym(j).lt.cfg%y(cfg%jmax+1)-Lbu) then
            jmax=j
            exit
         end if
      end do

      ! Find the z bounds
      kmin=cfg%kmin
      do k=cfg%kmin,cfg%kmax
         if (cfg%zm(k).gt.cfg%z(cfg%kmin)+Lbu) then
            kmin=k
            exit
         end if
      end do
      kmax=cfg%kmax
      do k=cfg%kmax,cfg%kmin,-1
         if (cfg%zm(k).lt.cfg%z(cfg%kmax+1)-Lbu) then
            kmax=k
            exit
         end if
      end do
   end subroutine get_borders


   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read,param_getsize,param_exists
      use messager, only: die
      use fcmech, only: nspec
      implicit none
      real(WP) :: L_buffer


      ! Read-in the inputs
      call param_read('Buffer region length',L_buffer)
      n_Y=param_getsize('Ensight output species')

      
      ! Create a low-Mach flow solver with bconds
      create_velocity_solver: block
         use hypre_str_class, only: pcg_pfmg,smg
         use lowmach_class,   only: dirichlet,neumann
         real(WP) :: visc
         ! Create flow solver
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Boundary conditions
         call fs%add_bcond(name='ym_outflow',type=neumann,face='y',dir=-1,canCorrect=.True.,locator=ym_locator)
         call fs%add_bcond(name='yp_outflow',type=neumann,face='y',dir=+1,canCorrect=.True.,locator=yp_locator)
         call fs%add_bcond(name='xm_outflow',type=neumann,face='x',dir=-1,canCorrect=.True.,locator=xm_locator)
         call fs%add_bcond(name='xp_outflow',type=neumann,face='x',dir=+1,canCorrect=.True.,locator=xp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg,nst=7)
         ps%maxlevel=12
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Setup the flow solver
         call fs%setup(pressure_solver=ps)
      end block create_velocity_solver


      ! Create neural networks
      create_ann: block
         use string, only: str_medium
         character(len=str_medium) :: netmodel,aenfname,csnfname,trnfname,netpath
         ! Read-in the data file names
         netpath='networks/'
         if (param_exists('Network model')) then
            call param_read('Network model',netmodel)
            netpath=trim(netpath)//trim(netmodel)//'/'
         end if
         call param_read('Auto encoder'        ,aenfname)
         call param_read('Chemical source'     ,csnfname)
         call param_read('Transport properties',trnfname)
         ! The auto encoder network
         aen=ann(cfg=cfg,fdata=trim(netpath)//trim(aenfname),name='Auto encoder network')
         call aen%print()
         ! Species sub-array and latent space variables sizes
         nY_sub=size(aen%inp_sub_ind)
         nlatent=size(aen%encoder_weight,dim=2)
         ! The source terms network
         csn=ann(cfg=cfg,fdata=trim(netpath)//trim(csnfname),name='Chemical source network')
         call csn%print()
         ! The transport properties network
         trn=ann(cfg=cfg,fdata=trim(netpath)//trim(trnfname),name='Transport properties network')
         call trn%print()
      end block create_ann


      ! Create a multi scalar solver
      create_scalar_solver: block
         use multivdscalar_class, only: dirichlet,neumann,quick,bquick
         ! Create scalar solver
         sc=multivdscalar(cfg=cfg,scheme=bquick,nscalar=nlatent,name='Variable density multi scalar')
         ! Boundary conditions
         call sc%add_bcond(name='ym_outflow',type=neumann,locator=ym_locator_sc,dir='-y')
         call sc%add_bcond(name='yp_outflow',type=neumann,locator=yp_locator_sc,dir='+y')
         call sc%add_bcond(name='xm_outflow',type=neumann,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp_outflow',type=neumann,locator=xp_locator_sc,dir='+x')
         ! Setup the solver
         call sc%setup()
      end block create_scalar_solver


      ! Allocate work arrays
      allocate_work_arrays: block
         ! Flow solver
         allocate(resU  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resV  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resW  (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(resRHO(fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Ui    (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Vi    (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         allocate(Wi    (fs%cfg%imino_:fs%cfg%imaxo_,fs%cfg%jmino_:fs%cfg%jmaxo_,fs%cfg%kmino_:fs%cfg%kmaxo_))
         ! Scalar solver
         allocate(resSC (sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_,sc%nscalar))
         allocate(SCtmp (sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_,sc%nscalar))
         allocate(bqflag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,sc%nscalar))
         allocate(SC_src(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_,nlatent)); SC_src=0.0_WP
         ! Temporary arrays
         allocate(SC_init (sc%cfg%imin:sc%cfg%imax,sc%cfg%jmin:sc%cfg%jmax,sc%cfg%kmin:sc%cfg%kmax))
         allocate(SC_init_(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
         allocate(U_init  (fs%cfg%imin:fs%cfg%imax,fs%cfg%jmin:fs%cfg%jmax,fs%cfg%kmin:fs%cfg%kmax))
         allocate(V_init  (fs%cfg%imin:fs%cfg%imax,fs%cfg%jmin:fs%cfg%jmax,fs%cfg%kmin:fs%cfg%kmax))
         allocate(tmpfield(sc%cfg%imin:sc%cfg%imax,sc%cfg%jmin:sc%cfg%jmax,sc%cfg%kmin:sc%cfg%kmax)); tmpfield=0.0_WP
         ! Combustion
         allocate(T(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
         allocate(Y(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_,n_Y))
         allocate(Y_sub(nY_sub))  ; Y_sub=0.0_WP
         allocate(TYS(2*nY_sub+1)); TYS=0.0_WP
         allocate(iY_in_sub(n_Y)) ; iY_in_sub=0
         allocate(Y_init(nspec))  ; Y_init=0.0_WP
         allocate(spec_name(nspec))
         allocate(Y_name(n_Y))
      end block allocate_work_arrays


      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=fs%cfg%amRoot,name='reactor2d_ann')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Inner iterations',time%itmax)
         time%dt=time%dtmax
      end block initialize_timetracker


      ! Initialize species
      initialize_species: block
         use string, only: str_long
         use, intrinsic :: iso_fortran_env, only: output_unit
         integer  :: iY_sub,ispec
         character(len=str_long) :: errmsg

         ! Get all the species names in the mechanism
         call fcmech_get_speciesnames(spec_name)
         ! Read-in the names of post processed species
         call param_read('Ensight output species',Y_name)

         ! Print species information
         if (cfg%amRoot) then

            ! Species info
            write(output_unit,'(" >  Species information")')

            ! Sub-array used in ANN
            write(output_unit,'(" Sub-array species used in the ann:")',advance='No')
            do iY_sub=1,nY_sub
               write(output_unit,'(a)',advance='No') ' '//trim(spec_name(aen%inp_sub_ind(iY_sub)))
            end do
            write(output_unit,'(" ")')

            ! Post processed species
            write(output_unit,'(" Post processed species:")',advance='No')
            do iY=1,n_Y
               write(output_unit,'(a)',advance='No') ' '//trim(Y_name(iY))
            end do
            write(output_unit,'(" ")')

         end if

         ! Get the initial mass fractions
         do ispec=1,nspec
            ! Initial values
            if (param_exists('Initial '//trim(spec_name(ispec)))) then
               call param_read('Initial '//trim(spec_name(ispec)),Y_init(ispec))
            end if
         end do

         ! Process the species indices and set initial mass fractions
         if (cfg%amRoot) write(output_unit,'(" Initial mass fractions:")')
         do iY_sub=1,nY_sub
            ! Global species index
            ispec=aen%inp_sub_ind(iY_sub)
            Y_sub(iY_sub)=Y_init(ispec)
            if (cfg%amRoot) then
               write(output_unit,'(a10,": "es12.5)') trim(spec_name(ispec)),Y_sub(iY_sub)
            end if
         end do

         ! Localize post processed species inside the ann species sub-array
         do iY=1,n_Y
            do iY_sub=1,nY_sub
               ! Global species index
               ispec=aen%inp_sub_ind(iY_sub)
               ! Local species index in the sub-array species
               if (trim(spec_name(ispec)).eq.trim(Y_name(iY))) then
                  iY_in_sub(iY)=iY_sub
                  exit
               end if
            end do
         end do

         ! Check if we have found all the required species
         do iY=1,n_Y
            if (iY_in_sub(iY).eq.0) then
               errmsg='Could not find '//trim(Y_name(iY))//' in the sub-species list. Make sure that '//trim(Y_name(iY))//' exists in the sub-array species used by the ANN.'
               call die(trim(errmsg))
            end if
         end do
      end block initialize_species


      ! Initialize our scalar fields
      initialize_scalar: block
         use messager, only: die
         use parallel, only: MPI_REAL_WP
         use fcmech,   only: hsp
         integer  :: i,j,k,ierr,ispec
         integer  :: imin,imax,jmin,jmax,kmin,kmax
         real(WP) :: h_init,T_buf,T_range,T_min
         real(WP) :: SC_init_min,SC_init_max

         ! Read-in inputs
         call param_read('Buffer temperature',T_buf)
         call param_read('Temperature range',T_range)
         call param_read('Temperature lower bound',T_min)
         ! Get bounds of the region to be initialized
         call get_borders(L_buffer,imin,imax,jmin,jmax,kmin,kmax)
         ! Initialize the global double delta scalar field
         if (sc%cfg%amRoot) call doubledelta_SCinit(L_buffer,imin,imax,jmin,jmax,kmin,kmax)
         ! Communicate information
         call MPI_BCAST(SC_init,sc%cfg%nx*sc%cfg%ny*sc%cfg%nz,MPI_REAL_WP,0,sc%cfg%comm,ierr)
         ! Get the min and max
         SC_init_min=minval(SC_init)
         SC_init_max=maxval(SC_init)
         SC_init_=SC_init(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_)

         ! Set initial conditions
         do k=sc%cfg%kmino_, sc%cfg%kmaxo_
            do j=sc%cfg%jmino_, sc%cfg%jmaxo_
               do i=sc%cfg%imino_, sc%cfg%imaxo_

                  ! Temperature and enthalpy
                  if ((i.ge.imin).and.(i.le.imax).and.(j.ge.jmin).and.(j.le.jmax).and.(k.ge.kmin).and.(k.le.kmax)) then
                     T(i,j,k)=(SC_init(i,j,k)-SC_init_min)/(SC_init_max-SC_init_min)*T_range+T_min
                  else
                     T(i,j,k)=T_buf
                  end if
                  call fcmech_thermodata(T(i,j,k))
                  h_init=0.0_WP
                  do ispec=1, nspec
                     h_init=h_init+hsp(ispec)*Y_init(ispec)
                  end do
                  
                  ! Map Y and h to the neural network scalars
                  call aen%encode([h_init,Y_sub],sc%SC(i,j,k,:))

                  ! Get transport properties
                  call trn%forward_pass(sc%SC(i,j,k,:),trnprop)
                  sc%rho(i,j,k)   =exp(trnprop(2))
                  fs%visc(i,j,k)  =trnprop(3)
                  sc%diff(i,j,k,:)=trnprop(4)*sc%rho(i,j,k)

                  ! Map the neural network scalars to Y (for visualization purposes at t=0)
                  call aen%forward_pass(sc%SC(i,j,k,:),TYS)
                  do iY=1,n_Y
                     Y(i,j,k,iY)=TYS(iY_in_sub(iY)+1)
                  end do

               end do
            end do
         end do

         ! Sync fields
         call sc%cfg%sync(T)
         call sc%cfg%sync(sc%rho)
         call sc%cfg%sync(fs%visc)
         do isc=1,sc%nscalar
            call sc%cfg%sync(sc%SC  (:,:,:,isc))
            call sc%cfg%sync(sc%diff(:,:,:,isc))
         end do
         do iY=1,n_Y
            call sc%cfg%sync(Y(:,:,:,iY))
         end do

         ! Release unused memory
         deallocate(SC_init)
         deallocate(Y_init)
      end block initialize_scalar


      ! Initialize our velocity field
      initialize_velocity: block
         use lowmach_class, only: bcond
         use parallel,      only: MPI_REAL_WP
         integer  :: n,i,j,k,ierr
         real(WP) :: Ut,le
         type(bcond), pointer :: mybc
         ! Read-in the inputs
         call param_read('Velocity fluctuation',Ut)
         call param_read('Energetic scale',le)
         ! Initialize the global velocity field
         if (fs%cfg%amRoot) call PP_spectrum(L_buffer,0.0_WP,Ut,le)
         ! Communicate information
         call MPI_BCAST(U_init,fs%cfg%nx*fs%cfg%ny*fs%cfg%nz,MPI_REAL_WP,0,fs%cfg%comm,ierr)
         call MPI_BCAST(V_init,fs%cfg%nx*fs%cfg%ny*fs%cfg%nz,MPI_REAL_WP,0,fs%cfg%comm,ierr)
         ! Set the local velocity field
         fs%U(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)=U_init(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)
         fs%V(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)=V_init(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)
         ! Sync it
         call fs%cfg%sync(fs%U)
         call fs%cfg%sync(fs%V)
         ! Release memory
         deallocate(U_init)
         deallocate(V_init)
         ! Set density from scalar
         fs%rho=sc%rho
         ! Form momentum
         call fs%rho_multiply()
         ! Apply all other boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         call fs%interp_vel(Ui,Vi,Wi)
         resRHO=0.0_WP
         call fs%get_div(drhodt=resRHO)
         ! Compute MFR through all boundary conditions
         call fs%get_mfr()
      end block initialize_velocity


      ! Initialize turbulence statistics
      initialize_stats: block
         ! Mean density
         call cfg%integrate(sc%rho,rhomean)
         rhomean=rhomean/cfg%vol_total
         ! Favre-averaged temperature
         tmpfield=sc%rho*T
         call cfg%integrate(tmpfield,Tmean)
         Tmean=Tmean/(rhomean*cfg%vol_total)
         ! Favre-averaged temperature fluctuations
         tmpfield=sc%rho*(T-Tmean)**2
         call cfg%integrate(tmpfield,Trms)
         Trms=sqrt(Trms/(rhomean*cfg%vol_total))
         ! Maximum temperature
         call cfg%maximum(T,Tmax)
      end block initialize_stats

      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='reactor2d_ann')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_scalar('pressure'  ,fs%P)
         call ens_out%add_vector('velocity'  ,Ui,Vi,Wi)
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('density'   ,sc%rho)
         call ens_out%add_scalar('viscosity' ,fs%visc)
         call ens_out%add_scalar('T',T)
         call ens_out%add_scalar('thermal_diff',sc%diff(:,:,:,1))
         call ens_out%add_scalar('sc',SC_init_(sc%cfg%imino_:sc%cfg%imaxo_,sc%cfg%jmino_:sc%cfg%jmaxo_,sc%cfg%kmino_:sc%cfg%kmaxo_))
         do iY=1,n_Y
            call ens_out%add_scalar('Y_'//Y_name(iY),Y(:,:,:,iY))
         end do
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         ! Create simulation monitor
         mfile=monitor(fs%cfg%amRoot,'simulation')
         call mfile%add_column(time%n,'Timestep number')
         call mfile%add_column(time%t,'Time')
         call mfile%add_column(time%dt,'Timestep size')
         call mfile%add_column(time%cfl,'Maximum CFL')
         call mfile%add_column(fs%Umax,'Umax')
         call mfile%add_column(fs%Vmax,'Vmax')
         call mfile%add_column(fs%Wmax,'Wmax')
         call mfile%add_column(fs%Pmax,'Pmax')
         call mfile%add_column(sc%rhomax,'RHOmax')
         call mfile%add_column(sc%rhomin,'RHOmin')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Creat statistics monitor
         statfile=monitor(fs%cfg%amRoot,'stats')
         call statfile%add_column(time%t,'Time')
         call statfile%add_column(Tmean,'Tmean')
         call statfile%add_column(Trms,'Trms')
         call statfile%add_column(Tmax,'Tmax')
         call statfile%write()
      end block create_monitor

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      integer :: i,j,k


      ! Perform time integration
      do while (.not. time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old scalar
         sc%rhoold=sc%rho
         sc%SCold =sc%SC

         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW

         ! Get the source terms
         do k=sc%cfg%kmino_,sc%cfg%kmaxo_
            do j=sc%cfg%jmino_,sc%cfg%jmaxo_
               do i=sc%cfg%imino_,sc%cfg%imaxo_
                  call csn%forward_pass(sc%SC(i,j,k,:),SC_src(i,j,k,:))
               end do
            end do
         end do

         ! Perform sub-iterations
         do while (time%it.le.time%itmax)

            ! =============SCALAR SOLVER =============!
            
            ! Reset metric for bquick
            call sc%metric_reset()

            ! Build mid-time scalar
            sc%SC=0.5_WP*(sc%SC+sc%SCold)

            ! Explicit calculation of drhoSC/dt from scalar equation
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Assemble explicit residual
            do isc=1,sc%nscalar
               resSC(:,:,:,isc)=time%dt*resSC(:,:,:,isc)-2.0_WP*sc%rho*sc%SC(:,:,:,isc)+(sc%rho+sc%rhoold)*sc%SCold(:,:,:,isc)+time%dt*sc%rho*SC_src(:,:,:,isc)
               SCtmp(:,:,:,isc)=2.0_WP*sc%SC(:,:,:,isc)-sc%SCold(:,:,:,isc)+resSC(:,:,:,isc)/sc%rho
            end do

            ! Apply it to get explicit scalar prediction
            do isc=1,sc%nscalar
               do k=sc%cfg%kmino_,sc%cfg%kmaxo_
                  do j=sc%cfg%jmino_,sc%cfg%jmaxo_
                     do i=sc%cfg%imino_,sc%cfg%imaxo_
                        if ((SCtmp(i,j,k,isc).le.0.0_WP).or.(SCtmp(i,j,k,isc).ge.1.0_WP)) then
                           bqflag(i,j,k,isc)=.true.
                        else
                           bqflag(i,j,k,isc)=.false.
                        end if
                     end do
                  end do
               end do
            end do

            ! Adjust metrics
            call sc%metric_adjust(SCtmp,bqflag)

            ! Recompute drhoSC/dt
            call sc%get_drhoSCdt(resSC,fs%rhoU,fs%rhoV,fs%rhoW)

            ! Assemble explicit residual
            do isc=1,sc%nscalar
               resSC(:,:,:,isc)=(time%dt*resSC(:,:,:,isc)+sc%rhoold*sc%SCold(:,:,:,isc))/sc%rho-2.0_WP*sc%SC(:,:,:,isc)+sc%SCold(:,:,:,isc)+time%dt*SC_src(:,:,:,isc)
            end do

            ! Apply these residuals
            sc%SC=2.0_WP*sc%SC-sc%SCold+resSC

            ! Apply boundary conditions on the resulting field
            call sc%apply_bcond(time%t,time%dt)

            ! =============================================

            ! ============UPDATE PROPERTIES ====================

            ! Get transport properties
            do k=sc%cfg%kmino_,sc%cfg%kmaxo_
               do j=sc%cfg%jmino_,sc%cfg%jmaxo_
                  do i=sc%cfg%imino_,sc%cfg%imaxo_
                     call trn%forward_pass(sc%SC(i,j,k,:),trnprop)
                     sc%rho(i,j,k)   =exp(trnprop(2))
                     fs%visc(i,j,k)  =trnprop(3)
                     sc%diff(i,j,k,:)=trnprop(4)*sc%rho(i,j,k)
                  end do
               end do
            end do

            ! ===================================================

            ! ============VELOCITY SOLVER ======================

            ! Build n+1 density
            fs%rho=0.5_WP*(sc%rho+sc%rhoold)

            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)

            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)

            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! Divide by density to get velocity residuals
            do k=cfg%kmin_,cfg%kmax_+1
               do j=cfg%jmin_,cfg%jmax_+1
                  do i=cfg%imin_,cfg%imax_+1
                     resU(i,j,k)=resU(i,j,k)/sum(fs%itpr_x(:,i,j,k)*fs%rho(i-1:i,j,k))
                     resV(i,j,k)=resV(i,j,k)/sum(fs%itpr_y(:,i,j,k)*fs%rho(i,j-1:j,k))
                     resW(i,j,k)=resW(i,j,k)/sum(fs%itpr_z(:,i,j,k)*fs%rho(i,j,k-1:k))
                  end do
               end do
            end do
            call fs%cfg%sync(resU)
            call fs%cfg%sync(resV)
            call fs%cfg%sync(resW)

            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW

            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            call fs%rho_multiply()

            ! Solve Poisson equation
            call sc%get_drhodt(dt=time%dt, drhodt=resRHO)
            call fs%correct_mfr(drhodt=resRHO)
            call fs%get_div(drhodt=resRHO)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)

            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            call fs%rho_divide

            ! ===================================================

            ! Increment sub-iteration counter
            time%it=time%it+1

         end do

         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call sc%get_drhodt(dt=time%dt,drhodt=resRHO)
         call fs%get_div(drhodt=resRHO)

         ! Post process
         if (ens_evt%occurs()) then
            ! Map the neural network scalars to T and Y
            do k=sc%cfg%kmino_,sc%cfg%kmaxo_
               do j=sc%cfg%jmino_,sc%cfg%jmaxo_
                  do i=sc%cfg%imino_,sc%cfg%imaxo_
                     call aen%forward_pass(sc%SC(i,j,k,:),TYS)
                     T(i,j,k)=TYS(1)
                     do iY=1,n_Y
                        Y(i,j,k,iY)=TYS(iY_in_sub(iY)+1)
                     end do
                  end do
               end do
            end do
            ! Mean density
            call cfg%integrate(sc%rho,rhomean)
            rhomean=rhomean/cfg%vol_total
            ! Favre-averaged temperature
            tmpfield=sc%rho*T
            call cfg%integrate(tmpfield,Tmean)
            Tmean=Tmean/(rhomean*cfg%vol_total)
            ! Favre-averaged temperature fluctuations
            tmpfield=sc%rho*(T-Tmean)**2
            call cfg%integrate(tmpfield,Trms)
            Trms=sqrt(Trms/(rhomean*cfg%vol_total))
            ! Maximum temperature
            call cfg%maximum(T,Tmax)
            ! Output to ensight
            call ens_out%write_data(time%t)
            ! Output monitoring
            call statfile%write()
         end if

         ! Perform and output monitoring
         call fs%get_max()
         call sc%get_max()
         call sc%get_int()
         call mfile%write()
         call cflfile%write()

      end do


   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none


      ! Get rid of all objects-need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker

      ! Deallocate work arrays
      deallocate(spec_name,Y_sub,TYS,iY_in_sub,Y_name)
      deallocate(resU,resV,resW,resRHO,Ui,Vi,Wi,T)
      deallocate(resSC,SCtmp,SC_src,bqflag,Y)


   end subroutine simulation_final


end module simulation
