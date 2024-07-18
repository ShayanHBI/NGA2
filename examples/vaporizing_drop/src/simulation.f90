!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs
   use tpscalar_class,    only: tpscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   implicit none
   private
   
   !> Get a couple linear solvers, a two-phase flow solver and volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps
   !type(ddadi),       public :: vs
   type(ddadi),       public :: ss
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(timetracker), public :: time,pseudo_time
   
   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out,ens_mflux_out
   type(event)    :: ens_evt,ens_mflux_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,mfluxfile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: VFgradX,VFgradY,VFgradZ
   real(WP), dimension(:,:,:),   allocatable :: mflux
   real(WP), dimension(:,:,:),   allocatable :: mfluxL,mfluxL_old,resmfluxL,mfluxL_err_field
   real(WP), dimension(:,:,:),   allocatable :: mfluxG,mfluxG_old,resmfluxG
   
   !> Problem definition
   real(WP), dimension(3) :: center
   real(WP) :: radius,depth
   integer  :: iTl,iTg
   real(WP) :: mflux_int,mflux_err,mflux_tol
   real(WP) :: mfluxL_int,mfluxL_err,mfluxL_int_err
   real(WP) :: mfluxG_int,mfluxG_err,mfluxG_int_err
   real(WP) :: mflux_ens_time
   real(WP) :: evp_mass_flux
   
contains


   !> Function that defines a level set function for a drop problem
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=radius-sqrt(sum((xyz-center)**2))
      ! Add the pool
      ! G=max(G,depth-xyz(2))
   end function levelset_drop
   
   
   !> Initialization of problem solver
   subroutine simulation_init
      use param, only: param_read
      implicit none
      
      
      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:2))
         allocate(VFgradX   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradX=0.0_WP
         allocate(VFgradY   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradY=0.0_WP
         allocate(VFgradZ   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); VFgradZ=0.0_WP
         allocate(mflux     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mflux =0.0_WP
         allocate(mfluxL    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mfluxL=0.0_WP
         allocate(mfluxG    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mfluxG=0.0_WP
         allocate(mfluxL_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mfluxL_old=0.0_WP
         allocate(mfluxG_old(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mfluxG_old=0.0_WP
         allocate(resmfluxL (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resmfluxL=0.0_WP
         allocate(resmfluxG (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resmfluxG=0.0_WP
         allocate(mfluxL_err_field(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); mfluxL_err_field=0.0_WP
         allocate(resU(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block allocate_work_arrays
      
      
      ! Initialize time tracker with 2 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      
      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,VFhi,VFlo,remap,flux_storage,neumann
         integer :: i,j,k,n,si,sj,sk
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3) :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,name='VOF')
         ! call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=remap,name='VOF')
         vf%cons_correct=.false.
         ! Initialize to a droplet and a pool
         !center=[0.0_WP,0.05_WP,0.0_WP]
         call param_read('Droplet center',center)
         call param_read('Droplet diameter',radius); radius=radius/2.0_WP
         depth =0.02_WP
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  ! Set cube vertices
                  n=0
                  do sk=0,1
                     do sj=0,1
                        do si=0,1
                           n=n+1; cube_vertex(:,n)=[vf%cfg%x(i+si),vf%cfg%y(j+sj),vf%cfg%z(k+sk)]
                        end do
                     end do
                  end do
                  ! Call adaptive refinement code to get volume and barycenters recursively
                  vol=0.0_WP; area=0.0_WP; v_cent=0.0_WP; a_cent=0.0_WP
                  call cube_refine_vol(cube_vertex,vol,area,v_cent,a_cent,levelset_drop,0.0_WP,amr_ref_lvl)
                  vf%VF(i,j,k)=vol/vf%cfg%vol(i,j,k)
                  if (vf%VF(i,j,k).ge.VFlo.and.vf%VF(i,j,k).le.VFhi) then
                     vf%Lbary(:,i,j,k)=v_cent
                     vf%Gbary(:,i,j,k)=([vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]-vf%VF(i,j,k)*vf%Lbary(:,i,j,k))/(1.0_WP-vf%VF(i,j,k))
                  else
                     vf%Lbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                     vf%Gbary(:,i,j,k)=[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)]
                  end if
               end do
            end do
         end do
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Create discontinuous polygon mesh from IRL interface
         call vf%polygonalize_interface()
         ! Calculate distance from polygons
         call vf%distance_from_polygon()
         ! Calculate subcell phasic volumes
         call vf%subcell_vol()
         ! Calculate curvature
         call vf%get_curvature()
         ! Reset moments to guarantee compatibility with interface reconstruction
         call vf%reset_volume_moments()
      end block create_and_initialize_vof
      
      
      ! Create a two-phase flow solver without bconds
      create_and_initialize_flow_solver: block
         use hypre_str_class, only: pcg_pfmg2
         use mathtools,       only: Pi
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         call param_read('Liquid density',fs%rho_l)
         call param_read('Gas density',fs%rho_g)
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         call param_read('Static contact angle',fs%contact_angle)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_pfmg2,nst=7)
         ps%maxlevel=10
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         ! Configure implicit velocity solver
         !vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! Setup the solver
         call fs%setup(pressure_solver=ps)!,implicit_solver=vs)
         ! Zero initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver
      
      
      ! Create a one-sided scalar solver
      create_scalar: block
         use param, only: param_read
         use tpscalar_class, only: Lphase,Gphase
         use hypre_str_class, only: pcg_pfmg2,pfmg,gmres_pfmg
         integer :: i,j,k
         real(WP) :: Ldiff,Gdiff
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=2,name='tpscalar')
         ! Initialize the phase specific VOF
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Initialize the phase specific density
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         ! Temperature on the liquid and gas sides
         sc%SCname=[  'Tl',  'Tg']; iTl=1; iTg=2
         sc%phase =[Lphase,Gphase]
         ! Read diffusivity
         call param_read('Liquid thermal diffusivity',Ldiff)
         sc%diff(:,:,:,iTl)=Ldiff
         call param_read('Gas thermal diffusivity',Gdiff)
         sc%diff(:,:,:,iTg)=Gdiff
         ! Configure implicit scalar solver
         ss=ddadi(cfg=cfg,name='Scalar',nst=7)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize scalar fields
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  ! Liquid scalar
                  if (vf%VF(i,j,k).gt.0.0_WP) then
                     ! We are in the liquid
                     if (cfg%ym(j).gt.depth+cfg%dy(j)) then
                        ! We are above the pool
                        sc%SC(i,j,k,iTl)=1.0_WP
                     else
                        ! We are in the pool
                        sc%SC(i,j,k,iTl)=2.0_WP
                     end if
                  end if
                  ! Gas scalar
                  if (vf%VF(i,j,k).lt.1.0_WP) then
                     ! We are in the gas
                     sc%SC(i,j,k,iTg)=(cfg%ym(j)-depth)/(cfg%yL-depth)
                  end if
               end do
            end do
         end do
      end block create_scalar
      

      ! Create a framework for shifting the evaporation source term
      initialize_mflux: block
         ! Initialize a pseudo time tracker
         pseudo_time=timetracker(amRoot=cfg%amRoot,name='Pseudo',print_info=.false.)
         call param_read('Max pseudo timestep size',pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',pseudo_time%cflmax)
         call param_read('Max pseudo time steps',pseudo_time%nmax)
         call param_read('Tolerence',mflux_tol)
         call param_read('Evaporation mass flux',evp_mass_flux)
         pseudo_time%dt=pseudo_time%dtmax
         ! Initialize the evaporation source term and errors
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP)); mflux=evp_mass_flux; else where; mflux=0.0_WP; end where
         mfluxL=mflux; mfluxG=mflux
         mflux_err=0.0_WP; mfluxL_err=0.0_WP; mfluxG_err=0.0_WP
         mfluxL_int_err=0.0_WP; mfluxG_int_err=0.0_WP
         ! Integral
         call cfg%integrate(mflux,mflux_int)
         call cfg%integrate(mfluxL,mfluxL_int)
         call cfg%integrate(mfluxG,mfluxG_int)
      end block initialize_mflux
      

      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: nsc
         logical :: ens_for_mflux

         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='VaporizingDrop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_scalar('curvature',vf%curv)
         call ens_out%add_surface('plic',smesh)
         do nsc=1,sc%nscalar
           call ens_out%add_scalar(trim(sc%SCname(nsc)),sc%SC(:,:,:,nsc))
         end do
         call ens_out%add_scalar('mflux' ,mflux)
         call ens_out%add_scalar('mfluxL',mfluxL)
         call ens_out%add_scalar('mfluxG',mfluxG)
         call ens_out%add_scalar('mfluxL_err',mfluxL_err_field)
         call ens_out%add_vector('VFgrad',VFgradx,VFgrady,VFgradz)
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)

         ! Creat ensight for mflux evolution
         call param_read('Ensight for mflux evolution',ens_for_mflux)
         if (ens_for_mflux) then
            call param_read('When to write mflux ensight',mflux_ens_time)
            ens_mflux_out=ensight(cfg=cfg,name='mfluxEvolution')
            ens_mflux_evt=event(time=pseudo_time,name='mflux ensight')
            call param_read('Ensight output period for mflux',ens_mflux_evt%tper)
            call ens_mflux_out%add_scalar('VOF'  ,vf%VF)
            call ens_mflux_out%add_scalar('mflux' ,mflux)
            call ens_mflux_out%add_scalar('mfluxL',mfluxL)
            call ens_mflux_out%add_scalar('mfluxG',mfluxG)
            call ens_mflux_out%add_scalar('mfluxL_err',mfluxL_err_field)
            if (ens_mflux_evt%occurs()) call ens_mflux_out%write_data(pseudo_time%t)
         else
            mflux_ens_time=-1.0_WP
         end if
      end block create_ensight
      
      
      ! Create a monitor file
      create_monitor: block
         integer :: nsc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call sc%get_max(VF=vf%VF)
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
         call mfile%add_column(vf%VFmax,'VOF maximum')
         call mfile%add_column(vf%VFmin,'VOF minimum')
         call mfile%add_column(vf%VFint,'VOF integral')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%write()
         ! Create CFL monitor
         cflfile=monitor(fs%cfg%amRoot,'cfl')
         call cflfile%add_column(time%n,'Timestep number')
         call cflfile%add_column(time%t,'Time')
         call cflfile%add_column(fs%CFLst,'STension CFL')
         call cflfile%add_column(fs%CFLc_x,'Convective xCFL')
         call cflfile%add_column(fs%CFLc_y,'Convective yCFL')
         call cflfile%add_column(fs%CFLc_z,'Convective zCFL')
         call cflfile%add_column(fs%CFLv_x,'Viscous xCFL')
         call cflfile%add_column(fs%CFLv_y,'Viscous yCFL')
         call cflfile%add_column(fs%CFLv_z,'Viscous zCFL')
         call cflfile%write()
         ! Create scalar monitor
         scfile=monitor(sc%cfg%amRoot,'scalar')
         call scfile%add_column(time%n,'Timestep number')
         call scfile%add_column(time%t,'Time')
         do nsc=1,sc%nscalar
           call scfile%add_column(sc%SCmin(nsc),trim(sc%SCname(nsc))//'_min')
           call scfile%add_column(sc%SCmax(nsc),trim(sc%SCname(nsc))//'_max')
           call scfile%add_column(sc%SCint(nsc),trim(sc%SCname(nsc))//'_int')
         end do
         call scfile%write()
         ! Create mflux monitor
         mfluxfile=monitor(sc%cfg%amRoot,'mflux')
         call mfluxfile%add_column(time%n,'Timestep number')
         call mfluxfile%add_column(time%t,'Time')
         call mfluxfile%add_column(pseudo_time%t,'Pseudo time')
         call mfluxfile%add_column(pseudo_time%dt,'Pseudo time step')
         call mfluxfile%add_column(pseudo_time%cfl,'Maximum pseudo CFL')
         call mfluxfile%add_column(pseudo_time%n,'No. pseudo steps')
         call mfluxfile%add_column(mflux_int,'mflux int')
         call mfluxfile%add_column(mfluxL_int,'shifted mfluxL int')
         call mfluxfile%add_column(mfluxG_int,'shifted mfluxG int')
         call mfluxfile%add_column(mfluxL_int_err,'mfluxL int err')
         call mfluxfile%add_column(mfluxG_int_err,'mfluxG int err')
         call mfluxfile%add_column(mfluxL_err,'max mfluxL err')
         call mfluxfile%add_column(mfluxG_err,'max mfluxG err')
         call mfluxfile%write()
      end block create_monitor
      
      
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation - this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use tpns_class, only: static_contact,harmonic_visc
      use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
      implicit none
      real(WP), dimension(:,:,:), allocatable :: U_pc,V_pc,W_pc
      real(WP) :: VFm,VFp
      real(WP), dimension(3) :: normalm,normalp,normal_tmp
      logical :: is_interfacial_m,is_interfacial_p
      integer :: i,j,k

      ! Allocate the phase change velocity components
      allocate(U_pc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(V_pc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      allocate(W_pc(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      
      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF
         
         ! Remember old SC
         sc%SCold =sc%SC
         sc%PVFold=sc%PVF
         
         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         
         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)
         
         ! Update the evaporation source term
         where ((vf%VF.gt.0.0_WP).and.(vf%VF.lt.1.0_WP)); mflux=evp_mass_flux; else where; mflux=0.0_WP; end where
         
         ! Prepare a phase-change velocity field
         ! X-direction
         U_pc=fs%U
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_
               do i=cfg%imin_,cfg%imax_+1
                  VFm=vf%VF(i-1,j,k)
                  VFp=vf%VF(i  ,j,k)
                  is_interfacial_m=VFm.gt.0.0_WP.and.VFm.lt.1.0_WP
                  is_interfacial_p=VFp.gt.0.0_WP.and.VFp.lt.1.0_WP
                  if (is_interfacial_m) then
                     normalm=calculateNormal(vf%interface_polygon(1,i-1,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i-1,j,k)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i-1,j,k))
                        normalm=0.5_WP*(normalm+normal_tmp)
                     end if
                     if (is_interfacial_p) then
                        normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                        if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                           normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                           normalp=0.5_WP*(normalp+normal_tmp)
                        end if
                        U_pc(i,j,k)=U_pc(i,j,k)+sum(sc%itp_x(:,i,j,k)*mflux(i-1:i,j,k)*[normalm(1),normalp(1)]/fs%rho_l)
                     else
                        U_pc(i,j,k)=U_pc(i,j,k)+mflux(i-1,j,k)*normalm(1)/fs%rho_l
                     end if
                  else if (is_interfacial_p) then
                     normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                        normalp=0.5_WP*(normalp+normal_tmp)
                     end if
                     U_pc(i,j,k)=U_pc(i,j,k)+mflux(i,j,k)*normalp(1)/fs%rho_l
                  end if
               end do
            end do
         end do
         call cfg%sync(U_pc)
         ! Y-direction
         V_pc=fs%V
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_+1
               do i=cfg%imin_,cfg%imax_
                  VFm=vf%VF(i,j-1,k)
                  VFp=vf%VF(i,j  ,k)
                  is_interfacial_m=VFm.gt.0.0_WP.and.VFm.lt.1.0_WP
                  is_interfacial_p=VFp.gt.0.0_WP.and.VFp.lt.1.0_WP
                  if (is_interfacial_m) then
                     normalm=calculateNormal(vf%interface_polygon(1,i,j-1,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j-1,k)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i,j-1,k))
                        normalm=0.5_WP*(normalm+normal_tmp)
                     end if
                     if (is_interfacial_p) then
                        normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                        if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                           normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                           normalp=0.5_WP*(normalp+normal_tmp)
                        end if
                        V_pc(i,j,k)=V_pc(i,j,k)+sum(sc%itp_y(:,i,j,k)*mflux(i,j-1:j,k)*[normalm(2),normalp(2)]/fs%rho_l)
                     else
                        V_pc(i,j,k)=V_pc(i,j,k)+mflux(i,j-1,k)*normalm(2)/fs%rho_l
                     end if
                  else if (is_interfacial_p) then
                     normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                        normalp=0.5_WP*(normalp+normal_tmp)
                     end if
                     V_pc(i,j,k)=V_pc(i,j,k)+mflux(i,j,k)*normalp(2)/fs%rho_l
                  end if
               end do
            end do
         end do
         call cfg%sync(V_pc)
         ! Z-direction
         W_pc=fs%W
         do k=cfg%kmin_,cfg%kmax_
            do j=cfg%jmin_,cfg%jmax_+1
               do i=cfg%imin_,cfg%imax_
                  VFm=vf%VF(i,j,k-1)
                  VFp=vf%VF(i,j,k  )
                  is_interfacial_m=VFm.gt.0.0_WP.and.VFm.lt.1.0_WP
                  is_interfacial_p=VFp.gt.0.0_WP.and.VFp.lt.1.0_WP
                  if (is_interfacial_m) then
                     normalm=calculateNormal(vf%interface_polygon(1,i,j,k-1))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k-1)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k-1))
                        normalm=0.5_WP*(normalm+normal_tmp)
                     end if
                     if (is_interfacial_p) then
                        normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                        if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                           normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                           normalp=0.5_WP*(normalp+normal_tmp)
                        end if
                        W_pc(i,j,k)=W_pc(i,j,k)+sum(sc%itp_z(:,i,j,k)*mflux(i,j,k-1:k)*[normalm(3),normalp(3)]/fs%rho_l)
                     else
                        W_pc(i,j,k)=W_pc(i,j,k)+mflux(i,j,k-1)*normalm(3)/fs%rho_l
                     end if
                  else if (is_interfacial_p) then
                     normalp=calculateNormal(vf%interface_polygon(1,i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                        normal_tmp=calculateNormal(vf%interface_polygon(2,i,j,k))
                        normalp=0.5_WP*(normalp+normal_tmp)
                     end if
                     W_pc(i,j,k)=W_pc(i,j,k)+mflux(i,j,k)*normalp(3)/fs%rho_l
                  end if
               end do
            end do
         end do
         call cfg%sync(W_pc)
         
         ! VOF solver step
         call vf%advance(dt=time%dt,U=U_pc,V=V_pc,W=W_pc)

         ! ! Now transport our phase-specific scalars
         ! advance_scalar: block
         !    use tpscalar_class, only: Lphase,Gphase
         !    integer :: nsc
            
         !    ! Get the phas-specific VOF
         !    sc%PVF(:,:,:,Lphase)=vf%VF
         !    sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
            
         !    ! Explicit calculation of dSC/dt from advective term
         !    call sc%get_dSCdt(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,VFold=vf%VFold,VF=vf%VF,detailed_face_flux=vf%detailed_face_flux,dt=time%dt)
            
         !    ! Advance advection
         !    do nsc=1,sc%nscalar
         !       where (sc%mask.eq.0.and.sc%PVF(:,:,:,sc%phase(nsc)).gt.0.0_WP) sc%SC(:,:,:,nsc)=((sc%PVFold(:,:,:,sc%phase(nsc)))*sc%SCold(:,:,:,nsc)+time%dt*resSC(:,:,:,nsc))/(sc%PVF(:,:,:,sc%phase(nsc)))
         !       where (sc%PVF(:,:,:,sc%phase(nsc)).eq.0.0_WP) sc%SC(:,:,:,nsc)=0.0_WP
         !    end do

         !    ! Apply the mass/energy transfer source term for the species/temperature
         !    do nsc=1,sc%nscalar
         !       where (sc%PVF(:,:,:,sc%phase(nsc)).gt.0.0_WP.and.sc%PVF(:,:,:,sc%phase(nsc)).lt.1.0_WP) sc%SC(:,:,:,nsc)=sc%SC(:,:,:,nsc)-mflux/sc%Prho(sc%phase(nsc))*vf%SD*time%dt
         !    end do
               
         !    ! Advance diffusion
         !    call sc%solve_implicit(time%dt,sc%SC)
            
         !    ! Apply boundary conditions
         !    call sc%apply_bcond(time%t,time%dt)
               
         ! end block advance_scalar         
               
         ! Shift the evaporation source term away from the interface
         shift_mflux: block
            use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX
            use parallel, only: MPI_REAL_WP
            use irl_fortran_interface, only: calculateNormal,getNumberOfVertices
            real(WP), dimension(:,:,:), allocatable :: ccVFgradX,ccVFgradY,ccVFgradZ
            integer  :: ierr,i,j,k
            real(WP) :: my_mflux_err
            real(WP), dimension(3) :: n1,n2
            
            ! Allocate memory for the cell-centered VOF gradient
            allocate(ccVFgradX(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            allocate(ccVFgradY(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            allocate(ccVFgradZ(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
            
            ! Get the scaled gradient of VOF
            do k=vf%cfg%kmino_,vf%cfg%kmaxo_
               do j=vf%cfg%jmino_,vf%cfg%jmaxo_
                  do i=vf%cfg%imino_,vf%cfg%imaxo_
                     n1=calculateNormal(vf%interface_polygon(1,i,j,k))
                     if (getNumberOfVertices(vf%interface_polygon(2,i,j,k)).gt.0) then
                        n2=calculateNormal(vf%interface_polygon(2,i,j,k))
                        n1=0.5_WP*(n1+n2)
                     end if
                     ccVFgradX(i,j,k)=-n1(1)
                     ccVFgradY(i,j,k)=-n1(2)
                     ccVFgradZ(i,j,k)=-n1(3)
                  end do
               end do
            end do
            call sc%cellVec_to_face(ccf_x=ccVFgradX,ccf_y=ccVFgradY,ccf_z=ccVFgradZ,fcf_x=VFgradX,fcf_y=VFgradY,fcf_z=VFgradZ)
            
            ! Deallocate the unused
            deallocate(ccVFgradX,ccVFgradY,ccVFgradZ)
            
            ! Get the CFL based on the gradient of the VOF
            call sc%get_cfl(VFgradX,VFgradY,VFgradZ,pseudo_time%dt,pseudo_time%cfl)
            
            ! Reset the pseudo time
            call pseudo_time%reset()
            
            ! Adjust the pseudo time step
            call pseudo_time%adjust_dt()
            
            ! Initialize the evaporation source terms on the liquid and gas sides
            mfluxL=mflux; mfluxG=mflux

            ! Move the evaporation source term away from the interface
            do while (.not.pseudo_time%done())
               
               ! Remember old mflux
               mfluxL_old=mfluxL
               mfluxG_old=mfluxG
               
               ! Increment pseudo time
               call pseudo_time%increment()
               
               ! Assemble explicit residual
               call sc%get_dmfluxdtau( VFgradX, VFgradY, VFgradZ,mfluxL_old,resmfluxL)
               call sc%get_dmfluxdtau(-VFgradX,-VFgradY,-VFgradZ,mfluxG_old,resmfluxG)
               
               ! Apply these residuals
               mfluxL=mfluxL_old+pseudo_time%dt*resmfluxL
               mfluxG=mfluxG_old+pseudo_time%dt*resmfluxG

               ! Output to ensight
               mfluxL_err_field=mfluxL-mfluxL_old
               if (mflux_ens_time.eq.time%t.and.ens_mflux_evt%occurs()) then
                  call ens_mflux_out%write_data(pseudo_time%t)
               end if

               ! Calculate the error on the liquid side
               my_mflux_err=maxval(abs(mfluxL_err_field))
               call MPI_ALLREDUCE(my_mflux_err,mfluxL_err,1,MPI_REAL_WP,MPI_Max,sc%cfg%comm,ierr)
               ! Calculate the error on the gas side
               my_mflux_err=maxval(abs(mfluxG-mfluxG_old))
               call MPI_ALLREDUCE(my_mflux_err,mfluxG_err,1,MPI_REAL_WP,MPI_Max,sc%cfg%comm,ierr)
               ! Check convergence
               mflux_err=max(mfluxL_err,mfluxG_err)
               if (mflux_err.lt.mflux_tol) exit

            end do

            ! Integral of mflux
            call cfg%integrate(mflux,mflux_int)
            call cfg%integrate(mfluxL,mfluxL_int)
            call cfg%integrate(mfluxG,mfluxG_int)
            mfluxL_int_err=abs(mfluxL_int-mflux_int)
            mfluxG_int_err=abs(mfluxG_int-mflux_int)

         end block shift_mflux

         ! Prepare new staggered viscosity (at n+1)
         call fs%get_viscosity(vf=vf,strat=harmonic_visc)
         
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! Build mid-time velocity
            fs%U=0.5_WP*(fs%U+fs%Uold)
            fs%V=0.5_WP*(fs%V+fs%Vold)
            fs%W=0.5_WP*(fs%W+fs%Wold)
            
            ! Preliminary mass and momentum transport step at the interface
            call fs%prepare_advection_upwind(dt=time%dt)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            
            ! Add momentum source terms
            call fs%addsrc_gravity(resU,resV,resW)
            
            ! Assemble explicit residual
            resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
            resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
            resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW
            
            ! Form implicit residuals
            !call fs%solve_implicit(time%dt,resU,resV,resW)
            
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU/fs%rho_U
            fs%V=2.0_WP*fs%V-fs%Vold+resV/fs%rho_V
            fs%W=2.0_WP*fs%W-fs%Wold+resW/fs%rho_W
            
            ! Apply other boundary conditions
            call fs%apply_bcond(time%t,time%dt)
            
            ! Solve Poisson equation
            call fs%update_laplacian()
            call fs%correct_mfr()
            call fs%get_div()
            ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
            call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
            fs%psolv%rhs=-fs%cfg%vol*(fs%div+(mfluxG/fs%rho_g-mfluxL/fs%rho_l)*vf%SD)/time%dt ! Evaporation mass source term is taken into account here
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            
            ! Correct velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            fs%P=fs%P+fs%psolv%sol
            fs%U=fs%U-time%dt*resU/fs%rho_U
            fs%V=fs%V-time%dt*resV/fs%rho_V
            fs%W=fs%W-time%dt*resW/fs%rho_W
            
            ! Increment sub-iteration counter
            time%it=time%it+1
            
         end do
         
         ! Recompute interpolated velocity and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
         
         ! Output to ensight
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if
         
         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max(VF=vf%VF)
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call mfluxfile%write()
         
      end do
      
      deallocate(U_pc,V_pc,W_pc)

   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      deallocate(resU,resV,resW,Ui,Vi,Wi,resSC,mflux,mfluxL,mfluxL_old,mfluxG,mfluxG_old,resmfluxL,resmfluxG)
   end subroutine simulation_final
   
   
end module simulation