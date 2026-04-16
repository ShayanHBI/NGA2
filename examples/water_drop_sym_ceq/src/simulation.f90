!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use hypre_str_class,   only: hypre_str
   use ddadi_class,       only: ddadi
   use tpns_class,        only: tpns
   use vfs_class,         only: vfs,VFlo,VFhi
   use tpscalar_class,    only: tpscalar,Lphase,Gphase
   use lgpc_class,        only: lgpc
   use string,            only: str_short,str_medium
   use YAMLRead,          only: YAMLElement
   use chem_sys_class,    only: chem_sys
   use chem_state_class,  only: chem_state,fixed_PH
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use surfmesh_class,    only: surfmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use mathtools,         only: Pi
   implicit none
   private

   !> Get a couple linear solvers, a two-phase flow solver, a volume fraction solver and corresponding time tracker
   type(hypre_str),   public :: ps,ss
   type(ddadi),       public :: vs
   type(tpns),        public :: fs
   type(vfs),         public :: vf
   type(tpscalar),    public :: sc
   type(lgpc),        public :: lg
   type(timetracker), public :: time,timeSC

   !> The array of the species. Eeach stored as a YAMLElement object
   type(YAMLElement), dimension(:), allocatable :: species

   !> Species names
   character(len=str_medium), dimension(:), allocatable :: sp_names

   !> Chemical system and state
   type(chem_sys)   :: sys
   type(chem_state) :: state

   !> Ensight postprocessing
   type(surfmesh) :: smesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt,TYv_evt

   !> Simulation monitor file
   type(monitor) :: mfile,cflfile,scfile,lgfile

   public :: simulation_init,simulation_run,simulation_final

   !> Private work arrays
   real(WP), dimension(:,:,:,:), allocatable :: resSC
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:),   allocatable :: T
   real(WP), dimension(:),       allocatable :: MM
   ! Debug
   real(WP), dimension(:,:,:),   allocatable :: cluster_map

   !> Problem definition
   real(WP) :: R0,T_liq,T_amb,T_g,pressure,center(3)
   integer  :: iWv,iWl,iO2,iN2,iTl,iTg
   real(WP) :: rho_l,rho_g,k_l,k_g,Cp_l,Cp_g,alpha_l,alpha_g,h_lg
   real(WP) :: wv2air_rat,N2O_rat,YO2_exit,YN2_exit,Ywv_exit
   integer  :: ns,np=2
   real(WP) :: R_drp,V_drp,T_itf
   real(WP) :: ceq_liq_frac_thld=5.0e-2_WP !< Threshold for clustering trigger (configurable)
   ! Debug
   real(WP) :: prhs_int
   real(WP) :: mfr_err


contains


   !> Function that defines a level set function
   function levelset_drop(xyz,t) result(G)
      implicit none
      real(WP), dimension(3),intent(in) :: xyz
      real(WP), intent(in) :: t
      real(WP) :: G
      ! Create the drop
      G=R0-norm2(xyz-center)
   end function levelset_drop


   !> Function that returns the index of an input species name
   function get_sp_ind(name)
      use messager, only: die
      implicit none
      character(len=*), intent(in) :: name
      integer :: isc,get_sp_ind
      do isc=1,ns
         if (trim(name).eq.trim(sp_names(isc))) then
            get_sp_ind=isc
            return
         end if
      end do
      call die('[water_drop get_sp_ind] Unknown species')
   end function get_sp_ind


   !> Function that localizes the x- boundary
   function xm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin) isIn=.true.
   end function xm_locator


   !> Function that localizes the x- boundary for scalar fields
   function xm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imin-1) isIn=.true.
   end function xm_locator_sc


   !> Function that localizes the x+ boundary
   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (i.eq.pg%imax+1) isIn=.true.
   end function xp_locator


   !> Function that localizes y- boundary
   function ym_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin) isIn=.true.
   end function ym_locator


   !> Function that localizes y- boundary for scalar fields
   function ym_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmin-1) isIn=.true.
   end function ym_locator_sc


   !> Function that localizes y+ boundary
   function yp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (j.eq.pg%jmax+1) isIn=.true.
   end function yp_locator


   !> Function that localizes z- boundary
   function zm_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin) isIn=.true.
   end function zm_locator


   !> Function that localizes z- boundary for scalar fields
   function zm_locator_sc(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmin-1) isIn=.true.
   end function zm_locator_sc


   !> Function that localizes z+ boundary
   function zp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn=.false.
      if (k.eq.pg%kmax+1) isIn=.true.
   end function zp_locator


   subroutine apply_dirichlet()
      use tpns_class, only: bcond,dirichlet
      use mathtools,  only: Pi
      type(bcond), pointer :: my_bc
      real(WP) :: Ub,Ux,Uy,Uz,vfr,myR
      integer  :: i,j,k,n,stag
      call cfg%integrate(lg%div_vel,vfr)
      my_bc=>fs%first_bc
      do while (associated(my_bc))
         if (my_bc%type.ne.dirichlet) then
            my_bc=>my_bc%next
            cycle
         end if
         if (my_bc%itr%amIn) then
            select case (my_bc%face)
             case ('x')
               stag=min(my_bc%dir,0)
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
                  Ub=vfr/(4.0_WP*Pi*myR**2)
                  Ux=cfg%xm(i)/myR*Ub
                  Uy=cfg%ym(j)/myR*Ub
                  Uz=cfg%zm(j)/myR*Ub
                  fs%U(i     ,j    ,k    )=Ux
                  fs%V(i+stag,j:j+1,k    )=Uy
                  fs%W(i+stag,j    ,k:k+1)=Uz
               end do
             case ('y')
               stag=min(my_bc%dir,0)
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
                  Ub=vfr/(4.0_WP*Pi*myR**2)
                  Ux=cfg%xm(i)/myR*Ub
                  Uy=cfg%ym(j)/myR*Ub
                  Uz=cfg%zm(j)/myR*Ub
                  fs%U(i:i+1,j+stag,k    )=Ux
                  fs%V(i    ,j     ,k    )=Uy
                  fs%W(i    ,j+stag,k:k+1)=Uz
               end do
             case ('z')
               stag=min(my_bc%dir,0)
               do n=1,my_bc%itr%n_
                  i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
                  myR=sqrt(cfg%xm(i)**2+cfg%ym(j)**2+cfg%zm(k)**2)
                  Ub=vfr/(4.0_WP*Pi*myR**2)
                  Ux=cfg%xm(i)/myR*Ub
                  Uy=cfg%ym(j)/myR*Ub
                  Uz=cfg%zm(j)/myR*Ub
                  fs%U(i:i+1,j    ,k+stag)=Ux
                  fs%V(i    ,j:j+1,k+stag)=Uy
                  fs%W(i    ,j    ,k     )=Uz
               end do
            end select
         end if
         my_bc=>my_bc%next
      end do
   end subroutine apply_dirichlet


   subroutine sym_irl()
      use irl_fortran_interface, only: getPlane,new,construct_2pt,RectCub_type,&
      &                                setNumberOfPlanes,setPlane,matchVolumeFraction
      real(WP), dimension(1:4) :: plane
      type(RectCub_type) :: cell
      integer :: i,j,k
      integer :: ii,jj,kk
      call new(cell)
      if (vf%cfg%iproc.eq.1) then
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imin_-1
                  ii=(vf%cfg%imin_-i)+vf%cfg%imin_-1
                  plane=getPlane(vf%liquid_gas_interface(ii,j,k),0)
                  call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                  &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
                  plane(4)=dot_product([-plane(1),plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                  call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                  call setPlane(vf%liquid_gas_interface(i,j,k),0,[-plane(1),plane(2),plane(3)],plane(4))
                  call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
               end do
            end do
         end do
      end if
      if (vf%cfg%jproc.eq.1) then
         do k=vf%cfg%kmino_,vf%cfg%kmaxo_
            do j=vf%cfg%jmino_,vf%cfg%jmin_-1
               jj=(vf%cfg%jmin_-j)+vf%cfg%jmin_-1
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  plane=getPlane(vf%liquid_gas_interface(i,jj,k),0)
                  call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                  &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
                  plane(4)=dot_product([plane(1),-plane(2),plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                  call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                  call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),-plane(2),plane(3)],plane(4))
                  call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
               end do
            end do
         end do
      end if
      if (vf%cfg%kproc.eq.1) then
         do k=vf%cfg%kmino_,vf%cfg%kmin_-1
            kk=(vf%cfg%kmin_-k)+vf%cfg%kmin_-1
            do j=vf%cfg%jmino_,vf%cfg%jmaxo_
               do i=vf%cfg%imino_,vf%cfg%imaxo_
                  plane=getPlane(vf%liquid_gas_interface(i,j,kk),0)
                  call construct_2pt(cell,[vf%cfg%x(i  ),vf%cfg%y(j  ),vf%cfg%z(k  )],&
                  &                       [vf%cfg%x(i+1),vf%cfg%y(j+1),vf%cfg%z(k+1)])
                  plane(4)=dot_product([plane(1),plane(2),-plane(3)],[vf%cfg%xm(i),vf%cfg%ym(j),vf%cfg%zm(k)])
                  call setNumberOfPlanes(vf%liquid_gas_interface(i,j,k),1)
                  call setPlane(vf%liquid_gas_interface(i,j,k),0,[plane(1),plane(2),-plane(3)],plane(4))
                  call matchVolumeFraction(cell,vf%VF(i,j,k),vf%liquid_gas_interface(i,j,k))
               end do
            end do
         end do
      end if
   end subroutine sym_irl


   !> Calculate the spatially averaged temperature and vapor mass fraction versus radius
   subroutine get_T_Yv()
      use mathtools, only: Pi
      use string,    only: str_medium
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM,MPI_INT
      use parallel,  only: MPI_REAL_WP
      type(monitor) :: TYvfile
      real(WP) :: theta_s,theta_e,dtheta,theta
      real(WP) :: phi_s,phi_e,dphi,phi
      real(WP) :: sinTheta,cosTheta,sinPhi,cosPhi
      real(WP) :: l,dR,R_s,R_e,x,y,z
      real(WP), target :: r_out,T_out,Yv_out
      real(WP), dimension(:), allocatable :: r_,T_,T_avg,Yv_,Yv_avg
      integer  :: ntheta,itheta,nphi,iphi,nr,r_ind,nunsumble,ind(3),i,j,k,ierr
      character(len=str_medium) :: tstr
      integer :: found_local,found_global
      ! Create the monitors
      write(tstr,'(ES12.5)') time%t
      TYvfile=monitor(fs%cfg%amRoot,'T_Yv_'//trim(adjustl(tstr)))
      call TYvfile%add_column(r_out ,'r')
      call TYvfile%add_column(T_out ,'T')
      call TYvfile%add_column(Yv_out,'Yv')
      ! Grid for theta
      theta_s=0.0_WP
      theta_e=0.5_WP*Pi
      ntheta=5
      dtheta=(theta_e-theta_s)/real(ntheta-1,WP)
      ! Grid for phi
      phi_s=0.0_WP
      phi_e=0.5_WP*Pi
      nphi=5
      dphi=(phi_e-phi_s)/real(nphi-1,WP)
      ! Grid for radius
      R_s=0.0_WP
      R_e=590e-6
      dR=cfg%dx(1)
      nr=int((R_e-R_s)/dR)+1
      allocate(r_(nr))
      allocate(T_(nr)); T_=0.0_WP
      allocate(T_avg(nr)); T_avg=0.0_WP
      allocate(Yv_(nr)); Yv_=0.0_WP
      allocate(Yv_avg(nr)); Yv_avg=0.0_WP
      r_(1)=R_s
      do r_ind=2,nr
         r_(r_ind)=r_(r_ind-1)+dR
      end do
      ! Loop over theta
      nunsumble=0
      do itheta=1,ntheta
         ! Advance theta
         theta=theta_s+real(itheta-1,WP)*dtheta
         cosTheta=cos(theta)
         sinTheta=sin(theta)
         ! Loop over phi
         do iphi=1,nphi
            nunsumble=nunsumble+1
            ! Advance phi
            phi=phi_s+real(iphi-1,WP)*dphi
            sinPhi=sin(phi)
            cosPhi=cos(phi)
            ! Initilize the line marcher
            i=cfg%imin_
            j=cfg%jmin_
            k=cfg%kmin_
            ! March the line
            do r_ind=1,nr
               ! Get the point coordinates
               x=r_(r_ind)*sinTheta
               y=r_(r_ind)*cosTheta*cosPhi
               z=r_(r_ind)*cosTheta*sinPhi
               if (cfg%is_in_subdomain([x,y,z])) then
                  ! Get the corresponding indices of the cell containing the line segment end point
                  ind=cfg%get_ijk_local(pos=[x,y,z],ind_guess=[i,j,k])
                  i=ind(1)
                  j=ind(2)
                  k=ind(3)
                  ! Get temperature
                  T_(r_ind)=T_(r_ind)+cfg%get_scalar([x,y,z],i,j,k,T,'N')
                  ! Get vapor mass fraction
                  Yv_(r_ind)=Yv_(r_ind)+cfg%get_scalar([x,y,z],i,j,k,sc%SC(:,:,:,iWv),'N')
               end if
            end do
         end do
      end do
      ! Calculate the mean quantities
      call MPI_ALLREDUCE(T_,T_avg,nr,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      T_avg=T_avg/real(nunsumble,WP)
      call MPI_ALLREDUCE(Yv_,Yv_avg,nr,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      Yv_avg=Yv_avg/real(nunsumble,WP)
      do r_ind=1,nr
         r_out =r_    (r_ind)
         T_out =T_avg (r_ind); if (T_out .lt.1e-20) T_out =0.0_WP
         Yv_out=Yv_avg(r_ind); if (Yv_out.lt.1e-20) Yv_out=0.0_WP
         call TYvfile%write()
      end do
      ! Deallocate the arrays
      deallocate(r_,T_,T_avg,Yv_,Yv_avg)
      ! Finalize the monitor object
      call TYvfile%finalize()
   end subroutine get_T_Yv


   !> Calculate the spatially averaged interface temperature
   subroutine get_T_itf()
      use mpi_f08,   only: MPI_ALLREDUCE,MPI_SUM
      use parallel,  only: MPI_REAL_WP
      integer  :: index,i,j,k,ierr
      real(WP) :: A_I,area_itf_,area_itf,T_itf_
      ! Initialize
      area_itf_=0.0_WP
      T_itf_=0.0_WP
      ! Loop over the interfacial cells
      do index=1,vf%band_count(0)
         ! Get the interfacial cell indices
         i=vf%band_map(1,index)
         j=vf%band_map(2,index)
         k=vf%band_map(3,index)
         ! Get the interface area
         A_I=cfg%vol(i,j,k)*vf%SD(i,j,k)
         ! A_I=1.0_WP
         area_itf_=area_itf_+A_I
         ! Get the interface temperature
         T_itf_=T_itf_+A_I*T(i,j,k)
      end do
      ! Calculate the averaged interface temperature
      call MPI_ALLREDUCE(area_itf_,area_itf,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      call MPI_ALLREDUCE(T_itf_,T_itf,1,MPI_REAL_WP,MPI_SUM,cfg%comm,ierr)
      T_itf=T_itf/area_itf
   end subroutine get_T_itf


   !> Calculate the droplet radius
   subroutine get_R_drp()
      use mathtools, only: Pi
      call cfg%integrate(vf%VF,V_drp)
      V_drp=8.0_WP*V_drp
      R_drp=(0.75_WP*V_drp/Pi)**(1.0_WP/3.0_WP)
   end subroutine get_R_drp


   !
   subroutine interface_jump()
      use messager, only: die
      implicit none
      real(WP), dimension(:),     allocatable :: vol_new,vol_old,mp,N,phasicHoR,Y
      logical,  dimension(:,:,:), allocatable :: processed
      integer,  dimension(:,:,:), allocatable :: needs_clustering_flag
      logical,  dimension(:),     allocatable :: active
      integer,  dimension(:,:),   allocatable :: cell_indices
      real(WP), dimension(:),     allocatable :: Vscaled,vof_old,vof_new,w
      real(WP) :: Vnew,Vold,Nsum,vof,itf_area,Tl,Tg,Tln,Tgn
      integer  :: i,j,k,index,isc,p,n_clustered,m,cluster_id,nf
      integer  :: in,jn,kn
      integer  :: stx,sty,stz
      real(WP) :: mdot2p
      real(WP), parameter :: wmin=1.0e-16_WP,dVlmin=1.0e-16_WP
      integer, parameter :: nc_max=7
      integer :: nc_cap,bfs_head
      integer, dimension(:,:), allocatable :: cell_indices_tmp
      real(WP) :: dVl,dVl_i,dVl_rem,Vref,vof_tmp,interfaceness,wsum
      real(WP) :: cluster_liq_frac
      logical  :: cluster_done
      ! Gradient-preserving: per-cell original temperatures and species
      real(WP), dimension(:),   allocatable :: Tl_orig,Tg_orig
      real(WP), dimension(:,:), allocatable :: Y_orig
      real(WP) :: Tl_old_mean,Tg_old_mean,delta_Tl,delta_Tg
      real(WP), dimension(:), allocatable :: Y_eq,delta_Y
      real(WP) :: Y_sum

      ! Debug
      cluster_map=0.0_WP
      cluster_id=0

      Vref=minval(cfg%vol)

      ! Allocate arrays
      allocate(vol_new(Lphase:Gphase))
      allocate(vol_old(Lphase:Gphase))
      allocate(mp(Lphase:Gphase))
      allocate(N(ns))
      allocate(phasicHoR(Lphase:Gphase))
      allocate(Y(ns))
      allocate(processed(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); processed=.false.
      allocate(needs_clustering_flag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); needs_clustering_flag=0

      ! Initial capacity for BFS queue / cluster cell list (will grow as needed)
      nc_cap=64
      allocate(cell_indices(3,nc_cap))

      ! Clustering stencil
      if (cfg%nx.gt.1) then
         stx=1
      else
         stx=0
      end if
      if (cfg%ny.gt.1) then
         sty=1
      else
         sty=0
      end if
      if (cfg%nz.gt.1) then
         stz=1
      else
         stz=0
      end if

      ! ========================================================================
      ! Pass 1: Identify all cells that need clustering
      ! Try single-cell CEQ on each interfacial cell; process successful
      ! single-cell cases immediately, mark failures
      ! ========================================================================
      do index=1,vf%band_count(0)

         ! Get the interfacial cell indices
         i=vf%band_map(1,index)
         j=vf%band_map(2,index)
         k=vf%band_map(3,index)

         ! Initialize single-cell quantities
         n_clustered=1
         cell_indices(:,1)=[i,j,k]
         itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
         vol_old=sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         Y =sc%SC(i,j,k,1:ns)
         Tl=sc%SC(i,j,k,iTl)
         Tg=sc%SC(i,j,k,iTg)

         ! Pre-evaluate the equilibrium
         call get_equilibrium()

         ! Check if clustering is needed
         if ((.not.state%success).or.(state%N(iWl)/(state%N(iWl)+state%N(iWv)).lt.ceq_liq_frac_thld)) then
            needs_clustering_flag(i,j,k)=1
         else
            call apply_single_cell_equilibrium_result()
            processed(i,j,k)=.true.
         end if

      end do

      ! Sync needs_clustering_flag across MPI boundaries so ghost cells are aware of neighbors' failures
      call cfg%sync(needs_clustering_flag)

      ! ========================================================================
      ! Pass 2: Process all interfacial cells
      ! Skip cells already processed in pass 1
      ! Cells needing clustering: BFS flood-fill to form clusters, then process
      ! ========================================================================
      do index=1,vf%band_count(0)

         ! Get the interfacial cell indices
         i=vf%band_map(1,index)
         j=vf%band_map(2,index)
         k=vf%band_map(3,index)

         ! Skip if already processed
         if (processed(i,j,k)) cycle

         ! Initialize cluster with the seed cell
         n_clustered=1
         cell_indices(:,1)=[i,j,k]

         ! Initialize quantities
         itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
         vol_old=sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
         Y =sc%SC(i,j,k,1:ns)
         Tl=sc%SC(i,j,k,iTl)
         Tg=sc%SC(i,j,k,iTg)

         if (needs_clustering_flag(i,j,k).eq.1) then

            ! Increment cluster ID and mark seed cell
            cluster_id=cluster_id+1
            processed(i,j,k)=.true.
            cluster_map(i,j,k)=real(cluster_id,WP)

            ! Initialize the species mass
            do isc=1,ns
               p=sc%phase(isc)
               Y(isc)=sc%Prho(p)*sc%PVF(i,j,k,p)*cfg%vol(i,j,k)*sc%SC(i,j,k,isc)
            end do

            ! Initialize the mass-averaged temperatures (mass*T for accumulation)
            Tl=sc%Prho(Lphase)*sc%PVF(i,j,k,Lphase)*cfg%vol(i,j,k)*Tl
            Tg=sc%Prho(Gphase)*sc%PVF(i,j,k,Gphase)*cfg%vol(i,j,k)*Tg

            ! BFS flood-fill: explore connected interfacial neighbors
            ! Stop when cluster liquid fraction is sufficient or max size is reached
            bfs_head=1
            cluster_done=.false.
            do while (bfs_head.le.n_clustered.and..not.cluster_done)

               ! Pop the next cell from the BFS queue
               i=cell_indices(1,bfs_head)
               j=cell_indices(2,bfs_head)
               k=cell_indices(3,bfs_head)
               bfs_head=bfs_head+1

               ! Explore face-connected neighbors (6-connectivity, no diagonal bias)
               do nf=1,6
                  in=i; jn=j; kn=k
                  select case (nf)
                   case (1); in=i-stx
                   case (2); in=i+stx
                   case (3); jn=j-sty
                   case (4); jn=j+sty
                   case (5); kn=k-stz
                   case (6); kn=k+stz
                  end select

                  ! Skip if offset is zero (2D case)
                  if (in.eq.i.and.jn.eq.j.and.kn.eq.k) cycle

                  ! Check stopping criteria before adding more cells
                  if (n_clustered.ge.nc_max) then
                     cluster_done=.true.
                     exit
                  end if

                  ! Neighbor must be interfacial, not clustered yet, and also needs clustering
                  if (vf%VF(in,jn,kn).gt.VFlo.and.vf%VF(in,jn,kn).lt.VFhi.and..not.processed(in,jn,kn).and.needs_clustering_flag(in,jn,kn).eq.1) then

                     ! Mark it as clustered
                     n_clustered=n_clustered+1

                     ! Grow the array if needed
                     if (n_clustered.gt.nc_cap) then
                        allocate(cell_indices_tmp(3,nc_cap*2))
                        cell_indices_tmp(:,1:nc_cap)=cell_indices
                        nc_cap=nc_cap*2
                        call move_alloc(cell_indices_tmp,cell_indices)
                     end if

                     cell_indices(:,n_clustered)=[in,jn,kn]
                     processed(in,jn,kn)=.true.
                     cluster_map(in,jn,kn)=real(cluster_id,WP)

                     ! Accumulate old volumes
                     vol_old=vol_old+sc%PVF(in,jn,kn,:)*cfg%vol(in,jn,kn)

                     ! Accumulate mass and mass*temperature
                     do isc=1,ns
                        p=sc%phase(isc)
                        Y(isc)=Y(isc)+sc%Prho(p)*sc%PVF(in,jn,kn,p)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,isc)
                     end do
                     Tln=sc%SC(in,jn,kn,iTl)
                     Tgn=sc%SC(in,jn,kn,iTg)
                     Tl=Tl+sc%Prho(Lphase)*sc%PVF(in,jn,kn,Lphase)*cfg%vol(in,jn,kn)*Tln
                     Tg=Tg+sc%Prho(Gphase)*sc%PVF(in,jn,kn,Gphase)*cfg%vol(in,jn,kn)*Tgn

                     ! Accumulate interface area
                     itf_area=itf_area+cfg%vol(in,jn,kn)*vf%SD(in,jn,kn)

                     ! Check if cluster now has sufficient liquid fraction
                     cluster_liq_frac=vol_old(Lphase)/sum(vol_old)
                     if (cluster_liq_frac.ge.ceq_liq_frac_thld) then
                        cluster_done=.true.
                        exit
                     end if
                  end if

               end do

            end do

            ! If only the seed cell remains, this is not a real cluster
            if (n_clustered.eq.1) then
               cluster_map(cell_indices(1,1),cell_indices(2,1),cell_indices(3,1))=0.0_WP
               cluster_id=cluster_id-1
            end if

            ! Cluster-level phase masses
            mp=sc%Prho*vol_old

            ! Cluster-averaged mass fractions and temperatures
            do isc=1,ns
               Y(isc)=Y(isc)/mp(sc%phase(isc))
            end do
            Tl=Tl/mp(Lphase)
            Tg=Tg/mp(Gphase)

         end if

         ! Get the equilibrium state
         call get_equilibrium()

         if (.not.state%success) cycle

         call apply_equilibrium_result()

      end do

      ! Sync VOF
      call cfg%sync(vf%VF)

      ! Remove flotsams and thin structures if needed
      call vf%remove_flotsams()
      call vf%remove_thinstruct()

      ! Synchronize and clean-up barycenter fields
      call vf%sync_and_clean_barycenters()

      ! Update the interface band
      call vf%update_band()

      ! Perform interface reconstruction from transported moments
      call vf%build_interface()

      ! Create discontinuous polygon mesh from IRL interface
      call vf%polygonalize_interface()

      ! Calculate curvature
      call vf%get_curvature()

      ! Reset moments to guarantee compatibility with interface reconstruction
      call vf%reset_moments()

      ! Sync fields
      do isc=1,sc%nscalar
         call cfg%sync(sc%SC(:,:,:,isc))
      end do
      call cfg%sync(sc%PVF(:,:,:,Lphase))
      call cfg%sync(sc%PVF(:,:,:,Gphase))
      call cfg%sync(lg%mdot2p)

      ! Apply boundary conditions
      call sc%apply_bcond(time%t,time%dt)
      call vf%apply_bcond(time%t,time%dt)

      ! Sync cluster map (contains per-cluster IDs from the loop above)
      call cfg%sync(cluster_map)

      ! Deallocate arrays
      deallocate(vol_new,vol_old,mp,N,phasicHoR,Y,processed,needs_clustering_flag,cell_indices)


   contains

      subroutine apply_single_cell_equilibrium_result()
         implicit none
         real(WP) :: vf_new

         ! Store the old total volume
         Vold=sum(vol_old)

         ! Calculate the old cluster VOF
         vof=vol_old(Lphase)/Vold

         ! Update the phase masses
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N(isc)*MM(isc)
         end do

         ! Get the phase volumes
         vol_new=mp/sc%Prho
         Vnew=sum(vol_new)

         ! Get the phase change mass flux
         mdot2p=(Vnew-Vold)/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*itf_area)

         ! Direct single-cell VOF update
         vf_new=vol_new(Lphase)/cfg%vol(i,j,k)

         ! Assign VOF
         if (vf_new.lt.VFlo) then
            vf%VF(i,j,k)=0.0_WP
         else if (vf_new.gt.VFhi) then
            vf%VF(i,j,k)=1.0_WP
         else
            vf%VF(i,j,k)=vf_new
         end if
         sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
         sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)

         ! Assign equilibrium species mass fractions directly
         do isc=1,ns
            p=sc%phase(isc)
            if (sc%PVF(i,j,k,p).gt.0.0_WP.and.mp(p).gt.0.0_WP) then
               sc%SC(i,j,k,isc)=MM(isc)*N(isc)/mp(p)
            else
               sc%SC(i,j,k,isc)=0.0_WP
            end if
         end do

         ! Renormalize species per phase to ensure mass conservation
         do p=Lphase,Gphase
            Y_sum=0.0_WP
            do isc=1,ns
               if (sc%phase(isc).eq.p) Y_sum=Y_sum+sc%SC(i,j,k,isc)
            end do
            if (Y_sum.gt.0.0_WP) then
               do isc=1,ns
                  if (sc%phase(isc).eq.p) sc%SC(i,j,k,isc)=sc%SC(i,j,k,isc)/Y_sum
               end do
            end if
         end do

         ! Assign equilibrium temperature directly
         if (vf%VF(i,j,k).eq.1.0_WP) then
            sc%SC(i,j,k,iTl)=state%T
            sc%SC(i,j,k,iTg)=0.0_WP
            lg%mdot2p(i,j,k)=0.0_WP
         else if (vf%VF(i,j,k).eq.0.0_WP) then
            sc%SC(i,j,k,iTl)=0.0_WP
            sc%SC(i,j,k,iTg)=state%T
            lg%mdot2p(i,j,k)=0.0_WP
         else
            sc%SC(i,j,k,iTl)=state%T
            sc%SC(i,j,k,iTg)=state%T
            lg%mdot2p(i,j,k)=mdot2p
         end if

      end subroutine apply_single_cell_equilibrium_result

      subroutine apply_equilibrium_result()
         implicit none

         ! Store the old total volume
         Vold=sum(vol_old)

         ! Calculate the cluster VOF
         vof=vol_old(Lphase)/Vold

         ! Update the phase masses
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N(isc)*MM(isc)
         end do

         ! Get the phase volumes
         vol_new=mp/sc%Prho
         Vnew=sum(vol_new)

         ! Get the phase change mass flux
         mdot2p=(Vnew-Vold)/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*itf_area)

         ! Allocate per-cluster work arrays
         allocate(Vscaled(n_clustered)); Vscaled=0.0_WP
         allocate(vof_old(n_clustered)); vof_old=0.0_WP
         allocate(vof_new(n_clustered)); vof_new=0.0_WP
         allocate(w(n_clustered)); w=0.0_WP
         allocate(active(n_clustered)); active=.false.

         ! --- Gradient-preserving: store original per-cell T and Y ---
         allocate(Tl_orig(n_clustered))
         allocate(Tg_orig(n_clustered))
         allocate(Y_orig(ns,n_clustered))
         allocate(Y_eq(ns))
         allocate(delta_Y(ns))
         do m=1,n_clustered
            i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
            Tl_orig(m)=sc%SC(i,j,k,iTl)
            Tg_orig(m)=sc%SC(i,j,k,iTg)
            do isc=1,ns
               Y_orig(isc,m)=sc%SC(i,j,k,isc)
            end do
         end do

         ! Compute old cluster-mean temperature (mass-weighted) and species
         Tl_old_mean=Tl
         Tg_old_mean=Tg

         ! Compute deltas: new equilibrium vs old cluster mean
         delta_Tl=state%T-Tl_old_mean
         delta_Tg=state%T-Tg_old_mean

         ! Compute equilibrium species mass fractions
         do isc=1,ns
            p=sc%phase(isc)
            if (mp(p).gt.0.0_WP) then
               Y_eq(isc)=MM(isc)*N(isc)/mp(p)
            else
               Y_eq(isc)=0.0_WP
            end if
         end do

         ! Delta for species mass fractions
         delta_Y=Y_eq-Y

         ! Gather geometry and current VOF per clustered cell
         do m=1,n_clustered
            i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
            active(m)=.true.
            Vscaled(m)=cfg%vol(i,j,k)/Vref ! Scale it for more accurate calculations
            vof_old(m)=vf%VF(i,j,k)
            vof_new(m)=vof_old(m)
         end do

         ! Total liquid volume change ( > 0 condensation, < 0 vaporization)
         dVl=(vol_new(Lphase)-vol_old(Lphase))/Vref ! Scale it for more accurate calculations
         dVl_rem=dVl

         if (abs(dVl).gt.dVlmin) then

            ! Build weights
            do m=1,n_clustered
               i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
               itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
               interfaceness=minval(sc%PVF(i,j,k,:))
               w(m)=max(itf_area*interfaceness,wmin)
            end do

            ! Iteratively redistribute liquid
            do

               ! Update weights sum
               wsum=0.0_WP
               do m=1,n_clustered
                  if (active(m)) wsum=wsum+w(m)
               end do

               ! Terminate if succeeded
               if (abs(dVl_rem).le.dVlmin.or.wsum.le.0.0_WP) exit

               ! Distribute
               dVl=dVl_rem
               do m=1,n_clustered
                  if (.not.active(m)) cycle

                  ! Estimate VOF
                  dVl_i =(w(m)/wsum)*dVl
                  vof_tmp=vof_new(m)+dVl_i/Vscaled(m)

                  ! Clip it
                  if (vof_tmp.gt.1.0_WP) then
                     dVl_i=(1.0_WP-vof_new(m))*Vscaled(m)
                     vof_new(m)=1.0_WP
                     active(m)=.false.
                  else if (vof_tmp.lt.0.0_WP) then
                     dVl_i=(0.0_WP-vof_new(m))*Vscaled(m)
                     vof_new(m)=0.0_WP
                     active(m)=.false.
                  else
                     vof_new(m)=vof_tmp
                  end if

                  ! Correct the liquid volume change
                  dVl_rem=dVl_rem-dVl_i

               end do

            end do

         end if

         ! Assign per-cell fields (Need to treat cells with VOF=0 and 1, differently)
         do m=1,n_clustered

            ! Get the cell indices
            i=cell_indices(1,m)
            j=cell_indices(2,m)
            k=cell_indices(3,m)

            ! Skip ghost cells (only write to owned cells; ghosts will be overwritten by sync)
            if (i.lt.cfg%imin_.or.i.gt.cfg%imax_) cycle
            if (j.lt.cfg%jmin_.or.j.gt.cfg%jmax_) cycle
            if (k.lt.cfg%kmin_.or.k.gt.cfg%kmax_) cycle

            ! Assign VOF
            if (vof_new(m).lt.VFlo) then
               vf%VF(i,j,k)=0.0_WP
            else if (vof_new(m).gt.VFhi) then
               vf%VF(i,j,k)=1.0_WP
            else
               vf%VF(i,j,k)=vof_new(m)
            end if
            sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
            sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)

            ! Composition: gradient-preserving with additive delta
            do isc=1,ns
               p=sc%phase(isc)
               if(sc%PVF(i,j,k,p).gt.0.0_WP) then
                  sc%SC(i,j,k,isc)=Y_orig(isc,m)+delta_Y(isc)
               else
                  sc%SC(i,j,k,isc)=0.0_WP
               end if
            end do
            ! Renormalize species per phase to ensure mass conservation
            do p=Lphase,Gphase
               Y_sum=0.0_WP
               do isc=1,ns
                  if (sc%phase(isc).eq.p) Y_sum=Y_sum+sc%SC(i,j,k,isc)
               end do
               if (Y_sum.gt.0.0_WP) then
                  do isc=1,ns
                     if (sc%phase(isc).eq.p) sc%SC(i,j,k,isc)=sc%SC(i,j,k,isc)/Y_sum
                  end do
               end if
            end do

            ! Temperature and phase change mass flux
            if (vf%VF(i,j,k).eq.1.0_WP) then
               sc%SC(i,j,k,iTl)=Tl_orig(m)+delta_Tl
               sc%SC(i,j,k,iTg)=0.0_WP
               lg%mdot2p(i,j,k)=0.0_WP
            else if (vf%VF(i,j,k).eq.0.0_WP) then
               sc%SC(i,j,k,iTl)=0.0_WP
               sc%SC(i,j,k,iTg)=Tg_orig(m)+delta_Tg
               lg%mdot2p(i,j,k)=0.0_WP
            else
               sc%SC(i,j,k,iTl)=Tl_orig(m)+delta_Tl
               sc%SC(i,j,k,iTg)=Tg_orig(m)+delta_Tg
               lg%mdot2p(i,j,k)=mdot2p
            end if

         end do

         ! Free per-cluster work arrays
         deallocate(Vscaled,vof_old,vof_new,w,active)
         deallocate(Tl_orig,Tg_orig,Y_orig,Y_eq,delta_Y)

      end subroutine apply_equilibrium_result

      subroutine get_equilibrium()
         implicit none

         ! Calculate and normalize the mole numbers
         do isc=1,ns
            N(isc)=Y(isc)*mp(sc%phase(isc))/MM(isc)
         end do
         Nsum=sum(N)
         if (Nsum.gt.0.0_WP) N=N/Nsum

         ! Get the phasic enthalpies
         call state%get_phasic_HoR(Lphase,N,Tl,phasicHoR(Lphase))
         call state%get_phasic_HoR(Gphase,N,Tg,phasicHoR(Gphase))

         ! Reinitialize the mole numbers
         call state%N_init(N=N,HoR=sum(phasicHoR),T_g=T_g)
         if (.not.state%success) then
            print*,'Cluster VOF = ',vof
            print*,'N = ',N
            print*,'N*Nsum = ',N*Nsum
            print*,'HoR = ',sum(phasicHoR)
            print*,'T_g = ',T_g
            print*,'Clustered cells info:'
            print*,'n_clustered = ',n_clustered
            do m=1,n_clustered
               i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
               print*,'i,j,k = ',i,j,k
               print*,'VOF = ',vf%VF(i,j,k)
               print*,'SD = ',vf%SD(i,j,k)
            end do
            call die('interface_jump: N_init failed')
         end if

         ! Get the chemical equilibrium
         call state%equilibrate()

         ! Re-scale the mole numbers
         N=state%N*Nsum

      end subroutine get_equilibrium

   end subroutine interface_jump 


   ! Claude version
   ! subroutine interface_jump()
   !    use messager, only: die
   !    implicit none
   !    real(WP), dimension(:),     allocatable :: vol_new,vol_old,mp,N,phasicHoR,Y
   !    logical,  dimension(:,:,:), allocatable :: clustered
   !    integer,  dimension(:,:,:), allocatable :: needs_clustering_flag
   !    logical,  dimension(:),     allocatable :: active
   !    integer,  dimension(:,:),   allocatable :: cell_indices
   !    real(WP), dimension(:),     allocatable :: Vscaled,vof_old,vof_new,w
   !    real(WP) :: Vnew,Vold,Nsum,vof,itf_area,Tl,Tg,Tln,Tgn
   !    integer  :: i,j,k,index,isc,p,n_clustered,m,cluster_id
   !    integer  :: in,jn,kn
   !    integer  :: stx,sty,stz
   !    real(WP) :: mdot2p
   !    real(WP), parameter :: wmin=1.0e-16_WP,dVlmin=1.0e-16_WP
   !    integer, parameter :: nc_max=27
   !    integer :: nc_cap,bfs_head
   !    integer, dimension(:,:), allocatable :: cell_indices_tmp
   !    real(WP) :: dVl,dVl_i,dVl_rem,Vref,vof_tmp,interfaceness,wsum
   !    real(WP) :: cluster_liq_frac
   !    logical  :: cluster_done
   !    ! Gradient-preserving: per-cell original temperatures and species
   !    real(WP), dimension(:),   allocatable :: Tl_orig,Tg_orig
   !    real(WP), dimension(:,:), allocatable :: Y_orig
   !    real(WP) :: Tl_old_mean,Tg_old_mean,delta_Tl,delta_Tg
   !    real(WP), dimension(:), allocatable :: Y_eq,delta_Y
   !    real(WP) :: Y_sum

   !    ! Debug
   !    dbg_flg=0.0_WP
   !    cluster_id=0

   !    Vref=minval(cfg%vol)

   !    ! Allocate arrays
   !    allocate(vol_new(Lphase:Gphase))
   !    allocate(vol_old(Lphase:Gphase))
   !    allocate(mp(Lphase:Gphase))
   !    allocate(N(ns))
   !    allocate(phasicHoR(Lphase:Gphase))
   !    allocate(Y(ns))
   !    allocate(clustered(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); clustered=.false.
   !    allocate(needs_clustering_flag(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); needs_clustering_flag=0

   !    ! Initial capacity for BFS queue / cluster cell list (will grow as needed)
   !    nc_cap=64
   !    allocate(cell_indices(3,nc_cap))

   !    ! Clustering stencil
   !    if (cfg%nx.gt.1) then
   !       stx=1
   !    else
   !       stx=0
   !    end if
   !    if (cfg%ny.gt.1) then
   !       sty=1
   !    else
   !       sty=0
   !    end if
   !    if (cfg%nz.gt.1) then
   !       stz=1
   !    else
   !       stz=0
   !    end if

   !    ! ========================================================================
   !    ! Pass 1: Identify all cells that need clustering
   !    ! Try single-cell CEQ on each interfacial cell; mark failures
   !    ! ========================================================================
   !    do index=1,vf%band_count(0)

   !       ! Get the interfacial cell indices
   !       i=vf%band_map(1,index)
   !       j=vf%band_map(2,index)
   !       k=vf%band_map(3,index)

   !       ! Pre-evaluate the equilibrium
   !       mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
   !       Y =sc%SC(i,j,k,1:ns)
   !       Tl=sc%SC(i,j,k,iTl)
   !       Tg=sc%SC(i,j,k,iTg)
   !       call get_equilibrium()

   !       ! Check if clustering is needed
   !       if ((.not.state%success).or.(state%N(iWl)/(state%N(iWl)+state%N(iWv)).lt.ceq_liq_frac_thld)) then
   !          needs_clustering_flag(i,j,k)=1
   !       end if

   !    end do

   !    ! Sync needs_clustering_flag across MPI boundaries so ghost cells are aware of neighbors' failures
   !    call cfg%sync(needs_clustering_flag)

   !    ! ========================================================================
   !    ! Pass 2: Process all interfacial cells
   !    ! Cells not needing clustering: do single-cell CEQ and apply results
   !    ! Cells needing clustering: BFS flood-fill to form clusters, then process
   !    ! ========================================================================
   !    do index=1,vf%band_count(0)

   !       ! Get the interfacial cell indices
   !       i=vf%band_map(1,index)
   !       j=vf%band_map(2,index)
   !       k=vf%band_map(3,index)

   !       ! Skip if already processed as part of a cluster
   !       if (clustered(i,j,k)) cycle

   !       ! Initialize cluster with the seed cell
   !       n_clustered=1
   !       cell_indices(:,1)=[i,j,k]

   !       ! Initialize quantities
   !       itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
   !       vol_old=sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
   !       mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
   !       Y =sc%SC(i,j,k,1:ns)
   !       Tl=sc%SC(i,j,k,iTl)
   !       Tg=sc%SC(i,j,k,iTg)

   !       if (needs_clustering_flag(i,j,k).eq.1) then

   !          ! Increment cluster ID and mark seed cell
   !          cluster_id=cluster_id+1
   !          clustered(i,j,k)=.true.
   !          dbg_flg(i,j,k)=real(cluster_id,WP)

   !          ! Initialize the species mass
   !          do isc=1,ns
   !             p=sc%phase(isc)
   !             Y(isc)=sc%Prho(p)*sc%PVF(i,j,k,p)*cfg%vol(i,j,k)*sc%SC(i,j,k,isc)
   !          end do

   !          ! Initialize the mass-averaged temperatures (mass*T for accumulation)
   !          Tl=sc%Prho(Lphase)*sc%PVF(i,j,k,Lphase)*cfg%vol(i,j,k)*Tl
   !          Tg=sc%Prho(Gphase)*sc%PVF(i,j,k,Gphase)*cfg%vol(i,j,k)*Tg

   !          ! BFS flood-fill: explore connected interfacial neighbors
   !          ! Stop when cluster liquid fraction is sufficient or max size is reached
   !          bfs_head=1
   !          cluster_done=.false.
   !          do while (bfs_head.le.n_clustered.and..not.cluster_done)

   !             ! Pop the next cell from the BFS queue
   !             i=cell_indices(1,bfs_head)
   !             j=cell_indices(2,bfs_head)
   !             k=cell_indices(3,bfs_head)
   !             bfs_head=bfs_head+1

   !             ! Explore all neighbors in the stencil
   !             z_loop: do kn=k-stz,k+stz
   !                y_loop: do jn=j-sty,j+sty
   !                   x_loop: do in=i-stx,i+stx

   !                      ! Check stopping criteria before adding more cells
   !                      if (n_clustered.ge.nc_max) then
   !                         cluster_done=.true.
   !                         exit z_loop
   !                      end if

   !                      ! Neighbor must be interfacial, not clustered yet, and also needs clustering
   !                      if (vf%VF(in,jn,kn).gt.VFlo.and.vf%VF(in,jn,kn).lt.VFhi.and..not.clustered(in,jn,kn).and.needs_clustering_flag(in,jn,kn).eq.1) then

   !                         ! Mark it as clustered
   !                         n_clustered=n_clustered+1

   !                         ! Grow the array if needed
   !                         if (n_clustered.gt.nc_cap) then
   !                            allocate(cell_indices_tmp(3,nc_cap*2))
   !                            cell_indices_tmp(:,1:nc_cap)=cell_indices
   !                            nc_cap=nc_cap*2
   !                            call move_alloc(cell_indices_tmp,cell_indices)
   !                         end if

   !                         cell_indices(:,n_clustered)=[in,jn,kn]
   !                         clustered(in,jn,kn)=.true.
   !                         dbg_flg(in,jn,kn)=real(cluster_id,WP)

   !                         ! Accumulate old volumes
   !                         vol_old=vol_old+sc%PVF(in,jn,kn,:)*cfg%vol(in,jn,kn)

   !                         ! Accumulate mass and mass*temperature
   !                         do isc=1,ns
   !                            p=sc%phase(isc)
   !                            Y(isc)=Y(isc)+sc%Prho(p)*sc%PVF(in,jn,kn,p)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,isc)
   !                         end do
   !                         Tln=sc%SC(in,jn,kn,iTl)
   !                         Tgn=sc%SC(in,jn,kn,iTg)
   !                         Tl=Tl+sc%Prho(Lphase)*sc%PVF(in,jn,kn,Lphase)*cfg%vol(in,jn,kn)*Tln
   !                         Tg=Tg+sc%Prho(Gphase)*sc%PVF(in,jn,kn,Gphase)*cfg%vol(in,jn,kn)*Tgn

   !                         ! Accumulate interface area
   !                         itf_area=itf_area+cfg%vol(in,jn,kn)*vf%SD(in,jn,kn)

   !                         ! Check if cluster now has sufficient liquid fraction
   !                         cluster_liq_frac=vol_old(Lphase)/sum(vol_old)
   !                         if (cluster_liq_frac.ge.ceq_liq_frac_thld) then
   !                            cluster_done=.true.
   !                            exit z_loop
   !                         end if
   !                      end if

   !                   end do x_loop
   !                end do y_loop
   !             end do z_loop

   !          end do

   !          ! If only the seed cell remains, this is not a real cluster
   !          if (n_clustered.eq.1) then
   !             dbg_flg(cell_indices(1,1),cell_indices(2,1),cell_indices(3,1))=0.0_WP
   !             cluster_id=cluster_id-1
   !          end if

   !          ! Cluster-level phase masses
   !          mp=sc%Prho*vol_old

   !          ! Cluster-averaged mass fractions and temperatures
   !          do isc=1,ns
   !             Y(isc)=Y(isc)/mp(sc%phase(isc))
   !          end do
   !          Tl=Tl/mp(Lphase)
   !          Tg=Tg/mp(Gphase)

   !       end if

   !       ! Get the equilibrium state
   !       call get_equilibrium()

   !       ! Need to change the following line
   !       if (.not.state%success) cycle

   !       ! Store the old total volume
   !       Vold=sum(vol_old)

   !       ! Calculate the cluster VOF
   !       vof=vol_old(Lphase)/Vold

   !       ! Update the phase masses
   !       mp=0.0_WP
   !       do isc=1,ns
   !          p=sc%phase(isc)
   !          mp(p)=mp(p)+N(isc)*MM(isc)
   !       end do

   !       ! Get the phase volumes
   !       vol_new=mp/sc%Prho
   !       Vnew=sum(vol_new)

   !       ! Get the phase change mass flux
   !       mdot2p=(Vnew-Vold)/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*itf_area)

   !       ! Allocate per-cluster work arrays
   !       allocate(Vscaled(n_clustered)); Vscaled=0.0_WP
   !       allocate(vof_old(n_clustered)); vof_old=0.0_WP
   !       allocate(vof_new(n_clustered)); vof_new=0.0_WP
   !       allocate(w(n_clustered)); w=0.0_WP
   !       allocate(active(n_clustered)); active=.false.

   !       ! --- Gradient-preserving: store original per-cell T and Y ---
   !       allocate(Tl_orig(n_clustered))
   !       allocate(Tg_orig(n_clustered))
   !       allocate(Y_orig(ns,n_clustered))
   !       allocate(Y_eq(ns))
   !       allocate(delta_Y(ns))
   !       do m=1,n_clustered
   !          i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !          Tl_orig(m)=sc%SC(i,j,k,iTl)
   !          Tg_orig(m)=sc%SC(i,j,k,iTg)
   !          do isc=1,ns
   !             Y_orig(isc,m)=sc%SC(i,j,k,isc)
   !          end do
   !       end do

   !       ! Compute old cluster-mean temperature (mass-weighted) and species
   !       Tl_old_mean=Tl   ! These are already the cluster-averaged values from above
   !       Tg_old_mean=Tg
   !       ! Compute deltas: new equilibrium vs old cluster mean
   !       delta_Tl=state%T-Tl_old_mean
   !       delta_Tg=state%T-Tg_old_mean
   !       ! Compute equilibrium species mass fractions
   !       do isc=1,ns
   !          p=sc%phase(isc)
   !          if (mp(p).gt.0.0_WP) then
   !             Y_eq(isc)=MM(isc)*N(isc)/mp(p)
   !          else
   !             Y_eq(isc)=0.0_WP
   !          end if
   !       end do
   !       ! Delta for species mass fractions
   !       delta_Y=Y_eq-Y

   !       ! Gather geometry and current VOF per clustered cell
   !       do m=1,n_clustered
   !          i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !          active(m)=.true.
   !          Vscaled(m)=cfg%vol(i,j,k)/Vref ! Scale it for more accurate calculations
   !          vof_old(m)=vf%VF(i,j,k)
   !          vof_new(m)=vof_old(m)
   !       end do

   !       ! Total liquid volume change ( > 0 condensation, < 0 vaporization)
   !       dVl=(vol_new(Lphase)-vol_old(Lphase))/Vref ! Scale it for more accurate calculations
   !       dVl_rem=dVl

   !       if (abs(dVl).gt.dVlmin) then

   !          ! Build weights
   !          do m=1,n_clustered
   !             i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !             itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
   !             interfaceness=minval(sc%PVF(i,j,k,:))
   !             w(m)=max(itf_area*interfaceness,wmin)
   !          end do

   !          ! Iteratively redistribute liquid
   !          do

   !             ! Update weights sum
   !             wsum=0.0_WP
   !             do m=1,n_clustered
   !                if (active(m)) wsum=wsum+w(m)
   !             end do

   !             ! Terminate if succeeded
   !             if (abs(dVl_rem).le.dVlmin.or.wsum.le.0.0_WP) exit

   !             ! Distribute
   !             dVl=dVl_rem
   !             do m=1,n_clustered
   !                if (.not.active(m)) cycle

   !                ! Estimate VOF
   !                dVl_i =(w(m)/wsum)*dVl
   !                vof_tmp=vof_new(m)+dVl_i/Vscaled(m)

   !                ! Clip it
   !                if (vof_tmp.gt.1.0_WP) then
   !                   dVl_i=(1.0_WP-vof_new(m))*Vscaled(m)
   !                   vof_new(m)=1.0_WP
   !                   active(m)=.false.
   !                else if (vof_tmp.lt.0.0_WP) then
   !                   dVl_i=(0.0_WP-vof_new(m))*Vscaled(m)
   !                   vof_new(m)=0.0_WP
   !                   active(m)=.false.
   !                else
   !                   vof_new(m)=vof_tmp
   !                end if

   !                ! Correct the liquid volume change
   !                dVl_rem=dVl_rem-dVl_i

   !             end do

   !          end do

   !       end if

   !       ! Assign per-cell fields (Need to treat cells with VOF=0 and 1, differently)
   !       do m=1,n_clustered

   !          ! Get the cell indices
   !          i=cell_indices(1,m)
   !          j=cell_indices(2,m)
   !          k=cell_indices(3,m)

   !          ! Skip ghost cells (only write to owned cells; ghosts will be overwritten by sync)
   !          if (i.lt.cfg%imin_.or.i.gt.cfg%imax_) cycle
   !          if (j.lt.cfg%jmin_.or.j.gt.cfg%jmax_) cycle
   !          if (k.lt.cfg%kmin_.or.k.gt.cfg%kmax_) cycle

   !          ! Assign VOF
   !          if (vof_new(m).lt.VFlo) then
   !             vf%VF(i,j,k)=0.0_WP
   !          else if (vof_new(m).gt.VFhi) then
   !             vf%VF(i,j,k)=1.0_WP
   !          else
   !             vf%VF(i,j,k)=vof_new(m)
   !          end if
   !          sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
   !          sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)

   !          ! Composition: gradient-preserving with additive delta
   !          do isc=1,ns
   !             p=sc%phase(isc)
   !             if(sc%PVF(i,j,k,p).gt.0.0_WP) then
   !                sc%SC(i,j,k,isc)=Y_orig(isc,m)+delta_Y(isc)
   !             else
   !                sc%SC(i,j,k,isc)=0.0_WP
   !             end if
   !          end do
   !          ! Renormalize species per phase to ensure mass conservation
   !          do p=Lphase,Gphase
   !             Y_sum=0.0_WP
   !             do isc=1,ns
   !                if (sc%phase(isc).eq.p) Y_sum=Y_sum+sc%SC(i,j,k,isc)
   !             end do
   !             if (Y_sum.gt.0.0_WP) then
   !                do isc=1,ns
   !                   if (sc%phase(isc).eq.p) sc%SC(i,j,k,isc)=sc%SC(i,j,k,isc)/Y_sum
   !                end do
   !             end if
   !          end do

   !          ! Temperature and phase change mass flux
   !          if (vf%VF(i,j,k).eq.1.0_WP) then
   !             sc%SC(i,j,k,iTl)=Tl_orig(m)+delta_Tl
   !             sc%SC(i,j,k,iTg)=0.0_WP
   !             lg%mdot2p(i,j,k)=0.0_WP
   !          else if (vf%VF(i,j,k).eq.0.0_WP) then
   !             sc%SC(i,j,k,iTl)=0.0_WP
   !             sc%SC(i,j,k,iTg)=Tg_orig(m)+delta_Tg
   !             lg%mdot2p(i,j,k)=0.0_WP
   !          else
   !             sc%SC(i,j,k,iTl)=Tl_orig(m)+delta_Tl
   !             sc%SC(i,j,k,iTg)=Tg_orig(m)+delta_Tg
   !             lg%mdot2p(i,j,k)=mdot2p
   !          end if

   !       end do

   !       ! Free per-cluster work arrays
   !       deallocate(Vscaled,vof_old,vof_new,w,active)
   !       deallocate(Tl_orig,Tg_orig,Y_orig,Y_eq,delta_Y)

   !    end do

   !    ! Sync VOF
   !    call cfg%sync(vf%VF)

   !    ! Remove flotsams and thin structures if needed
   !    call vf%remove_flotsams()
   !    call vf%remove_thinstruct()

   !    ! Synchronize and clean-up barycenter fields
   !    call vf%sync_and_clean_barycenters()

   !    ! Update the interface band
   !    call vf%update_band()

   !    ! Perform interface reconstruction from transported moments
   !    call vf%build_interface()

   !    ! Create discontinuous polygon mesh from IRL interface
   !    call vf%polygonalize_interface()

   !    ! Calculate curvature
   !    call vf%get_curvature()

   !    ! Reset moments to guarantee compatibility with interface reconstruction
   !    call vf%reset_moments()

   !    ! Sync fields
   !    do isc=1,sc%nscalar
   !       call cfg%sync(sc%SC(:,:,:,isc))
   !    end do
   !    call cfg%sync(sc%PVF(:,:,:,Lphase))
   !    call cfg%sync(sc%PVF(:,:,:,Gphase))
   !    call cfg%sync(lg%mdot2p)

   !    ! Apply boundary conditions
   !    call sc%apply_bcond(time%t,time%dt)
   !    call vf%apply_bcond(time%t,time%dt)

   !    ! Debug: dbg_flg already contains per-cluster IDs from the loop above
   !    call cfg%sync(dbg_flg)

   !    ! Deallocate arrays
   !    deallocate(vol_new,vol_old,mp,N,phasicHoR,Y,clustered,needs_clustering_flag,cell_indices)


   ! contains

   !    subroutine get_equilibrium()
   !       implicit none

   !       ! Calculate and normalize the mole numbers
   !       do isc=1,ns
   !          N(isc)=Y(isc)*mp(sc%phase(isc))/MM(isc)
   !       end do
   !       Nsum=sum(N)
   !       if (Nsum.gt.0.0_WP) N=N/Nsum

   !       ! Get the phasic enthalpies
   !       call state%get_phasic_HoR(Lphase,N,Tl,phasicHoR(Lphase))
   !       call state%get_phasic_HoR(Gphase,N,Tg,phasicHoR(Gphase))

   !       ! Reinitialize the mole numbers
   !       call state%N_init(N=N,HoR=sum(phasicHoR),T_g=T_g)
   !       if (.not.state%success) then
   !          print*,'Cluster VOF = ',vof
   !          print*,'N = ',N
   !          print*,'N*Nsum = ',N*Nsum
   !          print*,'HoR = ',sum(phasicHoR)
   !          print*,'T_g = ',T_g
   !          print*,'Clustered cells info:'
   !          print*,'n_clustered = ',n_clustered
   !          do m=1,n_clustered
   !             i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !             print*,'i,j,k = ',i,j,k
   !             print*,'VOF = ',vf%VF(i,j,k)
   !             print*,'SD = ',vf%SD(i,j,k)
   !          end do
   !          call die('interface_jump: N_init failed')
   !       end if

   !       ! Get the chemical equilibrium
   !       call state%equilibrate()

   !       ! Re-scale the mole numbers
   !       N=state%N*Nsum

   !    end subroutine get_equilibrium

   ! end subroutine interface_jump


   ! Old version
   ! subroutine interface_jump()
   !    use messager, only: die
   !    implicit none
   !    real(WP), dimension(:),     allocatable :: vol_new,vol_old,mp,N,phasicHoR,Y
   !    logical,  dimension(:,:,:), allocatable :: clustered
   !    logical,  dimension(:),     allocatable :: active
   !    integer,  dimension(:,:),   allocatable :: cell_indices
   !    real(WP), dimension(:),     allocatable :: Vscaled,vof_old,vof_new,w
   !    real(WP) :: Vnew,Vold,Nsum,vof,itf_area,Tl,Tg,Tln,Tgn
   !    integer  :: i,j,k,index,isc,p,n_clustered,m
   !    integer  :: in,jn,kn
   !    integer  :: stx,sty,stz
   !    real(WP) :: mdot2p
   !    real(WP), parameter :: wmin=1.0e-16_WP,dVlmin=1.0e-16_WP
   !    integer,  parameter :: nc_max=27
   !    real(WP) :: dVl,dVl_i,dVl_rem,Vref,vof_tmp,interfaceness,wsum

   !    ! Debug
   !    dbg_flg=0.0_WP

   !    Vref=minval(cfg%vol)

   !    ! Allocate arrays
   !    allocate(vol_new(Lphase:Gphase))
   !    allocate(vol_old(Lphase:Gphase))
   !    allocate(mp(Lphase:Gphase))
   !    allocate(N(ns))
   !    allocate(phasicHoR(Lphase:Gphase))
   !    allocate(Y(ns))
   !    allocate(clustered(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); clustered=.false.
   !    allocate(cell_indices(3,nc_max)); cell_indices=0
   !    allocate(Vscaled(nc_max)); Vscaled=0.0_WP
   !    allocate(vof_old(nc_max)); vof_old=0.0_WP
   !    allocate(vof_new(nc_max)); vof_new=0.0_WP
   !    allocate(w(nc_max)); w=0.0_WP
   !    allocate(active(nc_max)); active=.false.

   !    ! Clustering stencil
   !    if (cfg%nx.gt.1) then
   !       stx=1
   !    else
   !       stx=0
   !    end if
   !    if (cfg%ny.gt.1) then
   !       sty=1
   !    else
   !       sty=0
   !    end if
   !    if (cfg%nz.gt.1) then
   !       stz=1
   !    else
   !       stz=0
   !    end if

   !    ! if (cfg%amRoot) print*,'beginning of jump: VOF(31,12,1) = ',vf%VF(31,12,1)

   !    ! Loop over the interfacial cells
   !    do index=1,vf%band_count(0)

   !       ! Get the interfacial cell indices
   !       i=vf%band_map(1,index)
   !       j=vf%band_map(2,index)
   !       k=vf%band_map(3,index)

   !       ! Skip if already clustered
   !       if (clustered(i,j,k)) cycle

   !       ! Add current cell to the potential cluster
   !       n_clustered=1
   !       cell_indices(:,1)=[i,j,k]

   !       ! Initialize the interfacial area and old volumes
   !       itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
   !       vol_old=sc%PVF(i,j,k,:)*cfg%vol(i,j,k)

   !       ! Pre-evaluate the equilibrium
   !       mp=sc%Prho*sc%PVF(i,j,k,:)*cfg%vol(i,j,k)
   !       Y =sc%SC(i,j,k,1:ns)
   !       Tl=sc%SC(i,j,k,iTl)
   !       Tg=sc%SC(i,j,k,iTg)

   !       ! if (i.eq.31.and.j.eq.12.and.k.eq.1) then
   !       !    print*,'*******'
   !       !    print*,'Before CEQ:'
   !       !    print*,'VOF = ',vf%VF(i,j,k)
   !       !    print*,'Tg = ',Tg
   !       !    print*,'Tl = ',Tl
   !       !    print*,'Y  = ',Y
   !       !    print*,'mp = ',mp
   !       !    print*,'vol_old = ',vol_old
   !       !    print*,'Vold = ',sum(vol_old)
   !       ! end if
   !       call get_equilibrium()
   !       ! Debug
   !       ! if (.not.state%Nming_success) then
   !       !       print*,'*******'
   !       !       print*,'This cell failed ceq get_Nming:'
   !       !       print*,'i,j,k = ',i,j,k
   !       !       print*,'x, y, z = ',cfg%xm(i),cfg%ym(j),cfg%zm(k)
   !       !       print*,'Tg = ',Tg
   !       !       print*,'Tl = ',Tl
   !       !       print*,'Y = ',Y
   !       !       print*,'mp = ',mp
   !       !       print*,'PVF = ',sc%PVF(i,j,k,:)
   !       !       print*,'VF = ',vf%VF(i,j,k)
   !       !       print*,'N = ',N
   !       !       print*,'state%N = ',state%N
   !       !       print*,'state%HoR0 = ',state%HoR0
   !       !       print*,'*******'
   !       !       call die('')
   !       ! end if

   !       ! Cluster cells if the equilibrium failed
   !       if ((.not.state%success).or.(state%N(iWl)/(state%N(iWl)+state%N(iWv)).lt.5e-2)) then

   !          ! if (vf%VF(i,j,k).gt.0.1_WP) then
   !          !    print*,'*******'
   !          !    print*,'This high VOF cell failed ceq:'
   !          !    print*,'i,j,k = ',i,j,k
   !          !    print*,'Tg = ',Tg
   !          !    print*,'Tl = ',Tl
   !          !    print*,'Y = ',Y
   !          !    print*,'mp = ',mp
   !          !    print*,'PVF = ',sc%PVF(i,j,k,:)
   !          !    print*,'VF = ',vf%VF(i,j,k)
   !          !    print*,'N = ',N
   !          !    print*,'state%N = ',state%N
   !          !    print*,'state%HoR0 = ',state%HoR0
   !          !    print*,'*******'
   !          ! end if

   !          ! Mark it as clustered
   !          clustered(i,j,k)=.true.

   !          ! Initialize the species mass
   !          do isc=1,ns
   !             p=sc%phase(isc)
   !             Y(isc)=sc%Prho(p)*sc%PVF(i,j,k,p)*cfg%vol(i,j,k)*sc%SC(i,j,k,isc)
   !          end do

   !          ! Initialize the mass-averaged temperatures
   !          Tl=sc%Prho(Lphase)*sc%PVF(i,j,k,Lphase)*cfg%vol(i,j,k)*Tl
   !          Tg=sc%Prho(Gphase)*sc%PVF(i,j,k,Gphase)*cfg%vol(i,j,k)*Tg

   !          ! Loop over the cluster stencil skipping the ghost cells
   !          z_loop: do kn=k-stz,k+stz
   !             if (kn.lt.cfg%kmin_.or.kn.gt.cfg%kmax_) cycle
   !             y_loop: do jn=j-sty,j+sty
   !                if (jn.lt.cfg%jmin_.or.jn.gt.cfg%jmax_) cycle
   !                x_loop: do in=i-stx,i+stx
   !                   if (in.lt.cfg%imin_.or.in.gt.cfg%imax_) cycle

   !                   ! Neighbor must be interfacial and not clustered yet
   !                   if (vf%VF(in,jn,kn).gt.VFlo.and.vf%VF(in,jn,kn).lt.VFhi.and..not.clustered(in,jn,kn)) then

   !                      ! Mark it as clustered
   !                      n_clustered=n_clustered+1
   !                      cell_indices(:,n_clustered)=[in,jn,kn]
   !                      clustered(in,jn,kn)=.true.

   !                      ! Accumulate old volumes
   !                      vol_old=vol_old+sc%PVF(in,jn,kn,:)*cfg%vol(in,jn,kn)

   !                      ! Accumulate mass*SC and mass*temperature
   !                      do isc=1,ns
   !                         p=sc%phase(isc)
   !                         Y(isc)=Y(isc)+sc%Prho(p)*sc%PVF(in,jn,kn,p)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,isc)
   !                      end do
   !                      ! Tl=Tl+sc%Prho(Lphase)*sc%PVF(in,jn,kn,Lphase)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,iTl)
   !                      ! Tg=Tg+sc%Prho(Gphase)*sc%PVF(in,jn,kn,Gphase)*cfg%vol(in,jn,kn)*sc%SC(in,jn,kn,iTg)
   !                      Tln=sc%SC(in,jn,kn,iTl)
   !                      Tgn=sc%SC(in,jn,kn,iTg)
   !                      Tl=Tl+sc%Prho(Lphase)*sc%PVF(in,jn,kn,Lphase)*cfg%vol(in,jn,kn)*Tln
   !                      Tg=Tg+sc%Prho(Gphase)*sc%PVF(in,jn,kn,Gphase)*cfg%vol(in,jn,kn)*Tgn

   !                      ! Accumulate interface area
   !                      itf_area=itf_area+cfg%vol(in,jn,kn)*vf%SD(in,jn,kn)
   !                   end if

   !                end do x_loop
   !             end do y_loop
   !          end do z_loop

   !          ! Cluster-level phase masses
   !          mp=sc%Prho*vol_old

   !          ! Cluster-averaged mass fractions and temperatures
   !          do isc=1,ns
   !             Y(isc)=Y(isc)/mp(sc%phase(isc))
   !          end do
   !          Tl=Tl/mp(Lphase)
   !          Tg=Tg/mp(Gphase)

   !          ! Get the equilibrium state of the cluster
   !          call get_equilibrium()

   !       end if

   !       ! Calculate the cluster VOF
   !       vof=vol_old(Lphase)/sum(vol_old)

   !       ! Store the old total volume
   !       Vold=sum(vol_old)

   !       ! Debug
   !       if (.not.state%success) then
   !          print*,'Cluster VOF = ',vof
   !          print*,'N initial scaled and fed into ceq = ',N
   !          print*,'N initial actual = ',N*Nsum
   !          print*,'Y = ',Y
   !          print*,'HoR = ',sum(phasicHoR)
   !          print*,'T_g = ',T_g
   !          print*,'Tg = ',Tg
   !          print*,'Tl = ',Tl
   !          print*,'Clustered cells info:'
   !          print*,'n_clustered = ',n_clustered
   !          do m=1,n_clustered
   !             i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !             print*,'i,j,k = ',i,j,k
   !             print*,'VOF = ',vf%VF(i,j,k)
   !             print*,'SD = ',vf%SD(i,j,k)
   !          end do
   !          call die('line 685')
   !       end if

   !       ! Update the phase masses
   !       mp=0.0_WP
   !       do isc=1,ns
   !          p=sc%phase(isc)
   !          mp(p)=mp(p)+N(isc)*MM(isc)
   !       end do

   !       ! Get the phase volumes
   !       vol_new=mp/sc%Prho
   !       Vnew=sum(vol_new)

   !       ! Get the phase change mass flux
   !       mdot2p=(Vnew-Vold)/(time%dt*(1.0_WP/sc%Prho(Gphase)-1.0_WP/sc%Prho(Lphase))*itf_area)

   !       ! if (i.eq.16.and.j.eq.27.and.k.eq.12) then
   !       !    print*,'mp = ',mp
   !       !    print*,'HoR = ',phasicHoR
   !       !    print*,'sum(HoR) = ',sum(phasicHoR)
   !       !    print*,'vol_new = ',vol_new
   !       !    print*,'Vnew = ',Vnew
   !       !    print*,'itf_area = ',itf_area
   !       !    print*,'dV = ',Vnew-Vold
   !       !    print*,'mdot2p = ',mdot2p
   !       ! end if

   !       ! Gather geometry and current VOF per clustered cell
   !       do m=1,n_clustered
   !          i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !          active(m)=.true.
   !          Vscaled(m)=cfg%vol(i,j,k)/Vref ! Scale it for more accurate calculations
   !          vof_old(m)=vf%VF(i,j,k)
   !          vof_new(m)=vof_old(m)
   !       end do

   !       ! Total liquid volume change ( > 0 condensation, < 0 vaporization)
   !       dVl=(vol_new(Lphase)-vol_old(Lphase))/Vref ! Scale it for more accurate calculations
   !       dVl_rem=dVl

   !       if (abs(dVl).gt.dVlmin) then

   !          ! Build weights
   !          do m=1,n_clustered
   !             i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !             itf_area=cfg%vol(i,j,k)*vf%SD(i,j,k)
   !             interfaceness=minval(sc%PVF(i,j,k,:))
   !             w(m)=max(itf_area*interfaceness,wmin)
   !          end do

   !          ! Iteratively redistribute liquid
   !          do

   !             ! Update weights sum
   !             wsum=0.0_WP
   !             do m=1,n_clustered
   !                if (active(m)) wsum=wsum+w(m)
   !             end do

   !             ! Terminate if succssesd
   !             if (abs(dVl_rem).le.dVlmin.or.wsum.le.0.0_WP) exit

   !             ! Distribute
   !             dVl=dVl_rem
   !             do m=1,n_clustered
   !                if (.not.active(m)) cycle

   !                ! Estimate VOF
   !                dVl_i =(w(m)/wsum)*dVl
   !                vof_tmp=vof_new(m)+dVl_i/Vscaled(m)

   !                ! Clip it
   !                if (vof_tmp.gt.1.0_WP) then
   !                   dVl_i=(1.0_WP-vof_new(m))*Vscaled(m)
   !                   vof_new(m)=1.0_WP
   !                   active(m)=.false.
   !                else if (vof_tmp.lt.0.0_WP) then
   !                   dVl_i=(0.0_WP-vof_new(m))*Vscaled(m)
   !                   vof_new(m)=0.0_WP
   !                   active(m)=.false.
   !                else
   !                   vof_new(m)=vof_tmp
   !                end if

   !                ! Correct the liquid volume change
   !                dVl_rem=dVl_rem-dVl_i

   !             end do

   !          end do

   !       end if

   !       ! Assign per-cell fields (Need to treat cells with VOF=0 and 1, differently)
   !       do m=1,n_clustered

   !          ! Get the cell indices
   !          i=cell_indices(1,m)
   !          j=cell_indices(2,m)
   !          k=cell_indices(3,m)

   !          ! Assign VOF
   !          if (vof_new(m).lt.VFlo) then
   !             vf%VF(i,j,k)=0.0_WP
   !          else if (vof_new(m).gt.VFhi) then
   !             vf%VF(i,j,k)=1.0_WP
   !          else
   !             vf%VF(i,j,k)=vof_new(m)
   !          end if
   !          sc%PVF(i,j,k,Lphase)=vf%VF(i,j,k)
   !          sc%PVF(i,j,k,Gphase)=1.0_WP-vf%VF(i,j,k)

   !          ! Composition (Assuming the same mass fraction for all non-empty the cells in the cluster)
   !          do isc=1,ns
   !             p=sc%phase(isc)
   !             if(sc%PVF(i,j,k,p).gt.0.0_WP) then
   !                sc%SC(i,j,k,isc)=MM(isc)*N(isc)/mp(sc%phase(isc))
   !             else
   !                sc%SC(i,j,k,isc)=0.0_WP
   !             end if
   !          end do

   !          ! Temperature and phase change mass flux
   !          if (vf%VF(i,j,k).eq.1.0_WP) then
   !             sc%SC(i,j,k,iTl)=state%T ! Not sure if this is good enough.
   !             sc%SC(i,j,k,iTg)=0.0_WP
   !             lg%mdot2p(i,j,k)=0.0_WP
   !          else if (vf%VF(i,j,k).eq.0.0_WP) then
   !             sc%SC(i,j,k,iTl)=0.0_WP
   !             sc%SC(i,j,k,iTg)=state%T ! Not sure if this is good enough.
   !             lg%mdot2p(i,j,k)=0.0_WP
   !          else
   !             sc%SC(i,j,k,iTl)=state%T
   !             sc%SC(i,j,k,iTg)=state%T
   !             lg%mdot2p(i,j,k)=mdot2p
   !          end if

   !          ! if (i.eq.16.and.j.eq.27.and.k.eq.12) then
   !          !    print*,'VOF = ',vf%VF(i,j,k)
   !          ! end if

   !       end do

   !    end do

   !    ! Sync VOF
   !    call cfg%sync(vf%VF)

   !    ! if (cfg%amRoot) print*,'after cluster operations and before the IRL stuff: VOF(31,12,1) = ',vf%VF(31,12,1)

   !    ! Update the interface (Do I need it? I don't think so)
   !    ! call vf%advect_interface(0.0_WP,fs%U,fs%V,fs%W)

   !    ! Remove flotsams and thin structures if needed
   !    call vf%remove_flotsams()
   !    call vf%remove_thinstruct()

   !    ! Synchronize and clean-up barycenter fields
   !    call vf%sync_and_clean_barycenters()

   !    ! Update the interface band
   !    call vf%update_band()

   !    ! Perform interface reconstruction from transported moments
   !    call vf%build_interface()

   !    ! Create discontinuous polygon mesh from IRL interface
   !    call vf%polygonalize_interface()

   !    ! Perform interface sensing (Do I need it?)
   !    ! if (vf%two_planes) call vf%sense_interface()

   !    ! Calculate distance from polygons (I don't think it's needed anywhere)
   !    ! call vf%distance_from_polygon()

   !    ! Calculate subcell phasic volumes (I don't think it's needed anywhere)
   !    ! call vf%subcell_vol()

   !    ! Calculate curvature
   !    call vf%get_curvature()

   !    ! Reset moments to guarantee compatibility with interface reconstruction
   !    call vf%reset_moments()

   !    ! if (cfg%amRoot) print*,'after the IRL stuff: VOF(31,12,1) = ',vf%VF(31,12,1)

   !    ! Sync fields
   !    do isc=1,sc%nscalar
   !       call cfg%sync(sc%SC(:,:,:,isc))
   !    end do
   !    call cfg%sync(vf%VF)
   !    call cfg%sync(lg%mdot2p)

   !    ! Apply boundary conditions
   !    call sc%apply_bcond(time%t,time%dt)
   !    call vf%apply_bcond(time%t,time%dt)

   !    ! Debug
   !    where (clustered) dbg_flg=1.0_WP
   !    call cfg%sync(dbg_flg)

   !    ! Deallocate arrays
   !    deallocate(vol_new,vol_old,mp,N,phasicHoR,Y,clustered,cell_indices,Vscaled,vof_old,vof_new,w,active)

   !    ! Debug
   !    ! debug: block
   !    !    use messager, only: die
   !    !    use irl_fortran_interface, only: getNumberOfVertices
   !    !    if (cfg%amRoot) print*,'NumberOfVertices for cell 31,12,1',getNumberOfVertices(vf%interface_polygon(1,31,12,1))
   !    !    do k=cfg%kmin_,cfg%kmax_
   !    !       do j=cfg%jmin_,cfg%jmax_
   !    !          do i=cfg%imin_,cfg%imax_
   !    !             if (vf%VF(i,j,k) > 0.0_WP .and. vf%VF(i,j,k) < 1.0_WP) then
   !    !                if (getNumberOfVertices(vf%interface_polygon(1,i,j,k)) == 0) then
   !    !                   print *, "empty poly:", i,j,k, "VF=",vf%VF(i,j,k), "mask=",vf%mask(i,j,k)
   !    !                   call die('')
   !    !                end if
   !    !             end if
   !    !          end do
   !    !       end do
   !    !    end do
   !    ! end block debug

   !    ! if (cfg%amRoot) print*,'End of the interface jump: VOF(31,12,1) = ',vf%VF(31,12,1)

   ! contains

   !    subroutine get_equilibrium()
   !       implicit none

   !       ! Calculate and normalize the mole numbers
   !       do isc=1,ns
   !          N(isc)=Y(isc)*mp(sc%phase(isc))/MM(isc)
   !       end do
   !       Nsum=sum(N)
   !       if (Nsum.gt.0.0_WP) N=N/Nsum

   !       ! Get the phasic enthalpies
   !       call state%get_phasic_HoR(Lphase,N,Tl,phasicHoR(Lphase))
   !       call state%get_phasic_HoR(Gphase,N,Tg,phasicHoR(Gphase))

   !       ! Reinitialize the mole numbers
   !       call state%N_init(N=N,HoR=sum(phasicHoR),T_g=T_g)
   !       if (.not.state%success) then
   !          print*,'Cluster VOF = ',vof
   !          print*,'N = ',N
   !          print*,'N*Nsum = ',N*Nsum
   !          print*,'HoR = ',sum(phasicHoR)
   !          print*,'T_g = ',T_g
   !          print*,'Clustered cells info:'
   !          print*,'n_clustered = ',n_clustered
   !          do m=1,n_clustered
   !             i=cell_indices(1,m); j=cell_indices(2,m); k=cell_indices(3,m)
   !             print*,'i,j,k = ',i,j,k
   !             print*,'VOF = ',vf%VF(i,j,k)
   !             print*,'SD = ',vf%SD(i,j,k)
   !          end do
   !          call die('line 943')
   !       end if

   !       ! Get the chemical equilibrium
   !       call state%equilibrate()

   !       ! if (i.eq.16.and.j.eq.27.and.k.eq.12) then
   !       !    print*,'N = ',N
   !       !    print*,'After CEQ:'
   !       !    print*,'N eq = ',state%N
   !       !    print*,'T eq = ',state%T
   !       !    print*,'T iter = ',state%iter_T
   !       ! end if

   !       ! Re-scale the mole numbers
   !       ! N=state%N*Nsum
   !       ! Debug: If not successful, don't assign N so I know what initial moles caused this
   !       if (state%success) N=state%N*Nsum

   !    end subroutine get_equilibrium

   ! end subroutine interface_jump


   !> Initialization of problem solver
   subroutine simulation_init
      use param,    only: param_exists,param_read,param_getsize
      use messager, only: die
      implicit none
      integer :: ne,ncs
      character(len=str_short), dimension(:), allocatable :: e_names
      real(WP), dimension(:,:), allocatable :: elem_mat
      real(WP), dimension(:,:), allocatable :: phse_mat
      real(WP), allocatable :: nasa_coef(:,:)
      character(len=str_medium), dimension(:), allocatable :: const_sp
      integer,  dimension(:), allocatable :: CS


      ! Read problem inputs
      read_inputs: block
         call param_read('Liquid density',rho_l)
         call param_read('Gas density',rho_g)
         call param_read('Latent heat',h_lg)
         call param_read('Liquid thermal conductivity',k_l)
         call param_read('Liquid specific heat capacity',Cp_l)
         call param_read('Gas thermal conductivity',k_g)
         call param_read('Gas specific heat capacity',Cp_g)
         call param_read('Liquid temperature',T_liq)
         call param_read('Ambient temperature',T_amb)
         call param_read('Pressure',Pressure)
         alpha_l=k_l/(rho_l*Cp_l)
         alpha_g=k_g/(rho_g*Cp_g)
         call param_read('CEQ liquid fraction threshold',ceq_liq_frac_thld,default=5.0e-2_WP)
      end block read_inputs


      ! Parse the mechanism file
      parse_mech: block
         use chem_sys_class, only: ncof
         use YAMLRead,       only: YAMLHandler,YAMLSequence,YAMLMap,yaml_open_file,yaml_start_from_sequence,yaml_close_file
         character(len=str_medium) :: mch_file
         character(len=str_short), dimension(:), allocatable :: sp_names_copy,const_sp_copy
         type(YAMLHandler)  :: domain
         type(YAMLSequence) :: sp_list,phases,elements
         type(YAMLElement)  :: sp,gas
         type(YAMLMap)      :: thermo,comp
         integer :: isc,nn,i,j,k,e,code
         character(len=:), allocatable :: name_arr(:)
         character(len=:), allocatable :: name
         real(WP), allocatable :: T_range(:)
         real(WP), dimension(:,:), allocatable :: a
         logical :: new_elem
         ! Get the target species from input
         ns=param_getsize('Species')
         if (param_exists('Constrained species')) then
            ncs=param_getsize('Constrained species')
            allocate(const_sp(1:ncs))
            allocate(const_sp_copy(1:ncs))
            call param_read('Constrained species',const_sp)
            const_sp_copy=const_sp
         else
            ncs=0
         end if
         allocate(sp_names(1:ns))
         allocate(sp_names_copy(1:ns))
         allocate(CS(ncs))
         call param_read('Species',sp_names)
         sp_names_copy=sp_names
         ! Read the mechanism file path
         call param_read('Mechanism file',mch_file)
         ! Open the mechanism
         domain=yaml_open_file(trim(mch_file))
         ! Get the list of all species
         sp_list=yaml_start_from_sequence(domain,'species')
         ! Extract the target species from the mechanism
         allocate(species(1:ns))
         nn=0
         k=0
         do isc=0,sp_list%size-1 ! Index in YAMLSequence starts from 0
            sp=sp_list%element(isc)
            name_arr=sp%value_str('name',code)
            name=''
            do i=1,size(name_arr)
               name=trim(name//name_arr(i))
            end do
            do i=1,ns
               if (sp_names_copy(i).eq.name) then
                  nn=nn+1
                  species(nn)=sp
                  sp_names(nn)=name
                  do j=1,ncs
                     if (const_sp_copy(j).eq.name) then
                        k=k+1
                        const_sp(k)=name
                        CS(k)=nn
                     end if
                  end do
               end if
            end do
            call sp%destroy()
         end do
         if(nn.ne.ns) call die('Some species are missing in the mechanism file.')
         ! Get the elements that exist in the target species
         phases=yaml_start_from_sequence(domain,'phases')
         gas=phases%element(0)
         elements=gas%value_sequence('elements',code) ! For some reason I couldn't directly get the element names using yaml-fortran
         allocate(e_names(elements%size)); e_names=''
         ne=0
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do e=1,size(comp%labels)
               name=''
               do i=1,size(comp%labels(e)%str)
                  name=name//trim(comp%labels(e)%str(i))
               end do
               new_elem=.true.
               do i=1,ne
                  if (trim(e_names(i)).eq.trim(name)) then
                     new_elem=.false.
                     exit
                  end if
               end do
               if (new_elem) then
                  ne=ne+1
                  e_names(ne)=trim(name)
               end if
            end do
         end do
         e_names=e_names(1:ne)
         ! Form the element matrix
         allocate(elem_mat(ns,ne)); elem_mat=0.0_WP
         do isc=1,ns
            sp=species(isc)
            comp=sp%value_map('composition')
            do i=1,size(comp%labels)
               name=''
               do nn=1,size(comp%labels(i)%str)
                  name=name//trim(comp%labels(i)%str(nn))
               end do
               do e=1,ne
                  if (trim(e_names(e)).eq.trim(name)) elem_mat(isc,e)=real(comp%value_int(name,code),WP)
               end do
            end do
         end do
         ! Read the NASA-7 polynomials
         allocate(nasa_coef(1:ns,2*ncof+1)); nasa_coef=0.0_WP
         allocate(a(1:2,1:ncof)); a=0.0_WP
         do isc=1,ns
            sp=species(isc)
            thermo=sp%value_map('thermo')
            T_range=thermo%value_double_1d('temperature-ranges',code)
            select case (size(T_range))
             case (3)
               a=thermo%value_double_2d('data',code)
             case (2)
               a(1,:)=thermo%value_double_1d('data',code)
             case default
               call die('Invalid temperature range')
            end select
            nasa_coef(isc,1)=T_range(2)
            nasa_coef(isc,2:  ncof+1)=a(1,:)
            nasa_coef(isc,9:2*ncof+1)=a(2,:)
         end do
         ! Form the phase summation matrix
         allocate(phse_mat(ns,Lphase+1:Gphase+1)); phse_mat(:,Lphase+1)=0.0_WP; phse_mat(:,Gphase+1)=1.0_WP
         do isc=1,ns
            if (len_trim(sp_names(isc)).ge.3) then
               if (sp_names(isc)(len_trim(sp_names(isc))-2:len_trim(sp_names(isc))).eq.'(L)') then
                  phse_mat(isc,Lphase+1)=1.0_WP
                  phse_mat(isc,Gphase+1)=0.0_WP
               end if
            end if
         end do
         ! Close the mechanism file and clean up
         call yaml_close_file(domain)
         call sp_list%destroy()
         call sp%destroy()
         call comp%destroy()
         call thermo%destroy()
         ! Destroy and deallocate species array (no longer needed after data extraction)
         do isc=1,ns
            call species(isc)%destroy()
         end do
         deallocate(species)
         deallocate(sp_names_copy,const_sp_copy)
         if (allocated(const_sp_copy)) deallocate(const_sp_copy)
      end block parse_mech


      ! Allocate work arrays
      allocate_work_arrays: block
         allocate(resSC (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_,1:ns+2))
         allocate(resU  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi    (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(T     (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(MM(ns)); MM=[32.0_WP,18.0_WP,18.0_WP,28.0_WP]; MM=0.001_WP*MM
         ! Debug
         allocate(cluster_map(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); cluster_map=0.0_WP
      end block allocate_work_arrays


      ! Initialize the chemical equilibrium framework
      ceq_init: block
         use chem_state_class, only: NR,LS
         integer :: ng=1
         real(WP), dimension(:,:), allocatable :: Bg
         character(len=2) :: eq_cond
         integer :: isc
         ! Allocate arrays
         allocate(Bg(ns,ng));  Bg=0.0_WP
         ! Create the general constraints
         do isc=1,ns
            if (sp_names(isc).eq.'H2O')    Bg(isc,1)=1.0_WP
            if (sp_names(isc).eq.'H2O(L)') Bg(isc,1)=1.0_WP
         end do
         ! Initialize the chemical system
         call sys%initialize(np=np,ns=ns,ne=ne,ncs=ncs,ng=ng,P=phse_mat,Ein=elem_mat,CS=CS,Bg=Bg,thermo_in=nasa_coef,diag=5)
         ! Initialize the chemical state
         call state%initialize(sys=sys,cond=fixed_PH,PH_method=NR,dNdT_method=LS,p=pressure)
         call param_read('Newton tolerance',state%tol_N)
         call param_read('Newton max iterations',state%iter_N_max)
         call param_read('T tolerance',state%tol_T)
         call param_read('T max iterations',state%iter_T_max)
         ! call param_read('Temperature initial guess',T_g)
         ! T_g=0.5_WP*(T_amb+T_liq)
         T_g=T_liq
         ! Deallocate arrays
         deallocate(Bg)
      end block ceq_init


      ! Initialize time tracker
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot,name='Main')
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         call param_read('Sub-iterations',time%itmax)
         time%dt=time%dtmax
         timeSC=timetracker(amRoot=cfg%amRoot,name='SC Time')
         ! call param_read('Scalar time step',timeSC%dtmax)
         timeSC%t=time%t
         timeSC%tmax=time%tmax
         timeSC%dtmax=time%dtmax
         timeSC%dt=timeSC%dtmax
      end block initialize_timetracker


      ! Initialize our VOF solver and field
      create_and_initialize_vof: block
         use mms_geom,  only: cube_refine_vol
         use vfs_class, only: lvira,flux_storage,neumann
         integer :: i,j,k,n,si,sj,sk,ierr
         real(WP), dimension(3,8) :: cube_vertex
         real(WP), dimension(3)   :: v_cent,a_cent
         real(WP) :: vol,area
         integer, parameter :: amr_ref_lvl=4
         ! Create a VOF solver
         call vf%initialize(cfg=cfg,reconstruction_method=lvira,transport_method=flux_storage,nband=6,name='VOF')
         ! Boundary conditinos
         call vf%add_bcond(name='xm',type=neumann,locator=xm_locator_sc,dir='-x')
         call vf%add_bcond(name='xp',type=neumann,locator=xp_locator   ,dir='+x')
         call vf%add_bcond(name='ym',type=neumann,locator=ym_locator_sc,dir='-y')
         call vf%add_bcond(name='yp',type=neumann,locator=yp_locator   ,dir='+y')
         call vf%add_bcond(name='zm',type=neumann,locator=zm_locator_sc,dir='-z')
         call vf%add_bcond(name='zp',type=neumann,locator=zp_locator   ,dir='+z')
         ! Initialize the VOF field
         call param_read('Drop center',center)
         call param_read('Drop radius',R0)
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
         ! Correct mask outside physical domain
         if (cfg%iproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imin_-1
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%iproc.eq.cfg%npx) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imax_+1,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.1) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmin_-1
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%jproc.eq.cfg%npy) then
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmax_+1,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.1) then
            do k=cfg%kmino_,cfg%kmin_-1
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         if (cfg%kproc.eq.cfg%npz) then
            do k=cfg%kmax_+1,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     vf%mask(i,j,k)=1
                  end do
               end do
            end do
         end if
         ! Apply boundary conditions
         call vf%apply_bcond(time%t,time%dt)
         ! Update the band
         call vf%update_band()
         ! Perform interface reconstruction from VOF field
         call vf%build_interface()
         ! Set interface planes at the boundaries
         call vf%set_full_bcond()
         ! Apply IRL boundary conditions
         call sym_irl()
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
         use hypre_str_class, only: gmres_smg,pcg_pfmg2,pcg_smg
         use tpns_class,      only: dirichlet,slip
         ! Create flow solver
         fs=tpns(cfg=cfg,name='Two-phase NS')
         ! Assign constant viscosity to each phase
         call param_read('Liquid dynamic viscosity',fs%visc_l)
         call param_read('Gas dynamic viscosity',fs%visc_g)
         ! Assign constant density to each phase
         fs%rho_l=rho_l
         fs%rho_g=rho_g
         ! Read in surface tension coefficient
         call param_read('Surface tension coefficient',fs%sigma)
         fs%contact_angle=fs%contact_angle*Pi/180.0_WP
         ! Assign acceleration of gravity
         call param_read('Gravity',fs%gravity)
         ! Boundary conditions
         call fs%add_bcond(name='xm',type=slip     ,face='x',dir=-1,canCorrect=.false.,locator=xm_locator)
         call fs%add_bcond(name='xp',type=dirichlet,face='x',dir=+1,canCorrect=.true. ,locator=xp_locator)
         call fs%add_bcond(name='ym',type=slip     ,face='y',dir=-1,canCorrect=.false.,locator=ym_locator)
         call fs%add_bcond(name='yp',type=dirichlet,face='y',dir=+1,canCorrect=.true. ,locator=yp_locator)
         call fs%add_bcond(name='zm',type=slip     ,face='z',dir=-1,canCorrect=.false.,locator=zm_locator)
         call fs%add_bcond(name='zp',type=dirichlet,face='z',dir=+1,canCorrect=.true. ,locator=zp_locator)
         ! Configure pressure solver
         ps=hypre_str(cfg=cfg,name='Pressure',method=pcg_smg,nst=7)
         call param_read('Pressure iteration',ps%maxit)
         call param_read('Pressure tolerance',ps%rcvg)
         call param_read('Max coarsening levels',ps%maxlevel)
         ! Implicit velocity solver
         vs=ddadi(cfg=cfg,name='Velocity',nst=7)
         ! vs=hypre_str(cfg=cfg,name='Velocity',method=pcg_smg,nst=7)
         ! call param_read('Scalar iteration',vs%maxit)
         ! call param_read('Scalar tolerance',vs%rcvg)
         ! Setup the solver
         call fs%setup(pressure_solver=ps,implicit_solver=vs)
         ! Initial field
         fs%U=0.0_WP; fs%V=0.0_WP; fs%W=0.0_WP
         ! Apply boundary conditions
         call fs%apply_bcond(time%t,time%dt)
         ! Calculate cell-centered velocities and divergence
         call fs%interp_vel(Ui,Vi,Wi)
         call fs%get_div()
      end block create_and_initialize_flow_solver


      ! Create a one-sided scalar solver
      create_scalar: block
         use param,           only: param_read
         use tpscalar_class,  only: bcond,dirichlet,neumann
         use mpi_f08,         only: MPI_ALLREDUCE,MPI_MAX
         use parallel,        only: MPI_REAL_WP
         use hypre_str_class, only: gmres_pfmg2
         type(bcond), pointer :: my_bc
         real(WP) :: mp(Lphase:Gphase),N_init(ns),spDiff,my_Y
         integer  :: n,i,j,k,isc,p,ierr
         integer  :: pos_open,pos_close
         ! Read-in inputs
         call param_read('Water vapor to air mole ratio',wv2air_rat)
         call param_read('Nitrogen to oxygen mole ratio',N2O_rat)
         call param_read('Species diffusivity',spDiff)
         ! Create scalar solver
         call sc%initialize(cfg=cfg,nscalar=ns+2,name='tpscalar')
         sc%skip(get_sp_ind('H2O(L)'))=.true.
         ! Boundary conditinos
         call sc%add_bcond(name='xm',type=neumann  ,locator=xm_locator_sc,dir='-x')
         call sc%add_bcond(name='xp',type=dirichlet,locator=xp_locator   ,dir='+x')
         call sc%add_bcond(name='ym',type=neumann  ,locator=ym_locator_sc,dir='-y')
         call sc%add_bcond(name='yp',type=dirichlet,locator=yp_locator   ,dir='+y')
         call sc%add_bcond(name='zm',type=neumann  ,locator=zm_locator_sc,dir='-z')
         call sc%add_bcond(name='zp',type=dirichlet,locator=zp_locator   ,dir='+z')
         ! Assign scalar names and phases
         sc%SCname=[sp_names,'Tl','Tg']; iTl=ns+1; iTg=ns+2
         do isc=1,ns
            pos_open =index(sc%SCname(isc),'(')
            pos_close=index(sc%SCname(isc),')')
            if (pos_open.gt.0.and.pos_close.gt.pos_open) then
               sc%SCname(isc)(pos_open:pos_open)  ='_'
               sc%SCname(isc)(pos_close:pos_close)=' '
               sc%SCname(isc)=adjustl(trim(sc%SCname(isc)))
            end if
            sc%phase(isc)=sys%get_pind(isc)
         end do
         sc%SCname(get_sp_ind('H2O'))='H2O_g'
         sc%phase(iTl)=Lphase
         sc%phase(iTg)=Gphase
         ! Initialize the phasic density and VOF
         sc%Prho(Lphase)=fs%rho_l
         sc%Prho(Gphase)=fs%rho_g
         sc%PVF(:,:,:,Lphase)=vf%VF
         sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
         ! Assign diffusivities
         sc%diff(:,:,:,1:ns)=spDiff
         sc%diff(:,:,:,iTl)=alpha_l
         sc%diff(:,:,:,iTg)=alpha_g
         ! do isc=1,sc%nscalar
         !    p=sc%phase(isc)
         !    if (p.eq.Lphase) then
         !       sc%diff(:,:,:,isc)=alpha_l
         !    else
         !       sc%diff(:,:,:,isc)=alpha_g
         !    end if
         ! end do
         ! Initialize the linear solver
         ss=hypre_str(cfg=cfg,name='Scalar',method=gmres_pfmg2,nst=7)
         call param_read('Scalar iteration',ss%maxit)
         call param_read('Scalar tolerance',ss%rcvg)
         ! Setup the solver
         call sc%setup(implicit_solver=ss)
         ! Initialize mole numbers
         iWv=get_sp_ind('H2O')
         iWl=get_sp_ind('H2O(L)')
         iO2=get_sp_ind('O2')
         iN2=get_sp_ind('N2')
         N_init(iWl)=1.0_WP
         N_init(iO2)=1.0_WP
         N_init(iN2)=N2O_rat*N_init(iO2)
         N_init(iWv)=wv2air_rat*sum(N_init(iO2:iN2))
         ! Get the phase mass
         mp=0.0_WP
         do isc=1,ns
            p=sc%phase(isc)
            mp(p)=mp(p)+N_init(isc)*MM(isc)
         end do
         ! Initialize scalars
         do isc=1,ns
            p=sc%phase(isc)
            where (sc%PVF(:,:,:,p).gt.VFlo)
               sc%SC(:,:,:,isc)=MM(isc)*N_init(isc)/mp(p)
            else where
               sc%SC(:,:,:,isc)=0.0_WP
            end where
         end do
         do i=sc%cfg%imino_,sc%cfg%imaxo_
            do j=sc%cfg%jmino_,sc%cfg%jmaxo_
               do k=sc%cfg%kmino_,sc%cfg%kmaxo_
                  if (vf%VF(i,j,k).gt.VFlo) then
                     sc%SC(i,j,k,iTl)=T_liq
                  end if
                  if (vf%VF(i,j,k).lt.VFhi) then
                     sc%SC(i,j,k,iTg)=T_amb
                  end if
               end do
            end do
         end do
         my_Y=maxval(sc%SC(:,:,:,iO2))
         call MPI_ALLREDUCE(my_Y,YO2_exit,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         my_Y=maxval(sc%SC(:,:,:,iN2))
         call MPI_ALLREDUCE(my_Y,YN2_exit,1,MPI_REAL_WP,MPI_MAX,cfg%comm,ierr)
         Ywv_exit=1.0_WP-YO2_exit-YN2_exit
         ! Apply boundary conditions
         call sc%apply_bcond(timeSC%t,time%dt)
         call sc%get_bcond('xp',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i-1,j,k,iWv)
            sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i-1,j,k,iO2)
            sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i-1,j,k,iN2)
            sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i-1,j,k,iTg)
         end do
         call sc%get_bcond('yp',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i,j-1,k,iWv)
            sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i,j-1,k,iO2)
            sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i,j-1,k,iN2)
            sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i,j-1,k,iTg)
         end do
         call sc%get_bcond('zp',my_bc)
         do n=1,my_bc%itr%no_
            i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
            sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i,j,k-1,iWv)
            sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i,j,k-1,iO2)
            sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i,j,k-1,iN2)
            sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i,j,k-1,iTg)
         end do
         ! Get the phasic face apertures
         call sc%get_face_apt()
         ! One-field temperature
         T=sc%PVF(:,:,:,Lphase)*sc%SC(:,:,:,iTl)+sc%PVF(:,:,:,Gphase)*sc%SC(:,:,:,iTg)
      end block create_scalar


      ! Create and initialize an lgpc object
      create_lgpc: block
         use lgpc_class, only: symmetry
         integer :: i,j,k
         ! Create the object
         call lg%initialize(cfg=cfg,vf=vf,sc=sc%SC,iTl=iTl,iTg=iTg,itp_x=fs%itpr_x,itp_y=fs%itpr_y,itp_z=fs%itpr_z,div_x=fs%divp_x,div_y=fs%divp_y,div_z=fs%divp_z,name='liquid gas pc')
         call param_read('Mass flux tolerence',     lg%mdot3p_tol)
         call param_read('Max pseudo timestep size',lg%pseudo_time%dtmax)
         call param_read('Max pseudo cfl number',   lg%pseudo_time%cflmax)
         call param_read('Max pseudo time steps',   lg%pseudo_time%nmax)
         lg%pseudo_time%dt=lg%pseudo_time%dtmax
         ! Boundary conditions
         call lg%add_bcond(name='xm',type=symmetry,face='x',dir=-1,locator=xm_locator_sc)
         call lg%add_bcond(name='ym',type=symmetry,face='y',dir=-1,locator=ym_locator_sc)
         call lg%add_bcond(name='zm',type=symmetry,face='z',dir=-1,locator=zm_locator_sc)
         ! Get densities from the flow solver
         lg%rho_l=fs%rho_l
         lg%rho_g=fs%rho_g
      end block create_lgpc


      ! Apply the interface jump conditions
      ! call interface_jump()
      ! Get the volumetric lgpc mass flux
      call lg%get_mdot3p()
      ! Initialize the liquid and gas mass fluxes
      call lg%init_mdot3pLG()
      ! Get the interface normal
      call lg%get_normal()


      ! Initialize outputs
      call get_T_Yv()
      call get_T_itf()
      call get_R_drp()


      ! Create surfmesh object for interface polygon output
      create_smesh: block
         smesh=surfmesh(nvar=0,name='plic')
         call vf%update_surfmesh(smesh)
      end block create_smesh


      ! Add Ensight output
      create_ensight: block
         integer :: isc
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='water_drop')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',Ui,Vi,Wi)
         call ens_out%add_scalar('VOF',vf%VF)
         call ens_out%add_scalar('pressure',fs%P)
         call ens_out%add_surface('plic',smesh)
         do isc=1,ns
            call ens_out%add_scalar('Y_'//trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
         end do
         do isc=ns+1,sc%nscalar
            call ens_out%add_scalar(trim(sc%SCname(isc)),sc%SC(:,:,:,isc))
         end do
         call ens_out%add_scalar('mdot3p',lg%mdot3p)
         call ens_out%add_scalar('evp_div',lg%div_vel)
         call ens_out%add_scalar('mdot2p',lg%mdot2p)
         call ens_out%add_scalar('mdot3pL',lg%mdot3pLG(:,:,:,Lphase))
         call ens_out%add_scalar('mdot3pG',lg%mdot3pLG(:,:,:,Gphase))
         call ens_out%add_scalar('divergence',fs%div)
         call ens_out%add_scalar('Temperature',T)
         call ens_out%add_vector('normal',lg%normal(:,:,:,1),lg%normal(:,:,:,2),lg%normal(:,:,:,3))
         ! Debug
         call ens_out%add_scalar('cluster_map',cluster_map)
         call ens_out%add_scalar('PVFL',sc%PVF(:,:,:,Lphase))
         call ens_out%add_scalar('PVFG',sc%PVF(:,:,:,Gphase))
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
         ! Create an event for T and Yv output
         TYv_evt=event(time=time,name='T_Yv output')
         call param_read('T and Yv file output period',TYv_evt%tper)
      end block create_ensight


      ! Create a monitor file
      create_monitor: block
         integer :: isc
         ! Prepare some info about fields
         call fs%get_cfl(time%dt,time%cfl)
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
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
         call mfile%add_column(R_drp,'Droplet radius')
         call mfile%add_column(T_itf,'Interface temperature')
         call mfile%add_column(fs%divmax,'Maximum divergence')
         call mfile%add_column(fs%psolv%it,'Pressure iteration')
         call mfile%add_column(fs%psolv%rerr,'Pressure error')
         call mfile%add_column(prhs_int,'prhs_int')
         ! Debug
         call fs%get_mfr()
         call lg%cfg%integrate(lg%div_vel,mfr_err)
         mfr_err=abs(mfr_err-sum(fs%mfr))
         call mfile%add_column(mfr_err,'mfr_err')
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
         do isc=1,sc%nscalar
            call scfile%add_column(sc%SCmin(isc),trim(sc%SCname(isc))//' min')
            call scfile%add_column(sc%SCmax(isc),trim(sc%SCname(isc))//' max')
            call scfile%add_column(sc%SCint(isc),trim(sc%SCname(isc))//' int')
         end do
         call scfile%write()
         ! Create lgpc monitor
         lgfile=monitor(lg%cfg%amRoot,'lgpc')
         call lgfile%add_column(time%n,'Timestep number')
         call lgfile%add_column(time%t,'Time')
         call lgfile%add_column(lg%pseudo_time%dt,'Pseudo time step')
         call lgfile%add_column(lg%pseudo_time%cfl,'Maximum pseudo CFL')
         call lgfile%add_column(lg%pseudo_time%n,'No. pseudo steps')
         call lgfile%add_column(lg%mdot3p_int,'mdot3p int')
         call lgfile%add_column(lg%mdot3pL_int,'shifted mdot3pL int')
         call lgfile%add_column(lg%mdot3pG_int,'shifted mdot3pG int')
         call lgfile%add_column(lg%mdot3pL_int_err,'mdot3pL int err')
         call lgfile%add_column(lg%mdot3pG_int_err,'mdot3pG int err')
         call lgfile%add_column(lg%mdot3pL_err,'max mdot3pL err')
         call lgfile%add_column(lg%mdot3pG_err,'max mdot3pG err')
         call lgfile%write()
      end block create_monitor

      ! debug
      ! find_cell: block
      !    integer :: ind(3)
      !    ind=cfg%get_ijk_local(pos=[0.000118_WP,0.0_WP,0.000097_WP],ind_guess=[20,20,1])
      !    print*,'ind = ',ind
      !    print*,'x = ',cfg%xm(ind(1)),'y = ',cfg%ym(ind(2)),'z = ',cfg%zm(ind(3))
      !    print*,'cfg%rank = ',cfg%rank
      !    print*,'cfg%imin_ = ',cfg%imin_,'cfg%imax_ = ',cfg%imax_
      !    print*,'cfg%jmin_ = ',cfg%jmin_,'cfg%jmax_ = ',cfg%jmax_
      !    print*,'cfg%kmin_ = ',cfg%kmin_,'cfg%kmax_ = ',cfg%kmax_
      ! end block find_cell


   end subroutine simulation_init


   !> Perform an NGA2 simulation-this mimicks NGA's old time integration for multiphase
   subroutine simulation_run
      use messager, only: die
      use tpns_class, only: static_contact,harmonic_visc
      use mathtools,  only: Pi
      implicit none
      integer  :: i,j,k
      ! Debug
      logical :: flg

      ! Perform time integration
      do while (.not.time%done())

         ! Increment time
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Remember old VOF
         vf%VFold=vf%VF

         ! Remember old lgpc divergence
         lg%div_vel_old=lg%div_vel

         ! Remember old velocity
         fs%Uold=fs%U
         fs%Vold=fs%V
         fs%Wold=fs%W

         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced

         ! Prepare old staggered density (at n)
         call fs%get_olddensity(vf=vf)

         ! if (cfg%amRoot) print*,'Beggining of time step, before advancing VOF: VOF(31,12,1) = ',vf%VF(31,12,1)

         ! VOF solver step
         call vf%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W)
         call vf%apply_bcond(time%t,time%dt)

         ! if (cfg%amRoot) print*,'After advancing VOF: VOF(31,12,1) = ',vf%VF(31,12,1)


         ! debug_vof_div: block
         !    integer, parameter :: i0=26, j0=1, k0=21
         !    real(WP) :: vfold, vfnew, dvf
         !    real(WP) :: net, dtdiv_face
         !    real(WP) :: dt, vol

         !    if (cfg%amRoot) then
         !       vfold = vf%VFold(i0,j0,k0)
         !       vfnew = vf%VF   (i0,j0,k0)

         !       ! Only print when it actually becomes nonzero (your case: ~2.6e-9)
         !       if (vfnew > 1.0e-12_WP .or. vfold > 1.0e-12_WP) then
         !          dt  = time%dt
         !          vol = cfg%vol(i0,j0,k0)
         !          dvf = vfnew - vfold

         !          ! Net "dt * div" computed exactly like the crude face volumes (dt*U*A etc.)
         !          net = (-dt*fs%U(i0+1,j0,k0)*cfg%dy(j0)*cfg%dz(k0) + dt*fs%U(i0,j0,k0)*cfg%dy(j0)*cfg%dz(k0)) &
         !             + (-dt*fs%V(i0,j0+1,k0)*cfg%dz(k0)*cfg%dx(i0) + dt*fs%V(i0,j0,k0)*cfg%dz(k0)*cfg%dx(i0)) &
         !             + (-dt*fs%W(i0,j0,k0+1)*cfg%dx(i0)*cfg%dy(j0) + dt*fs%W(i0,j0,k0)*cfg%dx(i0)*cfg%dy(j0))

         !          dtdiv_face = net / vol

         !          print*, 'VOFDBG step=',time%n,' cell=',i0,j0,k0, &
         !                   ' band=',vf%band(i0,j0,k0), &
         !                   ' VFold=',vfold,' VF=',vfnew,' dVF=',dvf, &
         !                   ' dtDiv_face=',dtdiv_face
         !       end if
         !    end if
         ! end block debug_vof_div


         ! Debug
         ! debug1: block
         !    use irl_fortran_interface, only: getNumberOfVertices
         !    if (cfg%amRoot) print*,'NumberOfVertices for cell 31,12,1',getNumberOfVertices(vf%interface_polygon(1,31,12,1))
         ! end block debug1

         ! ================== SCALAR ================== !

         advance_scalar: block
            use tpscalar_class, only: bcond
            type(bcond), pointer :: my_bc
            integer  :: isc,p,n
            real(WP) :: dt_sc

            ! Increment scalar time step
            if (timeSC%t+timeSC%dt.gt.time%t) then
               dt_sc=timeSC%dt
               timeSC%dt=time%t-timeSC%t
               call timeSC%increment()
               timeSC%dt=dt_sc
            else
               call timeSC%increment()
            end if

            flg=.false.

            ! Remember old SC
            sc%SCold =sc%SC
            sc%PVFold=sc%PVF

            ! Update the phasic VOF and face apertures
            sc%PVF(:,:,:,Lphase)=vf%VF
            sc%PVF(:,:,:,Gphase)=1.0_WP-vf%VF
            call sc%get_face_apt()

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    ! print*,'PVF(9 :11)',sc%PVF(9:11,27,15,Gphase)
            !    ! print*,'face_apt_x(10:11) = ',sc%face_apt_x(10:11,27,15,Gphase)
            !    ! print*,'PVF(26:28)',sc%PVF(10,26:28,15,Gphase)
            !    ! print*,'face_apt_y(27:28) = ',sc%face_apt_y(10,27:28,15,Gphase)
            !    ! print*,'PVF(14:16)',sc%PVF(10,27,14:16,Gphase)
            !    ! print*,'face_apt_z(15:16) = ',sc%face_apt_z(10,27,15:16,Gphase)
            !    print*,'before extrapolation Tl = ',sc%SC(31,12,1,iTl)
            !    print*,'VOF = ',vf%VF(31,12,1)
            ! end if

            ! Correct for the emptied out and new interfacial cells
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               do k=cfg%kmino_,cfg%kmaxo_
                  do j=cfg%jmino_,cfg%jmaxo_
                     do i=cfg%imino_,cfg%imaxo_
                        if (sc%PVF(i,j,k,p).eq.0.0_WP) sc%SC(i,j,k,isc)=0.0_WP
                     end do
                  end do
               end do
               call lg%pure_zero_interfacial_extp(p,sc%SC(:,:,:,isc))
            end do

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'after extrapolation Tl = ',sc%SC(31,12,1,iTl)
            ! end if

            ! Explicit calculation of dVOFSC/dt from scalar advection
            call sc%get_dSCdt_adv(dSCdt=resSC,U=fs%U,V=fs%V,W=fs%W,detailed_face_flux=vf%detailed_face_flux,dt=timeSC%dt)

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'advection res Tl = ',resSC(31,12,1,iTl)
            ! end if

            ! Advance scalar advection
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               where (sc%mask.eq.0.and.sc%PVF(:,:,:,p).gt.0.0_WP) sc%SC(:,:,:,isc)=(sc%PVFold(:,:,:,p)*sc%SCold(:,:,:,isc)+timeSC%dt*(resSC(:,:,:,isc)+lg%div_vel_old(:,:,:)*sc%SCold(:,:,:,isc)))/sc%PVF(:,:,:,p)
               where (sc%PVF(:,:,:,p).eq.0.0_WP) sc%SC(:,:,:,isc)=0.0_WP
            end do

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'VOFold*TlOld/VOF = ',sc%PVFold(31,12,1,Lphase)*sc%SCold(31,12,1,iTl)/sc%PVF(31,12,1,Lphase)
            !    print*,'dt*div(VOF*u*Tl)/VOF = ',timeSC%dt*resSC(31,12,1,iTl)/sc%PVF(31,12,1,Lphase)
            !    print*,'dt*div(u)*TlOld/VOF = ',timeSC%dt*lg%div_vel_old(31,12,1)*sc%SCold(31,12,1,iTl)/sc%PVF(31,12,1,Lphase)
            !    print*,'after advection Tl = ',sc%SC(31,12,1,iTl)
            ! end if

            ! Explicit calculation of dVOFSC/dt from scalar diffusion
            call sc%get_dSCdt_dff(dSCdt=resSC)
            do isc=1,sc%nscalar
               p=sc%phase(isc)
               where (sc%mask.eq.0.and.sc%PVF(:,:,:,p).gt.0.0_WP) resSC(:,:,:,isc)=timeSC%dt*resSC(:,:,:,isc)/sc%PVF(:,:,:,p)
               ! where (sc%PVF(:,:,:,p).eq.0.0_WP) resSC(:,:,:,isc)=0.0_WP
            end do
            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'rhs Tl for linear solver = ',resSC(31,12,1,iTl)
            ! end if

            ! Form implicit diffusive residual
            call sc%solve_implicit_dff(timeSC%dt,resSC)

            ! Advance scalar diffusion
            sc%SC=sc%SC+resSC

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'after diffusion Tl = ',sc%SC(31,12,1,iTl)
            ! end if

            !
            where (vf%VF.gt.0.0_WP) sc%SC(:,:,:,iWl)=1.0_WP

            ! do k=cfg%kmin_,cfg%kmax_
            !    do j=cfg%jmin_,cfg%jmax_
            !       do i=cfg%imin_,cfg%imax_
            !          if (sc%SC(i,j,k,iTl).gt.353.01_WP) then
            !             ! print*,'-------------------------'
            !             ! print*,'Tl = ',sc%SC(i,j,k,iTl)
            !             ! print*,'i,j,k = ',i,j,k
            !             ! print*,'VOF = ',vf%VF(i,j,k)
            !             ! print*,'-------------------------'
            !             flg=.true.
            !          end if
            !       end do
            !    end do
            ! end do
            ! if (flg) call die('')

            ! Apply boundary conditions
            call sc%apply_bcond(timeSC%t,timeSC%dt)
            call sc%get_bcond('xp',my_bc)
            do n=1,my_bc%itr%no_
               i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
               sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i-1,j,k,iWv)
               sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i-1,j,k,iO2)
               sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i-1,j,k,iN2)
               sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i-1,j,k,iTg)
            end do
            call sc%get_bcond('yp',my_bc)
            do n=1,my_bc%itr%no_
               i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
               sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i,j-1,k,iWv)
               sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i,j-1,k,iO2)
               sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i,j-1,k,iN2)
               sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i,j-1,k,iTg)
            end do
            call sc%get_bcond('zp',my_bc)
            do n=1,my_bc%itr%no_
               i=my_bc%itr%map(1,n); j=my_bc%itr%map(2,n); k=my_bc%itr%map(3,n)
               sc%SC(i,j,k,iWv)=2.0_WP*Ywv_exit-sc%SC(i,j,k-1,iWv)
               sc%SC(i,j,k,iO2)=2.0_WP*YO2_exit-sc%SC(i,j,k-1,iO2)
               sc%SC(i,j,k,iN2)=2.0_WP*YN2_exit-sc%SC(i,j,k-1,iN2)
               sc%SC(i,j,k,iTg)=2.0_WP*T_amb   -sc%SC(i,j,k-1,iTg)
            end do

            ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
            !    print*,'after bc Tl = ',sc%SC(31,12,1,iTl)
            ! end if

         end block advance_scalar

         ! ================== PHASE CHANGE ================== !

         ! Apply the interface jump conditions
         call interface_jump()
         ! if (cfg%amRoot) print*,'After jumping'
         ! if (16.ge.cfg%imin_.and.16.le.cfg%imax_.and.27.ge.cfg%jmin_.and.27.le.cfg%jmax_.and.12.ge.cfg%kmin_.and.12.le.cfg%kmax_) print*,'VOF = ',vf%VF(31,12,1)
         ! if (cfg%amRoot) print*,'*******'

         ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
         !    print*,'Tl = ',sc%SC(31,12,1,iTl)
         !    ! print*,'x,y,z = ',cfg%xm(16),cfg%ym(27),cfg%zm(12)
         ! end if

         ! Get the volumetric lgpc mass flux
         call lg%get_mdot3p()

         ! Shift the lgpc mass flux
         call lg%shift_mdot3p()

         ! Get the phase-change induced divergence
         call lg%get_div()

         ! if (cfg%iproc.eq.1.and.cfg%jproc.eq.1.and.cfg%kproc.eq.1) then
         !    print*,'after lgpc shift Tl = ',sc%SC(31,12,1,iTl)
         ! end if

         ! ================== VELOCITY ================== !

         advance_flow: block

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

               ! Add momentum mass fluxes
               call fs%addsrc_gravity(resU,resV,resW)

               ! Assemble explicit residual
               resU=-2.0_WP*fs%rho_U*fs%U+(fs%rho_Uold+fs%rho_U)*fs%Uold+time%dt*resU
               resV=-2.0_WP*fs%rho_V*fs%V+(fs%rho_Vold+fs%rho_V)*fs%Vold+time%dt*resV
               resW=-2.0_WP*fs%rho_W*fs%W+(fs%rho_Wold+fs%rho_W)*fs%Wold+time%dt*resW

               ! Form implicit residuals
               call fs%solve_implicit(time%dt,resU,resV,resW)

               ! Apply these residuals
               fs%U=2.0_WP*fs%U-fs%Uold+resU!/fs%rho_U
               fs%V=2.0_WP*fs%V-fs%Vold+resV!/fs%rho_V
               fs%W=2.0_WP*fs%W-fs%Wold+resW!/fs%rho_W

               ! Apply boundary conditions
               call fs%apply_bcond(time%t,time%dt)
               call apply_dirichlet()

               ! Solve Poisson equation
               call fs%update_laplacian()
               call fs%correct_mfr(src=lg%div_vel)
               call fs%get_div(src=lg%div_vel)
               ! call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf,contact_model=static_contact)
               call fs%add_surface_tension_jump(dt=time%dt,div=fs%div,vf=vf)
               fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dt
               fs%psolv%sol=0.0_WP
               call fs%psolv%solve()
               call fs%shift_p(fs%psolv%sol)

               ! Correct velocity
               call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
               call cfg%integrate(fs%psolv%rhs,prhs_int)
               fs%P=fs%P+fs%psolv%sol
               fs%U=fs%U-time%dt*resU/fs%rho_U
               fs%V=fs%V-time%dt*resV/fs%rho_V
               fs%W=fs%W-time%dt*resW/fs%rho_W

               ! Increment sub-iteration counter
               time%it=time%it+1

            end do

            ! Recompute interpolated velocity and divergence
            call fs%interp_vel(Ui,Vi,Wi)
            call fs%get_div(src=lg%div_vel)

            ! Debug
            call fs%get_mfr()
            call lg%cfg%integrate(lg%div_vel,mfr_err)
            mfr_err=abs(mfr_err-sum(fs%mfr))

         end block advance_flow

         ! Output to ensight
         T=vf%VF*sc%SC(:,:,:,iTl)+(1.0_WP-vf%VF)*sc%SC(:,:,:,iTg)
         if (ens_evt%occurs()) then
            call vf%update_surfmesh(smesh)
            call ens_out%write_data(time%t)
         end if

         ! Get T and Yv profiles
         if (TYv_evt%occurs()) call get_T_Yv()

         ! Get interface temperature
         call get_T_itf()

         ! Get droplet radius
         call get_R_drp()

         ! Perform and output monitoring
         call fs%get_max()
         call vf%get_max()
         call sc%get_max()
         call mfile%write()
         call cflfile%write()
         call scfile%write()
         call lgfile%write()

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
      deallocate(resU,resV,resW,Ui,Vi,Wi,resSC,T)

   end subroutine simulation_final


end module simulation
