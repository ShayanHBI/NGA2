!> Finite chem class
!> Extends multivdscalar class for finite rate chemistry calculations
module finitechem_class
   use precision,           only: WP
   use config_class,        only: config
   use multivdscalar_class, only: multivdscalar
   use fcmech
   use dvode_f90_m
   implicit none
   private

   ! Expose type/constructor/methods
   public :: finitechem

   ! Clipping values for temperature
   real(WP), parameter :: T_min=50.0_WP
   real(WP), parameter :: T_max=5000.0_WP

   ! Solver options
   logical :: use_jacanal=.false.
   logical :: use_lewis=.false.

   !> Finite chemistry solver object definition
   type, extends(multivdscalar) :: finitechem
      
      ! Thermodynamic pressure and density
      real(WP) :: Pthermo,Pthermo_old
      real(WP) :: RHOmean,RHO_0
      
      ! Scalar variables
      real(WP), dimension(:,:,:),   pointer     :: T                          !< Temperature
      real(WP), dimension(:,:,:,:), pointer     :: Y                          !< Composition
      real(WP), dimension(:,:,:),   allocatable :: visc                       !< Viscosity field
      real(WP), dimension(:,:,:,:), allocatable :: SRCchem                    !< Chemical source terms
      real(WP), dimension(:,:,:,:), allocatable :: SRC                        !< Total source terms for scalar equations
      real(WP), dimension(:,:,:),   allocatable :: h                          !< Enthalpy
      real(WP), dimension(:,:,:),   allocatable :: lambda                     !< Thermal conductivity
      real(WP), dimension(:,:,:),   allocatable :: Cp                         !< Heat capacity
      
      ! Metrics
      real(WP), dimension(:,:,:,:), allocatable :: grdsc_xm,grdsc_ym,grdsc_zm !< Scalar gradient for SC at centers
      logical :: use_explicit_try=.false.
      
      ! Monitoring quantities
      real(WP) :: visc_min,visc_max,diff_min,diff_max                         !< Maximum and minimum
   
      !> Scheduler variables
      logical :: use_scheduler
      ! Identities
      integer  :: imaster,inewmaster
      integer  :: nbundles,nbundles_max,nproc_waiting
      real(WP) :: bundleref
      integer, dimension(:), allocatable :: imaster_list,imaster_aware,iproc_waiting
      ! Buffers
      integer :: nwhere_buf
      real(WP), dimension(:), allocatable :: bufferR,bufferS
      ! Locations
      integer, dimension(:),   allocatable :: nwhere
      integer, dimension(:),   allocatable :: iwhere_buf,jwhere_buf,kwhere_buf
      integer, dimension(:,:), allocatable :: iwhere,jwhere,kwhere

   contains

      procedure :: scheduler_init
      procedure :: clip
      procedure :: react
      procedure :: get_density
      procedure :: get_viscosity
      procedure :: get_diffusivity
      procedure :: get_cpmix
      procedure :: get_enthalpy
      procedure :: get_Wmix
      procedure :: diffusive_source
      procedure :: pressure_source
      procedure :: update_pressure
      procedure :: get_src
      procedure :: get_max=>fc_get_max

   end type finitechem

   !> Declare finitechem constructor
   interface finitechem
      procedure constructor
   end interface finitechem

contains


   !> Default constructor for finitchem object
   function constructor(cfg,scheme,name) result(self)
      implicit none
      type(finitechem) :: self
      class(config),target,intent(in) :: cfg
      integer,intent(in) :: scheme
      character(len=*),optional :: name
      character(len=str_medium), dimension(nspec) :: names
      integer :: i,j,k

      ! Create a multivdscalar object
      self%multivdscalar=multivdscalar(cfg=cfg,scheme=scheme,nscalar=nspec+1,name=name)

      ! Get species names
      call fcmech_get_speciesnames(names)

      ! Scalar names
      do i=1,nspec
         self%SCname(i)=names(i)
      end do
      self%SCname(nspec+1)='T'

      ! Allocate variables
      allocate(self%visc(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));            self%visc=0.0_WP
      allocate(self%lambda(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));          self%lambda=0.0_WP
      allocate(self%Cp(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));              self%Cp=0.0_WP
      allocate(self%h(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));               self%h=0.0_WP
      allocate(self%SRCchem(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,nspec+1)); self%SRCchem=0.0_WP
      allocate(self%SRC(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,nspec+1));     self%SRC=0.0_WP

      ! Allocate finite difference gradient operators
      allocate (self%grdsc_xm(0:+1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
      allocate (self%grdsc_ym(0:+1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))
      allocate (self%grdsc_zm(0:+1,self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_))

      ! Create gradient coefficients to cell faces
      do k=self%cfg%kmin_,self%cfg%kmax_+1
         do j=self%cfg%jmin_,self%cfg%jmax_+1
            do i=self%cfg%imin_,self%cfg%imax_+1
               self%grdsc_xm(:,i,j,k)=self%cfg%dxi(i)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in x from [xm,ym,zm] to [x,ym,zm]
               self%grdsc_ym(:,i,j,k)=self%cfg%dyi(j)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in y from [xm,ym,zm] to [xm,y,zm]
               self%grdsc_zm(:,i,j,k)=self%cfg%dzi(k)*[-1.0_WP,+1.0_WP] !< FD gradient of SC in z from [xm,ym,zm] to [xm,ym,z]
            end do
         end do
      end do

      ! Pointers to composition and temperature
      self%Y=>self%SC(:,:,:,1:nspec)
      self%T=>self%SC(:,:,:,nspec+1)

   end function constructor


   !> Initializer for the scheduler
   subroutine scheduler_init(this,bundleref)
      use mpi_f08, only: MPI_MAX,MPI_INTEGER
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP),intent(in) :: bundleref
      integer :: ierr
      this%bundleref=bundleref
      this%nbundles=int((this%cfg%imax_-this%cfg%imin_+1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)/this%bundleref)
      call MPI_ALLREDUCE(this%nbundles,this%nbundles_max,1,MPI_INTEGER,MPI_MAX,this%cfg%comm,ierr)
      print*,'Number of cells per bundle:',this%nbundles,'[proc: ',this%cfg%rank,']'
      ! allocate(this%imaster_list(0:this%cfg%nproc-1),this%imaster_aware(0:this%cfg%nproc-1))
      allocate(this%imaster_list(this%cfg%nproc),this%imaster_aware(this%cfg%nproc))
      allocate(this%bufferS(this%nbundles_max*(nspec+1)),this%bufferR(this%nbundles_max*(nspec+1)))
      allocate(this%iwhere_buf(this%nbundles_max),this%jwhere_buf(this%nbundles_max),this%kwhere_buf(this%nbundles_max))
      ! allocate(this%iwhere(0:this%cfg%nproc-1,this%nbundles_max),this%jwhere(0:this%cfg%nproc-1,this%nbundles_max),this%kwhere(0:this%cfg%nproc-1,this%nbundles_max))
      ! allocate(this%nwhere(0:this%cfg%nproc-1))
      ! allocate(this%iproc_waiting(0:this%cfg%nproc-1))
      allocate(this%iwhere(this%cfg%nproc,this%nbundles_max),this%jwhere(this%cfg%nproc,this%nbundles_max),this%kwhere(this%cfg%nproc,this%nbundles_max))
      allocate(this%nwhere(this%cfg%nproc))
      allocate(this%iproc_waiting(this%cfg%nproc))

   end subroutine scheduler_init


   !> Clip and rescale mass fractions
   subroutine clip(this,myY)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP),intent(inout), dimension(nspec) :: myY
      myY=min(max(myY,0.0_WP),1.0_WP)
      myY=myY/sum(myY)
   end subroutine clip


   !> Calculate reaction source terms
   subroutine react(this,dt)
      use parallel, only: MPI_REAL_WP
      use mpi_f08, only: MPI_TAG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_SOURCE,MPI_STATUS_SIZE,MPI_INTEGER,MPI_LOGICAL
      use messager, only: die
      use,intrinsic :: ISO_C_BINDING, only: c_loc
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP),intent(in) :: dt  !< Timestep size over which to advance
      integer :: nsc,i,j,k,myi,myj,myk
      real(WP), dimension(nspec+1) :: sol,solold

      ! Local scheduler variables
      integer :: ndata,ibuf,istatus
      ! Tags
      integer :: itag_ndata,itag_data,itag_ihead,itag_idle,itag_done,itag_imaster
      integer :: idata,icount_,ncount_
      logical :: ldone,all_cells
      ! Misc
      integer :: islave,itag,ip
      ! MPI
      integer, dimension(MPI_STATUS_SIZE) :: status
      integer :: ierr,iexit_master
      ! List of particles to send out for direct integration
      integer :: nDI
      integer, dimension(:), pointer :: iDI
      ! Temp variables to update particle properties
      real(WP) :: rbuf,myw,mysv,myh
      CHARACTER(LEN=30) :: Format

      ! Initialize with zeros
      this%SRCchem=0.0_WP
      ! If only one processor,or if beginning of simulation,just do the work for all particles
      if (.not.this%use_scheduler.or.this%cfg%nproc.eq.1) then

         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  sol(1:nspec)=min(max(this%SC(i,j,k,1:nspec),0.0_WP),1.0_WP)
                  sol(1:nspec)=sol(1:nspec)/sum(sol(1:nspec))
                  sol(nspec+1)=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
                  solold=sol
                  if (solold(sN2).gt.0.8_WP) cycle                        ! Package initial solution vector
                  ! call this%clip(sol(1:nspec))
                  ! Advance the chemical equations for each particle
                  call fc_reaction_compute_sol(sol,this%pthermo,dt)
                  this%SRCchem(i,j,k,:)=sol-solold
               end do
            end do
         end do

         ! Sync
         do nsc=1,nspec+1
            call this%cfg%sync(this%SRCchem(:,:,:,nsc))
         end do

         ! Stop there,no dynamic scheduling in this case
         return

      end if

      ! ! ------------------------------------------- !
      ! ! Dynamic scheduler here for load balancing
      ! ! ------------------------------------------- !

      ! ! Tag meanings:
      ! ! TAG=1: sending number of data
      ! itag_ndata=1
      ! ! TAG=2: sending data
      ! itag_data=2
      ! ! TAG=3: identity of master
      ! itag_imaster=3
      ! ! TAG=4: idle
      ! itag_idle=4
      ! ! TAG=5: quit signal
      ! itag_done=5

      ! ! Initialize master/slave identities
      ! ! Initial master id
      ! this%imaster=1
      ! ! flag that a new master has been promoted
      ! this%inewmaster=-1
      ! ! Processors that have been master already
      ! this%imaster_list=0
      ! this%imaster_list(this%imaster)=1
      ! ! Flag indicating that data need to be sent back to master
      ! idata=0
      ! ! Extra integer buffer
      ! ibuf=1
      ! ! Not done to start with
      ! ldone=.false.

      ! ! Loop until all work is done
      ! scheduler_loop: do while (.not.ldone)


      !    ! Master loop
      !    if (this%cfg%rank+1.eq.this%imaster) then

      !       ! Initializing processor roles and buffers
      !       this%nproc_waiting=0
      !       this%iproc_waiting=0
      !       this%bufferR=0.0_WP
      !       this%bufferS=0.0_WP
      !       ! Initialize ndata
      !       ndata=0
      !       ! Initialize loop condition
      !       iexit_master=0
      !       !   nretrieves=0
      !       ! Set up isat mode for master
      !       !   call fcsubs_retrieve_start
      !       ! Initialize counter
      !       icount_=1
      !       ncount_=(this%cfg%imax_-this%cfg%imin_+1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)

      !       ! Let's send out all my data
      !       master_loop: do while (iexit_master.eq.0) ! No explicit exit conditions here,automatically handled below

      !          ! ------------------------------------------- !
      !          ! Gather compositions until there is a bundle of ndata particles ready to send
      !          ! Skip locally processed compositions
      !          do while (ndata.lt.this%nbundles.and.icount_.le.ncount_)

      !             ! Figure out which composition to consider next to get more data
      !             i=int(icount_/((this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)))+1
      !             j=int((icount_-(i-1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1))/(this%cfg%kmax_-this%cfg%kmin_+1))+1
      !             k=icount_-(i-1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)-(j-1)*(this%cfg%kmax_-this%cfg%kmin_+1)+1
      !             i=i+this%cfg%imin_-1
      !             j=j+this%cfg%jmin_-1
      !             k=k+this%cfg%kmin_-1

      !             ! Skip walls
      !             if (this%mask(i,j,k).eq.0) then

      !                ! Initialize solution vector
      !                sol(1:nspec)=min(max(this%SC(i,j,k,1:nspec),0.0_WP),1.0_WP)
      !                sol(1:nspec)=sol(1:nspec)/sum(sol(1:nspec))
      !                sol(nspec+1)=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
      !                ! Try PLP with ISAT to get solution at t+deltat
      !                !   call fcsubs_retrieve(sol,Pthermo,dt,istatus)
      !                istatus=0
      !                ! Check status
      !                if (istatus.eq.1) then
      !                   ! If retrieve was successful,transfer solution back!
      !                   !   this%SC(i,j,k,isc_1:isc_1+nspec-1)=sol(1:nspec)
      !                   !   this%SC(i,j,k,isc_1:isc_T)=sol(nspec+1)
      !                   !   nretrieves=nretrieves+1
      !                else
      !                   ! Buffer composition to send out to other processors
      !                   ndata=ndata+1
      !                   this%bufferS((ndata-1)*(nspec+1)+1:ndata*(nspec+1))=sol
      !                   this%iwhere_buf(ndata)=icount_ ! link between bundle and pmc array
      !                   this%nwhere_buf=ndata ! Keep track of actual number of compos in bundle
      !                end if

      !             end if

      !             ! Increment index of next composition to consider
      !             icount_=icount_+1

      !          end do

      !          ! ------------------------------------------- !
      !          ! Listen to the slaves until one raises its hand
      !          call MPI_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,this%cfg%comm,status,ierr)
      !          ! Got something: which slave is talking?
      !          islave=status(MPI_SOURCE)+1
      !          ! What message is it sending?
      !          itag=status(MPI_TAG)

      !          ! ------------------------------------------- !
      !          ! Does this slave have results to send?
      !          if (itag.eq.itag_data) then
      !             ! itag=1: slave has data to send back,receive them!
      !             call MPI_recv(this%bufferR(1:this%nwhere(islave)*(nspec+1)),this%nwhere(islave)*(nspec+1),MPI_REAL_WP,islave-1,itag_data,this%cfg%comm,status,ierr)
      !             ! Store each bufferR compo at correct location in pmc array
      !             do i=1,this%nwhere(islave)
      !                ! Figure out where to store the result
      !                myi=int(this%iwhere(islave,i)/((this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)))+1
      !                myj=int((this%iwhere(islave,i)-(myi-1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1))/(this%cfg%kmax_-this%cfg%kmin_+1))+1
      !                myk=this%iwhere(islave,i)-(myi-1)*(this%cfg%jmax_-this%cfg%jmin_+1)*(this%cfg%kmax_-this%cfg%kmin_+1)-(myj-1)*(this%cfg%kmax_-this%cfg%kmin_+1)+1
      !                myi=myi+this%cfg%imin_-1
      !                myj=myj+this%cfg%jmin_-1
      !                myk=myk+this%cfg%kmin_-1
      !                ! this%SC(myi,myj,myk,1:nspec)=this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1)-1)
      !                ! this%SC(myi,myj,myk,nspec+1)=this%bufferR(i*(nspec+1))
      !                this%SRCchem(myi,myj,myk,1:nspec)=this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1)-1)
      !                this%SRCchem(myi,myj,myk,nspec+1)=this%bufferR(i*(nspec+1))

      !             end do
      !          else
      !             ! Just receive its empty message
      !             call MPI_recv(ibuf,1,MPI_INTEGER,islave-1,itag_idle,this%cfg%comm,status,ierr)
      !          end if

      !          ! ------------------------------------------- !
      !          ! Do I have more work to do?
      !          if (ndata.gt.0) then ! ndata not 0: last bundle has not been sent yet,more work to do!
      !             ! Keep a note of what was sent to islave
      !             this%iwhere(islave,:)=this%iwhere_buf
      !             this%nwhere(islave)=this%nwhere_buf
      !             ! If yes,send the number of data that will be sent
      !             call MPI_send(this%nwhere(islave),1,MPI_INTEGER,islave-1,itag_ndata,this%cfg%comm,ierr)
      !             ! Send the next chunk to this slave
      !             call MPI_send(this%bufferS(1:this%nwhere(islave)*(nspec+1)),ndata*(nspec+1),MPI_REAL_WP,islave-1,itag_data,this%cfg%comm,ierr)
      !             ! Reset buffer to start accumulating more composition for next call
      !             ndata=0
      !             this%iwhere_buf=0
      !             this%nwhere_buf=0
      !             ! Go back to beginning of loop
      !             cycle master_loop
      !          end if

      !          ! ------------------------------------------- !
      !          ! If reaching here,no more work to do,
      !          ! ready to finish work as master and switch to slave status

      !          ! ------------------------------------------- !
      !          ! Has a new master being promoted?
      !          if (this%inewmaster.eq.-1) then !no,no new master yet

      !             ! Reset this%imaster_aware to 0
      !             this%imaster_aware=0
      !             ! But I am aware already,that counts for someting
      !             this%imaster_aware(this%cfg%rank+1)=1

      !             ! Has this slave been a master yet?
      !             if (this%imaster_list(islave).eq.0) then ! no,not yet
      !                ! Promote the slave to master status
      !                this%inewmaster=islave
      !                ! Tell the current slave that it is the new master and update aware list
      !                call MPI_send(this%inewmaster,1,MPI_INTEGER,islave-1,itag_imaster,this%cfg%comm,ierr)
      !                this%imaster_aware(islave)=1
      !                ! Send id of new master to all waiting slaves and update aware list
      !                do ip=1,this%nproc_waiting
      !                   call MPI_send(this%inewmaster,1,MPI_INTEGER,this%iproc_waiting(ip)-1,itag_imaster,this%cfg%comm,ierr)
      !                   this%imaster_aware(this%iproc_waiting(ip))=1
      !                end do

      !             else ! Yes,this proc has been master already
      !                ! Add identity of this slave to the waiting list
      !                this%nproc_waiting=this%nproc_waiting+1
      !                this%iproc_waiting(this%nproc_waiting)=islave
      !                ! Check if we are fully done (all procs have been masters already)
      !                if (this%nproc_waiting.eq.this%cfg%nproc-1) then
      !                   exit master_loop
      !                end if
      !             end if

      !             ! ------------------------------------------- !
      !          else ! Yes,a new master has been promoted
      !             ! Send the id of the new master to the current slave and update aware list
      !             call MPI_send(this%inewmaster,1,MPI_INTEGER,islave-1,itag_imaster,this%cfg%comm,ierr)
      !             this%imaster_aware(islave)=1
      !          end if

      !          ! ------------------------------------------- !
      !          ! Continue waiting for signals till all slaves have received new master id
      !          if (sum(this%imaster_aware).eq.this%cfg%nproc) iexit_master=1

      !       end do master_loop

      !       ! ------------------------------------------- !
      !       ! End isat mode for master
      !       !   call fcsubs_retrieve_stop

      !       ! Print number of retrieves done by the master
      !       !print*,'master',this%cfg%rank,': ',nretrieves,' out of ',npmc_,'[',real(nretrieves,WP)/real(npmc_),']'

      !       ! ------------------------------------------- !
      !       ! Check if I am the last master
      !       if (sum(this%imaster_list).eq.this%cfg%nproc) then
      !          ! If so,update done
      !          ldone=.true.
      !          ! Send a quit signal to every body
      !          do ip=1,this%cfg%nproc
      !             if (ip.eq.this%cfg%rank+1) cycle
      !             call MPI_send(ldone,1,MPI_LOGICAL,ip-1,itag_done,this%cfg%comm,ierr)
      !          end do
      !       else ! I am not the last master
      !          ! Switching the master id to somebody else
      !          this%imaster=this%inewmaster
      !          ! Keep it rolling
      !          cycle scheduler_loop
      !       end if

      !       ! ------------------------------------------- !
      !       ! ------------------------------------------- !
      !       ! Slave loop
      !    else

      !       slave_loop: do while (.true.) ! No explicit exit conditions here,automatically handled below
      !          ! ------------------------------------------- !
      !          ! Do I have results to send to the master?
      !          if (idata.eq.1) then ! yes,I have data to send
      !             ! Send them to master
      !             call MPI_send(this%bufferR(1:ndata*(nspec+1)),ndata*(nspec+1),MPI_REAL_WP,this%imaster-1,itag_data,this%cfg%comm,ierr)
      !             ! No more data to send for now
      !             idata=0
      !          else ! No,no data to send
      !             ! Just tell the master I am available
      !             call MPI_send(ibuf,1,MPI_INTEGER,this%imaster-1,itag_idle,this%cfg%comm,ierr)
      !          end if

      !          ! ------------------------------------------- !
      !          ! Listening to the master to see what is coming next
      !          call MPI_probe(this%imaster-1,MPI_ANY_TAG,this%cfg%comm,status,ierr)
      !          ! What message is it sending?
      !          itag=status(MPI_TAG)

      !          ! ------------------------------------------- !
      !          ! What is the message?
      !          if (itag.eq.itag_done) then ! Quit message
      !             ! Receive it (we just tested the tag here,still need to actually receive the integer)
      !             call MPI_recv(ldone,1,MPI_LOGICAL,this%imaster-1,itag_done,this%cfg%comm,status,ierr)
      !             cycle scheduler_loop

      !             ! ------------------------------------------- !
      !          elseif (itag.eq.itag_ndata) then ! Work message
      !             ! Receiving number of data expected
      !             call MPI_recv(ndata,1,MPI_integer,this%imaster-1,itag_ndata,this%cfg%comm,status,ierr)
      !             ! Receiving chunk of data to process
      !             call MPI_recv(this%bufferR(1:ndata*(nspec+1)),ndata*(nspec+1),MPI_REAL_WP,this%imaster-1,itag_data,this%cfg%comm,status,ierr)
      !             ! Do the work
      !             do i=1,ndata
      !                solold=this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1))
      !                call fc_reaction_compute_sol(this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1)),this%Pthermo,dt)
      !                this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1))=this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1))-solold
      !             end do
      !             ! Now I have data to send
      !             idata=1
      !             ! Cycle the slave loop again
      !             cycle slave_loop

      !             ! ------------------------------------------- !
      !          elseif (itag.eq.itag_imaster) then ! Master has changed
      !             ! Receive the id of the new master
      !             call MPI_recv(ibuf,1,MPI_integer,this%imaster-1,itag_imaster,this%cfg%comm,status,ierr)
      !             ! Update the id of master
      !             this%imaster=ibuf
      !             this%imaster_list(this%imaster)=1
      !             ! Cycle the scheduler loop
      !             cycle scheduler_loop

      !          end if
      !       end do slave_loop
      !    end if
      ! end do scheduler_loop

      ! Sync
      do nsc=1,nspec+1
         call this%cfg%sync(this%SRCchem(:,:,:,nsc))
      end do

   contains

      ! Compute solution after time step delta t
      subroutine fc_reaction_compute_sol(sol,mypressure,dt)
         use random
         use messager, only: die
         implicit none
         real(WP), dimension(nspec+1) :: sol,dsol,solcheck
         real(WP) :: dt,t_chem
         real(WP) :: mypressure
         real(WP) :: tstop,tstart
         type(vode_opts), save :: opts
         integer  :: istate,itask
         real(WP) :: atol=1.0e-15_WP
         real(WP) :: rtol=1.0e-12_WP
         integer :: n
         ! external :: fc_reaction_compute_rhs
         ! external :: fc_reaction_compute_jac
         ! external :: pdfrmap1
         real(WP), dimension(22) :: rstats
         integer, dimension(31) :: istats
         ! ISAT parameters
         real(WP) :: xx(nspec+2),f(nspec+1),dfdx(nspec+1,nspec+2),hvar(1),stats(100)
         integer :: iusr(1)
         ! Save solution
         solcheck=sol
         ! Set values for DVODE
         itask=1; istate=1; tstart=0.0_WP; tstop=dt
         ! Direct integration
         call fc_reaction_compute_rhs(nspec+1,0.0_WP,solcheck,dsol)
         t_chem=sol(nspec+1)/(abs(dsol(nspec+1))+epsilon(1.0_WP))
         if (dt.lt.0.0001_WP*t_chem.and.this%use_explicit_try) then
            sol=sol+dsol*dt
         else
            ! Integrate using DVODE
            if (use_jacanal) then
               continue
               ! opts=set_opts(dense_j=.true.,mxstep=500000,abserr=atol,relerr=rtol,tcrit=tstop,user_supplied_jacobian=.true.)
               !  call dvode_f90(fc_reaction_compute_rhs,nspec+1,sol,tstart,tstop,itask,istate,opts,j_fcn=fc_reaction_compute_jac)
            else
               opts=set_opts(method_flag=22,mxstep=500000,abserr=atol,relerr=rtol,tcrit=tstop)
               call dvode_f90(fc_reaction_compute_rhs,nspec+1,sol,tstart,tstop,itask,istate,opts)
            end if
            call get_stats(rstats,istats)
            call release_opts_arrays(opts)
            ! Error handling
            if (istate.ne.2) then
               print *,'NCF: ',istats(21) ! No. of convergence failures of the nonlinear solver so far.
               print *,'NEF: ',istats(22) ! No. of error test failures of the integrator so far.
               print *,'tstop: ',tstop,' tstart:',tstart,' itask: ',itask,' istate: ',istate
               print *,'solold---------'
               do n=1,nspec
                  print *,'sol(',n,')=',solcheck(n),'_WP !'
               end do
               print *,'sol(NT)=',solcheck(nspec+1),'_WP !','T'
               print *,'-------------------------------'
               print *,'sol',sol
               print *,'-------------------------------'
               do n=1,nspec
                  print *,'dsol(',n,')=',dsol(n),'_WP !'
               end do
               print *,'T','dsol(n)',dsol(nspec+1)
               call die('fc_reaction_source: Direct integration-DVODE failed to converge.')
            end if
         end if
         ! Clip and renormalize,just in case
         call this%clip(sol(1:nspec))
      end subroutine fc_reaction_compute_sol

      ! Computes the chemical source term of the system (called by solver)
      subroutine fc_reaction_compute_rhs(n_,t_,sol,rhs)
         implicit none
         integer,intent(in) :: n_
         real(WP),intent(in) :: t_
         real(WP), dimension(n_),intent(in)  :: sol
         real(WP), dimension(n_),intent(out) :: rhs
         real(WP) :: Cp_mix,Wmix,RHOmix
         real(WP), dimension(nspec) :: wdot
         ! Reset rhs
         rhs=0.0_WP
         ! Get W of mixture
         call this%get_Wmix(sol(1:nspec),Wmix)
         ! Get Cp of mixture and update hsp
         call this%get_cpmix(sol(1:nspec),sol(nspec+1),Cp_mix)
         ! Get RHO of mixture
         RHOmix=this%Pthermo*Wmix/(Rcst*sol(nspec+1))
         ! Get the reaction source terms
         call fcmech_get_wdot(this%Pthermo,sol(nspec+1),sol(1:nspec),wdot)
         ! Transform concentration into mass fraction
         rhs(1:nspec)=wdot*W_sp/RHOmix
         ! Temperature rhs from change in concentration
         rhs(nspec+1)=-sum(hsp*rhs(1:nspec))/Cp_mix
      end subroutine fc_reaction_compute_rhs

      ! ! ---------------------------------------------------------------------------------- !
      ! ! ========================================================== !
      ! ! Computes the approximate analytical jacobian of the system !
      ! ! ========================================================== !
      ! subroutine fc_reaction_compute_jac(n_,t_,sol,ml,mu,jac,nrpd)
      !     implicit none

      !     ! ! Input
      !     integer :: n_,ml,mu,nrpd
      !     real(WP) :: t_
      !     real(WP), dimension(n_) :: sol
      !     real(WP), dimension(nrpd,n_) :: jac

      !     ! Local variables
      !     real(WP) :: Cp_mix,Wmix,RHOmix
      !     real(WP), dimension(nspec) :: C,Cdot
      !     real(WP), dimension(nTB+nFO) :: M
      !     real(WP), dimension(nreac+nreac_reverse) :: W,K

      !     ! Get W of mixture
      !     call this%get_Wmix(sol(1:nspec),Wmix)

      !     ! Get Cp of mixture and update hsp
      !     call this%get_cpmix(sol(1:nspec),sol(nspec+1),Cp_mix)

      !     ! Get RHO of mixture
      !     RHOmix=this%Pthermo*Wmix/(Rcst*sol(nspec+1))

      !     call get_thirdbodies(M,c)

      !     call get_rate_coefficients(k,M,sol(nspec+1),this%Pthermo)

      !     call get_reaction_rates(w,k,M,c)

      !     call get_production_rates(Cdot,w)

      !     ! Get analytical Jacobian from mechanism file
      !     ! call fc_reaction_getjacobian(sol(1:nspec),sol(nspec+1),M,Cdot,K,RHOmix,Cp_mix,Wmix,jac)

      !     return
      ! end subroutine fc_reaction_compute_jac

      ! ------------------------------------------------------------------------------------------------ !
      ! ================================================ !
      ! Used by finitechem_getjacobian
      ! to compute one of the more involved terms
      ! dealing with pressure dependent rate coefficients
      ! ================================================ !
      subroutine fc_reaction_compute_dlnFCdT(fca_i,fcta_i,fcb_i,fctb_i,fcc_i,fctc_i,Tloc,FC,dlnFCdT)
         implicit none
         real(WP) :: fca_i,fcta_i,fcb_i,fctb_i,fcc_i,fctc_i
         real(WP) :: Tloc,FC,dlnFCdT
         real(WP) :: tmp
         real(WP), parameter :: eps=2.0_WP*epsilon(1.0_WP)
         ! fca_i=1-alpha,fcta_i=T***,fcb_i=alpha,fctb_i=T*,fcc_i=beta,fctc_i=T**
         FC=0.0_WP
         dlnFCdT=0.0_WP
         if (abs(fcta_i).gt.eps) then
            tmp=fca_i*exp(-Tloc/fcta_i)
            FC=FC+tmp
            dlnFCdT=dlnFCdT-tmp/fcta_i
         end if
         if (abs(fctb_i).gt.eps) then
            tmp=fcb_i*exp(-Tloc/fctb_i)
            FC=FC+tmp
            dlnFCdT=dlnFCdT-tmp/fctb_i
         end if
         tmp=fcc_i*exp(-fctc_i/Tloc)
         FC=FC+tmp
         dlnFCdT=dlnFCdT+tmp*fctc_i/(Tloc*Tloc)
         dlnFCdT=dlnFCdT/FC
      end subroutine fc_reaction_compute_dlnFCdT

      ! ------------------------------------------------------------------------------------------------ !
      ! ================================================ !
      ! Used by finitechem_getjacobian
      ! to compute one of the more involved terms
      ! dealing with pressure dependent rate coefficients
      ! ================================================ !
      subroutine fc_reaction_compute_dGdT(redP,kc_oInf,Tloc,FC,dlnFCdT,dGdT)
         implicit none
         real(WP),intent(in) :: redP,kc_oInf,Tloc,FC
         real(WP) :: J,H
         real(WP) :: ctmp,ntmp,dtmp
         real(WP),intent(out) :: dGdT,dlnFCdT
         ntmp=0.75_WP-1.27_WP*dlog10(FC)
         ctmp=-0.4_WP-0.67_WP*dlog10(FC)
         dtmp=0.14_WP
         J=log(redP)+ctmp
         H=J/(ntmp-dtmp*J)
         dlnFCdT=dlnFCdT/(1+H*H)
         dGdT=-2*H/((1+H*H)**2)*(ntmp/((ntmp-dtmp*J)**2))*(kc_oInf-1/Tloc)
      end subroutine fc_reaction_compute_dGdT

   end subroutine react

   
   !> Calculate mixture density
   subroutine get_density(this)
      implicit none
      class(finitechem),intent(inout) :: this
      integer :: i,j,k
      real(WP):: Tmix,Wmix
      real(WP), dimension(nspec) :: Ys
      ! Original version
      ! Compute the new density from the equation of state
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) then
                  this%rho(i,j,k)=1000.0_WP
               else
                  Tmix=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
                  Ys=this%SC(i,j,k,1:nspec)
                  Ys=min(max(Ys,0.0_WP),1.0_WP)
                  Ys=Ys/sum(Ys)
                  call this%get_Wmix(Ys,Wmix)
                  this%rho(i,j,k)=this%Pthermo*Wmix/(Rcst*Tmix)
               end if
            end do
         end do
      end do
   end subroutine get_density


   !> Calculate mixture viscosity
   subroutine get_viscosity(this)
      implicit none
      class(finitechem),intent(inout) :: this
      integer  :: i,j,k,sc1,sc2
      real(WP) :: Tmix,buf
      real(WP), dimension(nspec) :: eta
      real(WP), dimension(nspec,nspec) :: phi
      ! Compute the new viscosity from Wilke's method
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) cycle
               ! Pure compounds viscosity
               Tmix=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
               call fcmech_get_viscosity(eta,Tmix)
               ! Mixing coefficients
               do sc2=1,nspec
                  do sc1=1,nspec
                     if (sc1.eq.sc2) then
                        phi(sc1,sc2)=1.0_WP
                     else
                        buf=sqrt(eta(sc1)/eta(sc2))*(W_sp(sc2)/W_sp(sc1))**0.25_WP
                        phi(sc1,sc2)=(1.0_WP+buf)**2/sqrt(8.0_WP+8.0_WP*W_sp(sc1)/W_sp(sc2))
                     end if
                  end do
               end do
               ! Mixing rule
               this%visc(i,j,k)=0.0_WP
               do sc1=1,nspec
                  if (this%SC(i,j,k,1+sc1-1).le.0.0_WP) cycle
                  buf=sum(this%SC(i,j,k,1:nspec)*phi(sc1,:)/W_sp)
                  this%visc(i,j,k)=this%visc(i,j,k)+this%SC(i,j,k,1+sc1-1)*eta(sc1)/(W_sp(sc1)*buf)
               end do
            end do
         end do
      end do
   end subroutine get_viscosity


   !> Calculate species and thermal diffusivities
   subroutine get_diffusivity(this)
      implicit none
      class(finitechem),intent(inout) :: this
      integer  :: i,j,k,n
      real(WP) :: Wmix,Tmix
      real(WP) :: sum1,sum2,sumY,sumDiff
      real(WP), dimension(nspec) :: eta,cond,YoverW,Ys
      real(WP), dimension(nspec,nspec) :: invDij
      ! Compute the new diffusivity
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) cycle
               ! ---- Thermal diffusivity ---- !
               ! Mixture molar mass and temperature
               Ys=this%SC(i,j,k,1:nspec)
               Ys=min(max(Ys,0.0_WP),1.0_WP)
               Ys=Ys/sum(Ys)
               call this%get_Wmix(Ys,Wmix)
               Tmix=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
               ! Individual compounds viscosity
               call fcmech_get_viscosity(eta,Tmix)
               ! Individual compounds viscosity
               call fcmech_get_conductivity(cond,Tmix,eta)
               ! Mixture averaged thermal conductivity
               sum1=Wmix*sum(Ys/(cond*W_sp))
               sum2=Wmix*sum(Ys*cond/W_sp)
               this%lambda(i,j,k)=0.5_WP*(sum2+1.0_WP/sum1)
               ! Average Cp based on scalar field
               call this%get_cpmix(Ys,Tmix,this%Cp(i,j,k))
               ! Thermal diffusivity for enthalpy
               !  this%diff(i,j,k,isc_ENTH)=this%lambda(i,j,k)/this%Cp(i,j,k)
               ! Thermal diffusivity for temperature
               ! this%diff(i,j,k,nspec+1)=this%lambda(i,j,k)/this%Cp(i,j,k)
               this%diff(i,j,k,:)=this%lambda(i,j,k)/this%Cp(i,j,k)
               ! ! ---- Species diffusivity ---- !
               ! ! Inverse of binary diffusion coefficients
               ! call fcmech_get_invDij(invDij,Tmix,this%Pthermo)
               ! ! Constant terms
               ! Ys=this%SC(i,j,k,1:nspec)
               ! YOverW=Ys/W_sp
               ! sumY=sum(Ys)
               ! ! Compute mixture-average diffusion coefficient for each species
               ! do n=1,nspec
               !    ! Denominator
               !    sumDiff=sum(YOverW*invDij(n,:))
               !    if (sumDiff.gt.1.0e-15_WP) then
               !       ! Diffusion is well defined
               !       this%diff(i,j,k,n)=(sumY-Ys(n))/(Wmix*sumDiff);
               !    else
               !       ! Diffusion is ill defined
               !       sumDiff=sum(invDij(n,:)/W_sp)
               !       this%diff(i,j,k,n)=(real(nspec,WP))/(Wmix*sumDiff)
               !    end if
               !    this%diff(i,j,k,n)=this%rho(i,j,k)*this%diff(i,j,k,n)
               ! end do
            end do
         end do
      end do
   end subroutine get_diffusivity

   
   !> Calculate mixture specific heat capacity at constant pressure
   subroutine get_cpmix(this,scalar,Tmix,Cp_mix)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP), dimension(nspec),intent(in) :: scalar
      real(WP),intent(in) :: Tmix
      real(WP),intent(out) :: Cp_mix
      call fcmech_thermodata(Tmix)
      Cp_mix=sum(scalar*Cpsp)
   end subroutine get_cpmix


   !> Calculate mixture enthalpy
   subroutine get_enthalpy(this)
      implicit none
      integer :: i,j,k,n
      class(finitechem),intent(inout) :: this
      this%h=0.0_WP
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               call fcmech_thermodata(this%T(i,j,k))
               do n=1,nspec
                  this%h(i,j,k)=this%h(i,j,k)+hsp(n)*this%Y(i,j,k,n)
               end do
            end do
         end do
      end do
   end subroutine get_enthalpy


   !> Calculate mixture molar mass
   subroutine get_Wmix(this,scalar,Wmix)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP), dimension(nspec),intent(in) :: scalar
      real(WP),intent(out) :: Wmix
      real(WP), dimension(nspec) :: scalar_clip
      scalar_clip=scalar!min(max(scalar,0.0_WP),1.0_WP)
      Wmix=1.0_WP/sum(scalar_clip/W_sp)
   end subroutine get_Wmix


   !> Calculate diffusion source terms
   subroutine diffusive_source(this,dt)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP), dimension(:,:,:), allocatable :: Wmix,Cpmix,DFX_SUM,DFY_SUM,DFZ_SUM
      real(WP), dimension(:,:,:), allocatable :: FX,FY,FZ
      real(WP), dimension(:,:,:,:), allocatable :: DFX,DFY,DFZ
      real(WP),intent(in) :: dt
      real(WP) :: Tmix,Ttmp,df1,df2,df3
      real(WP), dimension(nspec) :: Ys,hs,Cps
      integer :: i,j,k,nsc

      ! Allocate arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(DFY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(DFZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(Wmix(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(Cpmix(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFX_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFY_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFZ_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Initialize diffusive fluxes
      DFX=0.0_WP
      DFY=0.0_WP
      DFZ=0.0_WP

      ! Get Wmix and Cpmix fields
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) then
                  Wmix(i,j,k)=1.0_WP
                  Cpmix(i,j,k)=1.0_WP
               else
                  Tmix=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
                  Ys=this%SC(i,j,k,1:nspec)
                  call this%get_Wmix(Ys,Wmix(i,j,k))
                  call this%get_cpmix(Ys,Tmix,Cpmix(i,j,k))
               end if
            end do
         end do
      end do

      ! Form species diffusive fluxes
      do nsc=1,nspec
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Molar diffusion correction-DIFF/Wmix*Yi*grad(Wmix)
                  FX(i,j,k)=sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k,nsc)/Wmix(i-1:i,j,k))*sum(this%itp_x(:,i,j,k)*this%SC(i-1:i,j,k,nsc))*sum(this%grdsc_x(:,i,j,k)*Wmix(i-1:i,j,k))
                  FY(i,j,k)=sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k,nsc)/Wmix(i,j-1:j,k))*sum(this%itp_y(:,i,j,k)*this%SC(i,j-1:j,k,nsc))*sum(this%grdsc_y(:,i,j,k)*Wmix(i,j-1:j,k))
                  FZ(i,j,k)=sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k,nsc)/Wmix(i,j,k-1:k))*sum(this%itp_z(:,i,j,k)*this%SC(i,j,k-1:k,nsc))*sum(this%grdsc_z(:,i,j,k)*Wmix(i,j,k-1:k))
                  ! Store full diffusive flux
                  DFX(i,j,k,nsc)=FX(i,j,k)+sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k,nsc))*sum(this%grdsc_x(:,i,j,k)*this%SC(i-1:i,j,k,nsc))
                  DFY(i,j,k,nsc)=FY(i,j,k)+sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k,nsc))*sum(this%grdsc_y(:,i,j,k)*this%SC(i,j-1:j,k,nsc))
                  DFZ(i,j,k,nsc)=FZ(i,j,k)+sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k,nsc))*sum(this%grdsc_z(:,i,j,k)*this%SC(i,j,k-1:k,nsc))
               end do
            end do
         end do
         ! Update species source term
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%SRC(i,j,k,nsc)=this%SRC(i,j,k,nsc)+dt*(sum(this%divsc_x(i,j,k,:)*FX(i:i+1,j,k)) &
                  &                                          +sum(this%divsc_y(i,j,k,:)*FY(i,j:j+1,k)) &
                  &                                          +sum(this%divsc_z(i,j,k,:)*FZ(i,j,k:k+1)))
               end do
            end do
         end do
      end do

      ! Correct species diffusion to ensure sum(DFX)=0, sum(DFY)=0, sum(DFZ)=0
      do k=this%cfg%kmin_,this%cfg%kmax_+1
         do j=this%cfg%jmin_,this%cfg%jmax_+1
            do i=this%cfg%imin_,this%cfg%imax_+1
               DFX_SUM(i,j,k)=sum(DFX(i,j,k,:))
               DFY_SUM(i,j,k)=sum(DFY(i,j,k,:))
               DFZ_SUM(i,j,k)=sum(DFZ(i,j,k,:))
            end do
         end do
      end do

      do nsc=1,nspec
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Diffusion correction: -Yi*sum(DFX)
                  FX(i,j,k)=-sum(this%itp_x(:,i,j,k)*this%SC(i-1:i,j,k,nsc))*DFX_SUM(i,j,k)
                  FZ(i,j,k)=-sum(this%itp_y(:,i,j,k)*this%SC(i,j-1:j,k,nsc))*DFY_SUM(i,j,k)
                  FY(i,j,k)=-sum(this%itp_z(:,i,j,k)*this%SC(i,j,k-1:k,nsc))*DFZ_SUM(i,j,k)
                  ! Update full diffusive flux
                  DFX(i,j,k,nsc)=DFX(i,j,k,nsc)+FX(i,j,k)
                  DFY(i,j,k,nsc)=DFY(i,j,k,nsc)+FY(i,j,k)
                  DFZ(i,j,k,nsc)=DFZ(i,j,k,nsc)+FZ(i,j,k)
               end do
            end do
         end do
         ! Update species source term
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  this%SRC(i,j,k,nsc)=this%SRC(i,j,k,nsc)+dt*(sum(this%divsc_x(i,j,k,:)*FX(i:i+1,j,k)) &
                  &                                          +sum(this%divsc_y(i,j,k,:)*FY(i,j:j+1,k)) &
                  &                                          +sum(this%divsc_z(i,j,k,:)*FZ(i,j,k:k+1)))
               end do
            end do
         end do
      end do

      ! Form thermal diffusion correction-lambda/Cp^2*grad(Cpmix).grad(T)
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+dt*this%diff(i,j,k,nspec+1)/Cpmix(i,j,k)*(sum(this%grdsc_xm(:,i,j,k)*Cpmix(i:i+1,j,k))*sum(this%grdsc_xm(:,i,j,k)*this%SC(i:i+1,j,k,nspec+1)) &
               &                                                                                        +sum(this%grdsc_ym(:,i,j,k)*Cpmix(i,j:j+1,k))*sum(this%grdsc_ym(:,i,j,k)*this%SC(i,j:j+1,k,nspec+1)) &
               &                                                                                        +sum(this%grdsc_zm(:,i,j,k)*Cpmix(i,j,k:k+1))*sum(this%grdsc_zm(:,i,j,k)*this%SC(i,j,k:k+1,nspec+1)))
            end do
         end do
      end do

      ! Update temperature source terms
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).eq.1) cycle
               ! Update Cp of species
               call fcmech_get_thermodata(hs,Cps,this%SC(i,j,k,nspec+1))
               ! Loop over species and compute temperature flux-sum(Cpsp*DF.grad(T))/Cpmix
               df1=0.0_WP; df2=0.0_WP; df3=0.0_WP
               do nsc=1,nspec
                  df1=df1+Cps(nsc)*sum(this%itp_x(:,i,j,k)*DFX(i-1:i,j,k,nsc))
                  df2=df2+Cps(nsc)*sum(this%itp_y(:,i,j,k)*DFY(i,j-1:j,k,nsc))
                  df3=df3+Cps(nsc)*sum(this%itp_z(:,i,j,k)*DFZ(i,j,k-1:k,nsc))
               end do
               ! Temperature source
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+dt/Cpmix(i,j,k)*(df1*sum(this%grdsc_xm(:,i,j,k)*this%SC(i:i+1,j,k,nspec+1)) &
               &                                                               +df2*sum(this%grdsc_ym(:,i,j,k)*this%SC(i,j:j+1,k,nspec+1)) &
               &                                                               +df3*sum(this%grdsc_zm(:,i,j,k)*this%SC(i,j,k:k+1,nspec+1)))
            end do
         end do
      end do

      ! Release memory
      deallocate (FX,FY,FZ,DFX,DFY,DFZ,Cpmix,Wmix,DFX_SUM,DFY_SUM,DFZ_SUM)
   end subroutine diffusive_source


   !> Calculate pressure source term
   subroutine pressure_source(this)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP) :: Cp_mix,Tmix
      integer :: i,j,k
      real(WP), dimension(nspec) :: Ys
      ! Compute pressure source term
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).eq.1) cycle
               Tmix=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
               Ys=this%SC(i,j,k,1:nspec)
               call this%get_cpmix(Ys,Tmix,Cp_mix)
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+(this%Pthermo-this%Pthermo_old)/Cp_mix
            end do
         end do
      end do
   end subroutine pressure_source


   !> Update thermodynamic pressure and density
   subroutine update_pressure(this)
      implicit none
      class(finitechem),intent(inout) :: this
      integer :: i,j,k
      ! Save the old background pressure
      this%Pthermo_old=this%Pthermo
      ! Recompute mean density
      call this%cfg%integrate(this%rho,integral=this%rhoint)
      this%RHOmean=this%rhoint/this%cfg%vol_total
      ! Update Pthermo
      this%Pthermo=this%Pthermo*this%RHO_0/this%RHOmean
      ! Update density
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) cycle
               this%RHO(i,j,k)=this%RHO(i,j,k)*(this%RHO_0/this%RHOmean)
            end do
         end do
      end do
   end subroutine update_pressure


   !> Calculate all the source terms
   subroutine get_src(this,dt)
      implicit none
      class(finitechem),intent(inout) :: this
      real(WP),intent(in) :: dt
      integer :: nsc
      ! Add chemical source terms
      do nsc=1,nspec+1
         this%SRC(:,:,:,nsc)=this%rho*this%SRCchem(:,:,:,nsc)
      end do
      ! Get pressure source term
      ! call this%pressure_source()
      ! Get diffusion source terms
      call this%diffusive_source(dt)
   end subroutine get_src


   !> Calculate the min and max of SC fields
   subroutine fc_get_max(this)
      use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(finitechem),intent(inout) :: this
      integer :: ierr,i,j,k,nsc
      real(WP) :: my_visc_max,my_visc_min,my_rhomax,my_rhomin,my_SCmax,my_SCmin,my_rhoSCmax,my_rhoSCmin,my_diff_max,my_diff_min
      my_SCmax=-huge(1.0_WP)
      my_SCmin=+huge(1.0_WP)
      my_rhomax=-huge(1.0_WP)
      my_rhomin=+huge(1.0_WP)
      my_rhoSCmax=-huge(1.0_WP)
      my_rhoSCmin=+huge(1.0_WP)
      my_rhomax=maxval(this%rho(:,:,:)); call MPI_ALLREDUCE(my_rhomax,this%rhomax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_rhomin=minval(this%rho(:,:,:)); call MPI_ALLREDUCE(my_rhomin,this%rhomin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      my_visc_max=-huge(1.0_WP)
      my_visc_min=+huge(1.0_WP)
      my_diff_max=-huge(1.0_WP)
      my_diff_min=+huge(1.0_WP)
      do nsc=1,this%nscalar
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  ! Skip only walls
                  if (this%mask(i,j,k).ne.1) then
                     my_SCmax=max(this%SC(i,j,k,nsc),my_SCmax)
                     my_SCmin=min(this%SC(i,j,k,nsc),my_SCmin)
                     my_rhoSCmax=max(this%rhoSC(i,j,k,nsc),my_rhoSCmax)
                     my_rhoSCmin=min(this%rhoSC(i,j,k,nsc),my_rhoSCmin)
                     my_diff_max=max(this%diff(i,j,k,nsc),my_diff_max)
                     my_diff_min=min(this%diff(i,j,k,nsc),my_diff_min)
                  end if
               end do
            end do
         end do
         call MPI_ALLREDUCE(my_SCmax,this%SCmax(nsc),1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(my_SCmin,this%SCmin(nsc),1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(my_rhoSCmax,this%rhoSCmax(nsc),1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
         call MPI_ALLREDUCE(my_rhoSCmin,this%rhoSCmin(nsc),1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      end do
      call MPI_ALLREDUCE(my_diff_max,this%diff_max,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      call MPI_ALLREDUCE(my_diff_min,this%diff_min,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      my_visc_max=maxval(this%visc); call MPI_ALLREDUCE(my_visc_max,this%visc_max,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_visc_min=minval(this%visc); call MPI_ALLREDUCE(my_visc_min,this%visc_min,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
      my_rhomax=maxval(this%rho(:,:,:)); call MPI_ALLREDUCE(my_rhomax,this%rhomax,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr)
      my_rhomin=minval(this%rho(:,:,:)); call MPI_ALLREDUCE(my_rhomin,this%rhomin,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr)
   end subroutine fc_get_max


end module finitechem_class
