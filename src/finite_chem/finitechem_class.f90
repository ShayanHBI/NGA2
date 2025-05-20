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

   ! Inverse molar mass of the species
   real(WP), dimension(:), allocatable :: Winv

   !> Finite chemistry solver object definition
   type, extends(multivdscalar) :: finitechem
      
      ! Thermodynamic pressure and density
      real(WP) :: Pthermo,Pthermo_old
      real(WP) :: RHOmean,RHO_0
      
      ! Scalar variables
      real(WP), dimension(:,:,:),   allocatable :: visc                       !< Viscosity field
      real(WP), dimension(:,:,:,:), allocatable :: SRCchem                    !< Chemical source terms
      real(WP), dimension(:,:,:,:), allocatable :: SRC                        !< Total source terms for scalar equations
      real(WP), dimension(:,:,:),   allocatable :: Cp                         !< Mixture heat capacity
      real(WP), dimension(:,:,:),   allocatable :: W                          !< Mixture molar mass
      
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
      procedure :: get_molarMass
      procedure :: get_Cp
      procedure :: get_density
      procedure :: get_visc_diff
      procedure :: diffusive_src
      procedure :: pressure_src
      procedure :: update_pressure
      procedure :: get_src
      procedure :: mixture_avg
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
      class(config),target, intent(in) :: cfg
      integer, intent(in) :: scheme
      character(len=*), optional :: name
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
      self%SCname(nspec+1)='temperature'

      ! Allocate variables
      allocate(self%visc(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));            self%visc=0.0_WP
      allocate(self%Cp(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));              self%Cp=0.0_WP
      allocate(self%W(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_));               self%W=0.0_WP
      allocate(self%SRCchem(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,nspec+1)); self%SRCchem=0.0_WP
      allocate(self%SRC(self%cfg%imino_:self%cfg%imaxo_,self%cfg%jmino_:self%cfg%jmaxo_,self%cfg%kmino_:self%cfg%kmaxo_,nspec+1));     self%SRC=0.0_WP
      allocate(Winv(nspec)); Winv=1.0_WP/W_sp

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

   end function constructor


   !> Initializer for the scheduler
   subroutine scheduler_init(this,bundleref)
      use mpi_f08, only: MPI_MAX,MPI_INTEGER
      implicit none
      class(finitechem), intent(inout) :: this
      real(WP), intent(in) :: bundleref
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
      class(finitechem), intent(inout) :: this
      real(WP), dimension(nspec), intent(inout) :: myY
      myY=min(max(myY,0.0_WP),1.0_WP)
      myY=myY/sum(myY)
   end subroutine clip


   !> Calculate reaction source terms
   subroutine react(this,dt)
      use parallel, only: MPI_REAL_WP
      use mpi_f08,  only: MPI_TAG,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_SOURCE,MPI_STATUS_SIZE,MPI_INTEGER,MPI_LOGICAL
      use messager, only: die
      implicit none
      class(finitechem), intent(inout) :: this
      real(WP), intent(in) :: dt
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

      ! Initialize with zeros
      this%SRCchem=0.0_WP

      ! If only one processor or if beginning of simulation
      if (.not.this%use_scheduler.or.this%cfg%nproc.eq.1) then
         do k=this%cfg%kmin_,this%cfg%kmax_
            do j=this%cfg%jmin_,this%cfg%jmax_
               do i=this%cfg%imin_,this%cfg%imax_
                  sol(1:nspec)=this%SC(i,j,k,1:nspec)
                  ! Clip and renormalize
                  call this%clip(sol(1:nspec))
                  sol(nspec+1)=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
                  ! Remember old solution
                  solold=sol
                  if (sol(sN2).gt.0.8_WP) cycle
                  ! Advance the reactions
                  call get_sol(sol)
                  ! Calculate the scalar chemical source terms
                  this%SRCchem(i,j,k,:)=sol-solold
               end do
            end do
         end do
         ! Sync
         do nsc=1,nspec+1
            call this%cfg%sync(this%SRCchem(:,:,:,nsc))
         end do
         ! Stop there, no dynamic scheduling in this case
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
      !                call get_sol(this%bufferR((i-1)*(nspec+1)+1:i*(nspec+1)))
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
      subroutine get_sol(mysol)
         use random
         use messager, only: die
         implicit none
         real(WP), dimension(nspec+1), intent(inout) :: mysol
         real(WP), dimension(nspec+1) :: dsol,solcheck
         real(WP) :: t_chem
         real(WP) :: tstop,tstart
         type(vode_opts), save :: opts
         integer  :: istate,itask
         real(WP) :: atol=1.0e-15_WP
         real(WP) :: rtol=1.0e-12_WP
         integer :: n
         real(WP), dimension(22) :: rstats
         integer, dimension(31) :: istats
         ! ISAT parameters
         real(WP) :: xx(nspec+2),f(nspec+1),dfdx(nspec+1,nspec+2),hvar(1),stats(100)
         integer :: iusr(1)
         ! Save solution
         solcheck=mysol
         ! Set values for DVODE
         itask=1; istate=1; tstart=0.0_WP; tstop=dt
         ! Direct integration
         call get_rhs(nspec+1,0.0_WP,solcheck,dsol)
         t_chem=mysol(nspec+1)/(abs(dsol(nspec+1))+epsilon(1.0_WP))
         if (dt.lt.0.0001_WP*t_chem.and.this%use_explicit_try) then
            mysol=mysol+dsol*dt
         else
            ! Integrate using DVODE
            if (use_jacanal) then
               continue
               ! opts=set_opts(dense_j=.true.,mxstep=500000,abserr=atol,relerr=rtol,tcrit=tstop,user_supplied_jacobian=.true.)
               !  call dvode_f90(get_rhs,nspec+1,mysol,tstart,tstop,itask,istate,opts,j_fcn=fc_reaction_compute_jac)
            else
               opts=set_opts(method_flag=22,mxstep=500000,abserr=atol,relerr=rtol,tcrit=tstop)
               call dvode_f90(get_rhs,nspec+1,mysol,tstart,tstop,itask,istate,opts)
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
               print *,'sol',mysol
               print *,'-------------------------------'
               do n=1,nspec
                  print *,'dsol(',n,')=',dsol(n),'_WP !'
               end do
               print *,'T','dsol(n)',dsol(nspec+1)
               call die('fc_reaction_source: Direct integration-DVODE failed to converge.')
            end if
         end if
         ! Clip and renormalize
         call this%clip(mysol(1:nspec))
      end subroutine get_sol

      ! Computes the chemical source term of the system (called by solver)
      subroutine get_rhs(n_,t_,mysol,rhs)
         implicit none
         integer, intent(in) :: n_
         real(WP), intent(in) :: t_
         real(WP), dimension(n_), intent(in)  :: mysol
         real(WP), dimension(n_), intent(out) :: rhs
         real(WP), dimension(nspec) :: wdot
         ! Reset rhs
         rhs=0.0_WP
         ! Get the reaction source terms
         call fcmech_get_wdot(this%Pthermo,mysol(nspec+1),mysol(1:nspec),wdot)
         ! Transform concentration into mass fraction
         rhs(1:nspec)=wdot*W_sp/this%rho(i,j,k)
         ! Temperature rhs from change in concentration
         rhs(nspec+1)=-sum(hsp*rhs(1:nspec))/this%Cp(i,j,k)
      end subroutine get_rhs

   end subroutine react


   !> Calculate mixture molar mass
   subroutine get_molarMass(this)
      implicit none
      class(finitechem), intent(inout) :: this
      integer :: i,j,k
      real(WP), dimension(nspec) :: Y
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               Y=this%SC(i,j,k,1:nspec)
               call this%clip(Y)
               this%W(i,j,k)=1.0_WP/this%mixture_avg(Winv,Y)
            end do
         end do
      end do
   end subroutine get_molarMass


   !> Calculate mixture heat capacity at constant pressure
   subroutine get_Cp(this)
      implicit none
      class(finitechem), intent(inout) :: this
      integer :: i,j,k
      real(WP):: T
      real(WP), dimension(nspec) :: Y
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               Y=this%SC(i,j,k,1:nspec)
               call this%clip(Y)
               T=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
               call fcmech_thermodata(T)
               this%Cp(i,j,k)=this%mixture_avg(Cpsp,Y)
            end do
         end do
      end do
   end subroutine get_Cp

   
   !> Calculate mixture density
   subroutine get_density(this)
      implicit none
      class(finitechem), intent(inout) :: this
      integer :: i,j,k
      real(WP):: T
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) then
                  this%rho(i,j,k)=1000.0_WP
               else
                  T=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
                  this%rho(i,j,k)=this%Pthermo*this%W(i,j,k)/(Rcst*T)
               end if
            end do
         end do
      end do
   end subroutine get_density


   !> Calculate mixture viscosity and diffusivity
   subroutine get_visc_diff(this)
      implicit none
      class(finitechem), intent(inout) :: this
      integer  :: i,j,k,nsc1,nsc2
      real(WP) :: T,buf
      real(WP) :: lambda1,lambda2
      real(WP), dimension(nspec) :: Y,eta,cond
      real(WP), dimension(nspec,nspec) :: phi
      do k=this%cfg%kmino_,this%cfg%kmaxo_
         do j=this%cfg%jmino_,this%cfg%jmaxo_
            do i=this%cfg%imino_,this%cfg%imaxo_
               if (this%mask(i,j,k).eq.1) cycle
               ! Composition and temperature
               Y=this%SC(i,j,k,1:nspec)
               call this%clip(Y)
               T=min(max(this%SC(i,j,k,nspec+1),T_min),T_max)
               ! Pure compounds viscosity
               call fcmech_get_viscosity(eta,T)
               ! Mixing coefficients
               do nsc2=1,nspec
                  do nsc1=1,nspec
                     if (nsc1.eq.nsc2) then
                        phi(nsc1,nsc2)=1.0_WP
                     else
                        buf=sqrt(eta(nsc1)/eta(nsc2))*(W_sp(nsc2)/W_sp(nsc1))**0.25_WP
                        phi(nsc1,nsc2)=(1.0_WP+buf)**2/sqrt(8.0_WP+8.0_WP*W_sp(nsc1)/W_sp(nsc2))
                     end if
                  end do
               end do
               ! Get viscosity from Wilke's method
               this%visc(i,j,k)=0.0_WP
               do nsc1=1,nspec
                  if (this%SC(i,j,k,1+nsc1-1).le.0.0_WP) cycle
                  buf=sum(this%SC(i,j,k,1:nspec)*phi(nsc1,:)/W_sp)
                  this%visc(i,j,k)=this%visc(i,j,k)+this%SC(i,j,k,1+nsc1-1)*eta(nsc1)/(W_sp(nsc1)*buf)
               end do
               ! Individual compounds conductivity
               call fcmech_get_conductivity(cond,T,eta)
               ! Mixture averaged thermal conductivity
               lambda1=this%W(i,j,k)*sum(Y/(cond*W_sp))
               lambda2=this%W(i,j,k)*sum(Y*cond/W_sp)
               ! Unity Lewis number
               this%diff(i,j,k,:)=0.5_WP*(lambda2+1.0_WP/lambda1)/this%Cp(i,j,k)
            end do
         end do
      end do
   end subroutine get_visc_diff


   !> Calculate diffusion source terms
   subroutine diffusive_src(this,dt)
      implicit none
      class(finitechem), intent(inout) :: this
      real(WP), dimension(:,:,:),   allocatable :: DFX_SUM,DFY_SUM,DFZ_SUM
      real(WP), dimension(:,:,:),   allocatable :: FX,FY,FZ
      real(WP), dimension(:,:,:,:), allocatable :: DFX,DFY,DFZ
      real(WP), intent(in) :: dt
      real(WP) :: df1,df2,df3
      integer :: i,j,k,nsc

      ! Allocate flux arrays
      allocate(FX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(FZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFX(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(DFY(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(DFZ(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_,1:nspec))
      allocate(DFX_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFY_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))
      allocate(DFZ_SUM(this%cfg%imino_:this%cfg%imaxo_,this%cfg%jmino_:this%cfg%jmaxo_,this%cfg%kmino_:this%cfg%kmaxo_))

      ! Initialize diffusive fluxes
      DFX=0.0_WP
      DFY=0.0_WP
      DFZ=0.0_WP

      ! Form species diffusive fluxes
      do nsc=1,nspec
         do k=this%cfg%kmin_,this%cfg%kmax_+1
            do j=this%cfg%jmin_,this%cfg%jmax_+1
               do i=this%cfg%imin_,this%cfg%imax_+1
                  ! Molar diffusion correction-DIFF/Wmix*Yi*grad(Wmix)
                  FX(i,j,k)=sum(this%itp_x(:,i,j,k)*this%diff(i-1:i,j,k,nsc)/this%W(i-1:i,j,k))*sum(this%itp_x(:,i,j,k)*this%SC(i-1:i,j,k,nsc))*sum(this%grdsc_x(:,i,j,k)*this%W(i-1:i,j,k))
                  FY(i,j,k)=sum(this%itp_y(:,i,j,k)*this%diff(i,j-1:j,k,nsc)/this%W(i,j-1:j,k))*sum(this%itp_y(:,i,j,k)*this%SC(i,j-1:j,k,nsc))*sum(this%grdsc_y(:,i,j,k)*this%W(i,j-1:j,k))
                  FZ(i,j,k)=sum(this%itp_z(:,i,j,k)*this%diff(i,j,k-1:k,nsc)/this%W(i,j,k-1:k))*sum(this%itp_z(:,i,j,k)*this%SC(i,j,k-1:k,nsc))*sum(this%grdsc_z(:,i,j,k)*this%W(i,j,k-1:k))
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
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+dt*this%diff(i,j,k,nspec+1)/this%Cp(i,j,k)*(sum(this%grdsc_xm(:,i,j,k)*this%Cp(i:i+1,j,k))*sum(this%grdsc_xm(:,i,j,k)*this%SC(i:i+1,j,k,nspec+1)) &
               &                                                                                          +sum(this%grdsc_ym(:,i,j,k)*this%Cp(i,j:j+1,k))*sum(this%grdsc_ym(:,i,j,k)*this%SC(i,j:j+1,k,nspec+1)) &
               &                                                                                          +sum(this%grdsc_zm(:,i,j,k)*this%Cp(i,j,k:k+1))*sum(this%grdsc_zm(:,i,j,k)*this%SC(i,j,k:k+1,nspec+1)))
            end do
         end do
      end do

      ! Update temperature source term
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).eq.1) cycle
               ! Update Cp of species
               call fcmech_thermodata(this%SC(i,j,k,nspec+1))
               ! Loop over species and compute temperature flux-sum(Cpsp*DF.grad(T))/Cpmix
               df1=0.0_WP; df2=0.0_WP; df3=0.0_WP
               do nsc=1,nspec
                  df1=df1+Cpsp(nsc)*sum(this%itp_x(:,i,j,k)*DFX(i-1:i,j,k,nsc))
                  df2=df2+Cpsp(nsc)*sum(this%itp_y(:,i,j,k)*DFY(i,j-1:j,k,nsc))
                  df3=df3+Cpsp(nsc)*sum(this%itp_z(:,i,j,k)*DFZ(i,j,k-1:k,nsc))
               end do
               ! Temperature source
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+dt/this%Cp(i,j,k)*(df1*sum(this%grdsc_xm(:,i,j,k)*this%SC(i:i+1,j,k,nspec+1)) &
               &                                                                 +df2*sum(this%grdsc_ym(:,i,j,k)*this%SC(i,j:j+1,k,nspec+1)) &
               &                                                                 +df3*sum(this%grdsc_zm(:,i,j,k)*this%SC(i,j,k:k+1,nspec+1)))
            end do
         end do
      end do

      ! Deallocate flux arrays
      deallocate (FX,FY,FZ,DFX,DFY,DFZ,DFX_SUM,DFY_SUM,DFZ_SUM)
   end subroutine diffusive_src


   !> Calculate pressure source term
   subroutine pressure_src(this)
      implicit none
      class(finitechem), intent(inout) :: this
      integer :: i,j,k
      ! Compute pressure source term
      do k=this%cfg%kmin_,this%cfg%kmax_
         do j=this%cfg%jmin_,this%cfg%jmax_
            do i=this%cfg%imin_,this%cfg%imax_
               if (this%mask(i,j,k).eq.1) cycle
               this%SRC(i,j,k,nspec+1)=this%SRC(i,j,k,nspec+1)+(this%Pthermo-this%Pthermo_old)/this%Cp(i,j,k)
            end do
         end do
      end do
   end subroutine pressure_src


   !> Update thermodynamic pressure and density
   subroutine update_pressure(this)
      implicit none
      class(finitechem), intent(inout) :: this
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
      class(finitechem), intent(inout) :: this
      real(WP), intent(in) :: dt
      integer :: nsc
      ! Add chemical source terms
      do nsc=1,nspec+1
         this%SRC(:,:,:,nsc)=this%rho*this%SRCchem(:,:,:,nsc)
      end do
      ! Get pressure source term
      ! call this%pressure_src()
      ! Get diffusion source terms
      call this%diffusive_src(dt)
   end subroutine get_src


   !> Calculate the mixture-averaged of a given quantity
   function mixture_avg(this,input,Y)
      implicit none
      class(finitechem), intent(in) :: this
      real(WP), dimension(nspec), intent(in) :: input,Y
      real(WP) :: mixture_avg
      mixture_avg=sum(Y*input)
   end function mixture_avg


   !> Calculate the min and max of SC fields
   subroutine fc_get_max(this)
      use mpi_f08, only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN
      use parallel, only: MPI_REAL_WP
      implicit none
      class(finitechem), intent(inout) :: this
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
