!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use timetracker_class, only: timetracker
   use shockslab_class,   only: shockslab
   use event_class,       only: event
   use monitor_class,     only: monitor
   use timer_class,       only: timer
   implicit none
   private; public :: simulation_init,simulation_run,simulation_final

   !> Track time from here
   type(timetracker) :: time

   !> Shock-slab simulation
   type(shockslab) :: ss

   !> Ensight output event
   type(event) :: ens_evt

   !> Timing
   type(monitor) :: timefile !< Timing monitoring
   type(timer)   :: tstep    !< Timer for step

   !> Equations of state
   real(WP) :: PinfL,GammaL,CvL
   real(WP) :: PinfG,GammaG,CvG

   !> Flow parameters
   real(WP) :: Ms,Xs
   real(WP) :: rho1,p1,u1,M1
   real(WP) :: rho2,p2,u2,M2
   real(WP) :: rho_ratio,c_ratio
   real(WP) :: rhoL,ML

   !> Slab location
   real(WP) :: slab_left,slab_right

contains


   !> Function that returns a smooth Heaviside of thickness delta
   real(WP) function Hshock(x,delta)
      real(WP), intent(in) :: x,delta
      ! Goes from 0 to 1 as x goes from negative to positive
      Hshock=1.0_WP/(1.0_WP+exp(-x/delta))
   end function Hshock


   !> P=EOS(RHO,I) for liquid
   pure real(WP) function get_PL(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PL=RHO*I*(GammaL-1.0_WP)-GammaL*PinfL
   end function get_PL
   !> T=f(RHO,P) for liquid
   pure real(WP) function get_TL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TL=(P+PinfL)/(CvL*RHO*(GammaL-1.0_WP))
   end function get_TL
   !> C=f(RHO,P) for liquid
   pure real(WP) function get_CL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CL=sqrt(GammaL*(P+PinfL)/RHO)
   end function get_CL
   !> S=f(RHO,P) for liquid
   pure real(WP) function get_SL(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SL=CvL*log((P+PinfL)/RHO**GammaL)
   end function get_SL


   !> P=EOS(RHO,I) for gas
   pure real(WP) function get_PG(RHO,I)
      implicit none
      real(WP), intent(in) :: RHO,I
      get_PG=RHO*I*(GammaG-1.0_WP)-GammaG*PinfG
   end function get_PG
   !> T=f(RHO,P) for gas
   pure real(WP) function get_TG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_TG=(P+PinfG)/(CvG*RHO*(GammaG-1.0_WP))
   end function get_TG
   !> C=f(RHO,P) for gas
   pure real(WP) function get_CG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_CG=sqrt(GammaG*(P+PinfG)/RHO)
   end function get_CG
   !> S=f(RHO,P) for gas
   pure real(WP) function get_SG(RHO,P)
      implicit none
      real(WP), intent(in) :: RHO,P
      get_SG=CvG*log((P+PinfG)/RHO**GammaG)
   end function get_SG


   !> Mechanical relaxation model (implicit)
   subroutine P_relax_implicit(VF,Q)
      implicit none
      real(WP),                intent(inout) :: VF
      real(WP), dimension(1:), intent(inout) :: Q
      real(WP) :: a,b,d,d1,d0,Peq,VFeq,invG1G,invG1L,facG,facL
      real(WP), parameter :: RHOGmin=1.0e-3_WP
      ! Handle gas flotsams
      if (Q(2)/(1.0_WP-VF).lt.RHOGmin) return
      ! Setup quadratic problem
      invG1G=1.0_WP/(GammaG-1.0_WP); invG1L=1.0_WP/(GammaL-1.0_WP)
      d0=PinfL*GammaL*invG1L; d1=1.0_WP+invG1L
      facG=GammaG*PinfG*invG1G; facL=invG1G+VF
      a=d1*facL-VF*(invG1G+1.0_WP)
      b=d1*(facG-Q(4))-VF*facG+d0*facL-Q(3)*(invG1G+1.0_WP)
      d=d0*(facG-Q(4))-Q(3)*facG
      ! Get equilibrium pressure
      Peq=(-b+sqrt(b**2-4.0_WP*a*d))/(2.0_WP*a)
      ! Check if pressure is sound
      if (Peq.le.max(-PinfG,-PinfL)) return
      ! Get equilibrium volume fraction
      VFeq=(VF*Peq+Q(3))/(d1*Peq+d0)
      ! Adjust conserved quantities
      Q(3)=Q(3)-Peq*(VFeq-VF)
      Q(4)=Q(4)+Peq*(VFeq-VF)
      VF=VFeq
   end subroutine P_relax_implicit


   !> Solver initialization
   subroutine simulation_init
      implicit none

      ! Initialize eos and flow parameters - all cores
      initialize_parameters: block
         use string,   only: str_long
         use messager, only: log
         use parallel, only: amRoot
         use param,    only: param_read
         character(str_long) :: message
         ! Set PinfG to zero
         PinfG=0.0_WP
         ! Read in Gammas
         call param_read('Liquid gamma',GammaL)
         call param_read('Gas gamma'   ,GammaG)
         ! Read in shock Mach number and location
         call param_read('Shock Mach number',Ms)
         call param_read('Shock location',Xs)
         ! Read in slab location
         call param_read('Slab left',slab_left)
         call param_read('Slab right',slab_right)
         ! First generate static shock with normalized pre-shock conditions
         M1=Ms
         rho1=1.0_WP
         rho2=rho1*(GammaG+1.0_WP)*M1**2/((GammaG-1.0_WP)*M1**2+2.0_WP)
         p1=0.25_WP*rho1/GammaG*((GammaG+1.0_WP)*M1/(M1**2-1.0_WP))**2 ! Ensures that |u2-u1|=1
         p2=p1*(2.0_WP*GammaG/(GammaG+1.0_WP)*(M1**2-1.0_WP)+1.0_WP)
         u1=M1*sqrt(GammaG*p1/rho1)
         u2=u1*rho1/rho2
         ! Now shift frame of reference to obtain moving shock
         u2=abs(u2-u1); M2=u2/sqrt(GammaG*p2/rho2); u1=0.0_WP; M1=u1/sqrt(GammaG*p1/rho1)
         ! Read in density ratio and use it to set liquid density
         call param_read('Density ratio',rho_ratio); rhoL=rho_ratio*rho1
         ! Read in sound speed ratio and use it to set PinfL
         call param_read('Sound speed ratio',c_ratio)
         PinfL=p1*(rho_ratio*c_ratio**2*GammaG/GammaL-1.0_WP)
         ML=u2/sqrt(GammaL*(p1+PinfL)/rhoL)
         ! Set heat capacities corresponding to a normalized pre-shock and liquid temperature
         CvL=(p1+PinfL)/(rhoL*(GammaL-1.0_WP))
         CvG=(p1+PinfG)/(rho1*(GammaG-1.0_WP))
         ! Output case info
         if (amRoot) then
            write(message,'("[Liquid EOS] => Gamma=",es12.5)') GammaL; call log(message)
            write(message,'("[Liquid EOS] =>  Pinf=",es12.5)')  PinfL; call log(message)
            write(message,'("[Liquid EOS] =>    Cv=",es12.5)')    CvL; call log(message)
            write(message,'("[Gas EOS]    => Gamma=",es12.5)') GammaG; call log(message)
            write(message,'("[Gas EOS]    =>    Cv=",es12.5)')    CvG; call log(message)
            write(message,'("[Shock Mach number]     =>     Ms=",es12.5)')     Ms; call log(message)
            write(message,'("[Pre -shock conditions] =>   rho1=",es12.5)')   rho1; call log(message)
            write(message,'("[Pre -shock conditions] =>     p1=",es12.5)')     p1; call log(message)
            write(message,'("[Pre -shock conditions] =>     u1=",es12.5)')     u1; call log(message)
            write(message,'("[Pre -shock conditions] =>     M1=",es12.5)')     M1; call log(message)
            write(message,'("[Post-shock conditions] =>   rho2=",es12.5)')   rho2; call log(message)
            write(message,'("[Post-shock conditions] =>     p2=",es12.5)')     p2; call log(message)
            write(message,'("[Post-shock conditions] =>     u2=",es12.5)')     u2; call log(message)
            write(message,'("[Post-shock conditions] =>     M2=",es12.5)')     M2; call log(message)
            write(message,'("[Liquid Mach number] =>        ML=",es12.5)')     ML; call log(message)
            write(message,'("[Density ratio]      => rhoL/rho1=",es12.5)') rho_ratio; call log(message)
            write(message,'("[Sound speed ratio]  =>     cl/c1=",es12.5)')   c_ratio; call log(message)
            write(message,'("[Slab left]          =>          =",es12.5)') slab_left;  call log(message)
            write(message,'("[Slab right]         =>          =",es12.5)') slab_right; call log(message)
         end if
      end block initialize_parameters

      ! Initialize time tracker
      initialize_timetracker: block
         use parallel, only: amRoot
         use param,    only: param_read
         ! Create time tracker object
         time=timetracker(amRoot=amRoot,name='Global')
         ! Set time integration parameters
         call param_read('Max dt',time%dtmax); time%dt=time%dtmax
         call param_read('Max CFL',time%cflmax)
         call param_read('Max time',time%tmax)
      end block initialize_timetracker

      ! Setup shock-slab simulation - all cores
      setup_ss: block
         use param,    only: param_read
         use parallel, only: group
         real(WP), dimension(3) :: X0
         integer , dimension(3) :: meshsize,partition
         real(WP) :: dx
         ! Read in mesh parameters
         call param_read('dx',dx)
         call param_read('nx',meshsize(1))
         call param_read('ny',meshsize(2))
         call param_read('nz',meshsize(3))
         call param_read('X0',X0(1))
         X0(2)=-0.5_WP*real(meshsize(2),WP)*dx
         X0(3)=-0.5_WP*real(meshsize(3),WP)*dx
         call param_read('Partition',partition)
         ! Allocate and initialize the shock-slab solver
         call ss%initialize(dx=dx,meshsize=meshsize,startloc=X0,group=group,partition=partition,continue_monitor=.false.)
         ! Provide relaxation and thermodynamic models
         ss%fs%relax=>P_relax_implicit
         ss%fs%getPL=>get_PL; ss%fs%getCL=>get_CL; ss%fs%getSL=>get_SL; ss%fs%getTL=>get_TL
         ss%fs%getPG=>get_PG; ss%fs%getCG=>get_CG; ss%fs%getSG=>get_SG; ss%fs%getTG=>get_TG
         ! Set viscosities to zero (inviscid)
         ss%cst_viscL=0.0_WP; ss%cst_viscG=0.0_WP
      end block setup_ss

      ! Generate initial conditions for shock-slab problem
      initialize_ss: block
         use irl_fortran_interface, only: setNumberOfPlanes,setPlane
         integer :: i,j,k
         real(WP) :: xloc
         ! Initialize primary variables
         do k=ss%cfg%kmino_,ss%cfg%kmaxo_; do j=ss%cfg%jmino_,ss%cfg%jmaxo_; do i=ss%cfg%imino_,ss%cfg%imaxo_
                  xloc=ss%fs%cfg%xm(i)
                  ! Initialize VOF for a liquid slab between slab_left and slab_right
                  if (xloc.ge.slab_left.and.xloc.le.slab_right) then
                     ss%fs%VF(i,j,k)=1.0_WP
                  else
                     ss%fs%VF(i,j,k)=0.0_WP
                  end if
                  ! Set volume moment barycenters
                  ss%fs%BL(:,i,j,k)=[ss%fs%cfg%xm(i),ss%fs%cfg%ym(j),ss%fs%cfg%zm(k)]
                  ss%fs%BG(:,i,j,k)=[ss%fs%cfg%xm(i),ss%fs%cfg%ym(j),ss%fs%cfg%zm(k)]
                  ! Set PLIC interface
                  call setNumberOfPlanes(ss%fs%PLIC(i,j,k),1)
                  if (abs(xloc-slab_left).lt.0.5_WP*ss%fs%dx) then
                     ! Left interface of slab: normal pointing left (gas on left)
                     call setPlane(ss%fs%PLIC(i,j,k),0,[+1.0_WP,0.0_WP,0.0_WP],slab_left)
                  else if (abs(xloc-slab_right).lt.0.5_WP*ss%fs%dx) then
                     ! Right interface of slab: normal pointing right (gas on right)
                     call setPlane(ss%fs%PLIC(i,j,k),0,[-1.0_WP,0.0_WP,0.0_WP],-slab_right)
                  else
                     ! Pure gas or pure liquid cell
                     call setPlane(ss%fs%PLIC(i,j,k),0,[0.0_WP,0.0_WP,0.0_WP],sign(1.0_WP,ss%fs%VF(i,j,k)-0.5_WP))
                  end if
                  ! Initialize mixture velocity to normal shock
                  ss%fs%U(i,j,k)=u2*Hshock(Xs-ss%fs%cfg%x(i),delta=0.5_WP*ss%fs%dx)
                  ss%fs%V(i,j,k)=0.0_WP
                  ss%fs%W(i,j,k)=0.0_WP
                  ! Gas variables
                  if (ss%fs%VF(i,j,k).lt.1.0_WP) then
                     ss%fs%RHOG(i,j,k)=rho1+(rho2-rho1)*Hshock(Xs-ss%fs%cfg%xm(i),delta=0.5_WP*ss%fs%dx)
                     ss%fs%PG  (i,j,k)=p1  +(p2  -p1  )*Hshock(Xs-ss%fs%cfg%xm(i),delta=0.5_WP*ss%fs%dx)
                     ss%fs%IG  (i,j,k)=(ss%fs%PG(i,j,k)+GammaG*PinfG)/(ss%fs%RHOG(i,j,k)*(GammaG-1.0_WP))
                  end if
                  ! Liquid variables
                  if (ss%fs%VF(i,j,k).gt.0.0_WP) then
                     ss%fs%RHOL(i,j,k)=rhoL
                     ss%fs%PL  (i,j,k)=p1
                     ss%fs%IL  (i,j,k)=(ss%fs%PL(i,j,k)+GammaL*PinfL)/(ss%fs%RHOL(i,j,k)*(GammaL-1.0_WP))
                  end if
               end do; end do; end do
         ! Build PLIC interface
         call ss%fs%build_interface()
         ! Initialize conserved variables
         ss%fs%Q(:,:,:,1)=        ss%fs%VF *ss%fs%RHOL
         ss%fs%Q(:,:,:,2)=(1.0_WP-ss%fs%VF)*ss%fs%RHOG
         ss%fs%Q(:,:,:,3)= ss%fs%Q(:,:,:,1)*ss%fs%IL
         ss%fs%Q(:,:,:,4)= ss%fs%Q(:,:,:,2)*ss%fs%IG
         call ss%fs%get_momentum()
         ! Communicate conserved variables
         do i=1,ss%fs%nQ; call ss%fs%cfg%sync(ss%fs%Q(:,:,:,i)); end do
         ! Rebuild primitive variables
         call ss%fs%get_primitive()
         ! Interpolate velocity
         call ss%fs%interp_vel(ss%Ui,ss%Vi,ss%Wi)
         ! Compute local Mach number
         ss%Ma=sqrt(ss%Ui**2+ss%Vi**2+ss%Wi**2)/ss%fs%C
         ! Perform monitoring
         call ss%output_monitor()
      end block initialize_ss

      ! Initialize Ensight output event and perform initial Ensight output
      initialize_ensight: block
         use param, only: param_read
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         if (ens_evt%occurs()) then
            call ss%output_ensight(t=time%t)
         end if
      end block initialize_ensight

      ! Initialize timers
      initialize_timers: block
         use parallel, only: comm,amRoot
         ! Create timers
         tstep=timer(comm=comm,name='Timestep')
         ! Create corresponding monitor file
         timefile=monitor(amRoot,'timing')
         call timefile%add_column(time%n,'Timestep number')
         call timefile%add_column(time%t,'Time')
         call timefile%add_column(tstep%time,trim(tstep%name))
      end block initialize_timers

   end subroutine simulation_init


   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none

      ! Overall time integration
      do while (.not.time%done())

         ! Reset timer and start timestep timer
         call tstep%reset()
         call tstep%start()

         ! Adjust time step size using CFL info
         call ss%fs%get_cfl(dt=time%dt,cfl=time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Advance shock-slab simulation
         call ss%step(dt=time%dt)

         ! Perform monitoring
         call ss%output_monitor()

         ! Perform Ensight output
         if (ens_evt%occurs()) then
            call ss%output_ensight(t=time%t)
         end if

         ! Stop timestep timer
         call tstep%stop()

         ! Output timing info
         call timefile%write()

      end do

   end subroutine simulation_run


   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      call ss%finalize()
      call ens_evt%finalize()
      call timefile%finalize()
      call tstep%finalize()
   end subroutine simulation_final


end module simulation
