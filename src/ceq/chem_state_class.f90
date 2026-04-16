!> Chemical state class
module chem_state_class
   use precision,      only: WP
   use chem_sys_class, only: chem_sys,Gphase,ncof
   implicit none
   private
   
   !> Expose type/constructor/methods
   public :: chem_state

   !> List of available equilibrium conditions
   integer, parameter, public :: fixed_PT=1                        !< Fixed pressure and temperature
   integer, parameter, public :: fixed_PH=2                        !< Fixed pressure and enthalpy
   integer, parameter, public :: fixed_UV=3                        !< Fixed internal energy and volume

   !> List of available fixed_PH algorithms
   integer, parameter, public :: NR=1                              !< Newton-Raphson
   integer, parameter, public :: BS=2                              !< Bi-section

   !> List of available methods for dN/dT calculation
   integer, parameter, public :: FD=1                              !< Finite Difference
   integer, parameter, public :: LS=2                              !< Least Squares

   !> Fraction of Nm used in initial guess
   real(WP), parameter :: frac_Nm=0.1_WP

   !> Temperature bounds
   real(WP), parameter :: T_low =250.0_WP
   real(WP), parameter :: T_high=5000.0_WP

   !> Reference pressure (Pa) and the universal gas constant
   real(WP), parameter :: p0=101325_WP
   real(WP), parameter :: gas_cnst=8.31446261815324_WP

   !> Chemical statete object definition
   type :: chem_state

      ! This is our chemical system
      class(chem_sys), pointer :: sys                              !< This is the chemical system the solver is build for

      ! Equilibrium condition
      integer :: cond

      ! Thermochemical quantities
      real(WP) :: p                                                !< Pressure (Pa)
      real(WP) :: T                                                !< Temperature (K)
      real(WP) :: HoR,HoR0                                         !< Enthalpy over ideal gas constant R
      real(WP) :: UoR,UoR0                                         !< Internal energy over ideal gas constant R
      real(WP) :: V,V0                                             !< Volume of the mixture
      real(WP), dimension(:),   allocatable :: N                   !< Moles of species (ns). Follows the same order of initial moles
      real(WP), dimension(:),   allocatable :: Ndu                 !< Moles of species (ns) reordered into ceq format: [Nd, Nu]
      real(WP), dimension(:),   allocatable :: Nbar                !< Moles of phases (np)
      real(WP), dimension(:),   allocatable :: Nd                  !< Moles of determined species (nsd)
      real(WP), dimension(:),   allocatable :: Nu,Nuold,Nm         !< Moles of undetermined species (nsu)
      real(WP), dimension(:),   allocatable :: gu                  !< Gibbs functions of undetermined species (nsu)
      real(WP), dimension(:),   allocatable :: vmolar              !< Molar volumes (ns) reordered into ceq format: [Nd, Nu]
      real(WP), dimension(:),   allocatable :: lam                 !< Lagrange multipliers (nrc)
      real(WP), dimension(:),   allocatable :: cr                  !< Reduced constraint vector
      real(WP), dimension(:),   allocatable :: RC,RCd              !< Constraints residual arrays
      real(WP), dimension(:),   allocatable :: sol                 !< Chemical state solution vector: [lambdas, ln(Nbars)]
      real(WP), dimension(:,:), allocatable :: Btilde,Ptilde       !< Coefficient matrices
      real(WP), dimension(:,:), allocatable :: BtildeT,PtildeT     !< Coefficient matrices transposed

      ! Pointer to the avaiable chemical equilibrium procedures
      procedure(get_ceq_interface), pointer :: get_ceq=>NULL()     !< Get the chemical equilibrium

      ! Pointer to the avaiable dN/dpar calculation procedures
      procedure(get_dNdpar_interface), pointer :: get_dNdT=>NULL() !< Get dN/dT
      procedure(get_dNdpar_interface), pointer :: get_dNdp=>NULL() !< Get dN/dp

      ! Numerical parameters
      real(WP) :: Tlo                                              !< Lowest temperature at which h has been evaluated
      real(WP) :: Thi                                              !< Highest temperature at which h has been evaluated
      real(WP) :: tol_N                                            !< Tolerance for the residual norm
      real(WP) :: tol_T                                            !< Tolerance for the temperature
      real(WP) :: tol_p                                            !< Tolerance for the pressure
      real(WP) :: tol_H                                            !< Tolerance for the enthalpy residual
      real(WP) :: tol_U                                            !< Tolerance for the internal-energy residual
      real(WP) :: tol_V                                            !< Tolerance for the volume residual
      real(WP) :: dT                                               !< Residual error for the temperature
      real(WP) :: dp                                               !< Residual error for the pressure
      real(WP) :: RH                                               !< Residual error for the enthalpy
      real(WP) :: RU,RV                                            !< Residual errors for the internal energy and the volume
      integer  :: iter_N                                           !< Number of Newton-Raphson iterations
      integer  :: iter_T                                           !< Number of temperature iterations
      integer  :: iter_N_max                                       !< Maximum number of Newton-Raphson iterations
      integer  :: iter_T_max                                       !< Maximum number of temperature iterations
      integer  :: PH_method                                        !< Fixed PH algorithm
      integer  :: dNdT_method                                      !< dNdT calculation method
      integer  :: dNdp_method                                      !< dNdp calculation method
      logical  :: success                                          !< Flag for successful equilibrium calculations

   contains

      procedure :: initialize                                      !< Object initializer

      procedure :: N_init                                          !< Initialize the number of moles of species
      procedure :: N_re_init                                       !< Re-initialize the number of moles of species

      procedure :: get_hort                                        !< Get the normalized enthalpy
      procedure :: get_cpor                                        !< Get the normalized Cp
      procedure :: get_phasic_HoR                                  !< Get the phasic H/R
      procedure :: get_gort                                        !< Get the normalized Gibbs free energy
      procedure :: hor2T                                           !< Convert enthalpy to temperature
      procedure :: get_dgdT                                        !< Get the temperature derivative of the gibbs function
      procedure :: get_dgdp                                        !< Get the pressure derivative of the gibbs function

      procedure :: perturb                                         !< Perturb the chemical equilibrium problem
      procedure :: get_Nming                                       !< Get the composition that minimized G and satisfies the constraints
      procedure :: min_pert                                        !< Get the purturbed maxmin composition
      procedure :: maxmin_comp                                     !< Get the minmax composition
      procedure :: solve_linprog                                   !< Solve the linear programming problem

      procedure :: get_Nusqrt                                      !< Get the square root of the mole numbers
      procedure :: get_RC                                          !< Get the constraints residual vector
      procedure :: get_RH                                          !< Get the enthalpy residual
      procedure :: get_RUV                                         !< Get the internal energy and volume residual
      procedure :: sol_init                                        !< Initialize the chemical state solution vector
      procedure :: equilibrate                                     !< Obtain the chemical equilibrium state of the system
      procedure :: get_Cp_eff                                      !< Get the effective Cp (dH/dT)_p,c
      procedure :: get_dsoldpar                                    !< Get the derivative of the solution vector with respect to given parameter
      procedure :: get_BP                                          !< Get the coefficient matrices for constraints and phase summation

      procedure, private :: get_ceq_PT                             !< Get the chemical equilibrium state at constant pressure and temperature
      procedure, private :: get_ceq_PH_NR,get_ceq_PH_BS            !< Get the chemical equilibrium state at constant pressure and emthalpy
      procedure, private :: get_ceq_UV                             !< Get the chemical equilibrium state at constant internal energy and volume
      procedure, private :: get_dNdT_FD,get_dNdT_LS                !< Get the temperature derivative of the mole numbers
      procedure, private :: get_dNdp_FD,get_dNdp_LS                !< Get the pressure derivative of the mole numbers
      procedure, private :: get_dNdpar_FD,get_dNdpar_LS            !< Get the derivative of the mole numbers with resprect to generic variable par

   end type chem_state

   !> Interface for get_ceq
   interface
      subroutine get_ceq_interface(this)
         use precision, only: WP
         import chem_state
         class(chem_state), intent(inout) :: this
      end subroutine get_ceq_interface
   end interface

   !> Interface for dN/dpar evaluation procedures
   interface
      subroutine get_dNdpar_interface(this,dNdpar)
         use precision, only: WP
         import chem_state
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(out) :: dNdpar
      end subroutine get_dNdpar_interface
   end interface


   contains


      !> Chemical state initializer
      subroutine initialize(this,sys,cond,PH_method,dNdT_method,dNdp_method,p,vmolar)
         use messager,  only: die
         use mathtools, only: reorder_rows
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         class(chem_sys), target, intent(in) :: sys
         integer,  intent(in) :: cond
         integer,  intent(in), optional :: PH_method,dNdT_method,dNdp_method
         real(WP), intent(in) :: p
         real(WP), intent(in), optional :: vmolar(sys%ns)
         integer  :: np,nb,nc,ns,nsd,nsu,nrc

         ! Point to chemical system
         this%sys=>sys

         ! Set the eqiuilibrium condition
         select case (cond)
            ! Constant pressure and temperature
            case (fixed_PT)
               this%get_ceq=>get_ceq_PT
            ! Constant pressure and enthalpy
            case (fixed_PH)
               ! Select the solver
               if (present(PH_method)) then
                  select case (PH_method)
                     case (NR)
                        this%get_ceq=>get_ceq_PH_NR
                     case (BS)
                        this%get_ceq=>get_ceq_PH_BS
                     case default
                        call die('[chem_state initialize] Unknown fixed p and H algorithm')
                  end select
                  this%PH_method=PH_method
               else
                  call die('[chem_state initialize] Fixed p and H algorithm requires an input for the algorithm')
               end if
               ! Select dN/dT method
               if (present(dNdT_method)) then
                  select case (dNdT_method)
                  case (FD)
                     this%get_dNdT=>get_dNdT_FD
                  case (LS)
                     this%get_dNdT=>get_dNdT_LS
                  case default
                     call die('[chem_state initialize] Unknown dN/dT evaluation method')
                  end select
                  this%dNdT_method=dNdT_method
               else
                  call die('[chem_state initialize] Fixed p and H algorithm requires an input for the calculation method of dN/dT')
               end if
            ! Constant internal energy and volume
            case (fixed_UV)
               this%get_ceq=>get_ceq_UV
               ! Select dN/dT method
               if (present(dNdT_method)) then
                  select case (dNdT_method)
                  case (FD)
                     this%get_dNdT=>get_dNdT_FD
                  case (LS)
                     this%get_dNdT=>get_dNdT_LS
                  case default
                     call die('[chem_state initialize] Unknown dN/dT evaluation method')
                  end select
                  this%dNdT_method=dNdT_method
               else
                  call die('[chem_state initialize] Fixed U and V algorithm requires an input for the calculation method of dN/dT')
               end if
               ! Select dN/dp method
               if (present(dNdp_method)) then
                  select case (dNdp_method)
                  case (FD)
                     this%get_dNdp=>get_dNdp_FD
                  case (LS)
                     this%get_dNdp=>get_dNdp_LS
                  case default
                     call die('[chem_state initialize] Unknown dN/dp evaluation method')
                  end select
               else
                  call die('[chem_state initialize] Fixed U and V algorithm requires an input for the calculation method of dN/dp')
               end if
               this%dNdp_method=dNdp_method
            case default
               call die('[chem_state initialize] The chemical state must be at either constant temperature and pressure, constant enthalpy and pressure, or constant internal energy and volume')
         end select
         this%cond=cond

         ! Determine pressure
         if (p.le.0.0_WP) call die('[chem_state initialize] Pressure must be strictly positive')
         this%p=p

         ! Obtain indexes
         np =sys%np
         nb =sys%nb
         nc =sys%nc
         nrc=sys%nrc
         ns =sys%ns
         nsd=sys%nsd
         nsu=sys%nsu

         ! Allocate arrays
         allocate(this%N     (ns));         this%N      =0.0_WP
         allocate(this%Ndu   (ns));         this%Ndu    =0.0_WP
         allocate(this%Nbar  (np));         this%Nbar   =0.0_WP
         allocate(this%Nd    (nsd));        this%Nd     =0.0_WP
         allocate(this%Nu    (nsu));        this%Nu     =0.0_WP
         allocate(this%Nuold (nsu));        this%Nuold  =0.0_WP
         allocate(this%Nm    (nsu));        this%Nm     =0.0_WP
         allocate(this%gu    (nsu));        this%gu     =0.0_WP
         allocate(this%vmolar(ns));         this%vmolar =0.0_WP
         allocate(this%lam   (nrc));        this%lam    =0.0_WP
         allocate(this%cr    (nrc));        this%cr     =0.0_WP
         allocate(this%RC    (nrc+np));     this%RC     =0.0_WP
         allocate(this%RCd   (np));         this%RCd    =0.0_WP
         allocate(this%sol   (nrc+np));     this%sol    =0.0_WP
         allocate(this%Btilde(nsu,nrc));    this%Btilde =0.0_WP
         allocate(this%Ptilde(nsu,np));     this%Ptilde =0.0_WP
         allocate(this%BtildeT(nrc,nsu));   this%BtildeT=0.0_WP
         allocate(this%PtildeT(np,nsu));    this%PtildeT=0.0_WP

         ! Store and reorder the molar volumes
         if (present(vmolar)) call reorder_rows(vmolar,this%sys%sp_order,this%vmolar)

      end subroutine initialize


      !> Initialize the mole numbers
      subroutine N_init(this,T,c,N,HoR,UoR,V,N_h,T_h,N_g,T_g,p_g)

         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation. 

         ! Initializes the constrained equilibrium state of 
         ! an ideal liquid-gas mixture consisting of ns species,either at fixed pressure and
         ! temperature (p,T),or at fixed pressure and enthalpy (p,H).

         ! The nc equality constraints are written:  B'*N=c,where B is the
         ! ns x nc basic constraint matrix,N is the ns-vector of species moles,
         ! and c is the nc-vector of constraint values.

         ! Input:
         !  sys -object type chem_sys created by subroutine.
         !         Required first argument

         !  cond -Equilibrium condition: Either T_fixed or H_fixed

         !  (Either c or N must be specified.)
         !  c   -values of the nc basic constraints (real(nc))
         !  N   -moles of species used to calculate c as : c=B'*N (real(ns))

         !  p   -pressure

         ! (For a fixed (p,T) equilibrium calculation, specify T: do not specify HoR, N_h or T_h.)
         ! T-temperature (K) for fixed-temperature problem

         ! (For a fixed (p,H) equilibrium calculation, specify either HoR or N_h and T_h: 
         !  do not specify T.)
         ! HoR the fixed value of H/R [moles K], where H=enthalpy, R=universal gas constant.
         ! N_h species moles used to calculate H  (real(ns))
         ! T_h temperature used to calculate H as:  H/R=sum(N_h h(T_h)/R), where
         !     h(T_h)/R [which has dimensions K] is the molar specific species enthalpy.

         ! (For a fixed (p,H) equilibrium calculation, specify UoR and V)

         ! (Initial guesses are not needed, and should not be specified unless they
         !  are good guesses.)
         ! N_g  -initial guess for species moles
         ! T_g  -initial guess for temperature (for fixed (p,H) and fixed (U,V) only)
         ! p_g  -initial guess for pressure (for fixed (U,V) only)

         !  To diagnose an error condition,the case can be repeated with diagnostics turned on,
         !  by: call param_set(sys,diag=5).

         use, intrinsic :: iso_fortran_env, only: output_unit
         use messager,  only: die
         use mathtools, only: reorder_rows
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), intent(in), optional :: c(this%sys%nc),N(this%sys%ns),T,HoR,UoR,V,N_h(this%sys%ns),T_h,N_g(this%sys%ns),T_g,p_g
         integer  :: np,nb,nc,ns,nsd,nsu,nrc,npert,iret,i
         real(WP) :: max_pert,Numin,cb(this%sys%nc),cmod(this%sys%nb),Nd(this%sys%nsd),cr_norm,cb_norm,res,N_low,res_tol=1e-9
         real(WP), dimension(this%sys%ns)  :: N0,N1,h
         real(WP), dimension(this%sys%nsu) :: Nu,Nu0,Nm,Nupper,Ng,gu
         real(WP), dimension(this%sys%nrc) :: cr
         logical :: fail,diag,use_mmg=.true.
         this%success=.true.
         ! Get the neccesary inputs
         select case (this%cond)
            case (fixed_PT)
               if (present(T)) then
                  if (T.lt.T_low.or.T.gt.T_high) call die('[chem_state N_init] Temperature out of range')
                  this%T=T
               else
                  call die('[chem_state N_init] Temperature is required for the fixed temperature condition')
               end if
            case (fixed_PH)
               if (present(HoR)) then
                  this%HoR=HoR
               elseif(present(N_h).and.present(T_h)) then
                  call reorder_rows(N_h,this%sys%sp_order,N0)
                  call this%get_hort(this%sys%ns,T_h,this%sys%thermo,h)
                  this%HoR=sum(N0*h)*T_h
               else
                  call die('[chem_state N_init] Both N_h and T_h are required for the fixed enthalpy and pressure condition')
               end if
               if (present(T_g)) then
                  this%T=T_g
                  if (T_g.lt.T_low.or.T_g.gt.T_high) call die('[chem_state N_init] Guessed temperature out of range')
               else
                  this%T=sqrt(T_low*T_high)
                  this%T=max(this%T,0.1_WP*T_high)
               endif
            case (fixed_UV)
               if (present(UoR).and.present(V)) then
                  this%UoR=UoR
                  this%V=V
               else
                  call die('[chem_state N_init] Both UoR and V are required for the fixed internal energy and volume condition')
               end if
               if (present(T_g)) then
                  if (T_g.lt.T_low.or.T_g.gt.T_high) call die('[chem_state N_init] Guessed temperature out of range')
                  this%T=T_g
               else
                  this%T=sqrt(T_low*T_high)
                  this%T=max(this%T,0.1_WP*T_high)
               endif
               if (present(p_g)) then
                  if (p_g.le.0.0_WP) call die('[chem_state N_init] Guessed pressure must be strictly positive')
                  this%p=p_g
               else
                  this%p=p0
               endif
         end select

         ! Obtain indexes
         np =this%sys%np
         nb =this%sys%nb
         nc =this%sys%nc
         nrc=this%sys%nrc
         ns =this%sys%ns
         nsd=this%sys%nsd
         nsu=this%sys%nsu

         ! Initialize the Gibbs function
         call this%get_gort(nsu,this%T,this%p,this%sys%thermo(nsd+1:ns,:),this%sys%P(nsd+1:ns,Gphase),gu)

         ! Form the basic constraint vector
         if (present(c)) then
            cb(1:nc)=c
         elseif(present(N)) then
            cb(1:nc)=matmul(N,this%sys%B)
         else
            call die('[chem_state initialize] Neither c nor N specified')
         endif

         ! Form modified and reduced constraints
         cmod(1:nb)=matmul(this%sys%A(1:nb,1:nc),cb)
         Nd(1:nsd) =cmod(1:nsd)

         ! Treat the special case of no undetermined species
         if (nsu.eq.0) then
            this%Ndu=Nd
            if (this%cond.eq.fixed_PH) call this%hor2T(ns,this%Ndu,this%HoR,this%sys%thermo,this%T)
         else

            ! Reduced constraint vector
            cr(1:nrc) =cmod(nsd+1:nb)
            cr_norm   =norm2(cr)

            if (cr_norm.le.0.0_WP) then
            ! SBP added 4/9/2009
               if (cr_norm.eq.0.0_WP.and.nsd.gt.0.and.sum(Nd(1:nsd)).gt.0.0_WP) then
                  !  only determined species
                  this%Ndu=0.0_WP
                  this%Ndu(1:nsd)=Nd(1:nsd)
                  if (this%cond.eq.fixed_PH) call this%hor2T(ns,this%Ndu,this%HoR,this%sys%thermo,this%T)
               endif
               ! SBP end of added
               this%success=.false.
               ! write(output_unit,'(" >   [chem_state initialize] All zero composition")')
               return
               ! call die('[chem_state initialize] All zero composition')
            endif

            ! Use initial guess N_g if provided
            if (present(N_g)) then
               use_mmg=.false.
               call reorder_rows(N_g,this%sys%sp_order,N0)
               ! Guessed undetermined species
               Nu0(1:nsu)=N0(nsd+1:ns)
               ! Reduced c.v. based on N_g
               cb(1:nrc)=matmul(Nu0(1:nsu),this%sys%BR)
               cb_norm  =norm2(cb)
               if (cb_norm.eq.0.0_WP) then
                  use_mmg=.true.
               else
                  res=norm2((cb(1:nrc)/cb_norm-cr/cr_norm))
                  ! Adjust initial guess Nu0,store in Nm
                  if (res.gt.res_tol) then
                     N_low=sum(cr(1:this%sys%neu))*1e-15
                     call this%min_pert(nsu,nrc,this%sys%BR,cr,Nu0,N_low,Nm,iret)
                     ! min_pert failed
                     if (iret.ne.0) then
                        print*,'Warning: min_pert failed. Could not use N_g'
                        use_mmg=.true.
                     else
                        Nu0=Nm
                     endif
                  endif
               end if
            endif

            ! Use max-min and min_g if needed
            if (use_mmg) then
               ! Perturb if necessary
               call this%perturb(ns,nsd,nsu,this%sys%ne,this%sys%ned,this%sys%neu,nrc,Nd,cr,this%sys%BR,this%sys%E,this%sys%diag,this%sys%eps_el,this%sys%eps_sp, &
               &                 this%sys%pert_tol,this%sys%pert_skip,this%Nd,Nm,Nupper,this%cr,npert,max_pert,iret)
               this%Nm=Nm
               if (iret.eq.-1) then 
                  write(output_unit,'(" >   chem_state perturb: non-realizable constraint = ")')
               elseif(iret.eq.-2) then
                  this%success=.false.
                  write(output_unit,'(" >   [chem_state initialize] Perturb failed")')
                  return
               endif
               ! Determine min_g composition
               call this%get_Nming(nsu,nrc,this%sys%BR,this%cr,gu,Ng,iret)
               if (iret.lt.0) then
                  this%success=.false.
                  return
               end if
               ! Form initial guess Nu0
               Nu0=Ng+frac_Nm*(Nm-Ng)
            else
               this%Nd=Nd
               this%cr=cr
            end if

            ! Re-estimate T0 and re-evaluate gu if required
            if ((this%cond.eq.fixed_PH).and.(.not.present(T_g))) then
               N1(1:nsd)   =this%Nd
               N1(nsd+1:ns)=Nu0
               call this%hor2T(ns,N1,this%HoR,this%sys%thermo,this%T)
               ! Set gu based on T0
               call this%get_gort(nsu,this%T,this%p,this%sys%thermo(nsd+1:ns,:),this%sys%P(nsd+1:ns,Gphase),gu)
            endif

            ! Set the Gibbs functin and the undetermined species moles
            this%gu=gu
            this%Nu=Nu0

            ! Determine required output
            this%Ndu=[this%Nd,this%Nu]

         endif

         ! Update thermodynamic quantities
         if (this%cond.eq.fixed_PT) then
            call this%get_hort(ns,this%T,this%sys%thermo,h)
            this%HoR=sum(this%Ndu*h)*this%T
         elseif (this%cond.eq.fixed_UV) then
            this%UoR0=this%UoR
            this%V0=this%V
         endif

         ! Calculate the phase moles
         this%Nbar=matmul(transpose(this%sys%P),this%Ndu)

         ! Re-order species
         do i=1,ns
            this%N(this%sys%sp_order(i))=this%Ndu(i)
         end do

      end subroutine N_init


      !> Re-initialize mole numbers using min-g composition
      subroutine N_re_init(this)
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu) :: Ng
         integer  :: ns,nsd,nsu,nrc,i,iret
         ! Obtain indexes
         ns =this%sys%ns
         nsd=this%sys%nsd
         nsu=this%sys%nsu
         nrc=this%sys%nrc
         ! Get the gibbs of undetermined species
         call this%get_gort(nsu,this%T,this%p,this%sys%thermo(nsd+1:ns,:),this%sys%P(nsd+1:ns,Gphase),this%gu)
         ! Determine min_g composition
         call this%get_Nming(nsu,nrc,this%sys%BR,this%cr,this%gu,Ng,iret)
         this%success=.true.
         if (iret.lt.0) then
            this%success=.false.
            write(output_unit,'(" >   [chem_state N_re_init] get_Nming failed")')
            return
         end if
         ! Update the moles
         this%Nu=Ng+frac_Nm*(this%Nm-Ng)
         this%Ndu=[this%Nd,this%Nu]
         this%Nbar=matmul(transpose(this%sys%P),this%Ndu)
         ! Re-order species
         do i=1,ns
            this%N(this%sys%sp_order(i))=this%Ndu(i)
         end do
         ! Update the gibbs of undetermined species
         call this%get_gort(nsu,this%T,this%p,this%sys%thermo(nsd+1:ns,:),this%sys%P(nsd+1:ns,Gphase),this%gu)
      end subroutine N_re_init


      !> Get normalized enthalpies (Neglecting pressure dependence for liquid)
      subroutine get_hort(this,ns,T,thermo,hort)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in)   :: ns
         real(WP), intent(in)  :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: hort(ns)
         ! input:
         !	ns	     -number of species
         !  T       -temperature (K)
         !  thermo  -thermo data for all species
         ! output:
         !  hort    -h_j/(RT) -normalized enthalpies
         ! S. B. Pope 9/26/02
         real(WP) :: th(6),Tpnm1
         integer :: k,n
         th(1)=1.0_WP  ! coefficient multipliers for enthalpy
         th(6)=1./T
         Tpnm1=1.0_WP
         do n=2,5
            Tpnm1=Tpnm1*T            ! =T.^(n-1)
            th(n)=Tpnm1/float(n)     ! =T.^(n-1) ./ n
         end do
         do k=1,ns
            if (T<thermo(k,1)) then
               hort(k)=dot_product(thermo(k,2:7),th)  ! coefficients in lower temperature range
            else
               hort(k)=dot_product(thermo(k,9:14),th) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_hort


      !> Get the phasic enthalpy (H/R) (Neglecting pressure dependence for liquid)
      subroutine get_phasic_HoR(this,phase,N,T,HoR)
         use mathtools, only: reorder_rows
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in)  :: phase ! Follows IRL convention; 0 if liquid, 1 if gas
         real(WP), intent(in) :: T,N(this%sys%ns)
         real(WP), intent(out) :: HoR
         real(WP) :: th(6),Tpnm1,Nro(this%sys%ns)
         integer :: k,m
         ! Reorder mole numbers
         call reorder_rows(N,this%sys%sp_order,Nro)
         ! Initialize with zero
         HoR=0.0_WP
         ! Use NASA polynomials
         th(1)=1.0_WP  ! coefficient multipliers for enthalpy
         th(6)=1./T
         Tpnm1=1.0_WP
         do m=2,5
            Tpnm1=Tpnm1*T            ! =T.^(m-1)
            th(m)=Tpnm1/float(m)     ! =T.^(m-1) ./ m
         end do
         do k=1,this%sys%ns
            if (int(this%sys%P(k,Gphase)).eq.phase) then
               if (T<this%sys%thermo(k,1)) then
                  HoR=HoR+Nro(k)*dot_product(this%sys%thermo(k,2:7),th)  ! coefficients in lower temperature range
               else
                  HoR=HoR+Nro(k)*dot_product(this%sys%thermo(k,9:14),th) ! coefficients in upper temperature range
               endif
            end if
         end do
         HoR=T*HoR
      end subroutine get_phasic_HoR


      !> Get normalized Gibbs functions (Neglecting pressure dependence for liquid enthalpy)
      subroutine get_gort(this,ns,T,p,thermo,isGas,gort)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(WP), intent(in) :: T,p,thermo(ns,2*ncof+1),isGas(ns)
         real(WP), intent(out) :: gort(ns)
         ! input:
         !   ns     -number of species
         !   T      -temperature (K)
         !   p      -pressure (Pa)
         !   thermo -thermo data for all species
         !   isGas  -1 if gas,0 if liquid
         ! output:
         !   gort  -g_j/(RT) -normalized Gibbs functions
         ! S. B. Pope 9/26/02
         real(WP) :: tc(ncof),th(ncof),ts(ncof),tg(ncof)
         integer :: k,n
         if (ns.le.0) return
         tc=0.0_WP  ! coefficient multipliers for specific heats
         th=0.0_WP  ! coefficient multipliers for enthalpy
         ts=0.0_WP  ! coefficient multipliers for entropy
         tc(1)=1.0_WP
         th(1)=1.0_WP
         th(6)=1./T
         ts(1)=log(T)
         ts(7)=1.0_WP
         do n=2,5
            tc(n)=T*tc(n-1)   ! =T.^(n-1)
            th(n)=tc(n)/float(n)     ! =T.^(n-1) ./ n
            ts(n)=tc(n)/float((n-1)) ! =T.^(n-1) ./ (n-1)
         end do
         tg=th-ts
         do k=1,ns
            if (T<thermo(k,1)) then
               gort(k)=dot_product(thermo(k,2:8),tg)  ! coefficients in lower temperature range
            else
               gort(k)=dot_product(thermo(k,9:15),tg) ! coefficients in upper temperature range
            endif
         end do
         gort=gort+isGas*log(p/p0)
      end subroutine get_gort


      !> Determine temperature given enthalpy
      subroutine hor2T(this,ns,z,hin,thermo,T)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use messager, only: die
         implicit none
         class(chem_state),  intent(in)  :: this
         integer,            intent(in)  :: ns
         real(WP), intent(in)  :: z(ns),hin,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: T
         ! input:
         !	ns		  - number of species
         !   z      -moles of species
         !   hin    -enthalpy/R (K)=z'*h
         !   thermo -thermo data
         ! output:
         !   T      -temperature (K)

         ! Notes:  if the temperature is outside the range [T_low T_high]
         !   then T is returned as the closest of these bounds.
         !   If iteration fails,T is returned as T=-1.

         ! S. B. Pope 9/26/02

         integer :: itmax,it
         real(WP) :: T_tol,T0,hort(ns),hor,h_a,T_a,h_b,T_b,dT,&
            cpor(ns),hh,cpp

         itmax=100     ! maximum number of Newton iterations (usually only 3 required)
         T_tol=1e-6    ! error tolerance
         T0=1500.0_WP  ! initial guess

         !  determine if T>T0 and bracket T in [T_a T_b]
         call this%get_hort(ns,T0,thermo,hort)
         hor=dot_product(z,hort)*T0

         if (hin>hor) then	! T > T0=T_a
            h_a=hor
            T_a=T0
            call this%get_hort(ns,T_high,thermo,hort)
            h_b=dot_product(z,hort)*T_high
            if (hin.ge.h_b) then
               T=T_high   ! T > T_high (return T=T_high)
               return
            endif
            T_b=T_high
         else
            h_b=hor	! T < T0=T_b
            T_b=T0
            call this%get_hort(ns,T_low,thermo,hort)
            h_a=dot_product(z,hort)*T_low
            if (hin.le.h_a) then
               T=T_low    ! T < T_low (return T=T_low)
               return
            endif
            T_a=T_low
         endif

         !  estimate of T based on linear interpolation
         T=T_a+(hin-h_a)*(T_b-T_a)/(h_b-h_a)
            
         !  Newton iterations
         do it=1,itmax
            call this%get_cpor(ns,T,thermo,cpor)
            call this%get_hort(ns,T,thermo,hort)
            hh=dot_product(z,hort)*T
            cpp=dot_product(z,cpor)
            dT=(hin-hh)/cpp
            T=T+dT
            if (abs(dT).lt.T_tol) return  ! success
         end do

         ! Failure
         call die('[chem_state hor2T] Iterations failed')

      end subroutine hor2T


      !> Return d/dT of the normalized Gibbs functions
      subroutine get_dgdT(this,ns,T,thermo,dgdT)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer,  intent(in)  :: ns
         real(WP), intent(in)  :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: dgdT(ns)
         ! input:
         !   T      - temperature (K)
         !   thermo - thermo data for all species
         ! output:
         !   dgdT   - d/dT (g_j/(RT))
         ! S. B. Pope 7/1/03
         real(WP) :: tc(ncof),th(ncof),ts(ncof),tg(ncof)
         integer  :: k,n
         tc=0.d0  ! coefficient multipliers for specific heats
         th=0.d0  ! coefficient multipliers for enthalpy
         ts=0.d0  ! coefficient multipliers for entropy
         tc(1)=1.d0/T
         th(1)=0.d0
         th(6)=-1.d0/T**2
         ts(1)=1.d0/T
         ts(7)=0.d0
         do n=2,5
            tc(n)=T*tc(n-1)                  ! =T.^(n-2)
            th(n)=tc(n)*float(n-1)/float(n)  ! =T.^(n-2) * (n-1) / n
            ts(n)=tc(n)                      ! =T.^(n-2) 
         end do
         tg=th-ts
         do k=1,ns
            if (T<thermo(k,1)) then
               dgdT(k)=dot_product(thermo(k,2:8),tg)  ! coefficients in lower temperature range
            else
               dgdT(k)=dot_product(thermo(k,9:15),tg) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_dgdT


      !> Return d/dp of the normalized Gibbs functions (Neglecting pressure dependence for liquid)
      subroutine get_dgdp(this,ns,T,p,isGas,vmolar,dgdp)
         implicit none
         class(chem_state), intent(in) :: this
         integer,  intent(in)  :: ns
         real(WP), intent(in)  :: T,p,isGas(ns),vmolar(ns)
         real(WP), intent(out) :: dgdp(ns)
         ! input:
         !   T      - temperature (K)
         !   p      - pressure (Pa)
         !   isGas  -1 if gas,0 if liquid
         !   vmolar - molar volume (m^3/mol)
         ! output:
         !   dgdp   - d/dp (g_j/(RT))
         dgdp=isGas/p
         ! dgdp=isGas/p+(1.0_WP-isGas)*vmolar/(gas_cnst*T)
      end subroutine get_dgdp


      !> Generate (possibly) perturbed CE problem
      subroutine perturb(this,ns,nsd,nsu,ne,ned,neu,nrc,Nd,cr,BR,E,ifop,eps_el,eps_sp,pert_tol,pert_skip,zdp,zup,Nupper,crp,npert,max_pert,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns,nsd,nsu,ne,ned,neu,nrc,ifop,pert_skip
         real(WP), intent(in) :: Nd(nsd),cr(nrc),BR(nsu,nrc),E(ns,ne),eps_el,eps_sp,pert_tol
         integer, intent(out) :: npert,iret
         real(WP), intent(out) :: zdp(nsd),zup(nsu),Nupper(nsu),crp(nrc),max_pert

         !  Input:
         !	ns		   - number of species
         !	nsd		- number of determined species
         !	nsu		- number of undetermined species
         !	ne		   - number of elements
         !	ned		- number of determined elements
         !	neu		- number of undetermined elements
         !	nrc		- number of reduced constraints
         !  Nd       -moles of determined species
         !  cr       -reduced constraint vector
         !  BR       -reduced constraint matrix
         !  E        -element matrix
         !  ifop     >0 for output
         !  eps_el   -relative lower bound on element moles
         !  eps_sp   -relative lower bound on species moles
         !  pert_tol - largest allowed perturbation (moles/moles of atoms)
         !  pert_skip>0 to skip perturbing undetermined species

         !  Output:
         !   zdp    -perturbed moles of determined species
         !   zup    -min-max solution for undetermined species
         !   Nupper -upper bound on undetermined species moles
         !   crp    -perturbed reduced constraint vector
         !   npert  -number of perturbations made
         !   max_pert- largest normalized perturbation made
         !   iret   = 0  successful operation
         !          =-1,non-realizable large perturbation made
         !          =-2,failed to determine max-min composition
         !          =-3,zero atoms

         integer :: j,k
         real(WP) :: zdlim,zatoms,&
         cre(neu),sumcre,cref,zed(ned),zeu(neu),ze(ne),zeu_in(neu),zau,zelow,&
         zemax,zumm(nsu),Numin,zulow,zeumax

         iret=0         ! anticipate success
         npert=0        ! number of perturbations
         max_pert=0.0_WP  ! largest normalized perturbation

         zdp=Nd       ! check that determined species are non-negative
         zatoms=0.0_WP  ! estimate of moles of atoms
         if (nsd.gt.0) then
         do j=1,ne
            zatoms=zatoms+abs(dot_product(Nd,E(1:nsd,j)))
         end do
         endif

         do j=1,neu
         zatoms=zatoms+abs(cr(j))
         end do

         if (zatoms.le.0.0_WP) then
            ! write(output_unit,'(" >   chem_state perturb: no atoms")')
            iret=-3
            return
         endif

         zdlim=zatoms*pert_tol
         do k=1,nsd
            if (zdp(k)<0.0_WP) then
               if (abs(zdp(k))>zdlim ) then ! significantly negative
                     max_pert=max(max_pert,abs(zdp(k))/zatoms)
                     npert=npert+1
                     ! if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: negative determined species")')
               endif
               zdp(k)=0.0_WP
            endif
         end do

         crp=cr
         cre=cr(1:neu) ! check that undetermined elements are positive
         sumcre=0.0_WP
         do j=1,neu
         sumcre=sumcre+abs(cre(j))
         end do
         cref=sumcre*pert_tol

         do k=1,neu
            if (cre(k)<0.0_WP) then
               if (abs(cre(k))>cref) then ! significantly negative
                     max_pert=max(max_pert,abs(cre(k))/zatoms)
                     npert=npert+1
                     ! if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: negative undetermined element")')
               endif
               cre(k)=0.0_WP
            endif
         end do
         crp(1:neu)=cre
                     
         zed=matmul(zdp,E(1:nsd,1:ned))                !  moles of determined elements
         zeu=matmul(zdp,E(1:nsd,ned+1:ne))+crp(1:neu)  ! moles of undetermined elements
         ze(1:ned)=zed
         ze(ned+1:ne)=zeu    ! moles of elements
         zeu_in=zeu
         zemax=maxval(ze)
         zeumax=maxval(zeu)
         zau=sum(zeu)       ! moles of atoms in undetermined species

         !  impose lower bound on moles of undetermined elements
         zelow=max(eps_el*zeumax,eps_el**2*zemax)
         do j=1,neu      
            if (zeu(j)<zelow) then
               npert=npert+1
               max_pert=max(max_pert,(zelow-zeu(j))/zatoms)
               zeu(j)=zelow
            endif
         end do

         Nupper=0.0_WP    ! determine upper bound on undetermined species
         do j=1,nsu
            Nupper(j)=1/maxval(E(nsd+j,ned+1:ne)/zeu)
         end do

         !  determine max-min moles of undetermined species

         call this%maxmin_comp(nsu,nrc,BR,crp,zumm,Numin,iret)

         if (iret<0) then
            ! if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: maxmin_comp failed")')
            iret=-2 
            return 
         endif

         zup=zumm
         do j=1,nsu     ! impose lower limit on undetermined species
            zulow=eps_sp*Nupper(j)
            if (zumm(j)<zulow) then
               zup(j)=zulow
               npert=npert+1
               max_pert=max(max_pert,(zup(j)-zumm(j))/zatoms)
            endif
         end do

         if (pert_skip.gt.0) return  ! do not modify constraints

         !  modify constraints according to the perturbation in undetermined species
         crp=crp+matmul(zup-zumm,BR)

         if (max_pert.gt.pert_tol) then
            ! if (ifop.ge.1) write(output_unit,'(" >   chem_state perturb: large perturbation made")')
            iret=-1 
         endif

      end subroutine perturb


      !> Get Ng,the value of N which minimizes g'N,subject to B'*N=c,N(i)>=0.
      subroutine get_Nming(this,nz,nc,B,c,g,Ng,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(WP), intent(in) :: B(nz,nc),c(nc),g(nz)
         integer, intent(out) :: iret
         real(WP), intent(out) :: Ng(nz)
         !  Exit flag from linprog: iret<0 for failure.
         !  S.B. Pope 10/1/02
         integer :: i,iftest  
         real(WP) :: BT(nc,nz)
         BT=transpose(B)
         call this%solve_linprog(nz,nc,g,BT,c,Ng,iret)
         do i=1,nz  ! guard against small negative values due to round-off
         Ng(i)=max(Ng(i),0.0_WP)
         end do
      end subroutine get_Nming


      !> Determine N=N0+d which satisfies:
      !>    1) N(i).ge.eps
      !>    2) B'*N=c
      !>    3) t=max_i(|d(i)|) is minimized.
      subroutine min_pert(this,nz,nc,B,c,N0,eps,z,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         ! Input:
         !  nz -length of N
         !  nc -length of c
         !  B  -nz x nc equality constraint matrix
         !  c  -constraint vector
         !  N0 -nz-vector,N0
         !  eps-positive threshold
         ! Output:
         !  z   -solution
         !  iret=0 for success
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in)           :: nz,nc
         real(WP), intent(in)  :: B(nz,nc),c(nc),N0(nz),eps
         real(WP), intent(out) :: z(nz)
         integer, intent(out)          :: iret
         real(WP) :: A(nc+nz,3*nz),x(3*nz),f(3*nz),r(nc+nz),d(nz),tp,tm,res,Bsc,csc,zsc
         integer :: i,nx,nr
         logical :: linear=.false.
         ! x=[ u v w ]=[ z-eps (t+d)/2   (t-d)/2  ]
         ! f=[ 0 1 1 1   [ 0 0 0     1 1 1    1 1 1   ] 
         ! minimize sum(f*x)=nz*t subject to  A x=r
         nx=3*nz
         nr=nc+nz
         Bsc=max(maxval(B),-minval(B))  !  scale factors
         csc=max(maxval(c),-minval(c))
         zsc=csc/Bsc
         A(1:nz+nc,1:3*nz)=0.0_WP
         if (linear) then  !  min. sum of dz
            do i=1,nz
               A(i,i)     = 1.0_WP
               A(i,nz+i)  =-1.0_WP
               A(i,2*nz+i)= 1.0_WP
            end do
         else  !  min. sum of dz/N0
            do i=1,nz
               A(i,i)     = 1.0_WP
               A(i,nz+i)  =-(max(N0(i),eps)/zsc)**0.0_WP
               A(i,2*nz+i)= (max(N0(i),eps)/zsc)**0.0_WP
            end do
         endif
         A(nz+1:nz+nc,1:nz)  =transpose(B)/Bsc
         r(1:nz)      =(N0(1:nz) -eps) /zsc
         do i=1,nc
            r(nz+i)=(c(i)-eps*sum(B(:,i)))/csc
         end do
         f(1:nz)     =0.0_WP
         do i=1,nz
            f(nz+i)  =(eps/(eps+N0(i)))**0.0_WP
            f(2*nz+i)=f(nz+i)
         end do
         call this%solve_linprog(nx,nr,f,A,r,x,iret)
         if (iret.eq.0) then
            do i=1,nz
               z(i)=max(x(i),0.0_WP)*zsc+eps
            end do
         endif
         if (.false.) return
      end subroutine min_pert


      !> Get the normalized Cp's at temperature T
      subroutine get_cpor(this,ns,T,thermo,cpor)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: ns
         real(WP), intent(in) :: T,thermo(ns,2*ncof+1)
         real(WP), intent(out) :: cpor(ns)
         ! input:
         !	ns	  -number of species
         !   T     -temperature (K)
         !   thermo-thermo data for all species
         ! output:
         !   cpor  -Cp_j/R   -normalized constant-pressure specific heats
         ! S. B. Pope 9/26/02
         real(WP) :: tc(5)
         integer :: k,n
         cpor=0.0_WP
         tc(1)=1.0_WP  ! coefficient multipliers for specific heats
         do n=2,5
            tc(n)=T*tc(n-1)   ! =T.^(n-1) 
         end do
         do k=1,ns
            if (T.lt.thermo(k,1)) then
               cpor(k)=dot_product(thermo(k,2:6),tc)  ! coefficients in lower temperature range
            else
               cpor(k)=dot_product(thermo(k,9:13),tc) ! coefficients in upper temperature range
            endif
         end do
      end subroutine get_cpor


      !> Determine the max-min composition.
      subroutine maxmin_comp(this,nz,nc,B,c,Nm,zmin,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nz,nc
         real(WP), intent(in) :: B(nz,nc),c(nc)
         integer, intent(out) :: iret
         real(WP), intent(out) :: Nm(nz),zmin
         !  Find Nm which maximizes zmin=min_i(z(i)) subject to B'*z=c.
         !  iret<0 indicates failure
         !  Method: 
         !   Initially assume zmin>=0.
         !   Define: x=[ (z-zmin)' zmin]'
         !   Maximize zmin (i.e.,minimize -x(n)) subject to x(i)>=0
         !       and B'*z=c,which is equivalent to A*x=c,
         !       where A=[ B' -sum(B)'].
         !   If feasible solution not found,zmin<0,and
         !   Define: x=[ (z-zmin)' -zmin]'
         !   Maximize zmin (i.e.,minimize x(n)) subject to x(i)>=0
         !       and B'*z=c,which is equivalent to A*x=c,
         !       where A=[ B' sum(B)'].
         !   S.B. Pope 10/1/02
         integer :: nx,tries, j,jj
         real(WP) :: A(nc,nz+1),bsum(nc),f(nz+1),x(nz+1)
         nx=nz+1
         bsum=sum(B,dim=1)
         f=0
         ! First assume zmin>0,x=[z'-zmin zmin]'
         f(nx)=-1.0_WP  ! minimize -zmin
         A(1:nc,1:nz)=transpose(B)
         A(1:nc,nx)=bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if (iret.eq.0) then  !  success,zmin>=0
            zmin=x(nx)
            Nm=x(1:nz)+zmin
            return
         elseif(iret.ne.-2) then  !  failure
            return
         endif
         ! zmin<0,re-define x=[z'-zmin -zmin]'
         f(nx)=1.0_WP  ! minimize -zmin
         A(1:nc,nx)=-bsum
         call this%solve_linprog(nx,nc,f,A,c,x,iret)
         if (iret.ne.0) return  ! failure
         zmin=-x(nx)
         Nm=x(1:nz)+zmin
      end subroutine maxmin_comp


      !> Determine x=xm which minimizes g=f'*x subject to x(i)>=0 and A*x=b,where A has full rank.
      subroutine solve_linprog(this,nx,nb,f,A,b,xm,iret)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         ! Ideal Gas Mixtures using Gibbs Function Continuation.
         use linprog, only: lp
         implicit none
         class(chem_state), intent(in) :: this
         integer, intent(in) :: nx,nb
         real(WP), intent(in)  :: f(nx),A(nb,nx),b(nb)
         real(WP), intent(out) :: xm(nx)
         integer, intent(out) :: iret
         ! Input:
         !	nx	- number of components of x
         !	nb	- number of components of b
         !	f	- nx-vector f
         !	A	- nx x nb matrix A
         ! Output:
         !	xm	- solution	
         !	iret= 0 if solution is found 
         !	iret=-1 if g is unbounded
         !	iret=-2 if there is no feasible solution
         !	iret=-3 if A is rank deficient
         !  S.B. Pope 10/1/06
         real(WP) :: eps=1e-9,ale(1,1),age(1,1),ble(1),bge(1)
         call lp(nx,0,0,nb,ale,age,A,ble,bge,b,f,xm,iret,toler=eps)
         if (iret<0) then
            !XXX write(0,*)'solve_linprog,iret=',iret  ! SBP XXX
         endif
      end subroutine solve_linprog


      !> Get the square root of moles
      function get_Nusqrt(this) result(Nusqrt)
         class(chem_state), intent(inout)  :: this
         real(WP), dimension(this%sys%nsu) :: Nusqrt
         Nusqrt=exp(0.5_WP*(-this%gu+matmul(this%sys%BR,this%sol(1:this%sys%nrc))+matmul(this%sys%P(this%sys%nsd+1:this%sys%ns,:),this%sol(this%sys%nrc+1:this%sys%nrc+this%sys%np))))
      end function get_Nusqrt


      !> Get the constraints residual vector
      subroutine get_RC(this,Nusqrt)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: Nusqrt
         this%RC(1:this%sys%nrc)=matmul(this%BtildeT,Nusqrt)-this%cr
         this%RC(this%sys%nrc+1:this%sys%nrc+this%sys%np)=matmul(this%PtildeT,Nusqrt)-this%Nbar+this%RCd
      end subroutine get_RC


      !> Initialize the solution unknowns
      subroutine sol_init(this)
         use messager,  only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(:), allocatable :: rhs,lam
         integer :: info
         this%success=.true.
         ! Allocate intermediate arrays
         allocate(rhs(this%sys%nrc))
         allocate(lam(this%sys%nrc))
         ! Calculate the Lagrange multipliers
         call this%get_gort(this%sys%nsu,this%T,this%p,this%sys%thermo(this%sys%nsd+1:this%sys%ns,:),this%sys%P(this%sys%nsd+1:this%sys%ns,Gphase),this%gu)
         rhs=log(this%Nu)-matmul(this%sys%P(this%sys%nsd+1:this%sys%ns,:),log(this%Nbar))+this%gu
         call lss(this%sys%nsu,this%sys%nrc,this%sys%BR,rhs,lam,info)
         if (info.ne.0) then
            this%success=.false.
            write(output_unit,'(" >   [chem_state sol_init] Least squares solver for lambda initialization failed")')
            return
            ! call die('[chem_state sol_init] Least squares solver for lambda initialization failed.')
         end if
         ! Set the initial solution vector
         this%sol(1:this%sys%nrc)=lam
         this%sol(this%sys%nrc+1:this%sys%nrc+this%sys%np)=log(this%Nbar)
         ! Deallocate intermediate arrays
         deallocate(rhs,lam)
      end subroutine sol_init


      !> Find the chemical equilibium state
      subroutine equilibrate(this)
         use messager, only: die
         implicit none
         class(chem_state), intent(inout) :: this
         ! Calculate the contribution of Nd in the residual
         if (this%sys%nsd.gt.0) then
            this%RCd=matmul(transpose(this%sys%P(1:this%sys%nsd,:)),this%Nd)
         else
            this%RCd=0.0_WP
         end if
         ! Get the chemical equilibrium state
         call this%get_ceq()
      end subroutine equilibrate


      !> Find the chemical equilibium state at constant pressure and temperature
      subroutine get_ceq_PT(this)
         use messager, only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         integer :: i,j,iJ,jJ,info
         real(WP), dimension(:,:), allocatable :: Jac
         real(WP), dimension(:),   allocatable :: dsol,Nusqrt
         real(WP), dimension(:,:), allocatable :: BTB,PTP,BTP
         real(WP), dimension(:),   allocatable :: S,work
         real(WP) :: Rnorm,rcond
         integer  :: rank,lwork,iter
         ! Allocate arrays
         allocate(Jac   (this%sys%nrc+this%sys%np,this%sys%nrc+this%sys%np)); Jac=0.0_WP
         allocate(dsol  (this%sys%nrc+this%sys%np));  dsol  =0.0_WP
         allocate(BTB   (this%sys%nrc,this%sys%nrc)); BTB   =0.0_WP
         allocate(BTP   (this%sys%nrc,this%sys%np));  BTP   =0.0_WP
         allocate(PTP   (this%sys%np,this%sys%np));   PTP   =0.0_WP
         allocate(Nusqrt(this%sys%nsu));              Nusqrt=0.0_WP
         allocate(S     (this%sys%nrc+this%sys%np))
         lwork=10*(this%sys%nrc+this%sys%np)
         allocate(work(lwork))
         rcond=-1.0_WP
         ! Initialize the solution vector
         call this%sol_init()
         if (.not.this%success) return
         Nusqrt=this%get_Nusqrt()
         ! Newton-Raphson
         this%iter_N=0
         do iter=1,this%iter_N_max
            ! Build the coefficient matrices
            call this%get_BP(Nusqrt)
            ! Get the residual error
            call this%get_RC(Nusqrt)
            Rnorm=norm2(this%RC)
            ! Evaluate the error
            if (Rnorm.lt.this%tol_N) then
               this%iter_N=iter
               exit
            end if
            ! Build the Jacobian matrix
            BTB=matmul(this%BtildeT,this%Btilde)
            PTP=matmul(this%PtildeT,this%Ptilde)
            BTP=matmul(this%BtildeT,this%Ptilde)
            do j=1,this%sys%nrc
               jJ=j
               do i=1,j
                  iJ=i
                  Jac(iJ,jJ)=BTB(i,j)
               end do
            end do
            do j=1,this%sys%np
               jJ=this%sys%nrc+j
               do i=1,this%sys%nrc
                  iJ=i
                  Jac(iJ,jJ)=BTP(i,j)
               end do
            end do
            do j=1,this%sys%np
               jJ=this%sys%nrc+j
               do i=1,j
                  iJ=this%sys%nrc+i
                  Jac(iJ,jJ)=PTP(i,j)
               end do
            end do
            do j=1,this%sys%nrc+this%sys%np-1
               do i=j+1,this%sys%nrc+this%sys%np
                  Jac(i,j)=Jac(j,i)
               end do
            end do
            do i=1,this%sys%np
               iJ=this%sys%nrc+i
               jJ=this%sys%nrc+i
               Jac(iJ,jJ)=Jac(iJ,jJ)-this%Nbar(i)
            end do
            ! Solve for dsol
            dsol=-this%RC
            call dgelss(this%sys%nrc+this%sys%np,this%sys%nrc+this%sys%np,1,Jac,this%sys%nrc+this%sys%np,dsol,this%sys%nrc+this%sys%np,S,rcond,rank,work,lwork,info)
            ! if (rank.ne.this%sys%nrc+this%sys%np) call die('[chem_state get_ceq_PT]: Jacobian is not full rank')
            if (rank.ne.this%sys%nrc+this%sys%np) then
               this%success=.false.
               return
            end if
            ! if (info.ne.0) call die('[chem_state get_ceq_PT]: Least-squares solver failed')
            if (info.ne.0) then
               this%success=.false.
               return
            end if
            ! Update the solution
            this%sol=this%sol+dsol
            ! Get the species and phase moles
            Nusqrt=this%get_Nusqrt()
            this%Nu=Nusqrt*Nusqrt
            this%Nbar=exp(this%sol(this%sys%nrc+1:this%sys%nrc+this%sys%np))
         end do
         if (iter.gt.this%iter_N_max) then
            this%iter_N=iter-1
            this%success=.false.
            return
         end if
         ! Assemble the composition
         this%Ndu=[this%Nd,this%Nu]
         ! Reorder the composition
         do i=1,this%sys%ns
            this%N(this%sys%sp_order(i))=this%Ndu(i)
         end do
         ! Deallocate arrays
         deallocate(Jac,dsol,BTB,BTP,PTP,Nusqrt,S,work)
      end subroutine get_ceq_PT


      !> Find the chemical equilibium state at constant pressure and enthalpy (Newton-Raphson method)
      subroutine get_ceq_PH_NR(this)
         ! Extracted from Pope, Stephen. (2003). The Computation of Constrained and Unconstrained Equilibrium Compositions of 
         use messager, only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(:), allocatable :: hort
         real(WP) :: Tn
         real(WP) :: hlo,hhi,Cp_eff
         integer :: i
         ! Allocate arrays
         allocate(hort(this%sys%ns))
         ! Initialize the constant PH iterations
         this%HoR0=this%HoR
         this%dT=1e5*this%tol_T*this%T
         this%iter_T=0
         this%Tlo=-1e30
         this%Thi= 1e30
         ! Iterate over temperature
         do while (abs(this%dT/this%T).ge.this%tol_T)
            ! Increment the iterations
            this%iter_T=this%iter_T+1
            if (this%iter_T.gt.this%iter_T_max) then
               this%iter_T=this%iter_T-1
               this%success=.false.
               return
            end if
            ! Store the old mole numbers
            this%Nuold=this%Nu
            ! Get the enthalpy residual
            call this%get_RH(T=this%T,RH=this%RH)
            if (.not.this%success) then
               return
            end if
            ! Get the effective Cp
            call this%get_Cp_eff(Cp_eff)
            ! Predict dT
            this%dT=-this%RH/Cp_eff
            ! Check that T is within limits
            if (this%T.eq.T_high.and.this%dT.gt.0.0_WP) then
               this%success=.false.
               ! write(output_unit,'(" >   [chem_state get_ceq_PH] T > T_high")')
               return
            end if
            if (this%T.eq.T_low .and.this%dT.lt.0.0_WP) then
               this%success=.false.
               ! write(output_unit,'(" >   [chem_state get_ceq_PH] T < T_low")')
               return
            end if
            ! Ensure that Tn is within limits
            Tn=this%T+this%dT
            Tn=max(min(Tn,T_high),T_low)
            ! Use linear interpolation instead if Tn is closer to known bound
            if (this%dT.gt.0.0_WP) then
               this%Tlo=this%T
               hlo=this%HoR
               if (Tn.gt.0.5_WP*(this%Tlo+this%Thi)) then
                  Tn=this%Tlo+(this%Thi-this%Tlo)*(this%HoR0-hlo)/(hhi-hlo)
               endif
            else
               this%Thi=this%T
               hhi=this%HoR
               if (Tn.lt.0.5_WP*(this%Tlo+this%Thi)) then
                  Tn=this%Tlo+(this%Thi-this%Tlo)*(this%HoR0-hlo)/(hhi-hlo)
               endif
            endif
            ! Update temperature increment
            this%dT=Tn-this%T
            ! Update temperature
            this%T=Tn
         end do
         ! Assemble the composition
         this%Ndu=[this%Nd,this%Nu]
         ! Reorder the composition
         do i=1,this%sys%ns
            this%N(this%sys%sp_order(i))=this%Ndu(i)
         end do
         ! Dellocate arrays
         deallocate(hort)
      end subroutine get_ceq_PH_NR


      !> Find the chemical equilibium state at constant pressure and enthalpy (Bi-section method)
      subroutine get_ceq_PH_BS(this)
         use messager, only: die
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP) :: Tm
         real(WP) :: Cp_eff,Rlo,Rhi
         integer  :: i,iter
         ! Initialize
         this%HoR0=this%HoR
         ! Iterate
         do iter=1,this%iter_T_max
            ! Get equilibrium for Tm
            Tm=0.5_WP*(this%Tlo+this%Thi)
            call this%get_RH(T=Tm,RH=this%RH)
            if (.not.this%success) then
               return
            end if
            ! Evaluate the enthalpy residual
            if (abs(this%RH/this%HoR0).lt.this%tol_H) then
               ! Assign the iteration number
               this%iter_T=iter
               ! Assemble the composition
               this%Ndu=[this%Nd,this%Nu]
               ! Reorder the composition
               do i=1,this%sys%ns
                  this%N(this%sys%sp_order(i))=this%Ndu(i)
               end do
               return
            end if
            ! Get the equilibrium for Tlo
            call this%get_RH(T=this%Tlo,RH=Rlo)
            if (.not.this%success) then
               return
            end if
            ! Get the equilibrium for Thi
            call this%get_RH(T=this%Thi,RH=Rhi)
            if (.not.this%success) then
               return
            end if
            ! Adjust the temperature bounds
            if (Rlo*this%RH.lt.0.0_WP) then
               this%Thi=Tm
            else
               this%Tlo=Tm
            endif
         end do
         this%success=.false.
         this%iter_T=iter-1
      end subroutine get_ceq_PH_BS


      !> Find the chemical equilibium state at constant internal energy and volume
      subroutine get_ceq_UV(this)
         use, intrinsic :: iso_fortran_env, only: output_unit
         implicit none
         class(chem_state), intent(inout) :: this
         real(WP), dimension(:,:), allocatable :: Jac
         real(WP), dimension(:),   allocatable :: rhs
         real(WP) :: alpha,pn,Tn
         integer  :: info,i,ipiv
         ! Allocate the intermediate arrays
         allocate(Jac(1:2,1:2))
         allocate(rhs(1:2))
         ! Initialize the constant UV iterations
         this%UoR0=this%UoR
         this%V0  =this%V
         this%dT=1e5_WP*this%tol_T*this%T
         this%dp=1e5_WP*this%tol_p*this%p
         this%iter_T=0
         ! Iterate over temperature and pressure
         do while ((abs(this%dT/this%T).ge.this%tol_T).or.(abs(this%dp/this%p).ge.this%tol_p))
            ! Increment the iterations
            this%iter_T=this%iter_T+1
            if (this%iter_T.gt.this%iter_T_max) then
               this%iter_T=this%iter_T-1
               this%success=.false.
               return
            end if
            ! Store the old mole numbers
            this%Nuold=this%Nu
            ! Get the residuals and Jacobian
            call this%get_RUV(T=this%T,p=this%p,RU=this%RU,RV=this%RV,Jac=Jac)
            if (.not.this%success) then
               return
            end if
            ! Solve for the residuals using LU decomposition
            rhs=[-this%RU,-this%RV]
            call dgesv(2,1,Jac,2,ipiv,rhs,2,info)
            this%dT=rhs(1)
            this%dp=rhs(2)
            if (info.ne.0) then
               this%success=.false.
               return
            end if
            ! Relax the residuals
            alpha=1.0_WP
            do
               Tn=this%T+alpha*this%dT
               pn=this%p+alpha*this%dp
               if (Tn.ge.T_low.and.Tn.le.T_high.and.pn.gt.0.0_WP) exit
               alpha=0.5_WP*alpha
               if (alpha.lt.1e-8_WP) then
                  this%success=.false.
                  return
               end if
            end do
            this%dT=alpha*this%dT
            this%dp=alpha*this%dp
            this%T=Tn
            this%p=pn
         end do
         this%Ndu=[this%Nd,this%Nu]
         do i=1,this%sys%ns
            this%N(this%sys%sp_order(i))=this%Ndu(i)
         end do
         ! Deallocate the intermediate arrays
         deallocate(Jac,rhs)
      end subroutine get_ceq_UV


      !> Get the enthalpy residual for a constant PH iteration
      subroutine get_RH(this,T,RH)
         class(chem_state), intent(inout) :: this
         real(WP), intent(in)  :: T
         real(WP), intent(out) :: RH
         real(WP), dimension(:), allocatable :: hort
         ! Allocate arrays
         allocate(hort(this%sys%ns))
         ! Assign temperature
         this%T=T
         ! Re-initialize mole numbers using current temperature
         call this%N_re_init()
         if (.not.this%success) then
            return
         end if
         ! Determine equilibrium composition at current temperature
         call this%get_ceq_PT()
         if (.not.this%success) then
            return
         end if
         ! Obtain species molar h/(RT)
         call this%get_hort(this%sys%ns,this%T,this%sys%thermo,hort)
         ! Mixture H/R
         this%HoR=this%T*sum(this%Ndu*hort)
         RH=this%HoR-this%HoR0
      end subroutine get_RH


      !> Get the internal energy and volume residuals for a constant UV iteration
      subroutine get_RUV(this,T,p,RU,RV,Jac)
         class(chem_state), intent(inout) :: this
         real(WP), intent(in)  :: T,p
         real(WP), intent(out) :: RU,RV
         real(WP), dimension(:,:), optional :: Jac
         real(WP), dimension(:), allocatable :: isGas,isLiq,hort,cpor,dNdT,dNdp,uor,vmolar
         real(WP) :: Cp_eff
         ! Allocate arrays
         allocate(isGas (this%sys%ns))
         allocate(isLiq (this%sys%ns))
         allocate(hort  (this%sys%ns))
         allocate(cpor  (this%sys%ns))
         allocate(dNdT  (this%sys%ns))
         allocate(dNdp  (this%sys%ns))
         allocate(uor   (this%sys%ns))
         allocate(vmolar(this%sys%ns))
         ! Get the gas and liquid indices
         isGas=this%sys%P(:,Gphase)
         isLiq=1.0_WP-isGas
         ! Assign temperature and pressure
         this%T=T
         this%p=p
         ! Re-initialize mole numbers using current temperature and pressure
         call this%N_re_init()
         if (.not.this%success) then
            return
         end if
         ! Determine equilibrium composition at current temperature and pressure
         call this%get_ceq_PT()
         if (.not.this%success) then
            return
         end if
         ! Obtain species molar h/(RT)
         call this%get_hort(this%sys%ns,this%T,this%sys%thermo,hort)
         ! Get the internal energy: u/R = (h - pv) / R = T * (h/(R*T) - pv/(RT)); pv/(RT) = 1 for ideal gas
         uor=this%T*(hort-isGas)-isLiq*this%p*this%vmolar/gas_cnst
         ! Get the molar volumes
         vmolar=isGas*gas_cnst*this%T/this%p+isLiq*this%vmolar
         ! Update the internal energy and volume
         this%UoR=sum(this%Ndu*uor)
         this%V  =sum(this%Ndu*vmolar)
         ! Get the residuals
         RU=this%UoR-this%UoR0
         RV=this%V  -this%V0
         ! Get the Jacobian if needed
         if (present(Jac)) then
            ! Obtain the temperature and pressure derivatives of the mole numbers
            call this%get_dNdT(dNdT)
            call this%get_dNdp(dNdp)
            ! Obtain species molar Cp/R
            call this%get_cpor(this%sys%ns,this%T,this%sys%thermo,cpor)
            ! Get the constrained effective specific heat at constant pressure
            Cp_eff=sum(cpor*this%Ndu)+this%T*sum(hort*dNdT)
            ! Form the Jacobian
            Jac(1,1)=sum(uor*dNdT)+sum(this%Ndu*(cpor-isGas))
            Jac(1,2)=sum(uor*dNdp)
            Jac(2,1)=gas_cnst*(Cp_eff-Jac(1,1))/this%p
            Jac(2,2)=sum(vmolar*dNdp)+sum(this%Ndu*isGas*(-gas_cnst*this%T/this%p**2))
         end if
         deallocate(isGas,isLiq,hort,cpor,dNdT,dNdp,vmolar)
      end subroutine get_RUV


      !> Evaluate the effective heat capacity (extensive)
      !> Cp_eff = (dH/dT)/R at constant constraints and p
      subroutine get_Cp_eff(this,Cp_eff)
         class(chem_state), intent(inout) :: this
         real(WP), intent(out) :: Cp_eff
         real(WP), dimension(:), allocatable :: dgudT,dxdT,cpor,hort,N,dNddT,dNudT,dNdT
         ! Allocate arrays
         allocate(cpor(this%sys%ns))
         allocate(hort(this%sys%ns))
         allocate(dNdT(this%sys%ns))
         ! Get the specific heat and enthalpy
         call this%get_cpor(this%sys%ns,this%T,this%sys%thermo,cpor)
         call this%get_hort(this%sys%ns,this%T,this%sys%thermo,hort)
         ! Rates of change of moles
         call this%get_dNdT(dNdT)
         Cp_eff=sum(cpor*this%Ndu)+this%T*sum(hort*dNdT)
         ! Deallocate arrays
         deallocate(cpor,hort,dNdT)
      end subroutine get_Cp_eff


      !> Get the derivative of the solution vector with respect to generic parameter (par)
      subroutine get_dsoldpar(this,dgudpar,dxdpar)
         use messager,  only: die
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: dgudpar
         real(WP), dimension(this%sys%nrc+this%sys%np), intent(out) :: dxdpar
         real(WP), dimension(:),   allocatable :: lamdotg,Sig,Sinv,work,Nusqrt,Ygdot,dlnNbardpar,rhs
         real(WP), dimension(:,:), allocatable :: Btildeinv,lamdoty,U,VT,M,Btilde_cp
         real(WP) :: srlim=1e-9
         integer :: info,lwork,i,n_small
         ! Allocate arrays
         allocate(Btildeinv (this%sys%nrc,this%sys%nsu))
         allocate(lamdotg   (this%sys%nrc))
         allocate(Sig       (this%sys%nrc))
         allocate(Sinv      (this%sys%nrc))
         allocate(work      (20*(this%sys%nsu+this%sys%nrc)))
         allocate(Nusqrt    (this%sys%nsu))
         allocate(Ygdot     (this%sys%nsu))
         allocate(dlnNbardpar(this%sys%np))
         allocate(rhs       (this%sys%np))
         allocate(lamdoty   (this%sys%nrc,this%sys%np))
         allocate(U         (this%sys%nsu,this%sys%nrc))
         allocate(VT        (this%sys%nrc,this%sys%nrc))
         allocate(M         (this%sys%np,this%sys%np))
         allocate(Btilde_cp (this%sys%nsu,this%sys%nrc))
         lwork=size(work)
         ! Get the Nusqrt vector and update the coefficient matrices
         Nusqrt=this%get_Nusqrt()
         call this%get_BP(Nusqrt)
         Btilde_cp=this%Btilde(1:this%sys%nsu,1:this%sys%nrc)
         ! Get the SVD of Btilde
         call dgesvd('S','A',this%sys%nsu,this%sys%nrc,Btilde_cp,this%sys%nsu,Sig(1:this%sys%nrc),U(1:this%sys%nsu,1:this%sys%nrc),this%sys%nsu,VT(1:this%sys%nrc,1:this%sys%nrc),this%sys%nrc,work(1:lwork),lwork,info)
         if (info.ne.0) call die('[chem_state get_dsoldpar] SVD of B tilde failed')
         ! Get the inverse of Sigma
         call get_Sinv(this%sys%nrc,Sig,Sinv,srlim,n_small)
         ! Store Sinv * V' in VT
         do i=1,this%sys%nrc
            VT(i,:)=Sinv(i)*VT(i,:)
         end do
         ! Btilde^-1 = V * Sinv * U'
         Btildeinv=transpose(matmul(U,VT))
         ! Solve for dsol/dpar
         Ygdot=Nusqrt*dgudpar
         lamdotg=matmul(Btildeinv,Ygdot)
         lamdoty=matmul(Btildeinv,this%Ptilde)
         M=matmul(this%PtildeT,matmul(this%Btilde,lamdoty))
         rhs=matmul(this%PtildeT,matmul(this%Btilde,lamdotg)-Ygdot)
         call lss(this%sys%np,this%sys%np,M,rhs,dlnNbardpar,info)
         if (info.ne.0) call die('[chem_state get_dsoldpar] Least squares solver failed')
         dxdpar(1:this%sys%nrc)=lamdotg-matmul(lamdoty,dlnNbardpar)
         dxdpar(this%sys%nrc+1:this%sys%nrc+this%sys%np)=dlnNbardpar
         ! Deallocate arrays
         deallocate(Btildeinv,lamdotg,lamdoty,Sig,Sinv,work,Nusqrt,Ygdot,dlnNbardpar,rhs,U,VT,M,Btilde_cp)
         contains
            ! Get the inverse of Sigma
            subroutine get_Sinv(n,S,Si,srat_lim,n_s)
               ! Given the n-vector of singular values, S, which are in decreasing order,
               ! return the n-vector of pseudo-inverses, Si, and the number n_s of
               ! small singular values.  The j-th singular value is deemed to be small
               ! if S(j)/S(1) < srat_lim.
               integer,  intent(in)  :: n
               real(WP), intent(in)  :: S(n),srat_lim
               integer,  intent(out) :: n_s
               real(WP), intent(out) :: Si(n)
               integer  :: j
               real(WP) :: slim
               slim=srat_lim*S(1)
               do j=1,n
                  if (S(j).gt.slim) then
                     Si(j)=1.0_WP/S(j)
                     n_s=n-j
                  else
                     Si(j)=0.0_WP
                  endif
               end do
            end subroutine get_Sinv
      end subroutine get_dsoldpar


      !> Get dN/dT using first order FD
      subroutine get_dNdT_FD(this,dNdT)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(out) :: dNdT
         call this%get_dNdpar_FD(this%dT,dNdT)
      end subroutine get_dNdT_FD


      !> Get dN/dp using first order FD
      subroutine get_dNdp_FD(this,dNdp)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(out) :: dNdp
         call this%get_dNdpar_FD(this%dp,dNdp)
      end subroutine get_dNdp_FD


      !> Get dN/dT using the least squares approach
      subroutine get_dNdT_LS(this,dNdT)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(out) :: dNdT
         real(WP), dimension(:), allocatable :: dgudT
         ! Allocate intermediate arrays
         allocate(dgudT(this%sys%nsu))
         ! Get d/dT of the normalized Gibbs functions of the undetermined species
         call this%get_dgdT(this%sys%nsu,this%T,this%sys%thermo(this%sys%nsd+1:this%sys%ns,:),dgudT)
         ! Get dN/dT
         call this%get_dNdpar_LS(dgudT,dNdT)
         ! Deallocate intermediate arrays
         deallocate(dgudT)
      end subroutine get_dNdT_LS


      !> Get dN/dp using the least squares approach
      subroutine get_dNdp_LS(this,dNdp)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%ns), intent(out) :: dNdp
         real(WP), dimension(:), allocatable :: dgudp
         ! Allocate intermediate arrays
         allocate(dgudp(this%sys%nsu))
         ! Get d/dp of the normalized Gibbs functions of the undetermined species
         call this%get_dgdp(this%sys%nsu,this%T,this%p,this%sys%P(this%sys%nsd+1:this%sys%ns,Gphase),this%vmolar,dgudp)
         ! Get dN/dp
         call this%get_dNdpar_LS(dgudp,dNdp)
         ! Deallocate intermediate arrays
         deallocate(dgudp)
      end subroutine get_dNdp_LS


      !> Get the mole numbers derivative with respect to a generic paramenter (par) using first order FD
      subroutine get_dNdpar_FD(this,dpar,dNdpar)
         class(chem_state), intent(inout) :: this
         real(WP), intent(in) :: dpar
         real(WP), dimension(this%sys%ns), intent(out) :: dNdpar
         real(WP), dimension(:), allocatable :: dNddpar,dNudpar
         ! Allocate arrays
         allocate(dNddpar(this%sys%nsd))
         allocate(dNudpar(this%sys%nsu))
         dNddpar=0.0_WP
         dNudpar=0.0_WP
         if (this%iter_T.gt.1) dNudpar=(this%Nu-this%Nuold)/dpar
         dNdpar=[dNddpar,dNudpar]
         ! Deallocate arrays
         deallocate(dNddpar,dNudpar)
      end subroutine get_dNdpar_FD


      !> Get the mole numbers derivative with respect to a generic paramenter (par) using the least squares approach
      subroutine get_dNdpar_LS(this,dgudpar,dNdpar)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: dgudpar
         real(WP), dimension(this%sys%ns), intent(out) :: dNdpar
         real(WP), dimension(:), allocatable :: dxdpar,dNddpar,dNudpar
         ! Allocate arrays
         allocate(dxdpar (this%sys%nrc+this%sys%np))
         allocate(dNddpar(this%sys%nsd))
         allocate(dNudpar(this%sys%nsu))
         ! Get dlambda/dpar and dlnNbar/dpar
         call this%get_dsoldpar(dgudpar,dxdpar)
         ! Get dN/dpar
         dNddpar=0.0_WP
         dNudpar=this%Nu*(-dgudpar+matmul(this%sys%BR,dxdpar(1:this%sys%nrc))+dxdpar(this%sys%nrc+1:this%sys%nrc+this%sys%np))
         dNdpar=[dNddpar,dNudpar]
         ! Deallocate arrays
         deallocate(dxdpar,dNddpar,dNudpar)
      end subroutine get_dNdpar_LS


      !> Get Btilde and Ptilde
      subroutine get_BP(this,Nusqrt)
         class(chem_state), intent(inout) :: this
         real(WP), dimension(this%sys%nsu), intent(in) :: Nusqrt
         integer :: j
         do j=1,this%sys%nrc
            this%Btilde(:,j)=Nusqrt*this%sys%BR(:,j)
         end do
         do j=1,this%sys%np
            this%Ptilde(:,j)=Nusqrt*this%sys%P(this%sys%nsd+1:this%sys%ns,j)
         end do
         this%BtildeT=transpose(this%Btilde)
         this%PtildeT=transpose(this%Ptilde)
      end subroutine get_BP


      !> Determine the least-squares/minimum-norm solution x to the linear equation Ax = b.
      subroutine lss(nb,nx,A,b,x,info)
         !	S.B. Pope 10/2/02
         implicit none
         integer,  intent(in)  :: nb,nx
         real(WP), intent(in)  :: A(nb,nx),b(nb)
         real(WP), intent(out) :: x(nx)
         integer,  intent(out) :: info
         !  Input:
         !	nb	- number of rows in b
         !	nx	- number of rows in A and x
         !	A	- the nb x nx matrix A
         !	b	- the nb-vector b
         !  Output:
         !	x	- the solution nx-vector
         !	info=0 for successful solution
         integer :: lwork,rank
         real(WP) :: tol=1.d-9,aa(nb,nx),bb(nb+nx),sv(nb+nx),work(4*(nb+nx+1)*(nb+nx+1))
         lwork= size(work)
         aa=A
         bb=0.d0
         bb(1:nb)=b
         call dgelss(nb,nx,1,aa(1:nb,1:nx),nb,bb(1:nb+nx),nb+nx,sv(1:nb+nx),tol,rank,work(1:lwork),lwork,info)
         x=bb(1:nx)
      end subroutine lss


end module chem_state_class