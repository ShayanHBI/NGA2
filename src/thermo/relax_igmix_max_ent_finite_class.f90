!> Finite-rate maximum-entropy relaxation (relaxation_ent_max.tex, Finite-rate relaxation): rate law dx~/dt=R*F^,
!> R=diag(lambda,vartheta,nu). The transfers faster than the finite one are instantaneous, the slower ones inactive;
!> the finite stage moves dx~=omega*K^-1*F^0, omega=1-exp(-r*K*dt) (eq. finite_solution), with K its pivot at the step start.
module relax_igmix_max_ent_finite_class
   use precision,                 only: WP
   use relax_igmix_max_ent_class, only: relax_igmix_max_ent,mestate,max_ent_apply
   use relax_igmix_sg_class,      only: Prelax,PTrelax,PTgrelax
   use thermorelax_class,         only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS
   implicit none
   private

   public :: relax_igmix_max_ent_finite
   public :: closure_restore,closure_pelanti

   !> Closure after the finite step: restore the faster equilibria by their stages (Activation), or Pelanti's
   !> mixture-energy relation at the new V_l and masses (Pelanti 2022 eqs. 59, 70, 82)
   integer, parameter :: closure_restore=1,closure_pelanti=2

   type, extends(relax_igmix_max_ent) :: relax_igmix_max_ent_finite
      real(WP) :: MRR =huge(1.0_WP)   !< Mechanical rate lambda [1/(Pa s)]; huge = instantaneous
      real(WP) :: TRR =huge(1.0_WP)   !< Thermal rate vartheta [W/(m^3 K)]; huge = instantaneous
      real(WP) :: CRR =huge(1.0_WP)   !< Chemical rate nu [kg s/m^5]; huge = instantaneous
      real(WP) :: tauC=-1.0_WP        !< Chemical relaxation time; >0 sets nu=1/(K^C*tauC) per cell
      integer  :: closure=closure_restore
   contains
      procedure :: apply
      procedure :: p_relax
      procedure :: pT_relax
      procedure :: pTg_relax
      procedure :: finite_step
      procedure :: pelanti_closure
   end type relax_igmix_max_ent_finite

contains

   !> Rate-separation checks, then the instantaneous class's apply on this, so that it dispatches to the finite workers below
   subroutine apply(this,dt,VF,Q,Pjump,ierr)
      use messager, only: die
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      logical :: nuc0
      if (this%model.ne.Prelax.and.this%MRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pT and pTg require MRR=huge')
      if (this%model.eq.PTgrelax.and.this%TRR.lt.huge(1.0_WP)) call die('[relax_igmix_max_ent_finite apply] pTg requires TRR=huge')
      ! No nucleation at zero chemical rate: the seed is itself a phase change
      nuc0=this%do_nucleate
      if (this%CRR.le.0.0_WP.and.this%tauC.le.0.0_WP) this%do_nucleate=.false.
      call max_ent_apply(this,dt,VF,Q,Pjump,ierr)
      this%do_nucleate=nuc0
   end subroutine apply

   !> Mechanical relaxation at rate lambda
   subroutine p_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      if (this%MRR.ge.huge(1.0_WP)) then
         call this%stage(1,VF,Q,Pjump,ier)
      else
         call this%finite_step(1,dt,VF,Q,Pjump,ier)
      end if
      if (present(ierr)) ierr=ier
   end subroutine p_relax

   !> Instantaneous mechanical stage, then thermal relaxation at rate vartheta
   subroutine pT_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      call this%stage(1,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) then
         if (this%TRR.ge.huge(1.0_WP)) then
            call this%stage(2,VF,Q,Pjump,ier)
         else
            call this%finite_step(2,dt,VF,Q,Pjump,ier)
         end if
      end if
      if (present(ierr)) ierr=ier
   end subroutine pT_relax

   !> Instantaneous mechanical and thermal stages, then chemical relaxation at rate nu
   subroutine pTg_relax(this,dt,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,  optional,                intent(out)   :: ierr
      integer :: ier
      call this%stage(1,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) call this%stage(2,VF,Q,Pjump,ier)
      if (ier.eq.RELAX_OK) then
         if (this%CRR.ge.huge(1.0_WP).and.this%tauC.le.0.0_WP) then
            call this%stage(3,VF,Q,Pjump,ier)
         else
            call this%finite_step(3,dt,VF,Q,Pjump,ier)
         end if
      end if
      if (present(ierr)) ierr=ier
   end subroutine pTg_relax

   !> Finite-rate step of stage ns from a state where the faster forces vanish (eq. finite_solution), then the closure
   subroutine finite_step(this,ns,dt,VF,Q,Pjump,ierr)
      use relax_igmix_nasg_finite_class, only: log_E1,inv_log_E1
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      integer,                           intent(in)    :: ns
      real(WP),                          intent(in)    :: dt
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      type(mestate) :: st
      real(WP), dimension(3,3) :: D,J
      real(WP), dimension(3) :: dx,K
      real(WP) :: rate,omega,Rv,ms,u0,u1,tgt
      call this%get_state(VF,Q,Pjump,st,ierr)
      if (ierr.ne.RELAX_OK) return
      D=this%get_D(st,Pjump)
      J=this%get_jac(st,D,Pjump)
      Rv=this%gas%R(this%indV)
      if (ns.eq.3.and.(st%dry.or.st%xv*st%pg.lt.this%pv_min)) then
         ! Departure (as relax_igmix_nasg_finite): dry gas keeps only the mixing force F_C=Rv*T*ln(m_s/m_v) and integrates
         ! it exactly, E1(ln(m_s/m_v))=E1(ln(m_s/m_v^0))+nu*Rv*T*dt/m_s, with m_s the dilute saturation vapor mass
         ms=st%ma*this%gas%R(this%indA)/Rv*exp((st%gl-(this%gas%cp(this%indV)*st%Tg+this%gas%q(this%indV) &
         &  -st%Tg*(this%gas%cp(this%indV)*log(st%Tg)+this%gas%qp(this%indV))))/(Rv*st%Tg)-log(st%pg))
         u0=st%mv/ms
         if (u0.ge.1.0_WP) return
         if (this%tauC.gt.0.0_WP) then
            this%CRR_loc=ms/(Rv*st%Tg*this%tauC)
         else
            this%CRR_loc=this%CRR
         end if
         u1=u0
         if (this%CRR_loc.gt.0.0_WP) then
            this%tauC_loc=ms/(this%CRR_loc*Rv*st%Tg)
            u1=1.0_WP
            if (log(this%CRR_loc)+log(Rv*st%Tg*dt/ms).lt.log(690.0_WP)) then
               tgt=this%CRR_loc*Rv*st%Tg*dt/ms
               if (u0.gt.0.0_WP) tgt=tgt+exp(log_E1(-log(u0)))
               if (tgt.lt.690.0_WP) u1=exp(-inv_log_E1(log(tgt)))
            end if
         end if
         call this%eliminate(ns,J,st%F,st%mv,.false.,1.0_WP,dx,K,ierr,dxs=(max(u1,u0)-u0)*ms)
         if (ierr.ne.RELAX_OK) return
      else
         ! Pivot K of the stage at the step start (eq. stiff_def)
         call this%eliminate(ns,J,st%F,st%mv,.false.,1.0_WP,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
         select case (ns)
         case (1); rate=this%MRR; this%KM=K(1)
         case (2); rate=this%TRR; this%KT=K(2)
         case (3)
            rate=this%CRR; this%KC=K(3)
            if (this%tauC.gt.0.0_WP) rate=1.0_WP/(K(3)*this%tauC)
            this%CRR_loc=rate
            if (rate.gt.0.0_WP) this%tauC_loc=1.0_WP/(rate*K(3))
         end select
         ! omega=1-exp(-r*K*dt), 1 once r*K*dt>=50
         omega=0.0_WP
         if (rate.gt.0.0_WP) then
            omega=1.0_WP
            if (log(rate)+log(K(ns)*dt).lt.log(50.0_WP)) omega=1.0_WP-exp(-rate*K(ns)*dt)
         end if
         call this%eliminate(ns,J,st%F,st%mv,.false.,omega,dx,K,ierr)
         if (ierr.ne.RELAX_OK) return
      end if
      call this%transfer(D,dx,Pjump,VF,Q)
      ! Closure of the faster equilibria
      if (ns.eq.1) return
      select case (this%closure)
      case (closure_pelanti)
         call this%pelanti_closure(VF,Q,Pjump,ierr)
      case default
         call this%stage(ns-1,VF,Q,Pjump,ierr)
      end select
   end subroutine finite_step

   !> Pelanti 2022 eqs. 59, 70, 82: V_l and the masses kept, energies split so that p_l-p_g=p_j with E_l+E_g fixed
   subroutine pelanti_closure(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: gG,qG,Vlb,p
      y=0.0_WP; y(this%indV)=min(Q(this%iVQ())/Q(2),1.0_WP); y(this%indA)=1.0_WP-y(this%indV)
      gG=sum(y*this%gas%cp(1:this%gas%ns))/sum(y*this%gas%cv(1:this%gas%ns))
      qG=sum(y*this%gas%q(1:this%gas%ns))
      Vlb=max(VF-this%bL*Q(1),(1.0_WP-this%brhomax)*VF)
      p=(Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*qG-Vlb*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)*Pjump/(gG-1.0_WP)) &
      & /(Vlb/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)/(gG-1.0_WP))
      ierr=RELAX_BAD_LIQUID
      if (p+this%liq%pinf.le.0.0_WP) return
      ierr=RELAX_BAD_GAS
      if (p-Pjump.le.0.0_WP) return
      Q(3)=VF*this%liq%get_rhoe_from_p_rho(p=p,rho=Q(1)/VF,y=[1.0_WP])
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=p-Pjump,rho=Q(2)/(1.0_WP-VF),y=y)
      ierr=RELAX_OK
   end subroutine pelanti_closure

end module relax_igmix_max_ent_finite_class
