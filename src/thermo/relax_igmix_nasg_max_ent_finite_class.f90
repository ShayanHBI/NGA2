!> Finite-rate maximum-entropy relaxation (relaxation_ent_max.tex, Finite-rate relaxation) for a NASG (or stiffened-gas, b=0)
!> liquid: relax_igmix_max_ent_finite with the closed forms of the mechanical stage with p_I=p* and of the pT-relaxed state,
!> the liquid state from the NASG relations, and Pelanti's closed-form mixture-energy closure
module relax_igmix_nasg_max_ent_finite_class
   use precision,                           only: WP
   use relax_igmix_max_ent_inst_class,      only: relax_igmix_max_ent_inst,phasic_state,pI_end
   use relax_igmix_max_ent_finite_class,    only: relax_igmix_max_ent_finite
   use relax_igmix_nasg_max_ent_inst_class, only: nasg_covolume,nasg_liq_state,nasg_liq_v_from_p_T,nasg_relax_M,nasg_relax_MT,nasg_admissible
   use stiffened_gas_class,                 only: stiffened_gas
   use igmix_class,                         only: igmix
   use thermorelax_class,                   only: RELAX_OK,RELAX_BAD_LIQUID,RELAX_BAD_GAS
   implicit none
   private

   public :: relax_igmix_nasg_max_ent_finite

   type, extends(relax_igmix_max_ent_finite) :: relax_igmix_nasg_max_ent_finite
      !> Liquid co-volume (0 for a plain stiffened gas); NASG is used without the brhomax clamp of nasg_class
      real(WP) :: bL=0.0_WP
   contains
      procedure :: initialize
      procedure :: admissible
      procedure :: liq_state
      procedure :: liq_v_from_p_T
      procedure :: relax_M
      procedure :: relax_MT
      procedure :: mixture_energy_closure
   end type relax_igmix_nasg_max_ent_finite

contains

   !> Parent initialization, then the NASG co-volume and p_I=p* for the instantaneous mechanical stage
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(inout) :: this
      class(stiffened_gas), target, intent(in) :: liq
      class(igmix),         target, intent(in) :: gas
      integer,                      intent(in) :: indV,indA
      call this%relax_igmix_max_ent_finite%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      this%bL=nasg_covolume(liq)
      !> Default p_I=p* for the instantaneous mechanical stage: its end state is the pressure equilibrium and the maximum of S
      !> on its path (tex: Choice of p_I); a finite mechanical step (MRR<huge) holds the impedance p_I of its start state
      this%pI_mode=pI_end
   end subroutine initialize

   !> Start state of a relaxation call from the NASG and ideal-gas relations
   subroutine admissible(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(in)  :: this
      real(WP),                               intent(in)  :: VF
      real(WP), dimension(:),                 intent(in)  :: Q
      real(WP),                               intent(in)  :: Pjump
      integer,                                intent(out) :: ierr
      call nasg_admissible(this,this%bL,VF,Q,ierr)
   end subroutine admissible

   !> Liquid state from the NASG relations
   subroutine liq_state(this,st,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(in)    :: this
      type(phasic_state),                          intent(inout) :: st
      integer,                                intent(out)   :: ierr
      call nasg_liq_state(this,this%bL,st,ierr)
   end subroutine liq_state

   !> NASG liquid specific volume at (p,T)
   real(WP) function liq_v_from_p_T(this,p,T) result(v)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(in) :: this
      real(WP),                               intent(in) :: p,T
      v=nasg_liq_v_from_p_T(this,this%bL,p,T)
   end function liq_v_from_p_T

   !> Instantaneous mechanical stage: closed form with p_I=p*, Newton with p_I held at the start state for the other closures
   subroutine relax_M(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(inout) :: this
      real(WP),                               intent(inout) :: VF
      real(WP), dimension(:),                 intent(inout) :: Q
      real(WP),                               intent(in)    :: Pjump
      integer,                                intent(out)   :: ierr
      if (this%pI_mode.eq.pI_end) then
         call nasg_relax_M(this,this%bL,VF,Q,Pjump,ierr)
      else
         call this%iterate(1,VF,Q,Pjump,ierr)
      end if
   end subroutine relax_M

   !> pT-relaxed state in closed form
   subroutine relax_MT(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(inout) :: this
      real(WP),                               intent(inout) :: VF
      real(WP), dimension(:),                 intent(inout) :: Q
      real(WP),                               intent(in)    :: Pjump
      integer,                                intent(out)   :: ierr
      call nasg_relax_MT(this,this%bL,VF,Q,Pjump,ierr)
   end subroutine relax_MT

   !> Pelanti 2022 eqs. 59, 70, 82: V_l and the masses kept, energies split so that p_l-p_g=p_j with E_l+E_g fixed
   subroutine mixture_energy_closure(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_finite), intent(inout) :: this
      real(WP),                          intent(inout) :: VF
      real(WP), dimension(:),            intent(inout) :: Q
      real(WP),                          intent(in)    :: Pjump
      integer,                           intent(out)   :: ierr
      real(WP), dimension(this%gas%ns) :: y
      real(WP) :: gG,qG,Vlb,p
      y=0.0_WP; y(this%indV)=min(Q(this%iVQ())/Q(2),1.0_WP); y(this%indA)=1.0_WP-y(this%indV)
      gG=sum(y*this%gas%cp(1:this%gas%ns))/sum(y*this%gas%cv(1:this%gas%ns))
      qG=sum(y*this%gas%q(1:this%gas%ns))
      Vlb=VF-this%bL*Q(1)
      ierr=RELAX_BAD_LIQUID
      if (Vlb.le.0.0_WP) return
      p=(Q(3)+Q(4)-Q(1)*this%liq%q-Q(2)*qG-Vlb*this%liq%gamma*this%liq%pinf/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)*Pjump/(gG-1.0_WP)) &
      & /(Vlb/(this%liq%gamma-1.0_WP)+(1.0_WP-VF)/(gG-1.0_WP))
      if (p+this%liq%pinf.le.0.0_WP) return
      ierr=RELAX_BAD_GAS
      if (p-Pjump.le.0.0_WP) return
      Q(3)=Vlb*(p+this%liq%gamma*this%liq%pinf)/(this%liq%gamma-1.0_WP)+Q(1)*this%liq%q
      Q(4)=(1.0_WP-VF)*this%gas%get_rhoe_from_p_rho(p=p-Pjump,rho=Q(2)/(1.0_WP-VF),y=y)
      ierr=RELAX_OK
   end subroutine mixture_energy_closure

end module relax_igmix_nasg_max_ent_finite_class
