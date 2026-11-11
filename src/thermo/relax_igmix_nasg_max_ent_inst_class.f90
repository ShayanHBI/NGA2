!> Maximum-entropy relaxation (relaxation_ent_max.tex) for a NASG (or stiffened-gas, b=0) liquid: relax_igmix_max_ent_inst with
!> the closed forms of the mechanical stage with p_I=p* and of the pT-relaxed state (tex: Closed form of the mechanical stage,
!> Closed form of the pT-relaxed state), and the liquid state from the NASG relations
module relax_igmix_nasg_max_ent_inst_class
   use precision,                      only: WP
   use relax_igmix_max_ent_inst_class, only: relax_igmix_max_ent_inst,phasic_state,pI_end
   use stiffened_gas_class,            only: stiffened_gas
   use igmix_class,                    only: igmix
   use thermorelax_class,              only: RELAX_OK,RELAX_FAILED,RELAX_BAD_LIQUID,RELAX_BAD_GAS,RELAX_DEGENERATE
   implicit none
   private

   public :: relax_igmix_nasg_max_ent_inst
   public :: nasg_covolume,nasg_liq_state,nasg_liq_v_from_p_T,nasg_relax_M,nasg_relax_MT,nasg_admissible

   type, extends(relax_igmix_max_ent_inst) :: relax_igmix_nasg_max_ent_inst
      !> Liquid co-volume (0 for a plain stiffened gas); NASG is used without the brhomax clamp of nasg_class
      real(WP) :: bL=0.0_WP
   contains
      procedure :: initialize
      procedure :: admissible
      procedure :: liq_state
      procedure :: liq_v_from_p_T
      procedure :: relax_M
      procedure :: relax_MT
   end type relax_igmix_nasg_max_ent_inst

contains

   !> Parent initialization, then the NASG co-volume and p_I=p* for the instantaneous mechanical stage
   subroutine initialize(this,liq,gas,indV,indA)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(inout) :: this
      class(stiffened_gas), target, intent(in) :: liq
      class(igmix),         target, intent(in) :: gas
      integer,                      intent(in) :: indV,indA
      call this%relax_igmix_max_ent_inst%initialize(liq=liq,gas=gas,indV=indV,indA=indA)
      this%bL=nasg_covolume(liq)
      !> Interface-pressure closure of d_M. Default p_I=p*: the only constant p_I whose end state is both the pressure
      !> equilibrium and the maximum of S on its path; a start-state p_I leaves dS/dV_l/=0 there (tex: Choice of p_I)
      this%pI_mode=pI_end
   end subroutine initialize

   !> Start state of a relaxation call from the NASG and ideal-gas relations
   subroutine admissible(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(in)  :: this
      real(WP),                             intent(in)  :: VF
      real(WP), dimension(:),               intent(in)  :: Q
      real(WP),                             intent(in)  :: Pjump
      integer,                              intent(out) :: ierr
      call nasg_admissible(this,this%bL,VF,Q,ierr)
   end subroutine admissible

   !> Liquid state from the NASG relations
   subroutine liq_state(this,st,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(in)    :: this
      type(phasic_state),                        intent(inout) :: st
      integer,                              intent(out)   :: ierr
      call nasg_liq_state(this,this%bL,st,ierr)
   end subroutine liq_state

   !> NASG liquid specific volume at (p,T)
   real(WP) function liq_v_from_p_T(this,p,T) result(v)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(in) :: this
      real(WP),                             intent(in) :: p,T
      v=nasg_liq_v_from_p_T(this,this%bL,p,T)
   end function liq_v_from_p_T

   !> Mechanical stage: closed form with p_I=p*, Newton with p_I held at the start state for the other closures
   subroutine relax_M(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(inout) :: this
      real(WP),                             intent(inout) :: VF
      real(WP), dimension(:),               intent(inout) :: Q
      real(WP),                             intent(in)    :: Pjump
      integer,                              intent(out)   :: ierr
      if (this%pI_mode.eq.pI_end) then
         call nasg_relax_M(this,this%bL,VF,Q,Pjump,ierr)
      else
         call this%iterate(1,VF,Q,Pjump,ierr)
      end if
   end subroutine relax_M

   !> pT-relaxed state in closed form
   subroutine relax_MT(this,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_nasg_max_ent_inst), intent(inout) :: this
      real(WP),                             intent(inout) :: VF
      real(WP), dimension(:),               intent(inout) :: Q
      real(WP),                             intent(in)    :: Pjump
      integer,                              intent(out)   :: ierr
      call nasg_relax_MT(this,this%bL,VF,Q,Pjump,ierr)
   end subroutine relax_MT

   !> Co-volume b of a NASG liquid, 0 for a stiffened gas
   real(WP) function nasg_covolume(liq) result(b)
      use nasg_class, only: nasg
      implicit none
      class(stiffened_gas), intent(in) :: liq
      b=0.0_WP
      select type (liq)
      type is (nasg)
         b=liq%b
      end select
   end function nasg_covolume

   !> Both phases have a state: masses, 0<VF<1, W0=V_l-b*m_l>0, T_l>0 (E_l-m_l*q>pinf*W0, NASG T=(e-q-pinf*(v-b))/c_v), T_g>0 (E_g>Q_g)
   subroutine nasg_admissible(this,bL,VF,Q,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)  :: this
      real(WP),                        intent(in)  :: bL,VF
      real(WP), dimension(:),          intent(in)  :: Q
      integer,                         intent(out) :: ierr
      real(WP) :: ml,mv,ma,W0
      integer  :: iV
      ierr=RELAX_DEGENERATE
      iV=this%iVQ()
      ml=Q(1); mv=min(Q(iV),Q(2)); ma=Q(2)-mv
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP.or.ml.le.0.0_WP.or.Q(2).le.0.0_WP.or.mv.lt.0.0_WP) return
      ierr=RELAX_BAD_LIQUID
      W0=VF-bL*ml
      if (W0.le.0.0_WP.or.Q(3)-ml*this%liq%q-this%liq%pinf*W0.le.0.0_WP) return
      ierr=RELAX_BAD_GAS
      if (Q(4)-mv*this%gas%q(this%indV)-ma*this%gas%q(this%indA).le.0.0_WP) return
      ierr=RELAX_OK
   end subroutine nasg_admissible

   !> Liquid p,T,s,g,c and Z_l from st%Vl,El,ml
   subroutine nasg_liq_state(this,bL,st,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in)    :: this
      real(WP),                        intent(in)    :: bL
      type(phasic_state),              intent(inout) :: st
      integer,                    intent(out)   :: ierr
      real(WP) :: vl,el,vb,dpdv,dpde,dTdv,dTde
      vl=st%Vl/st%ml; el=st%El/st%ml
      ! NASG p=(gamma-1)(e-q)/(v-b)-gamma*pinf, T=(e-q-pinf*(v-b))/cv, c^2=gamma*(p+pinf)*v^2/(v-b), admissible for v>b
      vb=vl-bL
      if (vb.le.0.0_WP) then; ierr=RELAX_BAD_LIQUID; return; end if
      st%pl=(this%liq%gamma-1.0_WP)*(el-this%liq%q)/vb-this%liq%gamma*this%liq%pinf
      st%Tl=(el-this%liq%q-this%liq%pinf*vb)/this%liq%cv
      if (st%pl+this%liq%pinf.le.0.0_WP.or.st%Tl.le.0.0_WP) then; ierr=RELAX_BAD_LIQUID; return; end if
      st%sl=this%liq%get_s_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%gl=this%liq%get_g_from_p_T(p=st%pl,T=st%Tl,y=[1.0_WP])
      st%cl=vl*sqrt(this%liq%gamma*(st%pl+this%liq%pinf)/vb)
      ! Liquid Z_l from p=(gamma-1)(e-q)/(v-b)-gamma*pinf, T=(e-q-pinf*(v-b))/cv, dg=v*dp-s*dT
      dpdv=-(st%pl+this%liq%gamma*this%liq%pinf)/vb
      dpde=(this%liq%gamma-1.0_WP)/vb
      dTdv=-this%liq%pinf/this%liq%cv
      dTde=1.0_WP/this%liq%cv
      st%Zl(1,:)=[dpdv,dpde,-(vl*dpdv+el*dpde)]/st%ml
      st%Zl(2,:)=[dTdv,dTde,-(vl*dTdv+el*dTde)]/st%ml
      st%Zl(3,:)=vl*st%Zl(1,:)-st%sl*st%Zl(2,:)
      ierr=RELAX_OK
   end subroutine nasg_liq_state

   !> NASG liquid specific volume v=R*T/(p+pinf)+b at (p,T); negative where the liquid has no state (p+pinf<=0, T<=0)
   real(WP) function nasg_liq_v_from_p_T(this,bL,p,T) result(v)
      implicit none
      class(relax_igmix_max_ent_inst), intent(in) :: this
      real(WP),                        intent(in) :: bL,p,T
      v=-1.0_WP
      if (p+this%liq%pinf.le.0.0_WP.or.T.le.0.0_WP) return
      v=this%liq%R*T/(p+this%liq%pinf)+bL
   end function nasg_liq_v_from_p_T

   !> Mechanical stage with p_I=p*, the end pressure (tex: Closed form of the mechanical stage, eqs. mech_gas, mech_liq, mech_quad):
   !> the larger root of eq. mech_quad, unique and admissible for an admissible start
   subroutine nasg_relax_M(this,bL,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                        intent(in)    :: bL
      real(WP),                        intent(inout) :: VF
      real(WP), dimension(:),          intent(inout) :: Q
      real(WP),                        intent(in)    :: Pjump
      integer,                         intent(out)   :: ierr
      real(WP) :: ml,mv,ma,Gg,gl,W0,Vg0,pl0,pg0,a,b,c,sq,z,y,dl
      integer  :: iV
      this%it_M=0
      ierr=RELAX_DEGENERATE
      iV=this%iVQ()
      ml=Q(1); mv=min(Q(iV),Q(2)); ma=Q(2)-mv
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP.or.ml.le.0.0_WP.or.Q(2).le.0.0_WP.or.mv.lt.0.0_WP) return
      ! Start state: NASG liquid (p_l+gamma_l*pinf)*W0=(gamma_l-1)(E_l-m_l*q_l), ideal-gas mixture p_g*V_g=Gamma_g*(E_g-Q_g)
      gl=this%liq%gamma
      W0=VF-bL*ml
      ierr=RELAX_BAD_LIQUID
      if (W0.le.0.0_WP) return
      pl0=(gl-1.0_WP)*(Q(3)-ml*this%liq%q)/W0-gl*this%liq%pinf
      if (pl0+this%liq%pinf.le.0.0_WP) return
      Gg=(mv*this%gas%R(this%indV)+ma*this%gas%R(this%indA))/(mv*this%gas%cv(this%indV)+ma*this%gas%cv(this%indA))
      Vg0=1.0_WP-VF
      pg0=Gg*(Q(4)-mv*this%gas%q(this%indV)-ma*this%gas%q(this%indA))/Vg0
      ierr=RELAX_BAD_GAS
      if (pg0.le.0.0_WP) return
      ierr=RELAX_FAILED
      if (Pjump+this%liq%pinf.le.0.0_WP) return
      ! Eq. mech_quad in z=y-p_g0, y=p*-p_j the gas pressure of the end state
      a=gl*Vg0+(1.0_WP+Gg)*W0
      b=a*pg0+gl*Vg0*(Pjump+this%liq%pinf)-(1.0_WP+Gg)*W0*(pl0-pg0-Pjump)
      c=-(1.0_WP+Gg)*W0*pg0*(pl0-pg0-Pjump)
      sq=sqrt(max(b**2-4.0_WP*a*c,0.0_WP))
      if (b.gt.0.0_WP) then
         z=-2.0_WP*c/(b+sq)
      else
         z=(-b+sq)/(2.0_WP*a)
      end if
      y=pg0+z
      if (y.le.0.0_WP) return
      ! Volume received by the liquid (eq. mech_gas); liquid works at p*=y+p_j, gas at y (eq. dirs_lg)
      dl=Vg0*z/((1.0_WP+Gg)*y)
      ! Pivot K^M at the start with p_I=p* (eq. K_M; NASG numerators of Algorithm step 2)
      this%KM=(pl0+gl*this%liq%pinf+(gl-1.0_WP)*(y+Pjump))/W0+(pg0+Gg*y)/Vg0
      VF=VF+dl
      Q(3)=Q(3)-(y+Pjump)*dl
      Q(4)=Q(4)+y*dl
      ierr=RELAX_OK
   end subroutine nasg_relax_M

   !> pT-relaxed state at fixed masses (tex: Closed form of the pT-relaxed state, eqs. pT_energy, pT_volume, pT_quad):
   !> the positive root y=p_g of eq. pT_quad, unique and admissible for E_s>0, omega>0, p_j>-pinf
   subroutine nasg_relax_MT(this,bL,VF,Q,Pjump,ierr)
      implicit none
      class(relax_igmix_max_ent_inst), intent(inout) :: this
      real(WP),                        intent(in)    :: bL
      real(WP),                        intent(inout) :: VF
      real(WP), dimension(:),          intent(inout) :: Q
      real(WP),                        intent(in)    :: Pjump
      integer,                         intent(out)   :: ierr
      real(WP) :: ml,mv,ma,Cvg,Qg,Al,Ag,Ct,om,pi,Es,Etot,a,b,c,sq,y,Vg,T
      real(WP) :: W0,Vg0,Tl0,Tg0,pl0,pg0,rl,rg,kl,kg
      integer  :: iV
      this%it_T=0
      ierr=RELAX_DEGENERATE
      iV=this%iVQ()
      ml=Q(1); mv=min(Q(iV),Q(2)); ma=Q(2)-mv
      if (VF.le.0.0_WP.or.VF.ge.1.0_WP.or.ml.le.0.0_WP.or.Q(2).le.0.0_WP.or.mv.lt.0.0_WP) return
      ierr=RELAX_BAD_LIQUID
      if (VF-bL*ml.le.0.0_WP) return
      ! Capacities and gas constants at fixed masses
      Cvg=mv*this%gas%cv(this%indV)+ma*this%gas%cv(this%indA)
      Qg =mv*this%gas%q (this%indV)+ma*this%gas%q (this%indA)
      Ag =mv*this%gas%R (this%indV)+ma*this%gas%R (this%indA)
      Al =ml*this%liq%R
      Ct =ml*this%liq%cv+Cvg
      om =1.0_WP-bL*ml
      pi =this%liq%pinf+Pjump
      ! Eq. pT_energy: C*T+pi*W=Es (C is Ct here: Fortran does not distinguish it from the coefficient c)
      Etot=Q(3)+Q(4)+Pjump*VF
      Es=Q(3)+Q(4)-ml*this%liq%q-Qg+Pjump*(VF-bL*ml)
      ierr=RELAX_FAILED
      if (pi.le.0.0_WP.or.Es.le.0.0_WP) return
      ! Pivot K^T at the start (eq. K_T) with Gamma_l/V_l=(gamma-1)/W0, K_s,l/V_l=gamma*(p_l+pinf)/W0, Gamma_g=A_g/C_v,g
      W0=VF-bL*ml; Vg0=1.0_WP-VF
      Tl0=(Q(3)-ml*this%liq%q-this%liq%pinf*W0)/(ml*this%liq%cv)
      Tg0=(Q(4)-Qg)/Cvg
      if (Tl0.gt.0.0_WP.and.Tg0.gt.0.0_WP) then
         pl0=Al*Tl0/W0-this%liq%pinf; pg0=Ag*Tg0/Vg0
         rl=(this%liq%gamma-1.0_WP)/W0; rg=Ag/(Cvg*Vg0)
         kl=this%liq%gamma*(pl0+this%liq%pinf)/W0; kg=(1.0_WP+Ag/Cvg)*pg0/Vg0
         this%KT=1.0_WP/(ml*this%liq%cv)+1.0_WP/Cvg-(rl*Tl0+rg*Tg0)*(rl+rg)/(kl+kg)
      end if
      ! Eq. pT_quad: larger root, without cancellation
      a=Ct*om
      b=om*pi*(Ct+Al)-Es*(Al+Ag)
      c=-Es*Ag*pi
      sq=sqrt(b**2-4.0_WP*a*c)
      if (b.ge.0.0_WP) then
         y=-2.0_WP*c/(b+sq)
      else
         y=(-b+sq)/(2.0_WP*a)
      end if
      ! Gas volume and temperature (eq. pT_volume), gas energy from its EOS, liquid energy from conservation
      Vg=om*Ag*(y+pi)/(Al*y+Ag*(y+pi))
      T=y*Vg/Ag
      VF=1.0_WP-Vg
      Q(4)=Cvg*T+Qg
      Q(3)=Etot-Q(4)-Pjump*VF
      ierr=RELAX_OK
   end subroutine nasg_relax_MT

end module relax_igmix_nasg_max_ent_inst_class
