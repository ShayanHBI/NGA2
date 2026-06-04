!> Single-cell test for the PTg relaxation model.
!> Reads EOS parameters directly from input_NASG and input_SG, then writes
!> saturation-curve CSVs for comparison with IAPWS-IF97 (NGA2_VS_IAPWS.py).
!>
!> Compile:
!>   gfortran -O2 \
!>     ../../src/libraries/precision_dp.f90 \
!>     ../../src/eos/eos_class.f90 ../../src/eos/ig_class.f90   \
!>     ../../src/eos/sg_class.f90  ../../src/eos/nasg_class.f90 \
!>     ../../src/eos/mix_class.f90 ../../src/eos/igmix_class.f90 \
!>     ../../src/amrmpsolvers/relax_class.f90 \
!>     ../../src/amrmpsolvers/relax_sg_ig_class.f90 \
!>     ../../src/amrmpsolvers/relax_nasg_ig_class.f90 \
!>     test_relax.f90 -o test_relax
!> Run:  ./test_relax   (must be run from the Sembian_blastwave directory)


!> make -f Makefile.test_relax        # build
!> make -f Makefile.test_relax run    # build + run
!> make -f Makefile.test_relax clean  # remove binary and .mod files


program test_relax
   use precision,           only: WP
   use eos_class,           only: eos
   use mix_class,           only: mix
   use relax_class,         only: relax
   use ig_class,            only: ig
   use igmix_class,         only: igmix
   use sg_class,            only: sg
   use nasg_class,          only: nasg
   use relax_sg_ig_class,   only: relax_sg_ig
   use relax_nasg_ig_class, only: relax_nasg_ig
   implicit none

   real(WP), parameter :: Mv=0.0180153_WP, Ma=0.02897_WP   ! molar masses [kg/mol]

   real(WP) :: VF,VF0,Q(8),Q0(8),PL,PG,rhoL,rhoG,TL,TG,Yv,hG

   ! ── NASG configuration ────────────────────────────────────────────────────
   nasg_block: block
      type(nasg),      target :: liq
      type(ig),        target :: gas_species(2)
      type(igmix),     target :: gas
      type(relax_nasg_ig),target :: rm
      real(WP) :: CvL,GammaL,PinfL,bL,qL,qpL
      real(WP) :: CvV,GammaV,      qV,qpV
      real(WP) :: CvA,GammaA,      qA,qpA

      call read_real('input_NASG','Liquid specific heat capacity at constant volume',CvL,   3610.0_WP)
      call read_real('input_NASG','Liquid specific heat capacity ratio',             GammaL, 1.19_WP)
      call read_real('input_NASG','Liquid reference energy shift',                   qL,    -1177788.0_WP)
      call read_real('input_NASG','Liquid reference entropy shift',                  qpL,   0.0_WP)
      call read_real('input_NASG','Liquid stiffening pressure',                      PinfL,  7.028e8_WP)
      call read_real('input_NASG','Liquid co-volume',                                bL,     6.61e-4_WP)
      call read_real('input_NASG','Vapor specific heat capacity at constant volume', CvV,    955.0_WP)
      call read_real('input_NASG','Vapor specific heat capacity ratio',              GammaV, 1.47_WP)
      call read_real('input_NASG','Vapor reference energy shift',                    qV,     2077616.0_WP)
      call read_real('input_NASG','Vapor reference entropy shift',                   qpV,    14317.0_WP)
      call read_real('input_NASG','Air specific heat capacity at constant volume',   CvA,    718.0_WP)
      call read_real('input_NASG','Air specific heat capacity ratio',                GammaA, 1.40_WP)
      call read_real('input_NASG','Air reference energy shift',                      qA,     0.0_WP)
      call read_real('input_NASG','Air reference entropy shift',                     qpA,    0.0_WP)
      call liq%initialize(pinf=PinfL,b=bL,gamma=GammaL,cv=CvL,q=qL,qp=qpL)
      call gas_species(1)%initialize(gamma=GammaV,cv=CvV,q=qV,qp=qpV)
      call gas_species(2)%initialize(gamma=GammaA,cv=CvA,q=qA,qp=qpA)
      call gas%initialize(ns=2); call gas%set_species(gas_species)
      call rm%initialize(liq=liq,gas=gas,indV=1,indA=2)

      ! Single-cell test
      call make_Q(liq,gas,p=1.0_WP,T=400.0_WP,VF_in=1.0_WP,Yv_in=0.0_WP,VF=VF,Q=Q0)
      VF0=VF; VF=VF0; Q=Q0; call rm%relax_pTg(VF=VF,Q=Q,Pjump=0.0_WP)
      call get_thermo(liq,gas,VF0,Q0,PL,PG,rhoL,rhoG,TL,TG,Yv,hG)
      print '(/,A)','── NASG  initial ───────────────────────────────────────────'
      print '(3(A,ES12.4,3X))','VF=',VF0,'p=',PL,'T=',TL
      call get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv,hG)
      print '(A)',  '── NASG  post-relax ────────────────────────────────────────'
      print '(3(A,ES12.4,3X))','VF=',VF,'pL=',PL,'pG=',PG,'TL=',TL,'TG=',TG,'Yv=',Yv
      print '(2(A,ES8.1,3X))','Δρ/ρ=',(sum(Q(1:2))-sum(Q0(1:2)))/sum(Q0(1:2)),'ΔΕ/E=',(sum(Q(3:4))-sum(Q0(3:4)))/sum(Q0(3:4))

      call write_PTg_curve('test_relax_NASG.csv',liq,gas,rm,p0=1.0e5_WP,VF0=0.5_WP,Yv0=0.5_WP,Tmin=300.0_WP,Tmax=500.0_WP,nT=200)
      ! call write_PTg_curve('test_relax_NASG.csv',liq,gas,rm,p0=-8.0017e+08_WP,VF0=1.0_WP,Yv0=0.0_WP,Tmin=300.0_WP,Tmax=500.0_WP,nT=200)
   end block nasg_block

   ! ── SG configuration ──────────────────────────────────────────────────────
   sg_block: block
      type(sg),      target :: liq
      type(ig),      target :: gas_species(2)
      type(igmix),   target :: gas
      type(relax_sg_ig),target :: rm
      real(WP) :: CvL,GammaL,PinfL,qL,qpL
      real(WP) :: CvV,GammaV,      qV,qpV
      real(WP) :: CvA,GammaA,      qA,qpA

      call read_real('input_SG','Liquid specific heat capacity at constant volume',CvL,   1816.0_WP)
      call read_real('input_SG','Liquid specific heat capacity ratio',             GammaL, 2.35_WP)
      call read_real('input_SG','Liquid reference energy shift',                   qL,    -1167000.0_WP)
      call read_real('input_SG','Liquid reference entropy shift',                  qpL,   0.0_WP)
      call read_real('input_SG','Liquid stiffening pressure',                      PinfL,  1.0e9_WP)
      call read_real('input_SG','Vapor specific heat capacity at constant volume', CvV,    1040.0_WP)
      call read_real('input_SG','Vapor specific heat capacity ratio',              GammaV, 1.43_WP)
      call read_real('input_SG','Vapor reference energy shift',                    qV,     2030000.0_WP)
      call read_real('input_SG','Vapor reference entropy shift',                   qpV,   -23400.0_WP)
      call read_real('input_SG','Air specific heat capacity at constant volume',   CvA,    718.0_WP)
      call read_real('input_SG','Air specific heat capacity ratio',                GammaA, 1.40_WP)
      call read_real('input_SG','Air reference energy shift',                      qA,     0.0_WP)
      call read_real('input_SG','Air reference entropy shift',                     qpA,    0.0_WP)
      call liq%initialize(pinf=PinfL,gamma=GammaL,cv=CvL,q=qL,qp=qpL)
      call gas_species(1)%initialize(gamma=GammaV,cv=CvV,q=qV,qp=qpV)
      call gas_species(2)%initialize(gamma=GammaA,cv=CvA,q=qA,qp=qpA)
      call gas%initialize(ns=2); call gas%set_species(gas_species)
      call rm%initialize(liq=liq,gas=gas,indV=1,indA=2)

      call write_PTg_curve('test_relax_SG.csv',liq,gas,rm,p0=1.0e5_WP,VF0=0.5_WP,Yv0=0.5_WP,Tmin=300.0_WP,Tmax=500.0_WP,nT=200)
   end block sg_block

contains

   ! ── Read a real value from a NGA2-format input file ───────────────────────
   !> Scans for "Key : Value"; returns default if key is absent or file missing.
   subroutine read_real(file,key,val,default)
      character(len=*), intent(in)  :: file,key
      real(WP),         intent(out) :: val
      real(WP),         intent(in)  :: default
      character(len=512) :: line,kpart,vpart
      integer :: u,ic
      val=default
      open(newunit=u,file=file,status='old',action='read',err=99)
      do
         read(u,'(A)',end=99) line
         line=adjustl(line)
         if (line(1:1).eq.'#'.or.len_trim(line).eq.0) cycle
         ic=index(line,':'); if (ic.eq.0) cycle
         kpart=adjustl(line(1:ic-1)); vpart=adjustl(line(ic+1:))
         if (trim(kpart).eq.trim(key)) then
            read(vpart,*,err=99) val; close(u); return
         end if
      end do
      99 close(u)
   end subroutine read_real

   ! ── Build Q from (p,T,VF,Yv) ─────────────────────────────────────────────
   subroutine make_Q(liq,gas,p,T,VF_in,Yv_in,VF,Q)
      class(eos), intent(in)  :: liq
      class(mix), intent(in)  :: gas
      real(WP),   intent(in)  :: p,T,VF_in,Yv_in
      real(WP),   intent(out) :: VF,Q(8)
      real(WP) :: y(2)
      y=[Yv_in,1.0_WP-Yv_in]
      Q=0.0_WP
      Q(1)=VF_in*liq%get_rho_from_p_T(p,T)
      Q(3)=Q(1)*liq%get_e_from_p_T(p,T)
      if (VF_in.eq.1.0_WP) then
         Q(2)=0.0_WP
         Q(4)=0.0_WP
      else
         Q(2)=(1.0_WP-VF_in)*gas%get_rho_from_p_T(p,T,y)
         Q(4)=Q(2)*gas%get_e_from_p_T(p,T,y)
      end if
      Q(8)=Q(2)*Yv_in
      VF=VF_in
   end subroutine make_Q

   ! ── Extract primitives from (VF,Q) ───────────────────────────────────────
   subroutine get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv,hG)
      class(eos), intent(in)  :: liq
      class(mix), intent(in)  :: gas
      real(WP),   intent(in)  :: VF,Q(8)
      real(WP),   intent(out) :: PL,PG,rhoL,rhoG,TL,TG,Yv,hG
      real(WP) :: y(2)
      Yv  =Q(8)/max(Q(2),tiny(1.0_WP)); y=[Yv,1.0_WP-Yv]
      rhoL=Q(1)/max(VF,         tiny(1.0_WP))
      rhoG=Q(2)/max(1.0_WP-VF, tiny(1.0_WP))
      PL  =liq%get_p_from_rho_e(rhoL,Q(3)/max(Q(1),tiny(1.0_WP)))
      PG  =gas%get_p_from_rho_e(rhoG,Q(4)/max(Q(2),tiny(1.0_WP)),y)
      TL  =liq%get_T_from_p_rho(PL,rhoL)
      TG  =gas%get_T_from_p_rho(PG,rhoG,y)
      hG  =gas%get_h_from_p_T(PG,TG,y)
   end subroutine get_thermo

   ! ── Sweep T and write saturation-curve CSV ────────────────────────────────
   !> Columns: T [K], p [Pa], Yv [-], VF [-], rhoL [kg/m^3], rhoG [kg/m^3], hG [J/kg]
   subroutine write_PTg_curve(file,liq,gas,rm,p0,VF0,Yv0,Tmin,Tmax,nT)
      character(len=*),  intent(in)    :: file
      class(eos),        intent(in)    :: liq
      class(mix),        intent(in)    :: gas
      class(relax),      intent(inout) :: rm
      real(WP),          intent(in)    :: p0,VF0,Yv0,Tmin,Tmax
      integer,           intent(in)    :: nT
      integer  :: u,it
      real(WP) :: VF,Q(8),PL,PG,rhoL,rhoG,TL,TG,Yv,hG
      real(WP), parameter :: VFlo=1e-4_WP,VFhi=1.0_WP-1e-4_WP

      open(newunit=u,file=file,status='replace',action='write')
      write(u,'(A)') 'T,p,Yv,VF,rhoL,rhoG,hG'
      do it=1,nT
         call make_Q(liq,gas,p0,Tmin+real(it-1,WP)*(Tmax-Tmin)/real(max(nT-1,1),WP),VF0,Yv0,VF,Q)
         call rm%relax_pTg(VF=VF,Q=Q,Pjump=0.0_WP)
         if (VF.le.VFlo.or.VF.ge.VFhi.or.Q(1).le.0.0_WP.or.Q(2).le.0.0_WP) cycle
         call get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv,hG)
         write(u,'(*(G0.15,:,","))') TL,PL,Yv,VF,rhoL,rhoG,hG
      end do
      close(u)
      print '(A,A)','Written: ',file
   end subroutine write_PTg_curve

end program test_relax
