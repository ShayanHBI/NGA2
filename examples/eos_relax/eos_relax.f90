!> Single-cell test for the PTg relaxation model, standalone against the
!> thermo-branch EOS/relaxation library (no AMReX/MPI needed).
!> Reads EOS parameters from input_impact and input_SG, test/sweep conditions
!> from input, then writes saturation-curve CSVs for comparison with
!> IAPWS-IF97 (NGA2_VS_IAPWS.py).
!>
!> make       # build
!> make run   # build + run
!> make clean # remove binary, .mod, .o files

program eos_relax
   use precision,                only: WP
   use stiffened_gas_class,      only: stiffened_gas
   use nasg_class,               only: nasg
   use igmix_class,              only: igmix
   use thermorelax_class,        only: thermorelax
   use relax_igmix_sg_class,     only: relax_igmix_sg,PTgrelax
   use relax_igmix_nasg_class,   only: relax_igmix_nasg
   implicit none

   real(WP), parameter :: Mv=0.0180153_WP, Ma=0.02897_WP   ! molar masses [kg/mol]

   real(WP) :: VF,VF0,Q(8),Q0(8),PL,PG,rhoL,rhoG,TL,TG,Yv

   ! ── What to run, single-cell test and sweep conditions (shared by both blocks)
   logical  :: run_test,run_sweep,run_edge,raw_test
   real(WP) :: p_test,T_test,VF_in_test,Yv_in_test,Pjump_test
   real(WP) :: VF_raw,Q_raw(8)
   real(WP) :: p_sweep,VF0_sweep,Yv0_sweep,Tmin_sweep,Tmax_sweep,Pjump_sweep
   integer  :: nT_sweep

   ! ── Fixed edge-case matrix (mirrors the old cases.txt / edge_cases.txt by hand)
   integer, parameter :: n_edge=7
   character(len=90) :: edge_label(n_edge)
   real(WP) :: edge_p(n_edge),edge_T(n_edge),edge_VF(n_edge),edge_Yv(n_edge)

   !> Post-relax diagnostics for one EOS on one edge case (initial + post-relax state)
   type :: case_result
      real(WP) :: VF0,PL0,PG0,TL0,TG0,Yv0
      real(WP) :: VF,PL,PG,TL,TG,Yv
      real(WP) :: drho_rel,de_rel
   end type case_result
   type(case_result) :: nasg_res(n_edge),sg_res(n_edge)

   call read_logical('input','Run single-cell test',              run_test)
   call read_logical('input','Write saturation curve',            run_sweep)
   call read_logical('input','Write edge cases',                  run_edge)
   call read_real('input','Test pressure',                        p_test)
   call read_real('input','Test temperature',                     T_test)
   call read_real('input','Test initial liquid volume fraction',  VF_in_test)
   call read_real('input','Test initial vapor mass fraction',     Yv_in_test)
   call read_real('input','Test Pjump',                           Pjump_test)
   call read_logical('input','Use raw Q for single-cell test',    raw_test)
   call read_real('input','Raw VF',                                VF_raw)
   call read_real('input','Raw Q1',                                Q_raw(1))
   call read_real('input','Raw Q2',                                Q_raw(2))
   call read_real('input','Raw Q3',                                Q_raw(3))
   call read_real('input','Raw Q4',                                Q_raw(4))
   call read_real('input','Raw Q5',                                Q_raw(5))
   call read_real('input','Raw Q6',                                Q_raw(6))
   call read_real('input','Raw Q7',                                Q_raw(7))
   call read_real('input','Raw Q8',                                Q_raw(8))
   call read_real('input','Sweep pressure',                       p_sweep)
   call read_real('input','Sweep initial liquid volume fraction', VF0_sweep)
   call read_real('input','Sweep initial vapor mass fraction',    Yv0_sweep)
   call read_real('input','Sweep minimum temperature',            Tmin_sweep)
   call read_real('input','Sweep maximum temperature',            Tmax_sweep)
   call read_int ('input','Sweep number of points',               nT_sweep)
   call read_real('input','Sweep Pjump',                          Pjump_sweep)

   edge_label(1)='CASE 1.1  liquid only (VF=1, no vapor, no air)   (p=1e6, T=470, VF_in=1, Yv_in=0)'
   edge_p(1)=1.0e6_WP; edge_T(1)=470.0_WP; edge_VF(1)=1.0_WP; edge_Yv(1)=0.0_WP
   edge_label(2)='CASE 1.2  liquid+air, no vapor (VF interfacial)   (p=1e6, T=400, VF_in=0.7, Yv_in=0)'
   edge_p(2)=1.0e6_WP; edge_T(2)=400.0_WP; edge_VF(2)=0.7_WP; edge_Yv(2)=0.0_WP
   edge_label(3)='CASE 1.3  liquid+vapor, no air (VF interfacial)   (p=1e6, T=470, VF_in=0.5, Yv_in=1)'
   edge_p(3)=1.0e6_WP; edge_T(3)=470.0_WP; edge_VF(3)=0.5_WP; edge_Yv(3)=1.0_WP
   edge_label(4)='CASE 1.4  liquid+vapor+air, general baseline   (p=1e6, T=350, VF_in=0.3, Yv_in=0.6)'
   edge_p(4)=1.0e6_WP; edge_T(4)=350.0_WP; edge_VF(4)=0.3_WP; edge_Yv(4)=0.6_WP
   edge_label(5)='CASE 2.1  vapor only (VF=0, no liquid, no air)   (p=1e6, T=350, VF_in=0, Yv_in=1)'
   edge_p(5)=1.0e6_WP; edge_T(5)=350.0_WP; edge_VF(5)=0.0_WP; edge_Yv(5)=1.0_WP
   edge_label(6)='CASE 2.2  vapor+air, no liquid (VF=0)   (p=1e6, T=350, VF_in=0, Yv_in=0.6)'
   edge_p(6)=1.0e6_WP; edge_T(6)=350.0_WP; edge_VF(6)=0.0_WP; edge_Yv(6)=0.6_WP
   edge_label(7)='CASE 2.3  air only, no vapor, no liquid (VF=0)   (p=1e6, T=350, VF_in=0, Yv_in=0)'
   edge_p(7)=1.0e6_WP; edge_T(7)=350.0_WP; edge_VF(7)=0.0_WP; edge_Yv(7)=0.0_WP

   ! ── NASG configuration ────────────────────────────────────────────────────
   nasg_block: block
      type(nasg),            target :: liq
      type(igmix),           target :: gas
      type(relax_igmix_nasg),target :: rm
      real(WP) :: CvL,GammaL,PinfL,bL,qL,qpL
      real(WP) :: CvV,GammaV,      qV,qpV
      real(WP) :: CvA,GammaA,      qA,qpA
      real(WP) :: p_cav,Tctol

      call read_real('input_impact','Liquid specific heat capacity at constant volume',CvL)
      call read_real('input_impact','Liquid specific heat capacity ratio',             GammaL)
      call read_real('input_impact','Liquid reference energy shift',                   qL)
      call read_real('input_impact','Liquid reference entropy shift',                  qpL)
      call read_real('input_impact','Liquid stiffening pressure',                      PinfL)
      call read_real('input_impact','Liquid co-volume',                                bL)
      call read_real('input_impact','Vapor specific heat capacity at constant volume', CvV)
      call read_real('input_impact','Vapor specific heat capacity ratio',              GammaV)
      call read_real('input_impact','Vapor reference energy shift',                    qV)
      call read_real('input_impact','Vapor reference entropy shift',                   qpV)
      call read_real('input_impact','Air specific heat capacity at constant volume',   CvA)
      call read_real('input_impact','Air specific heat capacity ratio',                GammaA)
      call read_real('input_impact','Air reference energy shift',                      qA)
      call read_real('input_impact','Air reference entropy shift',                     qpA)
      call read_real('input_impact','Cavitation pressure threshold',                   p_cav)
      call read_real('input_impact','Condensation temperature tolerance',              Tctol)
      call liq%initialize(gamma=GammaL,pinf=PinfL,b=bL,cv=CvL,q=qL,qp=qpL,name='water')
      call gas%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,qA],qp=[qpV,qpA], &
      &                    species_names=['vapor','air  '],name='gas')
      call rm%initialize(liq=liq,gas=gas,indV=1,indA=2)
      rm%model=PTgrelax
      rm%p_cav=p_cav
      rm%Tctol=Tctol

      ! Single-cell test
      if (run_test) then
         if (raw_test) then
            VF=VF_raw; Q0=Q_raw
         else
            call make_Q(liq,gas,p=p_test,T=T_test,VF_in=VF_in_test,Yv_in=Yv_in_test,VF=VF,Q=Q0)
         end if
         VF0=VF; VF=VF0; Q=Q0; call rm%apply(dt=1.0_WP,VF=VF,Q=Q,Pjump=Pjump_test)
         call get_thermo(liq,gas,VF0,Q0,PL,PG,rhoL,rhoG,TL,TG,Yv)
         print '(/,A)','── NASG  initial ───────────────────────────────────────────'
         print '(3(A,ES12.4,3X))','VF=',VF0,'pL=',PL,'pG=',PG,'TL=',TL,'TG=',TG,'Yv=',Yv
         call get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv)
         print '(A)',  '── NASG  post-relax ────────────────────────────────────────'
         print '(3(A,ES12.4,3X))','VF=',VF,'pL=',PL,'pG=',PG,'TL=',TL,'TG=',TG,'Yv=',Yv
         print '(2(A,ES8.1,3X))','Δρ/ρ=',(sum(Q(1:2))-sum(Q0(1:2)))/sum(Q0(1:2)),'ΔΕ/E=',(sum(Q(3:4))-sum(Q0(3:4)))/sum(Q0(3:4))
         if (VF.gt.0.0_WP.and.VF.lt.1.0_WP) print '(A,ES12.4,A,ES12.4)','   (PL-PG)=',PL-PG,'   target Pjump=',Pjump_test
      end if

      if (run_sweep) then
         call write_PTg_curve('eos_relax_NASG.csv',liq,gas,rm,p0=p_sweep,VF0=VF0_sweep,Yv0=Yv0_sweep,Tmin=Tmin_sweep,Tmax=Tmax_sweep,nT=nT_sweep,Pjump=Pjump_sweep)
      end if

      if (run_edge) call compute_edge_cases(liq,gas,rm,nasg_res)
   end block nasg_block

   ! ── SG configuration ──────────────────────────────────────────────────────
   sg_block: block
      type(stiffened_gas), target :: liq
      type(igmix),         target :: gas
      type(relax_igmix_sg),target :: rm
      real(WP) :: CvL,GammaL,PinfL,qL,qpL
      real(WP) :: CvV,GammaV,      qV,qpV
      real(WP) :: CvA,GammaA,      qA,qpA
      real(WP) :: p_cav,Tctol

      call read_real('input_SG','Liquid specific heat capacity at constant volume',CvL)
      call read_real('input_SG','Liquid specific heat capacity ratio',             GammaL)
      call read_real('input_SG','Liquid reference energy shift',                   qL)
      call read_real('input_SG','Liquid reference entropy shift',                  qpL)
      call read_real('input_SG','Liquid stiffening pressure',                      PinfL)
      call read_real('input_SG','Vapor specific heat capacity at constant volume', CvV)
      call read_real('input_SG','Vapor specific heat capacity ratio',              GammaV)
      call read_real('input_SG','Vapor reference energy shift',                    qV)
      call read_real('input_SG','Vapor reference entropy shift',                   qpV)
      call read_real('input_SG','Air specific heat capacity at constant volume',   CvA)
      call read_real('input_SG','Air specific heat capacity ratio',                GammaA)
      call read_real('input_SG','Air reference energy shift',                      qA)
      call read_real('input_SG','Air reference entropy shift',                     qpA)
      call read_real('input_SG','Cavitation pressure threshold',                   p_cav)
      call read_real('input_SG','Condensation temperature tolerance',              Tctol)
      call liq%initialize(gamma=GammaL,pinf=PinfL,cv=CvL,q=qL,qp=qpL,name='water')
      call gas%initialize(gamma=[GammaV,GammaA],cv=[CvV,CvA],q=[qV,qA],qp=[qpV,qpA], &
      &                    species_names=['vapor','air  '],name='gas')
      call rm%initialize(liq=liq,gas=gas,indV=1,indA=2)
      rm%model=PTgrelax
      rm%p_cav=p_cav
      rm%Tctol=Tctol

      ! Single-cell test
      if (run_test) then
         if (raw_test) then
            VF=VF_raw; Q0=Q_raw
         else
            call make_Q(liq,gas,p=p_test,T=T_test,VF_in=VF_in_test,Yv_in=Yv_in_test,VF=VF,Q=Q0)
         end if
         VF0=VF; VF=VF0; Q=Q0; call rm%apply(dt=1.0_WP,VF=VF,Q=Q,Pjump=Pjump_test)
         call get_thermo(liq,gas,VF0,Q0,PL,PG,rhoL,rhoG,TL,TG,Yv)
         print '(/,A)','── SG  initial ─────────────────────────────────────────────'
         print '(3(A,ES12.4,3X))','VF=',VF0,'pL=',PL,'pG=',PG,'TL=',TL,'TG=',TG,'Yv=',Yv
         call get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv)
         print '(A)',  '── SG  post-relax ──────────────────────────────────────────'
         print '(3(A,ES12.4,3X))','VF=',VF,'pL=',PL,'pG=',PG,'TL=',TL,'TG=',TG,'Yv=',Yv
         print '(2(A,ES8.1,3X))','Δρ/ρ=',(sum(Q(1:2))-sum(Q0(1:2)))/sum(Q0(1:2)),'ΔΕ/E=',(sum(Q(3:4))-sum(Q0(3:4)))/sum(Q0(3:4))
         if (VF.gt.0.0_WP.and.VF.lt.1.0_WP) print '(A,ES12.4,A,ES12.4)','   (PL-PG)=',PL-PG,'   target Pjump=',Pjump_test
      end if

      if (run_sweep) then
         call write_PTg_curve('eos_relax_SG.csv',liq,gas,rm,p0=p_sweep,VF0=VF0_sweep,Yv0=Yv0_sweep,Tmin=Tmin_sweep,Tmax=Tmax_sweep,nT=nT_sweep,Pjump=Pjump_sweep)
      end if

      if (run_edge) call compute_edge_cases(liq,gas,rm,sg_res)
   end block sg_block

   if (run_edge) call write_edge_cases('edge_cases.txt',edge_label,nasg_res,sg_res)

contains

   ! ── Read a real value from a NGA2-format input file ───────────────────────
   !> Scans for "Key : Value"; dies if the file is missing or the key is absent.
   subroutine read_real(file,key,val)
      character(len=*), intent(in)  :: file,key
      real(WP),         intent(out) :: val
      character(len=512) :: line,kpart,vpart
      integer :: u,ic
      open(newunit=u,file=file,status='old',action='read',err=98)
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
      98 write(*,*) 'ERROR: could not open input file '''//trim(file)//''''
      error stop
      99 write(*,*) 'ERROR: key '''//trim(key)//''' not found in '''//trim(file)//''''
      error stop
   end subroutine read_real

   ! ── Read an integer value from a NGA2-format input file ───────────────────
   !> Scans for "Key : Value"; dies if the file is missing or the key is absent.
   subroutine read_int(file,key,val)
      character(len=*), intent(in)  :: file,key
      integer,           intent(out) :: val
      character(len=512) :: line,kpart,vpart
      integer :: u,ic
      open(newunit=u,file=file,status='old',action='read',err=98)
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
      98 write(*,*) 'ERROR: could not open input file '''//trim(file)//''''
      error stop
      99 write(*,*) 'ERROR: key '''//trim(key)//''' not found in '''//trim(file)//''''
      error stop
   end subroutine read_int

   ! ── Read a logical value from a NGA2-format input file ────────────────────
   !> Scans for "Key : Value"; dies if the file is missing or the key is absent.
   subroutine read_logical(file,key,val)
      character(len=*), intent(in)  :: file,key
      logical,           intent(out) :: val
      character(len=512) :: line,kpart,vpart
      integer :: u,ic
      open(newunit=u,file=file,status='old',action='read',err=98)
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
      98 write(*,*) 'ERROR: could not open input file '''//trim(file)//''''
      error stop
      99 write(*,*) 'ERROR: key '''//trim(key)//''' not found in '''//trim(file)//''''
      error stop
   end subroutine read_logical

   ! ── Build Q from (p,T,VF,Yv) ─────────────────────────────────────────────
   subroutine make_Q(liq,gas,p,T,VF_in,Yv_in,VF,Q)
      class(stiffened_gas), intent(in)  :: liq
      class(igmix),         intent(in)  :: gas
      real(WP),              intent(in)  :: p,T,VF_in,Yv_in
      real(WP),              intent(out) :: VF,Q(8)
      real(WP) :: yL(1),yG(2)
      yL=[1.0_WP]
      yG=[Yv_in,1.0_WP-Yv_in]
      Q=0.0_WP
      Q(1)=VF_in*liq%get_rho_from_p_T(p,T,yL)
      Q(3)=Q(1)*liq%get_e_from_p_T(p,T,yL)
      if (VF_in.eq.1.0_WP) then
         Q(2)=0.0_WP
         Q(4)=0.0_WP
      else
         Q(2)=(1.0_WP-VF_in)*gas%get_rho_from_p_T(p,T,yG)
         Q(4)=Q(2)*gas%get_e_from_p_T(p,T,yG)
      end if
      Q(8)=Q(2)*Yv_in
      VF=VF_in
   end subroutine make_Q

   ! ── Extract primitives from (VF,Q) ───────────────────────────────────────
   subroutine get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv,pV,rhoV,hV)
      class(stiffened_gas), intent(in)  :: liq
      class(igmix),         intent(in)  :: gas
      real(WP),              intent(in)  :: VF,Q(8)
      real(WP),              intent(out) :: PL,PG,rhoL,rhoG,TL,TG,Yv
      real(WP), optional,    intent(out) :: pV,rhoV,hV

      real(WP) :: yL(1),yG(2),yVpure(2),Xv

      yL=[1.0_WP]

      ! Gas composition
      if (Q(2).gt.0.0_WP) then
         Yv=Q(8)/Q(2)
         Yv=min(max(Yv,0.0_WP),1.0_WP)
      else
         Yv=0.0_WP
      end if

      yG=[Yv,1.0_WP-Yv]

      ! Liquid phase
      if (VF.gt.0.0_WP.and.Q(1).gt.0.0_WP) then
         rhoL=Q(1)/VF
         PL  =liq%get_p_from_rho_e(rhoL,Q(3)/Q(1),yL)
         TL  =liq%get_T_from_rho_e(rhoL,Q(3)/Q(1),yL)
      else
         rhoL=0.0_WP
         PL  =0.0_WP
         TL  =0.0_WP
      end if

      ! Gas phase
      if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
         rhoG=Q(2)/(1.0_WP-VF)
         PG  =gas%get_p_from_rho_e(rhoG,Q(4)/Q(2),yG)
         TG  =gas%get_T_from_rho_e(rhoG,Q(4)/Q(2),yG)
      else
         rhoG=0.0_WP
         PG  =0.0_WP
         TG  =0.0_WP
      end if

      ! Optional vapor partial pressure, vapor-only density, and vapor-only enthalpy.
      ! These are meaningful only when the gas phase exists.
      yVpure=[1.0_WP,0.0_WP]
      if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
         Xv=Yv*Ma/(Yv*Ma+(1.0_WP-Yv)*Mv)
      else
         Xv=0.0_WP
      end if

      if (present(pV)) then
         if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
            pV=Xv*PG
         else
            pV=0.0_WP
         end if
      end if

      if (present(rhoV)) then
         if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
            rhoV=gas%get_rho_from_p_T(Xv*PG,TG,yVpure)
         else
            rhoV=0.0_WP
         end if
      end if

      if (present(hV)) then
         if (VF.lt.1.0_WP.and.Q(2).gt.0.0_WP) then
            hV=gas%get_h_from_p_T(Xv*PG,TG,yVpure)
         else
            hV=0.0_WP
         end if
      end if

   end subroutine get_thermo

   ! ── Sweep T and write saturation-curve CSV ────────────────────────────────
   !> Columns: T [K], pV [Pa], Yv [-], VF [-], rhoL [kg/m^3], rhoV [kg/m^3], hV [J/kg]
   subroutine write_PTg_curve(file,liq,gas,rm,p0,VF0,Yv0,Tmin,Tmax,nT,Pjump)
      character(len=*),   intent(in)    :: file
      class(stiffened_gas),intent(in)    :: liq
      class(igmix),        intent(in)    :: gas
      class(thermorelax),  intent(inout) :: rm
      real(WP),            intent(in)    :: p0,VF0,Yv0,Tmin,Tmax,Pjump
      integer,             intent(in)    :: nT

      integer  :: u,it,nwrite
      real(WP) :: T0,VF,Q(8)
      real(WP) :: PL,PG,rhoL,rhoG,TL,TG,Yv,pV,rhoV,hV

      open(newunit=u,file=file,status='replace',action='write')
      write(u,'(A)') 'T,pV,Yv,VF,rhoL,rhoV,hV'

      nwrite=0

      do it=1,nT

         T0=Tmin+real(it-1,WP)*(Tmax-Tmin)/real(max(nT-1,1),WP)

         call make_Q(liq,gas,p0,T0,VF0,Yv0,VF,Q)
         call rm%apply(dt=1.0_WP,VF=VF,Q=Q,Pjump=Pjump)

         ! Keep any genuine two-phase state.
         ! Do not reject small liquid volume fractions: condensation can easily
         ! produce VF ~ 1e-5 while still being physically meaningful.
         if (VF.le.0.0_WP.or.VF.ge.1.0_WP) cycle
         if (Q(1).le.0.0_WP.or.Q(2).le.0.0_WP) cycle

         call get_thermo(liq,gas,VF,Q,PL,PG,rhoL,rhoG,TL,TG,Yv,pV=pV,rhoV=rhoV,hV=hV)

         ! Basic sanity filter for plotting.
         if (TL.le.0.0_WP.or.pV.le.0.0_WP) cycle
         if (rhoL.le.0.0_WP.or.rhoV.le.0.0_WP) cycle

         write(u,'(*(G0.15,:,","))') TL,pV,Yv,VF,rhoL,rhoV,hV
         nwrite=nwrite+1

      end do

      close(u)

      print '(A,A,A,I0,A)','Written: ',file,'  (',nwrite,' rows)'

   end subroutine write_PTg_curve

   ! ── Run the fixed edge-case matrix (edge_label/edge_p/edge_T/edge_VF/edge_Yv,
   !    host-associated from the program) through one EOS/relax model ──────────
   subroutine compute_edge_cases(liq,gas,rm,res)
      class(stiffened_gas), intent(in)    :: liq
      class(igmix),         intent(in)    :: gas
      class(thermorelax),   intent(inout) :: rm
      type(case_result),    intent(out)   :: res(n_edge)
      real(WP) :: VFl,VF0l,Ql(8),Q0l(8),rhoLl,rhoGl
      integer  :: iec
      do iec=1,n_edge
         call make_Q(liq,gas,p=edge_p(iec),T=edge_T(iec),VF_in=edge_VF(iec),Yv_in=edge_Yv(iec),VF=VFl,Q=Q0l)
         VF0l=VFl; Ql=Q0l
         call rm%apply(dt=1.0_WP,VF=VFl,Q=Ql,Pjump=0.0_WP)
         call get_thermo(liq,gas,VF0l,Q0l,res(iec)%PL0,res(iec)%PG0,rhoLl,rhoGl,res(iec)%TL0,res(iec)%TG0,res(iec)%Yv0)
         res(iec)%VF0=VF0l
         call get_thermo(liq,gas,VFl,Ql,res(iec)%PL,res(iec)%PG,rhoLl,rhoGl,res(iec)%TL,res(iec)%TG,res(iec)%Yv)
         res(iec)%VF=VFl
         res(iec)%drho_rel=(sum(Ql(1:2))-sum(Q0l(1:2)))/sum(Q0l(1:2))
         res(iec)%de_rel  =(sum(Ql(3:4))-sum(Q0l(3:4)))/sum(Q0l(3:4))
      end do
   end subroutine compute_edge_cases

   ! ── Write the edge-case matrix (both EOS side by side per case), same format
   !    as the original hand-assembled cases.txt ────────────────────────────────
   subroutine write_edge_cases(file,label,resA,resB)
      character(len=*),  intent(in) :: file
      character(len=*),  intent(in) :: label(:)
      type(case_result), intent(in) :: resA(:),resB(:)
      integer :: u,ic
      open(newunit=u,file=file,status='replace',action='write')
      do ic=1,size(label)
         write(u,'(A)') repeat('#',64)
         write(u,'(A)') trim(label(ic))
         write(u,'(A)') repeat('#',64)
         write(u,'(A)') ''
         write(u,'(A)') '── NASG  initial ───────────────────────────────────────────'
         write(u,'(3(A,ES12.4,3X))') 'VF=',resA(ic)%VF0,'pL=',resA(ic)%PL0,'pG=',resA(ic)%PG0,&
         &                           'TL=',resA(ic)%TL0,'TG=',resA(ic)%TG0,'Yv=',resA(ic)%Yv0
         write(u,'(A)') '── NASG  post-relax ────────────────────────────────────────'
         write(u,'(3(A,ES12.4,3X))') 'VF=',resA(ic)%VF,'pL=',resA(ic)%PL,'pG=',resA(ic)%PG,&
         &                           'TL=',resA(ic)%TL,'TG=',resA(ic)%TG,'Yv=',resA(ic)%Yv
         write(u,'(2(A,ES8.1,3X))') 'Δρ/ρ=',resA(ic)%drho_rel,'ΔΕ/E=',resA(ic)%de_rel
         write(u,'(A)') ''
         write(u,'(A)') '── SG  initial ─────────────────────────────────────────────'
         write(u,'(3(A,ES12.4,3X))') 'VF=',resB(ic)%VF0,'pL=',resB(ic)%PL0,'pG=',resB(ic)%PG0,&
         &                           'TL=',resB(ic)%TL0,'TG=',resB(ic)%TG0,'Yv=',resB(ic)%Yv0
         write(u,'(A)') '── SG  post-relax ──────────────────────────────────────────'
         write(u,'(3(A,ES12.4,3X))') 'VF=',resB(ic)%VF,'pL=',resB(ic)%PL,'pG=',resB(ic)%PG,&
         &                           'TL=',resB(ic)%TL,'TG=',resB(ic)%TG,'Yv=',resB(ic)%Yv
         write(u,'(2(A,ES8.1,3X))') 'Δρ/ρ=',resB(ic)%drho_rel,'ΔΕ/E=',resB(ic)%de_rel
         write(u,'(A)') ''
      end do
      close(u)
      print '(A,A,A,I0,A)','Written: ',file,'  (',size(label),' cases)'
   end subroutine write_edge_cases

end program eos_relax
