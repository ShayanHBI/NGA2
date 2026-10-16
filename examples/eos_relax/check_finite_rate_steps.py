import mpmath as mp
mp.mp.dps=40
# NASG water (Pelanti 2022 Table 3), vapor + air ideal-gas mixture
gL,piL,qL,qpL,cvL,bL=mp.mpf('1.187'),mp.mpf('7028e5'),mp.mpf(-1177788),mp.mpf(0),mp.mpf(3610),mp.mpf('6.61e-4')
gV,cvV,qV,qpV=mp.mpf('1.467'),mp.mpf(955),mp.mpf(2077616),mp.mpf(14317)
gA,cvA,qA,qpA=mp.mpf('1.4'),mp.mpf(718),mp.mpf(0),mp.mpf(0)
RV,RA=(gV-1)*cvV,(gA-1)*cvA
pj=mp.mpf(3000)

# ---------------- state S=[Vl,ml,mv,ma,El,Eg] -> primitives and coefficients (tex Sec. "Cell state") ----------------
def prim(S):
   Vl,ml,mv,ma,El,Eg=S; Vg=1-Vl; mg=mv+ma
   rl=ml/Vl; psi=1-bL*rl
   pl=(gL-1)*(El/Vl-rl*qL)/psi-gL*piL; Tl=(pl+piL)*psi/((gL-1)*cvL*rl)
   Rm=(mv*RV+ma*RA)/mg; cvm=(mv*cvV+ma*cvA)/mg; qm=(mv*qV+ma*qA)/mg
   Tg=(Eg-mg*qm)/(mg*cvm); pg=mg/Vg*Rm*Tg; xv=mv*RV/(mg*Rm); pv=xv*pg
   c=dict(Vl=Vl,Vg=Vg,ml=ml,mg=mg,mv=mv,ma=ma,rl=rl,rg=mg/Vg,psi=psi,pl=pl,Tl=Tl,pg=pg,Tg=Tg,xv=xv,pv=pv,Rm=Rm,cvm=cvm,qm=qm)
   c.update(Ksl=gL*(pl+piL)/psi,KTl=(pl+piL)/psi,apl=psi/Tl,Gl=(gL-1)/psi,Cpl=ml*gL*cvL,Cvl=ml*cvL,
            vl=1/rl,hl=gL*cvL*Tl+bL*pl+qL,sl=cvL*mp.log(Tl**gL/(pl+piL)**(gL-1))+qpL)
   c.update(Ksg=(cvm+Rm)/cvm*pg,KTg=pg,apg=1/Tg,Gg=Rm/cvm,Cpg=mg*(cvm+Rm),Cvg=mg*cvm)
   c.update(vbar=RV*Tg/pg,hbar=gV*cvV*Tg+qV,sbar=cvV*mp.log(Tg**gV/pv**(gV-1))+qpV,KCC=RV*Tg*(1-xv)/mv)
   c['gl']=c['hl']-Tl*c['sl']; c['gv']=c['hbar']-Tg*c['sbar']
   c['FM']=pl-pg-pj; c['FT']=Tl-Tg; c['FC']=c['gl']-c['gv']
   return c
def Esurf(S): return S[4]+S[5]+pj*S[0]
# Energy closure: liquid pressure p at (Vl,ml,mv,ma) such that E_l(p)+E_g(p-pj)=Et; returns El,Eg
def close_E(Vl,ml,mv,ma,Et):
   Vg=1-Vl; mg=mv+ma; Rm=(mv*RV+ma*RA)/mg; cvm=(mv*cvV+ma*cvA)/mg; qm=(mv*qV+ma*qA)/mg
   Et=Et-pj*Vl; Al=(Vl-bL*ml)/(gL-1); Ag=Vg*cvm/Rm
   p=(Et-ml*qL-mg*qm-Al*gL*piL+Ag*pj)/(Al+Ag)
   El=Al*(p+gL*piL)+ml*qL; return [Vl,ml,mv,ma,El,Et-El]
# (close_E takes the conserved total E_l+E_g+p_j V_l)
# Exact thermo-mechanical closure at fixed masses and Esurf: solve for (Vl,El)
def close_pT(S,Et):
   Vl,ml,mv,ma=S[:4]
   def f(Vl_,El_):
      c=prim([Vl_,ml,mv,ma,El_,Et-pj*Vl_-El_]); return [c['FM']/1e5,c['FT']]
   Vl_,El_=mp.findroot(f,(S[0],S[4])); return [Vl_,ml,mv,ma,El_,Et-pj*Vl_-El_]
def mk(pl,Tl,Tg,Vl,Yv):
   ml=Vl*((pl+piL)/((gL-1)*cvL*Tl+bL*(pl+piL))); pg=pl-pj; Rm=Yv*RV+(1-Yv)*RA; mg=(1-Vl)*pg/(Rm*Tg)
   mv,ma=Yv*mg,(1-Yv)*mg; El=ml*(cvL*Tl*(pl+gL*piL)/(pl+piL)+qL); Eg=mv*(cvV*Tg+qV)+ma*(cvA*Tg+qA)
   return [Vl,ml,mv,ma,El,Eg]

# ---------------- the three steps as listed in the tex ----------------
def step_M(S,lam,dt):
   c=prim(S); Zl=c['rl']*mp.sqrt(c['Ksl']/c['rl']); Zg=c['rg']*mp.sqrt(c['Ksg']/c['rg'])
   pI=(Zg*c['pl']+Zl*(c['pg']+pj))/(Zl+Zg)
   KM=(c['Ksl']+c['Gl']*(pI-c['pl']))/c['Vl']+(c['Ksg']+c['Gg']*(pI-pj-c['pg']))/c['Vg']
   om=1 if lam==mp.inf else 1-mp.exp(-lam*KM*dt)
   dV=om*c['FM']/KM; Vl,ml,mv,ma,El,Eg=S
   if lam==mp.inf: return close_E(Vl+dV,ml,mv,ma,Esurf(S)),KM
   return [Vl+dV,ml,mv,ma,El-pI*dV,Eg+(pI-pj)*dV],KM
def step_T(S,th,dt):
   c=prim(S); al,ag=c['Vl']*c['apl'],c['Vg']*c['apg']; WTl,WTg=c['Vl']/c['KTl'],c['Vg']/c['KTg']; WT=WTl+WTg
   ApT=c['Vl']*c['Tl']*c['apl']+c['Vg']*c['Tg']*c['apg']; Cp=c['Cpl']+c['Cpg']; D=Cp*WT-(al+ag)*ApT
   phl=(WT*c['Cpg']-ApT*ag)/D; phg=(WT*c['Cpl']-ApT*al)/D
   N=c['Cpl']*c['Cvg']*WTg+c['Cpg']*c['Cvl']*WTl; KT=D/N
   SV=KT*(WTg*al*phl+WTl*ag*phg)/WT
   om=1 if th==mp.inf else 1-mp.exp(-th*KT*dt)
   Vl,ml,mv,ma,El,Eg=S
   return close_E(Vl-om*SV*c['FT']/KT,ml,mv,ma,Esurf(S)),KT,SV
def coefC(c):
   al=c['Vl']*c['apl']; WTl,WTg=c['Vl']/c['KTl'],c['Vg']/c['KTg']; WT=WTl+WTg
   Ap=al+c['Vg']*c['apg']; ApT=c['Vl']*c['Tl']*c['apl']+c['Vg']*c['Tg']*c['apg']; Cp=c['Cpl']+c['Cpg']; D=Cp*WT-Ap*ApT
   Dv,Dh=c['vl']-c['vbar'],c['hl']-c['hbar']
   Pi=(Cp*Dv-Ap*Dh)/D; Th=(ApT*Dv-WT*Dh)/D
   xil=c['vl']*Pi-c['sl']*Th; xig=-(c['vbar']*Pi-c['sbar']*Th)+c['KCC']
   return dict(Pi=Pi,Th=Th,Sa=al*Th-WTl*Pi+c['vl'],xil=xil,xig=xig,KC=xil+xig)
def step_C(S,nu,dt,closure='A',logm=False):
   c=prim(S); k=coefC(c)
   om=1 if nu==mp.inf else 1-mp.exp(-nu*k['KC']*dt)
   if logm:
      du=om*c['FC']/(c['mv']*k['KC']); M=-c['mv']*(mp.exp(du)-1)
   else:
      M=-om*c['FC']/k['KC']
   Vl,ml,mv,ma,El,Eg=S; Vn=Vl+k['Sa']*M
   if closure=='A': return close_E(Vn,ml+M,mv-M,ma,Esurf(S)),k
   if closure=='ATL':
      M=-c['mv']*(mp.exp(om*c['FC']/(c['mv']*k['KC']))-1); Vn=Vl+k['Sa']*M
      return step_T(close_E(Vn,ml+M,mv-M,ma,Esurf(S)),mp.inf,0)[0],k
   if closure=='AT': return step_T(close_E(Vn,ml+M,mv-M,ma,Esurf(S)),mp.inf,0)[0],k
   return close_pT([Vn,ml+M,mv-M,ma,El,Eg],Esurf(S)),k

def rel(a,b): return mp.nstr(abs(a-b)/max(abs(b),mp.mpf('1e-300')),3)

# ---------------- Mechanical ----------------
print('=== mechanical ===')
S0=mk(mp.mpf('1.3e5'),mp.mpf(360),mp.mpf(380),mp.mpf('0.3'),mp.mpf('0.4')); S0[5]+=mp.mpf(2e4)
S1,KM=step_M(S0,mp.inf,0); c=prim(S1)
print('inst  : F_M*=',mp.nstr(c['FM'],3),' dEsurf=',mp.nstr(Esurf(S1)-Esurf(S0),3),' dm=',mp.nstr(S1[1]-S0[1],3))
for x in [mp.mpf('0.1'),mp.mpf(1),mp.mpf(3)]:
   S1,KM=step_M(S0,x/(KM*1),1); c0,c1=prim(S0),prim(S1)
   print(f'lam*K*dt={mp.nstr(x,2)}: F_M*/F_M0={mp.nstr(c1["FM"]/c0["FM"],8)}  exp={mp.nstr(mp.exp(-x),8)}  dEsurf={mp.nstr(Esurf(S1)-Esurf(S0),3)}')

# ---------------- Thermal (from a mechanically relaxed state) ----------------
print('=== thermal ===')
S0=close_E(*mk(mp.mpf('1.3e5'),mp.mpf(360),mp.mpf(380),mp.mpf('0.3'),mp.mpf('0.4'))[:4],Esurf(mk(mp.mpf('1.3e5'),mp.mpf(360),mp.mpf(380),mp.mpf('0.3'),mp.mpf('0.4'))))
c0=prim(S0); print('start : F_M=',mp.nstr(c0['FM'],3),' F_T=',mp.nstr(c0['FT'],6))
# S_V^T closed form against finite difference of V_l per unit heat along the p-manifold
_,KT,SV=step_T(S0,mp.mpf(0),1)
def heat_path(dV):
   S=close_E(S0[0]+dV,*S0[1:4],Esurf(S0)); return S[4]+prim(S0)['pl']*dV,prim(S)['FT']
h=mp.mpf('1e-14'); Qp,Fp=heat_path(h); Qm,Fm=heat_path(-h)
print('S_V^T closed form vs FD:',rel(SV,2*h/(Qp-Qm)),'   K_T vs FD:',rel(KT,(Fp-Fm)/(Qp-Qm)))
S1,KT,SV=step_T(S0,mp.inf,0); c1=prim(S1)
S2,_,_=step_T(S1,mp.inf,0); c2=prim(S2)
print('inst  : F_T after 1 step=',mp.nstr(c1['FT'],4),' after 2=',mp.nstr(c2['FT'],4),' dEsurf=',mp.nstr(Esurf(S1)-Esurf(S0),3))
for x in [mp.mpf('0.1'),mp.mpf(1),mp.mpf(3)]:
   S1,KT,_=step_T(S0,x/KT,1); c1=prim(S1)
   print(f'th*K*dt={mp.nstr(x,2)}: F_T**/F_T*={mp.nstr(c1["FT"]/c0["FT"],8)}  exp={mp.nstr(mp.exp(-x),8)}  F_M**={mp.nstr(c1["FM"],3)}')

# ---------------- Chemical (from a pT-relaxed state) ----------------
print('=== chemical ===')
for Yv in [mp.mpf('0.1'),mp.mpf('0.4'),mp.mpf('0.95')]:
   B=mk(mp.mpf('0.6e5'),mp.mpf(368),mp.mpf(368),mp.mpf('0.3'),Yv); S0=close_pT(B,Esurf(B)); c0=prim(S0); k0=coefC(c0)
   print(f'-- Y_v={mp.nstr(Yv,2)}: F_C**={mp.nstr(c0["FC"],6)} J/kg  K_CC/K_C={mp.nstr(c0["KCC"]/k0["KC"],3)}')
   for cl in ['A','B']:
      S1,_=step_C(S0,mp.inf,0,cl); c1=prim(S1)
      print(f'   inst closure {cl}: F_C={mp.nstr(c1["FC"],4)} F_T={mp.nstr(c1["FT"],3)} F_M={mp.nstr(c1["FM"],3)} '
            f'dm_tot={mp.nstr(S1[1]+S1[2]-S0[1]-S0[2],3)} dEsurf={mp.nstr(Esurf(S1)-Esurf(S0),3)}')
   # finite rate: one frozen step vs fine integration of the parent law m_l'=-nu F_C with exact pT closure each substep
   for x in [mp.mpf('0.1'),mp.mpf(1)]:
      nu=x/k0['KC']; dt=1
      S1,_=step_C(S0,nu,dt,'B'); F1=prim(S1)['FC']
      S1A,_=step_C(S0,nu,dt,'A'); F1A=prim(S1A)['FC']
      S1L,_=step_C(S0,nu,dt,'B',logm=True); F1L=prim(S1L)['FC']
      Sr=list(S0); n=400
      for i in range(n):  # RK2 on m_l with exact pT closure
         def adv(S,dm): return close_pT([S[0],S[1]+dm,S[2]-dm,S[3],S[4],S[5]],Esurf(S))
         k1=-nu*prim(Sr)['FC']; Sm=adv(Sr,k1*dt/n/2); k2=-nu*prim(Sm)['FC']; Sr=adv(Sr,k2*dt/n)
      Fr=prim(Sr)['FC']
      print(f'   nu*K*dt={mp.nstr(x,2)}: F_C/F_C** frozen(B)={mp.nstr(F1/c0["FC"],7)} frozen(A)={mp.nstr(F1A/c0["FC"],7)} '
            f'log-m={mp.nstr(F1L/c0["FC"],7)} exp={mp.nstr(mp.exp(-x),7)} reference={mp.nstr(Fr/c0["FC"],7)}')

print('=== chemical, substepped: n=ceil(omega*|F_C|/(eta*R_v*T)) ===')
def step_C_sub(S,nu,dt,eta):
   c=prim(S); k=coefC(c); om=1-mp.exp(-nu*k['KC']*dt)
   n=max(1,int(mp.ceil(om*abs(c['FC'])/(eta*RV*c['Tl']))))
   for i in range(n): S,_=step_C(S,nu,dt/n,'B')
   return S,n
for Yv in [mp.mpf('0.1'),mp.mpf('0.95')]:
   B=mk(mp.mpf('0.6e5'),mp.mpf(368),mp.mpf(368),mp.mpf('0.3'),Yv); S0=close_pT(B,Esurf(B)); c0=prim(S0); k0=coefC(c0)
   for x in [mp.mpf('0.1'),mp.mpf(1),mp.mpf(10)]:
      nu=x/k0['KC']; dt=1
      Sr=list(S0); n=2000 if x>1 else 400
      for i in range(n):
         def adv(S,dm): return close_pT([S[0],S[1]+dm,S[2]-dm,S[3],S[4],S[5]],Esurf(S))
         k1=-nu*prim(Sr)['FC']; Sm=adv(Sr,k1*dt/n/2); k2=-nu*prim(Sm)['FC']; Sr=adv(Sr,k2*dt/n)
      Fr=prim(Sr)['FC']/c0['FC']
      out=f'Y_v={mp.nstr(Yv,2)} nu*K*dt={mp.nstr(x,3)}: reference={mp.nstr(Fr,6)}'
      for eta in [mp.inf,mp.mpf('0.3'),mp.mpf('0.1')]:
         S1,n1=step_C_sub(S0,nu,dt,eta); out+=f' | eta={mp.nstr(eta,2)}: {mp.nstr(prim(S1)["FC"]/c0["FC"],6)} (n={n1})'
      print(out)

print('=== chemical closure A + one instantaneous thermal step (AT) vs exact pT closure (B) ===')
for Yv in [mp.mpf('0.1'),mp.mpf('0.95')]:
   B=mk(mp.mpf('0.6e5'),mp.mpf(368),mp.mpf(368),mp.mpf('0.3'),Yv); S0=close_pT(B,Esurf(B)); c0=prim(S0); k0=coefC(c0)
   for x in [mp.mpf('0.1'),mp.mpf(1),mp.inf]:
      nu=x/k0['KC'] if x<mp.inf else mp.inf
      SA,_=step_C(S0,nu,1,'A'); SAT,_=step_C(S0,nu,1,'AT'); SB,_=step_C(S0,nu,1,'B')
      cA,cAT,cB=prim(SA),prim(SAT),prim(SB)
      print(f'Y_v={mp.nstr(Yv,2)} x={mp.nstr(x,3)}: F_T after A={mp.nstr(cA["FT"],3)} K, after AT={mp.nstr(cAT["FT"],3)} K;  |V_l(AT)-V_l(B)|={mp.nstr(abs(SAT[0]-SB[0]),3)}  dEsurf(AT)={mp.nstr(Esurf(SAT)-Esurf(S0),3)}')

print('=== nu -> inf: repeated omega=1 steps (closure AT) ===')
for Yv in [mp.mpf('0.02'),mp.mpf('0.1'),mp.mpf('0.95')]:
   for T0 in [mp.mpf(368),mp.mpf(330)]:
      B=mk(mp.mpf('0.6e5'),T0,T0,mp.mpf('0.3'),Yv); S=close_pT(B,Esurf(B)); out=f'Y_v={mp.nstr(Yv,2)} T={mp.nstr(T0,3)}: F_C/(R_v T) ='
      for it in range(7):
         c=prim(S); out+=' '+mp.nstr(c['FC']/(RV*c['Tl']),3)
         S,_=step_C(S,mp.inf,0,'AT')
      print(out)

print('=== log-m update (closure ATL): nu -> inf iteration, and substepped finite rate vs reference ===')
for Yv,T0 in [(mp.mpf('0.95'),mp.mpf(330)),(mp.mpf('0.1'),mp.mpf(368)),(mp.mpf('0.02'),mp.mpf(368))]:
   B=mk(mp.mpf('0.6e5'),T0,T0,mp.mpf('0.3'),Yv); S=close_pT(B,Esurf(B)); out=f'Y_v={mp.nstr(Yv,2)} T={mp.nstr(T0,3)}: F_C/(R_v T) ='
   for it in range(8):
      c=prim(S); out+=' '+mp.nstr(c['FC']/(RV*c['Tl']),3)
      S,_=step_C(S,mp.inf,0,'ATL')
   print(out)
def step_C_subL(S,nu,dt,eta):
   c=prim(S); k=coefC(c); om=1-mp.exp(-nu*k['KC']*dt)
   n=max(1,int(mp.ceil(om*abs(c['FC'])/(eta*RV*c['Tl']))))
   for i in range(n): S,_=step_C(S,nu,dt/n,'ATL')
   return S,n
for Yv,T0 in [(mp.mpf('0.95'),mp.mpf(330)),(mp.mpf('0.1'),mp.mpf(368))]:
   B=mk(mp.mpf('0.6e5'),T0,T0,mp.mpf('0.3'),Yv); S0=close_pT(B,Esurf(B)); c0=prim(S0); k0=coefC(c0)
   for x in [mp.mpf('0.1'),mp.mpf(1),mp.mpf(10)]:
      nu=x/k0['KC']; dt=1
      Sr=list(S0); n=2000 if x>1 else 400
      for i in range(n):
         def adv(S,dm): return close_pT([S[0],S[1]+dm,S[2]-dm,S[3],S[4],S[5]],Esurf(S))
         k1=-nu*prim(Sr)['FC']; Sm=adv(Sr,k1*dt/n/2); k2=-nu*prim(Sm)['FC']; Sr=adv(Sr,k2*dt/n)
      Fr=prim(Sr)['FC']/c0['FC']
      out=f'Y_v={mp.nstr(Yv,2)} T={mp.nstr(T0,3)} nu*K*dt={mp.nstr(x,3)}: reference={mp.nstr(Fr,6)}'
      for eta in [mp.inf,mp.mpf('0.1')]:
         S1,n1=step_C_subL(S0,nu,dt,eta); out+=f' | eta={mp.nstr(eta,2)}: {mp.nstr(prim(S1)["FC"]/c0["FC"],6)} (n={n1})'
      print(out)

print('=== linear update, substepped (closure AT), condensation and evaporation ===')
def step_C_subA(S,nu,dt,eta):
   c=prim(S); k=coefC(c); om=1-mp.exp(-nu*k['KC']*dt)
   n=max(1,int(mp.ceil(om*abs(c['FC'])/(eta*RV*c['Tl'])))); mvmin=S[2]
   for i in range(n):
      S,_=step_C(S,nu,dt/n,'AT'); mvmin=min(mvmin,S[2])
   return S,n,mvmin
for Yv,T0 in [(mp.mpf('0.95'),mp.mpf(330)),(mp.mpf('0.3'),mp.mpf(330)),(mp.mpf('0.1'),mp.mpf(368))]:
   B=mk(mp.mpf('0.6e5'),T0,T0,mp.mpf('0.3'),Yv); S0=close_pT(B,Esurf(B)); c0=prim(S0); k0=coefC(c0)
   print(f'Y_v={mp.nstr(Yv,2)} T={mp.nstr(T0,3)}: F_C/(R_v T)={mp.nstr(c0["FC"]/(RV*T0),3)}')
   for x in [mp.mpf('0.1'),mp.mpf(1),mp.mpf(10)]:
      nu=x/k0['KC']; dt=1
      Sr=list(S0); n=2000 if x>1 else 400
      for i in range(n):
         def adv(S,dm): return close_pT([S[0],S[1]+dm,S[2]-dm,S[3],S[4],S[5]],Esurf(S))
         k1=-nu*prim(Sr)['FC']; Sm=adv(Sr,k1*dt/n/2); k2=-nu*prim(Sm)['FC']; Sr=adv(Sr,k2*dt/n)
      Fr=prim(Sr)['FC']/c0['FC']
      out=f'   nu*K*dt={mp.nstr(x,3)}: reference={mp.nstr(Fr,6)}'
      for eta in [mp.mpf('0.1'),mp.mpf('0.05')]:
         try:
            S1,n1,mvmin=step_C_subA(S0,nu,dt,eta); out+=f' | eta={mp.nstr(eta,2)}: {mp.nstr(prim(S1)["FC"]/c0["FC"],6)} (n={n1}, min m_v/m_v0={mp.nstr(mvmin/S0[2],3)})'
         except Exception as e: out+=f' | eta={mp.nstr(eta,2)}: FAILED ({type(e).__name__})'
      print(out)
