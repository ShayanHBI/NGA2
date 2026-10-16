import mpmath as mp
mp.mp.dps=40
# NASG water (Pelanti 2022 Table 3), vapor + air ideal-gas mixture
gL,piL,qL,qpL,cvL,bL=mp.mpf('1.187'),mp.mpf('7028e5'),mp.mpf(-1177788),mp.mpf(0),mp.mpf(3610),mp.mpf('6.61e-4')
gV,cvV,qV,qpV=mp.mpf('1.467'),mp.mpf(955),mp.mpf(2077616),mp.mpf(14317)
gA,cvA,qA,qpA=mp.mpf('1.4'),mp.mpf(718),mp.mpf(0),mp.mpf(0)
RV,RA=(gV-1)*cvV,(gA-1)*cvA
# Liquid EOS
lv=lambda p,T:(gL-1)*cvL*T/(p+piL)+bL
le=lambda p,T:cvL*T*(p+gL*piL)/(p+piL)+qL
lh=lambda p,T:gL*cvL*T+bL*p+qL
ls=lambda p,T:cvL*mp.log(T**gL/(p+piL)**(gL-1))+qpL
lg=lambda p,T:lh(p,T)-T*ls(p,T)
def l_pT(V,E,m):
   v,e=V/m,E/m; p=(gL-1)*(e-qL)/(v-bL)-gL*piL; return p,(v-bL)*(p+piL)/(cvL*(gL-1))
# Gas mixture EOS; vapor at its partial pressure
vs=lambda T,pv:cvV*mp.log(T**gV/pv**(gV-1))+qpV
vh=lambda T:gV*cvV*T+qV
vg=lambda T,pv:vh(T)-T*vs(T,pv)
def g_state(V,E,mv,ma):
   T=(E-mv*qV-ma*qA)/(mv*cvV+ma*cvA); return (mv*RV+ma*RA)*T/V,T,mv*RV*T/V
def forces(X,pj):
   Vl,El,ml,mv,ma,Eg=X; pl,Tl=l_pT(Vl,El,ml); pg,Tg,pv=g_state(1-Vl,Eg,mv,ma)
   return pl-pg-pj,Tl-Tg,lg(pl,Tl)-vg(Tg,pv),pl,Tl,pg,Tg,pv
def build(pl,Tl,Tg,Vl,Yv,pj):
   ml=Vl/lv(pl,Tl); pg=pl-pj; Rm=Yv*RV+(1-Yv)*RA; mg=(1-Vl)*pg/(Rm*Tg); mv,ma=Yv*mg,(1-Yv)*mg
   return [Vl,ml*le(pl,Tl),ml,mv,ma,mv*(cvV*Tg+qV)+ma*(cvA*Tg+qA)]
# main.tex response coefficients of each phase at (V,E,m)
def coeffs(X,pj):
   Vl,El,ml,mv,ma,Eg=X; FM,FT,FC,pl,Tl,pg,Tg,pv=forces(X,pj); Vg=1-Vl
   omb=1-bL*ml/Vl
   L=dict(V=Vl,Ks=gL*(pl+piL)/omb,KT=(pl+piL)/omb,ap=omb/Tl,Gam=(gL-1)/omb,Cp=ml*gL*cvL,Cv=ml*cvL,p=pl,T=Tl,m=ml,
          v=Vl/ml,h=lh(pl,Tl),s=ls(pl,Tl),rho=ml/Vl)
   cvm=(mv*cvV+ma*cvA)/(mv+ma); Rm=(mv*RV+ma*RA)/(mv+ma)
   G=dict(V=Vg,Ks=(cvm+Rm)/cvm*pg,KT=pg,ap=1/Tg,Gam=Rm/cvm,Cp=(mv+ma)*(cvm+Rm),Cv=(mv+ma)*cvm,p=pg,T=Tg,rho=(mv+ma)/Vg)
   return L,G
def rel(a,b): return mp.nstr(abs(a-b)/abs(b),3)
h=mp.mpf('1e-12')

for pj in [mp.mpf(0),mp.mpf(3000)]:
   print(f'================ p_j = {mp.nstr(pj,4)} Pa ================')
   # ---------------- Mechanical: fixed masses, dE_l=-p_I dV_l, dE_g=(p_I-p_j) dV_l ----------------
   X=build(mp.mpf('1.3e5'),mp.mpf(360),mp.mpf(380),mp.mpf('0.3'),mp.mpf('0.4'),pj)
   Vl,El,ml,mv,ma,Eg=X; X[5]=Eg+mp.mpf(2e4)    # push gas off mechanical equilibrium
   L,G=coeffs(X,pj); FM=forces(X,pj)[0]
   pI=(L['p']+G['p']+pj)/2; pIl,pIg=pI,pI-pj
   Wle=(L['V']/L['Ks'])/(1-L['Gam']*(L['p']-pIl)/L['Ks']); Wge=(G['V']/G['Ks'])/(1-G['Gam']*(G['p']-pIg)/G['Ks'])
   KM=1/Wle+1/Wge
   xiPl=(L['Gam']*(pIl-L['p'])+L['Ks'])/L['V']; xiPg=(G['Gam']*(pIg-G['p'])+G['Ks'])/G['V']   # Pelanti eq. 51
   def FMof(d):
      Y=list(X); Y[0]+=d; Y[1]-=pIl*d; Y[5]+=pIg*d; return forces(Y,pj)[0]
   dFM=(FMof(h)-FMof(-h))/(2*h)
   print('MECH   F_M/F_T at base       :',mp.nstr(FM,6),mp.nstr(forces(X,pj)[1],6))
   print('MECH   -dF_M/dV_l  FD vs K_M :',rel(KM,-dFM),'   Pelanti 1/xi^p sum vs K_M:',rel(xiPl+xiPg,KM))

   # ---------------- Thermal: fixed masses, on p_l-p_g=p_j, total E_l+E_g+p_j V_l fixed ----------------
   X=build(mp.mpf('1.3e5'),mp.mpf(360),mp.mpf(380),mp.mpf('0.3'),mp.mpf('0.4'),pj)
   L,G=coeffs(X,pj); Etot=X[1]+X[5]+pj*X[0]
   def onM(Vl_):
      def f(El_):
         Y=list(X); Y[0]=Vl_; Y[1]=El_; Y[5]=Etot-pj*Vl_-El_; return forces(Y,pj)[0]
      El_=mp.findroot(f,X[1]); Y=list(X); Y[0]=Vl_; Y[1]=El_; Y[5]=Etot-pj*Vl_-El_; return Y
   Yp,Ym=onM(X[0]+h),onM(X[0]-h)
   dEl=(Yp[1]-Ym[1])/(2*h); dQ=dEl+L['p']            # heat to liquid per dV_l: dQ=dE_l+p_l dV_l
   Fp,Fm=forces(Yp,pj),forces(Ym,pj)
   dFT=(Fp[1]-Fm[1])/(2*h)/dQ; dTl=(Fp[4]-Fm[4])/(2*h)/dQ; dTg=(Fp[6]-Fm[6])/(2*h)/dQ; dp=(Fp[3]-Fm[3])/(2*h)/dQ
   al,ag=L['V']*L['ap'],G['V']*G['ap']; bl,bg=L['V']*L['T']*L['ap'],G['V']*G['T']*G['ap']
   WTl,WTg=L['V']/L['KT'],G['V']/G['KT']; WT=WTl+WTg; Ap,ApT,Cp=al+ag,bl+bg,L['Cp']+G['Cp']
   D=Cp*WT-Ap*ApT; phl=(WT*G['Cp']-ApT*ag)/D; phg=(WT*L['Cp']-ApT*al)/D
   N=L['Cp']*G['Cv']*WTg+G['Cp']*L['Cv']*WTl; KT=D/N
   dVl_main=al*(phl*KT)-WTl*((al*phl*KT-ag*phg*KT)/WT)
   print('THERM  F_T at base           :',mp.nstr(forces(X,pj)[1],6))
   print('THERM  dF_T/dQ   FD vs K_T=D/N:',rel(KT,dFT))
   print('THERM  dT_l/dQ   FD vs phi_l K_T:',rel(phl*KT,dTl),'  dT_g/dQ FD vs -phi_g K_T:',rel(-phg*KT,dTg))
   print('THERM  dV_l/dQ   FD vs main.tex :',rel(dVl_main,1/dQ))
   if pj==0:
      # Pelanti eqs. 10, 61, 62 with c^2=Ks/rho, phi=-rho*ap, zeta=rho/KT
      a1,a2=L['V'],G['V']; r1,r2=L['rho'],G['rho']; Z=a2*L['Gam']+a1*G['Gam']; DP=a1*G['Ks']+a2*L['Ks']
      COF=L['Gam']*G['Ks']-G['Gam']*L['Ks']
      f1,f2=-r1*L['ap'],-r2*G['ap']; z1,z2=r1/L['KT'],r2/G['KT']
      xT1=-(r1/a1*Z+z1*COF)/(f1*DP); xT2=-(r2/a2*Z-z2*COF)/(f2*DP)
      print('THERM  Pelanti 1/xi_1^T vs phi_l K_T:',rel(xT1,phl*KT),'  1/xi_2^T vs phi_g K_T:',rel(xT2,phg*KT))
      print('THERM  Pelanti Z/D vs dV_l/dQ      :',rel(Z/DP,1/dQ))

   # ---------------- Chemical: on p_l-p_g=p_j, T_l=T_g, fixed V, E_l+E_g+p_j V_l, m_l+m_v, m_a ----------------
   for Yv in [mp.mpf('0.4'),mp.mpf(1)]:
      p0,T0=mp.mpf('0.8e5'),mp.mpf(360)
      X=build(p0,T0,T0,mp.mpf('0.3'),Yv,pj)
      L,G=coeffs(X,pj); Vl,El,ml,mv,ma,Eg=X; Etot=El+Eg+pj*Vl; mtot=ml+mv
      def onMT(ml_):
         mv_=mtot-ml_
         def f(p,T):
            Vl_=ml_*lv(p,T); Vg_=(mv_*RV+ma*RA)*T/(p-pj)
            return [Vl_+Vg_-1,(ml_*le(p,T)+mv_*(cvV*T+qV)+ma*(cvA*T+qA)+pj*Vl_-Etot)/Etot]
         p,T=mp.findroot(f,(p0,T0)); pv=mv_*RV*T/(1-ml_*lv(p,T))
         return dict(p=p,T=T,Vl=ml_*lv(p,T),FC=lg(p,T)-vg(T,pv),gl=lg(p,T),mu=vg(T,pv))
      A,B=onMT(ml+h),onMT(ml-h); d=lambda k:(A[k]-B[k])/(2*h)
      FM,FT,FC,pl,Tl,pg,Tg,pv=forces(X,pj)
      al,ag=L['V']*L['ap'],G['V']*G['ap']; bl,bg=L['V']*L['T']*L['ap'],G['V']*G['T']*G['ap']
      WTl,WTg=L['V']/L['KT'],G['V']/G['KT']; WT=WTl+WTg; Ap,ApT,Cp=al+ag,bl+bg,L['Cp']+G['Cp']; D=Cp*WT-Ap*ApT
      xv=mv*RV/(mv*RV+ma*RA); vbar,hbar,sbar=RV*T0/pg,vh(T0),vs(T0,pv)
      Dv,Dh,Ds=L['v']-vbar,L['h']-hbar,L['s']-sbar
      Pi=(Cp*Dv-Ap*Dh)/D; Th=(ApT*Dv-WT*Dh)/D; KCC=RV*T0*(1-xv)/mv if xv<1 else 0
      KC=Dv*Pi-Ds*Th+KCC
      xil=L['v']*Pi-L['s']*Th; xig=-(vbar*Pi-sbar*Th)+KCC; Sa=al*Th-WTl*Pi+L['v']
      print(f'CHEM   Y_v={mp.nstr(Yv,3)}  base F_M,F_T  :',mp.nstr(FM,3),mp.nstr(FT,3))
      print('CHEM   Pi      FD:',rel(Pi,d('p')),'  Theta FD:',rel(Th,d('T')),'  S_alpha FD:',rel(Sa,d('Vl')))
      print('CHEM   K_C     FD:',rel(KC,d('FC')),'  1/xi_l^g=dg_l/dm_l FD:',rel(xil,d('gl')),'  1/xi_g^g=-dmu_v/dm_l FD:',rel(xig,-d('mu')))
      print('CHEM   K_CC/K_C  :',mp.nstr(KCC/KC,4))
