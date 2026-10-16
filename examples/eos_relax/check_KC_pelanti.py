import numpy as np
from scipy.optimize import fsolve
# NASG water, Pelanti 2022 Table 3 (liquid 1, vapor 2)
g1,pi1,eta1,etp1,cv1,b1=1.187,7028e5,-1177788.,0.,3610.,6.61e-4
g2,pi2,eta2,etp2,cv2,b2=1.467,0.,2077616.,14317.,955.,0.
P=[(g1,pi1,eta1,etp1,cv1,b1),(g2,pi2,eta2,etp2,cv2,b2)]
def rho(k,p,T): g,pi,eta,etp,cv,b=P[k]; return (p+pi)/((g-1)*cv*T+b*(p+pi))
def eint(k,p,T): g,pi,eta,etp,cv,b=P[k]; return (p+g*pi)/(p+pi)*cv*T+eta
def s(k,p,T): g,pi,eta,etp,cv,b=P[k]; return cv*np.log(T**g/(p+pi)**(g-1))+etp
def gibbs(k,p,T): g,pi,eta,etp,cv,b=P[k]; return g*cv*T+b*p+eta-T*s(k,p,T)
# state ** : p,T equilibrium, alpha1
p0,T0,a1=0.5e5,360.,0.5
r1,r2=rho(0,p0,T0),rho(1,p0,T0); a2=1-a1
m1,m2=a1*r1,a2*r2; E=m1*eint(0,p0,T0)+m2*eint(1,p0,T0)
# ground truth: move dm liquid->..., solve pT equilibrium at fixed total V=1, E
def pT_at(m1_,m2_):
    f=lambda x:[m1_/rho(0,x[0]*1e5,x[1])+m2_/rho(1,x[0]*1e5,x[1])-1,
                (m1_*eint(0,x[0]*1e5,x[1])+m2_*eint(1,x[0]*1e5,x[1])-E)/E]
    x=fsolve(f,[p0/1e5,T0],xtol=1e-14); return x[0]*1e5,x[1]
h=1e-4*m2
pp,Tp=pT_at(m1+h,m2-h); pm,Tm=pT_at(m1-h,m2+h)
Sp_fd=(pp-pm)/(2*h); ST_fd=(Tp-Tm)/(2*h)
FC=lambda m1_,m2_:(lambda pT:gibbs(0,*pT)-gibbs(1,*pT))(pT_at(m1_,m2_))
KC_fd=(FC(m1+h,m2-h)-FC(m1-h,m2+h))/(2*h)
# Pelanti (23): Gamma, chi, c, phi, zeta per phase
def coeffs(k,p,T):
    g,pi,eta,etp,cv,b=P[k]; r=rho(k,p,T); ombr=1-b*r; Ec=r*eint(k,p,T)
    Gam=(g-1)/ombr; chi=(g-1)*(-eta*ombr+b*(Ec-eta*r))/ombr**2
    c2=g*(p+pi)/(r*ombr); phi=-r*ombr/T
    d=1e-3*(p+pi); zeta=(rho(k,p+d,T)-rho(k,p-d,T))/(2*d)
    return r,Gam,chi,c2,phi,zeta
r1,G1,X1,c1,f1,z1=coeffs(0,p0,T0); r2,G2,X2,c2,f2,z2=coeffs(1,p0,T0)
DT=a1*a2*(r1*c1/G1-r2*c2/G2)*(f1*z2-f2*z1)+(a1/G1+a2/G2)*(a1*f1*r2+a2*f2*r1)
Sp=((X1/G1-X2/G2)*(a1*f1*r2+a2*f2*r1)+(r2*c2/G2-r1*c1/G1)*(a1*f1+a2*f2))/DT
ST=((X2/G2-X1/G1)*(a1*z1*r2+a2*z2*r1)+(r1*c1/G1-r2*c2/G2)*(a1*z1+a2*z2)+(a1/G1+a2/G2)*(r2-r1))/DT
s1,s2=s(0,p0,T0),s(1,p0,T0)
KC_correct=Sp*(1/r1-1/r2)-ST*(s1-s2)
KC_printed=Sp*(1/r1-1/r2)-ST*(s1+s2)   # literal (-1)^(k-1)/rho_k*Sp - s_k*ST, summed
# Your max_ent step_C closed form (Pi,Theta,K_pT; K_CC=0 for pure vapor)
Vl,Vg=a1,a2; apl=(1-b1*r1)/T0; apg=1/T0; KTl=(p0+pi1)/(1-b1*r1); KTg=p0
Cpl=m1*g1*cv1; Cpg=m2*g2*cv2
Ap=Vl*apl+Vg*apg; ApT=Vl*T0*apl+Vg*T0*apg; WT=Vl/KTl+Vg/KTg; Cp=Cpl+Cpg; D=Cp*WT-Ap*ApT
hl=g1*cv1*T0+b1*p0+eta1; hv=g2*cv2*T0+eta2
Dv=1/r1-1/r2; Dh=hl-hv; Ds=s1-s2
Pi=(Cp*Dv-Ap*Dh)/D; Th=(ApT*Dv-WT*Dh)/D; KpT=Dv*Pi-Ds*Th
print(f"S_p : FD {Sp_fd:.6e}  Pelanti(23c) {Sp:.6e}  yours Pi    {Pi:.6e}")
print(f"S_T : FD {ST_fd:.6e}  Pelanti(23b) {ST:.6e}  yours Theta {Th:.6e}")
print(f"K_C : FD {KC_fd:.6e}  Pelanti fixed {KC_correct:.6e}  yours K_pT {KpT:.6e}")
print(f"K_C with eq.(72) read literally (s1+s2): {KC_printed:.6e}")
# Pelanti 2022 eq. (23a) S_alpha against a finite difference of alpha_1 along the pT-equilibrium path
alpha_after=lambda m1_,m2_:(lambda pT:m1_/rho(0,*pT))(pT_at(m1_,m2_))
Sa_fd=(alpha_after(m1+h,m2-h)-alpha_after(m1-h,m2+h))/(2*h)
Sa=((a1/G1+a2/G2)*(a1*f1+a2*f2)+a1*a2*(X1/G1-X2/G2)*(f1*z2-f2*z1))/DT
print(f"S_a : FD {Sa_fd:.6e}  Pelanti(23a) {Sa:.6e}")
