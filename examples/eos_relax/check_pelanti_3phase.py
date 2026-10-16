import numpy as np
from scipy.optimize import fsolve
# Pelanti 2026: SG water (Table 3) liquid=1, vapor=2; air ideal gas=3 (gamma 1.4, kv 718)
P=[dict(g=2.35,pi=1e9,eta=-1167e3,kv=1816.),dict(g=1.43,pi=0.,eta=2030e3,kv=1040.),dict(g=1.4,pi=0.,eta=0.,kv=718.)]
rho =lambda k,p,T:(p+P[k]['pi'])/((P[k]['g']-1)*P[k]['kv']*T)
eint=lambda k,p,T:(p+P[k]['g']*P[k]['pi'])/(p+P[k]['pi'])*P[k]['kv']*T+P[k]['eta']
s   =lambda k,p,T:P[k]['kv']*np.log(T**P[k]['g']/(p+P[k]['pi'])**(P[k]['g']-1))
p0,T0=1e5,354.; a=np.array([0.6,0.3,0.1])
m=np.array([a[k]*rho(k,p0,T0) for k in range(3)]); E=sum(m[k]*eint(k,p0,T0) for k in range(3)); s30=s(2,p0,T0)
def pT(mm):   # all three at one T
    f=lambda x:[sum(mm[k]/rho(k,x[0]*1e5,x[1]) for k in range(3))-1,(sum(mm[k]*eint(k,x[0]*1e5,x[1]) for k in range(3))-E)/E]
    x=fsolve(f,[p0/1e5,T0],xtol=1e-14); p,T=x[0]*1e5,x[1]; return np.array([mm[0]/rho(0,p,T),mm[2]/rho(2,p,T),p,T])
def pTt(mm):  # liquid-vapor at T, gas isentropic at T3
    def f(x):
        p,T,T3=x[0]*1e5,x[1],x[2]
        return [mm[0]/rho(0,p,T)+mm[1]/rho(1,p,T)+mm[2]/rho(2,p,T3)-1,
                (mm[0]*eint(0,p,T)+mm[1]*eint(1,p,T)+mm[2]*eint(2,p,T3)-E)/E, (s(2,p,T3)-s30)/1e3]
    x=fsolve(f,[p0/1e5,T0,T0],xtol=1e-14); p,T,T3=x[0]*1e5,x[1],x[2]; return np.array([mm[0]/rho(0,p,T),mm[2]/rho(2,p,T3),p,T])
h=1e-5*m[1]; dm=np.array([h,-h,0])
fd_pT =(pT(m+dm)-pT(m-dm))/(2*h); fd_pTt=(pTt(m+dm)-pTt(m-dm))/(2*h)
# Coefficients (SG): Gamma=g-1, chi=-(g-1)eta, c2=g(p+pi)/rho, phi=-rho/T, zeta=rho/(p+pi)
r=[rho(k,p0,T0) for k in range(3)]
G=[P[k]['g']-1 for k in range(3)]; X=[-(P[k]['g']-1)*P[k]['eta'] for k in range(3)]
C=[P[k]['g']*(p0+P[k]['pi'])/r[k] for k in range(3)]; F=[-r[k]/T0 for k in range(3)]; Z=[r[k]/(p0+P[k]['pi']) for k in range(3)]
a1,a2,a3=a; r1,r2,r3=r; G1,G2,G3=G; X1,X2,X3=X; c1,c2,c3=C; f1,f2,f3=F; z1,z2,z3=Z   # c_k here = c_k^2
DT=(r3*a1*a2*(r1*c1/G1-r2*c2/G2)*(f1*z2-f2*z1)+r3*(a1/G1+a2/G2)*(a1*f1*r2+a2*f2*r1)
   +r2*a1*a3*(r1*c1/G1-r3*c3/G3)*(f1*z3-f3*z1)+r2*(a1/G1+a3/G3)*(a1*f1*r3+a3*f3*r1)
   +r1*a3*a2*(r3*c3/G3-r2*c2/G2)*(f3*z2-f2*z3)+r1*(a3/G3+a2/G2)*(a3*f3*r2+a2*f2*r3)
   -r1*r2*a3**2*f3/G3-r1*r3*a2**2*f2/G2-r2*r3*a1**2*f1/G1)
Sa1=((a1/G1+a2/G2+a3/G3)*(a1*f1*r3+a2*f2*r3+a3*f3*r2)+a1*a2*r3*(X1/G1-X2/G2)*(z2*f1-f2*z1)
    +a1*a3*r2*(X1/G1-X3/G3)*(z3*f1-f3*z1)+a1*a3*r2*(X3/G3-X2/G2)*(z3*f1-f3*z1)
    +a3*(r2*c2/G2-r3*c3/G3)*(a1*(z3*f1-f3*z1)+a2*(f2*z3-f3*z2)))/DT
Sa3=(a3**2*f3/G3*(r1-r2)+(a1/G1+a2/G2)*(-a3*f3*r2+a3*f3*r1)+a3*a2*r1*(X1/G1-X2/G2)*(z2*f3-f2*z3)
    -a3*a1*r2*(X2/G2-X1/G1)*(-f1*z3+f3*z1)+a3*(r2*c2/G2-r1*c1/G1)*(a1*(-f1*z3+f3*z1)+a2*(-f2*z3+f3*z2)))/DT
Sp=(r3*(X1/G1-X2/G2)*(a1*f1*r2+a2*f2*r1)+r3*(a1*f1+a2*f2)*(r2*c2/G2-r1*c1/G1)-a3*f3*r3*c3/G3*(r1-r2)
   +a3*f3*r1*r2*((c2-X2)/G2-(c1-X1)/G1))/DT
ST=(r3*(a1*z1*r2+a2*z2*r1)*(X2/G2-X1/G1)+r3*(r1*c1/G1-r2*c2/G2)*(a1*z1+a2*z2)+r3*(a1/G1+a2/G2)*(r2-r1)
   +r3*a3/G3*(c3*z3-1)*(r1-r2)-r1*r2*a3*z3*((c2-X2)/G2-(c1-X1)/G1))/DT
print("pT model  (eqs 16-17)   [S_a1, S_a3, S_p, S_T]")
print("  formula:",np.array([Sa1,Sa3,Sp,ST])); print("  FD     :",fd_pT)
DTt=(f1*r2*a1/G1*(a1*c3*r3+a3*r1*c1)+f2*r1*a2/G2*(a2*c3*r3+a3*c2*r2)
     +a1*a2*r3*c3*(((-c2*z2+1)*r2/G2+r1/G1*z2*c1)*f1+(r2/G2*z1*c2+r1/G1*(-c1*z1+1))*f2))
Sta1=(a1**2/G1*c3*f1*r3+(-(((X2*z2-1)/G2-z2*X1/G1)*f1-(z1*X2/G2-(X1*z1-1)/G1)*f2)*c3*r3*a2
      +r2*f1*a3*((c2-X2)/G2+X1/G1))*a1+a2/G2*f2*(a2*c3*r3+a3*c2*r2))/DTt
Sta3_pr=((r1*c1/G1-((c2-X2)/G2+X1/G1)*r2)*a3*f1*a1+a2*f2*a3*((X2/G2+(c1-X1))*r1/G1-r2*c2/G2))/DTt
Sta3_sy=((r1*c1/G1-((c2-X2)/G2+X1/G1)*r2)*a3*f1*a1+a2*f2*a3*((X2/G2+(c1-X1)/G1)*r1-r2*c2/G2))/DTt
Stp_pr=(-(f1*(r1*c1/G1-((c2-X2)/G2+X1/G1)*r2)*a1+a2*f2*((X2/G2+(c1-X1))*r1/G1-r2*c2/G2))*c3*r3)/DTt
Stp_sy=(-(f1*(r1*c1/G1-((c2-X2)/G2+X1/G1)*r2)*a1+a2*f2*((X2/G2+(c1-X1)/G1)*r1-r2*c2/G2))*c3*r3)/DTt
StT=(c3*((r1/G1*(c1*z1-1)-r2*(z1*(c2-X2)/G2+(X1*z1-1)/G1))*a1+a2*(((X2*z2-1)/G2+z2*(c1-X1))*r1/G1
     -r2/G2*(c2*z2-1)))*r3-((c2-X2)/G2-(c1-X1)/G1)*r2*a3*r1)/DTt
StT_sy=(c3*((r1/G1*(c1*z1-1)-r2*(z1*(c2-X2)/G2+(X1*z1-1)/G1))*a1+a2*(((X2*z2-1)/G2+z2*(c1-X1)/G1)*r1
     -r2/G2*(c2*z2-1)))*r3-((c2-X2)/G2-(c1-X1)/G1)*r2*a3*r1)/DTt
print("pT~ model (eqs 24-25)   [S_a1, S_a3, S_p, S_T]")
print("  printed  :",np.array([Sta1,Sta3_pr,Stp_pr,StT]))
print("  symmetric:",np.array([Sta1,Sta3_sy,Stp_sy,StT_sy]))
print("  FD       :",fd_pTt)
