"""Half-volume reactor relaxation: split the dominant phase into an inert part a and a part b
that, with the whole minor phase, forms a VF=0.5 reactor relaxed by the solver's own apply().
EOS and relaxation settings are those of examples/cav_tube (relax_batch.f90)."""
import os
import subprocess
import tempfile
import numpy as np
from scipy.optimize import brentq

EXE = os.environ.get('RELAX_BATCH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'relax_batch'))
KIND = os.environ.get('RELAX_KIND', 'inst')
PRELAX, PTGRELAX = 1, 3
VFLO = 1.0e-12

# Pelanti 2022 Table 3 (cav_tube input)
gL, pinfL, bL, cvL, qL, qpL = 1.187, 7.028e8, 6.61e-4, 3610.0, -1177788.0, 0.0
RL = (gL-1.0)*cvL
cpL = gL*cvL
gam = np.array([1.467, 1.4])
cvs = np.array([955.0, 718.0])
qs = np.array([2077616.0, 0.0])
qps = np.array([14317.0, 0.0])
cps = gam*cvs
Rs = cps-cvs
# Saturation-curve coefficients (relax_igmix_sg initialize, NASG ES)
RV = Rs[0]
AS = (cpL-cps[0]+qps[0]-qpL)/RV
BS = (qL-qs[0])/RV
CS = (cps[0]-cpL)/RV
DS = (cpL-cvL)/RV
ES = bL/RV


# ---------------------------------------------------------------- EOS
def liq_rho_pT(p, T):
   return (p+pinfL)/(RL*T+bL*(p+pinfL))


def liq_rhoe_prho(p, rho):
   return (1.0-bL*rho)*(p+gL*pinfL)/(gL-1.0)+rho*qL


def liq_p(rho, e):
   return (gL-1.0)*rho*(e-qL)/(1.0-bL*rho)-gL*pinfL


def liq_T(rho, e):
   return (e-qL-pinfL*(1.0-bL*rho)/rho)/cvL


def liq_s(p, T):
   return cpL*np.log(T)-RL*np.log(p+pinfL)+qpL


def mix(Yv):
   y = np.stack([Yv, 1.0-Yv])
   return (y*cvs[:, None]).sum(0), (y*Rs[:, None]).sum(0), (y*qs[:, None]).sum(0)


def gas_p(rho, e, Yv):
   cvm, Rm, qm = mix(Yv)
   return Rm*rho*(e-qm)/cvm


def gas_T(rho, e, Yv):
   cvm, Rm, qm = mix(Yv)
   return (e-qm)/cvm


def gas_s(p, T, Yv):
   y = np.stack([Yv, 1.0-Yv])
   xR = (y*Rs[:, None]).sum(0)
   s = 0.0
   for n in range(2):
      xn = y[n]*Rs[n]/xR
      s = s+y[n]*(cps[n]*np.log(T)-Rs[n]*np.log(np.maximum(xn*p, 1e-300))+qps[n])
   return s


def pvsat(pl, T):
   return np.exp(AS+(BS+ES*pl)/T)*T**CS*(pl+pinfL)**DS


def psat(T):
   p = 1e4
   for _ in range(60):
      p = pvsat(p, T)
   return p


def Tsat(p):
   return brentq(lambda T: np.log(pvsat(p, T))-np.log(p), 250.0, 600.0)


def xv(Yv):
   Mv, Ma = 0.0180153, 0.02897
   return Yv*Ma/(Yv*Ma+(1.0-Yv)*Mv)


# ---------------------------------------------------------------- cell state <-> Q
def build_Q(VF, pL, TL, pG, TG, Yv):
   """Q(1:8) per unit cell volume from phasic primitives (arrays of equal length)."""
   n = len(VF)
   Q = np.zeros((n, 8))
   rL = liq_rho_pT(pL, TL)
   cvm, Rm, qm = mix(Yv)
   rG = pG/(Rm*TG)
   Q[:, 0] = VF*rL
   Q[:, 1] = (1.0-VF)*rG
   Q[:, 2] = VF*liq_rhoe_prho(pL, rL)
   Q[:, 3] = (1.0-VF)*rG*(cvm*TG+qm)
   Q[:, 7] = Q[:, 1]*Yv
   return Q


def prims(VF, Q):
   """Phasic primitives; NaN where a phase is absent."""
   with np.errstate(all='ignore'):
      hasL = (VF > VFLO) & (Q[:, 0] > 0)
      hasG = (VF < 1.0-VFLO) & (Q[:, 1] > 0)
      rL = np.where(hasL, Q[:, 0]/VF, np.nan)
      eL = np.where(hasL, Q[:, 2]/Q[:, 0], np.nan)
      rG = np.where(hasG, Q[:, 1]/(1.0-VF), np.nan)
      eG = np.where(hasG, Q[:, 3]/Q[:, 1], np.nan)
      Yv = np.where(hasG, np.clip(Q[:, 7]/Q[:, 1], 0.0, 1.0), np.nan)
      out = dict(VF=VF, rhoL=rL, pL=liq_p(rL, eL), TL=liq_T(rL, eL),
                 rhoG=rG, pG=gas_p(rG, eG, Yv), TG=gas_T(rG, eG, Yv), Yv=Yv)
      out['sL'] = liq_s(out['pL'], out['TL'])
      out['sG'] = gas_s(out['pG'], out['TG'], Yv)
      # Total entropy per unit cell volume (absent phase contributes 0)
      out['S'] = np.where(hasL, Q[:, 0]*out['sL'], 0.0)+np.where(hasG, Q[:, 1]*out['sG'], 0.0)
   return out


# ---------------------------------------------------------------- solver call
def apply(VF, Q, model):
   """Run relax_igmix_nasg_finite%apply on every cell; returns (VF, Q, ierr)."""
   n = len(VF)
   models = np.full(n, model) if np.isscalar(model) else np.asarray(model)
   with tempfile.TemporaryDirectory() as d:
      fin, fout = os.path.join(d, 'in.dat'), os.path.join(d, 'out.dat')
      with open(fin, 'w') as f:
         f.write(f'{n}\n')
         for i in range(n):
            f.write(f'{models[i]:d} {VF[i]:.17e} '+' '.join(f'{x:.17e}' for x in Q[i])+'\n')
      subprocess.run([EXE, fin, fout, KIND], check=True, capture_output=True)
      out = np.loadtxt(fout, ndmin=2)
   return out[:, 0], out[:, 1:9], out[:, 9].astype(int)


# ---------------------------------------------------------------- half-volume reactor
LIQ = [0, 2]        # liquid Q components
GAS = [1, 3, 7]     # gas Q components
MOM = [4, 5, 6]


def split_relax(VF, Q, mech_after=False):
   """Relax a VF=0.5 reactor (whole minor phase + equal volume of the dominant phase), then
   add back the inert part a of the dominant phase (isochoric, adiabatic, mass-additive mixing).
   mech_after: follow with a full-cell mechanical relaxation (apply with Prelax)."""
   liq_minor = VF < 0.5
   Vm = np.minimum(VF, 1.0-VF)          # minor-phase volume = part-b volume
   Vd = np.maximum(VF, 1.0-VF)          # dominant-phase volume
   Vr = 2.0*Vm                          # reactor volume
   f = Vm/Vd                            # fraction of the dominant phase sent to the reactor
   # Reactor extensive content (per unit cell volume), split into kept part a and reactor part
   Qa = np.zeros_like(Q)
   Qr = Q.copy()
   for c in LIQ:
      Qa[:, c] = np.where(liq_minor, 0.0, (1.0-f)*Q[:, c])
      Qr[:, c] = Q[:, c]-Qa[:, c]
   for c in GAS:
      Qa[:, c] = np.where(liq_minor, (1.0-f)*Q[:, c], 0.0)
      Qr[:, c] = Q[:, c]-Qa[:, c]
   Qr[:, MOM] = 0.0
   VFr, Qr_n, ierr = apply(np.full(len(VF), 0.5), Qr/Vr[:, None], PTGRELAX)
   # Recombine: everything is additive per unit cell volume
   Qn = Qa+Vr[:, None]*Qr_n
   Qn[:, MOM] = Q[:, MOM]
   VFn = VFr*Vr+np.where(liq_minor, 0.0, 1.0-Vr)
   out = dict(VF=VFn, Q=Qn, ierr=ierr, VFr=VFr, Qr=Qr_n, Vr=Vr)
   if mech_after:
      VFm, Qm, ierm = apply(VFn, Qn, PRELAX)
      out.update(VF=VFm, Q=Qm, ierr_mech=ierm)
   return out


# ---------------------------------------------------------------- piston-coupled reactor (Python solver)
def liq_e_pT(p, T):
   rho = liq_rho_pT(p, T)
   return cvL*T+qL+pinfL*(1.0-bL*rho)/rho, rho


def piston_cell(VF, Q, piston=True):
   """pTg equilibrium of the half-volume reactor with part a as an isentropic piston: common p
   everywhere, T and chemical equilibrium between the minor phase and part b only, part a exchanges
   p dV work only. piston=False with part a removed is the plain full-cell pTg (solver check)."""
   from scipy.optimize import fsolve
   Q = np.asarray(Q, dtype=float)
   liq_minor = VF < 0.5
   Vm, Vd = min(VF, 1.0-VF), max(VF, 1.0-VF)
   f = Vm/Vd if piston else 1.0
   P0 = prims(np.array([VF]), Q[None, :])
   # Extensive inventory (per unit cell volume)
   mL0, EL0 = Q[0], Q[2]
   mG0, EG0, mV0 = Q[1], Q[3], Q[7]
   mA = mG0-mV0
   if liq_minor:
      ma, sa, ra = (1.0-f)*mG0, P0['sG'][0], mV0/mG0
      mLr, ELr, mGr, EGr, mVr = mL0, EL0, f*mG0, f*EG0, f*mV0
      mAr = f*mA
      Ea0 = (1.0-f)*EG0
   else:
      ma, sa = (1.0-f)*mL0, P0['sL'][0]
      mLr, ELr, mGr, EGr, mVr = f*mL0, f*EL0, mG0, EG0, mV0
      mAr = mA
      Ea0 = (1.0-f)*EL0
   Er0 = ELr+EGr
   Va0 = 1.0-Vm*2.0 if piston else 0.0
   Vr0 = 1.0-Va0

   def part_a(p):
      """Isentropic part a at pressure p: (volume, energy)."""
      if ma <= 0.0:
         return 0.0, 0.0
      if liq_minor:
         cvm, Rm, qm = [float(v) for v in mix(np.array([ra]))]
         cpm = cvm+Rm
         # gas_s(p,T)=sa at frozen composition
         T = np.exp((sa-gas_s(np.array([p]), np.array([1.0]), np.array([ra]))[0])/cpm)
         return ma*Rm*T/p, ma*(cvm*T+qm)
      T = np.exp((sa+RL*np.log(p+pinfL)-qpL)/cpL)
      e, rho = liq_e_pT(p, T)
      return ma/rho, ma*e

   def resid(x):
      # Unknowns: ln p, T, evaporated mass as a fraction of the reactor liquid; residuals scaled by the reactor
      lp, T, z = x
      p = np.exp(lp)
      eps = z*mLr
      mL = mLr-eps
      mV = mVr+eps
      Yv = mV/(mV+mAr)
      cvm, Rm, qm = [float(v) for v in mix(np.array([Yv]))]
      eL, rL = liq_e_pT(p, T)
      Va, Ea = part_a(p)
      r1 = (mL/rL+(mV+mAr)*Rm*T/p-Vr0+(Va-Va0))/Vr0
      r2 = (mL*eL+(mV+mAr)*(cvm*T+qm)-Er0+(Ea-Ea0))/abs(Er0)
      r3 = np.log(xv(Yv)*p)-np.log(pvsat(p, T)) if mAr > 0.0 else np.log(p)-np.log(pvsat(p, T))
      return [r1, r2, r3]

   # Several initial guesses; keep the physical root with the smallest residual
   Tm, Tl, Tgs = np.nanmean([P0['TL'][0], P0['TG'][0]]), P0['TL'][0], P0['TG'][0]
   pd = P0['pG'][0] if liq_minor else P0['pL'][0]
   best = None
   for Tg in (Tl, Tm, Tgs):
      for pg in (psat(Tg), pd):
         if not np.isfinite(pg) or pg <= 0.0:
            continue
         sol, info, ier, msg = fsolve(resid, [np.log(pg), Tg, 0.0], full_output=True, xtol=1e-13)
         r = np.linalg.norm(resid(sol))
         phys = (-1.0 < sol[2] < 1.0) and (mVr+sol[2]*mLr > 0.0) and 200.0 < sol[1] < 1000.0
         if phys and (best is None or r < best[1]):
            best = (sol, r)
      if best is not None and best[1] < 1e-10:
         break
   if best is None:
      return dict(VF=VF, Q=Q.copy(), ok=False, p=np.nan, T=np.nan, eps=0.0)
   lp, T, z = best[0]
   ok = best[1] < 1e-8
   eps = z*mLr
   p = np.exp(lp)
   mL, mV = mLr-eps, mVr+eps
   Yv = mV/(mV+mAr)
   eL, rL = liq_e_pT(p, T)
   cvm, Rm, qm = [float(v) for v in mix(np.array([Yv]))]
   Va, Ea = part_a(p)
   # Recombine into one state per phase (additive mixing of parts a and b)
   Qn = Q.copy()
   VLr = mL/rL
   if liq_minor:
      Qn[0], Qn[2] = mL, mL*eL
      Qn[1], Qn[3], Qn[7] = mV+mAr+ma, (mV+mAr)*(cvm*T+qm)+Ea, mV+ma*ra
      VFn = VLr
   else:
      Qn[0], Qn[2] = mL+ma, mL*eL+Ea
      Qn[1], Qn[3], Qn[7] = mV+mAr, (mV+mAr)*(cvm*T+qm), mV
      VFn = VLr+Va
   return dict(VF=VFn, Q=Qn, ok=ok, p=p, T=T, eps=eps)


# ---------------------------------------------------------------- experiments
METHODS = ['full cell', 'split rigid + M', 'split piston']
COLORS = {'full cell': '#eb6834', 'split rigid + M': '#2a78d6', 'split piston': '#1baf7a'}
MARKERS = {'full cell': 'o', 'split rigid + M': 's', 'split piston': '^'}


def run_methods(VF, Q):
   """Post-relaxation primitives and Q for every method."""
   VFf, Qf, _ = apply(VF, Q, PTGRELAX)
   o = split_relax(VF, Q, mech_after=True)
   rr = [piston_cell(VF[i], Q[i]) for i in range(len(VF))]
   Qp = np.array([r['Q'] for r in rr])
   VFp = np.array([r['VF'] for r in rr])
   okp = np.array([r['ok'] for r in rr])
   res = {'full cell': (VFf, Qf, np.ones(len(VF), bool)),
          'split rigid + M': (o['VF'], o['Q'], np.ones(len(VF), bool)),
          'split piston': (VFp, Qp, okp)}
   return {m: dict(P=prims(v, q), Q=q, ok=ok) for m, (v, q, ok) in res.items()}


def scenarios():
   p0 = psat(353.0)
   Ts = 353.0
   return {'S1: liquid +20 K superheat, saturated vapor': (p0, Ts+20, p0, Ts, 1.0),
           'S2: saturated liquid, vapor +50 K superheat': (p0, Ts, p0, Ts+50, 1.0),
           'S3: liquid -20 K subcooled, saturated vapor': (p0, Ts-20, p0, Ts, 1.0),
           'S4: liquid 300 K, dry air 350 K (Yv=0.005), 1 bar': (1e5, 300.0, 1e5, 350.0, 0.005)}


def random_cells(n, seed=1):
   rng = np.random.default_rng(seed)
   minor = 10.0**rng.uniform(-6.0, np.log10(0.5), n)
   VF = np.where(rng.random(n) < 0.5, minor, 1.0-minor)
   pure = rng.random(n) < 0.5
   p0 = 10.0**rng.uniform(4.0, 5.3, n)
   Ts = np.array([Tsat(p) for p in p0])
   TL = np.where(pure, Ts+rng.uniform(-25.0, 25.0, n), rng.uniform(280.0, 360.0, n))
   TG = np.where(pure, Ts+rng.uniform(-10.0, 60.0, n), rng.uniform(280.0, 420.0, n))
   p0 = np.where(pure, p0, 1e5*10.0**rng.uniform(-0.3, 0.5, n))
   Yv = np.where(pure, 1.0, 10.0**rng.uniform(-3.0, -0.5, n))
   pG = p0*(1.0+rng.uniform(-0.1, 0.1, n))
   return VF, build_Q(VF, p0, TL, pG, TG, Yv), dict(pure=pure, TL=TL, TG=TG)


def excursion(P, T0lo, T0hi):
   """Largest distance of any present phase temperature outside the initial [min,max] range."""
   out = np.zeros_like(T0lo)
   for k in ('TL', 'TG'):
      T = P[k]
      e = np.where(np.isfinite(T), np.maximum(0.0, np.maximum(T0lo-T, T-T0hi)), 0.0)
      out = np.maximum(out, e)
   return out


def main(outdir):
   import matplotlib
   matplotlib.use('Agg')
   import matplotlib.pyplot as plt
   from matplotlib.backends.backend_pdf import PdfPages
   plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False,
                        'axes.grid': True, 'grid.color': '#e4e3df', 'grid.linewidth': 0.5,
                        'axes.edgecolor': '#8a8984', 'text.color': '#0b0b0b', 'axes.labelcolor': '#52514e',
                        'xtick.color': '#52514e', 'ytick.color': '#52514e', 'legend.frameon': False})
   os.makedirs(outdir, exist_ok=True)
   pdf = PdfPages(os.path.join(outdir, 'split_reactor.pdf'))

   # VF sweeps
   VF = 1.0/(1.0+np.exp(-np.linspace(np.log(1e-6/(1-1e-6)), np.log((1-1e-6)/1e-6), 41)))
   n = len(VF)
   sc = scenarios()
   fig, ax = plt.subplots(4, 4, figsize=(13, 11), sharex=True)
   for j, (name, (pL, TL, pG, TG, Yv)) in enumerate(sc.items()):
      Q = build_Q(VF, np.full(n, pL), np.full(n, TL), np.full(n, pG), np.full(n, TG), np.full(n, Yv))
      R = run_methods(VF, Q)
      for m in METHODS:
         P, ok = R[m]['P'], R[m]['ok']
         dm = (R[m]['Q'][:, 0]-Q[:, 0])/Q[:, 0]
         kw = dict(color=COLORS[m], marker=MARKERS[m], ms=3, lw=1.5, label=m)
         ax[0, j].plot(VF[ok], P['TL'][ok], **kw)
         ax[1, j].plot(VF[ok], P['TG'][ok], **kw)
         ax[2, j].plot(VF[ok], np.abs(dm[ok]), **kw)
         ax[3, j].plot(VF[ok], (P['Yv'] if Yv < 1 else P['pL']/1e3)[ok], **kw)
      for i in range(4):
         ax[i, j].set_xscale('logit')
         ax[i, j].axhline(np.nan)
      for i, T0 in ((0, TL), (1, TG)):
         ax[i, j].axhline(T0, color='#8a8984', lw=0.8, ls='--')
      ax[0, j].set_title(name, fontsize=8, loc='left')
      ax[2, j].set_yscale('log')
      ax[3, j].set_xlabel('liquid VF')
      ax[3, j].set_ylabel('gas Yv' if Yv < 1 else 'pressure [kPa]')
   ax[0, 0].set_ylabel('liquid T [K]')
   ax[1, 0].set_ylabel('gas T [K]')
   ax[2, 0].set_ylabel('|liquid mass change| / liquid mass')
   ax[0, 0].legend(loc='best')
   fig.suptitle('One relaxation call vs liquid VF (dashed: initial phase temperature; missing points: phase consumed or no root)', x=0.01, ha='left')
   fig.tight_layout()
   pdf.savefig(fig)
   plt.close(fig)

   # Random ensemble
   VF, Q, meta = random_cells(1500)
   P0 = prims(VF, Q)
   R = run_methods(VF, Q)
   lo = np.minimum(meta['TL'], meta['TG'])
   hi = np.maximum(meta['TL'], meta['TG'])
   mminor = np.minimum(VF, 1.0-VF)
   fig, ax = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey='row')
   summary = []
   for k, m in enumerate(METHODS):
      P, ok = R[m]['P'], R[m]['ok']
      vf = P['VF']
      consumed = ok & ((vf < VFLO) | (vf > 1.0-VFLO) | ~np.isfinite(P['TL']) | ~np.isfinite(P['TG']))
      ex = excursion(P, lo, hi)
      liqmin = VF < 0.5
      for row, sel in ((0, liqmin), (1, ~liqmin)):
         s = ok & sel & ~consumed
         ax[row, k].scatter(mminor[s & meta['pure']], np.maximum(ex[s & meta['pure']], 1e-6), s=6, color=COLORS[m], label='pure vapor')
         ax[row, k].scatter(mminor[s & ~meta['pure']], np.maximum(ex[s & ~meta['pure']], 1e-6), s=6, facecolor='none', edgecolor=COLORS[m], label='vapor + air')
         ax[row, k].set_xscale('log')
         ax[row, k].set_yscale('log')
      ax[0, k].set_title(m, loc='left')
      ax[1, k].set_xlabel('minor-phase volume fraction')
      dS = P['S']-P0['S']
      summary.append(f'{m:16s} no root {np.sum(~ok):4d}   minor phase consumed {np.sum(consumed):4d}   '
                     f'new T extremum >0.1 K {np.sum(ok&(ex>0.1)):4d}   >5 K {np.sum(ok&(ex>5)):4d}   '
                     f'max {np.max(ex[ok]):8.1f} K   dS<0 {np.sum(ok&(dS<-1e-9*np.abs(P0["S"]))):4d}')
   ax[0, 0].set_ylabel('cells with minor liquid: worst phase T\noutside initial [min,max] [K]')
   ax[1, 0].set_ylabel('cells with minor gas: worst phase T\noutside initial [min,max] [K]')
   ax[0, 0].legend(loc='upper left')
   fig.suptitle(f'Random ensemble, {len(VF)} cells: new temperature extrema after one call (floored at 1e-6 K)', x=0.01, ha='left')
   fig.tight_layout()
   pdf.savefig(fig)
   plt.close(fig)
   pdf.close()
   with open(os.path.join(outdir, 'split_reactor_summary.txt'), 'w') as f:
      f.write('\n'.join(summary)+'\n')
   print('\n'.join(summary))


if __name__ == '__main__':
   import sys
   import warnings
   warnings.filterwarnings('ignore')
   main(sys.argv[1] if len(sys.argv) > 1 else '.')
