#!/usr/bin/env python
"""Analyze relax_clustered's solo-vs-clustered diagnostic CSV.

Answers two questions from one single-timestep data-collection run:
  1. For each clustered cell, did clustering actually change the outcome (was it needed)?
  2. What input-state quantity predicts a relaxation overshoot -- VF, mass, energy,
     heat capacity, or pressure?

Usage:  python analyze_relaxdiag.py [csv] [plotfile]
"""
import sys
import numpy as np
import pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else 'impact_relax_pTg_datacollect_relaxdiag_rank0.csv'
PLT = sys.argv[2] if len(sys.argv) > 2 else 'amrviz/impact_relax_pTg_datacollect/plt.nga2.cell.000002'


def read_input_consts(path='input'):
    """cvL, cvV, cvA and gammaA from the case input (cvA from the T2=1 shock normalization)."""
    vals = {}
    for line in open(path):
        line = line.split('#')[0]
        if ':' not in line:
            continue
        k, v = line.split(':', 1)
        try:
            vals[k.strip()] = float(v.strip())
        except ValueError:
            pass
    gammaA = vals['GammaA']
    M2 = vals['Gas Mach number']
    pG2 = 1.0 / (gammaA * M2 ** 2)
    return dict(cvL=vals['Liquid cv'], cvV=vals['Vapor cv'], cvA=pG2 / (gammaA - 1.0))


def spearman(x, y):
    """Rank correlation, no scipy dependency."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 10:
        return np.nan
    rx = pd.Series(x[ok]).rank().to_numpy()
    ry = pd.Series(y[ok]).rank().to_numpy()
    rx -= rx.mean()
    ry -= ry.mean()
    den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / den) if den > 0 else np.nan


df = pd.read_csv(CSV, skipinitialspace=True)
c = read_input_consts()

# ---- Derived input-state quantities (the candidate predictors) ----
df['cvG0'] = df.Yv0 * c['cvV'] + (1.0 - df.Yv0) * c['cvA']
df['CL0'] = df.Q1_0 * c['cvL']                    # liquid heat capacity in the cell
df['CG0'] = df.Q2_0 * df.cvG0                     # gas    heat capacity in the cell
df['fG'] = df.CG0 / (df.CL0 + df.CG0)             # gas share of pool heat capacity
df['fmin'] = np.minimum(df.CL0, df.CG0) / (df.CL0 + df.CG0)   # minority-phase capacity share
df['dT0'] = df.TG0 - df.TL0                       # thermal driving gap
df['dP0'] = df.PL0 - df.PG0                       # mechanical driving gap
df['mass_frac_G'] = df.Q2_0 / (df.Q1_0 + df.Q2_0)
df['E_frac_G'] = df.Q4_0 / (df.Q3_0 + df.Q4_0)
df['pred_dTL'] = df.fG * df.dT0                   # closed-form single-cell prediction

# ---- Outcomes ----
df['dTL_solo'] = df.TLs - df.TL0                  # what relax alone does to this cell
df['dTL_act'] = df.TLa - df.TL0                   # what actually happened
df['dRHOG_solo'] = df.RHOGs - df.RHOG0
df['dRHOG_act'] = df.RHOGa - df.RHOG0
df['eff_TL'] = df.TLa - df.TLs                    # how much clustering changed the result
df['eff_RHOG'] = df.RHOGa - df.RHOGs
df['rel_TL_solo'] = df.dTL_solo / df.TL0.replace(0, np.nan)
df['rel_eff_TL'] = df.eff_TL / df.TL0.replace(0, np.nan)

clus = df[df.event.isin(['HOST', 'GUEST'])].copy()
strand = df[df.event == 'STRANDED'].copy()

print('=' * 78)
print(f'{len(clus)} clustered cells ({len(clus)//2} pairs), {len(strand)} stranded cells')
print('=' * 78)

# ---- Q1: was clustering needed? ----
print('\n--- Q1. Did clustering change the outcome? (|TLa-TLs| relative to TL0) ---')
r = clus.rel_eff_TL.abs()
for thr in [1e-6, 1e-4, 1e-3, 1e-2]:
    print(f'  |dTL| from clustering > {thr:8.0e} of TL0 : {(r > thr).sum():5d}  ({100*(r>thr).mean():5.1f}%)')
print(f'  median |relative TL change from clustering| : {r.median():.3e}')
rg = (clus.eff_RHOG / clus.RHOG0.replace(0, np.nan)).abs()
print(f'  median |relative RHOG change from clustering|: {rg.median():.3e}')
print(f'  RHOG made WORSE by clustering (|RHOGa-ref| > |RHOGs-ref|, ref=partner RHOG0): '
      f'{(clus.RHOGa.sub(clus.RHOG0p).abs() > clus.RHOGs.sub(clus.RHOG0p).abs()).mean()*100:.1f}%')

# ---- Q2: what predicts the size of the solo relaxation excursion? ----
print('\n--- Q2. Spearman rank correlation vs |relative TL jump under solo relax| ---')
target = df.rel_TL_solo.abs().to_numpy()
preds = ['VF0', 'Q1_0', 'Q2_0', 'Q3_0', 'Q4_0', 'RHOG0', 'RHOL0',
         'CL0', 'CG0', 'fG', 'fmin', 'mass_frac_G', 'E_frac_G', 'dT0', 'dP0', 'PL0', 'TL0']
rows = [(p, spearman(df[p].to_numpy(dtype=float), target)) for p in preds]
for p, s in sorted(rows, key=lambda t: -abs(t[1] if np.isfinite(t[1]) else 0)):
    print(f'  {p:14s} {s:+.3f}')

print('\n--- Q2b. Same, but vs |relative RHOG jump under solo relax| ---')
target_g = (df.dRHOG_solo / df.RHOG0.replace(0, np.nan)).abs().to_numpy()
rows = [(p, spearman(df[p].to_numpy(dtype=float), target_g)) for p in preds]
for p, s in sorted(rows, key=lambda t: -abs(t[1] if np.isfinite(t[1]) else 0)):
    print(f'  {p:14s} {s:+.3f}')

# ---- Q2c: VF and capacity share are collinear -- which one carries the information? ----
def partial(primary, secondary, target, nbins=8):
    """Spearman(secondary,target) computed WITHIN narrow bins of primary, so any
    correlation that survives is information the primary does not already carry."""
    d = df[[primary, secondary]].copy()
    d['t'] = target
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    try:
        d['bin'] = pd.qcut(d[primary], nbins, duplicates='drop')
    except ValueError:
        return []
    out = []
    for b, g in d.groupby('bin', observed=True):
        if len(g) >= 30:
            out.append(spearman(g[secondary].to_numpy(), g.t.to_numpy()))
    return [s for s in out if np.isfinite(s)]


print('\n--- Q2c. Partial check: does capacity share predict beyond VF (and vice versa)? ---')
for tname, tgt in [('|dTL/TL0|', target), ('|dRHOG/RHOG0|', target_g)]:
    a = partial('VF0', 'fG', tgt)
    b = partial('fG', 'VF0', tgt)
    fa = f'{np.median(a):+.3f} [{min(a):+.3f},{max(a):+.3f}]' if a else 'n/a'
    fb = f'{np.median(b):+.3f} [{min(b):+.3f},{max(b):+.3f}]' if b else 'n/a'
    print(f'  {tname:14s} fG within VF0 bins: {fa:28s} VF0 within fG bins: {fb}')

# ---- Q3: does the closed-form capacity-share prediction hold? ----
print('\n--- Q3. Closed-form check:  dTL ?= [CG/(CL+CG)]*(TG0-TL0)  under solo relax ---')
m = np.isfinite(df.pred_dTL) & np.isfinite(df.dTL_solo) & (df.dTL_solo.abs() > 1e-12)
print(f'  Spearman(pred, actual) = {spearman(df.pred_dTL[m].to_numpy(), df.dTL_solo[m].to_numpy()):+.3f}   (n={m.sum()})')
ratio = (df.dTL_solo[m] / df.pred_dTL[m]).replace([np.inf, -np.inf], np.nan).dropna()
print(f'  actual/predicted: median {ratio.median():.3f}, IQR [{ratio.quantile(.25):.3f}, {ratio.quantile(.75):.3f}]')

# ---- Q4: are stranded cells actually the problematic ones? ----
print('\n--- Q4. Stranded vs clustered populations (solo-relax excursion) ---')
for name, g in [('clustered', clus), ('stranded', strand)]:
    a = g.rel_TL_solo.abs()
    b = (g.dRHOG_solo / g.RHOG0.replace(0, np.nan)).abs()
    print(f'  {name:10s} |dTL/TL0|: med {a.median():.2e} p90 {a.quantile(.9):.2e} max {a.max():.2e} | '
          f'|dRHOG/RHOG0|: med {b.median():.2e} p90 {b.quantile(.9):.2e}')
print(f'  stranded VF0: min {strand.VF0.min():.4f} med {strand.VF0.median():.4f} max {strand.VF0.max():.4f}')
print(f'  stranded fG : med {strand.fG.median():.3e}   clustered fG: med {clus.fG.median():.3e}')

# ---- Q5: neighbor-relative outlier score from the end-of-step plotfile ----
print('\n--- Q5. Spatial outlier score at end of step (automates the ParaView eyeball) ---')
try:
    import yt
    yt.set_log_level(50)
    ds = yt.load(PLT)
    lvl = ds.max_level
    dims = ds.domain_dimensions * ds.refine_by ** lvl
    cg = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge, dims=dims)
    fields = {f: np.array(cg[('boxlib', f)][:, :, 0]) for f in ('RHOG', 'TL', 'VF', 'stranded')}

    def outlier(field, ii, jj):
        """value / median(neighbours in which the phase actually exists).

        A neighbour where the phase is absent reads exactly 0 (pure-gas cell has TL=0,
        pure-liquid cell has RHOG=0); including those in the median drags it toward 0 and
        manufactures a ~2x ratio out of nothing, so they are masked out.
        """
        nx, ny = field.shape
        ii = np.clip(ii, 1, nx - 2)
        jj = np.clip(jj, 1, ny - 2)
        nb = np.stack([field[ii - 1, jj], field[ii + 1, jj], field[ii, jj - 1], field[ii, jj + 1]])
        nb = np.where(nb > 1e-30, nb, np.nan)
        self_val = field[ii, jj]
        with np.errstate(invalid='ignore'):
            med = np.nanmedian(nb, axis=0)
        ok = np.isfinite(med) & (med > 1e-30) & (self_val > 1e-30)
        return np.where(ok, self_val / np.where(ok, med, 1.0), np.nan)

    for name, g in [('clustered', clus), ('stranded', strand)]:
        ii = g.i.to_numpy().astype(int)
        jj = g.j.to_numpy().astype(int)
        for f in ('RHOG', 'TL'):
            o = outlier(fields[f], ii, jj)
            o = o[np.isfinite(o)]
            print(f'  {name:10s} {f:5s} value/median(neighbours): med {np.median(o):.3f} '
                  f'p90 {np.percentile(o,90):.3f} max {o.max():.2f}  (>2x: {100*(o>2).mean():.1f}%)')
except Exception as e:  # noqa: BLE001
    print(f'  skipped ({type(e).__name__}: {e})')

df.to_csv('relaxdiag_derived.csv', index=False)
print('\nderived table written to relaxdiag_derived.csv')
