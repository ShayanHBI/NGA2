#!/usr/bin/env python
"""Profile which input characteristics separate overshooting cells from the bulk.

For each outcome channel, compares the extreme tail against the rest using a rank-based
AUC (0.5 = no signal, 1.0 = feature always higher in the tail, 0.0 = always lower), then
grows a small decision tree to expose rules -- distinct rules mean distinct mechanisms.
"""
import numpy as np
import pandas as pd

d = pd.read_csv('relax_dataset.csv')

# outcomes, each restricted to rows where the relevant phase survives the solve
d['rel_dTL'] = ((d.TLs - d.TL0) / d.TL0.replace(0, np.nan)).where(d.liq_survives)
d['rel_dTG'] = ((d.TGs - d.TG0) / d.TG0.replace(0, np.nan)).where(d.gas_survives)
d['rel_dRHOG'] = ((d.RHOGs - d.RHOG0) / d.RHOG0.replace(0, np.nan)).where(d.gas_survives)
d['rel_dRHOL'] = ((d.RHOLs - d.RHOL0) / d.RHOL0.replace(0, np.nan)).where(d.liq_survives)
# liquid pressure is measured against the stiffness scale, not its own tiny value
d['rel_dPL'] = ((d.PLs - d.PL0) / d.PL_plus_pinf).where(d.liq_survives)
d['rel_dPG'] = ((d.PGs - d.PG0) / d.PG0.abs().replace(0, np.nan)).where(d.gas_survives)
d['liq_lost'] = ((d.TL0 > 0) & (d.TLs <= 0)).astype(float)

FEATS = ['VF0', 'VF_dist_lo', 'VF_dist_hi', 'need_self',
         'm_liq', 'm_gas', 'mass_ratio',
         'C_liq', 'C_gas', 'fG', 'C_ratio',
         'e_liq', 'e_gas', 'E_ratio', 'eL_spec', 'eG_spec',
         'PL0', 'PG0', 'dP0', 'dP_rel', 'PL_plus_pinf',
         'TL0', 'TG0', 'dT0', 'dT_rel',
         'RHOL0', 'RHOG0', 'aL2', 'aG2', 'Z_ratio', 'covol_frac', 'Yv0', 'pv0']


def auc(x, mask):
    """P(feature higher in tail than in bulk), from rank sums. NaN-safe."""
    ok = np.isfinite(x)
    a, b = x[ok & mask], x[ok & ~mask]
    if len(a) < 20 or len(b) < 20:
        return np.nan, len(a)
    r = pd.Series(np.concatenate([a, b])).rank().to_numpy()
    return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b)), len(a)


def profile(name, target, q=0.99, absval=True):
    t = d[target].abs() if absval else d[target]
    v = t.replace([np.inf, -np.inf], np.nan)
    thr = v.quantile(q)
    mask = (v >= thr).to_numpy() & np.isfinite(v).to_numpy()
    print(f'\n{"="*90}\n{name}: tail = top {100*(1-q):.1f}% (|value| >= {thr:.3e}), n={mask.sum()}')
    rows = []
    for f in FEATS:
        a, n = auc(d[f].to_numpy(dtype=float), mask)
        if np.isfinite(a):
            med_t = d.loc[mask, f].median()
            med_b = d.loc[~mask, f].median()
            rows.append((f, a, med_t, med_b))
    rows.sort(key=lambda r: -abs(r[1] - 0.5))
    print(f'{"feature":14s}{"AUC":>7s}{"median in tail":>18s}{"median in bulk":>18s}')
    for f, a, mt, mb in rows[:12]:
        print(f'{f:14s}{a:7.3f}{mt:18.4g}{mb:18.4g}')
    return mask


def tree(mask, feats, depth=3, min_leaf=200, indent='', parent_rate=None):
    """Greedy binary splits maximising separation of `mask`; prints the rules found."""
    if depth == 0:
        return
    n = len(mask)
    base = mask.mean()
    best = None
    for f in feats:
        x = d[f].to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() < 2 * min_leaf:
            continue
        qs = np.nanquantile(x[ok], [.05, .1, .25, .5, .75, .9, .95, .99])
        for thr in np.unique(qs):
            L = ok & (x <= thr)
            R = ok & (x > thr)
            if L.sum() < min_leaf or R.sum() < min_leaf:
                continue
            # weighted gini decrease
            def gini(m):
                p = mask[m].mean()
                return p * (1 - p)
            g = (L.sum() * gini(L) + R.sum() * gini(R)) / n
            if best is None or g < best[0]:
                best = (g, f, thr, L, R)
    if best is None:
        return
    _, f, thr, L, R = best
    for lab, m in ((f'{f} <= {thr:.4g}', L), (f'{f} >  {thr:.4g}', R)):
        rate = mask[m].mean()
        lift = rate / base if base > 0 else np.nan
        print(f'{indent}{lab:34s} n={m.sum():6d}  overshoot rate={100*rate:6.2f}%  lift={lift:5.1f}x')
        if lift > 1.5 and m.sum() > min_leaf:
            sub = d[m]
            # recurse on the subset by re-indexing
            globals()['_d_stack'] = sub
            tree_sub(m, mask, feats, depth - 1, indent + '    ', min_leaf)


def tree_sub(subset_mask, mask, feats, depth, indent, min_leaf):
    if depth == 0:
        return
    idx = np.where(subset_mask)[0]
    sm = mask[idx]
    base = sm.mean()
    n = len(idx)
    best = None
    for f in feats:
        x = d[f].to_numpy(dtype=float)[idx]
        ok = np.isfinite(x)
        if ok.sum() < 2 * min_leaf:
            continue
        for thr in np.unique(np.nanquantile(x[ok], [.1, .25, .5, .75, .9])):
            L = ok & (x <= thr)
            R = ok & (x > thr)
            if L.sum() < min_leaf or R.sum() < min_leaf:
                continue
            g = (L.sum() * (sm[L].mean() * (1 - sm[L].mean())) +
                 R.sum() * (sm[R].mean() * (1 - sm[R].mean()))) / n
            if best is None or g < best[0]:
                best = (g, f, thr, L, R)
    if best is None:
        return
    _, f, thr, L, R = best
    for lab, m in ((f'{f} <= {thr:.4g}', L), (f'{f} >  {thr:.4g}', R)):
        rate = sm[m].mean()
        print(f'{indent}{lab:30s} n={m.sum():6d}  rate={100*rate:6.2f}%  lift={rate/base if base>0 else np.nan:5.1f}x')


for nm, tgt in [('TL overshoot', 'rel_dTL'), ('TG overshoot', 'rel_dTG'),
                ('RHOG overshoot', 'rel_dRHOG'), ('RHOL overshoot', 'rel_dRHOL'),
                ('PL jump (vs p+pinf)', 'rel_dPL'), ('PG overshoot', 'rel_dPG')]:
    m = profile(nm, tgt)
    print('  --- rules ---')
    tree(m, FEATS, depth=2)

print(f'\n{"="*90}\nLIQUID PHASE LOSS (TL0>0 -> TLs=0)')
m = (d.liq_lost > 0).to_numpy()
print(f'n={m.sum()}')
rows = []
for f in FEATS:
    a, n = auc(d[f].to_numpy(dtype=float), m)
    if np.isfinite(a):
        rows.append((f, a, d.loc[m, f].median(), d.loc[~m, f].median()))
rows.sort(key=lambda r: -abs(r[1] - 0.5))
print(f'{"feature":14s}{"AUC":>7s}{"median in tail":>18s}{"median in bulk":>18s}')
for f, a, mt, mb in rows[:12]:
    print(f'{f:14s}{a:7.3f}{mt:18.4g}{mb:18.4g}')
