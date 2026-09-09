#!/usr/bin/env python
"""Build the overshoot table: one row per clustered/stranded cell, with every candidate
input metric alongside the measured outcome and a spatial overshoot score.

Writes overshoot_table.csv.
"""
import numpy as np
import pandas as pd
import yt

CSV = 'impact_relax_pTg_datacollect_relaxdiag_rank0.csv'
PLT = 'amrviz/impact_relax_pTg_datacollect/plt.nga2.cell.000002'

# --- EOS constants (cvA from the T2=1 shock normalization in simulation.f90) ---
vals = {}
for line in open('input'):
    line = line.split('#')[0]
    if ':' in line:
        k, v = line.split(':', 1)
        try:
            vals[k.strip()] = float(v.strip())
        except ValueError:
            pass
cvL, cvV = vals['Liquid cv'], vals['Vapor cv']
cvA = (1.0 / (vals['GammaA'] * vals['Gas Mach number'] ** 2)) / (vals['GammaA'] - 1.0)

d = pd.read_csv(CSV, skipinitialspace=True)

# --- candidate input metrics ---
d['cvG0'] = d.Yv0 * cvV + (1 - d.Yv0) * cvA
d['C_liq'] = d.Q1_0 * cvL                        # liquid thermal capacity  (mass * cv)
d['C_gas'] = d.Q2_0 * d.cvG0                     # gas    thermal capacity
d['fG'] = d.C_gas / (d.C_liq + d.C_gas)          # gas share of total capacity
d['mass_ratio'] = d.Q2_0 / (d.Q1_0 + d.Q2_0)     # gas share of mass
d['E_ratio'] = d.Q4_0 / (d.Q3_0 + d.Q4_0)        # gas share of internal energy
d['dT0'] = d.TG0 - d.TL0                         # thermal driving gap
d['dP0'] = d.PL0 - d.PG0                         # mechanical driving gap
d['pred_dTL'] = d.fG * d.dT0                     # capacity-share prediction

# --- outcomes. A phase that disappears reads exactly 0, which is NOT an overshoot ---
d['liq_survives'] = (d.TL0 > 0) & (d.TLs > 0)
d['gas_survives'] = (d.RHOG0 > 0) & (d.RHOGs > 0)
d['dTL_solo'] = np.where(d.liq_survives, d.TLs - d.TL0, np.nan)
d['dTL_clust'] = np.where((d.TL0 > 0) & (d.TLa > 0), d.TLa - d.TL0, np.nan)
d['dRHOG_solo'] = np.where(d.gas_survives, d.RHOGs - d.RHOG0, np.nan)
d['dRHOG_clust'] = np.where((d.RHOG0 > 0) & (d.RHOGa > 0), d.RHOGa - d.RHOG0, np.nan)
d['clust_effect_TL'] = d.TLa - d.TLs             # what clustering changed vs relaxing alone
d['clust_effect_RHOG'] = d.RHOGa - d.RHOGs
d['pred_err'] = (d.dTL_solo - d.pred_dTL).abs() / d.dTL_solo.abs()

# --- spatial overshoot score from the end-of-step field (automates the ParaView eyeball) ---
yt.set_log_level(50)
ds = yt.load(PLT)
lvl = ds.max_level
cg = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge,
                      dims=ds.domain_dimensions * ds.refine_by ** lvl)
fld = {f: np.array(cg[('boxlib', f)][:, :, 0]) for f in ('TL', 'RHOG')}


def nbr_ratio(field, ii, jj):
    """cell value / median of the face neighbours in which the phase exists (0 = absent)."""
    nx, ny = field.shape
    ii = np.clip(ii, 1, nx - 2)
    jj = np.clip(jj, 1, ny - 2)
    nb = np.stack([field[ii - 1, jj], field[ii + 1, jj], field[ii, jj - 1], field[ii, jj + 1]])
    nb = np.where(nb > 1e-30, nb, np.nan)
    with np.errstate(invalid='ignore'):
        med = np.nanmedian(nb, axis=0)
    self_val = field[ii, jj]
    ok = np.isfinite(med) & (med > 1e-30) & (self_val > 1e-30)
    return np.where(ok, self_val / np.where(ok, med, 1.0), np.nan)


ii = d.i.to_numpy().astype(int)
jj = d.j.to_numpy().astype(int)
d['TL_vs_nbrs'] = nbr_ratio(fld['TL'], ii, jj)
d['RHOG_vs_nbrs'] = nbr_ratio(fld['RHOG'], ii, jj)

COLS = ['event', 'i', 'j', 'need_self',
        'VF0', 'Q1_0', 'Q2_0', 'Q3_0', 'Q4_0', 'RHOL0', 'RHOG0', 'PL0', 'PG0', 'TL0', 'TG0', 'Yv0',
        'C_liq', 'C_gas', 'fG', 'mass_ratio', 'E_ratio', 'dT0', 'dP0',
        'TLs', 'TLa', 'RHOGs', 'RHOGa',
        'dTL_solo', 'dTL_clust', 'dRHOG_solo', 'dRHOG_clust',
        'clust_effect_TL', 'clust_effect_RHOG',
        'pred_dTL', 'pred_err', 'TL_vs_nbrs', 'RHOG_vs_nbrs', 'liq_survives']
out = d[COLS].sort_values('dTL_solo', ascending=False, na_position='last')
out.to_csv('overshoot_table.csv', index=False)
print(f'overshoot_table.csv written: {len(out)} rows '
      f'({(d.event!="STRANDED").sum()} clustered, {(d.event=="STRANDED").sum()} stranded)\n')

pd.set_option('display.width', 250)
fmt = lambda v: f'{v:.4g}'  # noqa: E731
show = ['event', 'i', 'j', 'VF0', 'Q1_0', 'Q2_0', 'C_liq', 'C_gas', 'fG', 'dT0',
        'TL0', 'TLs', 'dTL_solo', 'pred_dTL', 'pred_err', 'TL_vs_nbrs']
print('=== TOP 20 BY TL OVERSHOOT (liquid surviving both sides) ===')
print(out[out.liq_survives].head(20)[show].to_string(index=False, float_format=fmt))

showg = ['event', 'i', 'j', 'VF0', 'Q1_0', 'Q2_0', 'fG', 'RHOG0', 'RHOGs', 'RHOGa',
         'dRHOG_solo', 'clust_effect_RHOG', 'RHOG_vs_nbrs']
print('\n=== TOP 20 BY SPATIAL RHOG OUTLIER AT END OF STEP ===')
print(out.nlargest(20, 'RHOG_vs_nbrs')[showg].to_string(index=False, float_format=fmt))

print('\n=== STRANDED CELLS, TOP 15 BY TL OVERSHOOT ===')
s = out[(out.event == 'STRANDED') & out.liq_survives]
print(s.head(15)[show].to_string(index=False, float_format=fmt))
