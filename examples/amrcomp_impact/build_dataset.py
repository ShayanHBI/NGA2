#!/usr/bin/env python
"""Merge every relaxdiag CSV in datacollect/ into one dataset with derived metrics.

Writes relax_dataset.csv. No conclusions drawn here -- this only assembles the columns
so the mechanisms can be looked for.
"""
import glob
import numpy as np
import pandas as pd

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
gammaL, pinfL, bL = vals['GammaL'], vals['Liquid pinf'], vals['Liquid covolume']
gammaA = vals['GammaA']
cvA = (1.0 / (gammaA * vals['Gas Mach number'] ** 2)) / (gammaA - 1.0)

files = sorted(glob.glob('datacollect/*_relaxdiag_n*_rank*.csv'))
d = pd.concat([pd.read_csv(f, skipinitialspace=True) for f in files], ignore_index=True)
print(f'{len(files)} files -> {len(d)} rows, steps {int(d.step.min())}..{int(d.step.max())}, '
      f't={d.time.min():.3f}..{d.time.max():.3f}')

# ---------- candidate metrics, grouped by mechanism family ----------
# volume
d['VF_dist_lo'] = np.maximum(0.05 - d.VF0, 0.0)      # how far into the thin-liquid band
d['VF_dist_hi'] = np.maximum(d.VF0 - 0.95, 0.0)      # how far into the thin-gas band
# mass
d['m_liq'] = d.Q1_0
d['m_gas'] = d.Q2_0
d['mass_ratio'] = d.Q2_0 / (d.Q1_0 + d.Q2_0)
# thermal capacity
d['cvG0'] = d.Yv0 * cvV + (1 - d.Yv0) * cvA
d['C_liq'] = d.Q1_0 * cvL
d['C_gas'] = d.Q2_0 * d.cvG0
d['fG'] = d.C_gas / (d.C_liq + d.C_gas)
d['C_ratio'] = d.C_gas / d.C_liq.replace(0, np.nan)
# energy
d['e_liq'] = d.Q3_0
d['e_gas'] = d.Q4_0
d['E_ratio'] = d.Q4_0 / (d.Q3_0 + d.Q4_0)
d['eL_spec'] = d.Q3_0 / d.Q1_0.replace(0, np.nan)
d['eG_spec'] = d.Q4_0 / d.Q2_0.replace(0, np.nan)
# pressure / mechanical
d['dP0'] = d.PL0 - d.PG0
d['dP_rel'] = (d.PL0 - d.PG0) / d.PG0.replace(0, np.nan)
d['PL_plus_pinf'] = d.PL0 + pinfL
# thermal driving
d['dT0'] = d.TG0 - d.TL0
d['dT_rel'] = (d.TG0 - d.TL0) / d.TL0.replace(0, np.nan)
# acoustic / stiffness
d['aL2'] = gammaL * (d.PL0 + pinfL) / d.RHOL0.replace(0, np.nan) / (1 - bL * d.RHOL0).clip(lower=1e-12)
d['aG2'] = (d.Yv0 * vals['GammaV'] + (1 - d.Yv0) * gammaA) * d.PG0 / d.RHOG0.replace(0, np.nan)
d['Z_ratio'] = np.sqrt(d.aL2 * d.RHOL0 ** 2) / np.sqrt(d.aG2 * d.RHOG0 ** 2)   # impedance ratio
d['covol_frac'] = bL * d.RHOL0                                                 # NASG packing proximity
# saturation proximity
d['pv0'] = d.Yv0 / (d.Yv0 + (1 - d.Yv0) * vals['GammaV'] / gammaA) * d.PG0     # crude vapour partial p

# ---------- outcomes ----------
d['liq_survives'] = (d.TL0 > 0) & (d.TLs > 0)
d['gas_survives'] = (d.RHOG0 > 0) & (d.RHOGs > 0)
d['dTL_solo'] = np.where(d.liq_survives, d.TLs - d.TL0, np.nan)
d['dTG_solo'] = np.where(d.gas_survives, d.TGs - d.TG0, np.nan)
d['dRHOG_solo'] = np.where(d.gas_survives, d.RHOGs - d.RHOG0, np.nan)
d['dRHOL_solo'] = np.where(d.liq_survives, d.RHOLs - d.RHOL0, np.nan)
d['dPL_solo'] = np.where(d.liq_survives, d.PLs - d.PL0, np.nan)
d['dPG_solo'] = np.where(d.gas_survives, d.PGs - d.PG0, np.nan)
d['dVF_solo'] = d.VFs - d.VF0
d['dYv_solo'] = d.Yvs - d.Yv0
d['rel_dTL'] = d.dTL_solo / d.TL0.replace(0, np.nan)
d['rel_dRHOG'] = d.dRHOG_solo / d.RHOG0.replace(0, np.nan)
# what clustering changed relative to relaxing alone
d['clust_dTL'] = d.TLa - d.TLs
d['clust_dRHOG'] = d.RHOGa - d.RHOGs
d['clust_dVF'] = d.VFa - d.VFs
# phase-change magnitude (mass transferred across the interface in the solo solve)
d['dm_phase'] = d.Q1_s - d.Q1_0

d.to_csv('relax_dataset.csv', index=False)
print(f'relax_dataset.csv written: {len(d)} rows, {len(d.columns)} columns')
print('\nrows per event type:'); print(d.event.value_counts().to_string())
print('\nrows per logged step:'); print(d.groupby('step').size().to_string())
