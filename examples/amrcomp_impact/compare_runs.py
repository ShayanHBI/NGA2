#!/usr/bin/env python
"""Compare end-of-step fields between the clustering-on / clustering-off / repeat runs.

The repeat run has the SAME config as the clustered run, so its difference from it is the
run-to-run noise floor -- the clustering-off difference only means something if it exceeds it.
"""
import numpy as np
import yt

yt.set_log_level(50)
RUNS = {
    'clustered': 'amrviz/impact_relax_pTg_datacollect/plt.nga2.cell.000002',
    'nocluster': 'amrviz/impact_relax_pTg_nocluster/plt.nga2.cell.000002',
    'repeat':    'amrviz/impact_relax_pTg_repeat/plt.nga2.cell.000002',
}
FIELDS = ('RHOG', 'RHOL', 'TL', 'TG', 'PL', 'PG', 'VF')

data = {}
for name, path in RUNS.items():
    ds = yt.load(path)
    lvl = ds.max_level
    dims = ds.domain_dimensions * ds.refine_by ** lvl
    cg = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge, dims=dims)
    data[name] = {f: np.array(cg[('boxlib', f)][:, :, 0]) for f in FIELDS}
    print(f'{name}: level {lvl}, grid {dims[0]}x{dims[1]}')

print('\nmax |A-B| over all cells, and how many cells differ at all')
print(f'{"field":6s} {"clustered vs nocluster":>30s} {"clustered vs repeat (noise)":>32s}')
for f in FIELDS:
    a, b, c = data['clustered'][f], data['nocluster'][f], data['repeat'][f]
    d1, d2 = np.abs(a - b), np.abs(a - c)
    print(f'{f:6s} {d1.max():12.4e} ({(d1>0).sum():6d} cells) {d2.max():12.4e} ({(d2>0).sum():6d} cells)')

print('\nglobal max of each field per run (is anything overshooting?)')
print(f'{"field":6s} {"clustered":>14s} {"nocluster":>14s} {"repeat":>14s}')
for f in FIELDS:
    print(f'{f:6s} ' + ''.join(f'{data[r][f].max():14.4e}' for r in RUNS))

# where the two configs disagree most on RHOG
d = np.abs(data['clustered']['RHOG'] - data['nocluster']['RHOG'])
if d.max() > 0:
    idx = np.dstack(np.unravel_index(np.argsort(d.ravel())[::-1][:8], d.shape))[0]
    print('\ntop-8 cells by |RHOG_clustered - RHOG_nocluster|:')
    print(f'{"i":>6s}{"j":>7s}{"clustered":>13s}{"nocluster":>13s}{"VF_clus":>10s}')
    for i, j in idx:
        print(f'{i:6d}{j:7d}{data["clustered"]["RHOG"][i,j]:13.4e}'
              f'{data["nocluster"]["RHOG"][i,j]:13.4e}{data["clustered"]["VF"][i,j]:10.4f}')
