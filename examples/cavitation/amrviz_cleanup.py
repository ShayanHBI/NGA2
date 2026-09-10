#!/usr/bin/env python3
"""Thin an amrviz/<case> directory to a uniform time grid within [tmin,tmax].

Usage:
  ./amrviz_cleanup.py amrviz/cavitation_wall_pTg --tmin 0 --tmax 2.45e-5 --dt 1e-6 --execute

Dry-run by default (reports what would be deleted); pass --execute to actually delete.
"""
import argparse
import glob
import os
import re
import shutil
import sys


def read_plt_time(header_path):
    with open(header_path) as f:
        lines = f.read().splitlines()
    nvars = int(lines[1])
    time_line = 2 + nvars + 1  # version,nvars,vars[nvars],dim,time
    return float(lines[time_line])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case_dir", help="path to amrviz/<case_name> directory")
    ap.add_argument("--tmin", type=float, required=True)
    ap.add_argument("--tmax", type=float, required=True)
    ap.add_argument("--dt", type=float, required=True, help="grid spacing to keep, e.g. 1e-6")
    ap.add_argument("--tol", type=float, default=5e-2, help="tolerance for on-grid match, as a fraction of dt")
    ap.add_argument("--execute", action="store_true", help="actually delete; default is dry-run")
    args = ap.parse_args()

    case_dir = args.case_dir
    if not os.path.isdir(case_dir):
        sys.exit(f"not a directory: {case_dir}")

    plt_dirs = glob.glob(os.path.join(case_dir, "plt.nga2.*.[0-9][0-9][0-9][0-9][0-9][0-9]"))
    by_index = {}
    for d in plt_dirs:
        m = re.search(r"\.(\d{6})$", d)
        if m:
            by_index.setdefault(m.group(1), []).append(d)
    if not by_index:
        sys.exit(f"no plt.nga2.*.NNNNNN directories found under {case_dir}")

    index_time = {idx: read_plt_time(os.path.join(sorted(dirs)[0], "Header"))
                  for idx, dirs in by_index.items()}

    # Dimensional times (seconds) need a tolerance relative to dt
    tol = args.tol * args.dt
    keep, drop = set(), set()
    for idx, t in index_time.items():
        nearest = round(t / args.dt) * args.dt
        on_grid = abs(t - nearest) < tol
        if on_grid and args.tmin - tol <= t <= args.tmax + tol:
            keep.add(idx)
        else:
            drop.add(idx)

    tvals = index_time.values()
    print(f"{len(index_time)} output steps found (t={min(tvals):.4e} to {max(tvals):.4e})")
    print(f"  keep: {len(keep)}  drop: {len(drop)}")

    srf_names = [os.path.splitext(os.path.basename(p))[0]
                 for p in glob.glob(os.path.join(case_dir, "*.pvd"))]

    drop_paths = []
    for idx in drop:
        drop_paths.extend(by_index[idx])
        for srf in srf_names:
            vtp = os.path.join(case_dir, f"{srf}_{idx}.vtp")
            if os.path.exists(vtp):
                drop_paths.append(vtp)

    total_bytes = 0
    for p in drop_paths:
        if os.path.isdir(p):
            for root, _, files in os.walk(p):
                for fn in files:
                    total_bytes += os.path.getsize(os.path.join(root, fn))
        else:
            total_bytes += os.path.getsize(p)
    print(f"  will free ~{total_bytes / 1e9:.2f} GB across {len(drop_paths)} paths")

    if not args.execute:
        print("dry-run only -- pass --execute to actually delete")
        return

    for p in drop_paths:
        if os.path.isdir(p):
            shutil.rmtree(p)
        else:
            os.remove(p)

    for srf in srf_names:
        entries = []
        for idx in sorted(keep):
            vtp = f"{srf}_{idx}.vtp"
            if os.path.exists(os.path.join(case_dir, vtp)):
                entries.append((index_time[idx], vtp))
        lines = ['<?xml version="1.0"?>',
                 '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">',
                 '  <Collection>']
        for t, fname in entries:
            lines.append(f'    <DataSet timestep="{t}" file="{fname}"/>')
        lines.append('  </Collection>')
        lines.append('</VTKFile>')
        pvd_path = os.path.join(case_dir, f"{srf}.pvd")
        with open(pvd_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        print(f"  regenerated {pvd_path} ({len(entries)} entries)")

    print("done")


if __name__ == "__main__":
    main()
