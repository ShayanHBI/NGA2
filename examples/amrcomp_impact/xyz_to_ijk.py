#!/usr/bin/env python3
"""Convert a physical point to the (i,j,k) indices of its containing cell
at the finest AMR level, for an NGA2 example case.

x-domain bounds come from the runtime 'input' file (xlo defaults to 0.0,
xhi = xlo + Domain length). y/z-domain bounds are usually hardcoded in
src/simulation.f90 (not in 'input'), so this script greps for
'ylo=', 'yhi=', 'zlo=', 'zhi=' assignments there; pass --ylo/--yhi/etc.
explicitly if that lookup fails or the case does it differently.

If Base nz==1, k is always 0 (z is never refined regardless of AMR level,
per amrgrid_class.f90 forcing rrefz=1 for a single base cell in z).

Usage:
    python3 xyz_to_ijk.py X Y [Z] [--case DIR] [--xlo X0] [--ylo Y0] [--yhi Y1] [--zlo Z0] [--zhi Z1] [--refratio R]
"""
import argparse
import re
import sys
from pathlib import Path


def read_input_params(input_path):
    params = {}
    for line in input_path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line or ':' not in line:
            continue
        key, val = line.split(':', 1)
        params[key.strip()] = val.strip()
    return params


def find_hardcoded(sim_path, name):
    if not sim_path.exists():
        return None
    m = re.search(rf'\b{name}\s*=\s*([-+0-9.eEdD_WP]+)', sim_path.read_text())
    if not m:
        return None
    try:
        return float(m.group(1).replace('_WP', '').replace('d', 'e').replace('D', 'e'))
    except ValueError:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('x', type=float)
    ap.add_argument('y', type=float)
    ap.add_argument('z', type=float, nargs='?', default=0.0)
    ap.add_argument('--case', type=Path, default=Path(__file__).resolve().parent,
                     help="example case directory (default: this script's directory)")
    ap.add_argument('--input', type=Path, default=None, help="input file (default: <case>/input)")
    ap.add_argument('--xlo', type=float, default=None)
    ap.add_argument('--ylo', type=float, default=None)
    ap.add_argument('--yhi', type=float, default=None)
    ap.add_argument('--zlo', type=float, default=None)
    ap.add_argument('--zhi', type=float, default=None)
    ap.add_argument('--refratio', type=int, default=2, help="AMR refinement ratio per level (default: 2)")
    args = ap.parse_args()

    examples_dir = Path(__file__).resolve().parent.parent
    if not (args.case / 'input').exists() and (examples_dir / args.case / 'input').exists():
        args.case = examples_dir / args.case  # bare example name, e.g. --case amrcomp_impact

    input_path = args.input or (args.case / 'input')
    if not input_path.exists():
        sys.exit(f"error: input file not found: {input_path}")
    params = read_input_params(input_path)

    def need(key):
        if key not in params:
            sys.exit(f"error: '{key}' not found in {input_path}")
        return params[key]

    Lx = float(need('Domain length'))
    nx = int(need('Base nx'))
    ny = int(need('Base ny'))
    nz = int(need('Base nz'))
    maxlvl = int(need('Max level'))

    sim_path = args.case / 'src' / 'simulation.f90'
    xlo = args.xlo if args.xlo is not None else (find_hardcoded(sim_path, 'xlo') or 0.0)
    xhi = xlo + Lx
    ylo = args.ylo if args.ylo is not None else find_hardcoded(sim_path, 'ylo')
    yhi = args.yhi if args.yhi is not None else find_hardcoded(sim_path, 'yhi')
    if ylo is None or yhi is None:
        sys.exit(f"error: could not find ylo/yhi in {sim_path}; pass --ylo/--yhi explicitly")
    zlo = args.zlo if args.zlo is not None else (find_hardcoded(sim_path, 'zlo') or 0.0)
    if nz == 1:
        zhi = zlo + Lx * nz / nx  # unused for k, kept only for a consistent dz report
    else:
        zhi = args.zhi if args.zhi is not None else find_hardcoded(sim_path, 'zhi')
        if zhi is None:
            sys.exit(f"error: could not find zhi in {sim_path}; pass --zhi explicitly")

    r = args.refratio ** maxlvl
    dx = (xhi - xlo) / (nx * r)
    dy = (yhi - ylo) / (ny * r)
    dz = (zhi - zlo) / (nz * r)

    i = int((args.x - xlo) // dx)
    j = int((args.y - ylo) // dy)
    k = 0 if nz == 1 else int((args.z - zlo) // dz)

    xcc = xlo + (i + 0.5) * dx
    ycc = ylo + (j + 0.5) * dy
    zcc = 0.0 if nz == 1 else zlo + (k + 0.5) * dz

    print(f"domain: x=[{xlo},{xhi}] y=[{ylo},{yhi}]" + ("" if nz == 1 else f" z=[{zlo},{zhi}]"))
    print(f"dx={dx} dy={dy}" + ("" if nz == 1 else f" dz={dz}") + f"  (maxlvl={maxlvl}, ref ratio={args.refratio})")
    print(f"i={i} j={j} k={k}")
    print(f"cell center: ({xcc}, {ycc}, {zcc})")
    print(f"offset from target: ({args.x - xcc}, {args.y - ycc}, {args.z - zcc})")
    if abs(args.x - xcc) > 0.5 * dx or abs(args.y - ycc) > 0.5 * dy:
        print("warning: target point is more than half a cell from the computed cell center", file=sys.stderr)


if __name__ == '__main__':
    main()
