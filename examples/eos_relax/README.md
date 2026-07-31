# `examples/eos_relax`

Standalone Fortran test harness for the `thermo`-branch EOS + PTg
phase-relaxation library (`src/thermo/*`, `src/libraries/material_class.f90`,
`src/libraries/thermorelax_class.f90`). No grid, no flow solver, no AMReX,
no MPI — just the EOS classes and the relaxation model exercised directly.

It is a port of `test_relax.f90` from the `lgpc` branch's
`examples/Sembian_blastwave/`, rewritten against `thermo`'s restructured EOS
API (`src/eos/*` → `src/thermo/*`, different type names, different
`initialize`/`apply` signatures). ENASG was dropped (no ENASG class exists on
`thermo`); only NASG and SG are exercised.

## What it does

For each of two liquid EOS models — NASG (`nasg_class`) and SG
(`stiffened_gas_class`) — paired with an ideal-gas mixture (`igmix_class`,
species 1 = vapor, species 2 = air) and its matching relaxation model
(`relax_igmix_nasg_class` / `relax_igmix_sg_class`):

1. **Single-cell smoke test** (optional, gated by `Run single-cell test`):
   builds a conserved state `Q` from a primitive state (p, T, VF, Yv), calls
   `rm%apply(...)` once, and prints before/after `VF, pL, pG, TL, TG, Yv`.
2. **Saturation-curve sweep** (optional, gated by `Write saturation curve`):
   sweeps temperature over a range, relaxes at each point, and writes a CSV
   (`eos_relax_NASG.csv`, `eos_relax_SG.csv`) of the resulting saturation
   curve for comparison against IAPWS-IF97.

Both phases are driven by the **same** conditions for NASG and SG (read once,
shared) — only the EOS coefficients differ per model.

## Files

| File | Purpose |
|---|---|
| `eos_relax.f90` | The program. `nasg_block` / `sg_block` set up EOS + relax model, run the test/sweep. Local `contains`: `read_real`/`read_int`/`read_logical` (input-file key scanner), `make_Q` (primitives → conserved `Q`), `get_thermo` (conserved `Q` → primitives, incl. vapor partial pressure/density/enthalpy), `write_PTg_curve` (T-sweep → CSV). |
| `input` | Run-control flags + single-cell test state + sweep state (shared by both blocks). |
| `input_NASG`, `input_SG` | Per-model EOS coefficients + cavitation/condensation thresholds. |
| `amrvof_class.f90` | Stub module providing just `VFlo`/`VFhi` (the real one lives in `src/amrmpsolvers/amrvof_class.f90` and pulls in AMReX; `relax_igmix_sg_class%apply`/`pTg_relax` only need these two parameters). |
| `messager_stub.f90` | Stub module providing `die`/`log` (the real `messager` module depends on MPI via `parallel`). |
| `Makefile` | Plain `gfortran` build, no AMReX/MPI/NGA2 build system involved. `make`, `make run`, `make clean`. |
| `NGA2_VS_IAPWS.py` | Plots `eos_relax_{NASG,SG}.csv` against IAPWS-IF97 (p, h_v, ρ_L, ρ_V vs T) → `NGA2_VS_IAPWS.pdf`. |
| `fit_eos.py` | Independent EOS coefficient fitter (SG/NASG/ENASG) against IAPWS-IF97 + Rice & Walsh Hugoniot data. No dependency on this Fortran program's output — self-contained. |

## Build & run

```bash
cd examples/eos_relax
make run                       # builds ./eos_relax, runs it
source ~/venv/combustion/bin/activate
python NGA2_VS_IAPWS.py        # -> NGA2_VS_IAPWS.pdf
python fit_eos.py --mode all   # coefficient fitter, independent of the above
```

`make clean` removes the binary, `.mod`, `.o` files (does not touch inputs
or previously-written CSVs — clean those manually if needed).

## `input` — run control + test/sweep state (shared by NASG and SG)

```
Run single-cell test   : true     ! print before/after relax diagnostics
Write saturation curve : false    ! sweep T and write eos_relax_{NASG,SG}.csv

Test pressure                        : 1e6    ! [Pa] single-cell test
Test temperature                     : 350    ! [K]
Test initial liquid volume fraction  : 0      ! VF_in
Test initial vapor mass fraction     : 0.6    ! Yv_in

Sweep pressure                       : 1e6    ! [Pa] held fixed along the sweep
Sweep initial liquid volume fraction : 0.3    ! VF0 fed to make_Q at every T
Sweep initial vapor mass fraction    : 0.6    ! Yv0 fed to make_Q at every T
Sweep minimum temperature            : 300    ! [K]
Sweep maximum temperature            : 500    ! [K]
Sweep number of points               : 200
```

(Values above are whatever is currently checked into `input` — edit freely to
explore other conditions; the two flags are independent, so you can run only
the smoke test, only the sweep, both, or neither.)

No defaults are hard-coded anywhere in `eos_relax.f90` — every key must be
present in its file or the program `error stop`s with a message naming the
missing file/key. This is deliberate: a missing key should be loud, not
silently masked by a fallback value.

## `input_NASG` / `input_SG` — per-model EOS coefficients

14 keys for NASG (12 for SG — no co-volume `b`): liquid/vapor/air
`cv, gamma, q, qp`, plus `Liquid stiffening pressure` (`pinf`) and
`Liquid co-volume` (`b`, NASG only). These are the same fitted coefficients
used in `examples/cavitation/input_{NASG,SG}`.

Also: `Cavitation pressure threshold` (→ `rm%p_cav`) and
`Condensation temperature tolerance` (→ `rm%Tctol`) — same keys/semantics as
`examples/cavitation`: `p_cav` delays liquid nucleation until the liquid
pressure drops below this threshold (a very negative value ⇒ deep tension
required before cavitating; `1e30`/huge ⇒ no delay, nucleate as soon as the
liquid is supersaturated). `Tctol` delays vapor nucleation until the gas is
subcooled by more than this below the saturation temperature (`0` ⇒ nucleate
right at `Tsat`, no subcooling required).

## Q / conserved-state layout

Matches production (`src/amrmpsolvers/amrmpcomp_class.f90`) exactly:

```
Q(1) = VF*rhoL              Q(2) = (1-VF)*rhoG
Q(3) = VF*rhoL*eL           Q(4) = (1-VF)*rhoG*eG
Q(5:7) = momentum (unused here, always 0)
Q(8) = (1-VF)*rhoG*Yv       (vapor partial density; liq%ns=1, indV=1)
```

## Relaxation call

`rm%model=PTgrelax` is set right after `initialize`, then relaxation is
invoked exactly as production code does it
(`src/amrmpsolvers/amrmpcomp_class.f90:2808`):

```fortran
call rm%apply(dt=1.0_WP, VF=VF, Q=Q, Pjump=0.0_WP)
```

`dt` is threaded through the API but not used by the mechanical/thermal/
chemical relaxation numerics themselves — `1.0_WP` is a harmless placeholder.

## CSV output schema

`eos_relax_NASG.csv` / `eos_relax_SG.csv`, header `T,pV,Yv,VF,rhoL,rhoV,hV`:

| Column | Meaning |
|---|---|
| `T` | post-relax **liquid** temperature `TL` [K] |
| `pV` | vapor partial pressure [Pa] |
| `Yv` | vapor mass fraction in the gas phase [-] |
| `VF` | liquid volume fraction [-] |
| `rhoL` | liquid density [kg/m³] |
| `rhoV` | vapor-only density (pure-vapor closure at `T`, `pV`) [kg/m³] |
| `hV` | vapor-only enthalpy (pure-vapor closure at `T`, `pV`) [J/kg] |

Rows with degenerate states (`VF` at 0/1, non-positive phase mass, or a
sanity-check failure post-relax) are skipped — see `write_PTg_curve`.

## Known non-bug: single-cell smoke test can show "no change" for `Yv_in=1`

With `VF_in=0` and `Yv_in=1` (a bit-exact pure-vapor cell, no air), `rm%apply`
can legitimately leave `VF`/`Q` untouched: `VF=0.0` exactly is outside
`relax_igmix_sg_class`'s `interfacial` window (`[VFlo,VFhi]=[1e-12,
1-1e-12]`), and `pTg_relax`'s nucleation Newton solve can fail to converge
from deep supersaturation, returning `RELAX_NUC_FAILED` with no fallback
available for a non-interfacial cell (by design — production code defers
pure cells to the solver's own PLIC pure-cell snap). With `Yv_in<1` (vapor +
air mixture, e.g. the current `input`'s `Yv_in=0.6`) nucleation converges
normally and condensation is visible in the printed diagnostics (`VF` goes
from `0` to a small positive value, `TL`/`TG` jump up to the new equilibrium
temperature). Either way this does **not** affect the saturation-curve
sweep, whose initial `VF0=0.3` is genuinely interfacial and engages
`pTg_relax` normally — confirmed by the sweep output tracking IAPWS-IF97
closely.

## Provenance

Ported from `lgpc:examples/Sembian_blastwave/test_relax.f90` +
`NGA2_VS_IAPWS.py` + `fit_eos.py`. Dropped: ENASG block (no `enasg_class` on
`thermo`), `plot_relax.py` (stale CSV filenames), `fit_eos_H.py`/
`fit_eos_backup.py` (different purpose / superseded).
