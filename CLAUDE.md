# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is NGA2

NGA2 is a high-performance Fortran CFD research library for solving fluid PDEs (incompressible/compressible Navier-Stokes, two-phase flows, Lagrangian particles, AMR, etc.) on Cartesian meshes using MPI parallelism. It is built on top of [AMReX](https://github.com/AMReX-Codes/amrex) for distributed data structures and AMR infrastructure. The AMReX source is the authoritative reference for distributed data structures and communication patterns — check `~/Repositories/amrex` before inventing solutions.

---

## Build System

NGA2 uses a GNUmake build system. Each example in `examples/` has its own `GNUmakefile`. There is no top-level build; you build from within a specific example directory.

```bash
cd examples/Sembian_blastwave
make -j8
make -j8 DEBUG=TRUE
make clean
make realclean   # removes all compiled objects
```

Key `GNUmakefile` flags:
| Flag | Default | Purpose |
|---|---|---|
| `USE_MPI` | TRUE | MPI parallelism |
| `USE_AMREX` | TRUE | AMReX AMR backend |
| `USE_HDF5` | TRUE | HDF5 I/O |
| `USE_HYPRE` | TRUE | HYPRE linear solvers |
| `DEBUG` | FALSE | Debug build |
| `COMP` | gnu | Compiler suite |
| `EXEBASE` | nga2 | Executable name prefix |

The build system is in `tools/GNUMake/`. New source files must be registered in `Make.package` in their directory (e.g., `f90EXE_sources += mynewfile.f90`).

---

## Running

```bash
mpiexec -n 4 ./nga2.*   # use mpiexec -n, not mpirun -np
./nga2.*                 # single-core
```

Runtime parameters are read from a plain-text `input` file in the working directory. Parameters are key-value pairs (`Key : Value`) parsed by `inputfile_class` / `param`.

---

## Code Architecture

### Driver Pattern

Every example follows the same driver pattern defined in `src/core/main.f90`:

```
parallel_init → messager_init → param_init → geometry_init → simulation_init → simulation_run → simulation_final → param_final → parallel_final
```

Each example provides two modules in its `src/` directory:
- **`geometry.f90`** — defines `geometry_init`: constructs the grid, sets domain size, BCs, periodicity.
- **`simulation.f90`** — defines `simulation_init/run/final`: instantiates solvers, EOS/relax objects, monitors, I/O, and drives the time loop.

### `src/` Library Modules

| Directory | Purpose |
|---|---|
| `core/` | Entry point (`main.f90`), `geometry.f90` stub |
| `libraries/` | Infrastructure: `precision`, `parallel`, `param`, `inputfile_class`, `monitor_class`, `timetracker_class`, `event_class`, `messager`, `timer_class`, `mathtools`, etc. |
| `grid/` | `sgrid_class`, `pgrid_class`, `coupler_class`, `iterator_class` |
| `amrbase/` | AMR grid (`amrgrid_class`), data containers (`amrdata_class`, `amrflux_class`), visualization (`amrviz_class`), I/O (`amrio_class`), solvers (`amrsolver_class`, `amrmg_class`) |
| `amrmpsolvers/` | AMR multiphase solvers: `amrvof_class` (VOF), `amrmpcomp_class` (compressible two-phase), `amrmpinc_class` (incompressible two-phase); **`relax_class`, `relax_sg_class`, `relax_nasg_class`** (PTg relaxation models) |
| `eos/` | EOS and mixture classes — see [EOS Library](#eos-library) below |
| `solver/` | Single-grid linear solvers: HYPRE, DDADI, FFT, diagonal, BBMG |
| `two_phase/` | Two-phase flow: `tpns_class`, `vfs_class`, `lgpc_class`, `ccl_class`, `stracker_class`, `plicnet.f90` |
| `compressible/` | Compressible solvers: `mast_class`, `matm_class` |

### Object-Oriented Conventions

- All objects use `initialize` / `finalize` type-bound procedures (no constructors).
- Classes live in `*_class.f90` files; the type is named without `_class`.
- Objects that will be stored by pointer require the `target` attribute at declaration.
- AMReX callbacks use module-level singleton pointers for callback routing — intentional, not a bug.

---

## EOS Library (`src/eos/`)

### Class Hierarchy

```
eos  (abstract interface — no data)           eos_class.f90
└── ig   (ideal gas: γ, cv, q, qp)            ig_class.f90
    └── sg   (stiffened gas: + pinf)          sg_class.f90
        └── nasg (Noble-Abel SG: + b)         nasg_class.f90

mix  (abstract interface — no data)           mix_class.f90
└── igmix (ideal-gas mixture of ig species)  igmix_class.f90
```

### Key Conventions

- **`eos`** has no data. Caloric parameters live in `ig` and are inherited by all subclasses. Access all EOS parameters as **public fields directly** — there are no getter methods.
- **`ig` public fields**: `gamma`, `cv`, `cp` (=γcv, stored at init), `R` (=(γ−1)cv, stored at init), `q`, `qp`. `sg` adds `pinf`; `nasg` adds `b`.
- **Argument ordering** follows `p, rho, T` precedence: `get_*(rho,e)`, `get_*(p,rho)`, `get_*(p,T)`, `get_*(rho,T)`.
- **`mix` closures** carry an additional trailing `y(:)` (gas-phase mass fractions, unnormalized; normalized internally). Interface otherwise mirrors `eos`.
- **Initialize chain**: each subclass initializes only its own new field, then delegates up:
  ```fortran
  ! ig:   initialize(gamma, cv [, q, qp])
  ! sg:   initialize(gamma, cv [, pinf, q, qp])    ! q,qp must come before pinf
  ! nasg: initialize(gamma, cv [, pinf, b, q, qp]) ! same ordering rule
  ```
  Arguments `pinf`, `b`, `q`, `qp` are all optional and default to zero.
- **`igmix%set_species(eosG(:))`** takes the full array of species in one call (`class(ig), target, intent(in) :: eosG(:)`). Species pointers are non-owning; the caller keeps the objects alive.
- **`igmix` mixture rules**: `cv`, `R`, and `q` are mass-fraction-weighted averages; entropy uses per-species partial pressures.
- **`igmix` species accessors**: `get_species_cv(is)`, `get_species_cp(is)`, `get_species_gamma(is)`, `get_species_q(is)` expose per-species scalars to external code.

---

## Relaxation Library (`src/amrmpsolvers/relax*.f90`)

### Class Hierarchy

```
relax  (abstract base — liq/gas pointers, saturation coefficients, tolerances)
├── relax_nasg  (Pelanti 2022 3-step PTg relaxation for NASG liquid)
└── relax_sg    (quadratic-solver 3-step PTg relaxation for SG liquid)
```

### Usage Pattern

```fortran
type(relax_nasg), target :: rm   ! or relax_sg

call rm%initialize(liq=eosL, gas=mixG, AS=AS, BS=BS, CS=CS, DS=DS, ES=ES)
! Then assign to the solver:
fs%relax_model => rm
```

The `initialize` signature is `(liq, gas, AS, BS, CS, DS, ES)` where AS–ES are the Clausius-Clapeyron saturation curve coefficients:
```
AS = (CpL - CpV + qpV - qpL) / (CpV - CvV)
BS = (qL - qV)                / (CpV - CvV)
CS = (CpV - CpL)              / (CpV - CvV)
DS = (CpL - CvL)              / (CpV - CvV)
ES = bL                       / (CpV - CvV)   ! ES=0 for SG
```

`initialize` extracts `PinfL`, `bL`, and gas species scalars from `liq` and `gas` at setup time. The full 3-step relaxation (mechanical → thermal → chemical) is invoked via `rm%apply(VF, Q, Pjump)`.

### Differences Between relax_nasg and relax_sg

| Feature | NASG | SG |
|---|---|---|
| Step 1 (mechanical) | Pelanti 2022 ODE (impedance + acoustic) | Direct quadratic in pressure |
| Step 2 (thermal) | ODE system with `zetaL = ρL(1−bρL)/(p+P∞)` | Quadratic with `etaL=qL`, `etaG=qG`, PinfG=0 |
| PTsat curve | `AS + (BS+ES·p)/T + CS·ln T + DS·ln(p+P∞) − ln pv` | `AS + BS/T + CS·ln T + DS·ln(p+P∞) − ln pv` |
| `get_T_lvg` numerator | `(1−Yv) − bL·(ρ0(1−Yv) − ρA0)` | `(1−Yv)` |
| `get_coeffs_lv` | Uses `(1−bL·Σρ)` co-volume factor | Uses `etaL=qL, etaG=qG, PinfG=0` |

### Standalone Compilation (no AMReX/MPI needed)

The relax classes depend only on `precision`, `eos`, and `mix` — they can be compiled standalone for testing:

```bash
gfortran -O2 \
  ../../src/libraries/precision_dp.f90 \
  ../../src/eos/eos_class.f90 ../../src/eos/ig_class.f90 \
  ../../src/eos/sg_class.f90  ../../src/eos/nasg_class.f90 \
  ../../src/eos/mix_class.f90 ../../src/eos/igmix_class.f90 \
  ../../src/amrmpsolvers/relax_class.f90 \
  ../../src/amrmpsolvers/relax_sg_class.f90 \
  ../../src/amrmpsolvers/relax_nasg_class.f90 \
  mytest.f90 -o mytest
```

---

## `amrmpcomp` — Updated API

### EOS and Relaxation Setup

The old function-pointer API (`fs%getPL`, `fs%getPG`, etc.) has been replaced:

```fortran
! Before initialize:
call fs%set_thermo(eosL, mixG)   ! polymorphic: class(eos), class(mix)
fs%relax_model => rm             ! polymorphic: class(relax)
call fs%initialize(amr=amr, name=trim(case_name))
```

`set_thermo` stores `fs%liq => eosL` and `fs%gas => mixG`. `fs%relax_model` is a `class(relax), pointer`.

### Conserved Variable Layout

`nQ = 7 + (ns−1)` — only `ns−1` species mass fractions are transported (one carrier species is reconstructed).

**Species ordering convention**: species 1…ns−1 are **transported** (stored in Q); species ns is the **carrier** (reconstructed as `1 − Σ Ys`). For Sembian (ns=2): species 1 = vapor (Q(8)), species 2 = air (carrier, not stored).

| Q index | Quantity |
|---|---|
| 1 | `VF·ρL` |
| 2 | `(1−VF)·ρG` |
| 3 | `VF·ρL·eL` |
| 4 | `(1−VF)·ρG·eG` |
| 5–7 | `ρu`, `ρv`, `ρw` |
| 8 : 6+ns | `(1−VF)·ρG·Ys` for transported species s=1…ns−1 |

For Sembian (ns=2): nQ=8, Q(8) = `(1−VF)·ρG·Yv`. In `get_primitive`, the carrier mass fraction is reconstructed as `y(ns) = 1 − Σ y(1:ns−1)` before EOS calls.

The `Yg` field (`amrmpcomp%Yg`) stores the `ns−1` transported species mass fractions as an amrdata object with ncomp=ns−1.

---

## Sembian Blastwave Example (`examples/Sembian_blastwave/`)

### Overview

Cylindrical shock wave impacting a water cylinder in air. AMR compressible two-phase simulation with PTg phase-change relaxation. Runtime selection between SG and NASG liquid EOS.

### Input Files

| File | EOS | Notes |
|---|---|---|
| `input_NASG` | NASG | Uses `Liquid EOS type : NASG`. Reference: Le Metayer & Saurel (2016), Table V |
| `input_SG`   | SG   | Uses `Liquid EOS type : SG`.   Reference: Flatten & Lund (2011) |

Copy or symlink the desired file to `input` before running:
```bash
ln -sf input_NASG input
mpiexec -n 8 ./nga2.*
```

### EOS Parameters

**NASG** (Le Metayer & Saurel 2016, Table V, range 300–500 K):

| Parameter | Liquid | Vapor |
|---|---|---|
| γ | 1.19 | 1.47 |
| cv [J/(kg·K)] | 3610 | 955 |
| P∞ [Pa] | 7.028×10⁸ | 0 |
| b [m³/kg] | 6.61×10⁻⁴ | 0 |
| q [J/kg] | −1 177 788 | 2 077 616 |
| q′ [J/(kg·K)] | 0 | 14 317 |

**SG** (Flatten & Lund 2011, range 300–500 K):

| Parameter | Liquid | Vapor |
|---|---|---|
| γ | 2.35 | 1.43 |
| cv [J/(kg·K)] | 1816 | 1040 |
| P∞ [Pa] | 10⁹ | 0 |
| q [J/kg] | −1 167 000 | 2 030 000 |
| q′ [J/(kg·K)] | 0 | −23 400 |

Air (both cases): γ=1.40, cv=718 J/(kg·K), q=0, q′=0.

### Gas Species Ordering

In `igmix`, species ordering is:
- **Species 1**: vapor (transported, Q(8) = `(1−VF)·ρG·Yv`)
- **Species 2**: air (carrier, reconstructed as `Y_air = 1 − Yv`)

```fortran
call eosG(1)%initialize(gamma=GammaV, cv=CvV, q=qV, qp=qpV)   ! vapor
call eosG(2)%initialize(gamma=GammaA, cv=CvA, q=qA, qp=qpA)   ! air (carrier)
call mixG%initialize(ns=2); call mixG%set_species(eosG)
```

### Case Name and Output

The case name is constructed at runtime: `Sembian_blastwave_NASG` or `Sembian_blastwave_SG`. This propagates to AMR grid name, solver name, visualization directory, and restart directory.

---

## Standalone Unit Test: `test_relax.f90`

Tests the EOS and relaxation classes without AMReX or MPI.

### Build and Run

```bash
cd examples/Sembian_blastwave

# Using the dedicated Makefile (recommended):
make -f Makefile.test_relax        # compile
make -f Makefile.test_relax run    # compile + run
make -f Makefile.test_relax clean  # remove binaries and .mod files
```

`test_relax` must be run from the `Sembian_blastwave/` directory because it reads `input_NASG` and `input_SG` directly using its own lightweight key-value file reader (`read_real` subroutine) — no NGA2 `param` library needed.

### What It Does

1. Reads NASG EOS parameters from `input_NASG` (with hard-coded defaults as fallback).
2. Reads SG EOS parameters from `input_SG`.
3. For each configuration: initializes EOS objects (`nasg`/`sg`, `igmix`) and a relaxation model (`relax_nasg`/`relax_sg`).
4. Runs a single-cell PTg relaxation test (T=490 K, p=10⁵ Pa, VF=0.5, Yv=1).
5. Writes saturation-curve CSVs by sweeping T from 300 to 500 K at fixed pressure.

### Output CSVs

| File | EOS | Contents |
|---|---|---|
| `test_relax_NASG.csv` | NASG | T, p, Yv, VF, ρL, ρG, hG |
| `test_relax_SG.csv` | SG | T, p, Yv, VF, ρL, ρG, hG |

CSV column units: T [K], p [Pa], Yv [–], VF [–], ρL [kg/m³], ρG [kg/m³], hG [J/kg].

### Key Subroutines

| Subroutine | Signature | Purpose |
|---|---|---|
| `read_real` | `(file, key, val, default)` | Parse `Key : Value` from NGA2 input format |
| `make_Q` | `(liq, gas, p, T, VF_in, Yv_in → VF, Q)` | Build conserved Q from (p,T,VF,Yv) |
| `get_thermo` | `(liq, gas, VF, Q → PL,PG,ρL,ρG,TL,TG,Yv,hG)` | Extract primitives from Q |
| `write_PTg_curve` | `(file, liq, gas, rm, p0, VF0, Yv0, Tmin, Tmax, nT)` | Sweep T and write CSV |

All subroutines are generic via `class(eos)`, `class(mix)`, `class(relax)` arguments — pass any EOS/relax combination.

---

## Post-Processing (`examples/Sembian_blastwave/`)

### `NGA2_VS_IAPWS.py`

Compares `test_relax` output against IAPWS-IF97 saturation data. Requires `pip install iapws`.

```bash
./test_relax               # generate CSVs first
python NGA2_VS_IAPWS.py    # produces NGA2_VS_IAPWS.pdf
```

**Configuration** (top of script):

```python
CURVES = {
    "NASG": "test_relax_NASG.csv",
    "SG":   "test_relax_SG.csv",
}
```

To add a new EOS variant, add one line to `CURVES`. The script produces a 2×2 figure:
- p vs T (log scale) — with IAPWS-IF97 saturation pressure
- hG vs T — with IAPWS-IF97 saturated vapor enthalpy
- ρL vs T — with IAPWS-IF97 saturated liquid density
- ρG vs T — with IAPWS-IF97 saturated vapor density

Uses LaTeX fonts (`text.usetex=True`). Uses IAPWS-IF97 via `IAPWS97(T=T, x=0/1)`.

---

## Coding Standards

- **Precision**: use `real(WP)` (from `use precision, only: WP`), literals as `1.0_WP`. Never `real(8)`, `1.0d0`, or `dble()`.
- **Indentation**: 3 spaces (not 2, not 4, not tabs).
- **Logical operators**: prefer `.eq.` / `.ne.` over `==` / `/=`.
- **`use` statements**: place close to where they are used (local `use` in subroutines) unless a symbol is used throughout the whole module.
- **Block statements**: use named `block`/`end block` to encapsulate sequential steps. Do not wrap an entire routine in a single block.
- **Explicit keyword arguments**: use for calls with many parameters.
- **No spaces around arithmetic operators** or after commas in array indices and call arguments (except in type/use declarations).
- **No documentation files in the repo** unless explicitly requested.

---

## Testing

No separate test suite. Tests are example cases and `test_relax`. Run an example and check monitor output files and visualization output.

```bash
cd examples/Sembian_blastwave
mpiexec -n 8 ./nga2.*          # full simulation
make -f Makefile.test_relax run  # standalone EOS/relax unit test
```
