# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is NGA2

NGA2 is a high-performance Fortran CFD research library for solving fluid PDEs (incompressible/compressible Navier-Stokes, two-phase flows, Lagrangian particles, AMR, etc.) on Cartesian meshes using MPI parallelism. It is built on top of [AMReX](https://github.com/AMReX-Codes/amrex) for distributed data structures and AMR infrastructure. The AMReX source is the authoritative reference for distributed data structures and communication patterns — check `~/Repositories/amrex` before inventing solutions.

---

## Build System

NGA2 uses a GNUmake build system. Each example in `examples/` has its own `GNUmakefile`. There is no top-level build; you build from within a specific example directory.

```bash
cd examples/amrcomp_impact
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
mpiexec -n 4 ./nga2.dp.gnu.opt.mpi.exe -i input   # use mpiexec -n, not mpirun -np
```

The `-i <file>` flag is **required**: `param_init` (`src/libraries/param.f90`) only parses a parameter file when `-i`/`--input=` is given; without it no file is read and the run aborts on the first required `param_read`. Parameters are key-value pairs (`Key : Value`) parsed by `inputfile_class` / `param`.

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
| `amrmpsolvers/` | AMR multiphase solvers: `amrvof_class` (VOF), `amrmpcomp_class` (compressible two-phase), `amrmpinc_class` (incompressible two-phase) |
| `thermo/` | EOS and relaxation classes — see [Thermodynamics Library](#thermodynamics-library-srcthermo) below |
| `solver/` | Single-grid linear solvers: HYPRE, DDADI, FFT, diagonal, BBMG |
| `two_phase/` | Two-phase flow: `tpns_class`, `vfs_class`, `lgpc_class`, `ccl_class`, `stracker_class`, `plicnet.f90` |
| `compressible/` | Compressible solvers: `mast_class`, `matm_class` |

### Object-Oriented Conventions

- All objects use `initialize` / `finalize` type-bound procedures (no constructors).
- Classes live in `*_class.f90` files; the type is named without `_class`.
- Objects that will be stored by pointer require the `target` attribute at declaration.
- AMReX callbacks use module-level singleton pointers for callback routing — intentional, not a bug.

---

## Thermodynamics Library (`src/thermo/`)

### Abstract Bases (`src/libraries/`)

- **`material`** (`material_class.f90`): fields `name`, `ns` (=1 for pure substances), `species_names(:)`. Every method takes a trailing composition `y(:)` (ignored by pure substances). Methods follow `p, rho, T` argument precedence: `get_p/T/c/cv/gruneisen_from_rho_e(rho,e,y)`, `get_T/c/e/rhoe_from_p_rho(p,rho,y)`, `get_e/rho/h/s/g/rhoe_from_p_T(p,T,y)`, `get_hk_from_p_T(p,T,y,hk)` (subroutine), `get_p_from_rho_T`, `get_cv_from_rho_T`. Prefer `*_from_rho_e` in hot loops.
- **`thermorelax`** (`thermorelax_class.f90`): no fields; deferred `apply(dt,VF,Q,Pjump,ierr)`; default no-op `get_chem_transfer(...)`. Exports error codes `RELAX_OK=0` … `RELAX_VACUUM_VAPOR=8`.

### EOS Classes

```
material
├── ideal_gas      (gamma, cv, cp, R, q, qp)        ideal_gas_class.f90
│   └── stiffened_gas  (+ pinf)                     stiffened_gas_class.f90
│       └── nasg       (+ b, brhomax=0.9 clamp)     nasg_class.f90
├── igmix          (per-species arrays gamma, cv, cp, R, q, qp)   igmix_class.f90
└── coolprop       (only with USE_COOLPROP=TRUE)    coolprop_class.f90
```

`initialize` signatures (all keyword-optional `name`):
```fortran
call g%initialize(gamma,cv,q,qp,name)                  ! ideal_gas
call s%initialize(gamma,pinf,cv,q,qp,name)             ! stiffened_gas
call w%initialize(gamma,pinf,b,cv,q,qp,name)           ! nasg
call m%initialize(gamma,cv,q,qp,species_names,name)    ! igmix: array args, ns=size(gamma)
```
Parameters are public fields; there are no getters.

### Relaxation Classes

All store EOS objects as `liq`/`gas` pointers (pass `target` objects) and share the base `apply(dt,VF,Q,Pjump,ierr)`.

| Type | Extends | `initialize` | Model |
|---|---|---|---|
| `relax_ig_sg` | `thermorelax` | `(gas,liq)` — gas first | IG/SG p, pT (`Prelax=1`, `PTrelax=2`) |
| `relax_ig_nasg` | `thermorelax` | `(gas,liq)` | IG/NASG; adds `PThybrid=3` |
| `relax_igmix_sg` | `thermorelax` | `(liq,gas,indV,indA)` | SG/NASG liquid + igmix gas; `model` = `Prelax=1`, `PTrelax=2`, `PTgrelax=3`, `PThybrid=4` |
| `relax_igmix_nasg` | `relax_igmix_sg` | same | Pelanti 2022 instantaneous relaxation |
| `relax_igmix_nasg_finite` | `relax_igmix_sg` | same | Pelanti 2022 finite-rate (Sec. 7.2.3, eqs. 51, 71–83); rates `lambda`, `theta`, `nu` (default `huge` = instantaneous) |
| `relax_igmix_max_ent` | `thermorelax` | `(liq,gas,indV,indA)` | Constrained entropy maximization, sequential M→T→C |
| `relax_numerical` | `thermorelax` | `(liq,gas)` | Mechanical-only Newton, any `material` |

`indV`/`indA` are the vapor/air species indices in the `igmix` gas. Model constants (`Prelax`, …) are exported by `relax_igmix_sg_class`.

---

## `amrmpcomp` Solver (`src/amrmpsolvers/amrmpcomp_class.f90`)

### Thermo Setup

The case sets public pointers directly before `initialize` (which dies if `liq`/`gas` are unset):
```fortran
fs%liq=>water; fs%gas=>gas
call fs%initialize(amr=amr,name=trim(case_name))
call relax_model%initialize(liq=water,gas=gas,indV=1,indA=2)
relax_model%model=PTgrelax
fs%relax=>relax_model          ! class(thermorelax), pointer
```

### Conserved Variable Layout

`nQ = 7 + (liq%ns−1) + (gas%ns−1)`. Of each phase's `ns` species, the first `ns−1` are transported; species `ns` is the carrier, reconstructed as `1 − Σ`.

| Q index | Quantity |
|---|---|
| 1 | `VF·ρL` |
| 2 | `(1−VF)·ρG` |
| 3 | `VF·ρL·eL` |
| 4 | `(1−VF)·ρG·eG` |
| 5–7 | `ρu`, `ρv`, `ρw` |
| `Yl_lo:Yl_hi` = 8 : 6+liq%ns | liquid partial densities |
| `Yg_lo:Yg_hi` = 7+liq%ns : 5+liq%ns+gas%ns | `(1−VF)·ρG·Ys` gas partial densities |

For a pure liquid and vapor+air gas (vapor=1, air=2): nQ=8, Q(8) = `(1−VF)·ρG·Yv`. Primitive fields `fs%Yl`, `fs%Yg` (amrdata, ncomp=ns−1) exist only when ns>1.

---

## Reference Case: `examples/amrcomp_impact/`

Shock–droplet impact with NASG water and vapor+air `igmix` gas, `relax_igmix_nasg` relaxation applied through `fs%apply_relax_splitremix` (the `relax_igmix_nasg_finite` lines are kept commented out in `simulation.f90`). `Relaxation type` (p, pT, pThybrid, pTg) selects `relax_model%model`; `Merge Q` switches `merge_Q`. Only `src/simulation.f90` is in its `Make.package`.

---

## Standalone Test Drivers (no AMReX/MPI)

| Directory | Build | Tests |
|---|---|---|
| `examples/eos_relax/` | `make` (targets `run`, `verify`, `residuals`, `bench`, `compare`, `finite`, …) | single-cell relaxation, max-ent verification, equilibrium residuals, finite-rate tests, `test_pelanti_ref` (literal Pelanti transcription) |
| `examples/amrcomp_venturi/eos_check/` | `make -f Makefile.test_eos run` | NASG water/vapor/air vs Liu et al. 2026 Table 4 |

These link `src/libraries/{precision,material_class,thermorelax_class}.f90` and `src/thermo/*.f90` with local stubs (`messager_stub.f90`). Many of these files are untracked; check `git status` before assuming they exist.

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
- **Comments**: new comments are 1 line, 2 at most. No doc-comment essays, no derivations or "why this replaces the old approach" in code; rationale goes in chat. Never reference a session, a debug trace or a specific cell in a comment.
- **Never delete existing comments**, even stale-looking ones.
- No semicolon packing (except a short obviously-related pair) and no unnecessary `&` line breaks.

---

## Working Rules for Claude

- **Say before you edit**: before any Edit/Write to code, state in 1–2 sentences what will change.
- **Never edit the user's `input` files** for test runs, not even temporarily. Copy to `input_test` and run with `-i input_test`.
- **Keep `input` in sync once a change is agreed**: when a code change adds, renames or removes a parameter, or a setting tried in `input_test` is agreed on, update `input` to match (new keys with their agreed values, removed keys deleted) so the user's next run uses the current version. Say which lines changed.
- **No background runs**: never launch simulations or long processes in the background (`&`, `nohup`, `run_in_background`) without asking first. Earlier approval does not carry forward. Check for orphaned `nga2` processes after a crash or timeout.
- **Follow papers literally**: when implementing a published method, use the paper's variable names and equation order, comment blocks with equation numbers, and add no machinery the paper lacks (substeps, extra closures). Mark unavoidable departures in one line. If the paper's scheme has a weakness, report it with numbers and ask.
- **Rigorous derivations**: derive from first principles, state what every derivative/capacity holds fixed, and verify key identities numerically against the solver's own EOS routines (standalone Fortran drivers) before presenting. Label results derived / verified / hypothesis.
- **Fix, don't work around**: when a feature is blocked by other solver logic, analyse why and propose a structural fix with evidence. No debug-only bypass switches unless asked.
- **Theory documents** must read top to bottom: every term defined before use, every equation derived, no case-specific sample numbers in the framework (those go in verification drivers/notes).
- **Verify against purpose before handoff**: an implementation is done when it does its job on real data (e.g. the cells it exists for are actually handled), not when it compiles and components match formulas.

---

## Testing

No separate test suite. Tests are example cases (check `monitor/` files and visualization output) and the standalone drivers above.

```bash
cd examples/amrcomp_impact
mpiexec -n 8 ./nga2.dp.gnu.opt.mpi.exe -i input_test   # full simulation from a copy of input
cd ../eos_relax && make run                             # standalone EOS/relax test
```
