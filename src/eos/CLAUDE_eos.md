# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a collection of Fortran 90 modules implementing equation-of-state (EOS) models for compressible multi-phase flow simulations. The modules are designed to plug into a larger simulation framework (e.g., an AMReX-based solver) via polymorphic pointers; they are not standalone executables and have no build system of their own. The `precision` module (not in this directory) must provide the `WP` kind parameter.

## Module Hierarchy

```
eos (pure interface, eos_class.f90)
└── ig  (ideal gas + caloric params, ig_class.f90)
    └── sg   (stiffened gas, sg_class.f90)
        └── nasg (Noble-Abel stiffened gas, nasg_class.f90)

mix (pure interface, mix_class.f90)
└── igmix (ideal-gas mixture, igmix_class.f90)
```

`mix_class` also exports `eos_ptr` — a wrapper holding a `class(eos), pointer` — for storing heterogeneous EOS objects in arrays.

## Key Design Conventions

**`eos` is a pure interface with no data.** Caloric parameters (`gamma`, `cv`, `q`, `qp`) and their five accessors (`get_gamma`, `get_cv`, `get_cp`, `get_q`, `get_qp`) live in `ig`, the first concrete subclass. All subclasses of `ig` inherit them without redeclaring.

**Argument ordering for all thermodynamic closures follows `p, rho, T` precedence:**
- `get_*(rho, e)` — conserved-variable form
- `get_*(p, rho)` — pressure-density form
- `get_*(p, T)` — pressure-temperature form
- `get_*(rho, T)` — density-temperature form

`mix` closures carry an additional trailing `y(:)` (gas-phase mass fractions, unnormalized; internally normalized before use). The `mix` interface is otherwise structurally identical to `eos`.

**Each subclass `initialize` takes only the parameters it owns**, then chains up to its parent:

```fortran
! ig:   initialize(gamma, cv [, q, qp])
! sg:   initialize(gamma, cv [, pinf, q, qp])   — calls this%ig%initialize, then sets pinf
! nasg: initialize(gamma, cv [, pinf, b, q, qp]) — calls this%sg%initialize, then sets b
```

`q` and `qp` default to zero when omitted at every level. There are no positional-alias initializers (`initialize_sg`, `initialize_nasg` no longer exist).

**Subclass-specific accessors:** `sg` adds `get_pinf()`; `nasg` adds `get_b()`. These are not in the base interface.

**`get_rhoe_from_p_rho`, `get_rhoe_from_p_T`, and `get_g_from_p_T` are primary (deferred) closures**, not derived ones. Every concrete class implements them directly with the simplified formula for that EOS. For example, `ig`: `ρe(p,ρ) = p/(γ-1) + ρq`; `nasg`: `ρe(p,ρ) = (1-bρ)(p+γp∞)/(γ-1) + ρq`. The `igmix` exception: `g(p,T,y)` cannot be further simplified because entropy requires a per-species loop, so it calls `h - T*s` inline.

**`igmix` species pointers are polymorphic** (`class(ig), pointer`), so `sg` or `nasg` objects can be registered as gas-phase species. `set_species(ispecies, eos_model)` takes `class(ig), target`.

**`igmix` mixture rules:** `cv`, `R`, and `q` are mass-fraction-weighted averages of the species values; entropy uses partial pressures (mole fractions × total pressure) evaluated per-species.

## Integration with `amrmpcomp`

The solver-side component stores two polymorphic pointers:
- `eosL` — `class(eos)` — the liquid pure-substance EOS
- `mixG` — `class(mix)` — the gas mixture closure

Set both with `call fs%set_thermo(liquid_eos, gas_mix)` before `initialize`. This also fixes `ns` (number of gas species) and `nQ = 7 + ns`.

**Conserved variable layout `Q(1:nQ)`:**
| Index | Quantity |
|-------|----------|
| 1 | `VF·ρL` |
| 2 | `(1−VF)·ρG` |
| 3 | `VF·ρL·eL` |
| 4 | `(1−VF)·ρG·eG` |
| 5–7 | `ρu`, `ρv`, `ρw` |
| 8 : 7+ns | `(1−VF)·ρG·Ys` (one per species) |

Any array, loop range, or flux declaration that previously hardcoded index 8 or size 8 must be replaced with `this%nQ` / `7+ns` loops. See `amrmpcomp_generic_species_snippets.f90` for the canonical patterns (flux reconstruction, `clean_Q`, `update_Q`, stats).
