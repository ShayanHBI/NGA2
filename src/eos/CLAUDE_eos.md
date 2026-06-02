# CLAUDE.md — `src/eos/`

EOS and mixture classes for compressible multi-phase flow. Standalone Fortran modules; no build system of their own. Require only `precision` (from `src/libraries/precision_dp.f90`).

## Module Hierarchy

```
eos  (abstract interface, no data)              eos_class.f90
└── ig   (γ, cv, q, qp)                        ig_class.f90
    └── sg   (+ pinf)                          sg_class.f90
        └── nasg (+ b)                         nasg_class.f90

mix  (abstract interface, no data)             mix_class.f90
└── igmix (array of ig species)               igmix_class.f90
```

## Argument Ordering

All thermodynamic closures follow `p, rho, T` precedence:

| Form | Arguments |
|---|---|
| Conserved | `get_*(rho, e)` |
| Pressure-density | `get_*(p, rho)` |
| Pressure-temperature | `get_*(p, T)` |
| Density-temperature | `get_*(rho, T)` |

`mix` closures carry a trailing `y(:)` (mass fractions, unnormalized; normalized internally).

## Initialize Signatures

```fortran
! ig:   initialize(gamma, cv [, q, qp, pinf, b])
! sg:   initialize(gamma, cv [, q, qp, pinf, b])   ! pinf used, b ignored
! nasg: initialize(gamma, cv [, q, qp, pinf, b])   ! both pinf and b used
```

All optional args default to zero. **Important**: gfortran requires overriding procedures to have the same number of arguments as the parent — hence `ig/sg/nasg` all share the same 7-argument signature `(this, gamma, cv, q, qp, pinf, b)` with the extras optional.

## Public Fields

| Type | Public fields |
|---|---|
| `ig` | `gamma`, `cv`, `cp` (=γcv), `R` (=(γ−1)cv), `q`, `qp` |
| `sg` | inherits + `pinf` |
| `nasg` | inherits + `b` |

`cp` and `R` are computed once in `ig%initialize` and stored. Access all parameters directly as fields — there are no getter methods.

## igmix

`igmix` holds `ns` species as `ig_ptr` (non-owning pointers; caller keeps objects alive).

```fortran
type(ig), allocatable, target :: eosG(:)
type(igmix), target           :: mixG

allocate(eosG(ns))
call eosG(1)%initialize(...)   ! each species
call mixG%initialize(ns=ns)
call mixG%set_species(eosG)    ! takes full array in one call
```

**`set_species`** signature: `(this, eos_models(:))` where `eos_models` is `class(ig), target, intent(in)`.

**Species accessors** (for external code that can't reach the private `species` array directly):
- `get_species_cv(is)`, `get_species_cp(is)`, `get_species_gamma(is)`, `get_species_q(is)`

**Mixture rules**: cv, R, q are mass-fraction-weighted averages; entropy uses per-species partial pressures.

## Integration with amrmpcomp

```fortran
class(ig), allocatable, target :: eosL     ! or type(sg)/type(nasg)
type(ig),  allocatable, target :: eosG(:)
type(igmix),            target :: mixG

! Initialize all objects, then:
call fs%set_thermo(eosL, mixG)   ! before fs%initialize
```

`set_thermo` stores `fs%liq => eosL` and `fs%gas => mixG`. After `set_thermo`, `fs%nQ = 6 + mixG%ns` (i.e., 7 + (ns−1) transported species).

**Species ordering in the simulation**: species 1…ns−1 are transported (in Q); species ns is the carrier (reconstructed). For Sembian: species 1 = vapor (Q(8)), species 2 = air (carrier).
