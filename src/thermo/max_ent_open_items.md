# Max-ent relaxation: open items

Working list for the max-ent relaxation classes and `relaxation_ent_max.tex`, to go through one by one.

## HANDOFF 2026-11-11 REVIEW (read first; the class-by-class review of the four classes is complete)
The user reviewed the four classes one at a time (flow chart, concerns, code-vs-tex mismatches). Everything below is in the
code and the tex (PDF rebuilt), verified on the eos_relax drivers and scratch drivers; nothing is committed.
### Scope rules (user, also in memory)
- Never test at or past the critical point (inputs and end states: p_v<22.064 MPa, every T<647.096 K; a compressed liquid at
  GPa with T<T_c is in scope). Max-ent classes are used only through `apply_relax_splitremix`: test VF=0.5 reactors and
  cells that nucleate; extreme-VF inputs are out of scope.
### Changes, by class
- `relax_igmix_max_ent_inst` (EOS-agnostic):
  - `require_unclamped` at the top of every `pTg_relax` (base and finite): dies unless a NASG liquid has `brhomax>=1`. Drivers
    set `liq%brhomax=1`; the impact case already sets `water%brhomax=1e10`. Cases are changed as they are wired in.
  - `material%get_rho_max` (huge; 1/b in `nasg`) checked first in `liq_state`: the agnostic class used to accept p_l=+Inf past
    the co-volume (brhomax=1 divides by 0), which failed 20 bench states; the test is finite, so it holds under -ffast-math.
  - `max_ent_apply`: `admissible(VF,Q,Pjump,ierr)` start check before p, pT, pTg (base: `get_state`). pTg near-pure cells
    (VF<0.01 or >0.99) are relaxed only if seeded; unseeded two-phase -> `RELAX_DEGENERATE` (split-remix reactor takes them),
    unseeded pure -> `RELAX_OK` (single phase is the S maximum, eq. nuc_cond). `RELAX_NUC_FAILED` only when a seed was undone.
  - `relax_M` default `pI_end`: `NR_stage_pI_end`, damped Newton in (delta,p*) on R_l=p_l-p*, R_g=p_g+p_j-p* (tex: Mechanical
    stage with p_I=p* for a general equation of state, eqs. mech_path, mech_end, mech_newton); = NASG closed form to 5e-12,
    max 6 iterations, dS>=0. Other `pI_mode` keep `NR_stage_block(1)` with p_I held at the start (comment: no dS>=0 guarantee).
  - `iterate` renamed `NR_stage_block`, `iterate_pI_end` renamed `NR_stage_pI_end`, `closure_restore` renamed `restore_faster_eq`, `closure_pelanti` renamed `energy_split`, `chemical_transfer_exact` renamed `finite_stage_C` (user): damped Newton on the stand-alone leading ns x ns block; older blocks below keep the old name.
  - `eliminate` lost its unused `mv`, `logm` arguments (user); `bench_stages` stageA does its ln m_v mass step itself.
  - `force_tol(st,k,Pjump)` factored out of `converged` (eq. conv_tol, now in the tex with the pivot rule).
  - Comments fixed: `get_dir` (p_I=p_l only for relax_C), p_I option rationale.
- `relax_igmix_nasg_max_ent_inst`: closed `admissible` (`nasg_admissible`: masses, 0<VF<1, W0>0, T_l>0, T_g>0; before, the
  closed-form pT relaxed cells with T_g<0 or T_l<0); `nasg_relax_M` sets K^M (p_I=p*), `nasg_relax_MT` sets K^T (eq. K_T at the
  start). No formula mismatch with the tex.
- `relax_igmix_max_ent_finite`:
  - chi of the chemical path from `eliminate` with zero forces (eq. VE_C), not the residuals.
  - Root bracket symmetric: evaporation doubles m_v, condensation halves it, both with `path_admissible` and the phase-mass
    limit (tex updated).
  - Chemical switch: frozen exponential (eq. finite_solution) if |dm|<=`dmv_lin`*m_v** (0.1, error <=5%) and m_l kept; else, and
    for a dry gas, the exact integral (eq. dm_dt_exact) as a safeguard the user does not want reached; `chem_exact` flags it.
    Tex paragraph "Choice of the solution". Rejected after tests: consumption-only rule (nearly dry evaporation 97-100% wrong),
    seeding a dry gas then exponential (the seed, not nu, sets the evaporation).
  - Zero chemical rate returns before any root search. tauC is in the tex now (relaxation-time paragraph), not a departure.
  - `restore_faster_eq` default; `energy_split` kept as Pelanti's model, comment says it leaves T_l/=T_g (by design).
  - The frozen-path solution is meant for small rates (user): at nu*K*dt>=10 it can pass the equilibrium by ~1.3% of the
    transfer; accepted.
- `relax_igmix_nasg_max_ent_finite`: closed `admissible`; duplicate ierr line removed. Finite pT now agrees with the agnostic
  finite class to 1e-12 (both p_I=p*).
- Outside the classes: `apply_relax_splitremix` regular branch tests masses `pQ(1:2)` only (Q(3)<0 is legal NASG with q<0).
### Answers recorded
- No air: complete condensation (no two-phase root) needs (m_l+m_v)/V>=1/(b+R_l T/(pinf+p_j)); never reached by 400k
  subcritical VF=0.5 reactors -> past critical only, not treated (tex says so).
- End of pTg: all three F^ below `force_tol` (only exit with RELAX_OK); F^_M is floor-limited in ~18% of states (round-off of
  p_l from the cell totals), F^_T, F^_C at the set tolerances.
- The frozen-exponential chemical failures that led to the exact integral were subcritical (1e5 Pa reactor), not past-critical.
### Open
- Impact case compiled after the get_rho_max change but not after `admissible`, `pI_end`, pivots and the finite switch:
  compile-check it (back up the exe). Max-ent classes not wired into any simulation.
- Cost: finite pTg with the integral 4-5 us/cell vs ~1.3 us with the exponential; options offered (adaptive panels), not chosen.
- Old classes (not max-ent): test_pelanti_ref FAIL 7.6e-6 (`relax_igmix_nasg_finite`); `relax_igmix_nasg` pTg failures.
- Scratch drivers of this session (not in the repo) in
  /private/tmp/claude-501/-Users-shayanhbi-Repositories-NGA2/11bfe110-51f0-4b23-b1bb-9dcce07c4907/scratchpad/
  {noair,ftol,fin,pend,expo,sw,dry,nasgadm}/.

## HANDOFF 2026-11-11 morning (superseded by the review handoff above; read first; supersedes the 2026-11-10 handoff under item 10)
### Class structure (installed 2026-11-11)
| class | extends | liquid | p_relax / pT_relax | pTg_relax |
|---|---|---|---|---|
| `relax_igmix_max_ent_inst` | `relax_igmix_sg` | any EOS (material API; Z_l from Gamma, c, c_v, tex Jacobian in the transfer coordinates) | `relax_M` = `iterate(1)`, p_I held at the start state (`pI_mode`, default `pI_impedance`; `pI_end` stops with a message) / `relax_MT` = `iterate(2)` | `relax_MT` then `relax_C` |
| `relax_igmix_nasg_max_ent_inst` | `relax_igmix_max_ent_inst` | SG/NASG | closed forms: `relax_M` with p_I=p* (`pI_end` default; other modes `iterate(1)`), `relax_MT` (eq. pT_quad) | same |
| `relax_igmix_max_ent_finite` | `relax_igmix_max_ent_inst` | any EOS | rates `MRR`/`TRR`/`CRR`/`tauC`, exponential solution `finite_step`, `chemical_transfer_exact` | instantaneous: `relax_MT` then `relax_C` |
| `relax_igmix_nasg_max_ent_finite` | `relax_igmix_max_ent_finite` | SG/NASG | closed forms as the NASG inst class; `mixture_energy_closure` = Pelanti's closed form (eqs. 59/70/82) | same |
- The NASG closed forms are module procedures in `relax_igmix_nasg_max_ent_inst_class.f90` (`nasg_covolume`, `nasg_liq_state`,
  `nasg_liq_v_from_p_T`, `nasg_relax_M`, `nasg_relax_MT`), called by both NASG classes through thin type-bound wrappers.
- `liq_v_from_p_T` (binding): liquid v(p,T) for the chemical transfer path and the condensation seed (material API / NASG formula).
- EOS-agnostic `mixture_energy_closure` = tex option 3 (heat along d_T until F^_M=0, Newton, one step for SG/NASG; marked Departure).
- Nucleation stays in `max_ent_apply` (user), after the model check, pTg only; both finite classes get it through their apply.
- `full_condense` removed (user: we never work past the critical point). `iterate(3)` (direct 3x3) kept in the code, not used by pTg.
- Fortran pitfall found: `this%parent%proc()` passes the parent component, so every binding inside resolves to the parent type
  (the NASG finite class silently used the iterative pT). Call the stages on `this` instead.
### Names agreed with the user (2026-11-11)
`phasic_state` (was mestate), `st%dry_gas`, `get_dir` (get_D), `dampedNR` (damped_step), `iterate` (stage), `relax_M`/`relax_MT`/
`relax_C` (mech_relax/mt_relax/chem_stage), `xv_dry_seed`, `it_M`/`it_T`/`it_C` (nit_*), `chemical_transfer_exact` (chem_transfer),
`mixture_energy_closure` (pelanti_closure), `itmax_mv_eq`/`tol_mv_eq` (itmax_root/root_tol), and inside the chemical transfer
`chemical_rate`, `fc_path`, `path_admissible`, `vapor_entropy`, `chemical_force`, `chemical_stiffness`, `mv_of_theta`, `dt_dtheta`,
`gauss_panel`. Unchanged on purpose: parent/solver names (apply, p_relax, pT_relax, pTg_relax, model, do_nucleate, nucleated,
p_cav, Tctol, CRR_loc, tauC_loc), MRR/TRR/CRR/tauC (shared with relax_igmix_nasg_finite and case inputs), closure_pelanti,
itmax_time, time_tol, max_ent_apply, iVQ, get_state, liq_state, get_jac, eliminate, transfer, converged, KM/KT/KC, bL, tex-symbol locals.
### Convergence floor in `converged` (code only, not in the tex)
Force k is accepted at max(tol_k, 10*eps*(|z_k| + V*(|Z_l(k,1)|+|Z_g(k,1)|) + (|E_l|+|E_g|)*|Z_l(k,2)| + |E_g|*|Z_g(k,2)|)):
the round-off of the force from the split of the cell totals, V_g=1-VF and E_l from energy conservation. The old floor
10*eps*gamma*pinf missed the V/(V-b*m_l) amplification (1 bail at VF=0.5, cold water) and the trace-phase terms (nucleation seeds).
### Verification 2026-11-11 (user: test only VF=0.5 and cells that nucleate)
- VF=0.5, 7 sets x 200k subcritical (make-cmp generator with a bisection reference of F^_C on the pT manifold instead of the
  scheme-a filter; liquid 1e8-2e9 and 2e9-6e9 Pa pure NASG, near equilibrium 1e-3/1e-2, pure vapor near p_c): p, pT, pTg 0 failures
  for the inst and finite classes; NASG pTg 0.57-0.80 us, m_v to 1e-10; agnostic pTg 1.3-2.8 us, m_v to 4e-9; pT 0.016 vs 0.4-1.3 us.
- 200k cells that nucleate (cavitation and condensation seeds, 10% exactly pure), through apply: 0 failures for all four classes
  (code before 2026-11-11: 118 failures).
- make maxent (same ierr in every case), cmp (100%), ps14 PASS, bench (157777/157780: the 3 are b*rho_l>1 inputs that exist only
  on the nasg_class clamp), finite and run execute.
- Scratch drivers (not in the repo): ptgsplit.f90 (VF=0.5 sets, base vs NASG), ptgorder.f90 (pTg order study), nucdrv.F90
  (nucleating cells, -DNASG/-DFIN/-DNASGFIN) in /private/tmp/claude-501/-Users-shayanhbi-Repositories-NGA2/
  37e2b555-67cc-487a-be8b-772568144b46/scratchpad/{final,split,pTg}/ (build scripts build.sh, nuc/build_nuc.sh; /tmp is cleared on reboot).
### pTg order question (user, answered)
Direct 3x3 and pT-then-chemical reach the same state; the 3x3 failed 0.7-1% of far-from-equilibrium VF=0.5 states (ierr 1 step
leaves the admissible states, ierr 6 non-positive pivot), each then needing the pT-then-chemical fallback; the latter alone is
100% and about half the cost.
### Known limits / still open
- Liquid traces (VF~1e-9..1e-8): the trace volume is resolved only to ~4e-13 absolute (dp_l/dV_l~2e18 Pa near the co-volume, VF
  is one double); within tolerance, physically negligible.
- brhomax clamp still in nasg_class, relax_igmix_nasg(_finite), relax_igmix_sg floor_ok; the EOS-agnostic classes evaluate a NASG
  liquid through nasg_class, i.e. on the clamped EOS past b*rho=0.9 unless brhomax is set to 1 (test drivers do). User decision.
- Not built yet in amrcomp_impact (thermo Make.package registers the two new files); classes not wired into a simulation.
- test_pelanti_ref FAIL 7.6e-6 (relax_igmix_nasg_finite, untouched, not investigated); relax_igmix_nasg inherited pTg_relax
  (63.7k fails) to be replaced by the max-ent chemical step; finite chemical_transfer_exact frozen-path failure (vapor past p_c).
### User working rules reinforced this session
- Order: verify, then tex, then code; show evidence and wait for an explicit go before each code edit; build in scratch first.
- Never drop or bypass the user's code or theory: EOS-specific shortcuts go into a subclass, the general version stays.
- Name things after the physics; ask before renaming. Explanations must be concrete (file:line, numbers), not vague.
- Test sets: VF=0.5 reactors and cells that nucleate only; inputs and end states subcritical; no clamp in test inputs.

Context: the simulation uses split-remix (`apply_relax_splitremix`), whose reactor is always at VF=0.5,
so complete evaporation (red Open note at the end of the tex Nucleation section) cannot happen there
and is not on this list.

## 1. Finite-rate condensation into a pure or vapor-rich gas (tex Open iii) — DONE 2026-11-06
- Failure was: `test_max_ent` section 9, reactor "cold water / pure vapor 360 K" (VF=0.5), `tauC=dt` gave
  ierr=5; the linear F^_C with v_v=R_vT/p_g frozen condensed 3.1x the vapor present.
- Pelanti 2026 has no fix: her eq. 44 clips only volume fractions outside (0,1); here VF stays ~0.50.
- Fix (tex "Chemical step with exact vapor potential", eqs. pv_dalton, FC_path, FC_exact, KC_exact;
  `chem_transfer`): p_l, T, V_g linear in the transfer with Pi, Theta, chi_V frozen at **, and
  F^_C=g_l(p_l,T)-g_v(p_v,T) from the EOS with the Dalton partial pressure p_v=m_v*R_v*T/V_g. One
  formula for pure vapor and mixtures (the separate linear pure-vapor branch is gone); F^_C->+inf as m_v->0.
- Results: new test 5b (VF=0.5 condensation, Yv=1, 0.9, 0.5, nu*K*dt=0.1..10) vs RK4 rate law: rel.err
  <=4.3e-7 (frozen exponential gave -1.7..-3.1x the vapor present). Test 5 (vapor+air): rel.err <=1e-5,
  was 1e-4..5e-2. Test 4 (vs Pelanti's linear solution, small omega): unchanged to <=2e-9. All section-9
  reactors ierr=0.

## 2. Condensation seed can lower the entropy — DONE 2026-11-06
- Was: near T_sat (humid air, T_sat-T_g=0.1..10 K) `nucleate` removed ~3% of m_v with dS<0 (-3.4, -3.2,
  -1.4): the 5% heating cap let the latent heat carry the gas past saturation.
- Fix (user: unify, drop heating cap): one rule for both seeds, tex eqs. nuc_dA, nuc_cap:
  A=T_P*a=g_P-B_s, dA/ddm=-(Z_P3-v_s*Z_P1+s_s*Z_P2).sigma, dm<=A(0)/|dA/ddm|; nuc_cap_cav is now its
  stiff-liquid approximation. Kept: dm<=eps_n*rho_s, half-mass limits. Dropped: dT_nuc_max.
- Code: `liq_state`/`gas_state` split out of `get_state` (same arithmetic) so `nucleate` gets Z_P of a
  single phase.
- Verified (scratch nuccap.f90, solver EOS): analytic dA/ddm vs FD of exact A <=7e-8; at the cap
  a(dm)/a(0) in [-0.008,0] and dS/(a(0)dm)=0.49-0.50 (cavitation p_l=2e3..6.3e4, condensation Yv=0.5 and
  pure vapor, 0.1..10 K subcooling); new cavitation cap within 0.15% of the old one. Code nucleate on the
  former failures: dm/m_v=8e-5..8e-3, dS>0 at 0.50x first order. test_max_ent: end states unchanged.

## 3. Rate values (tex Open ii)
- lambda, vartheta, nu (or tauC) are physical input.
- Code-only departure: `tauC` sets nu=1/(K^C*tauC) with K^C at the end state of the step; Pelanti freezes
  nu at the start state (nu**). Not in the tex.
- Earlier plan: pivot to a local thermal rate (memory: local rate framework).
- Seen after item 1: in the pure-vapor reactor, K^C at the end state is ~23x K^C at the start (K^C ~ 1/m_v),
  so `tauC=dt` now gives a small nu and condenses little (p 8.3e4 -> 7.2e4 Pa; equilibrium 3.7e3 Pa).
  Which K^C defines tauC is part of this item.

## 4. Not wired into the simulation
- `amrcomp_impact` still runs `relax_igmix_nasg` through split-remix. The max-ent classes compile in the
  case but have never run there. On hold until the user says so.

## 5. Split-remix rhoG/pG overshoots
- Liquid-liquid squeezing, `src/amrmpsolvers/split_remix.md` §10. Planned next steps: equilibrate the gas
  against the interface liquid pressure; mechanical rate per unit time.

## 6. Seed energy conservation with p_j — DONE 2026-11-06
- Was: the seed conserved E_l+E_g, not E_l+E_g+p_j*V_l. Relevant: apply_relax_splitremix passes
  Pjump=sigma*curv to every cell, near-pure ones included.
- Fix: seed exchanged with dx_g=-A dx_l (seed state held, parent dx_P=-A_P*sigma*dm, A_P=A^-1 liquid,
  A for gas). Consistently, B_s uses the Laplace pressure p_P*=p_l-p_j (cavitation) / p_g+p_j
  (condensation): criteria p_vsat(p_l,T_l)>p_l-p_j and T_g<T_sat(p_g+p_j,p_v); liquid seed at (p_g+p_j,T_g);
  nuc_dA uses A_P*sigma; nuc_cap_cav numerator p_vsat-p_l+p_j; nuc_cons now includes p_j*dalpha_l.
- Verified (scratch nuccap2.f90, p_j=-2e4,0,2e4): analytic dA vs FD <=3e-8, dS/(a0*dm)=0.50 at the cap.
  test_max_ent section 8: p_l=7e4>p_vsat nucleates through p_j=+2e4 and ends at p_l-p_g=2e4; p_l=5e4 with
  p_j=-2e4 is left untouched (the old criterion seeded it with dS<0); E_l+E_g+p_j*V_l conserved to 1e-16.
  Condensation T_sat shifts only 0.003 K per 2e4 Pa (liquid Poynting), so there the energy term dominates.

## 7. Simultaneous transfers (tex Open i, theory only)
- Not needed by the sequential (rate-separated) scheme.

## Out of scope of the 2026-11-06/07 chat
Items 3, 4, 5, 7 need study and are parked. The next chat works only on item 8.

## Code housekeeping done 2026-11-07
- Rename: `relax_igmix_max_ent_class.f90` -> `relax_igmix_max_ent_inst_class.f90` (git mv), module
  `relax_igmix_max_ent_inst_class`, type `relax_igmix_max_ent_inst`. Updated Make.package, the finite class,
  eos_relax Makefile/drivers, CLAUDE.md. Not renamed: `finite_rate_relaxation.md` (historical, old step_C design).
- Comments: stale "Departure" notes (ln m_v Newton, dry-gas seed) now cite tex "Dilute vapor"; get_D cites
  "Matrix decomposition"; `tauC` marked "Departure (not in tex)".
- Finite class `chem_transfer` limits exposed as components: `itmax_bracket`=2000, `itmax_root`=200,
  `root_tol`=1e-14, `itmax_time`=60, `time_tol`=1e-13. Not converged -> RELAX_FAILED. A collapsed bracket
  (hi-lo below the tolerance) also counts as converged: round-off in the quadrature kept |dtheta|~1e-11 with
  the bracket at 2e-15. test_max_ent output identical before/after.
- Instantaneous class already has `FM_tol`, `FT_tol`, `FC_tol` (1e-10 relative, floored at round-off) and
  `itmax`=50 per stage (RELAX_FAILED when reached).
- State update (tex options): option 2 (EOS at every iterate) in the instantaneous class and the finite
  default `closure_restore`; option 3 = `closure_pelanti`; option 1 not implemented. Option 2 is the right one
  (test 7: Pelanti closure leaves T_g 333.9 K vs 300.2 K after one linear mechanical step).
- p_I options: `pI_mode` = pI_liquid, pI_gas (p_g+p_j), pI_mean, pI_impedance (default). Only the M stage
  uses it; T and C stages hold p_I=p_l (get_D with ns>=2).

## 9. relax_igmix_nasg p_relax / pT_relax: conservation and failures — DONE 2026-11-09
Context: cost comparison max_ent vs relax_igmix_nasg (item 10) showed that relax_igmix_nasg failed some states and did
not conserve. Only `p_relax` and `pT_relax` were fixed; its `pTg_relax` (inherited from relax_igmix_sg) is untouched
because the user will replace its chemical step by the max-ent one.
- Causes and fixes (`src/thermo/relax_igmix_nasg_class.f90`):
  1. Interface energy: both phases worked at one p_I, so E_l+E_g was conserved instead of E_l+E_g+p_j*V_l.
     p_relax: liquid works at p_I=phist*Peq+phi0*(Pint+cJ*p_j), gas at p_I-p_j (tex d_M); Pint+cJ*p_j is the
     impedance-weighted pressure between p_l and p_g+p_j; defaults phist=1, phi0=0 -> each phase works at its own final
     pressure. Quadratic: the old coefficients with the gas work pressure c -> c-p_j; update Q(4)+=(p_I-p_j)*dVF.
     pT_relax: from C*T+(pinf+p_j)*(V_l-b*m_l)=Es, Es=E_th+p_j*(V_l^0-b*m_l), and the volume constraint,
     omega*(P-p_j)*[C*(P+pinf)+(pinf+p_j)*A1]=Es*[A1*(P-p_j)+A2*(P+pinf)] (C=m_l*cv_l+m_g*cv_g, A1=m_l*R_l, A2=m_g*R_g,
     omega=1-b*m_l): the old quadratic with E_th->Es and pinf->pinf+p_j in the A1 terms; larger root (f(p_j)<0 for Es>0).
  2. Energy error even at p_j=0 (up to 1.6e-6 rel.): (-b+sqrt(b^2-4ad))/(2a) cancels when b~pinf and the root is ~1e2 Pa,
     and Q(3)=(VF-b*m_l)(P+gamma*pinf)/(gamma-1)+m_l*q_l is a difference of ~1e9 numbers. Fix: same root as 2d/(-b-sqrt)
     when b>=0 (both routines); pT_relax sets Q(4) from the gas EOS and Q(3) from conservation of E_l+E_g+p_j*V_l.
  3. 2106 pT states rejected: `any(Q(1:4)<=0)` treated a negative liquid energy (q_l<0, cold compressed liquid) as
     degenerate. Now: masses >0 and 0<VF<1 (both routines).
- Verified (`make cmp`, 200k random subcritical VF=0.5 reactors, RHOGmin=0 as in amrcomp_impact): p_relax and pT_relax
  converge all; masses exact, E_l+E_g+p_j*V_l to <=7e-16; p_relax |p_l-p_g-p_j|<=2.3e-9 rel.; pT end states = max_ent pT
  to 5e-12 for every p_j (was up to 5e-5); pT cost 0.12 us (was 0.11).
- Consequences: amrcomp_impact results change (split-remix calls these routines); not rebuilt or rerun yet.
  Inherited pTg_relax still fails 63.7k of the 200k states (ierr 1) and drifts the air mass by up to 6.6e-5.
- `make finite`: test_pelanti_ref prints FAIL (7.6e-6); it uses relax_igmix_nasg_finite (extends relax_igmix_sg), not
  this class; not investigated, unknown whether it failed before this session.

## 10. Cost of max_ent p and pT relaxation — DONE 2026-11-11
### Outcome 2026-11-11
- Classes: relax_igmix_max_ent_inst (EOS-agnostic, iterative; start-state p_I, default impedance), relax_igmix_nasg_max_ent_inst
  (closed-form relax_M with p_I=p*, relax_MT), relax_igmix_max_ent_finite (EOS-agnostic), relax_igmix_nasg_max_ent_finite.
  Renamed per user (stage->iterate, mech_relax/mt_relax/chem_stage->relax_M/relax_MT/relax_C, get_D->get_dir, mestate->phasic_state, ...).
- pTg = relax_MT then relax_C (direct 3x3 iterate(3) kept, not used: it failed 0.7-1% of far-from-equilibrium VF=0.5 states);
  full_condense removed (never past the critical point); nucleation stays in max_ent_apply, pTg only.
- converged(): round-off floor of each force from the split of V (V_g=1-VF) and E (E_l by conservation) between the phases.
- Verified: VF=0.5 sets (7 x 200k, subcritical, incl. 1e8-6e9 Pa pure NASG and near equilibrium) and 200k nucleating cells
  through apply: 0 failures for all four classes (old class: 118 nucleation failures); NASG pTg 0.8 us, base 1.3-2.8 us.
- Fortran pitfall: this%parent%proc() passes the parent component, so bindings inside resolve to the parent type.
### HANDOFF 2026-11-10 (superseded by the 2026-11-11 handoff at the top)
- Goal reached: closed-form p_relax (`mech_relax`, p_I=p*, default `pI_end`) and pT_relax (`mt_relax`, from `stage(2)`) for a
  NASG liquid; tex paragraphs "Choice of p_I", "Closed form of the mechanical stage", "Closed form of the pT-relaxed state".
  Cost 0.022/0.023 us vs 0.61 (Newton). Checks: make maxent, make cmp, make finite, make ps14 (Pelanti & Shyue 2014 PASS).
- brhomax clamp REMOVED from the max_ent classes (user: Claude's clamp was a huge mistake); still in nasg_class,
  relax_igmix_nasg(_finite), relax_igmix_sg floor_ok, amrcomp_impact default 0.9 -> user decision, not started.
- Dropped by the user: items 3 (rates), 4 (wiring), 5 (split-remix overshoots), 7 (simultaneous); near-pure random-VF failures.
- NEXT: pTg order. Build a clamp-free VF=0.5 subcritical compressed-liquid set (liquid up to ~1-2 GPa with pure NASG; inputs
  and end states subcritical) plus the make cmp set; compare (a) current pTg_relax (full_condense, direct stage(3), fallback
  stage(2)+chem_stage) vs (b) full_condense, mt_relax, chem_stage. Earlier (with the clamp, before mt_relax): (b) 0.85 us vs
  (a) 1.76 us with 0 failures on the make cmp set, but 29 failures on the clamped set. Keep (b) only if 100% robust.
- Scratch drivers of this session (not in repo): scratchpad item10/ (item10.f90 sets 1-5, closed.f90 is STALE: refers to
  brhomax). Fortran pitfall: names are case-insensitive (a2/A2, C/c, st/ST, C1/c1 all bit this session).
- Other leftovers: finite chem_transfer frozen-path failure (vapor past p_c + 5 extreme states); test_pelanti_ref FAIL 7.6e-6
  (relax_igmix_nasg_finite mu->inf vs Pelanti one-pass, threshold 1e-9 unsuitable; not investigated);
  relax_igmix_nasg inherited pTg_relax (63.7k fails) to be replaced by the max-ent chemical step.
### Question
relax_igmix_max_ent_inst is ~5x slower than relax_igmix_nasg for pT (and the mechanical stage alone), while max_ent
pTg is already cheaper than relax_igmix_nasg pTg. Make max_ent p and pT as cheap as possible without losing robustness
or exact conservation of E_l+E_g+p_j*V_l.
### Data (`make cmp`, M2 Max one core, case flags, 200k random VF=0.5 reactors without inputs or end states past the
critical point; liquid -1e8..1e9 Pa, 273-640 K; gas 1e2-1e8 Pa, 200-3000 K; Y_v 0..1; p_j 0 or +-1e2..1e6)
- pT: max_ent 0.65 us/cell (direct damped 2x2 Newton, ~4 iterations from far, ~2 near equilibrium), relax_igmix_nasg
  0.12 us (closed-form quadratic). Both converge all; end states agree to 5e-12. TOL=1e-5 for max_ent: 0.55 us.
- pTg: max_ent 1.77 us, relax_igmix_nasg 2.44 us (and fails 32%). Per cell, pTg = stage(3) direct, so a cheaper
  pT does not change pTg unless pTg's fallback or a pT-first path is used.
- `make bench` (old staged A vs class): VF=0.5 grid pT 0.73 us, near equilibrium 0.42 us.
### Candidate (derived, NOT yet verified)
For a NASG liquid (b*rho<brhomax) and an ideal-gas mixture of fixed composition both slow stages have closed forms:
- pT-relaxed state: the corrected relax_igmix_nasg quadratic above (unique, independent of p_I).
- Mechanical stage with p_I held at its start value (max_ent `stage(1)`): with E_l=E_l^0-p_I*dV, E_g=E_g^0+(p_I-p_j)*dV,
  p_l(V_l) and p_g(V_l) are linear-fractional, and p_l-p_g-p_j=0 times (V_l-b*m_l)*(1-V_l) is a quadratic in V_l.
Open points: which root and its admissibility; the brhomax clamp of nasg_class; whether stage(2)/stage(1) inside
pTg (chem_stage restores, finite closure_restore) should also use them; class structure (a NASG-specific path in an
SG/NASG-general class; the tex Iteration of the stages has a red Open note on this). User's order: verify, tex, code.
### Issue 1 (mechanical p_I) — DONE 2026-11-09
- With p_I held at a start value p_I0 the F^_M=0 end state has dS/dV_l=(p*-p_I0)(1/T_l-1/T_g)/=0 (tex FM_end): pressure
  equilibrium, but not the maximum of S on its path (200k set: 139k states peak before the end, up to 5% of the step,
  entropy given back median 3e-10, max 3e-3 of the produced; never S_end<S_0). eq. pI_w needs p_I between the pressures of
  the same state; p_I=p* meets it on the whole path.
- Decision (user): instantaneous mechanical stage uses p_I=p* (tex "Choice of p_I", "Closed form of the mechanical stage",
  eqs. FM_end, mech_gas, mech_liq, mech_quad incl. the branch continued beyond b*rho=brhomax; unique, always admissible).
- Code: `pI_end` (default) -> `mech_relax` from `stage(1)`; other pI_mode values keep the Newton iteration with p_I held;
  get_D maps pI_end to the impedance p_I (finite-rate step, D of the pT/pTg iteration).
- Verified: class stage(1) = relax_igmix_nasg p_relax default to 1e-15 (200k), S maximum on its path (end slope 9e-8 of
  start, round-off), 0.023 us vs 0.61 (Newton); clamped set (p_l 2-6 GPa): exactly one branch root in every state, 0.028 us.
  `make maxent`: only model 1 changes (0 M iterations, K^M diagnostic 0), end states unchanged.
### Issue 2 (pT non-uniqueness with the clamp) — DONE 2026-11-10
- False alarm: the 53k count was a bug in the scratch check (Fortran case-insensitivity, local a2 overwrote A2). With it
  fixed: 0 states with a root on both branches; pT state unique (V_l(y)+V_g(y) continuous and decreasing in the gas pressure).
- Found on the way: liq_state Z_l was wrong by 1/kappa=10 on the clamped branch (FD check), the cause of Newton failures there.
- Decision (user): the brhomax clamp of nasg_class (added by Claude in commit 1b534bf2, impact stability campaign) does not
  belong in the max-ent formulation: NASG bounds the compression itself (p_l->inf as V_l->b*m_l, closed forms stay
  admissible), and beyond the knee the clamped EOS is not thermodynamically consistent. Removed from the max_ent classes:
  liq_state, mech_relax (one branch), full_condense, nucleate, finite Kfun and pelanti_closure use plain NASG (V_l>b*m_l
  required, else RELAX_BAD_LIQUID); brhomax field gone; tex continued-branch paragraph removed. make maxent / make cmp:
  round-off changes only.
- Still in the library (user to decide): nasg_class accessors (brhomax=0.9 default), relax_igmix_nasg (get_p_eq, get_T_lvg,
  get_coeffs_lv), relax_igmix_nasg_finite, relax_igmix_sg floor_ok; cavitation case sets brhomax=1e10 (inert).
### pT closed form — DONE 2026-11-10
- tex "Closed form of the pT-relaxed state" (eqs. pT_energy, pT_volume, pT_quad: C*omega*y^2+[omega*pi(C+A_l)-E_s(A_l+A_g)]y
  -E_s*A_g*pi=0 in the gas pressure y; E_s>0 for an admissible start, so one positive root, always admissible).
- Code: `mt_relax`, called by `stage(2)` for every pI_mode (pT_relax, chem_stage start/trials/dry seed, finite pTg and
  closure_restore). make maxent: round-off changes only (model 2 now 0 iterations, K^T diagnostic 0); pT vs nasg_finite
  pT_exact 1e-13. make cmp: pT 0.023 us (was 0.61, relax_igmix_nasg 0.123), = relax_igmix_nasg to 3.6e-12; pTg 1.54 us (was
  1.77). Pelanti 2022 one-pass theta->inf (eqs. 68, 70) is first order: leaves T_l-T_g=-3.8e-4 K on her thermal test state
  (E_g differs by 1.5e-7); the closed form gives T_l-T_g=1e-13.
- Pelanti & Shyue 2014 (`make ps14`, test_pelanti_shyue_ref.f90, literal App. C.1/C.2, SG water of their Table 1, 200k
  states): pT = their eqs. C.3-C.5 to 1e-15 in alpha (E_1 within their own mixture energy error 2e-11), one admissible root
  of C.3 in every state, PASS; their mechanical step (p_I=(p_I0+p*)/2, eq. 37) differs from p_I=p* by up to 2.3e-2 in E_1.
  tex "Relation to Pelanti and Shyue" now says pT_quad reduces to their C.3.
- Remaining in item 10: pTg order (direct 3x3 first vs closed pT + chem_stage), to rerun on a clamp-free compressed set.
### Constraints from the user
- Test inputs at VF=0.5, never past the critical point (T_l<647.096 K, vapor partial pressure <22.064 MPa, end states too).
- Keep the fastest scheme only if it is 100% robust; end states must conserve mass and E_l+E_g+p_j*V_l and have forces
  below tolerance.

## 8. Staged vs direct iteration of pT/pTg (user Q6) — DONE 2026-11-09
- Kept: pT = direct damped 2x2 (`stage(2)`); pTg = complete-condensation test of a pure vapor (tex eq. full_cond),
  direct damped 3x3 (`stage(3)`), on failure `chem_stage` (Newton in ln m_v on F^_C along the pT-relaxed states,
  bounds from the sign of F''_C). tex "Iteration of the stages" and "Dilute vapor" rewritten.
- Causes fixed: m_v+m_v(e^d-1) rounded to 0 (now m_v*e^d set directly); exp overflow of the ln m_v step (halving acts
  on the Newton unknowns, m_l bound checked first); undamped C stage left the admissible states; air mass reassociated
  with m_v under -ffast-math (`ma=max(Q(2)-Q(iV),0)` in `transfer`).
- VF=0.5 grid without inputs past the critical point (157780 states; scratch study had 168070 incl. p_v>p_c): class
  converges all (pT, pTg), conservation <=6e-15, forces < tol, entropy never decreases, m_v vs bisection reference
  <=2e-10; old A failed 3062 pTg. Time class/A: pT 0.75, pTg 1.03 (both converged), near equilibrium pT 0.63, pTg 0.61.
- Left: finite `chem_transfer` fails where its frozen path leaves the admissible states before the root (only vapor
  past p_c, and 5 pure-vapor states at 1e7 Pa/150 K with liquid at 3e9 Pa). Random-VF set (minor phase 1e-4..1e-2,
  outside split-remix): ~5% pTg not relaxed (mostly trace liquid in gas), median 22 us, worst 20 ms per failed cell.
- p alone (test_max_ent section 6, VF=0.999, p_l=-2e7, p_I near p_g): no admissible mechanical end state on the fixed-p_I
  path (gas energy exhausted at T_g->0); was ierr 3, now ierr 1.

### Original notes
### Question
`pT_relax` = stage(1) then stage(2); `pTg_relax` = stage(1), stage(2), stage(3), each converged (sequential,
NOT nested: worst case 3*itmax Newton iterations for pTg, not itmax^3). The pT and pTg end states are unique
(independent of p_I and of the path), so the self-contained pT (2x2) / pTg (3x3) Newton system can be
iterated directly from the start state. Only p_relax must keep its own M iteration (its end state depends on p_I).

### Schemes
- A (current): stage(1), stage(2)[, stage(3)] converged in sequence.
- B (candidate): one Newton loop on the leading ns block from the start state, D held at the start with the
  impedance p_I for all columns (get_D(st,Pjump) without ns; p_I=p_l gives a negative/singular mechanical pivot
  for liquid in tension, e.g. p_l=-2e7, VF=0.999 -> ierr 6), dry-gas x_v seed as in stage(3), and step halving:
  halve dx~ (up to 40 times) until get_state of the trial state is OK. Same convergence test as stage.
  Without halving, B failed cold water / pure vapor at VF=0.5, 0.99 (ierr 3).
- Reference implementation of both: `examples/eos_relax/bench_stages.f90` (schemeA, schemeB, conv), build and
  run with `make bench` in examples/eos_relax (case flags -O3 -ffast-math ...). ~10 s.

### Results so far (M2 Max, one core, -O3 case flags, best of 3 alternating passes)
- 200k random states (p_l -1e7..1e8 incl. 20% tension, T_l 280-500 K, p_g 1e3-1e7, T_g 250-2000 K,
  25% dry / 25% pure vapor, VF log-distributed 1e-4 from either end, Pjump=0):
  - pT: both-converged per cell A 0.93 us, B 0.52 us (B/A 0.56); converged A 182881, B 197094, A-only 0.
  - pTg: both-converged per cell A 1.74 us, B 1.26 us (B/A 0.72); converged A 156384, B 184323;
    A-only 1453 (B ierr 6 singular pivot: 1006, ierr 5 degenerate: 440, ierr 1: 4, ierr 2: 3).
    A fails 22% (ierr 5: 34813, ierr 3: 5586, ierr 2: 3163). Whole set B slower (0.352 vs 0.307 s) because
    it solves 29k more states and ~4k hit itmax=50.
  - Hybrid B (itmax 15) then A on failure: pTg 184370 converged in 0.316 s (A: 156384 in 0.307 s).
  - End states agree to <=1.4e-10 where both converge.
- Near equilibrium (A's converged pTg end states perturbed by 1e-3 / 1e-2 relative in VF, E_l, E_g):
  pT B/A 0.56 / 0.60, pTg 0.47 / 0.52; mean iterations pTg A 6.6 / 7.7, B 2.9 / 3.7; B converges slightly more.

### Plan for the next chat (user)
1. Redo the test with VF=0.5 only (split-remix always calls relaxation at VF=0.5), covering all possible input
   sets: liquid in tension to strongly compressed, cold to hot, gas dry / humid / pure vapor, low to high p_g
   and T_g, Pjump=0 and nonzero, models pT and pTg (and p alone for reference).
2. List the failing cases of A and B, classify by ierr and by input region, find why they fail (replay single
   cells, print the Newton history) and fix the causes (not with bypasses).
3. Clean run-time comparison A vs B on the VF=0.5 set after the fixes; keep the faster scheme.
4. Then (user's order: verify, tex, code): rewrite tex "Iteration of the stages" (and the Algorithm text) for
   the kept scheme; implement in `stage`/`pT_relax`/`pTg_relax` of relax_igmix_max_ent_inst and in the finite
   class's pTg path (pT solve then the finite chemical step); re-run `make maxent` and `make bench`.
### Known failure seen earlier (staged A)
- Shocked water / hot dry air (p_l=1e8, T_l=320, p_g=5e6, T_g=2000): A fails at VF=0.99 (pT, pTg) and 0.01
  (pTg) with ierr 5 (undamped Newton step leaves the admissible states); damped B converges them in 6-9 its.
  At VF=0.5 A converges this state (10 its).
