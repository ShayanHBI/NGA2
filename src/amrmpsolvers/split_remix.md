# Split-remix relaxation: status and open items

Working notes for `amrmpcomp%apply_relax_splitremix` (sessions of 2026-10-22 and 2026-10-23). Section 0 is the session history and its decisions. Sections 1–6 are the technical record. Sections 7–9 cover what is untested, what has not been tried, and how to resume. Section 10 is the rhoG/pG overshoot investigation of the second 2026-10-23 session.

**Current status:**
- Split-remix runs both test cases to completion: amrcomp_impact (past the earlier failure at step 11787) and the cavitation example.
- The rhoG/pG overshoots in amrcomp_impact were traced (§10). **Split-remix cannot remove them:** repeated split-remix calls converge to the whole-cell equilibrium, and in the overshooting cells that equilibrium is the crushed gas. No remix variant changes this (§10.5).
- **Next task (§9):** explore two ideas that change the end state rather than the speed: (1) equilibrate the gas against the liquid pressure at the interface (Jr-like), (2) a mechanical relaxation rate per unit time.

## 0. Session history (2026-10-22 to 2026-10-23)

### Starting point
- **Setting.** Sharp-interface compressible two-phase solver (amrmpcomp) with interface relaxation (p, pT, pTg) borrowed from diffuse-interface methods. Fractional steps: transport, then relaxation of each interfacial cell as an isolated reactor at fixed volume, energy and mass.
- **Problem.** The reactor is "not nice" when the minor phase is small. Example: too little superheated liquid; equilibrium would need to evaporate more liquid than the cell holds, so the remainder ends too hot. The result is local over- and undershoots.
- **Tried before this session:**
  - clustering: expensive, and some needy cells never find a cluster;
  - cell-local rules (the `relax_igmix_nasg_finite` safety rules): did not fully fix it and are considered inaccurate.
- **User's proposal.** Split the dominant phase into parts a and b so that b plus the whole minor phase forms a VF=0.5 sub-cell. Relax the sub-cell, then remix b with a. The open question was the mixing rule.

### How it proceeded
1. **Assessment.** The idea amounts to an interface-centred reactor that is symmetric in volume. That bounds the phase ratio, so the minor phase cannot run out or be dragged.
   - The mixing rule is forced to be additive (section 1); the real choice is what part a does during relaxation.
   - Argument offered: pressure equalises acoustically across a cell within one step, while heat and species do not. So M is whole-cell, while T and C are local.
2. **Offline study** (section 2). A Fortran batch driver runs the solver's own `apply` on cells; Python does the split/remix.
   - The first runs used `relax_igmix_nasg_finite`, which turned out not to reach equilibrium at `huge` rates. The study switched to `relax_igmix_nasg`.
   - A Python pTg solver was verified against it, then extended with the piston variant.
   - Result: rigid throttles phase change; rigid+M creates bubble overshoots; piston is good for vapor + air but has two weak spots.
3. **User clarifications:**
   - split-remix must not apply to pure cells, which keep nucleation;
   - asked how a and b were mixed (answer: additive, section 1).
4. **User's implementation spec**, with the decisions:
   - New routine `apply_relax_splitremix`: a copy of `apply_relax` without clustering or dissolve; reference commit fd296a3.
   - The user swaps call sites by hand. My proposal of a dispatch flag inside `apply_relax` was rejected.
   - `VFrb` is a local parameter.
   - **Per-cell logic:**
     - pure cells: regular relaxation;
     - |VF−0.5| < VFrb: regular relaxation;
     - near-pure interfacial: try nucleation first; if it fails, split-remix;
     - no pure-phase stability test;
     - a pure cell whose nucleation fails gets no relaxation.
   - Scope: single-species liquid and vapor + air gas only. Kinetic energy noted as a possible mixing concern; it turns out not to matter (section 1).
   - Remix options now: rigid and rigid+M. Piston needs more discussion first.
   - Guards in the reactor. The user first said all safety rules off, then refined it: **only the dissolution switches off; keep the rest for p_relax and pT_relax**.
   - Ignore finite rates for now; run instantaneous `relax_igmix_nasg`, since Pelanti's finite-rate class is inaccurate in instantaneous mode.
   - `die` only on actual failures; a failed nucleation that later relaxes fine is not a failure.
5. **Process lesson.** I started editing before being told to and was stopped. Rule since then: design discussion or "do you have everything" gets questions or a plan only; edit only after an explicit implement or approve. A plan was written, approved and implemented (section 3).
6. **Verification:**
   - defaults bit-identical on 670 cells;
   - the nucleated flag test passes;
   - cav_tube compiles.
7. **amrcomp_impact switch.** User confirmed: instantaneous `relax_igmix_nasg` with split-remix. `input_test` was created and is now missing from the directory; recreate it before the next run.
8. **Failure F1** (section 4) was traced to the `Y_small` guard. Options:
   - (1) a switch to turn the guard off in the reactor;
   - (2) a distinct "chemistry declined" return code.

   **Chose 1** (`chem_guard`).
9. **Guard review**, requested by the user. It produced open items 5.1–5.3.
10. **Failure F2** (section 4) was traced to the pTg dry-gas starting guess plus negative vapor mass from transport. It produced open items 5.4–5.5, none implemented. These notes were written at this point.
11. **Cavitation example switched to split-remix** (`examples/cavitation`; not cav_tube). Only the two calls changed; the case already used `relax_igmix_nasg`.
12. **Failure F3** (cavitation, section 4) was traced to a convergence-test bug in `solve_lv` that appears only under `-ffast-math`. **Fixed**; see 3 and 4.
    - The first version also changed the convergence test in `solve_lv` and `solve_lvg`. That made 235 of 670 regression cells stop at the 1e-5 tolerance instead of about 1e-11, so it was narrowed to the failure case only.
    - The old-path upper-root problem found along the way is recorded as 5.6.
13. **Failure F4** (cavitation) was traced with a temporary diagnostic to the case's `Tmin_liq=300` relaxation floor acting on condensation seeds (section 4).
    - Option (a) was chosen: lower only `relax_model%Tmin_liq` to 0.1 K.
    - My first attempt lowered `fs%Tmin_liq`, which `clean_Q`'s rescue also uses, and broke initialisation; it was corrected.
    - With the fix the user's run reached t = 2.436e-5 of 2.45e-5.
14. **Failure F5** (cavitation, near the end of the run): the pure-vapor `solve_lv` path had no all-liquid/all-vapor boundary handling. Fixed as fix 7 (sections 3 and 4).
15. **Cavitation input study.** Gentler inputs were tried and rejected (section 5.8, context). The user kept the current input. Solver weaknesses 5.7 and 5.8 were recorded.
16. **amrcomp_impact F2 reproduced** at step 11787 in the user's rerun, same cell (43, 1191). Fixes 5.1–5.4 were implemented. The user restarts the case themselves (from the t = 11.7 checkpoint).
17. **amrcomp_impact passes F2.** The user's restart from t = 11.7 ran past step 11787 without failing. The user reports **overshoots** in this run that still need tracking; that is the next task (§9).

Second session of 2026-10-23 (overshoot tracking, details in §10):

18. **Cells (49–50, 997), t = 11.0–11.1.** pG/TG overshoot forms at 11.01–11.02 while the cells are pure or nearly pure gas: transport and the pressure step's p dV work, plus `merge_Q` heating the gas. Relaxation is not the origin there (§10.2). A `Merge Q` input key was added; with it off the overshoot remains.
19. **Box (44–52, 995–998), t = 10.80–10.90**, logged stage by stage with new debug prints (§3, §10.1). The ρG overshoot (up to 298) sits in cells squeezed by a liquid–liquid collision; relaxation moves the liquid's compression into the cell's gas (§10.3).
20. **User question: why did split-remix not fix it?** Answer, verified: repeated split-remix converges to the whole-cell equilibrium; the additive remix merges part a back into the liquid, so the next part b carries its compression (§10.5).
21. **Gas-strain guard** (skip relaxation of a gas-minor cell when one call would compress its gas by more than s_max = 0.02): implemented on request, tested, refined with a divU < 0 condition, then **rejected by the user and removed** (§10.6).
22. **Cell-local rules of branch `thermo_gpu_friendly`** reviewed (§10.7): only the sliver dissolve would help, and only for slivers.
23. **Remix variants** (rigid+M, piston, others) judged unable to help, by a fixed-point argument (§10.5). The user will explore the two ideas in §9 in a new chat.

### Decisions in force
| Topic | Options considered | Chosen for now |
|---|---|---|
| Mixing a+b | additive (forced) | additive |
| Part a during relaxation | rigid, rigid+M, isentropic piston | rigid (rigid+M available via `remix_mode`); piston deferred |
| Dispatch | flag inside `apply_relax`, or manual call swap | manual call swap |
| Near-pure interfacial cells | (a) always split; (b) regular call for VF<0.05 / >0.95; (c) nucleation first, then split | (c), with pure-phase stability off |
| Reactor guards | all off, or selective | dissolution off and `chem_guard` off; floors, `VFratmax`, `RHOGmin` kept |
| F1 fix | guard switch, or new return code | guard switch (`chem_guard`) |
| Relaxation class | finite (Pelanti), or instantaneous | instantaneous `relax_igmix_nasg` |
| Failure policy | die on any non-OK, or classify | die on any non-OK from the reactor or M step, except `RELAX_VACUUM_VAPOR` (accepted, 5.2). Cells with any Q(1:4) ≤ 0 never reach the reactor (5.1). |
| Cavitation example input | keep, gentler pulse, seeded bubble | keep the original strong pulse; condensation nucleation stays on (fog "dust" accepted as a modelling artefact) |
| Liquid relaxation floor (cavitation) | keep 300 K, lower it | `relax_model%Tmin_liq=0.1`; `fs%Tmin_liq=300` kept for the rescue |
| rhoG overshoot in squeezed cells | cell-local guard on relaxation, another remix, change the equilibrium target, rate per unit time | guard tried and removed; remix variants ruled out (§10.5); the last two to be explored |
| `merge_Q` in amrcomp_impact | on, off | on (`Merge Q : .true.`); its gas heating is a known side effect (§10.2) |

## 1. Idea

Relaxing a whole interfacial cell as one fixed-volume reactor lets the minor phase over- or undershoot. Examples: a small liquid evaporates completely, or the minor phase is pulled to the other phase's temperature.

Split-remix instead does three steps:
- **Split.** Form a sub-cell reactor at VF=0.5 from the whole minor phase plus an equal volume (part b) of the dominant phase.
- **Relax.** Relax that reactor with the usual `apply`.
- **Remix.** Add back the untouched rest of the dominant phase (part a).

Notation:
- Vm = min(VF, 1−VF), Vd = 1−Vm, reactor volume Vr = 2·Vm.
- Dominant fraction sent to the reactor: f = Vm/Vd.

The mixing rule is not a free choice. Conservation, one state per phase, and keeping the minor phase as the reactor left it together force **additive** mixing of parts a and b: mass, species, internal energy and volume all add. This is isochoric and adiabatic, and it never lowers entropy. Kinetic energy is unchanged, because momentum Q(5:7) never enters the reactor and total mass is conserved.

The real modelling choice is what part a does during the reactor relaxation:
- **Rigid wall** (implemented):
  - **remix_rigid:** additive remix only.
  - **remix_rigidM:** the same, then a whole-cell mechanical (`Prelax`) step.
- **Isentropic piston** (not implemented, needs discussion): part a exchanges only p dV work with the reactor. It needs a new constrained-equilibrium solver.

## 2. Offline study (Python + solver, cav_tube EOS, before implementation)

Tools in `examples/eos_relax/` (untracked):
- `relax_batch.f90`: runs `apply` on a file of cells. Argument 3 selects `inst` (`relax_igmix_nasg`) or `finite` (`relax_igmix_nasg_finite`).
- `split_reactor.py`: split/remix bookkeeping, a Python pTg solver and the piston variant, VF sweeps, and a random ensemble.

The Python pTg solver matches `relax_igmix_nasg` (p = 102295 Pa in S1 at VF=0.5; all three equilibrium conditions to 1e-15).

Build note: compile from a clean directory, not from `eos_relax/`, whose stale `.mod` files shadow the current ones. `relax_igmix_nasg_finite` needs a one-line stub module `parallel` with `integer :: rank=0`. In zsh, use a bash script for the compile line, since zsh does not split `$LIB`.

Findings, per relaxation call:
- **Full cell** fails at the VF extremes:
  - A small liquid (VF=1e-6) in dry air evaporates completely.
  - New temperature extremes of about 3 K appear: T = 350.2 below both initial temperatures, 353 and 403.
- **Split, rigid:**
  - The minor phase becomes VF-independent and bounded.
  - Phase change is throttled by confinement. Liquid 20 K superheated: 0.03% evaporates against 4.1% for the full cell, because the reactor pressure rises from 0.48 to 1.02 bar.
  - After remix, the phases are left at different pressures, up to a factor 2.4. A whole-cell `Prelax` (rigid+M) then drives tiny bubbles to 312 K and 473 K.
- **Split, piston:**
  - Pressure stays consistent and evaporation matches the full cell (3.8%). Vapor + air (the impact case) is clean.
  - Gas-dominant cells with strong evaporation reach about 408 K, because part a is compressed adiabatically.
  - Condensation onto a small liquid in *pure* vapor has no root. For pure vapor, part a should also be allowed to condense.
- **Separate issue:** `relax_igmix_nasg_finite` with all rates `huge` does not reach pTg equilibrium, although it returns ierr=0. TL−TG goes down to −7.5 K, and pvsat/pv reaches 2.3.

## 3. What is implemented

**`src/thermo/relax_igmix_sg_class.f90`** (inherited unchanged by `relax_igmix_nasg`). New fields, all defaults bit-identical to before (verified on 670 cells, inst and finite):
- `pure_stability` (default `.true.`): gates the `pure_phase_stability` block in `pTg_relax`.
- `nucleated` (result): `.true.` when the last `apply` created a nucleation seed and kept it. It is reset at `apply` entry, set after the nucleation block, and cleared in `undo_nucleation`.
- `chem_guard` (default `.true.`): gates the two "reducing phase below `Y_small`" checks after the pTg solve.
- `solve_lv` (pure-vapor saturation solve, inside `pTg_relax`): when the line search rejects every step but |F1| < `F1_tol`, the solve now reports convergence instead of failure. This is the F3 fix.
  - Before, F1 = 0 exactly made every trial fail the strict `|F1_try| < |F1|` test, so the solve exited unconverged.
  - Only that case changes: the 670-cell regression set stays bit-identical.
  - `solve_lvg` was left alone, because its line search only runs while |F1| or |F2| ≥ 0.3.
- `pTg_relax`, dry-gas starting guess (fix 5.4, F2): `xv=pvsat/(p+pvsat)` instead of `pvsat/p`, capped at halfway between Yv0 and Yv_hi. The old line is kept as a comment. The regression set stays bit-identical.
- `pTg_relax`, after the saturation solve (fix 7, F5): when the converged state implies VF > 1 or VF < 0, the cell is now resolved at the single-phase boundary using the existing `boundary_all_liquid` / `boundary_all_vapor` (as `solve_lvg` already does), and the densities and VF are recomputed.
  - Before, VF was clipped to [0, 1], which broke mass conservation, so `check_cons` returned `RELAX_FAILED`.
  - The old clip lines are kept as comments.
  - The 670-cell regression set stays bit-identical.

**`src/amrmpsolvers/amrmpcomp_class.f90`**, `apply_relax_splitremix`. It is built on the fd296a3 single loop: no clustering, no `dissolve_stranded_gas`, no Yv_ref.

Local parameters:
- `VFrb=0`
- `VF_nuc=0.01` (same as pTg_relax)
- `remix_mode=remix_rigid`

The routine requires `relax_igmix_nasg` and calls `die` otherwise. `pure_stability` is off for the whole routine and restored at the end.

Per cell:
1. **Pure cell, |VF−0.5| < VFrb, or any Q(1:4) ≤ 0** (the last is fix 5.1): regular `apply`. A failed nucleation in a pure cell leaves it unrelaxed and is not a failure.
2. **Interfacial with VF < VF_nuc or VF > 1−VF_nuc** (and `do_nucleate`): regular `apply` on a copy, kept only if `nucleated` and ierr=OK. Otherwise go to step 3.
3. **Split-remix:**
   - `split`: the reactor holds Q per reactor volume at VF=0.5. The dominant part is computed as Q/(2·Vd), directly, to avoid cancellation.
   - Reactor call through `relax_noguard`, which switches `diss_P`/`diss_RHO` and `chem_guard` off and restores them afterwards.
   - During the reactor call `rlx%vol` is set to Vr times the cell volume, so the ledger is in physical units (fix 5.3).
   - A reactor ierr of `RELAX_VACUUM_VAPOR` is accepted as OK (fix 5.2).
   - `remix`: Q = Qa + Vr·Qr; VF = Vr·VFr, plus 1−Vr if liquid is dominant.
   - rigid+M only: whole-cell `Prelax` through `relax_noguard`.
4. **Failure:** any non-OK ierr from the reactor or the M step prints the cell indices, VF, Pjump, cell Q and reactor Q, then calls `die`.
5. **Afterwards:** pure-cell snap, or PLIC/barycenter rebuild, as in fd296a3.

**Temporary diagnostic** (lines marked `! debug`, added to chase failure F4 in §4). For each cell it records the branch (1 regular, 2 nucleation kept, 3 split-remix, 4 split-remix+M) and whether a P/T floor fired. It prints the cell's state before and after relaxation, plus the reactor input and output, when:
- the cell enters with a phase T above `dbg_Tabs = 1e4` (SI; never triggers in impact's nondimensional units), or
- relaxation raises the energy density or the hottest phase T by more than `dbg_jump = 10`.

At most `dbg_nmax = 50` cells are printed per rank per call.

Kept on in the reactor, per decision: the P/T floors, `VFratmax`, `RHOGmin`, the conservation check, and the Yv clamps.

**Debug-box instrumentation** (second session of 2026-10-23, all lines marked `! debug`). Everything is gated by the `Debug ilo..khi` keys of `input` (commented out = off). Each line carries the step number `fs%relax_diag_step`.
- `apply_relax_splitremix`: for cells in the box, sets `dbg_cell`/`dbg_i,j,k` (so `relax_igmix_sg` prints its internal trace) and always prints the `[splitremix dbg]` summary (branch, floor, Pjump, pre/post Q, reactor in/out), now with `n=` and `t=`.
- `get_dQdt` (`[TR]`): per box band cell, VF before/after, liquid and gas volumes before and after transport (before renormalisation), the four face velocities, and per face (W,E,S,N) the liquid volume, gas volume and Q1–Q4 carried into the cell (per cell volume, positive inward).
- `add_phasic_pressure` (`[PP]`): VF, PL, PG, the face accelerations −2ΔPmix/(ρlo+ρhi)/dx on W,E,S,N, the face velocities after them, divU, and the p dV energy added to each phase.
- `merge_Q` (`[MQ]`): pool adoptions of box cells (gas/liquid, partner face, VF, mass before/after).
- `examples/amrcomp_impact/src/simulation.f90%probe_boxes`: probes the `Debug` box (old hard-coded box commented out), `probe_t0=0`, and writes `Vs,Vn` as two extra columns.

**`Merge Q` input key** (`examples/amrcomp_impact/src/simulation.f90`, read into `fs%merge_on`, default `.true.`); `input` has `Merge Q : .true.`.

**`examples/amrcomp_impact/src/simulation.f90`**. The old lines are commented out, not deleted:
- The relaxation model is `relax_igmix_nasg` (instantaneous).
- The finite rates, the four safety-rule reads, `vacmix_aG`, and the three sliver/chem-skip census reductions are commented out.
- Both calls go to `apply_relax_splitremix`.
- `input` is unchanged. Its rate and safety-rule keys are now unused.

## 4. Failures seen in amrcomp_impact and their resolution

**F1, step 123, cell (960, 1024, 0):** VF=0.908, gas Yv=9.8e-4, slightly supersaturated after the pT step.
- The pTg solve converged in 2 Newton iterations, condensing 1.8% of the vapor.
- The `Y_small` condensation guard (Yv0 < 0.001) then discarded it and returned `RELAX_FAILED`.
- The old `apply_relax` hit the same guard silently.
- Fixed with `chem_guard` off in the reactor; the replayed reactor then gives ierr=0.

**F2, step 11787, cell (43, 1191, 0):** a gas sliver, 1−VF = 2.4e-9, with **negative vapor mass** (Q(8) = −1.11e-9, Q(2) = 1.34e-9, so Yv = −0.83). The liquid is in tension (pL = −0.22), and the gas is hot (TG = 1.43, TL = 0.305).
- Yv is clamped to 0, so the gas is treated as dry air.
- The reactor VF is above 0.5, so `pTg_relax` uses its dry-gas starting guess `xv=pvsat(p,T)/p` (`relax_igmix_sg_class.f90` ~line 906).
- pvsat ≈ p gives a guess of Yv = 0.9962, above the largest reachable value Yv_hi = 0.99506. `solve_lvg` then exits on its first iteration (EXIT-E) and returns `RELAX_FAILED`.
- The old path also fails on this cell, and its energy swap leaves the sliver's gas at TG = 613.
- **Fixed (5.4).** The reactor converges with `xv=pvsat/(p+pvsat)`: ierr=0, p = 0.0668, T = 0.3058, Yv = 0.378.
- The user's restart from t = 11.7 then ran past step 11787.
- The underlying negative vapor mass (5.5) is still produced by transport.

**F3, cavitation example, step 43, cell (132, 116, 0):** VF = 0.9904, pure vapor (Yv = 1).
- The cell is extreme: liquid at pL = 3.46e9 Pa and TL = 638 K, next to vapor at pG = 4.6e5 Pa and TG = 823 K.
- The reactor's pT step converges (p = 7.38e5 Pa, T = 579.8 K). Then `solve_lv` reaches F1 = 0.0 exactly at iteration 8 (p = 6.5543e6 Pa).
- Its line search accepts only a strict decrease of |F1|, so it rejects, and `pTg_relax` returns `RELAX_FAILED`.
- The failure needs the case's `-O3 -ffast-math` build. Under `-O2`, F1 ends as roundoff and the solve converges to the same state.
- **Fixed** (section 3). The reactor now gives ierr=0 with both flag sets: VF = 0.7695, p = 6.554e6 Pa, T = 577.88 K, vapor pressure equal to saturation.
- The old `apply_relax` would have hit the same failure silently, keeping the pT state with no phase change.

**F4, cavitation example, step 195, cell (132, 131, 0):** VF = 6.06e-5. **Resolved for this case** by option (a) below. Option (b) is the general fix and is not done.
- The cell enters relaxation already non-physical: e_G ≈ 1.6e20 J/kg (about 1e17 K) and e_L ≈ 9.4e15 J/kg. At step 194, TGmax was 6.0e3 K, so the blow-up happens within step 195.
- Liquid and gas have equal energy per volume to 10 digits (9.258e18 J/m³). That is the signature of the p dV transport update (amrmpcomp ~859 and ~1016) acting at PL = PG with an enormous P. The likely sequence is that the stage-1 relaxation (or the wall fold after it) left an absurd pressure, which stage-2 transport turned into energy. This is inferred, not verified.
- Monitor: TGmax jumps from 1.2e3 to 1.31e4 K at step 187. In the same step the floors fire for the first time since step 120, and their count then grows (5 → 18 by step 194).
- Suspects:
  - (a) floors in the reactor (`brhomax=1e10` disables the packing-knee guard);
  - (b) the nucleation-kept whole-cell result;
  - (c) PL ≠ PG left by the rigid remix;
  - (d) unrelated to split-remix (no baseline run with `apply_relax` to step 195 yet).
- **Diagnosed with the diagnostic** (8-rank scratch run, which aborted at step 136):
  - **49 cells entered hot from transport** (steps 16–64, TG up to 6e4 K). All went through split-remix (branch 3) and left at ≤ 567 K, which the reactor handled correctly.
  - **6 cells blew up inside relaxation, every one with a floor firing,** in branch 1 (pure gas) or branch 2 (nucleation kept). None was in split-remix. TG reached up to 1.7e7 K.
  - **Mechanism:**
    1. Vapor expands and cools to about 270–290 K, which is more than `Tctol` below saturation, so `pTg_relax` nucleates a condensation seed of liquid at the gas temperature.
    2. The case's `Tmin_liq=300` floor (`apply` stage 3) raises the shared pressure until the stiff liquid reaches 300 K at fixed density.
    3. That pressure lands almost entirely in the gas energy.
    4. The same cells repeat this every step and grow. Transport then crushes a neighbour's liquid past 1/b, and its reactor returns `RELAX_BAD_LIQUID`.
  - The old `apply_relax` is exposed in the same way: it uses the same `apply` and the same floors.
- **Option (a) tried: lower the case's relaxation floor only.** In `examples/cavitation/src/simulation.f90`, `relax_model%Tmin_liq=0.1` is set after the floor copy. `fs%Tmin_liq=300` stays, because `clean_Q`'s liquid rescue uses it as its corner temperature.
  - A first attempt that lowered `fs%Tmin_liq` instead broke initialisation. With Tmin = 0.1 K the corner density is ρ(Pmin,Tmin) ≈ 1/b = 1380.6, so the rescue raised every liquid cell to that density (1.8e11 Pa at step 0) and the timestep collapsed.
- **Option (b), not done:** make the stage-3 floor safe for tiny or freshly nucleated phases, so it does not raise the shared pressure (applies to all cases).
- **Option (c), not done:** treat "nucleation kept, but a floor fired" as a failed nucleation. This would not cover pure cells.
- F4 cell data:
  - VF = 6.0600776807208824e-05
  - Q = [5.9626987237113287e-02, 5.8440899857121732e-02, 5.6105513770477719e+14, 9.2576558729965384e+18, 2.0495222625108991e+07, 4.0884914015628070e+07, −4.5766866085638651e-12, 5.8408516492753006e-02]

**F5, cavitation example, t = 2.436e-5 (step 1771), cell (127, 128, 0):** a vapor sliver, 1−VF = 3e-6. **Fixed (fix 7).**
- The cell is extreme: liquid at 3.1e9 Pa and 1194 K, next to dense vapor at 9.4e7 Pa and 2103 K.
- The reactor's pT step gives p = 1.69e8 Pa and T = 1133 K, far above psat (4.95e7 Pa), so all the vapor should condense.
- `solve_lv` has no boundary handling. It converges to a saturation state with VF = 1.0146, the clip breaks mass by 1.2%, and the result is `RELAX_FAILED`. This happens under `-O2` as well as fast-math.
- After the fix the reactor ends as pure compressed liquid: VF = 1, p = 6.307e7 Pa, T = 1110.2 K. That matches a hand computation of the all-liquid state (ρ = 622.3, stable since p > psat(T) = 4.94e7 Pa).
- The whole cell relaxed by the regular path lands on the spurious upper root again (p = 1.72e9 Pa), as in 5.6.
- **Note:** liquid at 1100–1200 K and 3e9 Pa is far beyond the water critical point and the NASG fit range. The late-time states of this case are suspect in their own right, presumably compression waves from the wall.
- F5 cell data:
  - VF = 9.9999705109486214e-01
  - Q = [1.1480747721175308e+03, 2.8476175562529011e-04, 3.9767372005702114e+09, 1.1668181485617006e+03, −1.1560525966761974e+06, −1.0790011098319679e+06, −2.4655417027978171e-13, 2.8476175562529011e-04]

Cavitation EOS (SI):
- Liquid (NASG): gamma 1.1173, pinf 5.0017e8, b 7.2341e-4, cv 3821.2, q −1.1711e6, qp 0, `brhomax` = 1e10.
- Vapor: gamma 1.4963, cv 929.06, q 2.1442e6, qp 1.9413e4.
- Air: gamma 1.4, cv 718, q = qp = 0.
- Relaxation settings: p_cav = −1e8, Tctol = 5, RHOGmin = 0.01, VFratmax = 10.
- Floors: Pmin_liq = −0.9·pinf, Tmin_liq = 300, Pmin_gas = 1e-4, Tmin_gas = 0.1.
- Pjump = 0.
- To reproduce the failure, the replay must be compiled with the case flags: `-O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer`.

F3 cell data:
- VF = 9.9041599391982282e-01
- Q = [1.2446977466304766e+03, 1.1607543601731204e-02, 1.6227547621057034e+09, 3.3764527263874777e+04, 7.6772257456853939e+05, −2.3033680141200111e+06, 7.8012121615202520e-14, 1.1607543601731204e-02]

Replay data. Impact EOS, all nondimensional:
- GammaA = 1.3; CvA = pG2/(GammaA−1) with pG2 = 1/(GammaA·Mach²), Mach = 1.95.
- Vapor: gamma 1.3115622, cv 1.0448401, q 1.5987106, qp 2.8517267.
- Liquid: gamma 1.4592522, pinf 54.360449, b 0.0036532648, cv 2.084228, q 0, qp 6.211819, `brhomax` = 1e10.
- Relaxation settings: p_cav = −0.8·pinf, Tctol = 0.1, RHOGmin = 0.
- Floors: Pmin_liq = −0.98·pinf, Tmin_liq = 0.1, Pmin_gas = 1e-4, Tmin_gas = 0.1.

Per-cell data:
- **F1:** cell VF = 9.0830017296304322e-01, Pjump = 7.1547178066394635e-05.
  - Q = [1.2414446358199905e+02, 2.2437799076370152e-01, 7.9581893485517199e+01, 3.2575651512548723e-02, −1.2232179246401974e+02, 4.7601840982517851e-01, 3.9523229889110235e-19, 2.1980103268593284e-04]
- **F2:** cell VF = 9.9999999755450442e-01, Pjump = 3.5575093821488507e-07.
  - Q = [1.1054344530819478e+02, 1.3420574648382725e-09, 1.0268502315941987e+02, 1.2977223204062312e-09, 3.4805016357932310e+01, 8.6800390406628878e+01, 2.1082562562271089e-18, −1.1137933085331929e-09]

## 5. Fixes and open items

Items 1–4 are **implemented** (2026-10-23, to clear F2 in amrcomp_impact). Items 5–8 are open.

1. **Zero- or negative-phase-mass interfacial cells.** *Implemented:* such cells go to the regular `apply` before splitting.
   - Problem: Q(1) ≤ 0 or Q(2) ≤ 0 (any of Q(1:4) ≤ 0) with VF interfacial. The reactor inherits it, `pT_relax` returns `RELAX_DEGENERATE`, and the routine dies. The old `apply` left such cells untouched and only counted them in `acc(7)`.
   - Proposal: route `any(Q(1:4).le.0)` to the regular `apply` before splitting.
2. **`RELAX_VACUUM_VAPOR` (8) from `activate_chem`.** *Implemented:* accepted as OK after the reactor call.
   - Problem: dry gas whose saturation pressure is below `pv_min` cannot hold vapor, so there is no phase change and the pT state is kept and consistent. It is not a failure, but the routine dies on it.
   - Proposal: accept it like OK.
3. **Floor-energy ledger.** *Implemented:* `rlx%vol=Vr*vol` around the reactor call.
   - Problem: the floors stay on in the reactor, but `acc(6)` uses `relax%vol` = cell volume while the reactor Q is per reactor volume. The rescue monitor's floor energy is overstated by 1/Vr.
   - Proposal: set `rlx%vol=Vr*cell_vol` around the reactor call.
4. **pTg dry-gas starting guess** (`relax_igmix_sg_class.f90` ~906, F2). *Implemented:* `xv=pvsat/(p+pvsat)` (Pjump still omitted, as the existing comment there requires), capped at halfway between Yv0 and Yv_hi, with the old line kept as a comment. The new F2 cell from the user's rerun now gives reactor ierr=0, Yv = 0.378. The 670-cell regression set stays bit-identical.
   - Problem: the current guess assumes the total pressure stays at p. In a liquid-dominated cell the stiff liquid holds the gas volume nearly fixed, so the dry gas keeps its partial pressure and the vapor adds pvsat on top.
   - Proposal: use `xv=pvsat/(p+pvsat)` and also cap the guess below Yv_hi.
   - This changes the guess for every caller, but only noticeably where pvsat is comparable to p.
5. **Negative vapor mass from transport** (F2 input). Options:
   - (a) leave it; the relaxation clamp already treats such gas as dry air;
   - (b) clip Q(8) into [0, Q(2)] before relaxation (not conservative);
   - (c) fix species positivity in transport, which is the root cause and a separate job.

6. **Old path: `solve_lv` converges to the spurious upper root** (found with F3; affects `apply_relax`, not split-remix).
   - Example: the F3 cell relaxed whole by the regular `apply`, with liquid pT pressure 3.07e9 Pa after pT relaxation. `solve_lv` starts from that pressure and converges to p = 1.49e9 Pa, T = 597.8 K. There p = pvsat(p,T) holds on the upper root; the physical saturation pressure at that T is 8.08e6 Pa.
   - The residual ln pvsat(p,T) − ln p has two roots in p, and only the lower one is physical. `get_psat` already brackets the lower one on purpose; see its comment.
   - The split reactor started from 7.4e5 Pa and found the physical root.
   - So `apply_relax` runs of strongly compressed near-pure liquid next to vapor may have produced cells at about 1.5e9 Pa in "saturation".
   - Possible fixes:
     - (a) start `solve_lv` from the lower root, e.g. `get_psat(T)` at the pT temperature, rather than from the pT pressure;
     - (b) reject a converged p above the residual's minimum `pmin`, which `get_psat` computes in closed form, and retry from below it;
     - (c) bracket in p on the monotone side, as `get_psat` does.
   - Not implemented; check first whether earlier cavitation runs contain such cells.

7. **Negative vapor mass in fresh cavitation seeds, then an absurd gas rescue** (found while tuning a gentler cavitation input; not split-remix specific).
   - **How it happens:** a seed nucleated under strong tension (`p_cav = −1e7`) is expanded to near-vacuum by mechanical relaxation. Its gas mass is then so small that the next transport step drives Q(8) negative, by up to 300 times |Q(2)|.
   - `clean_Q`'s gas rescue (Pmin_gas, Tmin_gas) then sets the gas to its corner state, pure carrier air at Tmin_gas = 0.1 K and Pmin_gas = 1e-4 Pa, i.e. ρG = 3.48e-6 kg/m³ and e = 71.8 J/kg. It leaves Q(8) negative.
   - The reactor cannot relax 0.1 K gas next to 350 K liquid under tension, so it returns `RELAX_FAILED`.
   - **Seen in** scratch trials B, D and E of the cavitation example. These used pulse (r0 = 1 mm, U0 = 2.33e5) or (r0 = 3 mm, U0 = 3.9e4 or 1.95e4) with `p_cav = −1e7`, and failed within 0.2 μs.
   - **Example cell (trial D):**
     - VF = 9.9446213426154995e-01
     - Q = [9.5874639110325199e+02, 1.9282262320508550e-08, 3.0986406189401692e+08, 1.3844664346125139e-06, −2.1667987760472392e+03, −6.5004797224548920e+03, −2.2130145706457889e-18, −6.0592599856273731e-07]
     - Pjump = 0
   - **Possible fixes:**
     - species positivity in the gas transport of tiny or near-vacuum gas;
     - a rescue that keeps Q(8) consistent (clip to [0, Q(2)]);
     - a rescue corner at a physical temperature rather than 0.1 K.
   - Related to 5.5.
8. **pTg Newton (`solve_lvg`) gives up when it starts from near-vacuum gas** (not split-remix specific).
   - **How it happens:** a gas sliver at about 2–10 Pa next to 350 K liquid (psat = 4.4e4 Pa) needs the pressure to rise about 20,000-fold, but the step limiter allows ×1.5 per iteration. Meanwhile Yv runs to its ceiling, `hi_lock ≥ 3` fires, neither single-phase boundary is valid, and the solve is rejected with `RELAX_FAILED`.
   - The true equilibrium is easy: the gas fills with saturated vapor at psat, which takes about 0.14 kg/m³ of the reactor's 482 kg/m³ of liquid.
   - **Seen in** scratch trials A (r0 = 1 mm, U0 = 1.17e5, `p_cav = 0`, died at 18.2 μs) and F (r0 = 3 mm, U0 = 1.95e4, `p_cav = 0`, died at 15.7 μs). Both hit stray seeds far from the centre, where reverberations put the liquid in tension.
   - **Example cell (trial A):**
     - VF = 9.9997431857407970e-01
     - Q = [9.6425100214054237e+02, 3.9420947612736675e-10, 3.1161505265459567e+08, 7.7511412844409242e-04, 9.4471073182581051e+01, 2.3776443938408867e+01, 1.6076518117920629e-18, 8.4730944668768291e-11]
     - Pjump = 0
     - The reactor starts the chemistry at p = 1.74 Pa after the pT step.
   - **Possible fixes:**
     - when p ≪ psat(T), start the Newton from p = psat(T) plus the dry-gas partial pressure (the same idea as 5.4);
     - only trust `hi_lock` once the pressure step is no longer limited.
   - **Context:** the cavitation example's input was left as is, the strong pulse with `p_cav = −1e8`. Gentler pulses either produce no visible cavity (an impulsive pulse radiates its energy acoustically within about r0/c) or a diffuse cavitation cloud, and most die on 7 or 8. A physical single-cavity demo would need an initial seeded bubble instead of nucleation from pure liquid. Not done.

Things checked and needing no guard:
- **`chem_guard` off is safe.**
  - The Newton limiter keeps Yv inside (Yvmin+fd_eps, Yv_hi), so vapor mass cannot go negative. Running out of liquid is resolved as the all-vapor boundary case.
  - The evaporation half of the guard cannot trigger at VF=0.5.
  - Residual risk: an equilibrium Yv below about 1e-7 cannot be reached, so the solve fails to converge and the run dies. Not realistic at impact conditions, where Ysat ≈ 1e-3.
- **Remixed states are always physical** (derived). T_mix is a heat-capacity-weighted average of the parts' temperatures, and v_mix is a mass-weighted average of their specific volumes. So T > 0 and v > b are preserved.
- **Rigid remix leaves PL−PG ≠ Pjump.** No downstream code assumes equality:
  - the pressure force uses the VF-weighted mixture pressure (amrmpcomp ~795);
  - the p dV work uses each phase's own pressure (~859).
- **Existing behavior, unchanged:**
  - the snap deletes a remaining minor phase below VFlo, a loss of about VFlo·ρ;
  - an unclamped Yv in part a survives the remix;
  - Pjump comes from cell curvature (a risk only for noisy slivers with large σ).

## 6. Later

- **Piston remix:** design the constrained equilibrium. In particular, for pure vapor, allow part a to condense while keeping it thermally insulated.
- **Repeated-step test:** check that repeated split-remix steps on an isolated cell converge to the whole-cell equilibrium, which is a fixed point.
- **Finite-rate class:** find out why `relax_igmix_nasg_finite` at `huge` rates misses equilibrium.
- **eos_relax Makefile:** it does not build the finite class (missing `parallel`).

## 7. Untested

- **rigid+M** (`remix_mode=remix_rigidM`) has never run in the flow solver. It was only studied offline, where it was poor.
- **`VFrb>0`** has never been used.
- **The near-pure "nucleation first" branch in the flow solver** has not been confirmed to trigger in a run, because there is no counter. The flag itself is unit-tested.
- **Pure cells with `pure_stability` off.** The outcome is argued equivalent (only the ierr code and ledger counts differ) but has not been checked in a run.
- **Physics quality.** No comparison yet of split-remix against the old `apply_relax` on real fields: Yv, T extrema, conservation monitors, overshoots at thin ligaments.
- **Run coverage.**
  - amrcomp_impact passes step 11787 after fixes 5.1–5.4 (user restart from t = 11.7). Its monitors have not been analysed yet.
  - The cavitation example completes.
  - Neither has been compared against the old `apply_relax`.
- **cav_tube** has never run with split-remix. Its case still uses the finite class and `apply_relax`.
- **Cost** has not been measured: the extra discarded whole-cell call in the near-pure branch, plus one reactor call per interfacial cell.
- **Conservation** has not been checked in a run: the floors (still on in the reactor) and the pure-cell snap of tiny remixed phases are both non-conservative.
- **Fixes 5.1–5.4 are in.** Only 5.4 has been exercised by a failing cell (F2). 5.1 and 5.2 are defensive: no run has hit them yet.
- **Items 5.5–5.8 are open.**

## 8. Options identified but not tried

- **Piston remix** (part a isentropic and mechanically coupled). Python prototype only; for pure vapor it needs condensation from part a.
- **Channel-wise split:** relax M over the whole cell, and T and C locally in the reactor. This follows the acoustic vs diffusive argument and the M→T→C channel order.
- **Isobaric reactor:** hold p fixed and absorb the volume change in part a.
- **Other reactor ratios** than equal volume:
  - equal heat capacity;
  - an effusivity-weighted ratio;
  - a knob φ* so that the reactor VF is max(VF, φ*), which interpolates between whole-cell and split.
- **dt-aware split fraction**, e.g. part b sized from a diffusion length √(D·dt). Today the per-step relaxed fraction f = Vm/Vd depends on the step count, not on time.
- **Finite rates with split.** Pelanti's rates are per unit volume, so the reactor cuts the effective rate by Vr. A scaling rule would be needed.
- **Repeated-step test** on an isolated cell, converging to the whole-cell equilibrium.
- **Distinct return codes** instead of `die`: classify "declined but consistent" outcomes, the alternative to 5.2 and to the F1 choice.

## 9. How to resume

**Uncommitted changes.**
- Already modified before this session: `relax_igmix_sg_class.f90`, `amrcomp_impact/src/simulation.f90`, `amrcomp_impact/input`, `thermorelax_class.f90`, `relax_igmix_nasg_finite_class.f90`.
- `amrmpcomp_class.f90` was clean at session start.

**Changed or created this session:**
- `src/thermo/relax_igmix_sg_class.f90`: `pure_stability`, `nucleated`, `chem_guard`, the `solve_lv` line-search convergence fix (F3), and the single-phase boundary resolution after the saturation solve (fix 7, F5).
- `src/amrmpsolvers/amrmpcomp_class.f90`: `apply_relax_splitremix` and its binding, fixes 5.1–5.3, and the `! debug` diagnostic (§3; strip it once no longer needed).
- `examples/amrcomp_impact/src/simulation.f90`: model switched to `relax_igmix_nasg`, the finite-only lines commented out, the calls swapped.
- `examples/cavitation/src/simulation.f90`:
  - Both `apply_relax` calls swapped to `apply_relax_splitremix`, with the old lines commented out. The case already used `relax_igmix_nasg`.
  - `relax_model%Tmin_liq=0.1` (F4, option a); `fs%Tmin_liq` stays 300 for the `clean_Q` rescue.
  - With fixes 6 and 7, it runs to completion with split-remix.
- New: `examples/eos_relax/relax_batch.f90`, `examples/eos_relax/split_reactor.py`, this file.

**Changed or created in the second 2026-10-23 session (overshoot tracking):**
- `src/amrmpsolvers/amrmpcomp_class.f90`: the debug-box instrumentation of §3 (`[TR]`, `[PP]`, `[MQ]`, `[splitremix dbg]` with `n=`). The gas-strain guard was added and then fully removed; nothing of it remains.
- `examples/amrcomp_impact/src/simulation.f90`: `Merge Q` key; `probe_boxes` on the `Debug` box with `Vs,Vn`.
- `examples/amrcomp_impact/input`: `Merge Q : .true.` added (after `Relax vacuum pressure ratio`). The user's run controls at the end of the session: `Max time 10.9`, `Restart from 1.08000E+01`, output and checkpoint every 0.01.
- New, untracked: `examples/amrcomp_impact/overshoot_tools/` (§10.8).

**Scratch tools are gone with the session; rebuild them as needed** (the overshoot tools of §10.8 are kept in the repo):
- **Replay driver** (`test_reactor.f90`-style): sets up the case EOS and relaxation settings (§4 has both cases), calls `relax_igmix_nasg%apply` on one cell or reactor with `dbg_cell=.true.`, and prints the phasic states. Section 4 has the cell data.
- **Build:** compile `precision_dp, string, messager_stub, parallel stub, eos_relax/amrvof_class, material, thermorelax, ideal_gas, stiffened_gas, nasg, igmix, relax_igmix_sg, relax_igmix_nasg` from a clean directory, with a bash script:
  - zsh does not split `$LIB`;
  - stale `.mod` files in `eos_relax/` shadow the current ones;
  - reproducing a failure from a run can need the case's flags `-O3 -ftree-vectorize -ffast-math -funroll-loops -fomit-frame-pointer` (F3 did).
- **Regression check:** 670 cells (pure, near-pure and interfacial VF, S1–S4 states plus a cavitating liquid, models pTg and p) through `relax_batch` (inst and finite), compared byte-for-byte against a pre-session baseline. Every relaxation-class change so far kept it bit-identical. Rebuild a baseline from `git stash` or HEAD if needed.
- **Frame analysis** (used for the cavitation dust study): yt 4.3 in `~/venv/combustion`.
  - Load each frame with `yt.load('amrviz/<case>/plt.nga2.cell.NNNNNN')` and build `ds.covering_grid(level=ds.index.max_level, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions*2**max_level)`.
  - Read fields as `('boxlib', name)`. Available: VF, RHOL, RHOG, PL, PG, TL, TG, Yv, IL, IG, U, V, W, Umag, Mach, C, beta, visc, cluster_idx, stranded.
- **Runs:**
  - Run cases in a scratch copy of the case directory, so your `monitor/`, `amrviz/` and `restart/` stay untouched.
  - Use `mpiexec -n 8 -x GFORTRAN_UNBUFFERED_ALL=1 ...`. Without it, an abort loses the failing cell's printout, because stdout is buffered.
  - The user prefers to start amrcomp_impact runs themselves.

**Next task: change the end state of the mechanical exchange, not its speed.** The overshoot tracking is done (§10). The user starts a new chat to explore two ideas:
1. **Equilibrate the gas against the liquid pressure at the interface**, not the cell-average liquid pressure. In a collision cell the cell average carries the water-hammer pressure of the liquid bulk (pL 8–40); the gas only touches the liquid at the interface. This changes the fixed point.
   - Prior art: the reconstruction jump Jr of branch `thermo_gpu_friendly` (commit 9cc125ec; `examples/amrcomp_impact/overshoot_rules_gpu_friendly.md` §5f). Per mixed cell, each phase's pressure is carried from its barycentre (CL/CG) to the PLIC polygon centroid with a minmod-limited gradient over face neighbours holding that phase; `Jr = dG − dL`, clamped to [0, PL−PG−σκ]; the relaxation solves with `Pjump+Jr` and gives `Jr·dVF` back to the liquid. Reported effect there: worst sliver rhoG 14.9× → 10.9× of its neighbours.
   - It reads neighbours but writes only its own cell.
2. **A mechanical relaxation rate per unit time**, e.g. set by the acoustic time of the gas pocket, so the flow can relieve the compression before the gas absorbs it. It does not change the fixed point, but the squeezing events last only 20–50 steps. Note: the user rejected the cell-local skip guard (§10.6), which acted similarly; the Pelanti finite class misses equilibrium at `huge` rates (§2, §6).
- **Test bed:** restart from the user's checkpoint t = 10.80 and run to 10.90 with the debug box i = 44–52, j = 995–998; baseline numbers in §10.3–10.4 (box max rhoG 298 at step 10900, (47, 998) 142, (46, 998) 252). Cell (43, 997) from t = 10.82 (§10.6) is the sliver case.
- Runs from the user's checkpoints reproduce the user's frames to round-off over this window (§10.1).

**Other open decisions:** fixes 5.5–5.8; the channel-wise split (§8; its mechanical part is whole-cell, so it cannot help the overshoot either, §10.5); stripping the `! debug` diagnostics; floor option (b) from F4; the sliver dissolve of §10.7; `merge_Q` heating the gas (§10.2).

## 10. rhoG/pG overshoot investigation (second session of 2026-10-23)

All numbers are nondimensional (impact case). K_L = ρL·cL² is the liquid bulk modulus (≈ 170–180 here). "Gas strain" is the relative loss of gas volume at fixed gas mass.

### 10.1 Method
- **Runs:** scratch copies of the case directory, restarted from the user's checkpoints (`restart/impact_relax_pTg/1.08000E+01`, `1.08200E+01`, `1.10000E+01`), `input_test` with the `Debug` box set. 8 ranks, `GFORTRAN_UNBUFFERED_ALL=1`. 100 steps take about 85 s.
- **Reproducibility:** `random_initialize` seeds from `/dev/urandom`, and every cell that turns mixed in relaxation gets a random PLIC normal, so two runs from one checkpoint drift apart once that happens. Locally the cells still matched the user's frames: the 10.80–10.90 box matched to ≤ 3e-11; the 11.00–11.04 run diverged elsewhere from step 11007 but cell (49, 997) still matched.
- **Stage accounting (RK2):** stage 2 restarts from the start-of-step state (`Q = Qold + dt·dQdt(Q*)`), so stage-1 relaxation and pressure updates affect a step only through the velocity. Per step: transport = S2a − A0, `build_plic` = S2b − S2a, relaxation = S2c − S2b, pressure step = S2d − S2c (probe labels of `probe_boxes`). The first analysis of cells 49/50 summed both stages, which counted relaxation about twice; its conclusions stand.

### 10.2 Cells (49–50, 997), t = 11.00–11.04
- pG and TG of both cells leave their face neighbours' range between 11.01 and 11.02 (pG 1.02/1.08 against at most 0.74; TG 0.83/0.86 against at most 0.64), while the cells are pure gas (liquid enters at steps 11015 and 11017).
- Drivers: converging flow (divU −35 to −53 for 20+ steps) through transport and the p dV work of the pressure step. Relaxation contributes little.
- `merge_Q` heats the gas. A liquid claimant (VF < `merge_VFlo`) adopts its partner's liquid velocity; the relative kinetic energy is thermalised and split by mass, so the gas gets ~85% ([amrmpcomp_class.f90] `merge_Q` pass 2, "single-velocity collapse"). It matched the measured gas-energy jump to all printed digits; over steps 11015–11019 it added ~15% to cell 50's gas energy (|uL − u| ≈ 0.67).
- With `Merge Q : .false.` (key added for this test): cell 50 pG 1.083 → 0.928 and TG 0.856 → 0.753 at step 11020, still above its neighbours; cell 49 unchanged (transport made up the difference).

### 10.3 The collision mechanism (box i = 44–52, j = 995–998, t = 10.80–10.90)
- The overshoot sits in row j = 998, under the liquid at j = 999: rhoG 138 in (47, 998) at step 10870, 252 in (46, 998), 298 in (45, 998) at step 10900 (pG 16.9), against 1–20 in gas cells two cells away.
- **(47, 998) is where two liquid streams meet:** east face u = −0.46 (liquid inflow, ρL ≈ 127), west face +0.25 → divU ≈ −65, the cell loses 6.5% of its volume per step. Transport strains both phases equally (VF is renormalised from the two phase volumes, `get_dQdt` divergence block), so pL jumps 1 → 11 per step (≈ 0.065·K_L); relaxation gives that compression to the gas.
- **The collision moves west:** (46, 998) is crushed next (rhoG 187 before it closes, VF → 1 at step 10893) and pushes gas at rhoG 100–250 through its west face into (45, 998), whose gas mass grows from 1.25 to 13.
- **Liquid pressures** around the collision reach 15–40 (pure liquid (44, 997) at pL 40), consistent with a water-hammer ρL·cL·Δu ≈ 132·1.1·0.25 ≈ 36.
- **Gas escape** (inference, not checked separately): the gas leaves only with the single, density-weighted face velocity (`get_face_velocity`), and its pressure enters the momentum only as (1−VF)·pG in the VF-weighted mixture pressure (`add_phasic_pressure`). Gas in a liquid-rich cell moves with the liquid.
- Cumulative Δ ln rhoG over the run (transport / plic / relaxation / total): (45, 998) +0.25 / +0.17 / **+4.18** / +4.61; (47, 998) −2.56 / 0 / **+6.26** / +3.70; (44, 997) −1.67 / −4.25 / **+8.72** / +2.80; (46, 998) +2.01 / −0.75 / +2.86 / +4.12.
- Some sliver cells show huge cancelling transport/plic pairs (±53 in (46, 997)): transport gives an absurd rhoG to a tiny gas volume and `merge_Q` repairs it. A separate sliver issue.

### 10.4 Relaxation strain law and what separates the cells
- **Per split-remix call** the gas strain is (pL − pG − Pjump)/K_L · min(VF, 1−VF)/(1−VF). Checked on all 2155 stage-2 split-remix calls with gas on both sides: median measured/predicted 1.05–1.29 in every VF bin (the excess is probably the thermal part of pTg).
- **Where the pre-relaxation over-pressure comes from** (medians): VF < 0.5: 1.12, of which −0.02 left over from the previous step and +1.18 added by transport; VF 0.5–0.7: 0.15 left over, +1.27 transport; VF 0.7–0.9: 0.93 / +0.31; VF > 0.9: 1.43 / +0.07. In gas-minor cells it is part a's stored compression.
- **High-VF cells (0.8 < VF < 1) split into two groups, same VF:**

| | Mean divU | Mean pL − pG before relaxation | Relaxation Δ ln rhoG per step | rhoG |
|---|---|---|---|---|
| No overshoot: (47–52, 997), (45, 995–996), (52, 996) | −1.6 to +1.3 | −1.6 to 1.7 | ≤ 0.01 | 1–6 |
| Overshoot: (44, 996–997), (45, 997–998), (46–47, 998) | −4.7 to −17.6 | 7 to 31 | +0.04 to +0.20 | up to 40–300 |

  VF is not the discriminant; squeezing and the liquid over-pressure are.

### 10.5 Why split-remix cannot remove it (verified, plus derivation)
- **Test** (`overshoot_tools/replay_repeat.f90`): freeze one cell (no transport) and apply split-remix 60 times; compare with one whole-cell `apply`. rhoG after 1 / 3 / 10 / 20 / 40 / 60 calls:
  - (45, 998), step 10880, VF 0.88, rhoG 42.4: 45.3 / 50.1 / 60.4 / 66.4 / 69.2 / 69.5; whole-cell 69.5.
  - (47, 998), step 10865, VF 0.84, rhoG 61.5: 68.3 / 78.3 / 96.2 / 103.5 / 105.5 / 105.6; whole-cell 105.6.
  - (43, 997), step 10830, VF 0.995, rhoG 11.8: 13.2 / 16.1 / 31.4 / 72.4 / 206 / 287; whole-cell 314.
  - pG and pL converge to the whole-cell values as well.
- **Why:** within one call only part b's liquid presses on the gas, but the additive remix merges a and b into one liquid state, and the next call's part b inherits part a's compression. Nothing remembers which liquid has relaxed, so after about VF/(1−VF) calls the gas has received all of it. Only stage 2's relaxation is kept each step, so the 20–50-step squeezing events give 20–50 calls.
- **Any remix has the same end state (derived):** if part b is cut from the cell's own liquid and the remix conserves mass, energy and volume and is the identity when the reactor does nothing, then the whole-cell equilibrium is a fixed point (the reactor does nothing there), and the only one (anywhere else the reactor moves part b and the remix carries the change). rigid+M and piston reach it faster; slower variants only delay it. Not tested numerically for rigid+M or piston.
- In these cells the whole-cell equilibrium *is* the overshoot: the gas at the liquid's collision pressure. A fix must change the end state (§9 next task).

### 10.6 Gas-strain guard (tried, rejected, removed)
- **v1:** skip the whole relaxation of a gas-minor (VF > 0.5) interfacial cell when the predicted strain s = (pL − pG − Pjump)/K_L > s_max = 0.02 (a priori, closed form). 10.80–10.90 against baseline: box max rhoG 298 → 31.5; divU in (45–47, 998) went to ≈ 0 within 5–10 steps (the liquid held pL 20–57 and slowed the collision); domain max rhoG at 10.9 298 → 65 (777 → 147 at step 10887); totals unchanged to 4–5 digits; ≈ 56 guarded cells per call.
  - Side effect 1: (41, 1051), guarded at VF 0.56 while still squeezed, reached pL 179, ρL 186 (pure liquid around it 36–77), decaying after ~25 steps.
  - Side effect 2: expanding cells with near-vacuum gas next to compressed liquid, e.g. (11, 1056), VF 0.996, pG 0.016, pL 8.3, divU ≈ +5, also tripped it; with no heat exchange the gas cooled (TG 0.119 against 0.285) and TGmin touched the 0.1 floor; gas rescues 124 → 153.
- **v2:** added divU < 0 (face divergence). Removed side effect 2 only partly (TG-floor steps 29 → 13), and missed cells that expand slightly while their liquid is still compressed: (43, 997) from t = 10.82 had divU +0.4 to +3 for 22 steps with s = 0.05–0.09, so relaxation pumped its gas (+0.07 to +0.13 in ln rhoG per step); (42, 1048) was crushed to rhoG 283 after the collision rebounded.
- The user rejected the guard; it is removed from the code and from `input`.
- Also seen in (43, 997): when VF passes 0.99 (step 10826), `merge_Q` pools its gas with a neighbour holding crushed gas, jumping rhoG 4.1 → 12.9.

### 10.7 Cell-local rules of branch `thermo_gpu_friendly`
Documented in `examples/amrcomp_impact/overshoot_rules_gpu_friendly.md`; code in commit 9cc125ec (only rules 1 and 2 are on `thermo`, in `relax_igmix_nasg_finite`). Against this mechanism:
- **5h, sliver dissolve** (gas fraction < 0.01 and pG < r·(pL − Pjump), r = 0.5 → gas folded into the liquid, VF = 1): would remove the sliver cases ((43, 997), (44, 997–998), late (46, 998), (42, 1048)); not the squeezed cells with more gas. Cost: a trace of air becomes liquid mass. The user may want it ported into `apply_relax_splitremix`.
- **Rule 2** (sliver gas along its isentrope): the gas still ends at the liquid pressure, only hotter; partial.
- **Jr**: lowers the equilibrium pressure; partial in its own tests; it is idea 1 of §9.
- **5g, rule 1, gas-density floor, positivity gate, rules 4/5/5e:** not applicable (TG ≈ TL already; dense gas; other quantities).
- Clustering and pooling rules were out of scope at the user's request.

### 10.8 Tools (`examples/amrcomp_impact/overshoot_tools/`, untracked)
- `trace.py`: `load(run_dir)` returns probe states per (i, j, step, stage) and the parsed `[TR]`, `[PP]`, `[MQ]` and `[splitremix dbg]` records. Run with `python -I` from the tools directory or add it to `sys.path`.
- `cellwin.py <amrviz case dir> i j [hi] [hj]`: field window around a cell for every frame. `pocket.py`: gas region connected to a cell, binned by VF. `neigh.py`: mixed cells against their pure-gas face neighbours, binned by VF.
- `replay_split.f90` (whole-cell vs one split-remix call) and `replay_repeat.f90` (60 repeated calls), case EOS and relaxation settings; input `cells.dat` lines `tag VF Pjump Q(1:8)`. Build with `bash build.sh` / `bash build_repeat.sh` (case flags, sources from `src/` and `examples/eos_relax/`).
