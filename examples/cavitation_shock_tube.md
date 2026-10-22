# Cavitation tube and dodecane shock tube: session handoff (2026-10-19 … 22)

Status: cav_tube is the active case (sections 3–5); dodecane_shocktube is parked (section 6).

## 1. Repository state

- Branch `thermo` at `b9fdeb80`. The GPU-friendly variant is `9cc125ec` on branch `thermo_gpu_friendly` (notes in `examples/amrcomp_impact/overshoot_rules_gpu_friendly.md`, untracked).
- Uncommitted changes on `thermo`:
  - `src/thermo/relax_igmix_nasg_finite_class.f90`:
    - VF-based chemical rate scaling (`nu_scale`, `nu_eff`) removed (Pelanti's nu per unit cell volume again);
    - cavitation seed cap (tension room, section 3.1);
    - condensation seed cap: latent heat may warm the cell by at most half its subcooling (fix 1, section 3.2);
    - `finite_chem`: liquid that no longer fits its volume (b·(Q1+dm) ≥ brhomax·VFeq) evaporates completely; a condensation step that would overpack the liquid is rejected; non-finite dm/VFeq are rejected (fixes 2 and 5);
    - `pTg_relax`: if pT_relax leaves VF outside [VFlo,VFhi], return before the chemical step (fix 4);
    - `dbg_cell` debug prints, and a `BADEXIT` check in `apply` around `sg_apply` (prints the entry state when a call returns NaN or >10× the internal energy), all marked `! debug`.
  - `src/thermo/relax_igmix_sg_class.f90`: Stage 3 floor skipped when the liquid is past the NASG packing knee (fix 3).
  - `src/amrmpsolvers/amrvof_class.f90`: face-slab field `subVF_width` removed again (2026-10-22); back to HEAD's half cells. Section 4.2.
  - `src/amrmpsolvers/amrmpcomp_class.f90`: unchanged (face-velocity memory was tried and removed, 2026-10-22; section 4.5).
  - `src/libraries/thermorelax_class.f90`: optional `nucleated` argument removed from `get_chem_transfer`.
  - `examples/amrcomp_impact`: `Chemical rate scaling` key removed.
  - Untracked: `examples/cav_tube/`, `examples/dodecane_shocktube/`, this file.
- `relax_igmix_sg_class` keeps its own copy of the nucleation block; fixes 1, 2, 4, 5 are only in the finite class.
- `examples/eos_relax` cannot currently build the finite class: nothing provides `module parallel` (`rank` is used in a debug block). A one-line stub works (used in the replay driver, section 3.3).
- `examples/dodecane_shocktube` was built before fixes 1–5; rebuild before running.

## 2. Common driver features (both `src/simulation.f90`)

1D SI cases (ny=nz=1, periodic y/z one cell thick, Max level 0), `relax_igmix_nasg_finite` with an `igmix` gas [vapor, air] (air is a zero-mass carrier, so with pure vapor the model is Pelanti 2022).

- **Time step:** `Max dt` is the initial dt; `adjust_dt` from step 1.
- **Artificial viscosity:** `Artificial bulk viscosity coefficient` (Cartif, 0 = off; both inputs use 5). `Artificial viscosity coefficient` is the shear-to-bulk ratio, active only when Cartif > 0.
- **Nucleation:** `Cavitation pressure threshold` [Pa] (pure liquid nucleates when pL < p_cav and pL < psat; 1e30 = at saturation) and `Condensation temperature tolerance` [K] (vapor nucleates when TG < Tsat − Tctol).
- **Safety rules and clustering:** keys present, all off. Cell merging (`merge_Q`, `merge_on=.true.` by default) is on and has no counter.
- **Debug line probe:** `Probe ilo/ihi`, `Probe start time`, `Probe stride` → `probe_r<rank>.csv` (VF, Q, face velocities, cell velocity, phasic thermo, `subVF_lo/hi`, PL/PG at 9 stage points per step: `A0_start`, `S1a_transport`, `S1b_plic`, `S1c_relax`, `S1d_end`, `S2a_transport`, …, `S2d_end`).
- **Relaxation debug box:** `Debug ilo..khi` prints the full relaxation trace for those cells.
- **Plot scripts** (`plot_pelanti_compare.py` in each folder): `python3 plot_pelanti_compare.py --p | --pT | --pTg [--case <amrviz dir>] [--paper <pdf>]`. Default case dir `amrviz/<case>_relax_<type>`; output `dodecane_shocktube_<type>.pdf` / `cav_tube_<type>.pdf`. Pelanti's curves come from the Zotero PDF's vector paths, calibrated on tick marks. Level-0 plotfiles only (no AMR).

## 3. cav_tube: relaxation failures found and fixed

### 3.1 Setup (current `input`)
- x in [0,1] m, 2000 cells (dx=5e-4), p=1e5 Pa, T=353 K, u=∓2 m/s, P22 Table 3 water, transmissive ends, dt=1.5e-7 (CFL≈0.47; use 6e-8 at 5000 cells); `Max time : 5e-4` (full case 3e-3); output every 5e-5 s.
- `Expansion location : 0.5` (velocity jump on the face between cells 2499 and 2500; see 4.1).
- `Cavitation pressure threshold : 0`, `Condensation temperature tolerance : 50` (user's choice for now: removes non-physical trace minority phase; to be revisited).
- `Artificial bulk viscosity coefficient : 5`. `Face slab width` removed (section 4.2).
- Pure-liquid start (user decision). Pelanti starts from α_vap=1e-2, whose slow mixture sound speed (her rarefaction head moves ≈130 m/s) makes it a different problem from pure water (c≈1.5 km/s, estimate); Fig. 4 comparisons are side by side, not a validation.
- Cost on 8 ranks to 5e-4 s: 2000 cells ≈21 s (≈7 ms/step), 5000 cells ≈135 s (≈16 ms/step). (The older ≈0.05 s/step figure included the debug probe.)

### 3.2 Failure 1 (deterministic, t≈0.104 ms, step 1727, cell 2501)
Chain in one relaxation call of a pure-vapor cell 0.92 K below Tsat (Tctol=0):
1. condensation seed 2.45e-3 kg/m³ (only the 5%·TG self-heating cap applied);
2. pT relaxation: the seed's latent heat warms the cell 17 K → strongly superheated;
3. linearized chemistry evaporates 98% of the seed, but VFeq=VF+S_α·dm leaves the rest at ρL=5004 kg/m³ (1/b=1513), TL=20.8 K, with ier=0;
4. Stage 3 floor (Tmin_liq=250 K) at that density raises p to 7.7e9 Pa and adds 1.66e10 J/m³ to the gas (= the `Floor dE` 0.132 J in `monitor/rescue`); the pressure step then produces ±3800 m/s faces.
Fixes 1–3 (section 1). Verified on the exact restart: step 1727 passes with Umax=2, PLmax=1e5, no floor.

### 3.3 Failure 2 (round-off sensitive, t≈0.21–0.35 ms; does not reproduce from restarts)
- Caught with the `BADEXIT` print and replayed offline (standalone driver in the session scratchpad, `replay/`, linking the thermo files with stubs).
- Entry: near-pure liquid (VF=1−1.2e-12) whose gas trace has **negative vapor mass** (Q8=−1.3e-16 vs Q2=1.7e-21); a tension-room-capped cavitation seed (2.45e-17); pT_relax returns VF=1.0 exactly with gas mass left → ρG=∞ → chemistry NaN, accepted by `finite_chem` (NaN passes all bounds tests), ier=0.
- Fixes 4–5. Verified: the replay returns finite values; the 0.28→0.5 ms stretch that died at 0.347 ms completes with no BADEXIT, floors or liquid rescues.
- Root cause not fixed: transport lets the vapor species mass leave [0, Q2] (same family as the "Yv overshoot" in the last commit). Fix would be in the flow solver (e.g. vapor flux = upwind Yv × gas mass flux).

### 3.4 Side effects seen
- With Tctol=0, cavity vapor just below Tsat re-seeds micro-droplets (VF≈1e-8) every call; their PL swings ±1e5 Pa → noisy PLmin/dPmax, harmless (mass ≈1.5e-5 kg/m³).
- With p_cav=1e30 (before the face fix) ≈700 two-phase cells formed across x≈0.23–0.81 by 0.21 ms.

## 4. cav_tube: flow-solver findings

### 4.1 Liquid sliver at the center (fixed by input)
**This case works much better when the velocity discontinuity sits exactly on a cell face. That is a workaround, not a fix: it points to an underlying numerical issue (a mid-cell jump should not freeze liquid or break the symmetry).**
- With `Expansion location : 0.5001` (mid-cell, from earlier this session) the cell holding x0 got the volume-weighted momentum u0·(1−2f_L) = 0: a full cell of water at rest on the symmetry plane. It evaporated while its faces opened (−19% by 0.05 ms), then froze at VF≈0.79 between two vapor cells at psat (no momentum, no superheat, no conduction).
- A single PLIC plane in that symmetric cell also put its vapor against one face, which broke the 2499↔2501 mirror symmetry from step 1.
- Fixed by `Expansion location : 0.5` (jump on a face; cavity now exactly mirror-symmetric). The volume-weighted momentum code in `cavtube_init` stays (inert for a face-aligned x0).

### 4.2 Velocity spikes from the half-cell face density (verified; open)
- Face velocity `u_F=(ρa·u_L+ρb·u_R)/(ρa+ρb)` and the pressure acceleration `−Δp/dx/((ρa+ρb)/2)` use the densities of the half cells next to the face (`subVF`).
- On the face between the interface cell and the cavity, the interface cell's half next to the face keeps some liquid for its whole draining (≈850 steps), so ρa (≈50 kg/m³ even at 5% liquid) ≫ ρb (vapor ≈0.3): the face moves with the liquid (≈−1.9 m/s) and the pressure step cannot change it. The cavity is emptied through that face while relaxation holds the interface cell at psat → a "dam": gap up to 1.9e4 Pa (cavity 2.87e4 vs 4.82e4).
- When that half drains, the weight flips to vapor and the stored gap is released at once (face +38 m/s, cell +58 m/s), then oscillates; repeats each time the interface enters a new cell (28–80 m/s spikes in the monitor).
- **Removed (2026-10-22), inconsistent:** `Face slab width` (cav_tube key → `amrvof%subVF_width`) = 1e-3: the face density is that of the phase touching the face, so the face sees vapor as soon as the interface enters the cell, before any gap builds. Momentum conservation is unchanged (each face's force is split ρ_slab/(ρa+ρb) between its two cells).
- Why removed: the half-cell density is the mass between the two cell centres that Δp/dx acts on (ghost-fluid density θρ1+(1−θ)ρ2). In a uniformly accelerated column it gives the exact face acceleration for any interface position; thin slabs are off by up to ≈1000× when vapor touches the face with liquid behind it (verified, `slab_check.py` in the session scratchpad). Hypothesis for the real cause: the one-velocity cell cannot hold the evaporation velocity jump m''(1/ρG−1/ρL) ≈ 2 m/s (m'' ≈ 0.58 kg/m²/s), which matches the 1.9 m/s dam face and ±1.9 m/s checkerboard.
- **Constraint (user, 2026-10-22):** the relaxation method exists to avoid a phase-change mass flux; no fix may compute m'' (rules out ghost-fluid velocity jumps from `dm`, Stefan/Hertz–Knudsen/CJ flux closures).
- **Test (0–0.3 ms, Tctol=50 in that test):** gap ≤ ±35 Pa (vs 19 kPa), no burst at step 2550 (face −0.02 vs +37.9 m/s), cavity at ≈4.80e4 Pa; brief 1 kPa gap at the next slab switch (step ≈3390), settling within a few hundred steps. Max Umax 7.5 vs 80.6 m/s, min PL 2.5e4 vs 9.0e3 Pa, TG 350–355 vs 335–397 K. Max dPmax 1.2e4 (vs 4.3e3), source not checked.

### 4.3 Vapor cell-velocity checkerboard (open)
- Without face memory (half cells, 5000 cells) the cavity cells alternate in sign at every output (e.g. 1.02, −1.52, 1.08, −1.28 m/s at 0.5 ms), with a hot centre (TG up to 370 K).
- With the slab fix, the first cavity cell settles at u≈+1.9 m/s (wrong sign) next to a face at ≈0; the mirror cell is −1.9, and so on. Face velocities (used by transport) are ≈0 and pressures balanced, so it no longer drives anything, but it sets Umax (2–7.5 m/s) and makes u look noisy.
- Why it persists: averaging cell velocities onto faces cancels an odd-even mode, and the cell momentum gets a central pressure gradient, which cannot see it.
- My suggestion (set pure-vapor cell velocities from the average of their face velocities) was rejected by the user as not right.

### 4.4 Earlier items (before this session's fixes)
- dPmax growth after ≈2.2 ms with Tctol=5 K (condensation droplets in the cavity, ±4 K vapor temperature scatter) went away with 50 K. Not rechecked with the face-slab fix and Tctol=0.
- Diffuse start (initial vapor fraction 1e-2 everywhere, Pelanti's setup) diverged from step 1. Reverted.

### 4.5 Face-velocity memory (tried 2026-10-22, removed by the user)
- Removed: an ad hoc change to a deliberately designed scheme whose consistency properties it was not shown to keep. Record of the test only.
- Diagnosis (derived): `get_face_velocity` rebuilds U = I(u) from cell velocities at every stage, so a face keeps its pressure impulse only through the two cell momenta. With uniform density one rebuild turns a face impulse F into (F₋+2F+F₊)/4: the compact Δp/dx cannot control the odd–even mode (4.3). Next to the interface cell the rebuild resets the cavity-side face to the liquid-weighted average every stage (hypothesis for the dam, 4.2). A mid-cell velocity jump is likewise averaged into one cell velocity and re-averaged at faces every stage (4.1).
- Method: increment form of time-step-consistent momentum interpolation (Choi 1999): U_f = U_f^n + I(u − u^n) + Δt·F_f, with I the existing half-cell density weighting, u the cell velocity after transport and relaxation, U^n = `Uold`, u^n = `Qold(5:7)/ρ_old`. Increment rather than carried face–cell difference, so a half-cell weight flip does not move the face by itself. Nothing changes U or ρu between steps (artificial viscosity only sets coefficients), so `Uold`/`Qold` are a consistent pair. Restart restores U.
- Boundary faces: a non-periodic domain-boundary face gets no pressure acceleration (extrapolated ghost, Δp=0), so with memory it froze at 2 m/s; when the rarefaction reached the outlets (≈3.2e-4 s) the outlet pressure fell to ≈600 Pa and a 1.6×/step instability followed. Fixed: those faces keep I(u) (their zero-gradient BC).
- Results to 5e-4 s (base = no memory; vapor mass per area, extremes over the run):

  | run | vapor mass | Umax | PLmin | PGmin | TG range | dPmax |
  |---|---|---|---|---|---|---|
  | base 2000 | 5.66e-4 | 58 | 3.1e3 | 2.5e4 | 288–564 | 49 |
  | memory 2000 | 5.92e-4 | 22 | 3.5e4 | 3.9e4 | 331–1529 (one event at step 109) | 28 |
  | base 5000 | 5.77e-4 | 85 | 3.1e3 | 2.5e4 | 288–396 | 101 |
  | memory 5000 | 5.86e-4 | 22 | 3.3e4 | 3.9e4 | 329–365 | 26 |

- Checkerboard gone at 5000 cells (cavity u smooth and monotone; TG 348–353 K, no hot centre); mirror symmetry 4e-6.
- Still there: the cavity pressure cycles with the cell-crossing period (≈1e-4 s): base dips to 2.9–4.0e4 Pa, memory to 4.4–4.6e4 (psat ≈4.72e4). With memory the cavity vapor sloshes coherently, up to −14.8 m/s next to the left interface at 0.5 ms (liquid −2 m/s), pressure higher at the interfaces than at the centre.
- At 2000 cells with memory the mirror symmetry is lost (|U+U_mirror| = 2.7 m/s; base 1.5e-2). Not investigated.
- Memory off reproduces the base run to round-off.

### 4.6 Grid study: vaporization is not mesh-dependent in cav_tube
- Runs to 5e-4 s without memory (1000/2000/4000/8000 cells; 8000 dropped from later studies as too many cells). Vapor mass per area at 5e-4 s: 5.92/5.66/5.76/5.79e-4 kg/m² (±2.5%, no trend); cavity width 1.965e-3 m on all grids; mean cavity vapor density 0.29–0.30 (saturation 0.29).
- Reason: the cavity width is set by kinematics (walls apart at 4 m/s) and relaxation fills it with vapor at ≈psat on demand; the latent heat (≈670 J/m² per side) comes from the liquid next to the interface.
- What does depend on the mesh: the dam noise (absent at 1000 cells; Umax 58–85 m/s and TG 288–564 K at 2000+) and, weakly, the interface liquid cooling (TLmin 352.38 → 352.10 K from 1000 to 8000 cells; cavity p 4.72e4 → 4.67e4). Intermediate times differ by up to 15% in vapor mass because of the dam cycles.
- The 0.160/0.108/0.085 trend in dodecane (section 6) is the supply-limited case: the evaporation front needs superheat it cannot reach without conduction.

## 5. Next steps (cav_tube)

1. Open, no accepted approach: the dam / cell-crossing pressure cycle (4.2) and the vapor checkerboard (4.3) with the original scheme.
3. Later: vapor-species boundedness in transport (3.3); remove the `BADEXIT` debug check once runs are clean; port fixes 1, 2, 4, 5 to `relax_igmix_sg_class`'s nucleation copy if that class is used; revisit the nucleation settings (p_cav=0, Tctol=50 suppress trace minority phase but should not be needed).

## 6. dodecane_shocktube (parked)

### Setup (current `input`)
2000 cells, liquid p=1e8 Pa, ρ=500 (642 K) left of x0, vapor 1e5 Pa, ρ=2 (1047 K) right; P22 Table 4 dodecane; dt=1e-7, t_end=473 µs. Level-set + volume-moment initial interface at `Interface location : 0.75025` (mid-cell). Bulk viscosity 5. `Relaxation type : p` currently. p_cav=1e30, Tctol=0. Keys `AMReX blocking factor` (default 8) and `Relaxation start time` (default 0) exist, both inert by default. Wall time on 8 ranks: p ≈25 s, pTg ≈100–150 s (cost ratio unexplained).

### What was established
- **p-relax validates the flow solver:** NGA2 matches Pelanti's p-relax and her exact 5-eq. solution (plateau 2.83e5 vs 2.81e5, u 240.3 vs 240.2, shock ≈0.93).
- **Bulk viscosity + mid-cell interface:** p-relax shock spike 297 → 251 m/s (exact 240), tail undershoot −5.8e5 → −3.5e5 Pa; pTg shock spike 322 → 270.
- **Tried and rejected:** AMR (no speed-up: 28 vs 25 s, no subcycling; liquid plateau wrong at level 0) and a 5 µs relaxation delay (no change in the final state; liquid stretched to −8.7e7 Pa during the delay).
- **Evaporation 3–4× too weak and grid-dependent** (0.160/0.108/0.085 kg/m² at 1000/2000/4000 cells; plateau 3.33e5 vs Pelanti 4.75e5). Derived: a sharp interface in equilibrium with no heat conduction has zero mass flux in the converged limit (energy jump m[(h_v+w_v²/2)−(h_l+w_l²/2)]=0 with h_v−h_l=L>0). Options: accept, or Pelanti's diffuse start (volumetric relaxation of the metastable liquid), or a relaxation layer of physical thickness at the interface. Mass-flux closures (interface conduction flux, Chapman–Jouguet front) are ruled out by the user: the relaxation method must not compute m''.
- **Vapor pinned on the grid inside moving liquid (verified):** a gas pocket between liquid neighbours is reconstructed against its upstream face (PLICnet normal from the stored barycenters, reset each step), so transport never moves it. The 3–5% vapor in the liquid plateau with p_cav=1e30 is a trail left by the moving interface. Fix would be in `amrvof` transport (needs a decision).
- **Rarefaction-tail undershoot** (≈0.9% of the pressure drop) also appears in amr_shocktube (water/air): flow-solver origin.

## 7. Practical notes

- Runs: `mpiexec -n 8 ./nga2.dp.gnu.opt.mpi.exe -i <input>` from the example folder. Test variants go in copies; to keep `amrviz/` untouched, run tests from a separate directory with the executable's full path. Foreground tool calls are limited to 10 min; no background runs without asking.
- Restart: `Restart from : <rtoa time>` (e.g. `1.00020E-04`; check `ls restart/<case>`). Restarts are not bit-identical; failure 2 moved between runs.
- To find a blow-up: restart from a checkpoint just before it with `Output period` ≈ dt, locate the cell, then probe its neighbourhood and use the debug box on the cell.
- Appending keys with `printf` needs a leading `\n` (inputs end without a newline). `param_read` reads the first token, so `0 ! 1e30` reads 0 (`!` is not a comment character in input files).
- `ny>1` must be a multiple of the blocking factor (8). For 1D AMR, ratio and blocking factor are 1 in single-cell directions.
- Python: `source ~/venv/combustion/bin/activate`; plot scripts use yt `covering_grid` at level 0.
- psat(T) for the dodecane pair: equal Gibbs energies g=(cp−q′)T − T(cp lnT − R ln(p+pinf)) + b·p + q; take the low root (bracket [1, 5e6] Pa).
- Session scratch (not in the repo): `evap_history.py`, `amr_profile.py`, `shocktube_exact.py`, `replay/` (offline relaxation replay), and the run folders `cav_blowup/`, `cav_center/`, `cav_seed/`, `cav_slab/`. 2026-10-22 scratch: `slab_check.py` (face-density consistency check), `grid/` (`base_nx2000`, `base_nx5000`, `mem_nx2000`, `mem_nx5000`, `off_nx2000`; `compare.py`, `load1d.py`).
- Shell: the tool shell is zsh, which does not word-split `$var` (use separate arguments or a function); macOS has no `timeout`. A foreground call over 10 min is moved to the background automatically, so keep each call to one run that fits.
- Conservation monitor: cross-section is dx² (1×1 cells in y,z), so per-area mass = monitor mass·nx²; gas volume per area = 1 − VFint·nx².
