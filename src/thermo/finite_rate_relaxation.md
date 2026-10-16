# Finite-rate relaxation (Pelanti) in NGA2 — working notes

Handoff for continuing the work in a new session. Status as of 2026-10-05: mechanical (`p_relax`) and
thermal (`pT_relax`) finite-rate relaxation are implemented and verified. Chemical (`pTg_relax`) is next.

## 1. Files

| File | Role |
|---|---|
| `src/thermo/relax_igmix_nasg_finite_class.f90` | Module `relax_igmix_nasg_finite_class`, type `relax_igmix_nasg_finite`. Registered in `src/thermo/Make.package`. |
| `examples/eos_relax/test_finite_prelax.f90` | Driver checking `p_relax` (see §6). |
| `examples/eos_relax/test_finite_pT.f90` | Driver checking `pT_relax` (see §6). |
| `examples/eos_relax/Makefile`, target `finite` | `make finite` builds and runs both drivers standalone (no AMReX/MPI). |
| `examples/eos_relax/check_KC_pelanti.py` | Numerical check that Pelanti's chemical stiffness equals the `step_C` stiffness of `relax_igmix_max_ent` (§7). Needs `source ~/venv/combustion/bin/activate`. |
| `~/Downloads/finite_rate_relaxation.{tex,pdf}` | Formal derivation of finite-rate relaxation within the entropy-maximization framework (§8). Draft subsubsection meant to follow "Matrix decomposition" in `~/Downloads/Multicomponent_evaporation_model-2/main.tex`. |
| `~/Downloads/relax_igmix_nasg_finite_class_draft.f90` | An earlier Claude draft of the same class (extends `relax_igmix_nasg`, has an `apply` override). Backup only; the repo file above is the user's version. |
| `~/Downloads/relax_nasg_ig_class_finite.f90` | The user's older finite-rate class (single-component gas). The repo class was ported from it. |

References: Pelanti 2022 (IJMF 153, 104097), two-phase, Sec. 7.2; Pelanti 2026 (Comput. Fluids 305,
106893), three-phase liquid–vapor–gas. Equation numbers below are Pelanti 2022 unless stated.

## 2. Class structure and conventions

- `relax_igmix_nasg_finite` extends `relax_igmix_sg` (it started as a copy of `relax_igmix_nasg`, so
  `initialize`, `get_p_eq`, `get_T_lvg`, `get_coeffs_lv` are duplicated from there). It overrides
  `p_relax` and `pT_relax`. `apply` and `pTg_relax` are inherited unchanged from `relax_igmix_sg`.
- Rates, all defaulting to `huge(1.0_WP)` = instantaneous:
  - `lambda` — mechanical, Pelanti's mu, `dalpha_l/dt = lambda*(p_l-p_g-Pjump)`, units 1/(Pa s)
  - `theta` — thermal, Pelanti's vartheta, heat to liquid `theta*(T_g-T_l)` per volume, units W/(m^3 K)
  - `nu` — chemical, Pelanti's nu, `dm_l/dt = nu*(g_v-g_l)`, units kg s/m^5 (not used yet)

  Because `Q` is per unit volume, Pelanti's values of mu, vartheta and nu apply directly.
- Pressure convention: `p` is the liquid pressure, and the gas sits at `p-Pjump` (Laplace jump).
- Energy convention, same as `relax_igmix_nasg`: `Q(3)+Q(4)` is conserved, and surface energy is not
  counted. `get_p_eq(VF,Q,qG,gammaG,Pjump)` returns the liquid pressure at which the phasic energies at
  fixed `VF` and fixed masses sum to `Q(3)+Q(4)`, with `p_l-p_g=Pjump`. This is the mixture energy
  relation (eqs. 59, 70).
- Gas composition is frozen during p and pT. The mixture-effective `cvG`, `cpG`, `qG`, `gammaG` come
  from `Yv=Q(iVQ)/Q(2)`, clamped to [0,1].
- Instantaneous test, used for both rates. Instantaneous means `rate=huge`, or `rate*K*dt >= 50` (so
  `exp(-50)` is below machine epsilon). The test is `log(rate)+log(K*dt) < log(50)`, so `rate*K*dt` is
  never formed. Otherwise the factor is `frac=1-exp(-rate*K*dt)`. `rate<=0` or `K*dt<=0` gives `frac=0`.
  NGA2's DEBUG build traps overflow on x86 (`-ffpe-trap=zero,overflow,invalid`), which is why the test
  is done in logs.

## 3. `p_relax` (Pelanti Sec. 7.2.1)

With phasic `p_k`, `rho_k`, `c_k`, Grüneisen `Gam_k`, and impedance-weighted interface pressure
`Pint=(Z_g p_l + Z_l p_g)/(Z_l+Z_g)`, `Z_k=rho_k c_k`:

- Inverse capacities (eq. 51): `1/xi_k = (Gam_k (Pint-p_k) + rho_k c_k^2)/alpha_k`. The same `Pint` is
  used for both phases, which is consistent with the shared-work energy convention.
- Exponential VF update (eqs. 53–56), with the Laplace jump in the force:
  `VFeq = VF - (p_g+Pjump-p_l)/(1/xi_l+1/xi_g) * frac`, where `K = 1/xi_l+1/xi_g`.
- The `VFratmax` clamp is applied before the energy update, so that energy is conserved at the final VF.
- Energy closure, which depends on the rate:
  - **Instantaneous:** `Peq=get_p_eq(VFeq,...)`, and the energies are rebuilt at `Peq` (liquid) and
    `Peq-Pjump` (gas). This is eq. 59.
  - **Finite:** the phases exchange only `Pint*dV` work: `Q(3)-=Pint*dVF`, `Q(4)+=Pint*dVF` (eqs. 1e–1f).
    The phases stay out of pressure equilibrium.

  The branch matters. Using `get_p_eq` for finite rates equalizes the pressures through the energy
  split even at `lambda=0`. That dumps liquid energy into the gas, which in a test showed up as roughly
  100 K of spurious gas heating.
- Failure codes follow `thermorelax_class` (`RELAX_FAILED`, `RELAX_BAD_LIQUID`, …). On failure the cell
  is left unmodified.

## 4. `pT_relax` (Pelanti Sec. 7.2.2)

- It first calls `this%p_relax`. The thermal ODE (eq. 60) is posed on the pressure-equilibrium
  manifold, so this mechanical step must be instantaneous. **User's decision:** the case sets
  `lambda=huge` (the default) whenever pT relaxation is used. `pT_relax` does not force it, and it
  does not propagate `p_relax`'s `ierr` (the original comment says "pT owns the final verdict").
- Post-mechanical state: `PL=get_p_eq(VF,Q,...)` and `PG=PL-Pjump`. Each phase is evaluated at its
  own pressure.
- `phi=(drho/dT)_p` and `zeta=(drho/dp)_T` (eqs. 16, 34). The EOS classes have no accessors for these,
  so they are coded explicitly:
  - NASG liquid: `phi=-rho(1-b rho)/T`, `zeta=rho(1-b rho)/(p+pinf)`. Past the `brhomax` knee, use
    `phi=-rho/T`, `zeta=rho/(p+pinf)`, matching the clamped branch of `nasg%get_rho_from_p_T`.
  - Ideal-gas mixture: `phi=-rho/T`, `zeta=rho/p_g`.
- Coefficients (eqs. 10, 61, 62):
  - `Z=(1-VF)Gam_l+VF Gam_g`
  - `D=VF rho_g c_g^2+(1-VF) rho_l c_l^2`
  - `COF=Gam_l rho_g c_g^2-Gam_g rho_l c_l^2`
  - `1/xiT_l=-(rho_l/VF Z+zeta_l COF)/(phi_l D)`
  - `1/xiT_g=-(rho_g/(1-VF) Z-zeta_g COF)/(phi_g D)`
- Update: `VFeq=VF+Z/D (T_g-T_l)/(1/xiT_l+1/xiT_g) * frac` (eq. 67), then the energies are rebuilt
  from `get_p_eq(VFeq)` (eq. 70). Unlike the mechanical case, the `get_p_eq` closure is exact here for
  any `theta`. On the pressure-equilibrium manifold the state is fixed by VF once the masses and total
  energy are given, so thermal disequilibrium is represented through VF alone.

## 5. Open items in the current code

1. **Inherited `apply` undoes finite `lambda`.** For `model=Prelax`, `relax_igmix_sg%apply` stage 2
   (around lines 338–361 of `relax_igmix_sg_class.f90`) re-splits the energy at fixed VF and forces
   `p_l-p_g=Pjump`. That erases finite-rate mechanical disequilibrium. It is harmless for
   `PTrelax`/`PThybrid`, because `pT_relax` already leaves the state in pressure equilibrium. A fix is
   an `apply` override that bypasses stage 2 for finite-`lambda` Prelax; one exists in the
   `~/Downloads` draft.
2. **Stale doc comments.** The `!>` headers above `p_relax` and `pT_relax` still describe the
   matm-style exact quadratic. The user should reword them (Claude was told never to remove the user's
   comments).
3. **Duplication.** The class could `extends(relax_igmix_nasg)` instead of copying its
   initialize/get_p_eq/get_T_lvg/get_coeffs_lv.
4. **`pTg_relax` doesn't yet respect finite rates.** It is inherited from `relax_igmix_sg`. It calls
   `this%pT_relax` (finite) and then the exact chemical solve (`solve_lv`/`solve_lvg`), which assumes
   the post-pT state is in thermal equilibrium. With finite `theta` that assumption is wrong, and `nu`
   is not used at all.

## 6. Verification (drivers in `examples/eos_relax`, `make finite`)

Both drivers use NASG water (Le Métayer–Saurel Table V) with a vapor+air ideal-gas mixture. Each runs
with `Pjump=0` and `Pjump≠0`. Build with
`FFLAGS="-O0 -g -ffree-line-length-none -Wall -ffpe-trap=zero,overflow,invalid"` to mirror the DEBUG
traps.

`test_finite_prelax`:
- Instantaneous (`lambda=huge` and `1e30`): `p_l-p_g-Pjump` at the round-off floor
  (eps·gamma·pinf ≈ 2e-7 Pa); mass and energy conserved to round-off.
- `lambda=0` leaves the cell bit-identical.
- Finite rate: the pressure gap after one step equals `exp(-Kp dt)` to 6 digits for `Kp dt` = 0.1, 1, 3.
  This checks the eq. 51 capacities together with the work closure.

`test_finite_pT`:
- Instantaneous `theta`: one step closes a thermal gap to O(gap²), and repeated steps converge. The
  result matches the exact quadratic of `relax_igmix_nasg%pT_relax` to about 1e-14 in VF.
- `theta=0`: no heat exchange.
- Finite `theta`: the one-step gap ratio follows `exp(-K_T dt)`. It also agrees to about 1e-5 with an
  independent reference that integrates Pelanti's actual source term: 20000 substeps of heat transfer
  `theta (T_g-T_l) h`, each followed by an instantaneous `p_relax`. This validates the eq. 62
  coefficients and the explicit `phi`/`zeta`.

Gotcha: Fortran is case-insensitive. `GL` (Grüneisen) collided with `gL` (liquid gamma) in `p_relax`;
the Grüneisen coefficients there are `GamL`/`GamG`.

## 7. Next step: finite-rate chemical relaxation (`pTg_relax`)

**Mapping, verified.** In Pelanti's chemical step (eqs. 71–83), mass transfer runs with mechanical and
thermal equilibrium held. Write `F_C = g_l - mu_v` (vapor at its partial pressure); `F_C>0` means
evaporation. Then:
- `dm_l/dt = -nu F_C`
- `dF_C/dt = K_C dm_l/dt`, where `K_C = Dv*Pi - Ds*Theta + K_CC` is exactly the stiffness in
  `relax_igmix_max_ent%step_C` (`src/thermo/relax_igmix_max_ent_class.f90`; coefficients from
  `get_state`).

Pelanti's `1/xi_1^g + 1/xi_2^g` equals the `Dv*Pi - Ds*Theta` part: her `S_p` and `S_T` are `Pi` and
`Theta`. `check_KC_pelanti.py` confirms this to 6–7 digits against a brute-force finite difference
along the pT-equilibrium manifold. `K_CC = R_v T (1-x_v)/m_v` is the mixing-entropy term. It is absent
in Pelanti 2022, where the gas is pure vapor.

With `K_C` frozen at the post-pT state `**`:
- `f = 1 - exp(-nu K_C dt)`
- `dm_l = -f F_C / K_C` (Pelanti eq. 77, with `Delta g = -F_C`)

As `nu -> infinity` this is one linearized `step_C`.

**Useful identity.** `F_C = +R_v T pTsat(p_l, p_v, T)`, where `pTsat` is the existing saturation
residual in `relax_igmix_sg`. It is only valid at thermal equilibrium. Note the sign.

**Proposed implementation.** Replace only the block in `relax_igmix_sg%pTg_relax` around lines 940–999,
"Solve chemical equilibrium … Update Q". Nucleation, `activate_chem`, the pure-phase snap and
`undo_nucleation` stay as they are.
1. At `**`, compute `F_C` and `K_C`. Factor the `K_C` computation out of `relax_igmix_max_ent` into a
   shared helper rather than copying it.
2. Move the mass in the log variable, as `step_C` does: `du=f*F_C/(K_C*m_v)`. Cap evaporation at the
   available liquid. Then `dm_l=-m_v(exp(du)-1)`. This keeps `m_v>0` automatically. Route the
   complete-evaporation cap to the existing pure-phase snap, not to Pelanti's epsilon clamp (eq. 83).
3. Update `Q(1)`, `Q(2)` and `Q(iVQ)`, and move energy with the donor phase's specific energy (this only
   seeds the split). Then call the exact pT closure at the new masses. That closure depends only on
   masses, composition and total energy, so VF, p and T land exactly on the pT-equilibrium manifold.
   This is stricter than Pelanti, who takes alpha from the linearization and leaves a small T_l≠T_g.
4. Optional refinement ("target-force" closure): solve the 1-D root `F_C(m_l) = F_C** exp(-nu K_C dt)`.
   Each evaluation is one pT closure plus one `pTsat`, and the derivative is `K_C`. `F_C` is monotone
   in `m_l`, so a bracketed Newton solve is safe. This reproduces the current exact equilibrium as
   `nu -> infinity`, whereas the linear step differs from it at O(F_C²).

**Pitfalls**
- Pelanti's eq. (72) (and eq. 39 in 2026), read literally, gives `s_1+s_2` in the sum. The correct
  form is `(-1)^(k-1) (S_p/rho_k - s_k S_T)`, which gives `s_1-s_2`. Read literally, the rate depends
  on the entropy reference constants.
- The 2026 three-phase paper does not apply directly. There, air is a separate phase with its own
  volume fraction (and optionally its own temperature), and vapor is evaluated at the total pressure.
  In NGA2, vapor and air share one gas phase (Dalton): `p_v` → `p_sat`, not `p`. Its `S_alpha` formulas
  (7–8) do not carry over, and its Fig. 1 should not be expected to match. The right template is
  Pelanti 2022 plus `K_CC`.
- Validation case: Pelanti 2022 §8.2.1, the water cavitation tube. NASG Table 3, `alpha_vap=1e-2`,
  `p=1 bar`, `T=353 K`, `u=±2 m/s`, 5000 cells, `nu=1e-4`, `t=3 ms`, chemical relaxation only when
  `T>T_sat`. With pure vapor, `K_CC=0`, so Figs. 4–5 should be reproducible.

## 8. Formal framework (summary of the PDF)

Entropy maximization fixes only the end state, `F=0`. The rate comes from the parent model's relaxation
laws.
- Pelanti's sources are exactly linear in the entropy forces `F=[F_V,F_E,F_m]`, with state-dependent
  coefficients. Linearized at equilibrium they give `xi_dot = L F` with `L = D Lambda D^T`.
- `D = [d_M,d_T,d_C]` with `d_M=[1,-p,0]`, `d_T=[0,1,0]`, `d_C=[0,g,1]`.
- `Lambda = diag(mu T, vartheta T^2, nu T)`.
- `L` is symmetric (Onsager) and diagonal in `D`, with no cross-coupling. The entropy production is
  `sum lambda_i f_i^2`, where `f=D^T F`.
- Freezing `K~=D^T K D` gives `X = K~^{-1}(I - exp(-K~ Lambda dt)) f`, the matrix generalization of
  Pelanti's scalar exponentials. Its `Lambda -> infinity` limit is the Newton step `J D X = -F`.
- Under rate separation `lambda_M >> lambda_T >> lambda_C`, the eigenvalues of `K~ Lambda` tend to
  `lambda_i s_i`, where `s_i` are the M→T→C elimination pivots (verified numerically). That recovers
  the sequential scheme with one scalar exponential per channel, including Pelanti's `K_p`, `K_T`, `K_g`.
- `d_C`: the matrix-decomposition section of `main.tex` uses `[v_l,e_l,1]`, the parent model `[0,g,1]`.
  This doesn't matter under instantaneous M and T, because the chemical pivot is invariant. It does
  matter for a fully coupled finite-rate solve.

---

## 9. Session 2026-10-06: `pTg_relax`, surface energy, check against Pelanti

Status at the end of this session: all three channels (`p_relax`, `pT_relax`, `pTg_relax`) are
implemented at finite rate in `relax_igmix_nasg_finite`. They conserve `E_l+E_g+Pjump*V_l`, and they
reproduce a literal transcription of Pelanti 2022 to round-off in her setting. Next: the 1D cavitation
tube (§9.6). Parts of §§5 and 7 are superseded; see §9.5.

### 9.1 Decisions made by the user

- **Model.** Pelanti 2022 eqs. 1a–1f, with three changes: `g_2 -> g_v(p_v,T_g)`, i.e. vapor at its
  partial pressure; an extra vapor equation `d((1-alpha)rho_g Y_v)/dt = nu (g_l - g_v)`; and the Laplace
  jump `Pjump`. Pelanti 2024/2026 add nothing usable, because their air is a separate phase and their
  vapor is at total pressure.
- **Energy convention: conserve `E_l+E_g+surface energy` = `Q(3)+Q(4)+Pjump*VF`.** This replaces the
  old `Q(3)+Q(4)` convention in the finite class. The user insisted on it after a version that kept
  `Q(3)+Q(4)`.
- **Follow Pelanti literally.** No substeps, no extra thermal step after the chemical closure, no
  ln(m_v) update. The user rejected a substepped version as "not Pelanti": the code must map
  equation by equation onto her paper.
- **Rate names:** `lambda` (Pelanti mu), `theta` (vartheta), `nu`.
- **Clustering.** When a phase runs out, the user's clustering framework can take over instead of
  Pelanti's eq. 83 clamp. Not wired in: the code currently uses eq. 83, see §9.3.

### 9.2 Formulation (`~/Downloads/finite_rate_relaxation.tex`, revision 4; rev 2 is backed up in the old session scratchpad)

The tex follows `main.tex` conventions: liquid minus gas, `F_M`, `F_T`, `F_C`, and the states
0, *, **, ⊗. It goes after "Sequential clustering". It is laid out like Pelanti eqs. 50–83 (ODEs,
coefficients, frozen solution, closure) and ends with an M1–C5 procedure list. Merge notes: drop the
recap paragraph and the three placeholder labels (`dp_M_Ks`, `dH_T`, `sys_T_1`), add `mathrsfs` and
`booktabs`, and add `Pelanti2022` and `Pelanti2026` to `ref.bib`.

- **Mechanical.** `1/W_k^eff` is Pelanti's `1/xi_k^p` (eq. 51), with the gas working pressure
  `p_I-Pjump`. The interface pressure is `p_I=(Z_g p_l+Z_l(p_g+Pjump))/(Z_l+Z_g)`.
- **Thermal.** Pelanti's eq. 62, with each phase at its own pressure, stays exact under
  surface-energy conservation for `Pjump≠0` (verified). It is equivalent to
  `K_T=D/N`, with `N=C_pl C_vg W_Tg + C_pg C_vl W_Tl` and `S_V^T=Z/D`.
- **Chemical (re-derived eq. 71 coefficients).** `[S_alpha,S_T,S_p]` solve Pelanti's eq. 97 matrix,
  with two changes:
  - `+Pjump` in the alpha_1 entry of the energy row;
  - the right-hand side `[1, -R_v/R_m, b3]`, with
    `b3=(1-Y_v)[-(R_v-R_a)e_g/R_m+(cv_v-cv_a)T+q_v-q_a]`. It reduces to her `[1,-1,0]` at `Y_v=1`.

  Also `1/xi_2^g=-(S_p R_v T/p_g - s_v(p_v,T) S_T)+K_CC` with `K_CC=R_v T(1-x_v)/m_v`, and
  `1/xi_1^g=S_p/rho_l-s_l S_T` (eq. 72, corrected). Plugging mixture properties into her eq. 23 is
  wrong, by up to 15× in `K_g`.
- **Chemical step only at thermal equilibrium.** When `T_l≠T_g`, `F_C` depends on the entropy datum.

### 9.3 Code (`src/thermo/relax_igmix_nasg_finite_class.f90`)

- **`p_relax`:**
  - `Pint` is shifted by `Pjump` on the gas side (the old line is kept, commented out).
  - `xiGinv` uses `Pint-Pjump-PG`.
  - The instantaneous closure is `get_p_eq(VFeq,[Q(1:3),Q(4)-Pjump*(VFeq-VF)],...)`.
  - The finite-rate gas energy is `Q(4)+dW-Pjump*(VFeq-VF)`.
- **`pT_relax`:** same closure target.
- **`pTg_relax`** is a verbatim copy of `relax_igmix_sg%pTg_relax`, with these changes:
  - **Surface-energy terms:** `-Pjump*dVF` wherever VF changes. That is the pure-gas and pure-liquid
    snaps (including the energies their stability tests use) and both nucleation seeds.
  - **Finite-rate branch:** the `finite_chemical` block goes before the exact solve when `nu<huge`.
    It calls the internal `finite_chem`, which maps onto Pelanti like this:
    - `chiL/G=c^2-Gamma*h`, `PHIL/G`, `zetaL/G` → eqs. 16, 34;
    - `CEp`, `CET` → eq. 101;
    - `M`, `rhs`, Cramer via `det3` → eqs. 97, 100, 23;
    - `xi1inv`, `xi2inv` → eq. 72;
    - `frac`, `dm`, `VFeq` → eqs. 76–78;
    - `abar` clamp → eq. 83;
    - `get_p_eq` → eq. 82.
  - **Fallbacks to the exact solve:** `nu=huge`, or vacuum vapor (`m_v<=0` or `pv<pv_min`). In the
    second case a dry cell still evaporates instantaneously.
  - **Rejected step:** `m_l<=0` or `m_v<=0` after the update, which is possible with air.
  - **Exact path:** the `surface_energy` loop iterates `rhoe0=Q0(3)+Q0(4)-Pjump*(VF-VF0)` (at most 5
    passes, 1 pass when `Pjump=0`).
  - **`p_v`** uses `get_xv` (molar masses), so `F_C=0` coincides with the exact solver's equilibrium.
  - **Assumes** `lambda=theta=huge` whenever `nu` is finite. Not enforced.
  - **Imports:** `Mv,Ma,dbg_cell,dbg_i,dbg_j` from `relax_igmix_sg_class`; `RELAX_NUC_FAILED`,
    `RELAX_VACUUM_VAPOR`.
- **Known limitation of the single step.** It is first order in the variation of the coefficients
  over `dt`, and with air `K_CC∝1/m_v`. Far from saturation (`|Delta g|≳R_vT`) one step is poor: at
  `nu*K*dt=10` evaporation removes 52% of the force where the reference removes 89%, and condensation
  overshoots. This is fine when relaxation runs every flow step (`|Delta g|≪R_vT`).

### 9.4 Verification (`examples/eos_relax`, `make finite` builds and runs all four)

Build with `FFLAGS="-O0 -g -ffree-line-length-none -fcheck=all -ffpe-trap=zero,overflow,invalid"`.

- `test_finite_prelax`, `test_finite_pT`: unchanged results. Their `dE` now measures
  `Q3+Q4+Pj*VF`, and it is at round-off.
- `test_finite_pTg` (new): evaporating and condensing cells, `Pjump` = 0 and 300.
  - `nu=0` leaves the cell unchanged.
  - `nu=huge` vs `relax_igmix_nasg`: ~1e-12 at `Pjump=0`; at `Pjump=300` it differs by ~1e-9,
    which is expected because that class conserves `Q3+Q4`.
  - Finite `nu` vs a fine integration of `m_l'=-nu F_C`: ~3e-4 at `nu*K*dt=0.1`.
  - Mass, air and `E_t` conserved to round-off.
- `test_pelanti_ref` (new) is the literal Pelanti 2022 transcription, with its own EOS from her eqs.
  28–34, her eq. 23 in closed form, and her Table 3 water with pure vapor at `Pjump=0`. It covers
  eqs. 51, 56–59, 10, 60–70, 23, 72 and 76–82. **Max relative difference 1.5e-11: PASS.** The 1e-11
  is round-off amplified in `Delta g`. Not covered: eq. 83, and the energy update for finite `mu`
  (which is not in her paper).
- `check_finite_rate_map.py` checks every mapping identity against finite differences at 40 digits.
  `check_finite_rate_steps.py` tests the step algorithm; it also holds this session's
  substep/ln(m_v) experiments, and it is slow (~15 min).

### 9.5 Status of earlier items

- §5.1 (inherited `apply` stage 2 undoes finite `lambda` for `Prelax`): **done** (2026-10-06). The
  finite class overrides `apply`: finite-`lambda` Prelax calls `p_relax` and skips stages 2–3 (and the
  dissolution snaps); everything else calls the parent's `apply`. The parent's `apply` is now `public`
  in `relax_igmix_sg_class` and is imported as `sg_apply`, so it can be called as `sg_apply(this,...)`.
  `this%relax_igmix_sg%apply` would lose dispatch, so `p_relax`/`get_p_eq` would resolve to the SG
  versions. `test_finite_prelax` covers the `apply` path. `PThybrid` with finite `lambda` is still
  undone by stage 2 (outside the `lambda=huge` rule).
- §5.2 (stale `!>` headers): still open. Also the module header line "Inherits pTg_relax unchanged
  from parent" is now stale. These are the user's comments: the user rewords them, Claude must not
  delete them.
- §5.3 (duplication of `relax_igmix_nasg`): still open.
- §5.4 (finite `pTg_relax`): **done**.
- §7 "Proposed implementation": **superseded.** No ln(m_v) update, no target-force variant, no
  substeps; Pelanti's eq. 82 closure and eq. 83 clamp.
- The inherited `apply`'s dissolution snaps (`diss_P`, `diss_RHO`) still conserve only `Q3+Q4`. They
  are off by default.

### 9.6 Next: 1D water cavitation tube (Pelanti 2022 §8.2.1)

- **Setup** (from `src/thermo/pelanti_papers_summary.md` §2.8):
  - domain [0,1] m;
  - `alpha_vap=1e-2`, `p=1e5` Pa, `T=353` K (below `T_sat`=372.3 K);
  - `u=-2` m/s for x<0.5 and `+2` m/s for x>0.5;
  - NASG **Table 3 of P22**: liquid γ=1.187, ϖ=7028e5, η=−1177788, η′=0, κv=3610, b=6.61e-4;
    vapor γ=1.467, η=2077616, η′=14317, κv=955. These are not the Le Métayer values in `input_NASG`
    (γ=1.19);
  - 5000 cells, CFL 0.5, t=3 ms.
- **Five runs:** p; pT(f) with `theta=2000` W/(m³K); pT; pTg(f) with `nu=1e-4` kg s/m⁵; pTg.
  Compare her Figs. 4–5 (p, T, α_vap, Y_vap, u), digitized e.g. with WebPlotDigitizer.
- **Pure vapor:** air species present with zero mass, `Q(8)=Q(2)`, so `K_CC=0` and the model is
  exactly P22. `sigma=0`, hence `Pjump=0`.
- **She relaxes chemically only where `T>T_sat` (evaporation only).** `activate_chem` also condenses,
  so a switch is needed.
- **Trap: `pTg_relax` pure-phase thresholds.** Liquid VF=0.99 is above `pure_VFhi=0.95`, so the
  pure-liquid snap fires and folds the gas into the liquid when `pL>=psat`, which is true at 1 bar
  and 353 K. VF=0.99 also sits at the `near_pure` boundary (`VF_nuc=0.01`), so nucleation logic may
  run. Both must be disabled or re-thresholded for this test, or the initial condition collapses.
  The same applies to the solver's pure-cell snap with `VFlo`/`VFhi`.
- **No case source yet:** `examples/cavitation1D` holds only build output and no source. Build it from
  `examples/cavitation` (`src/simulation.f90` already uses the relax classes). Relaxation must run
  every time step, with `lambda=huge` for all cases.
- **Expect differences:** Pelanti uses a 2nd-order wave-propagation HLLC solver, while NGA2 uses
  `amrmpcomp`, so profiles should agree up to discretization. The pT cases can also be checked against
  the exact pT-relaxed solution she mentions.

### 9.7 Session 2026-10-06 (later): no iterative solve, dry gas, amrcomp_impact

Supersedes the parts of §9.3 about the exact-solve fallback and the eq. 83 clamp.

- **`pTg_relax` has no Newton solve.** The chemical step is always `finite_chem`; `nu=huge` gives
  `frac=1`. The copied exact path (`solve_lv`/`solve_lvg`, boundary routines, `check_cons`), the Raoult
  starting guess and `activate_chem`'s dry-edge `Yv` reseeding are removed. The Raoult guess was a bug:
  it set `Yv` to saturation before `activate_chem`, so at `Pjump=0` liquid cells with `Yv<Y_small` never
  evaporated. **The same bug is still in `relax_igmix_sg%pTg_relax`** (untouched).
- **Pure-phase snap** skipped for finite `nu` (kept for `nu=huge`). **Nucleation** only for `nu>0`.
- **No eq. 83 clamp.** A phase that runs out (mass or linearized volume) transfers completely: all
  liquid evaporates (VF=0, `E_g=E_t`); all vapor condenses (air-free gas, `Q8>=(1-Y_small)Q2`: VF=1,
  `E_l=E_t-Pjump`; with air: `dm=Q8`, VF from eq. 78).
- **Dry gas** (`pv<pv_min`): eq. 77's limit is `dm=0`, but the ODE evaporates. Only the singular term
  `F_C=Rv T ln(m_s/m_v)` is kept and integrated exactly: `E1(ln(1/u))=E1(ln(1/u0))+nu Rv T dt/m_s`,
  `u=m_v/m_s`, `m_s=(Mv/Ma) m_a p_vs/p_g`, `p_vs=exp((g_l-g_v0)/(Rv T))`. `log_E1`/`inv_log_E1` are
  module procedures. Derivation: tex revision 5 (`~/Downloads/finite_rate_relaxation.tex`, rev 4 backed
  up as `_rev4.tex`). `test_finite_dry` (in `make finite`): matches a 1%-substep reference to 6 digits
  for `tau<=1e-5`, 3e-3 at `tau=1`. Gotcha hit again: local `mA` shadowed module `Ma`.
- **Guard** in `apply`: `die` unless `lambda=huge` for pT/pTg/pThybrid and `theta=huge` for pTg.
- **`test_finite_pTg`**: the `nu=huge` vs `relax_igmix_nasg` line now compares one linearized step with
  the exact equilibrium, so it differs (both test cells are far from saturation; case 2 overshoots,
  `F_C/F0=-1.10`). The label should be updated.
- **`amrcomp_impact`** uses `relax_igmix_nasg_finite`; input keys `Mechanical/Thermal/Chemical
  relaxation rate` (default instantaneous). Builds; not run yet.
- **Still open:** the `Y_small` reducing-phase check after the chemical step (blocks evaporation when
  liquid mass fraction <1e-3 and condensation when `Yv0<1e-3`); `dissolve_stranded_gas` is an unrated
  phase change; clustering Need is rate-unaware for finite `theta`; stale user comments at the
  `pTg_relax` header and the "solve_lv/solve_lvg" comment.

### 9.8 Session 2026-10-07: chemical and depletion clustering (amrmpcomp passes A3, A4)

- Observation (impact, t=0.01, `nu=1e-3`): unclustered VF≈1 cells reach `Yv` 5-25x their neighbours'.
  Cause: Pelanti's transfer is per unit cell volume, so `dYv ∝ 1/m_g`; a sliver of gas saturates within
  a few calls (replayed: `finite_chem` lines for `dm` and `Q(iVQ)`).
- `thermorelax%get_chem_transfer(dt,VF,Q,Pjump,dm,VFeq,Ysat,ok)`: default `ok=.false.`; the finite class
  implements it with the former coefficient block of `finite_chem` (eqs. 71-78 and the dry step), which
  now calls it. Refactor checked bit-for-bit on all five `make finite` drivers.
- `amrmpcomp`: `clust_chem` (A3) need `|dYv|/Ysat - cluster_dYtol` (gas side), `clust_depl` (A4) need
  `removed fraction of a phase - 1` (liquid on evaporation; gas on condensation only if air-free).
  Prediction at `Pjump=0`. Run after A1/A2 with the same VF/Q/cskip refresh. Flags `cluster_chem_on`,
  `cluster_depl_on` (default `.false.`); impact input keys `Chemical clustering`, `Depletion clustering`.
- Offline replay on plotfile t=0.01 (434 mixture cells, partner guards as in `score`): A3 28 claimants
  (all VF>=0.9857), A4 8 claimants (liquid slivers), 0 stranded; pooled needs below zero.

### 9.9 Session 2026-10-08: instantaneous closed forms, pass order, A3 merged

- Overshoot at (59,978,0), t=12.2: stage 2 A1 (then thermal) gave the only gas-rich neighbour (58,978)
  to (58,977); A2 (then mechanical, need 362) stranded the cell. Alone, Pelanti's linear `p_relax`
  failed (VFeq>=1) and `pT_relax` continued from the unrelaxed state: `get_p_eq` dumped liquid energy
  into the gas (`TG=1.22`). Even paired, eq. 57 + eq. 59 heated the gas to `TG=16.9` (residual liquid
  overpressure x `(V_l-b m_l)/(gamma_l-1)` goes to the gas).
- **Instantaneous rates now use closed forms**: `p_exact` (matm quadratic of `relax_igmix_nasg%p_relax`)
  and `pT_exact` (its p,T quadratic), modified to conserve `Q3+Q4+Pjump*VF` (gas works at `Pw-Pjump`;
  thermal target `Eth+Pjump(VF-b Q1)`). Used when `lambda`/`theta` = huge or rate*Kdt>=50; bitwise
  equal to `relax_igmix_nasg` at `Pjump=0`. The old eq. 59 `if (inst)` branch in `p_relax` is now
  unreachable. `test_pelanti_ref` now differs from Pelanti in its two infinite-rate lines (expected).
- With the closed form the solo cell is thermally correct (`TG=TL`, `RHOG=268`) but `PG=13.6` is its
  own mechanical equilibrium (liquid at 13.6, gas sliver 1.5e-4): it must be clustered.
- **Pass order = stage order**: A1 mechanical, A2 thermal, A3 chemical (`clust_mech=1`,
  `clust_thermal=2`, `clust_chem=3`). **A3 merges composition and depletion**: need =
  max(`|dYv|/Ysat - cluster_dYtol`, removed fraction - 1), side of the larger; one flag
  `Chemical clustering` (`Depletion clustering` and `cluster_depl_on` removed). Replay at t=0.01:
  36 claimants (28 gas, 8 liquid), 0 stranded.
