# Pelanti's arbitrary-rate relaxation papers — implementation summary

A self-contained digest of three papers by Marica Pelanti, meant to be read in place of the PDFs.
It keeps every equation, closure, parameter table and test setup needed to implement or validate the
methods. Equation numbers follow each paper.

- **[P22]** M. Pelanti, *Arbitrary-rate relaxation techniques for the numerical modeling of compressible
  two-phase flows with heat and mass transfer*, Int. J. Multiphase Flow 153 (2022) 104097.
  Two-phase, six-equation model; mechanical, thermal and chemical relaxation, all at arbitrary rate.
- **[P24]** M. Pelanti, *Numerical Modeling of Liquid-Vapor-Gas Flows with Arbitrary-Rate Mass Transfer*,
  in Hyperbolic Problems: Theory, Numerics, Applications Vol. II, SEMA SIMAI 35 (2024) 185–194.
  Short conference version of [P26].
- **[P26]** M. Pelanti, *Numerical relaxation techniques for mass transfer in three-phase liquid-vapor-gas
  flows*, Computers and Fluids 305 (2026) 106893. Three-phase model; arbitrary-rate chemical relaxation
  only, with mechanical and thermal relaxation instantaneous.

Section 8 is an errata list. Several printed equations are wrong; the corrected forms were checked
numerically (scripts in `examples/eos_relax/`).

---

## 1. Notation (all three papers)

Phase index: 1 = liquid, 2 = vapor, 3 = inert gas (three-phase papers only).

| Symbol | Meaning |
|---|---|
| $\alpha_k,\ \rho_k,\ p_k,\ T_k$ | volume fraction, density, pressure, temperature |
| $\mathcal{E}_k=\rho_k\varepsilon_k$ | internal energy per unit volume ($\varepsilon_k$ specific) |
| $E_k=\mathcal{E}_k+\rho_k|\vec u|^2/2$ | total energy per unit volume |
| $\rho=\sum\alpha_k\rho_k$, $Y_k=\alpha_k\rho_k/\rho$ | mixture density, mass fraction |
| $\mathcal{E}=\sum\alpha_k\mathcal{E}_k$, $E=\sum\alpha_kE_k$ | mixture internal and total energy per volume |
| $h_k=(\mathcal{E}_k+p_k)/\rho_k$, $s_k$, $g_k=h_k-T_ks_k$ | specific enthalpy, entropy, chemical potential |
| $c_k$ | phasic sound speed, $c_k^2=\Gamma_kh_k+\chi_k$ |
| $\Gamma_k=(\partial p_k/\partial\mathcal{E}_k)_{\rho_k}$ | Grüneisen coefficient |
| $\chi_k=(\partial p_k/\partial\rho_k)_{\mathcal{E}_k}$ | |
| $\phi_k=(\partial\rho_k/\partial T_k)_{p_k}=-\rho_k\beta_k$ | $\beta_k$: thermal expansion coefficient |
| $\zeta_k=(\partial\rho_k/\partial p_k)_{T_k}=\rho_k\mathcal{K}_{Tk}$ | $\mathcal{K}_{Tk}$: isothermal compressibility |
| $\kappa_{pk},\ \kappa_{vk}$ | specific heats at constant pressure and volume |
| $C_{pk}=\alpha_k\rho_k\kappa_{pk}$ | extensive heat capacity (per unit mixture volume) |
| $c_f=\sqrt{\sum Y_kc_k^2}$ | frozen (non-equilibrium) sound speed |

Useful relations [P22 eq. 17]:
- $\beta_k=\Gamma_k\kappa_{pk}/c_k^2=\Gamma_kC_{pk}/(c_k^2\alpha_k\rho_k)$
- $\mathcal{K}_{Tk}=1/(\rho_kc_k^2)+\beta_k^2T_k/(\rho_k\kappa_{pk})$

Mapping to NGA2 EOS parameters: $\varpi\to$ `pinf`, $\eta\to$ `q`, $\tilde\eta\to$ `qp`, $\kappa_v\to$ `cv`,
$\gamma\to$ `gamma`, $b\to$ `b`.

---

## 2. [P22] Two-phase model

### 2.1 Equations (eq. 1)
$$
\begin{aligned}
&\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mathcal{P},\\
&\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\mathcal{M},\qquad
\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\mathcal{M},\\
&\partial_t(\rho\vec u)+\nabla\cdot(\rho\vec u\otimes\vec u+(\alpha_1p_1+\alpha_2p_2)\mathbb I)=0,\\
&\partial_t(\alpha_1E_1)+\nabla\cdot(\alpha_1(E_1+p_1)\vec u)+\Sigma=-p_I\mathcal{P}+\mathcal{Q}+(g_I+|\vec u|^2/2)\mathcal{M},\\
&\partial_t(\alpha_2E_2)+\nabla\cdot(\alpha_2(E_2+p_2)\vec u)-\Sigma=p_I\mathcal{P}-\mathcal{Q}-(g_I+|\vec u|^2/2)\mathcal{M},
\end{aligned}
$$
with $\Sigma=-\vec u\cdot\vec\Xi$, $\vec\Xi=Y_2\nabla(\alpha_1p_1)-Y_1\nabla(\alpha_2p_2)$.

Relaxation terms (eq. 2):
$$\mathcal{P}=\mu(p_1-p_2),\qquad \mathcal{Q}=\vartheta(T_2-T_1),\qquad \mathcal{M}=\nu(g_2-g_1),\qquad \mu,\vartheta,\nu\ge0.$$
The rates may be functions of the state; in that case they are frozen at the start of the relaxation
step.

Interface closures. The second law holds if $p_I\in[\min p_k,\max p_k]$ and
$g_I\in[\min g_k,\max g_k]$ (Flåtten & Lund 2011). The tests use
$p_I=(Z_2p_1+Z_1p_2)/(Z_1+Z_2)$ with $Z_k=\rho_kc_k$. $g_I$ is never needed numerically.

Units of the rates (Q per unit volume):
- $\mu$: 1/(Pa s).
- $\vartheta$: W/(m³K). [P22] prints "Pa/(s K)", which is the same.
- $\nu$: kg s/m⁵. [P22] prints "Pa kg²/(s J²)", which is the same.

Mixture total energy is conserved: $\partial_tE+\nabla\cdot((E+\alpha_1p_1+\alpha_2p_2)\vec u)=0$ (eq. 7).

Phasic pressure equations (eq. 8):
$$\partial_tp_k+\vec u\cdot\nabla p_k+\rho_kc_k^2\nabla\cdot\vec u=\frac{(-1)^{k-1}}{\alpha_k}\Big\{-[\Gamma_k(\mathcal{E}_k+p_I)+\chi_k\rho_k]\mathcal{P}+\Gamma_k\mathcal{Q}+(\Gamma_kg_I+\chi_k)\mathcal{M}\Big\}.$$

[P22] drops the $\mathcal{M}/\rho_I$ term in the $\alpha_1$ equation that Pelanti & Shyue 2014 had. The
rigorous reduction (Appendix A) shows that mass transfer affects $\alpha$ only through pressure relaxation.

### 2.2 Relaxed-model hierarchy
Subcharacteristic ordering: $c_{pTg}\le c_{pT}\le c_p\le c_f$ (eq. 26).

**p-relaxed (Kapila five-equation) model** (eq. 9), the limit $\mu\to\infty$:
$$\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1-\frac{\alpha_1\alpha_2}{D}(\rho_2c_2^2-\rho_1c_1^2)\nabla\cdot\vec u=\frac1D(\alpha_2\Gamma_1+\alpha_1\Gamma_2)\mathcal Q+\frac1D(\alpha_2\omega_{I1}+\alpha_1\omega_{I2})\mathcal M,$$
- $D=\alpha_1\rho_2c_2^2+\alpha_2\rho_1c_1^2$ (eq. 10)
- $\omega_{Ik}=\Gamma_k(g_I-h_k)+c_k^2=\Gamma_kg_I+\chi_k$ (eq. 11)
- Mixture energy relation: $\mathcal E=\alpha_1\mathcal E_1(p,\rho_1)+\alpha_2\mathcal E_2(p,\rho_2)$ (eq. 12)
- Wood sound speed: $c_p^{-2}=\rho\sum\alpha_k/(\rho_kc_k^2)$ (eq. 13), equivalently $1/D=\rho c_p^2/(\rho_1c_1^2\rho_2c_2^2)$

Pressure (eq. 14):
$$\partial_tp+\vec u\cdot\nabla p+\rho c_p^2\nabla\cdot\vec u=\frac1D\big[(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\mathcal Q+(\rho_2c_2^2\omega_{I1}-\rho_1c_1^2\omega_{I2})\mathcal M\big].$$

Phasic temperatures (eq. 15):
$$\partial_tT_k+\vec u\cdot\nabla T_k+\frac{\rho c_p^2}{\phi_k}\Big(-\zeta_k+\frac1{c_k^2}\Big)\nabla\cdot\vec u=\frac{1}{\phi_kD}\Big[(-1)^k\frac{\rho_k}{\alpha_k}(\alpha_2\Gamma_1+\alpha_1\Gamma_2)-\zeta_k(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big]\mathcal Q$$
$$\qquad+\frac1{\phi_k}\Big\{\frac{(-1)^{k-1}}{\alpha_k}+\frac1D\Big[(-1)^k\frac{\rho_k}{\alpha_k}(\alpha_2\omega_{I1}+\alpha_1\omega_{I2})-\zeta_k(\rho_2c_2^2\omega_{I1}-\rho_1c_1^2\omega_{I2})\Big]\Big\}\mathcal M.$$

**pT-relaxed (four-equation) model** (eq. 18): the partial densities, momentum and $E$ are transported
with sources $\pm\mathcal M$. The mixture pressure law is eq. 12 plus $T_1(p,\rho_1)=T_2(p,\rho_2)$.

Sound speed (eq. 19):
$$\frac1{c_{pT}^2}=\frac1{c_p^2}+\frac{\rho TC_{p1}C_{p2}}{C_{p1}+C_{p2}}\Big(\frac{\Gamma_2}{\rho_2c_2^2}-\frac{\Gamma_1}{\rho_1c_1^2}\Big)^2.$$

Primitive equations (eqs. 20–22):
- $\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1+\rho c_{pT}^2[\alpha_1\alpha_2(\tfrac1{\rho_2c_2^2}-\tfrac1{\rho_1c_1^2})+\tfrac{TC_{p1}C_{p2}}{C_{p1}+C_{p2}}(\tfrac{\Gamma_2}{\rho_2c_2^2}-\tfrac{\Gamma_1}{\rho_1c_1^2})(\tfrac{\alpha_1\Gamma_2}{\rho_2c_2^2}+\tfrac{\alpha_2\Gamma_1}{\rho_1c_1^2})]\nabla\cdot\vec u=\mathcal MS_\alpha$
- $\partial_tT+\vec u\cdot\nabla T+\tfrac{\rho c_{pT}^2T}{C_{p1}+C_{p2}}(\tfrac{C_{p1}\Gamma_1}{\rho_1c_1^2}+\tfrac{C_{p2}\Gamma_2}{\rho_2c_2^2})\nabla\cdot\vec u=\mathcal MS_T$
- $\partial_tp+\vec u\cdot\nabla p+\rho c_{pT}^2\nabla\cdot\vec u=\mathcal MS_p$

with (eq. 23, **all four verified numerically**):
$$
\begin{aligned}
S_\alpha&=\tfrac1{D_T}\Big[\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(\alpha_1\phi_1+\alpha_2\phi_2)+\alpha_1\alpha_2\big(\tfrac{\chi_1}{\Gamma_1}-\tfrac{\chi_2}{\Gamma_2}\big)(\phi_1\zeta_2-\phi_2\zeta_1)\Big],\\
S_T&=\tfrac1{D_T}\Big[\big(\tfrac{\chi_2}{\Gamma_2}-\tfrac{\chi_1}{\Gamma_1}\big)(\alpha_1\zeta_1\rho_2+\alpha_2\zeta_2\rho_1)+\big(\tfrac{\rho_1c_1^2}{\Gamma_1}-\tfrac{\rho_2c_2^2}{\Gamma_2}\big)(\alpha_1\zeta_1+\alpha_2\zeta_2)+\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(\rho_2-\rho_1)\Big],\\
S_p&=\tfrac1{D_T}\Big[\big(\tfrac{\chi_1}{\Gamma_1}-\tfrac{\chi_2}{\Gamma_2}\big)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1)+\big(\tfrac{\rho_2c_2^2}{\Gamma_2}-\tfrac{\rho_1c_1^2}{\Gamma_1}\big)(\alpha_1\phi_1+\alpha_2\phi_2)\Big],\\
D_T&=\alpha_1\alpha_2\big(\tfrac{\rho_1c_1^2}{\Gamma_1}-\tfrac{\rho_2c_2^2}{\Gamma_2}\big)(\phi_1\zeta_2-\phi_2\zeta_1)+\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1).
\end{aligned}
$$
EOS-agnostic form (Appendix B, eq. 100). $S_\alpha$, $S_T$ and $S_p$ are $(d\tilde q/dw)^{-1}[1,-1,0]^T$ for
$\tilde q=[\alpha_1\rho_1,\alpha_2\rho_2,\mathcal E]$ and $w=[\alpha_1,T,p]$:
$$S_\alpha=\tfrac{C_{\mathcal ET}(\alpha_1\phi_1+\alpha_2\phi_2)-C_{\mathcal Ep}(\alpha_1\zeta_1+\alpha_2\zeta_2)}{D_T},\
S_T=\tfrac{C_{\mathcal ET}(\rho_2-\rho_1)+(\mathcal E_1-\mathcal E_2)(\alpha_1\zeta_1+\alpha_2\zeta_2)}{D_T},\
S_p=\tfrac{C_{\mathcal Ep}(\rho_1-\rho_2)-(\mathcal E_1-\mathcal E_2)(\alpha_1\phi_1+\alpha_2\phi_2)}{D_T},$$
where
- $C_{\mathcal Ep}=\sum\alpha_k(\partial\mathcal E_k/\partial T)_p$ and $C_{\mathcal ET}=\sum\alpha_k(\partial\mathcal E_k/\partial p)_T$
- $(\partial\mathcal E_k/\partial T_k)_{p_k}=-\chi_k\phi_k/\Gamma_k$ and $(\partial\mathcal E_k/\partial p_k)_{T_k}=(1-\chi_k\zeta_k)/\Gamma_k$ (eq. 101)
- $d\tilde q/dw=\begin{bmatrix}\rho_1&\alpha_1\phi_1&\alpha_1\zeta_1\\-\rho_2&\alpha_2\phi_2&\alpha_2\zeta_2\\\mathcal E_1-\mathcal E_2&C_{\mathcal Ep}&C_{\mathcal ET}\end{bmatrix}$ (eq. 97)

The same derivation holds at a constant temperature difference $T_2-T_1=\Delta T$, with phase-$k$
quantities taken at $(p,T_k)$.

**pTg-relaxed (HEM)** (eq. 24): $\rho$, $\rho\vec u$ and $E$ are transported, with $g_1(p,T)=g_2(p,T)$.
Sound speed (eq. 25):
$$\frac1{c_{pTg}^2}=\frac1{c_{pT}^2}+\frac{\rho T}{C_{p1}+C_{p2}}\Big[\frac{\Gamma_1C_{p1}}{\rho_1c_1^2}+\frac{\Gamma_2C_{p2}}{\rho_2c_2^2}-\frac1T\Big(\frac{dT}{dp}\Big)_{sat}(C_{p1}+C_{p2})\Big]^2.$$

### 2.3 Equations of state [P22 Sec. 4]
Both phases must cover the whole thermodynamic domain, because pure phases are modelled as
$\alpha_k=1-\epsilon$.

**Mie–Grüneisen** (eq. 27): $p=\Gamma(\rho)(\mathcal E-\rho\varepsilon_r(\rho))+p_r(\rho)$. The pressure-equilibrium
mixture pressure is explicit (eq. 28):
$$p=\frac{\mathcal E-(\alpha_1\rho_1\varepsilon_{r1}+\alpha_2\rho_2\varepsilon_{r2})+(\alpha_1p_{r1}/\Gamma_1+\alpha_2p_{r2}/\Gamma_2)}{\alpha_1/\Gamma_1+\alpha_2/\Gamma_2}.$$

**JWL** (eq. 29):
- $\Gamma=\Gamma_0$
- $\varepsilon_r=\frac{a}{r_1\rho_0}e^{-r_1\rho_0/\rho}+\frac{b}{r_2\rho_0}e^{-r_2\rho_0/\rho}-\varepsilon_0$
- $p_r=ae^{-r_1\rho_0/\rho}+be^{-r_2\rho_0/\rho}$

**NASG** (eqs. 30–34):
- $p=\frac{\gamma-1}{1-\rho b}(\mathcal E-\eta\rho)-\gamma\varpi$
- $T=\frac{(1-\rho b)(p+\varpi)}{\kappa_v\rho(\gamma-1)}$
- $\Gamma=\frac{\gamma-1}{1-\rho b}$, $\varepsilon_r=\eta$, $p_r=-\gamma\varpi$
- $s=\kappa_v\ln\frac{T^\gamma}{(p+\varpi)^{\gamma-1}}+\tilde\eta$
- $h=\kappa_pT+bp+\eta$ with $\kappa_p=\gamma\kappa_v$
- $g=(\gamma\kappa_v-\tilde\eta)T-\kappa_vT\ln\frac{T^\gamma}{(p+\varpi)^{\gamma-1}}+\eta+bp$
- $c=\sqrt{\gamma(p+\varpi)/(\rho(1-\rho b))}$
- $\phi=-\frac{\rho}{T}(1-\rho b)$
- $\zeta=\frac{1-\rho b}{T(\gamma-1)\kappa_v+b(p+\varpi)}\ \big(=\rho(1-\rho b)/(p+\varpi)\big)$

Setting $b=0$ gives the stiffened gas (SG).

**Saturation curve** (eqs. 35–36), from $p_1=p_2$, $T_1=T_2$, $g_1=g_2$, both phases NASG:
$$A_s+\frac{B_s}{T}+C_s\ln T+D_s\ln(p+\varpi_1)-\ln(p+\varpi_2)+\frac{pE_s}{T}=0,$$
$$A_s=\frac{\kappa_{p1}-\kappa_{p2}-\tilde\eta_1+\tilde\eta_2}{\kappa_{p2}-\kappa_{v2}},\ B_s=\frac{\eta_1-\eta_2}{\kappa_{p2}-\kappa_{v2}},\ C_s=\frac{\kappa_{p2}-\kappa_{p1}}{\kappa_{p2}-\kappa_{v2}},\ D_s=\frac{\kappa_{p1}-\kappa_{v1}}{\kappa_{p2}-\kappa_{v2}},\ E_s=\frac{b_1-b_2}{\kappa_{p2}-\kappa_{v2}}.$$

**Parameter tables** ($\varpi$ in Pa, $\eta$ in J/kg, $\tilde\eta$ in J/(kg K), $\kappa_v$ in J/(kg K), $b$ in m³/kg):

| [P22] Table | Fluid | Phase | $\gamma$ | $\varpi$ | $\eta$ | $\tilde\eta$ | $\kappa_v$ | $b$ |
|---|---|---|---|---|---|---|---|---|
| 3 (NASG, 300–500 K) | water | liquid | 1.187 | 7028e5 | −1 177 788 | 0 | 3610 | 6.61e-4 |
| | | vapor | 1.467 | 0 | 2 077 616 | 14 317 | 955 | 0 |
| 4 (NASG, 400–600 K) | dodecane | liquid | 1.206 | 1681e5 | −996 054 | 0 | 2532 | 7.51e-4 |
| | | vapor | 1.021 | 0 | −384 592 | −4301 | 2274 | 0 |
| 5 (NASG, 350–550 K) | water | liquid | 1.387 | 8899e5 | −1 244 191 | 0 | 3202 | 4.78e-4 |
| | | vapor | 1.954 | 0 | 2 287 484 | 6417 | 462 | 0 |

[P22] Table 2, JWL explosive: $\rho_0=1840$ kg/m³, $\Gamma_0=0.25$, $a=854.5\times10^9$ Pa,
$b=20.5\times10^9$ Pa, $r_1=4.6$, $r_2=1.35$, $\varepsilon_0=0$.

### 2.4 Fractional-step algorithm [P22 Sec. 5]
Write $q=[\alpha_1,\alpha_1\rho_1,\alpha_2\rho_2,\rho\vec u,\alpha_1E_1,\alpha_2E_2]$. Time scales are
assumed ordered: $\tau_\mu\ll\tau_\vartheta\ll\tau_\nu$.
1. Homogeneous hyperbolic step → state 0.
2. Relaxation steps:
   - (a) instantaneous mechanical relaxation, $\partial_tq=\psi_\mu$ → state $*$
   - (b) thermal, $\partial_tq=\psi_\mu+\psi_\vartheta$ with $\mu\to\infty$ → state $**$
   - (c) chemical, $\partial_tq=\psi_\mu+\psi_\vartheta+\psi_\nu$ with $\mu\to\infty$ → state $\otimes$

Step (a) always runs. Steps (b) and (c) are optional, and can be limited to interface cells
$\min_k\alpha_k>\epsilon$ (e.g. $10^{-6}$). When thermal and chemical relaxation run unconditionally, the
scheme approximates HEM.

**Mixture-energy consistency** (Sec. 5.1):
- (i) $E^0=E^{0,C}$: conservative total energy, guaranteed by HLLC.
- (ii) The relaxed state satisfies $\mathcal E^0=\alpha_1^\#\mathcal E_1(p^\#,\rho_1^\#)+\alpha_2^\#\mathcal E_2(p^\#,\rho_2^\#)$,
  with $\mathcal E^0=E^0-(\rho\vec u)^0\cdot(\rho\vec u)^0/(2\rho^0)$ (eq. 43).

**Invariants of every relaxation** (eqs. 47–49): $\rho$, $\vec u$, $E$ and $\mathcal E$ are constant. When
chemical relaxation is off, the $\alpha_k\rho_k$ are constant too. The mechanical and thermal steps
therefore have two unknowns ($\alpha_1$, $p$); the chemical step has three.

### 2.5 Instantaneous procedures (previous work, recalled in Sec. 7.1)
- **Mechanical:** invariants plus $p_1=p_2$. Integrate $\partial_t\mathcal E_k=(-1)^kp_I\partial_t\alpha_1$ with
  $p_I=\beta_pp^0+(1-\beta_p)p^*$, where $\beta_p=1/2$ (Pelanti & Shyue 2014) or $0$ (Pelanti 2017). For
  SG this gives a quadratic in $p^*$.
- **Thermal:** add $T_1=T_2$. For SG this gives a quadratic in $p^{**}$.
- **Chemical:** add $g_1=g_2$. The system has three unknowns and generally needs iteration, even for SG.

### 2.6 Arbitrary-rate relaxation [P22 Sec. 7.2] — the core of the method
General idea: write the relaxation ODEs of the parent model (for $\mathcal P$) or of the relaxed models
(for $\mathcal Q$ and $\mathcal M$), so that each process keeps the faster equilibria. Freeze the
coefficients at the start state; the potential difference then decays exponentially, and $\alpha$ (and
$\alpha_1\rho_1$) follow by quadrature. The pressure always comes from the mixture energy relation
(eq. 43), not from the exponential.

**Mechanical** (eqs. 50–59). From eq. 1a and eq. 8:
$$\partial_t\alpha_1=\mu(p_1-p_2),\qquad \partial_tp_1=\frac{\mu}{\xi_1^p}(p_2-p_1),\qquad \partial_tp_2=-\frac{\mu}{\xi_2^p}(p_2-p_1),$$
$$\frac1{\xi_k^p}=\frac1{\alpha_k}[\Gamma_k(\mathcal E_k+p_I)+\chi_k\rho_k]=\frac1{\alpha_k}[\Gamma_k(p_I-p_k)+\rho_kc_k^2]\quad\text{(eq. 51)}.$$
- Freeze $\xi_k^p=\xi_k^{p0}$ and let $\Delta p=p_2-p_1$.
- Then $\Delta p^*=\Delta p^0e^{-K_p\Delta t}$ with $K_p=\mu(1/\xi_1^{p0}+1/\xi_2^{p0})>0$.
- Volume fraction (eq. 56):
$$\alpha_1^*=\alpha_1^0-\frac{\Delta p^0}{1/\xi_1^{p0}+1/\xi_2^{p0}}\big(1-e^{-K_p\Delta t}\big).$$
- As $\mu\to\infty$, the exponential is dropped (eq. 57), and $p^*_{\mu\to\infty}=(\xi_1p_1^0+\xi_2p_2^0)/(\xi_1+\xi_2)$ (eq. 58).
- That $p^*$ is **not** used, because it is not energy-consistent. Instead, $p^*$ solves
  $\mathcal E^0=\alpha_1^*\mathcal E_1(p^*,(\alpha_1\rho_1)^0/\alpha_1^*)+\alpha_2^*\mathcal E_2(p^*,(\alpha_2\rho_2)^0/\alpha_2^*)$ (eq. 59),
  which is explicit for Mie–Grüneisen (eq. 28).
- Pelanti uses only $\mu\to\infty$.
- Appendix D proves that this procedure keeps $u$ and $p$ invariant across an isolated material
  interface, at first order for SG.

**Thermal** (eqs. 60–70), from the p-relaxed model (eqs. 9, 15), starting from state $*$:
$$\partial_t\alpha_1=\vartheta\frac ZD(T_2-T_1),\qquad
\partial_tT_1=\frac{\vartheta}{\phi_1D}\Big[-\frac{\rho_1}{\alpha_1}Z-\zeta_1(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big](T_2-T_1),$$
$$\partial_tT_2=\frac{\vartheta}{\phi_2D}\Big[\frac{\rho_2}{\alpha_2}Z-\zeta_2(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)\Big](T_2-T_1),$$
- $Z=\alpha_2\Gamma_1+\alpha_1\Gamma_2$ (eq. 61)
- Inverse capacities (eq. 62):
$$\frac1{\xi_k^T}=-\frac{\frac{\rho_k}{\alpha_k}Z+(-1)^{k-1}\zeta_k(\Gamma_1\rho_2c_2^2-\Gamma_2\rho_1c_1^2)}{\phi_kD}$$
  so $\partial_tT_1=\frac{\vartheta}{\xi_1^T}\Delta T$ and $\partial_tT_2=-\frac{\vartheta}{\xi_2^T}\Delta T$, with $\Delta T=T_2-T_1$.
- With frozen coefficients: $\Delta T^{**}=\Delta T^*e^{-K_T\Delta t}$, $K_T=\vartheta(1/\xi_1^{T*}+1/\xi_2^{T*})>0$.
- Volume fraction, freezing $Z/D$ as well (eq. 67):
$$\alpha_1^{**}=\alpha_1^*+\Big(\frac ZD\Big)^*\frac{\Delta T^*}{1/\xi_1^{T*}+1/\xi_2^{T*}}\big(1-e^{-K_T\Delta t}\big).$$
- Limit $\vartheta\to\infty$ (eq. 68): $T=(\xi_1T_1+\xi_2T_2)/(\xi_1+\xi_2)$ (eq. 69).
- $p^{**}$ from the energy relation at $\alpha_1^{**}$ (eq. 70).
- If $\vartheta$ depends on the state, it is frozen at $*$.

**Chemical** (eqs. 71–82), from the pT-relaxed model, starting from state $**$. Use
$\partial_tg_k=\frac1{\rho_k}\partial_tp-s_k\partial_tT$:
$$\partial_t\alpha_1=\nu S_\alpha\Delta g,\quad \partial_tg_1=\nu\Big(\frac{S_p}{\rho_1}-s_1S_T\Big)\Delta g,\quad \partial_tg_2=\nu\Big(\frac{S_p}{\rho_2}-s_2S_T\Big)\Delta g,\quad \partial_t(\alpha_1\rho_1)=\nu\Delta g,$$
with $\Delta g=g_2-g_1$; $g_I$ drops out.
- Define (corrected, see §8):
$$\frac1{\xi_1^g}=\frac{S_p}{\rho_1}-s_1S_T,\qquad \frac1{\xi_2^g}=-\Big(\frac{S_p}{\rho_2}-s_2S_T\Big),$$
  so that $\partial_t\Delta g=-\nu(1/\xi_1^g+1/\xi_2^g)\Delta g$, with
  $1/\xi_1^g+1/\xi_2^g=S_p(1/\rho_1-1/\rho_2)-S_T(s_1-s_2)$.
- Then $\Delta g^\otimes=\Delta g^{**}e^{-K_g\Delta t}$ with $K_g=\nu(1/\xi_1^{g**}+1/\xi_2^{g**})>0$.
- Partial density and volume fraction (eqs. 77–78):
$$(\alpha_1\rho_1)^\otimes=(\alpha_1\rho_1)^{**}+\frac{\Delta g^{**}}{1/\xi_1^{g**}+1/\xi_2^{g**}}\big(1-e^{-K_g\Delta t}\big),\qquad
\alpha_1^\otimes=\alpha_1^{**}+S_\alpha^{**}\frac{\Delta g^{**}}{1/\xi_1^{g**}+1/\xi_2^{g**}}\big(1-e^{-K_g\Delta t}\big).$$
- Limit $\nu\to\infty$: drop the exponential (eqs. 79–80), giving
  $g_\infty=(\xi_1g_1+\xi_2g_2)/(\xi_1+\xi_2)$ (eq. 81).
- $p^\otimes$ from the energy relation with the new partial densities (eq. 82).
- If $\nu$ depends on the state, it is frozen at $**$.
- The equations also hold for chemical relaxation at a constant $T_2-T_1$.

**Complete evaporation or condensation** (end of Sec. 7.2.3). If $\alpha_1^\otimes\notin(0,1)$:
- Set $\bar\alpha_1=\epsilon_\alpha$ when $s_2^{**}>s_1^{**}$ (phase 1 vanishes), else $1-\epsilon_\alpha$, with
  $\epsilon_\alpha=10^{-8}$.
- Solve for the rate (eq. 83):
$$\bar K_g=-\frac1{\Delta t}\ln\Big(1-\frac{\bar\alpha_1-\alpha_1^{**}}{S_\alpha^{**}\Delta g^{**}/(1/\xi_1^{g**}+1/\xi_2^{g**})}\Big),$$
- then update $(\alpha_1\rho_1)^\otimes$ with eq. 77 using $\bar K_g$.

The same clamp $\alpha_1=\bar\alpha_1$ is used if the mechanical or thermal steps give an unphysical $\alpha_1$.

### 2.7 Homogeneous solver (Sec. 6, Appendix C)
Second-order wave-propagation finite volume (LeVeque):
- Update (eq. 44): $Q_i^{n+1}=Q_i^n-\frac{\Delta t}{\Delta x}(\mathcal A^+\Delta Q_{i-1/2}+\mathcal A^-\Delta Q_{i+1/2})-\frac{\Delta t}{\Delta x}(F^h_{i+1/2}-F^h_{i-1/2})$
- Fluctuations: $\mathcal A^\pm\Delta Q=\sum_l(s^l)^\pm\mathcal W^l$
- Second-order correction: $F^h=\frac12\sum_l|s^l|(1-\frac{\Delta t}{\Delta x}|s^l|)\tilde{\mathcal W}^l$, with a minmod limiter.

HLLC-type solver: three waves at speeds $S_\ell,S^\star,S_r$, with $u$ and $p_m=\alpha_1p_1+\alpha_2p_2$
continuous across $S^\star$.
- Davis wave speeds (eq. 112): $S_\ell=\min(u_\ell-c_{f\ell},u_r-c_{fr})$, $S_r=\max(u_\ell+c_{f\ell},u_r+c_{fr})$.
- Contact speed (eq. 111):
$$S^\star=\frac{p_r-p_\ell+\rho_\ell u_\ell(S_\ell-u_\ell)-\rho_ru_r(S_r-u_r)}{\rho_\ell(S_\ell-u_\ell)-\rho_r(S_r-u_r)}.$$
- Middle states, $\iota=\ell,r$ (eq. 114). $\alpha_{1}$ is unchanged across the outer waves, and the rest are
$$\Big[(\alpha_k\rho_k)_\iota\tfrac{S_\iota-u_\iota}{S_\iota-S^\star},\ \rho_\iota\tfrac{S_\iota-u_\iota}{S_\iota-S^\star}S^\star,\ (\alpha_k\rho_k)_\iota\tfrac{S_\iota-u_\iota}{S_\iota-S^\star}\Big(\tfrac{E_{k,\iota}}{\rho_{k,\iota}}+(S^\star-u_\iota)\big(S^\star+\tfrac{p_{k,\iota}}{\rho_{k,\iota}(S_\iota-u_\iota)}\big)\Big)\Big].$$
- The non-conservative terms are neglected in the jump conditions across the outer waves. Across the
  middle wave the jump $-S^\star((\alpha_kp_k)^{\star r}-(\alpha_kp_k)^{\star\ell})$ approximates $-u\,\partial_x(\alpha_kp_k)$.
- [P22] corrects Pelanti & Shyue 2014: the Rankine–Hugoniot conditions for $\alpha_kE_k$ do not hold
  across the middle wave.

### 2.8 Test cases [P22 Sec. 8]
**8.1.1 Detonation gas–water shock tube**
- Domain $[-10,10]$ m, discontinuity at $x=0$.
- Left: JWL (Table 2), $\rho=2000$, $p=4.6406\times10^{10}$ Pa.
- Right: NASG water (Table 3), $\rho=1044$, $p=10^5$.
- Each region is nearly pure, with $\alpha=1-10^{-8}$; $u=0$.
- 500 cells, CFL 0.4, $t=0.9$ ms. Mechanical relaxation only. Compared against the exact p-relaxed solution.

**8.1.2 UNDEX near a rigid wall**
- Domain $[-1.35,1.35]\times[-1.5,0]$ m, wall at $y=0$.
- Bubble centred at $(0,-0.32)$, $r=0.05$: $p=4.6406\times10^{10}$, $\rho_{gas}=2000$, $\alpha_{gas}=1-10^{-6}$.
- Outside: $p=10^5$, $\rho_{liq}=1044$, $\alpha_{liq}=1-10^{-6}$.
- JWL + NASG Table 3; 721×400 cells, CFL 0.4.
- Outputs: pressure at 0, 0.03, 0.06, 0.08, 0.21, 0.3 ms, and wall-centre pressure history to 0.4 ms.

**8.2.1 Water cavitation tube** (two-phase)
- Domain $[0,1]$ m, $\alpha_{vap}=10^{-2}$, $p=10^5$, $T=353$ K (below $T_{sat}=372.3$ K).
- $u=-2$ m/s left of $x=0.5$, $+2$ m/s right.
- NASG Table 3; 5000 cells, CFL 0.5, $t=3$ ms.
- Five cases:
  - p-relax
  - pT(f), $\vartheta=2000$
  - pT
  - pTg(f), $\nu=10^{-4}$
  - pTg
- Chemical relaxation only when $T>T_{sat}(p)$, i.e. evaporation only (following Saurel et al. 2008).
- The pT cases match the exact pT-relaxed solution; with mass transfer, $p$ goes to $p_{sat}$ in the centre.

**8.2.2 Dodecane liquid–vapor shock tube**
- Unit length, discontinuity at $x=0.75$.
- Left: liquid, $p=10^8$, $\rho_{liq}=500$. Right: vapor, $p=10^5$, $\rho_{vap}=2$.
- $\alpha=10^{-8}$ of the opposite phase on each side; $u=0$.
- NASG Table 4; 2000 cells, CFL 0.5, $t=473\ \mu$s.
- Thermo-chemical relaxation only where $\min\alpha>10^{-4}$, and only when $T>T_{sat}$.
- Complete evaporation occurs; it is handled with eq. 83.

**8.2.3 Simões-Moreira & Shepherd dodecane evaporation fronts**
- SG dodecane with the parameters of Pelanti & Shyue 2014; instantaneous relaxation.
- Low-pressure chamber at 1 mbar; tube at $p_{sat}(T_{test})$; 8 values of $T_{test}$ between 180 and 300 °C.
- Front velocity $U=((\rho u)_a-(\rho u)_b)/(\rho_a-\rho_b)$, averaged over 400–500 μs.

**8.2.4 Barták blowdown**
- Water at $p=12.5$ MPa, $T=563.15$ K, $\alpha_{vap}=10^{-3}$; tube 1700 mm, closed on the right.
- Left boundary open to $10^5$ Pa. NASG Table 5; 1000 cells, CFL 0.5, $t=17.5$ ms; gauge at $x=48$ mm.
- $\vartheta\to\infty$, with a state-dependent $\nu$ (eq. 84):
  $\nu=C_r\alpha_2^{0.6}\big(\frac{p_{sat}-p}{p_{crit}-p_{sat}}\big)^{1.76}$, $C_r=1$, $p_{crit}=22.064$ MPa.
- Captures the metastable pressure undershoot; instantaneous $\nu$ cannot.

**8.2.5 High-pressure fuel injector (2-D)**
- Nozzle 10 cm long, 4 cm high, 1.2 cm throat; converging wall 45°, diverging wall 10°.
- Discontinuity at $x=0.8$ cm: liquid dodecane at $p=10^8$, $T=550$ K | vapor at $p=10^5$, $\rho_{vap}=5$.
- $\alpha_{vap}=10^{-4}$ in the liquid and $\alpha_{liq}=10^{-6}$ in the vapor; initially in thermal equilibrium.
- NASG Table 4. Activation where $\min\alpha>0.9\times10^{-4}$; $\vartheta\to\infty$; $\nu\in\{0,0.1,25,\infty\}$.
- 400×160 cells, CFL 0.4 for the $t=0.6$ ms fields; 200×80 cells for the $Y_{vap}$ snapshots at 0.27, 0.6, 1.2 ms.

---

## 3. [P26]/[P24] Three-phase liquid–vapor–inert gas model

### 3.1 Equations ([P26] eq. 1; [P24] eq. 1)
$$
\begin{aligned}
&\partial_t\alpha_1+\vec u\cdot\nabla\alpha_1=\mathcal P_{12}+\mathcal P_{13},\qquad
\partial_t\alpha_3+\vec u\cdot\nabla\alpha_3=\mathcal P_{31}+\mathcal P_{32},\\
&\partial_t(\alpha_1\rho_1)+\nabla\cdot(\alpha_1\rho_1\vec u)=\mathcal M,\quad
\partial_t(\alpha_2\rho_2)+\nabla\cdot(\alpha_2\rho_2\vec u)=-\mathcal M,\quad
\partial_t(\alpha_3\rho_3)+\nabla\cdot(\alpha_3\rho_3\vec u)=0,\\
&\partial_t(\rho\vec u)+\nabla\cdot\big(\rho\vec u\otimes\vec u+\textstyle\sum_k\alpha_kp_k\mathbb I\big)=0,\\
&\partial_t(\alpha_1E_1)+\nabla\cdot(\alpha_1(E_1+p_1)\vec u)+\Sigma_1=-p_{I12}\mathcal P_{12}-p_{I13}\mathcal P_{13}+\mathcal Q_{12}+\mathcal Q_{13}+(h_I+|\vec u|^2/2)\mathcal M,\\
&\partial_t(\alpha_2E_2)+\nabla\cdot(\alpha_2(E_2+p_2)\vec u)+\Sigma_2=-p_{I12}\mathcal P_{21}-p_{I23}\mathcal P_{23}+\mathcal Q_{21}+\mathcal Q_{23}-(h_I+|\vec u|^2/2)\mathcal M,\\
&\partial_t(\alpha_3E_3)+\nabla\cdot(\alpha_3(E_3+p_3)\vec u)+\Sigma_3=-p_{I13}\mathcal P_{31}-p_{I23}\mathcal P_{32}+\mathcal Q_{31}+\mathcal Q_{32},
\end{aligned}
$$
- $\Sigma_k=\vec u\cdot\vec\Upsilon_k$ with $\vec\Upsilon_k=Y_k\nabla(\sum_j\alpha_jp_j)-\nabla(\alpha_kp_k)$, and $\sum_k\Sigma_k=0$.
- [P26] writes the mass-transfer energy term with an interface **enthalpy** $h_I$; [P24] writes
  $g_I^t=g_I+|\vec u|^2/2$. Neither is needed numerically.
- Relaxation terms (eq. 2): $\mathcal P_{kj}=\mu_{kj}(p_k-p_j)$, $\mathcal Q_{kj}=\vartheta_{kj}(T_j-T_k)$, $\mathcal M=\nu(g_2-g_1)$.
- Mechanical relaxation is always instantaneous: $\mu_{kj}\equiv\mu\to\infty$.
- Thermal cases:
  - (i) all three phases: $\vartheta_{kj}\to\infty$
  - (ii) liquid–vapor only: $\vartheta_{12}\to\infty$, $\vartheta_{3k}=0$; reasonable when a resolved interface separates gas from liquid–vapor
  - (iii) none: $\vartheta_{kj}=0$
- $\nu\ge0$ is arbitrary, and mass transfer occurs at liquid–vapor thermal equilibrium.
- Mixture energy relation: $\mathcal E=\sum_k\alpha_k\mathcal E_k(p,\rho_k)$ (eq. 3).
- $c_f^2=\sum Y_kc_k^2$ (eq. 4).

**Physics difference from a Dalton gas mixture.** The inert gas is a separate phase with its own volume
fraction, and possibly its own temperature. The vapor chemical potential is evaluated at the full
pressure $p$, not at a partial pressure. In cavitation zones the **total** pressure therefore tends to
$p_{sat}$.

### 3.2 EOS used in [P26]: stiffened gas
- $p=(\gamma-1)(\mathcal E-\eta\rho)-\gamma\varpi$, $T=(p+\varpi)/(\kappa_v\rho(\gamma-1))$ (eq. 6)
- $s$, $h$, $g$ as for NASG with $b=0$ (eq. 7)
- Ideal gas: $\eta=\varpi=0$.
- Saturation curve: eqs. 8–9, i.e. eqs. 35–36 without $E_s$.
- All-SG mixture pressure (eq. 10):
$$p=\frac{\mathcal E-\sum_k\alpha_k\rho_k\eta_k-\sum_k\alpha_k\gamma_k\varpi_k/(\gamma_k-1)}{\sum_k\alpha_k/(\gamma_k-1)}.$$
- SG coefficients: $\Gamma=\gamma-1$, $\chi=-(\gamma-1)\eta$, $c^2=\gamma(p+\varpi)/\rho$, $\phi=-\rho/T$, $\zeta=\rho/(p+\varpi)$.

| [P26] Table | Fluid | Phase | $\gamma$ | $\varpi$ [Pa] | $\eta$ [J/kg] | $\tilde\eta$ [J/(kg K)] | $\kappa_v$ [J/(kg K)] |
|---|---|---|---|---|---|---|---|
| 3 | water | liquid | 2.35 | 1e9 | −1167e3 | 0 | 1816 |
| | | vapor | 1.43 | 0 | 2030e3 | −23.4e3 | 1040 |
| 4 | dodecane | liquid | 2.35 | 4e8 | −775.269e3 | 0 | 1077.7 |
| | | vapor | 1.025 | 0 | −237.547e3 | −24.4e3 | 1956.45 |

Air: ideal gas, $\gamma=1.4$, $\kappa_v=718$. [P24] uses $\gamma_g=2$ for the UNDEX combustion gas.

### 3.3 Relaxed models ([P26] Sec. 3, Appendix B)
The p-relaxed model (7 equations) is given in [P26] Appendix A. Below, $c_p$ is the three-phase Wood
speed, $c_p^{-2}=\rho\sum\alpha_k/(\rho_kc_k^2)$ (eq. 13).

**pT-relaxed model**, all phases at one $T$ (eqs. 11–17). Transported: $\alpha_k\rho_k$ (sources
$\pm\mathcal M$, 0), $\rho\vec u$ and $E$. Sound speed (eq. 12):
$$\frac1{c_{pT}^2}=\frac1{c_p^2}+\frac{\rho T}{\sum C_{pk}}\Big[C_{p1}C_{p2}\big(\tfrac{\Gamma_2}{\rho_2c_2^2}-\tfrac{\Gamma_1}{\rho_1c_1^2}\big)^2+C_{p1}C_{p3}\big(\tfrac{\Gamma_3}{\rho_3c_3^2}-\tfrac{\Gamma_1}{\rho_1c_1^2}\big)^2+C_{p2}C_{p3}\big(\tfrac{\Gamma_3}{\rho_3c_3^2}-\tfrac{\Gamma_2}{\rho_2c_2^2}\big)^2\Big].$$
Primitive form (eq. 14): $\partial_t\alpha_k+\vec u\cdot\nabla\alpha_k+\rho c_{pT}^2\mathcal A_{\alpha_k}\nabla\cdot\vec u=\mathcal S_{\alpha_k}\mathcal M$
for $k=1,3$; $\partial_tp+\dots+\rho c_{pT}^2\nabla\cdot\vec u=\mathcal S_p\mathcal M$;
$\partial_tT+\dots+\rho c_{pT}^2\mathcal A_T\nabla\cdot\vec u=\mathcal S_T\mathcal M$.

Mass-transfer coefficients (eqs. 16–17). **These are verified numerically as printed**, against a
finite difference of the three-phase pT equilibrium at fixed $\rho$, $\mathcal E$ and $\alpha_3\rho_3$.
Write $X_k\equiv\chi_k/\Gamma_k$ and $W_k\equiv\rho_kc_k^2/\Gamma_k$.
$$
\begin{aligned}
\mathcal S_{\alpha_1}=\tfrac1{D_T}\Big[&\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}+\tfrac{\alpha_3}{\Gamma_3}\big)(\alpha_1\phi_1\rho_3+\alpha_2\phi_2\rho_3+\alpha_3\phi_3\rho_2)
+\alpha_1\alpha_2\rho_3(X_1-X_2)(\zeta_2\phi_1-\phi_2\zeta_1)\\
&+\alpha_1\alpha_3\rho_2(X_1-X_3)(\zeta_3\phi_1-\phi_3\zeta_1)+\alpha_1\alpha_3\rho_2(X_3-X_2)(\zeta_3\phi_1-\phi_3\zeta_1)\\
&+\alpha_3(W_2-W_3)\big(\alpha_1(\zeta_3\phi_1-\phi_3\zeta_1)+\alpha_2(\phi_2\zeta_3-\phi_3\zeta_2)\big)\Big],\\
\mathcal S_{\alpha_3}=\tfrac1{D_T}\Big[&\tfrac{\alpha_3^2\phi_3}{\Gamma_3}(\rho_1-\rho_2)+\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(-\alpha_3\phi_3\rho_2+\alpha_3\phi_3\rho_1)
+\alpha_3\alpha_2\rho_1(X_1-X_2)(\zeta_2\phi_3-\phi_2\zeta_3)\\
&-\alpha_3\alpha_1\rho_2(X_2-X_1)(-\phi_1\zeta_3+\phi_3\zeta_1)+\alpha_3(W_2-W_1)\big(\alpha_1(-\phi_1\zeta_3+\phi_3\zeta_1)+\alpha_2(-\phi_2\zeta_3+\phi_3\zeta_2)\big)\Big],\\
\mathcal S_p=\tfrac1{D_T}\Big[&\rho_3(X_1-X_2)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1)+\rho_3(\alpha_1\phi_1+\alpha_2\phi_2)(W_2-W_1)
-\alpha_3\phi_3W_3(\rho_1-\rho_2)\\
&+\alpha_3\phi_3\rho_1\rho_2\Big(\tfrac{c_2^2-\chi_2}{\Gamma_2}-\tfrac{c_1^2-\chi_1}{\Gamma_1}\Big)\Big],\\
\mathcal S_T=\tfrac1{D_T}\Big[&\rho_3(\alpha_1\zeta_1\rho_2+\alpha_2\zeta_2\rho_1)(X_2-X_1)+\rho_3(W_1-W_2)(\alpha_1\zeta_1+\alpha_2\zeta_2)
+\rho_3\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(\rho_2-\rho_1)\\
&+\tfrac{\rho_3\alpha_3}{\Gamma_3}(c_3^2\zeta_3-1)(\rho_1-\rho_2)-\rho_1\rho_2\alpha_3\zeta_3\Big(\tfrac{c_2^2-\chi_2}{\Gamma_2}-\tfrac{c_1^2-\chi_1}{\Gamma_1}\Big)\Big],\\
D_T=\ &\rho_3\alpha_1\alpha_2(W_1-W_2)(\phi_1\zeta_2-\phi_2\zeta_1)+\rho_3\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_2}{\Gamma_2}\big)(\alpha_1\phi_1\rho_2+\alpha_2\phi_2\rho_1)\\
&+\rho_2\alpha_1\alpha_3(W_1-W_3)(\phi_1\zeta_3-\phi_3\zeta_1)+\rho_2\big(\tfrac{\alpha_1}{\Gamma_1}+\tfrac{\alpha_3}{\Gamma_3}\big)(\alpha_1\phi_1\rho_3+\alpha_3\phi_3\rho_1)\\
&+\rho_1\alpha_3\alpha_2(W_3-W_2)(\phi_3\zeta_2-\phi_2\zeta_3)+\rho_1\big(\tfrac{\alpha_3}{\Gamma_3}+\tfrac{\alpha_2}{\Gamma_2}\big)(\alpha_3\phi_3\rho_2+\alpha_2\phi_2\rho_3)\\
&-\rho_1\rho_2\alpha_3^2\tfrac{\phi_3}{\Gamma_3}-\rho_1\rho_3\alpha_2^2\tfrac{\phi_2}{\Gamma_2}-\rho_2\rho_3\alpha_1^2\tfrac{\phi_1}{\Gamma_1}.
\end{aligned}
$$
Checks:
- The two $\alpha_1\alpha_3\rho_2$ terms in $\mathcal S_{\alpha_1}$ combine to $\alpha_1\alpha_3\rho_2(X_1-X_2)(\zeta_3\phi_1-\phi_3\zeta_1)$.
- At $\alpha_3=0$: $D_T\to\rho_3D_T^{(2\text{-phase})}$, and $\mathcal S_p$, $\mathcal S_T$, $\mathcal S_{\alpha_1}$ reduce to the [P22] eq. 23 forms.
- [P24] eqs. 7a–e are the same formulas.

Convective coefficients (eq. 15), transcribed as printed and **not verified**; they are not needed for
relaxation:
- $\mathcal A_{\alpha_1}=\frac{\alpha_1\alpha_2}{\rho_2c_2^2}-\frac{\alpha_1(1-\alpha_1)}{\rho_1c_1^2}+\frac{\alpha_1\alpha_3}{\rho_3c_3^2}+\frac{T}{\sum C_{pk}}\big[C_{p1}(C_{p2}+C_{p3})(\alpha_1-1)\frac{\Gamma_1^2}{(\rho_1c_1^2)^2}-2C_{p1}C_{p3}(\alpha_1-\tfrac12)\frac{\Gamma_3\Gamma_1}{\rho_3c_3^2\rho_1c_1^2}+C_{p3}(C_{p1}+C_{p2})\alpha_1\frac{\Gamma_3^2}{(\rho_3c_3^2)^2}-2C_{p1}C_{p2}(\alpha_1-\tfrac12)\frac{\Gamma_2\Gamma_1}{\rho_2c_2^2\rho_1c_1^2}-2C_{p2}C_{p3}\frac{\Gamma_2\Gamma_3}{\rho_2c_2^2\rho_3c_3^2}\alpha_1+C_{p2}(C_{p1}+C_{p3})\alpha_1\frac{\Gamma_2^2}{(\rho_2c_2^2)^2}\big]$
- $\mathcal A_{\alpha_3}$: the same with $1\leftrightarrow3$
- $\mathcal A_T=\frac{T}{\sum C_{pk}}\sum_k\frac{C_{pk}\Gamma_k}{\rho_kc_k^2}$

**p T̃-relaxed model**, with $T_1=T_2=T$ and the gas at its own $T_3$ (eqs. 19–25; new in [P26]). The
gas has no heat exchange in the relaxation ($\mathcal Q_3=0$), so it evolves isentropically. Sound speed
(eq. 20):
$$\frac1{c_{p\tilde T}^2}=\frac1{c_p^2}+\frac{\rho T}{C_{p1}+C_{p2}}C_{p1}C_{p2}\Big(\frac{\Gamma_2}{\rho_2c_2^2}-\frac{\Gamma_1}{\rho_1c_1^2}\Big)^2.$$
Mass-transfer coefficients. These are the **corrected** forms; the printed eqs. 24b–d are wrong, see §8.
All of them, and $\tilde D_T$, were verified against a finite difference of the equilibrium
"liquid–vapor at $T$, gas isentropic at $T_3$, common $p$", at fixed $\rho$, $\mathcal E$ and $\alpha_3\rho_3$.
Define
$$\Lambda_1\equiv\frac{\rho_1c_1^2}{\Gamma_1}-\Big(\frac{c_2^2-\chi_2}{\Gamma_2}+\frac{\chi_1}{\Gamma_1}\Big)\rho_2,\qquad
\Lambda_2\equiv\Big(\frac{\chi_2}{\Gamma_2}+\frac{c_1^2-\chi_1}{\Gamma_1}\Big)\rho_1-\frac{\rho_2c_2^2}{\Gamma_2}.$$
Then
$$
\begin{aligned}
\tilde{\mathcal S}_{\alpha_1}&=\tfrac1{\tilde D_T}\Big[\tfrac{\alpha_1^2}{\Gamma_1}c_3^2\phi_1\rho_3+\Big(-\Big(\big(\tfrac{\chi_2\zeta_2-1}{\Gamma_2}-\tfrac{\zeta_2\chi_1}{\Gamma_1}\big)\phi_1-\big(\tfrac{\zeta_1\chi_2}{\Gamma_2}-\tfrac{\chi_1\zeta_1-1}{\Gamma_1}\big)\phi_2\Big)c_3^2\rho_3\alpha_2+\rho_2\phi_1\alpha_3\big(\tfrac{c_2^2-\chi_2}{\Gamma_2}+\tfrac{\chi_1}{\Gamma_1}\big)\Big)\alpha_1+\tfrac{\alpha_2}{\Gamma_2}\phi_2(\alpha_2c_3^2\rho_3+\alpha_3c_2^2\rho_2)\Big]\ \text{(as printed, correct)},\\
\tilde{\mathcal S}_{\alpha_3}&=\tfrac1{\tilde D_T}\big[\Lambda_1\alpha_3\phi_1\alpha_1+\Lambda_2\alpha_2\phi_2\alpha_3\big],\qquad
\tilde{\mathcal S}_p=-\tfrac{c_3^2\rho_3}{\tilde D_T}\big[\Lambda_1\phi_1\alpha_1+\Lambda_2\alpha_2\phi_2\big],\\
\tilde{\mathcal S}_T&=\tfrac1{\tilde D_T}\Big[c_3^2\rho_3\Big(\Big(\tfrac{\rho_1}{\Gamma_1}(c_1^2\zeta_1-1)-\rho_2\big(\tfrac{\zeta_1(c_2^2-\chi_2)}{\Gamma_2}+\tfrac{\chi_1\zeta_1-1}{\Gamma_1}\big)\Big)\alpha_1+\alpha_2\Big(\big(\tfrac{\chi_2\zeta_2-1}{\Gamma_2}+\tfrac{\zeta_2(c_1^2-\chi_1)}{\Gamma_1}\big)\rho_1-\tfrac{\rho_2}{\Gamma_2}(c_2^2\zeta_2-1)\Big)\Big)\\
&\qquad\ -\Big(\tfrac{c_2^2-\chi_2}{\Gamma_2}-\tfrac{c_1^2-\chi_1}{\Gamma_1}\Big)\rho_2\alpha_3\rho_1\Big],\\
\tilde D_T&=\phi_1\rho_2\tfrac{\alpha_1}{\Gamma_1}(\alpha_1c_3^2\rho_3+\alpha_3\rho_1c_1^2)+\phi_2\rho_1\tfrac{\alpha_2}{\Gamma_2}(\alpha_2c_3^2\rho_3+\alpha_3c_2^2\rho_2)\\
&\quad+\alpha_1\alpha_2\rho_3c_3^2\Big(\big((1-c_2^2\zeta_2)\tfrac{\rho_2}{\Gamma_2}+\tfrac{\rho_1}{\Gamma_1}\zeta_2c_1^2\big)\phi_1+\big(\tfrac{\rho_2}{\Gamma_2}\zeta_1c_2^2+\tfrac{\rho_1}{\Gamma_1}(1-c_1^2\zeta_1)\big)\phi_2\Big)\ \text{(as printed, correct)}.
\end{aligned}
$$
The heat-transfer coefficients $\mathcal B_{(\cdot)}$ (eqs. B.22–B.23) multiply $\mathcal Q_3$. They are
not needed when $\mathcal Q_3=0$, and are not reproduced here. The p T̃ convective coefficients
$\tilde{\mathcal A}$ (eq. 23) are also omitted: not needed for relaxation, and they look asymmetric as printed.

Derivation method for both models: the Chen–Levermore–Liu / Murrone–Guillard asymptotic reduction
(Appendix B).
- Write the p-relaxed model in $w=[\alpha_1,\alpha_3,u,p,T_1,T_2,T_3]$ and identify the stiff thermal part
  $\Psi$.
- Use the Maxwellian $\Xi$: $T_k=T$ for pT, or $T_1=T_2=T$ for p T̃.
- Take $P$ as the first $L$ rows of $C^{-1}$, with $C=[d\Xi,\ \text{range}(\Psi')]$; the reduced system is
  $\partial_tv+PA\,d\Xi\,\partial_xv=P\Phi$.

### 3.4 Numerical method in [P26]
Same fractional-step scheme as [P22], with $q\in\mathbb R^{8+d}$
($\alpha_1,\alpha_3,\alpha_k\rho_k,\rho\vec u,\alpha_kE_k$). The time scales satisfy $\tau_\mu\ll\tau_\vartheta\ll\tau_\nu$.
- **2(a), 2(b): instantaneous, by algebraic equilibrium** (from Pelanti & Shyue 2019).
  - Invariants: $\rho$, $\vec u$, $E$, $\mathcal E$, $\alpha_3\rho_3$, and $\alpha_k\rho_k$ (eqs. 35–36).
  - Mechanical: $p_k=p$, with $p^*_{Ikj}=p^*$ when integrating the phasic energies.
  - Thermal: $T_k=T$ (case i) or $T_1=T_2$ (case ii, plus the integrated gas energy equation).
  - For SG each reduces to a **cubic** in $p$, or a quadratic when two phases have $\varpi=0$.
- **2(c): arbitrary-rate chemical relaxation** (eqs. 37–44). Replace $\partial_tq=\psi_\mu+\psi_\vartheta+\psi_\nu$ by
  $\partial_tq=\hat\psi_\nu$:
$$\partial_t\alpha_1=\hat{\mathcal S}_{\alpha_1}\mathcal M,\quad \partial_t\alpha_3=\hat{\mathcal S}_{\alpha_3}\mathcal M,\quad \partial_t(\alpha_1\rho_1)=\mathcal M,\quad \partial_tT=\hat{\mathcal S}_T\mathcal M,\quad \partial_tp=\hat{\mathcal S}_p\mathcal M,$$
  - Case (i): $\hat{\mathcal S}=\mathcal S$ from the pT model (eqs. 16–17).
  - Case (ii): $\hat{\mathcal S}=\tilde{\mathcal S}$ from the p T̃ model with $\mathcal Q_3=0$ (eqs. 24–25).
  - Then $\partial_t\Delta g=-\nu(1/\xi_1^g+1/\xi_2^g)\Delta g$, with $1/\xi_k^g$ as in §2.6 (corrected) using $\hat{\mathcal S}_p$ and $\hat{\mathcal S}_T$.
  - Freeze at $**$: $\Delta g^\otimes=\Delta g^{**}e^{-K_g\Delta t}$ with $K_g=\nu(1/\xi_1^{g**}+1/\xi_2^{g**})$ (eq. 40).
  - Updates (eq. 41):
$$\alpha_k^\otimes=\alpha_k^{**}+\hat{\mathcal S}^{**}_{\alpha_k}\frac{\Delta g^{**}}{1/\xi_1^{g**}+1/\xi_2^{g**}}(1-e^{-K_g\Delta t})\ (k=1,3),\qquad
(\alpha_1\rho_1)^\otimes=(\alpha_1\rho_1)^{**}+\frac{\Delta g^{**}}{1/\xi_1^{g**}+1/\xi_2^{g**}}(1-e^{-K_g\Delta t}).$$
  - $\nu\to\infty$: drop the exponential (eq. 42).
  - $p^\otimes$ from $\mathcal E^0=\sum_{k=1,2}\alpha_k^\otimes\mathcal E_k(p^\otimes,(\alpha_k\rho_k)^\otimes/\alpha_k^\otimes)+\alpha_3^\otimes\mathcal E_3(p^\otimes,(\alpha_3\rho_3)^0/\alpha_3^\otimes)$ (eq. 43), explicit for SG via eq. 10.
- **Unphysical volume fractions** (eq. 44), e.g. complete evaporation or condensation, or very small $\alpha$:
  - Set $\alpha_1^\otimes=\epsilon_\alpha$ if $s_2^{**}>s_1^{**}$, or $\alpha_2^\otimes=\epsilon_\alpha$ if $s_1^{**}>s_2^{**}$, with $\epsilon_\alpha=10^{-8}$.
  - Compute $\bar K_g=-\frac1{\Delta t}\ln\big(1-\frac{\bar\alpha_k-\alpha_k^{**}}{\hat{\mathcal S}^{**}_{\alpha_k}\Delta g^{**}/(1/\xi_1^{g**}+1/\xi_2^{g**})}\big)$, with $\hat{\mathcal S}_{\alpha_2}=-(\hat{\mathcal S}_{\alpha_1}+\hat{\mathcal S}_{\alpha_3})$.
  - Update $\alpha_3$ and $\alpha_1\rho_1$ with $\bar K_g$; the last $\alpha$ comes from $\sum\alpha=1$.
  - The results are insensitive to $\epsilon_\alpha\in[10^{-12},10^{-8}]$.
- **Cost.** The new method is explicit. Against the earlier equilibrium-iteration method it gave CPU
  100.5 s vs 147.0 s (pTg) and 100.8 s vs 157.6 s (pT(lv)g) on the cavitation tube, and 2635 s vs 3641 s
  on the injector. The $\nu\to\infty$ results match the iterative method to ~1e-6–1e-7 relative.

HLLC solver: same structure as [P22] Appendix C, extended to three phases (eqs. C.1–C.11). Both $\alpha_1$
and $\alpha_3$ are unchanged across the outer waves.

### 3.5 Test cases
**[P26] 5.1 / [P24] 4.1: three-phase water cavitation tube**
- Domain $[0,1]$; $\alpha_v=10^{-2}$, $\alpha_g=10^{-1}$ (air); $p_0=10^5$; densities at $T_0=354$ K (below $T_{sat}$).
- $u_0=\mp20$ m/s about $x=0.5$; SG water Table 3; 3000 cells, CFL 0.5, $t=6$ ms.
- Six settings: pT, pT(lv), pTg, pTg(f), pT(lv)g, pT(lv)g(f), with $\nu=0.4\times10^{-3}$ kg s/m⁵ for the (f) cases.
- Chemical relaxation only when $T_l>T_{sat}(p)$, i.e. evaporation only.
- Plotted: $p$, $u$, $Y_v$, $\alpha_v+\alpha_g$, phasic temperatures and chemical potentials.

**[P26] 5.2 / [P24] 4.2: UNDEX near a flat wall**
- 481×280 cells over $[-0.6,0.6]\times[-0.7,0]$; bubble centred at $(0,-0.22)$, $r=0.05$.
- Inside the bubble: $p=8290\times10^5$, $\rho_g=1400$, $\alpha_l=\alpha_v=10^{-8}$.
- Outside: $p=10^5$, $T=303$ K, $\alpha_v=10^{-4}$, $\alpha_g=10^{-7}$.
- Combustion gas is ideal with $\gamma=2$; water SG Table 3.
- Thermal relaxation liquid–vapor only (p T̃ procedure); $\nu\in\{0,10^{-4},\infty\}$.
- At $t=0.375$ ms the wall minimum pressure is ≈ 4420 Pa ($\nu\to\infty$, = $p_{sat}$), ≈ 1342 Pa ($\nu=10^{-4}$), ≈ 76 Pa ($\nu=0$).
- Pressure fields at 0, 0.075, 0.2, 0.325, 0.375, 0.5 ms.

**[P26] 5.3: fuel injector**
- Nozzle as in [P22] 8.2.5.
- Liquid dodecane at $p=10^8$, $T=600$ K, $\alpha_v=10^{-4}$, $\alpha_g=10^{-6}$.
- Chamber: air+vapor at $p=10^5$, $\rho_v=4$, $\alpha_g=0.9$, $\alpha_l=10^{-6}$ (so $T=511.13$ K $>T_{sat}=475$ K, $Y_v=0.3946182$).
- All phases in thermal equilibrium (pT procedure); SG dodecane Table 4; air $\gamma=1.4$, $\kappa_v=718$.
- Activation where $\min(\alpha_l,\alpha_v)>10^{-3}$; 400×160 cells, CFL 0.5.
- $\nu\in\{0,1,10,\infty\}$ kg s/m⁵; snapshots at 0.12, 0.3, 0.6 ms.
- The unphysical-$\alpha$ correction is triggered for 2–6 cells per step, between $t\approx0.64\times10^{-4}$ and $1.12\times10^{-4}$ s.

---

## 4. What carries over to NGA2 and what doesn't
- NGA2's gas phase is a Dalton mixture of vapor and air at one temperature. That is neither [P26] case
  (i) nor (ii): the air is not a separate volume, and the vapor's chemical potential is taken at its
  partial pressure. The template for NGA2 is [P22] (two-phase) plus the mixture-entropy stiffness
  $K_{CC}=R_vT(1-x_v)/m_v$ added to $1/\xi_1^g+1/\xi_2^g$.
- See `src/thermo/finite_rate_relaxation.md` for the NGA2 implementation status, the mapping
  $K_C=\Delta v\,\Pi-\Delta s\,\Theta+K_{CC}$ (with Pelanti's $S_p=\Pi$ and $S_T=\Theta$), and the plan for the
  chemical step.
- [P26]'s three-phase $\mathcal S$ formulas are useful if NGA2 ever models the inert gas as a separate phase.

## 5. Key modelling assumptions, collected
1. Rates are constant or frozen at the start of each relaxation step; each relaxation ODE is solved
   with all coefficients frozen ($\xi$, $Z/D$, $S_\alpha$).
2. The pressure is always recomputed from the mixture energy relation, never from the exponential;
   this guarantees mixture-energy consistency.
3. Thermal relaxation is posed at pressure equilibrium. Chemical relaxation is posed at pressure and
   (full or liquid–vapor) temperature equilibrium.
4. Activation criteria are problem-dependent: interface cells $\min\alpha>\epsilon$; cavitation tests
   evaporate only, with $T>T_{sat}$.
5. Pure phases are $\alpha=1-\epsilon$; the EOS of each phase must stay valid everywhere.

## 6. Accuracy and relation to the Onsager formulation
The finite-rate formulas (exponential decay of the potential differences) are first-order accurate in
the variation of the frozen coefficients over $\Delta t$. They are exact for linear relaxation.
`~/Downloads/finite_rate_relaxation.pdf` shows that all of them follow from an Onsager form
$\dot\xi=L F$ with $L=D\Lambda D^T$ and $\Lambda=\mathrm{diag}(\mu T,\vartheta T^2,\nu T)$.

## 7. Verification scripts (in `examples/eos_relax/`, run with `source ~/venv/combustion/bin/activate`)
- `check_KC_pelanti.py` — two-phase NASG water ([P22] Table 3). [P22] eq. 23 ($S_\alpha$, $S_T$, $S_p$)
  agrees with a finite difference of the pT equilibrium; the corrected $1/\xi_1^g+1/\xi_2^g$ equals the
  NGA2 `step_C` stiffness; the literal eq. 72 does not.
- `check_pelanti_3phase.py` — three-phase SG water + air ([P26] Table 3). [P26] eqs. 16–17 agree as
  printed; eqs. 24a and 25 agree as printed; eqs. 24b–d agree only in the corrected form (§8).

## 8. Errata (verified numerically where marked)
1. **[P22] eq. 72, [P26] eq. 39, [P24] eq. 9.** Printed: $1/\xi_k^g=\frac{(-1)^{k-1}}{\rho_k}S_p-s_kS_T$.
   Correct: $1/\xi_k^g=(-1)^{k-1}\big(\frac{S_p}{\rho_k}-s_kS_T\big)$. Read literally, the printed form makes
   $K_g$ contain $s_1+s_2$, which depends on the entropy reference constants $\tilde\eta$. **Verified.**
2. **[P22] eqs. 79–80.** They print $\Delta g^*$ and $\alpha_1^*$; these should be $\Delta g^{**}$ and $\alpha_1^{**}$.
3. **[P26] eq. 42a.** Prints $\hat{\mathcal S}^{**}_k$; this should be $\hat{\mathcal S}^{**}_{\alpha_k}$.
4. **[P26] eqs. 24b, 24c, 24d (and [P24] eqs. 8b–d).** The vapor-side bracket prints
   $\big(\frac{\chi_2}{\Gamma_2}+(c_1^2-\chi_1)\big)\frac{\rho_1}{\Gamma_1}$ (24b, 24c) and
   $\big(\frac{\chi_2\zeta_2-1}{\Gamma_2}+\zeta_2(c_1^2-\chi_1)\big)\frac{\rho_1}{\Gamma_1}$ (24d).
   Correct: $\big(\frac{\chi_2}{\Gamma_2}+\frac{c_1^2-\chi_1}{\Gamma_1}\big)\rho_1$ and
   $\big(\frac{\chi_2\zeta_2-1}{\Gamma_2}+\frac{\zeta_2(c_1^2-\chi_1)}{\Gamma_1}\big)\rho_1$, which is the $1\leftrightarrow2$ mirror
   of the liquid-side bracket. With the printed forms $\tilde{\mathcal S}_T$ is wrong by ~27%, and
   $\tilde{\mathcal S}_{\alpha_3}$, $\tilde{\mathcal S}_p$ by ~1e-4 relative, in the test state. **Verified.**
5. **[P24] eq. 8a.** Parenthesization puts the $\alpha_1^2c_3^2\phi_1\rho_3/\Gamma_1$ term inside the factor
   multiplied by $\alpha_1$, which would make it cubic in $\alpha_1$. Use [P26] eq. 24a (term outside the
   bracket). **The [P26] form is verified**; the [P24] form was not evaluated, but it breaks the
   $1\leftrightarrow2$ symmetry with the $\alpha_2^2\phi_2c_3^2\rho_3/\Gamma_2$ term.
6. **[P22] Sec. 8.2.1.** The pT(f) rate is printed in "Pa/(s K)", which equals W/(m³K); and $\nu$ in
   "Pa kg²/(s J²)", which equals kg s/m⁵. Same units as [P24]/[P26].
