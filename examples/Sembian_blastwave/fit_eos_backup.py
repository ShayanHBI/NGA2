"""
Stiffened-Gas (SG), Noble-Abel Stiffened-Gas (NASG), and Extended
Noble-Abel Stiffened-Gas (ENASG) EOS coefficient fitter + plotter.

Implements the liquid SG/NASG and ENASG methodology described in:
    Le Metayer, O. and Saurel, R., "The Noble-Abel Stiffened-Gas equation of state",
    Physics of Fluids 28, 046102 (2016).
    and
    Chiapolino, A. and Saurel, R., "Extended Noble-Abel Stiffened-Gas Equation
    of State for Sub-and-Supercritical Liquid-Gas Systems Far from the Critical
    Point", Fluids 3, 48 (2018).

SG caloric EOS:     p = (gamma - 1) * (e - q) / v        -  gamma * p_inf
NASG caloric EOS:   p = (gamma - 1) * (e - q) / (v - b)  -  gamma * p_inf
ENASG correction to NASG:
    b(v)      = b1*v + b0
    p_inf(T)  = p_inf1*T + p_inf0
    p0_inf(T) = gamma*p_inf1*T + gamma*p_inf0*(1-b1)/(gamma-b1)

Vapor phase is kept as an ideal gas.

Usage:
    python fit_eos.py [--mode sg|nasg|enasg|all] [--constraint soundspeed|hugoniot]

Dependencies: numpy, scipy, matplotlib, iapws  (pip install iapws)
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.optimize import brentq, minimize_scalar, least_squares as _ls
from iapws import IAPWS97

plt.rcParams.update({
    "text.usetex":        True,
    "font.family":        "serif",
    "font.size":          24,
    "axes.labelsize":     34,
    "xtick.labelsize":    28,
    "ytick.labelsize":    28,
    "legend.fontsize":    22,
    "figure.dpi":         150,
    "axes.linewidth":     2.2,
    "xtick.major.size":   7,
    "ytick.major.size":   7,
    "xtick.major.width":  1.5,
    "ytick.major.width":  1.5,
    "xtick.minor.size":   4,
    "ytick.minor.size":   4,
})


# =============================================================================
# Reference data
# =============================================================================

# For SG and NASG
WATER_REFERENCE_STATE = {
    "rho0": 999.9747,  # kg/m^3
    "p0":   1.0100e5,  # Pa
    "c0":   1421.53,   # m/s
}

# Chiapolino & Saurel (2018), Appendix C, Tables A1-A2.
# Units are SI.  The ENASG liquid fit needs more reference information than the
# NASG fit because b(v) and p_inf(T) are variable.
WATER_ENASG_REFERENCE_STATE = {
    "Tc":          646.16,       # K
    "pc":          221.0e5,      # Pa
    "vc":          0.0025101,    # m^3/kg
    "p0_inf_c":    0.01,         # Pa, prescribed p0_inf at the critical point
    "bc":          1.0e-6,       # m^3/kg, prescribed covolume at the critical point
    "c0":          1552.1,       # m/s
    "p0":          1.0e5,        # Pa
    "v0":          0.0010182,    # m^3/kg
    "T_ref":       300.16,       # K
    "p_ref":       3570.2,       # Pa
    "v_ref":       0.0010035,    # m^3/kg
    "e_ref":       113.23e3,     # J/kg
    "b_ref":       0.0009125,    # m^3/kg
}

# Water vapor constants for the ideal-gas phase used by ENASG.
# R = R_univ / W for H2O.
WATER_VAPOR_IG_REFERENCE_STATE = {
    "R":           461.522,        # J/(kg.K)
    "p0":          1.0e5,          # Pa
    "T_range":     (372.76, 1275.0),
    "T_ref":       393.38,         # K, Appendix C.2
    "e_ref":       2537.7e3,       # J/kg, Appendix C.2
    "cv":          1500.0,         # J/(kg.K), Table 1
}

# Water principal Hugoniot, Rice & Walsh 1957 (J. Chem. Phys. 26, 824):
# tabulated (up, us) [m/s], ambient state rho0 = 1/V0 with V0 = 1.0018 cc/g.
RICE_WALSH = {
    "up": 1.0e3 * np.array([0.952, 1.392, 1.411, 1.655, 1.798, 1.806, 1.829,
                             2.370, 2.385, 4.13, 4.60, 4.72, 4.72, 4.81]),
    "us": 1.0e3 * np.array([3.354, 4.093, 4.126, 4.536, 4.757, 4.777, 4.811,
                             5.601, 5.626, 8.07, 8.45, 8.49, 8.58, 8.74]),
    "rho0": 998.2,   # kg/m^3
    "p0":   1.0e5,   # Pa
}

# Fitting window and resolution configuration.
FIT_CONFIG = {
    "T_fit":      (300.0, 500.0),  # K, saturation curve fitting window
    "T_plot":     (250.0, 650.0),  # K, saturation curve plot window
    "n_fit":      201,             # saturation table points in fitting window
    "n_plot":     374,             # saturation table points in plot window
    "n_vapor":    542,             # isobar table points for ENASG vapor fit
    "n_curve":    300,             # EOS curve points for saturation plots
    "n_hugoniot": 200,             # EOS curve points for Hugoniot plot
}

# NASG hugoniot isobar fit configuration.
NASG_HUGONIOT_CONFIG = {
    "T_pre":        [274., 290., 305., 320., 335., 348.],  # K, pre-shock isotherms
    "p_isobar_min": 1.0e5,                                 # Pa
    "p_isobar_max": 300.0e5,                               # Pa
    "n_isobar":     20,                                    # pressure points per isotherm
    "rho2_max":     2000.0,                                # kg/m^3, maximum post-shock density (caps b)
    "x0":           [1.4, 8.0e8, 4.0e-4, 3500.0],         # initial guess [gamma, p_inf, b, cv]
    "gamma_bounds": (1.01, 8.0),
    "p_inf_bounds": (1.0e5, 5.0e10),
    "cv_bounds":    (50.0, 1.0e4),
}


def iapws_saturation_table(T_min, T_max, n_points=101):
    """
    Build a saturation-curve table from IAPWS-IF97.

    Returns
    -------
    ndarray, shape (N, 12) with columns
        T, p_sat, h_l, h_g, v_l, v_g, e_l, e_g, s_l, s_g, w_l, w_g
    in SI units (w_l, w_g are the saturated-liquid/vapor sound speeds, m/s).

    Rows where IAPWS97 raises (e.g. above the critical point) are skipped.
    """
    T = np.linspace(T_min, T_max, n_points)
    rows = []
    for Ti in T:
        try:
            liq = IAPWS97(T=Ti, x=0)
            vap = IAPWS97(T=Ti, x=1)
            rows.append([Ti,
                         liq.P * 1e6,
                         liq.h * 1e3, vap.h * 1e3,
                         liq.v, vap.v,
                         liq.u * 1e3, vap.u * 1e3,
                         liq.s * 1e3, vap.s * 1e3,
                         liq.w, vap.w])
        except Exception:
            pass
    return np.array(rows)


def iapws_isobar_table(p, T_min, T_max, n_points=400):
    """
    Build an isobar table from IAPWS-IF97.

    Parameters
    ----------
    p : pressure in Pa

    Returns
    -------
    ndarray, shape (N, 5) with columns T, p, h, e, s in SI units.
    """
    T = np.linspace(T_min, T_max, n_points)
    rows = []
    for Ti in T:
        try:
            st = IAPWS97(P=p * 1e-6, T=Ti)
            rows.append([Ti, p, st.h * 1e3, st.u * 1e3, st.s * 1e3])
        except Exception:
            pass
    return np.array(rows)


# =============================================================================
# Vapor phase fit for SG/NASG: saturation-based ideal gas
# =============================================================================

def fit_vapor_phase(T, p, hg, vg):
    """
    Vapor-phase coefficients with b_g = 0 and p_inf_g = 0 (ideal-gas vapor).
    Implements Eqs. (50), (51), (55) of Le Metayer & Saurel (2016).
    """
    T_mean  = T.mean()
    hg_mean = hg.mean()

    # h_g(T) = cp_g * T + q_g  (least squares)
    cp_g = np.sum(T * (hg - hg_mean)) / np.sum(T * (T - T_mean))
    q_g  = hg_mean - cp_g * T_mean

    # v_g(T) = (cp_g - cv_g) * T / p_sat(T)  (least squares)
    r_g  = np.sum(vg * T / p) / np.sum((T / p) ** 2)
    cv_g = cp_g - r_g

    return {
        "eos":     "IG",
        "cp":      cp_g,
        "cv":      cv_g,
        "gamma":   cp_g / cv_g,
        "q":       q_g,
        "q_prime": np.nan,
        "b":       0.0,
        "b0":      0.0,
        "b1":      0.0,
        "p_inf":   0.0,
        "p_inf0":  0.0,
        "p_inf1":  0.0,
    }


# =============================================================================
# NASG liquid phase fit (b_l solved jointly with cp_l - cv_l from v_l data)
# =============================================================================

def _fit_liquid_nasg_for_given_pinf(T, p, hl, vl, p_inf_l):
    """
    Given a candidate p_inf,l, solve Eqs. (60)-(61), (64)-(65) for
    (cp_l, cv_l, q_l, b_l).  Both b_l and cp_l - cv_l are free.
    """
    T_mean  = T.mean()
    p_mean  = p.mean()
    hl_mean = hl.mean()
    vl_mean = vl.mean()

    Tp      = T / (p + p_inf_l)
    Tp_mean = Tp.mean()

    # Two-unknown least squares on v_l = (cp_l-cv_l)*T/(p+p_inf) + b_l
    r_l = np.sum(Tp * (vl - vl_mean)) / np.sum(Tp * (Tp - Tp_mean))
    b_l = vl_mean - r_l * Tp_mean

    # cp_l and q_l from h_l = cp_l*T + b_l*p + q_l
    num_cp = np.sum(T * (hl - hl_mean)) - b_l * np.sum(T * (p - p_mean))
    den_cp = np.sum(T * (T - T_mean))
    cp_l = num_cp / den_cp
    q_l  = hl_mean - cp_l * T_mean - b_l * p_mean

    return cp_l, cp_l - r_l, q_l, b_l


def _fit_nasg_hugoniot_isobar(T_sat, p_sat, hl_sat, hg_sat=None, vapor=None,
                               config=None):
    """
    Fit NASG liquid (gamma, p_inf, b, cv) jointly against IAPWS97 isobar data
    for rho(T,p) and c(T,p) plus, if vapor is provided, saturation p_sat shape
    and latent heat residuals.

    b is capped at 1/rho2_max so that rho2_max is the maximum post-shock density.
    q is derived from saturation enthalpy LS with the isobar-fitted (cp, b).
    """
    if config is None:
        config = NASG_HUGONIOT_CONFIG

    T_pre = np.array(config["T_pre"])
    p_arr = np.linspace(config["p_isobar_min"], config["p_isobar_max"],
                        config["n_isobar"])
    b_max = 1.0 / config["rho2_max"]

    t_list, p_list, rho_list, c_list = [], [], [], []
    for Ti in T_pre:
        for pi in p_arr:
            try:
                st = IAPWS97(P=pi * 1e-6, T=Ti)
                if (np.isfinite(st.rho) and np.isfinite(st.w)
                        and st.rho > 800.0 and st.w > 100.0):
                    t_list.append(Ti);     p_list.append(pi)
                    rho_list.append(st.rho);  c_list.append(st.w)
            except Exception:
                pass

    T_d   = np.array(t_list);    p_d   = np.array(p_list)
    rho_d = np.array(rho_list);  c_d   = np.array(c_list)
    if len(T_d) == 0:
        raise RuntimeError("No valid IAPWS97 isobar data")

    L_iapws = None
    if vapor is not None and hg_sat is not None:
        L_arr = hg_sat - hl_sat
        ok_L  = np.isfinite(L_arr) & (L_arr > 0.0)
        if ok_L.any():
            L_iapws  = L_arr[ok_L]
            T_sat_L  = T_sat[ok_L]
            p_sat_L  = p_sat[ok_L]

    def _rho(p, T, g, pi, b, cv):
        return (p + pi) / ((g - 1.0)*cv*T + b*(p + pi))

    def _c(p, T, g, pi, b, cv):
        rho = _rho(p, T, g, pi, b, cv)
        return np.sqrt(np.maximum(0.0, g*(p + pi) / (rho*(1.0 - b*rho))))

    def _q_ls(g, b, cv):
        cp_l = g * cv
        T_m  = T_sat.mean()
        return (hl_sat - b * p_sat).mean() - cp_l * T_m

    def res(x):
        g, pi, b, cv = x
        rr = (_rho(p_d, T_d, g, pi, b, cv) - rho_d) / rho_d
        rc = (_c(p_d,   T_d, g, pi, b, cv) - c_d  ) / c_d
        blocks = [rr, rc]

        if vapor is not None and L_iapws is not None:
            cp_v = vapor["cp"];  cv_v = vapor["cv"]
            q_v  = vapor["q"];   r_v  = cp_v - cv_v
            cp_l = g * cv
            q_l  = _q_ls(g, b, cv)
            # Latent heat: L = (cp_v - cp_l)*T + q_v - q_l - b*p_sat
            L_model = (cp_v - cp_l)*T_sat_L + q_v - q_l - b*p_sat_L
            rl = (L_model - L_iapws) / np.abs(L_iapws)
            # Saturation curve shape
            ds  = (cp_l - cv) / r_v
            es  = b / r_v
            bs  = (q_l - q_v) / r_v
            cs  = (cp_v - cp_l) / r_v
            lhs = (np.log(p_sat_L) - (bs + es*p_sat_L)/T_sat_L
                   - cs*np.log(T_sat_L) - ds*np.log(p_sat_L + pi))
            rp = lhs - lhs.mean()
            blocks.extend([rl, rp])

        return np.concatenate(blocks)

    g_lo,  g_hi  = config["gamma_bounds"]
    pi_lo, pi_hi = config["p_inf_bounds"]
    cv_lo, cv_hi = config["cv_bounds"]
    x0 = list(config["x0"])
    lo = [g_lo,  pi_lo, 1.0e-7, cv_lo]
    hi = [g_hi,  pi_hi, b_max,  cv_hi]
    fit = _ls(res, x0, bounds=(lo, hi),
              x_scale=[0.5, 1.0e9, 1.0e-4, 1000.0],
              xtol=1.0e-12, ftol=1.0e-12, max_nfev=30000)
    g, pi, b, cv = fit.x
    cp_l = g * cv
    q    = _q_ls(g, b, cv)

    return {
        "eos":     "NASG",
        "cp":      cp_l,
        "cv":      cv,
        "gamma":   g,
        "q":       q,
        "q_prime": np.nan,
        "b":       b,
        "b0":      b,
        "b1":      0.0,
        "p_inf":   pi,
        "p_inf0":  pi,
        "p_inf1":  0.0,
    }


def fit_liquid_phase_nasg(T, p, hl, vl, ref_state, p_inf_bounds=(1.0e7, 5.0e9),
                          constraint="soundspeed", hg=None, vapor=None):
    """
    Close the NASG liquid-phase system by fixing p_inf,l with one of:

    constraint="soundspeed" (default): the ambient reference sound speed,
    Le Metayer & Saurel (2016) Eq. (68). Brent's method finds p_inf,l such that:
        p0 + p_inf,l - (cv_l / cp_l) * rho0 * c0^2 * (1 - b_l * rho0) = 0

    constraint="hugoniot": least-squares match of the NASG principal Hugoniot
    to the Rice & Walsh (1957) shock data (see hugoniot_sse).
    """
    rho0 = ref_state["rho0"]
    p0   = ref_state["p0"]
    c0   = ref_state["c0"]

    def build_eos(p_inf_l):
        cp_l, cv_l, q_l, b_l = _fit_liquid_nasg_for_given_pinf(T, p, hl, vl, p_inf_l)
        return {
            "eos":     "NASG",
            "cp":      cp_l,
            "cv":      cv_l,
            "gamma":   cp_l / cv_l,
            "q":       q_l,
            "q_prime": np.nan,
            "b":       b_l,
            "b0":      b_l,
            "b1":      0.0,
            "p_inf":   p_inf_l,
            "p_inf0":  p_inf_l,
            "p_inf1":  0.0,
        }

    def solve_soundspeed():
        def residual(p_inf_l):
            cp_l, cv_l, _, b_l = _fit_liquid_nasg_for_given_pinf(T, p, hl, vl, p_inf_l)
            return p0 + p_inf_l - (cv_l / cp_l) * rho0 * c0**2 * (1.0 - b_l * rho0)

        a, b = p_inf_bounds
        fa, fb = residual(a), residual(b)
        if fa * fb > 0:
            for _ in range(8):
                b *= 5.0
                fb = residual(b)
                if fa * fb < 0:
                    break
            else:
                raise RuntimeError(
                    f"NASG: could not bracket root for p_inf_l. "
                    f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
                )
        return brentq(residual, a, b, xtol=1.0, rtol=1e-10)

    if constraint == "hugoniot":
        try:
            return _fit_nasg_hugoniot_isobar(T, p, hl, hg_sat=hg, vapor=vapor)
        except Exception as exc:
            print(f"  [NASG] Hugoniot (isobar) fit failed ({exc}); "
                  "falling back to soundspeed closure.")

    p_inf_l = solve_soundspeed()
    return build_eos(p_inf_l)


# =============================================================================
# SG liquid phase fit (b_l = 0 forced throughout)
# =============================================================================

def _fit_liquid_sg_for_given_pinf(T, p, hl, vl, p_inf_l):
    """
    Given a candidate p_inf,l, solve for (cp_l, cv_l, q_l) with b_l = 0.

    For SG:
      h_l(T) = cp_l * T + q_l             [b=0, no p term]
      v_l(T) = (cp_l - cv_l)*T/(p+p_inf)  [b=0, one-parameter fit]
    """
    T_mean  = T.mean()
    hl_mean = hl.mean()

    # cp_l and q_l from h_l = cp_l*T + q_l (independent of p_inf_l)
    cp_l = np.sum(T * (hl - hl_mean)) / np.sum(T * (T - T_mean))
    q_l  = hl_mean - cp_l * T_mean

    # cp_l - cv_l from v_l = r_l * T/(p+p_inf_l), one-parameter least squares
    Tp   = T / (p + p_inf_l)
    r_l  = np.sum(vl * Tp) / np.sum(Tp ** 2)
    cv_l = cp_l - r_l

    return cp_l, cv_l, q_l, 0.0


def fit_liquid_phase_sg(T, p, hl, vl, ref_state, p_inf_bounds=(1.0e7, 5.0e9),
                        constraint="soundspeed"):
    """
    Close the SG liquid-phase system by fixing p_inf,l with one of:

    constraint="soundspeed" (default): the ambient reference sound speed. With
    b_l = 0 the closure (Eq. 68) simplifies to, solved by Brent's method:
        p0 + p_inf,l - (cv_l / cp_l) * rho0 * c0^2 = 0

    constraint="hugoniot": least-squares match of the SG principal Hugoniot
    (NASG with b_l = 0) to the Rice & Walsh (1957) shock data.
    """
    rho0 = ref_state["rho0"]
    p0   = ref_state["p0"]
    c0   = ref_state["c0"]

    def build_eos(p_inf_l):
        cp_l, cv_l, q_l, b_l = _fit_liquid_sg_for_given_pinf(T, p, hl, vl, p_inf_l)
        return {
            "eos":     "SG",
            "cp":      cp_l,
            "cv":      cv_l,
            "gamma":   cp_l / cv_l,
            "q":       q_l,
            "q_prime": np.nan,
            "b":       b_l,
            "b0":      b_l,
            "b1":      0.0,
            "p_inf":   p_inf_l,
            "p_inf0":  p_inf_l,
            "p_inf1":  0.0,
        }

    def solve_soundspeed():
        def residual(p_inf_l):
            cp_l, cv_l, _, _ = _fit_liquid_sg_for_given_pinf(T, p, hl, vl, p_inf_l)
            return p0 + p_inf_l - (cv_l / cp_l) * rho0 * c0**2

        a, b = p_inf_bounds
        fa, fb = residual(a), residual(b)
        if fa * fb > 0:
            for _ in range(8):
                b *= 5.0
                fb = residual(b)
                if fa * fb < 0:
                    break
            else:
                raise RuntimeError(
                    f"SG: could not bracket root for p_inf_l. "
                    f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
                )
        return brentq(residual, a, b, xtol=1.0, rtol=1e-10)

    if constraint == "hugoniot":
        p_inf_l, feasible = _minimize_hugoniot(lambda pi: hugoniot_objective(build_eos(pi)),
                                                p_inf_bounds, xatol=1.0)
        if not feasible:
            print("  [SG] Hugoniot constraint infeasible (EOS admits no Rice & Walsh "
                  "shock states for any p_inf,l in bounds); falling back to soundspeed closure.")
            p_inf_l = solve_soundspeed()
    else:
        p_inf_l = solve_soundspeed()

    return build_eos(p_inf_l)


# =============================================================================
# ENASG liquid phase fit (Chiapolino & Saurel 2018, Appendix C.1)
# =============================================================================

def enasg_pinf(T, eos):
    """p_inf(T) = p_inf1*T + p_inf0."""
    return eos["p_inf1"] * T + eos["p_inf0"]


def enasg_p0inf(T, eos):
    """p0_inf(T) = gamma*p_inf1*T + gamma*p_inf0*(1-b1)/(gamma-b1)."""
    gamma = eos["gamma"]
    b1    = eos["b1"]
    return gamma * eos["p_inf1"] * T + gamma * eos["p_inf0"] * (1.0 - b1) / (gamma - b1)


def enasg_b_of_v(v, eos):
    """b(v) = b1*v + b0."""
    return eos["b1"] * v + eos["b0"]


def _enasg_auxiliary(A, C, b1, p, T):
    """
    D, E, F of Chiapolino & Saurel (2018), Appendix C, Eq. (A24).

    A = gamma*p_inf1 and C = gamma*p_inf0*(1-b1)/(gamma-b1).
    """
    denom = (1.0 - b1) * (p + A * T + C)
    D = C * (p + C) / denom
    E = -C * b1 / (1.0 - b1) + A * C * T / denom
    F = (p + A * T + C) / T - A
    return D, E, F


def _enasg_gamma_cv_for_given_C(T, p, v, e, ref_state, C, B_fixed=None):
    """
    Given C = gamma*p_inf0*(1-b1)/(gamma-b1), compute gamma and cv from
    the volume and internal-energy least-squares relations, Eqs. (A23)-(A29).

    If B_fixed is None, B = b0/(1-b1) is treated as a free parameter and
    fitted jointly with K = (gamma-1)*cv from the two-parameter volume LS
    (v = K*phi + B).  If B_fixed is provided it is used directly (one-param LS).

    Returns (gamma, cv, A, B).
    """
    Tc       = ref_state["Tc"]
    p0_inf_c = ref_state["p0_inf_c"]
    b1       = ref_state["b1"]
    A        = (p0_inf_c - C) / Tc

    denom = p + A * T + C
    if np.any(denom <= 0.0):
        raise ValueError("ENASG: p+p0_inf(T) became non-positive during the fit.")

    phi = T / ((1.0 - b1) * denom)

    if B_fixed is None:
        # Two-parameter LS: v = K*phi + B (both K and B free)
        A_ls = np.column_stack([phi, np.ones_like(phi)])
        K, B = np.linalg.lstsq(A_ls, v, rcond=None)[0]
    else:
        B = B_fixed
        K = np.sum((v - B) * phi) / np.sum(phi**2)

    D, E, F = _enasg_auxiliary(A, C, b1, p, T)
    Dr, Er, Fr = _enasg_auxiliary(A, C, b1, ref_state["p_ref"], ref_state["T_ref"])

    G  = (p + E) / F
    H  = D / F
    Gr = (ref_state["p_ref"] + Er) / Fr
    Hr = Dr / Fr

    de = e - ref_state["e_ref"]
    dG = G - Gr
    dH = H - Hr

    Se1 = np.sum(de * dG)
    Se2 = np.sum(dG**2)
    Se3 = np.sum(de * dH)
    Se4 = 2.0 * np.sum(dG * dH)
    Se5 = np.sum(dH**2)

    # A27 + A29 give a quadratic in gamma after eliminating cv.
    # cv = K/(gamma-1).
    qa = Se3 - K * Se5
    qb = Se1 - Se3 - K * Se4
    qc = -Se1 - K * Se2

    roots = np.roots([qa, qb, qc]) if abs(qa) > 1e-300 else np.roots([qb, qc])
    candidates = []
    for root in roots:
        if abs(root.imag) > 1e-7 * max(1.0, abs(root.real)):
            continue
        gamma = root.real
        if gamma <= 1.0:
            continue
        cv = K / (gamma - 1.0)
        if cv <= 0.0 or not np.isfinite(cv):
            continue
        # Convexity condition expects p_inf1 <= 0, p_inf0 >= 0, b1 < 1.
        p_inf1 = A / gamma
        p_inf0 = C * (gamma - b1) / (gamma * (1.0 - b1))
        if p_inf1 <= 0.0 and p_inf0 >= 0.0 and b1 < 1.0:
            candidates.append((gamma, cv))

    if not candidates:
        raise RuntimeError("ENASG: no physical gamma/cv root found for this C.")

    # The ENASG liquid root is the one just above unity for water/oxygen.
    gamma, cv = sorted(candidates, key=lambda pair: pair[0])[0]
    return gamma, cv, A, B


def enasg_c2_from_p_v(p, v, eos):
    """ENASG sound speed, Chiapolino & Saurel (2018), Eq. (42)."""
    gamma = eos["gamma"]
    cv    = eos["cv"]
    b1    = eos["b1"]
    p0bar = gamma * eos["p_inf0"] * (1.0 - b1) / (gamma - b1)
    A     = gamma * eos["p_inf1"]
    vb    = v - enasg_b_of_v(v, eos)
    den   = (gamma - 1.0) * cv - A * vb

    return (
        -v**2 * A * (p + p0bar) * ((gamma - 1.0) / den + 1.0 / cv)
        + ((p + p0bar) / den)
        * (v**2 * (gamma - b1) * (gamma - 1.0) * cv / vb)
    )


def fit_liquid_phase_enasg(T, p, e, v, s, ref_state, C_bounds=(1.0e7, 1.0e9),
                           constraint="soundspeed"):
    """
    Fit the liquid ENASG coefficients with Chiapolino & Saurel (2018),
    Appendix C.1.  gamma and cv come from the coupled least-squares
    volume/internal-energy relations for any candidate C; q comes from the
    reference internal energy and q_prime from the entropy fit.

    b0 and b1 are treated as free parameters:
      - b1 is fixed from the prescribed critical-point slope
        (b_ref + b1*(vc - v_ref) = bc).
      - For soundspeed closure, B = b0/(1-b1) is fitted jointly with
        K = (gamma-1)*cv from the two-parameter volume LS (v = K*phi + B).
      - For hugoniot closure, B is optimized over a 1-D search; for each
        candidate B the sound-speed condition determines C, and the best B
        is the one minimising the Rice & Walsh Hugoniot residual.

    The unknown C is closed with one of:

    constraint="soundspeed" (default): the reference sound-speed condition,
    Eq. (42), solved by Brent's method (B from two-parameter LS).

    constraint="hugoniot": 1-D search over B with C tied to the sound-speed
    condition for each B.  Minimises the Rice & Walsh (1957) Hugoniot residual.
    Falls back to soundspeed closure if no admissible B is found.
    """
    ref_state = dict(ref_state)
    ref_state["b1"] = (ref_state["bc"] - ref_state["b_ref"]) / (ref_state["vc"] - ref_state["v_ref"])

    b1 = ref_state["b1"]

    def build_eos_for_C(C, B_fixed=None):
        gamma, cv, A, B = _enasg_gamma_cv_for_given_C(T, p, v, e, ref_state, C, B_fixed=B_fixed)
        b0     = B * (1.0 - b1)
        p_inf1 = A / gamma
        p_inf0 = C * (gamma - b1) / (gamma * (1.0 - b1))

        # Reference energy q, Eq. (A25).  Computed here (not just after C is
        # fixed) so that eos_h is usable for any candidate C, e.g. inside the
        # Hugoniot objective.
        D_ref, E_ref, F_ref = _enasg_auxiliary(A, C, b1, ref_state["p_ref"], ref_state["T_ref"])
        q = ref_state["e_ref"] - cv * (ref_state["p_ref"] + gamma * D_ref + E_ref) / F_ref

        return {
            "eos":     "ENASG",
            "cp":      np.nan,   # cp is not constant for ENASG liquid.
            "cv":      cv,
            "gamma":   gamma,
            "q":       q,
            "q_prime": np.nan,
            "b":       np.nan,
            "b0":      b0,
            "b1":      b1,
            "p_inf":   np.nan,
            "p_inf0":  p_inf0,
            "p_inf1":  p_inf1,
        }

    def solve_soundspeed(B_fixed=None):
        def residual(C):
            eos = build_eos_for_C(C, B_fixed=B_fixed)
            return enasg_c2_from_p_v(ref_state["p0"], ref_state["v0"], eos) - ref_state["c0"]**2

        a, b = C_bounds
        fa, fb = residual(a), residual(b)
        if fa * fb > 0.0:
            # Fall back to a log scan.  This is more robust when a user changes the
            # fitting range or reference states.
            grid = np.logspace(5, 11, 300)
            vals = []
            for Ci in grid:
                try:
                    vals.append(residual(Ci))
                except Exception:
                    vals.append(np.nan)
            vals = np.array(vals)
            idx = None
            for i in range(len(grid) - 1):
                if not np.isfinite(vals[i]) or not np.isfinite(vals[i + 1]):
                    continue
                if vals[i] == 0.0 or vals[i] * vals[i + 1] < 0.0:
                    idx = i
                    break
            if idx is None:
                raise RuntimeError(
                    "ENASG: could not bracket C from the sound-speed condition. "
                    f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
                )
            a, b = grid[idx], grid[idx + 1]

        return brentq(residual, a, b, xtol=1.0e-2, rtol=1.0e-11)

    if constraint == "hugoniot":
        # Optimise B (and hence rho_max = 1/B) with C tied to the sound-speed
        # condition.  This keeps the thermal EOS consistent with the reference
        # sound speed at every B trial while letting the density ceiling move.
        # Upper bound: B must be below the value that makes the softest
        # Rice & Walsh point admissible (rho_max = 1/B > rho2_min).
        rho2_rw = RICE_WALSH["rho0"] * RICE_WALSH["us"] / (RICE_WALSH["us"] - RICE_WALSH["up"])
        rho2_min_rw = rho2_rw.min()
        B_hi = (1.0 / rho2_min_rw) * 0.9999   # just under the threshold
        B_lo = 1.0e-5

        def objective_B(B_val):
            try:
                C_val = solve_soundspeed(B_fixed=B_val)
                return hugoniot_objective(build_eos_for_C(C_val, B_fixed=B_val))
            except Exception:
                return np.inf, 0

        B_opt, feasible = _minimize_hugoniot(objective_B, (B_lo, B_hi), xatol=1.0e-8)
        if feasible:
            C = solve_soundspeed(B_fixed=B_opt)
            liquid = build_eos_for_C(C, B_fixed=B_opt)
        else:
            print("  [ENASG] Hugoniot constraint infeasible (no admissible B found); "
                  "falling back to soundspeed closure.")
            C = solve_soundspeed()
            liquid = build_eos_for_C(C)
    else:
        C = solve_soundspeed()
        liquid = build_eos_for_C(C)

    # Entropy reference q_prime, Eq. (A31).  The saturation liquid entropy is used.
    p0inf = enasg_p0inf(T, liquid)
    liquid["q_prime"] = np.mean(
        s
        - liquid["cv"] * (
            ((liquid["gamma"] - liquid["b1"]) / (1.0 - liquid["b1"])) * np.log(T)
            - ((liquid["gamma"] - 1.0) / (1.0 - liquid["b1"])) * np.log(p + p0inf)
        )
        + liquid["gamma"] * liquid["p_inf1"] * (liquid["gamma"] - 1.0)
        * liquid["cv"] * T / ((1.0 - liquid["b1"]) * (p + p0inf))
    )

    return liquid


def fit_vapor_phase_ideal_isobar(gas_data, gas_ref):
    """
    Fit the ideal-gas vapor coefficients used with ENASG.

    The 2018 paper chooses the gas coefficients away from saturation, on an
    isobar, and enforces Mayer's relation cp-cv=R.  Here cv and q are obtained
    from a least-squares e(T)=cv*T+q fit on that isobar, gamma=(cv+R)/cv,
    and q_prime follows from the ideal-gas entropy formula.
    """
    T = gas_data[:, 0]
    p = gas_data[:, 1]
    e = gas_data[:, 3]
    s = gas_data[:, 4]
    R = gas_ref["R"]

    if "cv" in gas_ref:
        # Chiapolino & Saurel choose a constant representative cv for the gas
        # phase rather than fitting it on the saturation curve.
        cv = gas_ref["cv"]
    else:
        T_mean = T.mean()
        e_mean = e.mean()
        cv = np.sum(T * (e - e_mean)) / np.sum(T * (T - T_mean))

    if "T_ref" in gas_ref and "e_ref" in gas_ref:
        q = gas_ref["e_ref"] - cv * gas_ref["T_ref"]
    else:
        e_mean = e.mean()
        q = e_mean - cv * T.mean()

    gamma = (cv + R) / cv
    cp    = cv + R

    q_prime = np.mean(s - cv * (gamma * np.log(T) - (gamma - 1.0) * np.log(p)))

    return {
        "eos":     "IG",
        "cp":      cp,
        "cv":      cv,
        "gamma":   gamma,
        "q":       q,
        "q_prime": q_prime,
        "b":       0.0,
        "b0":      0.0,
        "b1":      0.0,
        "p_inf":   0.0,
        "p_inf0":  0.0,
        "p_inf1":  0.0,
    }


# =============================================================================
# Entropy constants for SG/NASG
# Eq. (41)/(69-72).  E = (b_l - b_g)/Rg -> 0 automatically for SG.
# =============================================================================

def fit_entropy_constants(T, p, liquid, vapor):
    """
    Determine q'_l and q'_g via the saturated-vapor-pressure relation.
    Convention: q'_l = 0.
    For SG, b_l = b_g = 0 so E = 0 and the formula reduces to the SG form.
    """
    cp_l, cv_l, p_inf_l = liquid["cp"], liquid["cv"], liquid["p_inf"]
    cp_g, cv_g          = vapor["cp"],  vapor["cv"]
    b_l, b_g            = liquid["b"],  vapor["b"]
    q_l, q_g            = liquid["q"],  vapor["q"]

    r_g = cp_g - cv_g
    B   = (q_l - q_g) / r_g
    C   = (cp_g - cp_l) / r_g
    D   = (cp_l - cv_l) / r_g
    E   = (b_l - b_g) / r_g          # zero for SG

    A = np.mean(
        np.log(p)
        - (B + E * p) / T
        - C * np.log(T)
        - D * np.log(p + p_inf_l)
    )

    q_prime_l = 0.0
    q_prime_g = A * r_g - (cp_l - cp_g) + q_prime_l
    return q_prime_l, q_prime_g


# =============================================================================
# Top-level fitters
# =============================================================================

def compute_nasg_coefficients(sat_data, ref_state, T_range, constraint="soundspeed"):
    """Full NASG fit: vapor (ideal-gas) + liquid (with covolume b_l)."""
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 4:
        raise ValueError(f"Need >=4 points in [{Tmin},{Tmax}] K, got {len(d)}.")

    T, p, hl, hg, vl, vg = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

    vapor  = fit_vapor_phase(T, p, hg, vg)
    liquid = fit_liquid_phase_nasg(T, p, hl, vl, ref_state, constraint=constraint,
                                   hg=hg, vapor=vapor)
    qp_l, qp_g = fit_entropy_constants(T, p, liquid, vapor)
    liquid["q_prime"] = qp_l
    vapor["q_prime"]  = qp_g
    return liquid, vapor


def compute_sg_coefficients(sat_data, ref_state, T_range, constraint="soundspeed"):
    """Full SG fit: vapor (ideal-gas) + liquid (b_l = 0 forced)."""
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 4:
        raise ValueError(f"Need >=4 points in [{Tmin},{Tmax}] K, got {len(d)}.")

    T, p, hl, hg, vl, vg = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

    vapor  = fit_vapor_phase(T, p, hg, vg)
    liquid = fit_liquid_phase_sg(T, p, hl, vl, ref_state, constraint=constraint)
    qp_l, qp_g = fit_entropy_constants(T, p, liquid, vapor)
    liquid["q_prime"] = qp_l
    vapor["q_prime"]  = qp_g
    return liquid, vapor


def compute_enasg_coefficients(sat_data, liquid_ref, gas_ref, T_range, constraint="soundspeed"):
    """Full ENASG fit: liquid ENASG + vapor ideal gas."""
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 8:
        raise ValueError(f"Need >=8 points in [{Tmin},{Tmax}] K, got {len(d)}.")

    T, p, vl, el, sl = d[:, 0], d[:, 1], d[:, 4], d[:, 6], d[:, 8]

    gas_data = iapws_isobar_table(
        gas_ref["p0"], gas_ref["T_range"][0], gas_ref["T_range"][1],
        n_points=FIT_CONFIG["n_vapor"]
    )
    vapor  = fit_vapor_phase_ideal_isobar(gas_data, gas_ref)
    liquid = fit_liquid_phase_enasg(T, p, el, vl, sl, liquid_ref, constraint=constraint)
    return liquid, vapor


# =============================================================================
# Analytical saturation curves
# =============================================================================

def eos_v(T, p, phase):
    """Specific volume v(p,T) for SG/NASG/IG and ENASG liquid."""
    if phase["eos"] == "ENASG":
        return ((phase["gamma"] - 1.0) * phase["cv"] * T
                / ((1.0 - phase["b1"]) * (p + enasg_p0inf(T, phase)))
                + phase["b0"] / (1.0 - phase["b1"]))
    return (phase["cp"] - phase["cv"]) * T / (p + phase["p_inf"]) + phase["b"]


def eos_h(T, p, phase):
    """Specific enthalpy h(p,T) for SG/NASG/IG and ENASG liquid."""
    if phase["eos"] == "ENASG":
        gamma = phase["gamma"]
        b1    = phase["b1"]
        pinf  = enasg_pinf(T, phase)
        p0inf = enasg_p0inf(T, phase)
        return (
            phase["cv"] * T / (p + p0inf)
            * (gamma * (p + pinf) - p * b1 - gamma * b1 * pinf)
            / (1.0 - b1)
            + p * phase["b0"] / (1.0 - b1)
            + phase["q"]
        )
    return phase["cp"] * T + phase["b"] * p + phase["q"]


def eos_s(T, p, phase):
    """Specific entropy s(p,T) for SG/NASG/IG and ENASG liquid."""
    if phase["eos"] == "ENASG":
        gamma = phase["gamma"]
        b1    = phase["b1"]
        p0inf = enasg_p0inf(T, phase)
        return (
            phase["cv"] * (
                ((gamma - b1) / (1.0 - b1)) * np.log(T)
                - ((gamma - 1.0) / (1.0 - b1)) * np.log(p + p0inf)
            )
            - gamma * phase["p_inf1"] * (gamma - 1.0) * phase["cv"] * T
            / ((1.0 - b1) * (p + p0inf))
            + phase["q_prime"]
        )
    return phase["cv"] * (phase["gamma"] * np.log(T)
                          - (phase["gamma"] - 1.0) * np.log(p + phase["p_inf"])) + phase["q_prime"]


def eos_g(T, p, phase):
    """Gibbs free energy g(p,T)."""
    return eos_h(T, p, phase) - T * eos_s(T, p, phase)


def eos_psat(T, liquid, vapor):
    """
    Saturation pressure from Gibbs free-energy equality.

    For SG/NASG this uses the closed Le Metayer & Saurel saturation relation.
    For ENASG it solves g_l(p,T)=g_v(p,T) directly, which is equivalent to
    Chiapolino & Saurel's Eq. (52) but avoids hard-coding the long expression.
    """
    if liquid["eos"] == "ENASG" or vapor["eos"] == "ENASG":
        return eos_psat_gibbs(T, liquid, vapor)

    cp_l, cv_l, p_inf_l = liquid["cp"], liquid["cv"], liquid["p_inf"]
    cp_g, cv_g, p_inf_g = vapor["cp"],  vapor["cv"],  vapor["p_inf"]
    b_l, b_g             = liquid["b"],  vapor["b"]
    q_l, q_g             = liquid["q"],  vapor["q"]
    qp_l, qp_g           = liquid["q_prime"], vapor["q_prime"]

    r_g = cp_g - cv_g
    A   = (cp_l - cp_g + qp_g - qp_l) / r_g
    B   = (q_l - q_g) / r_g
    C   = (cp_g - cp_l) / r_g
    D   = (cp_l - cv_l) / r_g
    E   = (b_l - b_g) / r_g

    def residual(p, Tval):
        return (np.log(p + p_inf_g) - A - (B + E * p) / Tval
                - C * np.log(Tval) - D * np.log(p + p_inf_l))

    T = np.atleast_1d(T)
    out = np.empty_like(T, dtype=float)
    for i, Ti in enumerate(T):
        try:
            out[i] = brentq(residual, 1.0, 1.0e9, args=(Ti,),
                            xtol=1e-3, rtol=1e-10)
        except ValueError:
            out[i] = np.nan
    return out


def eos_psat_gibbs(T, liquid, vapor, p_bounds=(1.0, 1.0e9)):
    """Solve g_l(p,T)-g_v(p,T)=0 with a log-pressure bracket scan."""
    T = np.atleast_1d(T)
    out = np.empty_like(T, dtype=float)

    p_grid = np.logspace(np.log10(p_bounds[0]), np.log10(p_bounds[1]), 260)

    for i, Ti in enumerate(T):
        def residual(p):
            return eos_g(Ti, p, liquid) - eos_g(Ti, p, vapor)

        try:
            vals = np.array([residual(p) for p in p_grid])
            idx = None
            for j in range(len(p_grid) - 1):
                if not np.isfinite(vals[j]) or not np.isfinite(vals[j + 1]):
                    continue
                if vals[j] == 0.0 or vals[j] * vals[j + 1] < 0.0:
                    idx = j
                    break
            if idx is None:
                out[i] = np.nan
            else:
                out[i] = brentq(residual, p_grid[idx], p_grid[idx + 1],
                                xtol=1e-3, rtol=1e-10)
        except Exception:
            out[i] = np.nan
    return out


def eos_vl(T, p, liquid):
    return eos_v(T, p, liquid)


def eos_vg(T, p, vapor):
    return eos_v(T, p, vapor)


def eos_hl(T, p, liquid):
    return eos_h(T, p, liquid)


def eos_hg(T, p, vapor):
    return eos_h(T, p, vapor)


def eos_T_from_vp(v, p, phase):
    """
    Invert v(T,p) for T.  v(T,p) is affine in T at fixed p for every phase
    implemented here (SG/NASG/IG and ENASG liquid), so the inversion is
    closed-form.
    """
    if phase["eos"] == "ENASG":
        gamma, b1, cv = phase["gamma"], phase["b1"], phase["cv"]
        B    = phase["b0"] / (1.0 - b1)
        A    = gamma * phase["p_inf1"]
        Ctil = gamma * phase["p_inf0"] * (1.0 - b1) / (gamma - b1)
        return ((v - B) * (1.0 - b1) * (p + Ctil)
                / ((gamma - 1.0) * cv - (v - B) * (1.0 - b1) * A))
    return (v - phase["b"]) * (p + phase["p_inf"]) / (phase["cp"] - phase["cv"])


def eos_c(T, p, phase):
    """Sound speed c(p,T) for SG/NASG/IG and ENASG liquid."""
    if phase["eos"] == "ENASG":
        v = eos_v(T, p, phase)
        return np.sqrt(np.maximum(0.0, enasg_c2_from_p_v(p, v, phase)))
    gamma, b, p_inf = phase["gamma"], phase["b"], phase["p_inf"]
    rho = 1.0 / eos_v(T, p, phase)
    return np.sqrt(np.maximum(0.0, gamma * (p + p_inf) / (rho * (1.0 - b * rho))))


def eos_cl(T, p, liquid):
    return eos_c(T, p, liquid)


def eos_cg(T, p, vapor):
    return eos_c(T, p, vapor)


# =============================================================================
# Hugoniot constraint (Rice & Walsh 1957)
#
# Closes the Rankine-Hugoniot energy jump e2-e1 = 0.5*(p2+p1)*(v1-v2) with the
# phase's own caloric relation (via eos_h and eos_T_from_vp), so it applies to
# SG, NASG, and ENASG alike without a per-EOS closed-form Hugoniot.
# =============================================================================

def _hugoniot_residual_p2(p2, v2, v1, p1, e1, phase):
    T2 = eos_T_from_vp(v2, p2, phase)
    e2 = eos_h(T2, p2, phase) - p2 * v2
    return (e2 - e1) - 0.5 * (p2 + p1) * (v1 - v2)


def _hugoniot_p2_nasg_closed(v2, v1, p1, liquid):
    """
    Closed-form NASG principal Hugoniot pressure.

    Derived from the RH energy jump e2-e1 = 0.5*(p1+p2)*(v1-v2) using the
    NASG caloric EOS e(p,v) = (p + gamma*p_inf)*(v-b)/(gamma-1) + q.  The q
    terms cancel, giving the analytic formula.
    """
    g    = liquid["gamma"]
    pinf = liquid["p_inf"]
    b    = liquid["b"]
    if v2 <= b:
        return np.nan
    e0  = (p1 + g*pinf) * (v1 - b) / (g - 1.0)
    num = e0 + 0.5*p1*(v1 - v2) - g*pinf*(v2 - b) / (g - 1.0)
    den = (v2 - b) / (g - 1.0) - 0.5*(v1 - v2)
    if den <= 0.0:
        return np.nan
    p2 = num / den
    return p2 if (p2 > p1 and p2 < 1.0e12) else np.nan


def hugoniot_p2(v2, v1, p1, e1, phase, p_bounds=(1.0, 1.0e12)):
    """Shocked pressure p2 at specific volume v2 from the Rankine-Hugoniot jump."""
    if phase["eos"] == "NASG":
        return _hugoniot_p2_nasg_closed(v2, v1, p1, phase)
    grid = np.logspace(np.log10(p_bounds[0]), np.log10(p_bounds[1]), 200)
    vals = np.array([_hugoniot_residual_p2(p, v2, v1, p1, e1, phase) for p in grid])
    for j in range(len(grid) - 1):
        if not np.isfinite(vals[j]) or not np.isfinite(vals[j + 1]):
            continue
        if vals[j] == 0.0 or vals[j] * vals[j + 1] < 0.0:
            return brentq(_hugoniot_residual_p2, grid[j], grid[j + 1],
                          args=(v2, v1, p1, e1, phase), xtol=1.0, rtol=1e-10)
    return np.nan


def hugoniot_sse(phase, rw=RICE_WALSH):
    """
    Sum of squared relative pressure residuals between the EOS principal
    Hugoniot and the Rice & Walsh (1957) (up, us) shock data, and the number
    of points for which the EOS admits a shocked state (e.g. rho2 < 1/b_l
    for NASG).  Points outside the EOS's admissible range are skipped.

    Returns (sse, n_used).
    """
    rho0, p0 = rw["rho0"], rw["p0"]
    v1 = 1.0 / rho0
    p1 = p0
    T1 = eos_T_from_vp(v1, p1, phase)
    e1 = eos_h(T1, p1, phase) - p1 * v1

    rho2 = rho0 * rw["us"] / (rw["us"] - rw["up"])
    p2_rw = p0 + rho0 * rw["us"] * rw["up"]

    sse, n = 0.0, 0
    for v2, p2_ref in zip(1.0 / rho2, p2_rw):
        p2 = hugoniot_p2(v2, v1, p1, e1, phase)
        if not np.isfinite(p2):
            continue
        sse += ((p2 - p2_ref) / p2_ref) ** 2
        n += 1
    return sse, n


def hugoniot_rms(phase, rw=RICE_WALSH):
    """
    RMS relative pressure error [%] of the EOS Hugoniot vs Rice & Walsh (1957),
    restricted to points within the EOS's admissible density range.

    Returns (rms_percent, n_used, n_total).
    """
    sse, n = hugoniot_sse(phase, rw)
    n_total = len(rw["up"])
    if n == 0:
        return np.nan, 0, n_total
    return 100.0 * np.sqrt(sse / n), n, n_total


def hugoniot_objective(phase, rw=RICE_WALSH, penalty=1.0):
    """
    Penalized sum of squared relative pressure residuals, for use as a fit
    objective.  Points outside the EOS's admissible range (hugoniot_p2 -> nan,
    e.g. rho2 >= 1/b_l for NASG) are assigned a fixed relative residual of
    `penalty` instead of being dropped, so the objective cannot be minimized
    merely by shrinking the number of admissible points (hugoniot_sse's n).

    Returns (sse, n_used).
    """
    rho0, p0 = rw["rho0"], rw["p0"]
    v1 = 1.0 / rho0
    p1 = p0
    T1 = eos_T_from_vp(v1, p1, phase)
    e1 = eos_h(T1, p1, phase) - p1 * v1

    rho2  = rho0 * rw["us"] / (rw["us"] - rw["up"])
    p2_rw = p0 + rho0 * rw["us"] * rw["up"]

    sse, n = 0.0, 0
    for v2, p2_ref in zip(1.0 / rho2, p2_rw):
        p2 = hugoniot_p2(v2, v1, p1, e1, phase)
        if np.isfinite(p2):
            sse += ((p2 - p2_ref) / p2_ref) ** 2
            n += 1
        else:
            sse += penalty ** 2
    return sse, n


def _minimize_hugoniot(objective, bounds, n_grid=100, xatol=1.0):
    """
    Minimize a Hugoniot objective over `bounds`.

    `hugoniot_objective` has jump discontinuities where Rice & Walsh points
    enter/exit the EOS's admissible range (n_used changes), which can strand
    a Brent search at an arbitrary, possibly unphysical, bound. A coarse
    log-spaced grid scan locates the best branch first; Brent then refines
    within the (locally smooth) bracket around that grid point.

    `objective` must return a (sse, n_used) pair, as hugoniot_objective does.

    Returns (value, feasible). feasible=False means no point in `bounds`
    admits any Rice & Walsh shock state (n_used==0 at every grid point), so
    `value` is just the lower bound and should not be trusted as a fit.
    """
    grid = np.logspace(np.log10(bounds[0]), np.log10(bounds[1]), n_grid)
    results = [objective(p) for p in grid]
    vals   = np.array([r[0] for r in results])
    n_arr  = np.array([r[1] for r in results], dtype=int)
    i_best = int(np.argmin(vals))
    if n_arr[i_best] == 0:
        return grid[i_best], False
    lo = grid[max(i_best - 1, 0)]
    hi = grid[min(i_best + 1, len(grid) - 1)]
    if lo == hi:
        return grid[i_best], True
    result = minimize_scalar(lambda p: objective(p)[0], bounds=(lo, hi),
                              method="bounded", options={"xatol": xatol})
    return result.x, True


# =============================================================================
# Output: text table and figure
# =============================================================================

def _fmt_value(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return f"{'--':>18}"
    return f"{x:>18.4e}"


def print_coefficients(liquid, vapor, fluid_name, T_range, eos_name="NASG"):
    Tmin, Tmax = T_range
    header = f"{eos_name} coefficients for {fluid_name} in [{Tmin:.0f}-{Tmax:.0f}] K"
    print(header)
    print("=" * len(header))

    if liquid["eos"] == "ENASG":
        rows = [
            ("c_p    [J/(kg.K)]", None,                 vapor["cp"]),
            ("c_v    [J/(kg.K)]", liquid["cv"],         vapor["cv"]),
            ("gamma",             liquid["gamma"],      vapor["gamma"]),
            ("p_inf0 [Pa]",       liquid["p_inf0"],     vapor["p_inf0"]),
            ("p_inf1 [Pa/K]",     liquid["p_inf1"],     vapor["p_inf1"]),
            ("b0     [m^3/kg]",   liquid["b0"],         vapor["b0"]),
            ("b1",                liquid["b1"],         vapor["b1"]),
            ("q      [J/kg]",     liquid["q"],          vapor["q"]),
            ("q'     [J/(kg.K)]", liquid["q_prime"],    vapor["q_prime"]),
        ]
    else:
        rows = [
            ("c_p    [J/(kg.K)]", liquid["cp"],      vapor["cp"]),
            ("c_v    [J/(kg.K)]", liquid["cv"],      vapor["cv"]),
            ("gamma",             liquid["gamma"],   vapor["gamma"]),
            ("p_inf  [Pa]",       liquid["p_inf"],   vapor["p_inf"]),
            ("b      [m^3/kg]",   liquid["b"],       vapor["b"]),
            ("q      [J/kg]",     liquid["q"],       vapor["q"]),
            ("q'     [J/(kg.K)]", liquid["q_prime"], vapor["q_prime"]),
        ]

    print(f"{'Coefficient':<20}{'Liquid':>18}{'Vapor':>18}")
    print("-" * 56)
    for label, lv, gv in rows:
        print(f"{label:<20}{_fmt_value(lv)}{_fmt_value(gv)}")
    print()


def plot_saturation(results, sat_data,
                    T_plot_range=None, T_fit_range=None, n_curve=300,
                    savepath=None):
    """
    4x2 figure, arranged like Chiapolino & Saurel (2018), Fig. 1, with an
    extra row for the saturated sound speeds:
        p_sat, L_v
        h_l,   h_v
        rho_l, rho_v
        c_l,   c_v

    IAPWS-IF97 is plotted as a solid black line; EOS curves are plotted as
    colored markers.  The fitting window is shaded when T_fit_range is provided.

    Parameters
    ----------
    results : list of (label, liquid_dict, vapor_dict) tuples
    """
    T_ref    = sat_data[:, 0]
    p_ref    = sat_data[:, 1]
    hl_ref   = sat_data[:, 2]
    hg_ref   = sat_data[:, 3]
    vl_ref   = sat_data[:, 4]
    vg_ref   = sat_data[:, 5]
    wl_ref   = sat_data[:, 10]
    wg_ref   = sat_data[:, 11]
    Lv_ref   = hg_ref - hl_ref
    rhol_ref = 1.0 / vl_ref
    rhov_ref = 1.0 / vg_ref

    if T_plot_range is None:
        T_plot_range = (T_ref.min(), T_ref.max())

    T_curve   = np.r_[np.linspace(T_plot_range[0], T_plot_range[1], n_curve, endpoint=False),
                    T_plot_range[1]]
    colors     = ["b", "r", "c", "m"]
    linestyles = ["--", "-.", ":", (0, (7, 2, 1.5, 2))]

    fig, axes = plt.subplots(4, 2, figsize=(14.0, 20.5), sharex=True)

    panel_defs = [
        (axes[0, 0], r"$p\;[\mathrm{bar}]$"),
        (axes[0, 1], r"$L_v\;[\mathrm{kJ\,kg^{-1}}]$"),
        (axes[1, 0], r"$h_l\;[\mathrm{kJ\,kg^{-1}}]$"),
        (axes[1, 1], r"$h_v\;[\mathrm{kJ\,kg^{-1}}]$"),
        (axes[2, 0], r"$\rho_l\;[\mathrm{kg\,m^{-3}}]$"),
        (axes[2, 1], r"$\rho_v\;[\mathrm{kg\,m^{-3}}]$"),
        (axes[3, 0], r"$c_l\;[\mathrm{m\,s^{-1}}]$"),
        (axes[3, 1], r"$c_v\;[\mathrm{m\,s^{-1}}]$"),
    ]

    ref_curves = [
        p_ref * 1.0e-5,
        Lv_ref * 1.0e-3,
        hl_ref * 1.0e-3,
        hg_ref * 1.0e-3,
        rhol_ref,
        rhov_ref,
        wl_ref,
        wg_ref,
    ]

    x_ticks = [300, 400, 500, 600]
    x_tick_labels = ['300', '400', '500', '600']

    for (ax, ylabel), y_ref in zip(panel_defs, ref_curves):
        if T_fit_range is not None:
            ax.axvspan(T_fit_range[0], T_fit_range[1],
                       color="0.85", alpha=0.45, zorder=0)
        ax.plot(T_ref, y_ref, "k-", lw=4.4, label="IAPWS-IF97", zorder=3)
        ax.set_xlim(T_plot_range)
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels)
        ax.tick_params(axis="both", labelsize=28, width=4.0, length=8)
        ax.set_ylabel(ylabel)
        yfmt = ScalarFormatter(useMathText=True)
        yfmt.set_powerlimits((-2, 2))
        ax.yaxis.set_major_formatter(yfmt)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2), useMathText=True)
        for spine in ax.spines.values():
            spine.set_linewidth(2.2)
        ax.grid(True, which="both", alpha=0.3, lw=0.9)

    for ax in axes[:, 1]:
        ax.yaxis.set_label_position("right")
        ax.yaxis.tick_right()

    for ax in axes[-1, :]:
        ax.set_xlabel(r"$T\;[\mathrm{K}]$", fontsize=34)

    for idx, (label, liquid, vapor) in enumerate(results):
        p_c    = eos_psat(T_curve, liquid, vapor)
        hl_c   = eos_hl(T_curve, p_c, liquid)
        hg_c   = eos_hg(T_curve, p_c, vapor)
        vl_c   = eos_vl(T_curve, p_c, liquid)
        vg_c   = eos_vg(T_curve, p_c, vapor)
        Lv_c   = hg_c - hl_c
        rhol_c = 1.0 / vl_c
        rhov_c = 1.0 / vg_c
        cl_c   = eos_cl(T_curve, p_c, liquid)
        cv_c   = eos_cg(T_curve, p_c, vapor)

        eos_curves = [
            p_c * 1.0e-5,
            Lv_c * 1.0e-3,
            hl_c * 1.0e-3,
            hg_c * 1.0e-3,
            rhol_c,
            rhov_c,
            cl_c,
            cv_c,
        ]
        color = colors[idx % len(colors)]
        ls    = linestyles[idx % len(linestyles)]
        for (ax, _), y_c in zip(panel_defs, eos_curves):
            ax.plot(T_curve, y_c, linestyle=ls, lw=4.0, color=color,
                    label=label, zorder=4)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    shade_proxy = plt.Rectangle((0, 0), 1, 1, facecolor="0.85", alpha=0.45,
                                edgecolor="none", label="fit window")
    if T_fit_range is not None:
        handles = handles + [shade_proxy]
        labels  = labels  + ["fit window"]
    axes[0, 0].legend(handles, labels, markerscale=1.5, frameon=False, handlelength=2.8, fontsize=22)

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight")
        print(f"Saved {savepath}")
    return fig


def plot_hugoniot(results, rw=RICE_WALSH, n_curve=200, savepath=None):
    """
    1x3 figure: principal Hugoniot M_s, p_2, rho_2/rho_1 vs piston velocity u_p,
    compared against the Rice & Walsh (1957) water shock data.

    Each EOS's locus is traced by sweeping the shocked specific volume v2 from
    just below the ambient v1 down to the EOS's admissible limit, reusing the
    same eos_T_from_vp/eos_h/hugoniot_p2 machinery as the Hugoniot constraint.
    u_s is normalized by each EOS's own ambient sound speed a0=eos_c(T1,p0) to
    give the shock Mach number M_s=u_s/a0; the Rice & Walsh data are
    normalized by the IAPWS-IF97 sound speed of water at (rho0,p0).

    Parameters
    ----------
    results : list of (label, liquid_dict, vapor_dict) tuples
    """
    rho0, p0 = rw["rho0"], rw["p0"]
    v1 = 1.0 / rho0
    p1 = p0

    a0_rw = IAPWS97(T=293.15, P=p0 * 1.0e-6).w

    colors     = ["b", "r", "c", "m"]
    linestyles = ["--", "-.", ":", (0, (7, 2, 1.5, 2))]

    fig, axes = plt.subplots(1, 3, figsize=(19.0, 5.6))
    panel_defs = [
        (axes[0], r"$M_s$",              (1.0, 5.0)),
        (axes[1], r"$p_2\;[\mathrm{GPa}]$", (0.0, 12.0)),
        (axes[2], r"$\rho_2/\rho_1$",    (1.0, 2.0)),
    ]

    for idx, (label, liquid, _) in enumerate(results):
        T1 = eos_T_from_vp(v1, p1, liquid)
        e1 = eos_h(T1, p1, liquid) - p1 * v1
        a0 = eos_c(T1, p1, liquid)

        v2_grid = np.linspace(v1 * 0.999, v1 * 0.05, n_curve)
        p2_grid = np.array([hugoniot_p2(v2, v1, p1, e1, liquid) for v2 in v2_grid])
        ok = np.isfinite(p2_grid) & (p2_grid > p1)
        v2_ok, p2_ok = v2_grid[ok], p2_grid[ok]

        up = np.sqrt((p2_ok - p1) * (v1 - v2_ok))
        us = (p2_ok - p1) / (rho0 * np.maximum(up, 1.0e-12))
        rho_ratio = v1 / v2_ok

        color = colors[idx % len(colors)]
        ls    = linestyles[idx % len(linestyles)]
        for (ax, _, _), y in zip(panel_defs, [us / a0, p2_ok * 1.0e-9, rho_ratio]):
            ax.plot(up, y, linestyle=ls, lw=4.0, color=color, label=label, zorder=4)

    # Rice & Walsh (1957) reference data
    up_rw, us_rw = rw["up"], rw["us"]
    rho2_rw = rho0 * us_rw / (us_rw - up_rw)
    p2_rw   = p0 + rho0 * us_rw * up_rw
    for (ax, _, _), y in zip(panel_defs, [us_rw / a0_rw, p2_rw * 1.0e-9, rho2_rw / rho0]):
        ax.plot(up_rw, y, "ks", ms=12, mfc="none", mew=2.5,
                label=r"Rice \& Walsh (1957)", zorder=5)

    for ax, ylabel, ylim in panel_defs:
        ax.set_xlabel(r"$u_p\;[\mathrm{m\,s^{-1}}]$")
        ax.set_xlim(0.0, 3000.0)
        ax.set_ylim(*ylim)
        ax.tick_params(axis="both", labelsize=28, width=4.0, length=8)
        ax.set_ylabel(ylabel)
        for spine in ax.spines.values():
            spine.set_linewidth(2.2)
        ax.grid(True, which="both", alpha=0.3, lw=0.9)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, markerscale=1.5,
               frameon=False, handlelength=2.8, fontsize=22,
               bbox_to_anchor=(0.5, 1.12))

    fig.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight")
        print(f"Saved {savepath}")
    return fig


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit SG, NASG, and/or ENASG EOS coefficients to IAPWS data."
    )
    parser.add_argument(
        "--mode", choices=["nasg", "sg", "enasg", "all"], default="all",
        help="Which EOS(es) to fit and plot. Default: all.",
    )
    parser.add_argument(
        "--constraint", choices=["soundspeed", "hugoniot"], default="hugoniot",
        help=("Closure for the liquid EOS's free parameter (p_inf for SG/NASG, "
              "C for ENASG). 'soundspeed' uses the ambient reference sound "
              "speed (Le Metayer & Saurel 2016 Eq. 68 / Chiapolino & Saurel "
              "2018 Eq. 42). 'hugoniot' instead least-squares-fits the EOS "
              "principal Hugoniot to the Rice & Walsh (1957) shock data. "
              "Default: hugoniot."),
    )
    args = parser.parse_args()

    T_fit  = FIT_CONFIG["T_fit"]
    T_plot = FIT_CONFIG["T_plot"]

    sat_ref  = iapws_saturation_table(T_fit[0],  T_fit[1],  n_points=FIT_CONFIG["n_fit"])
    sat_wide = iapws_saturation_table(T_plot[0], T_plot[1], n_points=FIT_CONFIG["n_plot"])

    # Plot/legend order is intentionally kept as:
    #   IAPWS-IF97, SG, NASG, ENASG
    # where IAPWS-IF97 is added inside plot_saturation().
    plot_results = []

    if args.mode in ("sg", "all"):
        liq_sg, vap_sg = compute_sg_coefficients(sat_ref, WATER_REFERENCE_STATE, T_fit,
                                                  constraint=args.constraint)
        print_coefficients(liq_sg, vap_sg, "water", T_fit, eos_name="SG")
        print("Reference (Flatten & Lund 2011):")
        print("  Liquid: cp=4267  cv=1816  gamma=2.35  p_inf=1.0e9 Pa  b=0")
        print("          q=-1.167e6 J/kg   q'=0")
        print("  Vapor : cp=1487  cv=1040  gamma=1.43  p_inf=0  b=0")
        print("          q=2.030e6 J/kg   q'=-23400 J/(kg.K)")
        rms_sg, n_used_sg, n_tot_sg = hugoniot_rms(liq_sg)
        print(f"  Hugoniot RMS (Rice & Walsh 1957): {rms_sg:.2f}% ({n_used_sg}/{n_tot_sg} points in range)\n")
        plot_results.append(("SG", liq_sg, vap_sg))

    if args.mode in ("nasg", "all"):
        liq_nasg, vap_nasg = compute_nasg_coefficients(sat_ref, WATER_REFERENCE_STATE, T_fit,
                                                        constraint=args.constraint)
        print_coefficients(liq_nasg, vap_nasg, "water", T_fit, eos_name="NASG")
        print("Reference (Le Metayer & Saurel 2016, Table V):")
        print("  Liquid: cp=4285  cv=3610  gamma=1.19  p_inf=7.028e8 Pa  b=6.61e-4 m^3/kg")
        print("          q=-1.177788e6 J/kg   q'=0")
        print("  Vapor : cp=1401  cv=955   gamma=1.47  p_inf=0  b=0")
        print("          q=2.077616e6 J/kg   q'=14317 J/(kg.K)")
        rms_nasg, n_used_nasg, n_tot_nasg = hugoniot_rms(liq_nasg)
        print(f"  Hugoniot RMS (Rice & Walsh 1957): {rms_nasg:.2f}% ({n_used_nasg}/{n_tot_nasg} points in range)\n")
        plot_results.append(("NASG", liq_nasg, vap_nasg))

    if args.mode in ("enasg", "all"):
        liq_enasg, vap_enasg = compute_enasg_coefficients(
            sat_ref,
            WATER_ENASG_REFERENCE_STATE,
            WATER_VAPOR_IG_REFERENCE_STATE,
            T_fit,
            constraint=args.constraint,
        )
        print_coefficients(liq_enasg, vap_enasg, "water", T_fit, eos_name="ENASG")
        print("Reference (Chiapolino & Saurel 2018, Table 1):")
        print("  Liquid: gamma=1.0147  cv=4014  b1=-0.6050  b0=1.5196e-3 m^3/kg")
        print("          p_inf1=-471025 Pa/K  p_inf0=307078403 Pa")
        print("          q=-1.112426e6 J/kg   q'=-22049 J/(kg.K)")
        print("  Vapor : gamma=1.3079  cv=1500  p_inf0=p_inf1=b0=b1=0")
        print("          q=1.947630e6 J/kg    q'=1136 J/(kg.K)")
        rms_enasg, n_used_enasg, n_tot_enasg = hugoniot_rms(liq_enasg)
        print(f"  Hugoniot RMS (Rice & Walsh 1957): {rms_enasg:.2f}% ({n_used_enasg}/{n_tot_enasg} points in range)\n")
        plot_results.append(("ENASG", liq_enasg, vap_enasg))

    savepath   = f"fit_eos_{args.mode}_{args.constraint}_saturation.pdf"
    savepath_h = f"fit_eos_{args.mode}_{args.constraint}_hugoniot.pdf"

    plot_saturation(
        plot_results, sat_wide,
        T_plot_range=T_plot,
        T_fit_range=T_fit,
        n_curve=FIT_CONFIG["n_curve"],
        savepath=savepath,
    )

    plot_hugoniot(
        plot_results,
        n_curve=FIT_CONFIG["n_hugoniot"],
        savepath=savepath_h,
    )
