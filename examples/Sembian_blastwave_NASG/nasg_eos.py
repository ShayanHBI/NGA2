"""
Noble-Abel Stiffened-Gas (NASG) Equation of State: coefficient fitter + plotter.

Implements the methodology described in:
    Le Metayer, O. and Saurel, R., "The Noble-Abel Stiffened-Gas equation of state",
    Physics of Fluids 28, 046102 (2016).

The NASG caloric EOS is:
    P = (gamma - 1) * (e - q) / (v - b)  -  gamma * P_inf

This script:
  1. Pulls saturation data from IAPWS-IF97 (via the `iapws` package).
  2. Fits the NASG parameters (gamma, Cv, P_inf, b, q, q_prime) for liquid
     and vapor phases by least-squares on the saturation curves (Sec. IV).
  3. Reproduces Fig. 6 of the paper: 2x3 grid of (P_sat, L_v, h_l, h_g, v_l, v_g)
     versus T, with IAPWS data points overlaid on NASG analytical curves.

Default: liquid water / steam in the range 300-500 K (paper Table V).

Dependencies: numpy, scipy, matplotlib, iapws  (pip install iapws)

Usage:
    python nasg_eos.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from iapws import IAPWS97


# =============================================================================
# Reference data
# =============================================================================

# Liquid water reference state (paper Section IV.C)
WATER_REFERENCE_STATE = {
    # "rho0": 957.74,    # kg/m^3
    # "P0":   1.0453e5,  # Pa
    # "c0":   1542.0,    # m/s

    # "rho0": 997.0,     # kg/m^3
    # "P0":   3169.8,    # Pa
    # "c0":   1497.0,    # m/s

    "rho0": 999.9747,     # kg/m^3
    "P0":   1.0100e5,    # Pa
    "c0":   1421.53,    # m/s
}


def iapws_saturation_table(T_min, T_max, n_points=101):
    """
    Build a saturation-curve table from IAPWS-IF97.

    Returns
    -------
    ndarray, shape (N, 6) with columns (T, P_sat, h_l, h_g, v_l, v_g) in SI units.
    Rows where IAPWS97 raises (e.g. above the critical point) are skipped.
    """
    T = np.linspace(T_min, T_max, n_points)
    rows = []
    for Ti in T:
        try:
            liq = IAPWS97(T=Ti, x=0)
            vap = IAPWS97(T=Ti, x=1)
            rows.append([Ti, liq.P * 1e6, liq.h * 1e3, vap.h * 1e3,
                         liq.v, vap.v])
        except Exception:
            pass
    return np.array(rows)


# =============================================================================
# Phase fits (Section IV of the paper)
# =============================================================================

def fit_vapor_phase(T, P, hg, vg):
    """
    Vapor-phase coefficients with b_g = 0 and P_inf_g = 0 (ideal-gas vapor).
    Implements Eqs. (50), (51), (55).
    """
    T_mean  = T.mean()
    hg_mean = hg.mean()

    # h_g(T) = CP_g * T + q_g  (least squares)
    CP_g = np.sum(T * (hg - hg_mean)) / np.sum(T * (T - T_mean))
    q_g  = hg_mean - CP_g * T_mean

    # v_g(T) = (CP_g - Cv_g) * T / P_sat(T)  (least squares)
    R_g  = np.sum(vg * T / P) / np.sum((T / P) ** 2)
    Cv_g = CP_g - R_g

    return {
        "CP":    CP_g,
        "Cv":    Cv_g,
        "gamma": CP_g / Cv_g,
        "q":     q_g,
        "b":     0.0,
        "P_inf": 0.0,
    }


def _fit_liquid_for_given_Pinf(T, P, hl, vl, P_inf_l):
    """
    Given a candidate P_inf,l, solve Eqs. (60)-(61), (64)-(65) for
    (CP_l, Cv_l, q_l, b_l).
    """
    T_mean  = T.mean()
    P_mean  = P.mean()
    hl_mean = hl.mean()
    vl_mean = vl.mean()

    TP      = T / (P + P_inf_l)
    TP_mean = TP.mean()

    # (CP_l - Cv_l) and b_l from v_l data
    R_l = np.sum(TP * (vl - vl_mean)) / np.sum(TP * (TP - TP_mean))
    b_l = vl_mean - R_l * TP_mean

    # CP_l and q_l from h_l data (uses b_l already determined)
    num_CP = np.sum(T * (hl - hl_mean)) - b_l * np.sum(T * (P - P_mean))
    den_CP = np.sum(T * (T - T_mean))
    CP_l = num_CP / den_CP
    q_l  = hl_mean - CP_l * T_mean - b_l * P_mean

    return CP_l, CP_l - R_l, q_l, b_l


def fit_liquid_phase(T, P, hl, vl, ref_state, P_inf_bounds=(1.0e7, 5.0e9)):
    """
    Close the liquid-phase system with the sound-speed reference state (Eq. 68):

        P0 + P_inf_l - (Cv_l / CP_l) * rho0 * c0^2 * (1 - b_l * rho0) = 0

    All other liquid coefficients are determined as functions of P_inf_l.
    """
    rho0 = ref_state["rho0"]
    P0   = ref_state["P0"]
    c0   = ref_state["c0"]

    def residual(P_inf_l):
        CP_l, Cv_l, _, b_l = _fit_liquid_for_given_Pinf(T, P, hl, vl, P_inf_l)
        return P0 + P_inf_l - (Cv_l / CP_l) * rho0 * c0**2 * (1.0 - b_l * rho0)

    a, b = P_inf_bounds
    fa, fb = residual(a), residual(b)
    if fa * fb > 0:
        for _ in range(8):
            b *= 5.0
            fb = residual(b)
            if fa * fb < 0:
                break
        else:
            raise RuntimeError(
                f"Could not bracket root for P_inf_l. "
                f"f({a:.2e})={fa:.3e}, f({b:.2e})={fb:.3e}"
            )

    P_inf_l = brentq(residual, a, b, xtol=1.0, rtol=1e-10)
    CP_l, Cv_l, q_l, b_l = _fit_liquid_for_given_Pinf(T, P, hl, vl, P_inf_l)

    return {
        "CP":    CP_l,
        "Cv":    Cv_l,
        "gamma": CP_l / Cv_l,
        "q":     q_l,
        "b":     b_l,
        "P_inf": P_inf_l,
    }


def fit_entropy_constants(T, P, liquid, vapor):
    """
    Determine q'_l and q'_g via the saturated-vapor-pressure relation.
    Convention: q'_l = 0 (Eqs. 41, 69-72).
    """
    CP_l, Cv_l, P_inf_l = liquid["CP"], liquid["Cv"], liquid["P_inf"]
    CP_g, Cv_g          = vapor["CP"],  vapor["Cv"]
    b_l, b_g            = liquid["b"],  vapor["b"]
    q_l, q_g            = liquid["q"],  vapor["q"]

    Rg = CP_g - Cv_g
    B = (q_l - q_g) / Rg
    C = (CP_g - CP_l) / Rg
    D = (CP_l - Cv_l) / Rg
    E = (b_l - b_g)   / Rg

    A = np.mean(
        np.log(P)
        - (B + E * P) / T
        - C * np.log(T)
        - D * np.log(P + P_inf_l)
    )

    q_prime_l = 0.0
    q_prime_g = A * Rg - (CP_l - CP_g) + q_prime_l
    return q_prime_l, q_prime_g


def compute_nasg_coefficients(sat_data, ref_state, T_range):
    """
    Full NASG fit for a liquid/vapor pair.

    Parameters
    ----------
    sat_data : (N, 6) ndarray with columns (T, P_sat, h_l, h_g, v_l, v_g), SI units.
    ref_state : dict with keys 'rho0', 'P0', 'c0' (liquid sound-speed anchor).
    T_range : (T_min, T_max) tuple in K -- data is filtered to this window.

    Returns
    -------
    liquid, vapor : dicts with CP, Cv, gamma, P_inf, b, q, q_prime.
    """
    Tmin, Tmax = T_range
    mask = (sat_data[:, 0] >= Tmin) & (sat_data[:, 0] <= Tmax)
    d = sat_data[mask]
    if len(d) < 4:
        raise ValueError(
            f"Need at least 4 data points in [{Tmin}, {Tmax}] K, got {len(d)}."
        )

    T, P, hl, hg, vl, vg = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4], d[:, 5]

    vapor  = fit_vapor_phase(T, P, hg, vg)
    liquid = fit_liquid_phase(T, P, hl, vl, ref_state)
    qp_l, qp_g = fit_entropy_constants(T, P, liquid, vapor)
    liquid["q_prime"] = qp_l
    vapor["q_prime"]  = qp_g
    return liquid, vapor


# =============================================================================
# NASG analytical saturation curves (Eqs. 41, 43-46)
# =============================================================================

def nasg_psat(T, liquid, vapor):
    """
    Solve Eq. (41) for the NASG saturated vapor pressure at each T.

        ln(P + P_inf,g) = A + (B + E*P)/T + C*ln(T) + D*ln(P + P_inf,l)
    """
    CP_l, Cv_l, P_inf_l = liquid["CP"], liquid["Cv"], liquid["P_inf"]
    CP_g, Cv_g, P_inf_g = vapor["CP"],  vapor["Cv"],  vapor["P_inf"]
    b_l, b_g            = liquid["b"],  vapor["b"]
    q_l, q_g            = liquid["q"],  vapor["q"]
    qp_l, qp_g          = liquid["q_prime"], vapor["q_prime"]

    Rg = CP_g - Cv_g
    A = (CP_l - CP_g + qp_g - qp_l) / Rg
    B = (q_l - q_g) / Rg
    C = (CP_g - CP_l) / Rg
    D = (CP_l - Cv_l) / Rg
    E = (b_l - b_g)   / Rg

    def residual(P, Tval):
        return (np.log(P + P_inf_g) - A - (B + E * P) / Tval
                - C * np.log(Tval) - D * np.log(P + P_inf_l))

    T = np.atleast_1d(T)
    out = np.empty_like(T, dtype=float)
    P_lo, P_hi = 1.0, 1.0e9
    for i, Ti in enumerate(T):
        try:
            out[i] = brentq(residual, P_lo, P_hi, args=(Ti,),
                            xtol=1e-3, rtol=1e-10)
        except ValueError:
            out[i] = np.nan
    return out


def nasg_vl(T, P, liquid):
    """Eq. (46): saturated liquid specific volume."""
    return (liquid["CP"] - liquid["Cv"]) * T / (P + liquid["P_inf"]) + liquid["b"]


def nasg_vg(T, P, vapor):
    """Eq. (45): saturated vapor specific volume."""
    return (vapor["CP"] - vapor["Cv"]) * T / (P + vapor["P_inf"]) + vapor["b"]


def nasg_hl(T, P, liquid):
    """Eq. (44): saturated liquid specific enthalpy."""
    return liquid["CP"] * T + liquid["b"] * P + liquid["q"]


def nasg_hg(T, P, vapor):
    """Eq. (43): saturated vapor specific enthalpy."""
    return vapor["CP"] * T + vapor["b"] * P + vapor["q"]


# =============================================================================
# Output: text table and figure
# =============================================================================

def print_coefficients(liquid, vapor, fluid_name, T_range):
    Tmin, Tmax = T_range
    header = f"NASG coefficients for {fluid_name} in [{Tmin:.0f}-{Tmax:.0f}] K"
    print(header)
    print("=" * len(header))
    rows = [
        ("C_P    [J/(kg.K)]", liquid["CP"],      vapor["CP"]),
        ("C_v    [J/(kg.K)]", liquid["Cv"],      vapor["Cv"]),
        ("gamma  [-]",        liquid["gamma"],   vapor["gamma"]),
        ("P_inf  [Pa]",       liquid["P_inf"],   vapor["P_inf"]),
        ("b      [m^3/kg]",   liquid["b"],       vapor["b"]),
        ("q      [J/kg]",     liquid["q"],       vapor["q"]),
        ("q'     [J/(kg.K)]", liquid["q_prime"], vapor["q_prime"]),
    ]
    print(f"{'Coefficient':<20}{'Liquid':>18}{'Vapor':>18}")
    print("-" * 56)
    for label, lv, gv in rows:
        print(f"{label:<20}{lv:>18.4e}{gv:>18.4e}")
    print()


def plot_saturation(liquid, vapor, sat_data, T_fit_range,
                    T_plot_range=None, n_curve=400,
                    fluid_name="water/steam", savepath=None):
    """
    2x3 figure mirroring Fig. 6 of the paper:
    P_sat, L_v, h_l, h_g, v_l, v_g vs T. IAPWS data as dots, NASG as lines.
    The fit window is shaded.
    """
    T_ref  = sat_data[:, 0]
    P_ref  = sat_data[:, 1]
    hl_ref = sat_data[:, 2]
    hg_ref = sat_data[:, 3]
    vl_ref = sat_data[:, 4]
    vg_ref = sat_data[:, 5]
    Lv_ref = hg_ref - hl_ref

    if T_plot_range is None:
        T_plot_range = (T_ref.min(), T_ref.max())

    T_curve = np.linspace(*T_plot_range, n_curve)
    P_curve  = nasg_psat(T_curve, liquid, vapor)
    hl_curve = nasg_hl(T_curve, P_curve, liquid)
    hg_curve = nasg_hg(T_curve, P_curve, vapor)
    vl_curve = nasg_vl(T_curve, P_curve, liquid)
    vg_curve = nasg_vg(T_curve, P_curve, vapor)
    Lv_curve = hg_curve - hl_curve

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    fig.suptitle(
        f"NASG vs IAPWS-IF97 saturation curves for {fluid_name}  "
        f"(coefficients fitted on T \u2208 "
        f"[{T_fit_range[0]:.0f}, {T_fit_range[1]:.0f}] K)",
        fontsize=12,
    )

    panels = [
        (axes[0, 0], r"$P_{sat}$ [bar]",   P_ref  * 1e-5,  P_curve  * 1e-5,  False),
        (axes[0, 1], r"$L_v$ [kJ/kg]",     Lv_ref * 1e-3,  Lv_curve * 1e-3,  False),
        (axes[0, 2], r"$h_l$ [kJ/kg]",     hl_ref * 1e-3,  hl_curve * 1e-3,  False),
        (axes[1, 0], r"$h_g$ [kJ/kg]",     hg_ref * 1e-3,  hg_curve * 1e-3,  False),
        (axes[1, 1], r"$v_l$ [m$^3$/kg]",  vl_ref,         vl_curve,         False),
        (axes[1, 2], r"$v_g$ [m$^3$/kg]",  vg_ref,         vg_curve,         True ),
    ]
    for ax, label, y_ref, y_curve, log_y in panels:
        ax.plot(T_ref, y_ref, "o", color="tab:red", markersize=3.5,
                label="IAPWS-IF97", zorder=3)
        ax.plot(T_curve, y_curve, "-", color="tab:blue", linewidth=2,
                label="NASG", zorder=2)
        ax.axvspan(T_fit_range[0], T_fit_range[1], color="grey", alpha=0.10,
                   zorder=1, label="fit window")
        if log_y:
            ax.set_yscale("log")
        ax.set_xlabel("T [K]")
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if savepath:
        fig.savefig(savepath, dpi=140, bbox_inches="tight")
        print(f"Saved figure to {savepath}")
    return fig


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    # --- Configuration -------------------------------------------------------
    FLUID_NAME = "water/steam"
    T_FIT      = (270.0, 500.0)   # window used to fit NASG coefficients
    T_PLOT     = (275.0, 640.0)   # wider window to display reference data
    REF_STATE  = WATER_REFERENCE_STATE
    # ------------------------------------------------------------------------

    # IAPWS data: a dense table inside the fit window for least-squares,
    # plus a wider table for the plot.
    sat_fit  = iapws_saturation_table(T_FIT[0],  T_FIT[1],  n_points=101)
    sat_wide = iapws_saturation_table(T_PLOT[0], T_PLOT[1], n_points=51)

    liquid, vapor = compute_nasg_coefficients(sat_fit, REF_STATE, T_FIT)
    print_coefficients(liquid, vapor, FLUID_NAME, T_FIT)

    print("Values from Table V of Le Metayer & Saurel (2016):")
    print("  Liquid: CP=4285  Cv=3610  gamma=1.19  P_inf=7.028e8 Pa  b=6.61e-4 m^3/kg")
    print("          q=-1.177788e6 J/kg   q'=0")
    print("  Vapor : CP=1401  Cv=955   gamma=1.47  P_inf=0          b=0")
    print("          q=2.077616e6 J/kg    q'=14317 J/(kg.K)\n")

    plot_saturation(
        liquid, vapor, sat_wide, T_FIT,
        T_plot_range=T_PLOT,
        fluid_name=FLUID_NAME,
        savepath="nasg_water_saturation.pdf",
    )
    plt.show()
