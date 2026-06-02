"""
Compare NGA2 test_relax saturation outputs against IAPWS-IF97.

Reads the CSV files written by test_relax:
    test_relax_pure_water.csv   – Yv=1.0 (pure liquid-vapor, no inert gas)
    test_relax_air_water.csv    – Yv=0.5 (air-water mixture)

CSV columns:
    T    : post-relaxation liquid temperature [K]
    p    : post-relaxation equilibrium pressure [Pa]
    pv   : vapor partial pressure [Pa]
    Yv   : post-relaxation vapor mass fraction [-]
    VF   : post-relaxation liquid volume fraction [-]
    rhoL : post-relaxation liquid density [kg/m^3]

The IAPWS-IF97 reference curve is generated using the `iapws` package.
Install with:  pip install iapws

Usage:
    python plot_relax.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    from iapws import IAPWS97
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: iapws\n"
        "Install it with:\n"
        "    pip install iapws"
    ) from exc


# =============================================================================
# Paths
# =============================================================================

pure_csv  = Path("test_relax_pure_water.csv")
mix_csv   = Path("test_relax_air_water.csv")

pT_pdf    = Path("plot_relax_pv_vs_T.pdf")
rhoT_pdf  = Path("plot_relax_rhoL_vs_T.pdf")
Yv_pdf    = Path("plot_relax_Yv_vs_T.pdf")
VF_pdf    = Path("plot_relax_VF_vs_T.pdf")


# =============================================================================
# Helper functions
# =============================================================================

def read_csv(path):
    """Read a test_relax CSV and return a clean DataFrame."""
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna()
    if df.empty:
        raise ValueError(f"No valid rows in {path}")
    return df


def iapws_saturation(T_min, T_max, n=400):
    """Return IAPWS-IF97 saturated-liquid table over [T_min, T_max] K."""
    T_crit = 647.096
    T_lo = max(float(T_min), 273.16)
    T_hi = min(float(T_max), T_crit - 1e-6)
    rows = []
    for Ti in np.linspace(T_lo, T_hi, n):
        try:
            liq = IAPWS97(T=Ti, x=0)
            rows.append({"T": Ti, "p_sat": liq.P * 1e6, "rhoL": 1.0 / liq.v})
        except Exception:
            pass
    ref = pd.DataFrame(rows)
    if ref.empty:
        raise RuntimeError("IAPWS returned no valid states in the given range.")
    return ref


# =============================================================================
# Main
# =============================================================================

def main():
    pure = read_csv(pure_csv)
    mix  = read_csv(mix_csv)

    T_all = np.concatenate([pure["T"].values, mix["T"].values])
    pad   = 0.05 * max(T_all.max() - T_all.min(), 1.0)
    ref   = iapws_saturation(T_all.min() - pad, T_all.max() + pad)

    # -----------------------------------------------------------------------
    # 1. Vapor pressure vs Temperature
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    ax.semilogy(ref["T"], ref["p_sat"],
                linewidth=2, color="k", label="IAPWS-IF97")
    ax.semilogy(pure["T"], pure["pv"],
                "o", markersize=5, label="NGA2 NASG – pure water (Yv=1)")
    ax.semilogy(mix["T"], mix["pv"],
                "s", markersize=4, label="NGA2 NASG – air+water (Yv=0.5)")

    ax.set_xlabel("Temperature  T  [K]")
    ax.set_ylabel("Vapor partial pressure  pv  [Pa]")
    ax.grid(True, which="both", alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(pT_pdf)
    print(f"Saved: {pT_pdf}")

    # -----------------------------------------------------------------------
    # 2. Liquid density vs Temperature
    # -----------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 5.0))

    ax.plot(ref["T"], ref["rhoL"],
            linewidth=2, color="k", label="IAPWS-IF97 saturated liquid")
    ax.plot(pure["T"], pure["rhoL"],
            "o", markersize=5, label="NGA2 NASG – pure water")
    ax.plot(mix["T"], mix["rhoL"],
            "s", markersize=4, label="NGA2 NASG – air+water")

    ax.set_xlabel("Temperature  T  [K]")
    ax.set_ylabel(r"Liquid density  $\rho_L$  [kg/m$^3$]")
    ax.grid(True, alpha=0.35)
    ax.legend()
    fig.tight_layout()
    fig.savefig(rhoT_pdf)
    print(f"Saved: {rhoT_pdf}")

    # -----------------------------------------------------------------------
    # 3. Vapor mass fraction vs Temperature
    # -----------------------------------------------------------------------
    if "Yv" in pure.columns and "Yv" in mix.columns:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))

        ax.plot(pure["T"], pure["Yv"],
                "o", markersize=5, label="NGA2 NASG – pure water")
        ax.plot(mix["T"], mix["Yv"],
                "s", markersize=4, label="NGA2 NASG – air+water")

        ax.set_xlabel("Temperature  T  [K]")
        ax.set_ylabel("Vapor mass fraction  Yv  [-]")
        ax.grid(True, alpha=0.35)
        ax.legend()
        fig.tight_layout()
        fig.savefig(Yv_pdf)
        print(f"Saved: {Yv_pdf}")

    # -----------------------------------------------------------------------
    # 4. Liquid volume fraction vs Temperature
    # -----------------------------------------------------------------------
    if "VF" in pure.columns and "VF" in mix.columns:
        fig, ax = plt.subplots(figsize=(7.5, 5.0))

        ax.plot(pure["T"], pure["VF"],
                "o", markersize=5, label="NGA2 NASG – pure water")
        ax.plot(mix["T"], mix["VF"],
                "s", markersize=4, label="NGA2 NASG – air+water")

        ax.set_xlabel("Temperature  T  [K]")
        ax.set_ylabel("Liquid volume fraction  VF  [-]")
        ax.grid(True, alpha=0.35)
        ax.legend()
        fig.tight_layout()
        fig.savefig(VF_pdf)
        print(f"Saved: {VF_pdf}")


if __name__ == "__main__":
    main()
