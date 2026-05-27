"""
Compare NGA2 saturation outputs against IAPWS-IF97.

This version expects the simplified CSV files written by write_PTg_pT_curve:

    T,p,pv,Yv,VF,rhoL

where:
    T    : equilibrium temperature [K]
    p    : equilibrium mixture/mechanical pressure [Pa]
    pv   : vapor partial/saturation pressure [Pa]
    Yv   : vapor mass fraction in the gas phase [-]
    VF   : liquid volume fraction [-]
    rhoL : liquid density [kg/m^3]

The IAPWS reference curve is generated internally using the `iapws` Python
package; no separate IAPWS CSV file is read.

Install dependency if needed:
    pip install iapws
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
# User inputs
# =============================================================================

pure_path = Path("./NGA2_pT_pure_water_NASG.csv")
# mix_path = Path("./NGA2_pT_Yv0_0p5.csv")
mix_path = Path("./NGA2_pT_Yv0_0p5_NASG.csv")

pT_pdf_path = Path("./p_vs_T_nasg_pure_mix_comparison.pdf")
rhoT_pdf_path = Path("./rhoL_vs_T_nasg_pure_mix_comparison.pdf")


# =============================================================================
# Helpers
# =============================================================================

def read_ptg_csv(path):
    """Read a PTg_relax CSV and normalize old/new column names.

    New expected columns:
        T,p,pv,Yv,VF,rhoL

    Legacy columns with `_final` suffix are still accepted so older output files
    do not immediately break the plotting script.
    """
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = [col.strip() for col in df.columns]

    aliases = {
        "T":    ["T", "T_final"],
        "p":    ["p", "p_final", "Peq", "p_eq"],
        "pv":   ["pv", "pv_final", "ppv", "ppv_final"],
        "Yv":   ["Yv", "Yv_final"],
        "VF":   ["VF", "VF_final", "VFeq", "VF_eq"],
        "rhoL": ["rhoL", "rhoL_final", "rho_l", "rhoL_eq"],
    }

    normalized = pd.DataFrame(index=df.index)

    missing = []
    for clean_name, possible_names in aliases.items():
        for name in possible_names:
            if name in df.columns:
                normalized[clean_name] = pd.to_numeric(df[name], errors="coerce")
                break
        else:
            # Only T, pv, and rhoL are strictly needed for the current plots.
            if clean_name in ("T", "pv", "rhoL"):
                missing.append(clean_name)

    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {missing}. "
            "Expected at least T, pv, and rhoL. "
            "For the new test_relax output, the full header should be: "
            "T,p,pv,Yv,VF,rhoL"
        )

    if "valid_state" in df.columns:
        normalized = normalized[df["valid_state"] == 1]

    normalized = normalized.dropna(subset=["T", "pv", "rhoL"])

    if normalized.empty:
        raise ValueError(f"No valid rows found in {path}")

    return normalized


def iapws_saturation_table(T_min, T_max, n_points=300):
    """Build an IAPWS-IF97 saturation table.

    Returns a DataFrame with columns:
        T      : temperature [K]
        p_sat  : saturation pressure [Pa]
        rhoL   : saturated liquid density [kg/m^3]
    """
    T_crit = 647.096  # K
    T_low = max(float(T_min), 273.16)
    T_high = min(float(T_max), T_crit - 1.0e-6)

    if T_high <= T_low:
        raise ValueError(
            f"Invalid IAPWS temperature window: [{T_low}, {T_high}] K"
        )

    rows = []
    for Ti in np.linspace(T_low, T_high, n_points):
        try:
            liq = IAPWS97(T=Ti, x=0)
            rows.append({
                "T": Ti,
                "p_sat": liq.P * 1.0e6,   # IAPWS97 returns pressure in MPa
                "rhoL": 1.0 / liq.v,      # liq.v is m^3/kg
            })
        except Exception:
            # Skip states outside the package's valid region.
            pass

    ref = pd.DataFrame(rows)
    if ref.empty:
        raise RuntimeError(
            "IAPWS did not return any valid saturation states for the "
            f"temperature window [{T_low}, {T_high}] K."
        )

    return ref


def make_iapws_reference(pure, mix):
    """Use the PTg temperature range, padded slightly, for the reference curve."""
    Tmin = min(pure["T"].min(), mix["T"].min())
    Tmax = max(pure["T"].max(), mix["T"].max())

    # Pad the reference curve a bit so the line extends beyond the markers.
    pad = 0.05 * max(Tmax - Tmin, 1.0)
    return iapws_saturation_table(Tmin - pad, Tmax + pad)


# =============================================================================
# Main
# =============================================================================

def main():
    pure = read_ptg_csv(pure_path)
    mix = read_ptg_csv(mix_path)
    ref = make_iapws_reference(pure, mix)

    # -------------------------------------------------------------------------
    # Pressure-temperature plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))

    plt.semilogy(
        ref["T"],
        ref["p_sat"],
        linewidth=2,
        label="IAPWS-IF97"
    )

    plt.semilogy(
        pure["T"],
        pure["pv"],
        "o",
        markersize=6,
        label="NGA2 NASG pure water"
    )

    plt.semilogy(
        mix["T"],
        mix["pv"],
        "s",
        markersize=3,
        label="NGA2 NASG air-water"
    )

    plt.xlabel("Temperature, T [K]")
    plt.ylabel("Vapor pressure, pv [Pa]")
    plt.grid(True, which="both", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(pT_pdf_path)

    print("Saved:")
    print(pT_pdf_path)

    # -------------------------------------------------------------------------
    # Liquid-density-temperature plot
    # -------------------------------------------------------------------------
    plt.figure(figsize=(7.5, 5.0))

    plt.plot(
        ref["T"],
        ref["rhoL"],
        linewidth=2,
        label="IAPWS-IF97 saturated liquid"
    )

    plt.plot(
        pure["T"],
        pure["rhoL"],
        "o",
        markersize=6,
        label="NGA2 SG"
    )

    plt.plot(
        mix["T"],
        mix["rhoL"],
        "s",
        markersize=3,
        label="NGA2 NASG"
    )

    plt.xlabel("Temperature, T [K]")
    plt.ylabel("Liquid density, rhoL [kg/m^3]")
    plt.grid(True, which="both", alpha=0.35)
    plt.legend()
    plt.tight_layout()
    plt.savefig(rhoT_pdf_path)

    print("Saved:")
    print(rhoT_pdf_path)


if __name__ == "__main__":
    main()
