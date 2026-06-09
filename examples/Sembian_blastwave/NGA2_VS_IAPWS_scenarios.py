"""
Systematic verification of PTg relaxation across many initial conditions.
Produces 6 PDFs comparing NGA2 results against IAPWS-IF97.

1. Run: make -f Makefile.test_relax_scenarios run
2. Run: python NGA2_VS_IAPWS_scenarios.py

pip install iapws   (once, if missing)

PDFs produced:
  verification_VF_sweep.pdf      -- VF sweep (NASG and SG vs IAPWS)
  verification_Yv_sweep.pdf      -- Yv sweep (NASG and SG vs IAPWS)
  verification_liq_pressure.pdf  -- pure liquid nucleation across pressures (incl. cavitation)
  verification_gas_Yv.pdf        -- pure gas condensation, Yv sweep
  verification_gas_pressure.pdf  -- pure gas condensation, pressure sweep (Yv=0.5)
  verification_vap_pressure.pdf  -- pure vapor condensation, pressure sweep (Yv=1)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from iapws import IAPWS97

# ── Matplotlib defaults ────────────────────────────────────────────────────────

plt.rcParams.update({
    "text.usetex":       True,
    "font.family":       "serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "legend.fontsize":   8,
    "figure.dpi":        150,
    "axes.linewidth":    1.1,
    "xtick.major.size":  4,
    "ytick.major.size":  4,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.minor.size":  2,
    "ytick.minor.size":  2,
})

# (csv_col, iapws_col, y-label, log-scale)
PANELS = [
    ("pV",   "p_sat", r"$p_v\;[\mathrm{Pa}]$",              True ),
    ("hV",   "hg",    r"$h_v\;[\mathrm{J\,kg^{-1}}]$",      False),
    ("rhoL", "rhoL",  r"$\rho_L\;[\mathrm{kg\,m^{-3}}]$",   False),
    ("rhoV", "rhoG",  r"$\rho_v\;[\mathrm{kg\,m^{-3}}]$",   False),
]

# ── IAPWS-IF97 reference ──────────────────────────────────────────────────────

def iapws_sat(T_min=298.0, T_max=502.0, n=500):
    rows = []
    for T in np.linspace(max(T_min - 2, 273.16), min(T_max + 2, 647.09), n):
        try:
            liq = IAPWS97(T=T, x=0)
            vap = IAPWS97(T=T, x=1)
            rows.append({"T": T,
                          "p_sat": liq.P * 1e6,
                          "rhoL":  1.0 / liq.v,
                          "rhoG":  1.0 / vap.v,
                          "hg":    vap.h * 1e3})
        except Exception:
            pass
    return pd.DataFrame(rows)

# ── CSV loader ────────────────────────────────────────────────────────────────

def load(name):
    """Return DataFrame or None if file missing or empty."""
    if not os.path.exists(name):
        return None
    try:
        df = pd.read_csv(name)
        return df if len(df) > 0 else None
    except Exception:
        return None

# ── Core plotting helper ──────────────────────────────────────────────────────

def make_pdf(output, title, curves, ref, markevery=10):
    """
    curves: list of (csv_filename, label, color, marker, lw, ms, zorder)
    ref:    IAPWS DataFrame
    """
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(title, fontsize=11)

    for ax, (col, ref_col, ylabel, log) in zip(axes.flat, PANELS):
        ax.plot(ref["T"], ref[ref_col], "k-", lw=2.5, label="IAPWS-IF97", zorder=10)

        for (fname, label, color, marker, lw, ms, zo) in curves:
            df = load(fname)
            if df is None or col not in df.columns:
                continue
            plot_fn = ax.semilogy if log else ax.plot
            plot_fn(df["T"], df[col],
                    marker=marker, ms=ms, lw=lw,
                    color=color, markevery=markevery,
                    markerfacecolor=color, markeredgecolor=color,
                    label=label, zorder=zo)

        ax.set_xlabel(r"$T\;[\mathrm{K}]$")
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3, lw=0.6)
        ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(output)
    plt.close(fig)
    print(f"Saved {output}")

# ── Load reference ────────────────────────────────────────────────────────────

ref = iapws_sat()

# ── PDF 1: VF sweep ───────────────────────────────────────────────────────────
# All initial VF values should collapse onto the same saturation curve.

_vf_vals   = ["000","001","010","050","070","090","099"]
_vf_labels = [r"$\mathrm{VF}_0=0.00$",
               r"$\mathrm{VF}_0=0.01$",
               r"$\mathrm{VF}_0=0.10$",
               r"$\mathrm{VF}_0=0.50$",
               r"$\mathrm{VF}_0=0.70$",
               r"$\mathrm{VF}_0=0.90$",
               r"$\mathrm{VF}_0=0.99$"]
_vf_colors_n = [cm.Blues(v)   for v in np.linspace(0.35, 0.95, 7)]
_vf_colors_s = [cm.Reds(v)    for v in np.linspace(0.35, 0.95, 7)]
_vf_marks    = ["o","s","^","D","v","P","*"]

_vf_curves = []
for val, lbl, cn, cs, mk in zip(_vf_vals, _vf_labels, _vf_colors_n, _vf_colors_s, _vf_marks):
    _vf_curves.append((f"sc_NASG_VF{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _vf_curves.append((f"sc_SG_VF{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_VF_sweep.pdf",
    title=r"VF sweep ($Y_{v,0}=0.60$, $p_0=10^5\,\mathrm{Pa}$) --- all should overlap IAPWS",
    curves=_vf_curves,
    ref=ref,
)

# ── PDF 2: Yv sweep ───────────────────────────────────────────────────────────

_yv_vals   = ["000","001","030","090","099","100"]
_yv_labels = [r"$Y_{v,0}=0.00$",
               r"$Y_{v,0}=0.01$",
               r"$Y_{v,0}=0.30$",
               r"$Y_{v,0}=0.90$",
               r"$Y_{v,0}=0.99$",
               r"$Y_{v,0}=1.00$"]
_yv_colors_n = [cm.Blues(v)  for v in np.linspace(0.35, 0.95, 6)]
_yv_colors_s = [cm.Reds(v)   for v in np.linspace(0.35, 0.95, 6)]
_yv_marks    = ["o","s","^","D","v","P"]

_yv_curves = []
for val, lbl, cn, cs, mk in zip(_yv_vals, _yv_labels, _yv_colors_n, _yv_colors_s, _yv_marks):
    _yv_curves.append((f"sc_NASG_Yv{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _yv_curves.append((f"sc_SG_Yv{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_Yv_sweep.pdf",
    title=r"$Y_v$ sweep ($\mathrm{VF}_0=0.30$, $p_0=10^5\,\mathrm{Pa}$) --- all should overlap IAPWS",
    curves=_yv_curves,
    ref=ref,
)

# ── PDF 3: Pure liquid pressure sweep (nucleation + cavitation) ───────────────

_liq_vals   = ["pm1e6","pm1e5","pm1e4","pm1e2","p0","p1e4","p1e5","p1e6","p1e7"]
_liq_labels = [r"$p_0=-10^6\,\mathrm{Pa}$",
                r"$p_0=-10^5\,\mathrm{Pa}$",
                r"$p_0=-10^4\,\mathrm{Pa}$",
                r"$p_0=-10^2\,\mathrm{Pa}$",
                r"$p_0=0$",
                r"$p_0=10^4\,\mathrm{Pa}$",
                r"$p_0=10^5\,\mathrm{Pa}$",
                r"$p_0=10^6\,\mathrm{Pa}$",
                r"$p_0=10^7\,\mathrm{Pa}$"]
_liq_colors_n = [cm.Blues(v)  for v in np.linspace(0.25, 0.95, 9)]
_liq_colors_s = [cm.Reds(v)   for v in np.linspace(0.25, 0.95, 9)]
_liq_marks    = ["o","s","^","D","v","P","*","h","H"]

_liq_curves = []
for val, lbl, cn, cs, mk in zip(_liq_vals, _liq_labels, _liq_colors_n, _liq_colors_s, _liq_marks):
    _liq_curves.append((f"sc_NASG_liq_{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _liq_curves.append((f"sc_SG_liq_{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_liq_pressure.pdf",
    title=r"Pure liquid nucleation ($\mathrm{VF}_0=1$, $Y_{v,0}=0$) --- all should overlap IAPWS",
    curves=_liq_curves,
    ref=ref,
    markevery=20,
)

# ── PDF 4: Pure gas Yv sweep (condensation) ───────────────────────────────────

_gas_yv_vals   = ["000","001","010","050","090","099"]
_gas_yv_labels = [r"$Y_{v,0}=0.00$",
                   r"$Y_{v,0}=0.01$",
                   r"$Y_{v,0}=0.10$",
                   r"$Y_{v,0}=0.50$",
                   r"$Y_{v,0}=0.90$",
                   r"$Y_{v,0}=0.99$"]
_gas_yv_colors_n = [cm.Blues(v)  for v in np.linspace(0.35, 0.95, 6)]
_gas_yv_colors_s = [cm.Reds(v)   for v in np.linspace(0.35, 0.95, 6)]
_gas_yv_marks    = ["o","s","^","D","v","P"]

_gas_yv_curves = []
for val, lbl, cn, cs, mk in zip(_gas_yv_vals, _gas_yv_labels, _gas_yv_colors_n, _gas_yv_colors_s, _gas_yv_marks):
    _gas_yv_curves.append((f"sc_NASG_gas_Yv{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _gas_yv_curves.append((f"sc_SG_gas_Yv{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_gas_Yv.pdf",
    title=r"Pure gas condensation, $Y_v$ sweep ($\mathrm{VF}_0=0$, $p_0=10^5\,\mathrm{Pa}$)",
    curves=_gas_yv_curves,
    ref=ref,
)

# ── PDF 5: Pure gas pressure sweep (Yv=0.5) ───────────────────────────────────

_gas_p_vals   = ["1e3","1e4","1e5","1e6","1e7"]
_gas_p_labels = [r"$p_0=10^3\,\mathrm{Pa}$",
                  r"$p_0=10^4\,\mathrm{Pa}$",
                  r"$p_0=10^5\,\mathrm{Pa}$",
                  r"$p_0=10^6\,\mathrm{Pa}$",
                  r"$p_0=10^7\,\mathrm{Pa}$"]
_gas_p_colors_n = [cm.Blues(v) for v in np.linspace(0.35, 0.95, 5)]
_gas_p_colors_s = [cm.Reds(v)  for v in np.linspace(0.35, 0.95, 5)]
_gas_p_marks    = ["o","s","^","D","v"]

_gas_p_curves = []
for val, lbl, cn, cs, mk in zip(_gas_p_vals, _gas_p_labels, _gas_p_colors_n, _gas_p_colors_s, _gas_p_marks):
    _gas_p_curves.append((f"sc_NASG_gas_p{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _gas_p_curves.append((f"sc_SG_gas_p{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_gas_pressure.pdf",
    title=r"Pure gas condensation, pressure sweep ($\mathrm{VF}_0=0$, $Y_{v,0}=0.50$)",
    curves=_gas_p_curves,
    ref=ref,
)

# ── PDF 6: Pure vapor pressure sweep (Yv=1) ───────────────────────────────────

_vap_curves = []
for val, lbl, cn, cs, mk in zip(_gas_p_vals, _gas_p_labels, _gas_p_colors_n, _gas_p_colors_s, _gas_p_marks):
    _vap_curves.append((f"sc_NASG_vap_p{val}.csv", r"NASG " + lbl, cn, mk, 0, 5, 4))
    _vap_curves.append((f"sc_SG_vap_p{val}.csv",   r"SG "   + lbl, cs, mk, 0, 5, 3))

make_pdf(
    output="verification_vap_pressure.pdf",
    title=r"Pure vapor condensation, pressure sweep ($\mathrm{VF}_0=0$, $Y_{v,0}=1.00$)",
    curves=_vap_curves,
    ref=ref,
)
