"""
NGA2 vs IAPWS-IF97 saturation curves.

1. Run ./test_relax to generate the CSVs.
2. python NGA2_VS_IAPWS.py

pip install iapws   (once, if missing)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97

# ── Configuration ─────────────────────────────────────────────────────────────

CURVES = {
    "NASG": "test_relax_NASG.csv",
    "SG":   "test_relax_SG.csv",
}

PLOTS = [                                        # (csv_col, iapws_col, ylabel, log)
    ("p",    "p_sat", r"$p\;[\mathrm{Pa}]$",                     True ),
    ("hG",   "hg",    r"$h_g\;[\mathrm{J\,kg^{-1}}]$",           False),
    ("rhoL", "rhoL",  r"$\rho_l\;[\mathrm{kg\,m^{-3}}]$",        False),
    ("rhoG", "rhoG",  r"$\rho_g\;[\mathrm{kg\,m^{-3}}]$",        False),
]

OUTPUT = "NGA2_VS_IAPWS.pdf"

# ── Matplotlib defaults ────────────────────────────────────────────────────────

plt.rcParams.update({
    "text.usetex":        True,
    "font.family":        "serif",
    "font.size":          13,
    "axes.labelsize":     14,
    "legend.fontsize":    12,
    "figure.dpi":         150,
    "axes.linewidth":     1.2,
    "xtick.major.size":   5,
    "ytick.major.size":   5,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
    "xtick.minor.size":   3,
    "ytick.minor.size":   3,
})

# ── IAPWS-IF97 reference ──────────────────────────────────────────────────────

def iapws_sat(T_min, T_max, n=500):
    rows = []
    for T in np.linspace(max(T_min - 2, 273.16), min(T_max + 2, 647.09), n):
        try:
            liq = IAPWS97(T=T, x=0)
            vap = IAPWS97(T=T, x=1)
            rows.append({"T": T, "p_sat": liq.P * 1e6,
                          "rhoL": 1 / liq.v, "rhoG": 1 / vap.v,
                          "hg":   vap.h * 1e3})
        except Exception:
            pass
    return pd.DataFrame(rows)

# ── Load data ─────────────────────────────────────────────────────────────────

data = {label: pd.read_csv(path) for label, path in CURVES.items()}
T_all = np.concatenate([d["T"].values for d in data.values()])
ref   = iapws_sat(T_all.min(), T_all.max())

# ── Plot ──────────────────────────────────────────────────────────────────────

markers = ["o", "s", "^", "D"]
colors  = ['b', 'r', 'c', 'm']
filled  = {"NASG": True, "SG": False}

fig, axes = plt.subplots(2, 2, figsize=(10, 7))

for ax, (col, ref_col, ylabel, log) in zip(axes.flat, PLOTS):
    # IAPWS reference
    ax.plot(ref["T"], ref[ref_col], "k-", lw=2.5, label="IAPWS-IF97", zorder=3)

    # NGA2 curves
    for (label, df), mk, col_ in zip(data.items(), markers, colors):
        if col not in df.columns:
            continue
        mfc = col_ if filled.get(label, True) else "none"
        (ax.semilogy if log else ax.plot)(
            df["T"], df[col], markevery=5, marker=mk, ms=5, lw=0,
            markerfacecolor=mfc, markeredgecolor=col_, markeredgewidth=1.5,
            label=label, zorder=4)

    ax.set_xlabel(r"$T\;[\mathrm{K}]$")
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.3, lw=0.7)
    ax.legend(markerscale=1.5)

fig.tight_layout()
fig.savefig(OUTPUT)
print(f"Saved {OUTPUT}")