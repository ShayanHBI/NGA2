"""
Convergence of the new framework with the number of outer passes.

1. ./sweep_nouter [VF0] [DP] [DT] [PJ]      generates the CSVs
2. python NGA2_nouter.py

The sweep starts in mechanical, thermal AND chemical disequilibrium and runs with a
non-zero Laplace jump. That matters: if the cell starts at thermo-mechanical
equilibrium the mechanical and thermal steps have nothing to do on the first pass and
the nouter effect all but disappears.

Top row is the saturation curve itself. The solver error at k=1 is largely invisible
there, because p_v is on a log axis and because the NASG fit is itself 1-8% from IAPWS.
The bottom row is where the convergence actually shows: deviation of each pass count
from the new framework's own converged answer (nouter=12).

pip install iapws   (once, if missing)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97

# ── Configuration ─────────────────────────────────────────────────────────────

REFERENCE = "sweep_new_kref.csv"          # new framework, nouter=12

CURVES = [                                # (label, csv, colour, marker, filled)
    ("old framework",          "sweep_old.csv",    "k", "s", False),
    ("new framework, $k=1$",   "sweep_new_k1.csv", "r", "o", True ),
    ("new framework, $k=3$",   "sweep_new_k3.csv", "g", "^", True ),
    ("new framework, $k=5$",   "sweep_new_k5.csv", "b", "v", True ),
]

PLOTS = [                                 # (csv_col, iapws_col, ylabel, log)
    ("pV",   "p_sat", r"$p_v\;[\mathrm{Pa}]$",              True ),
    ("hV",   "hg",    r"$h_v\;[\mathrm{J\,kg^{-1}}]$",      False),
    ("rhoL", "rhoL",  r"$\rho_l\;[\mathrm{kg\,m^{-3}}]$",   False),
    ("rhoV", "rhoG",  r"$\rho_v\;[\mathrm{kg\,m^{-3}}]$",   False),
]

OUTPUT = "NGA2_nouter.pdf"

# ── Matplotlib defaults ───────────────────────────────────────────────────────

plt.rcParams.update({
    "text.usetex":        True,
    "font.family":        "serif",
    "font.size":          12,
    "axes.labelsize":     13,
    "legend.fontsize":    10,
    "figure.dpi":         150,
    "axes.linewidth":     1.2,
    "xtick.major.size":   5,
    "ytick.major.size":   5,
    "xtick.major.width":  1.2,
    "ytick.major.width":  1.2,
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

# ── Load ──────────────────────────────────────────────────────────────────────

conv = pd.read_csv(REFERENCE)
data = {label: pd.read_csv(path) for label, path, *_ in CURVES}
ref  = iapws_sat(conv["T"].min(), conv["T"].max())

# sweep_nouter writes only rows every curve relaxed, so the files are row-aligned
for label, df in data.items():
    if len(df) != len(conv):
        raise SystemExit(f"{label}: {len(df)} rows vs reference {len(conv)} -- rerun ./sweep_nouter")

FLOOR = 1e-16          # so that exact agreement still plots on a log axis

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))

for j, (col, ref_col, ylabel, log) in enumerate(PLOTS):

    # ── top: the curve itself ──
    ax = axes[0, j]
    ax.plot(ref["T"], ref[ref_col], "-", color="0.45", lw=3,
            label="IAPWS-IF97", zorder=2)
    for (label, _path, c, mk, filled) in CURVES:
        df = data[label]
        ax.plot(df["T"], df[col], marker=mk, ms=5, lw=0, markevery=9,
                markerfacecolor=(c if filled else "none"),
                markeredgecolor=c, markeredgewidth=1.4, label=label, zorder=4)
    if log:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r"$T\;[\mathrm{K}]$")
    ax.grid(True, which="both", alpha=0.3, lw=0.7)
    if j == 0:
        ax.legend(markerscale=1.4, loc="lower right")

    # ── bottom: deviation from the new framework's own converged answer ──
    ax = axes[1, j]
    base = conv[col].values
    for (label, _path, c, mk, filled) in CURVES:
        d = np.abs(data[label][col].values - base) / np.maximum(np.abs(base), 1e-300)
        ls = ":" if label == "old framework" else "-"
        ax.semilogy(conv["T"], np.maximum(d, FLOOR), ls, color=c, lw=1.8,
                    marker=mk, ms=4, markevery=9,
                    markerfacecolor=(c if filled else "none"),
                    markeredgecolor=c, label=label)
    ax.axhline(1e-2, color="0.6", lw=1, ls="--")
    ax.set_ylim(1e-16, 10)
    ax.set_ylabel(r"relative deviation from converged")
    ax.set_xlabel(r"$T\;[\mathrm{K}]$")
    ax.grid(True, which="both", alpha=0.3, lw=0.7)
    ax.set_title(ylabel, fontsize=12)
    if j == 0:
        ax.legend(markerscale=1.3, loc="lower right", framealpha=0.9)

# the figure labels itself with the sweep definition, since those are now CLI knobs
try:
    pr = pd.read_csv("sweep_params.csv").iloc[0]
    sub = (rf"$V\!F_0={pr.VF0:g}$,  "
           rf"$p_l/p_g={((1+pr.DP)/(1-pr.DP)):.2f}$,  "
           rf"$T_l/T_g={((1+pr.DT)/(1-pr.DT)):.2f}$,  "
           rf"$p_{{\mathrm{{jump}}}}={pr.PJ:.0f}$ Pa,  "
           rf"$Y_{{v,0}}={pr.Yv0:g}$")
    if pr.nskip:
        sub += rf";  {int(pr.nskip)} temperatures declined and are not shown"
except Exception:
    sub = r"cell started in mechanical, thermal and chemical disequilibrium"
fig.suptitle(r"Convergence with outer passes $k$" + "\n" + sub, fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(OUTPUT)
print(f"Saved {OUTPUT}")

# ── Numbers for the record ────────────────────────────────────────────────────

print(f"\n{'curve':<26}" + "".join(f"{c:>12}" for c, *_ in PLOTS))
for (label, _path, *_rest) in CURVES:
    plain = label.replace("$", "").replace("\\", "")
    line = f"{plain:<26}"
    for col, *_ in PLOTS:
        base = conv[col].values
        d = np.abs(data[label][col].values - base) / np.maximum(np.abs(base), 1e-300)
        line += f"{d.max():>12.2e}"
    print(line)
print("\n(max relative deviation from the new framework at nouter=12)")
