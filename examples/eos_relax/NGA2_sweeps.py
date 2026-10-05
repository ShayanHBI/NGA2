"""
Five saturation-curve sweeps, each in the exact format of NGA2_VS_IAPWS.pdf, each
showing the new framework at nouter = 1, 3 and 5.

    make sweeps          (or: ./sweep_nouter <tag> ... for each, then python NGA2_sweeps.py)

The five differ only in the entry state and the Laplace jump, and that is the point:
when the cell starts at thermo-mechanical equilibrium the mechanical and thermal steps
have nothing to do on the first pass and the three curves lie on top of one another.
The further the entry state is from equilibrium, the further k=1 falls behind.

pip install iapws   (once, if missing)
"""

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97

# ── The five sweeps ───────────────────────────────────────────────────────────
#   tag, title, VF0, DP, DT, PJ, p0, Yv0
#   entry state is  p_l=p0(1+DP), p_g=p0(1-DP), T_l=T(1+DT), T_g=T(1-DT)

SWEEPS = [
    ("sweepA", "A. equilibrium start, no Laplace jump",
     0.50, 0.00, 0.00, 0.0,   1.0e5, 0.5),
    ("sweepB", "B. mechanical + thermal disequilibrium, no Laplace jump",
     0.50, 0.25, 0.10, 0.0,   1.0e5, 0.5),
    ("sweepC", "C. thin liquid, disequilibrium + Laplace jump",
     0.05, 0.25, 0.15, 1.0e4, 1.0e5, 0.5),
    ("sweepD", "D. thin gas, disequilibrium + Laplace jump",
     0.95, 0.25, 0.15, 1.0e4, 1.0e5, 0.5),
    ("sweepE", "E. high pressure, disequilibrium + Laplace jump",
     0.50, 0.20, 0.12, 5.0e4, 1.0e6, 0.5),
]

KS = [1, 3, 5]

PLOTS = [                                        # (csv_col, iapws_col, ylabel, log)
    ("pV",   "p_sat", r"$p\;[\mathrm{Pa}]$",                     True ),
    ("hV",   "hg",    r"$h_v\;[\mathrm{J\,kg^{-1}}]$",           False),
    ("rhoL", "rhoL",  r"$\rho_l\;[\mathrm{kg\,m^{-3}}]$",        False),
    ("rhoV", "rhoG",  r"$\rho_v\;[\mathrm{kg\,m^{-3}}]$",        False),
]

# ── Matplotlib defaults: identical to NGA2_VS_IAPWS.py ────────────────────────

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

# ── Generate the CSVs ─────────────────────────────────────────────────────────

for tag, _title, VF0, DP, DT, PJ, p0, Yv0 in SWEEPS:
    cmd = ["./sweep_nouter", tag, f"{VF0}", f"{DP}", f"{DT}", f"{PJ}", f"{p0}", f"{Yv0}"]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"{tag}: sweep_nouter failed\n{out.stdout}\n{out.stderr}")
    nrow = sum(1 for ln in out.stdout.splitlines() if "rows" in ln)
    print(f"{tag}: {out.stdout.splitlines()[0]}")

# ── One figure per sweep, in the NGA2_VS_IAPWS layout ─────────────────────────

markers = ["o", "^", "s"]
colors  = ["r", "g", "b"]

summary = []

for tag, title, VF0, DP, DT, PJ, p0, Yv0 in SWEEPS:

    data = {rf"new framework, $k={k}$": pd.read_csv(f"{tag}_k{k}.csv") for k in KS}
    T_all = np.concatenate([d["T"].values for d in data.values()])
    ref   = iapws_sat(T_all.min(), T_all.max())

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    for ax, (col, ref_col, ylabel, log) in zip(axes.flat, PLOTS):
        ax.plot(ref["T"], ref[ref_col], "k-", lw=2.5, label="IAPWS-IF97", zorder=3)

        # k=3 and k=5 coincide wherever the pass has converged, so interleave the
        # marker positions rather than drawing one on top of the other
        for i, ((label, df), mk, c) in enumerate(zip(data.items(), markers, colors)):
            (ax.semilogy if log else ax.plot)(
                df["T"], df[col], markevery=(3 * i, 9), marker=mk, ms=5, lw=0,
                markerfacecolor=c, markeredgecolor=c, markeredgewidth=1.5,
                label=label, zorder=4)

        ax.set_xlabel(r"$T\;[\mathrm{K}]$")
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3, lw=0.7)
        ax.legend(markerscale=1.5)

    sub = (rf"$V\!F_0={VF0:g}$,  $p_l/p_g={(1+DP)/(1-DP):.2f}$,  "
           rf"$T_l/T_g={(1+DT)/(1-DT):.2f}$,  "
           rf"$p_{{\mathrm{{jump}}}}={PJ:.0f}$ Pa,  $p_0={p0:.0e}$ Pa")
    fig.suptitle(title + "\n" + sub, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    out = f"NGA2_{tag}.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"Saved {out}")

    # how far k=1 and k=3 sit from k=5, as a number to go with the picture
    base = data[rf"new framework, $k=5$"]
    row = [title.split(".")[0]]
    for k in (1, 3):
        d = data[rf"new framework, $k={k}$"]
        dev = np.abs(d["pV"].values - base["pV"].values) / np.abs(base["pV"].values)
        row.append(dev.max())
    summary.append(row)

# ── Numbers for the record ────────────────────────────────────────────────────

print(f"\n{'sweep':<8}{'max |dp_v| k=1 vs k=5':>24}{'k=3 vs k=5':>16}")
for tag, d1, d3 in summary:
    print(f"{tag:<8}{d1:>24.2e}{d3:>16.2e}")
