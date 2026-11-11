"""
NGA2 vs IAPWS-IF97 saturation curves.

1. Run `make run` (or ./eos_relax) to generate the CSVs.
2. python NGA2_VS_IAPWS.py

Saturation sweeps of both max-ent inst classes (one page per sweep):
1. make sat
2. python NGA2_VS_IAPWS.py sweeps

pip install iapws   (once, if missing)
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from iapws import IAPWS97

# ── Configuration ─────────────────────────────────────────────────────────────

CURVES = {
    "old framework":  "eos_relax_NASG.csv",
    "new framework": "eos_relax_NASG_max_ent.csv",
    # "iterative SG":   "eos_relax_SG.csv",
}

PLOTS = [                                        # (csv_col, iapws_col, ylabel, log)
    ("pV",    "p_sat", r"$p\;[\mathrm{Pa}]$",                     True ),
    ("hV",   "hg",    r"$h_v\;[\mathrm{J\,kg^{-1}}]$",           False),
    ("rhoL", "rhoL",  r"$\rho_l\;[\mathrm{kg\,m^{-3}}]$",        False),
    ("rhoV", "rhoG",  r"$\rho_v\;[\mathrm{kg\,m^{-3}}]$",        False),
]

OUTPUT = "NGA2_VS_IAPWS.pdf"

# Sweeps of sat_sweep.f90: {tag} is replaced by the sweep tag in sat_sweeps.csv
SWEEP_CURVES = {
    "max-ent inst":      "sat_{tag}_inst.csv",
    "max-ent NASG inst": "sat_{tag}_nasg_inst.csv",
}
SWEEP_TABLE  = "sat_sweeps.csv"
MODEL_CURVE  = "sat_nasg_model.csv"        # NASG g_l = g_v, pure water, flat interface
SWEEP_OUTPUT = "NGA2_VS_IAPWS_sweeps.pdf"

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

# ── Plot ──────────────────────────────────────────────────────────────────────

markers = ["o", "^", "s"]
colors  = ['b', 'r', 'g']
filled  = {"old framework": True, "new framework": False, "iterative SG": False,
           "max-ent inst": True, "max-ent NASG inst": False}

def plot_figure(data, model=None, stagger=False):
    """IAPWS reference, optional NASG model curve, NGA2 curves; stagger offsets the markers of coinciding curves"""
    T_all = np.concatenate([d["T"].values for d in data.values()])
    ref   = iapws_sat(T_all.min(), T_all.max())

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))

    for ax, (col, ref_col, ylabel, log) in zip(axes.flat, PLOTS):
        # IAPWS reference
        ax.plot(ref["T"], ref[ref_col], "k-", lw=2.5, label="IAPWS-IF97", zorder=3)

        # NASG model saturation curve, restricted to the data range
        if model is not None and col in model.columns:
            m = model[(model["T"] >= T_all.min() - 2) & (model["T"] <= T_all.max() + 2)]
            ax.plot(m["T"], m[col], color="0.5", ls="--", lw=1.5, label="NASG model", zorder=3)

        # NGA2 curves
        for i, ((label, df), mk, col_) in enumerate(zip(data.items(), markers, colors)):
            if col not in df.columns:
                continue
            mfc = col_ if filled.get(label, True) else "none"
            (ax.semilogy if log else ax.plot)(
                df["T"], df[col], markevery=(4 * i, 8) if stagger else 8, marker=mk, ms=5, lw=0,
                markerfacecolor=mfc, markeredgecolor=col_, markeredgewidth=1.5,
                label=label, zorder=4)

        ax.set_xlabel(r"$T\;[\mathrm{K}]$")
        ax.set_ylabel(ylabel)
        if log:
            ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.3, lw=0.7)
        ax.legend(markerscale=1.5)

    return fig

if len(sys.argv) > 1 and sys.argv[1] == "sweeps":
    from matplotlib.backends.backend_pdf import PdfPages

    sweeps = pd.read_csv(SWEEP_TABLE)
    model  = pd.read_csv(MODEL_CURVE)
    with PdfPages(SWEEP_OUTPUT) as pdf:
        for _, s in sweeps.iterrows():
            data = {label: pd.read_csv(path.format(tag=s["tag"])) for label, path in SWEEP_CURVES.items()}
            fig = plot_figure(data, model=model, stagger=True)
            sub = (rf"$V\!F_0={s.VF0:g}$,  $p_l/p_g={(1+s.DP)/(1-s.DP):.2f}$,  "
                   rf"$T_l/T_g={(1+s.DT)/(1-s.DT):.2f}$,  $p_j={s.PJ:.0f}$ Pa,  "
                   rf"$p_0={s.p0:.0e}$ Pa,  $Y_{{v,0}}={s.Yv0:g}$")
            ok = (rf"converged: inst {int(s.ok_inst)}/{int(s.ok_inst + s.fail_inst)},  "
                  rf"NASG inst {int(s.ok_nasg_inst)}/{int(s.ok_nasg_inst + s.fail_nasg_inst)}")
            fig.suptitle(f"{s['tag']}. {s['title']}\n{sub}\n{ok}", fontsize=12)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved {SWEEP_OUTPUT}")
else:
    data = {label: pd.read_csv(path) for label, path in CURVES.items()}
    fig = plot_figure(data)
    fig.tight_layout()
    fig.savefig(OUTPUT)
    print(f"Saved {OUTPUT}")
