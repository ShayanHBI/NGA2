"""Dodecane shock tube: NGA2 final-time profiles against Pelanti (2022) Fig. 6, whose curves are extracted
from the vector paths of the paper's page 15 (pdftocairo -svg).
Usage: python3 plot_pelanti_compare.py --p | --pT | --pTg [--case amrviz dir] [--paper pdf]
Reads amrviz/dodecane_relax_<type> by default and writes dodecane_shocktube_<type>.pdf"""
import argparse
import re
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import yt

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description="NGA2 dodecane shock tube against Pelanti (2022) Fig. 6")
mode = ap.add_mutually_exclusive_group(required=True)
for r in ["p", "pT", "pTg"]:
    mode.add_argument(f"--{r}", dest="ref", action="store_const", const=r, help=f"{r}-relax run")
ap.add_argument("--case", type=Path, default=None, help="amrviz case dir (default: amrviz/dodecane_relax_<type>)")
ap.add_argument("--paper", type=Path, default=Path.home()/"Zotero/storage/RNRA9IMX/Pelanti - 2022 - Arbitrary-rate relaxation techniques for the numer.pdf")
args = ap.parse_args()
# Fig. 6 has p and pTg; pT is shown against both
REF = args.ref
CASE = args.case if args.case is not None else HERE/"amrviz"/f"dodecane_relax_{REF}"
if not any(CASE.glob("plt.nga2.cell.*")): ap.error(f"no plotfiles in {CASE}")
PDF = args.paper
OUT = str(HERE/f"dodecane_shocktube_{REF}")

AXES_LABEL_FONTSIZE = 15
NUMBER_FONTSIZE     = 12
LEGEND_FONTSIZE     = 12
TITLE_FONTSIZE      = 15
plt.rcParams.update({
    "text.usetex":         True,
    "font.family":         "serif",
    "font.size":           12.0,
    "axes.labelsize":      AXES_LABEL_FONTSIZE,
    "xtick.labelsize":     NUMBER_FONTSIZE,
    "ytick.labelsize":     NUMBER_FONTSIZE,
    "axes.linewidth":      1.3,
    "xtick.major.width":   1.2,
    "ytick.major.width":   1.2,
    "xtick.minor.width":   0.9,
    "ytick.minor.width":   0.9,
    "xtick.major.size":    4.5,
    "ytick.major.size":    4.5,
    "xtick.minor.size":    2.6,
    "ytick.minor.size":    2.6,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}",
})


def tick_formatter(value, _):
    if abs(value-round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return rf"${text}$"


def exact_markers(seg, scale, ylim, n=40):
    # Uniform-x marker positions, keeping every third (and the last) marker inside flat runs
    s = seg[np.argsort(seg[:, 0])]; xm = np.linspace(s[0, 0], s[-1, 0], n); ym = np.interp(xm, s[:, 0], s[:, 1])
    f = np.log10 if scale == "log" else (lambda v: v)
    dy = np.abs(np.diff(f(np.maximum(ym, ylim[0]))))/(f(ylim[1])-f(ylim[0]))
    flat = np.r_[True, dy < 5.0e-3] & np.r_[dy < 5.0e-3, True]
    keep = ~flat
    for run in np.split(np.arange(n), np.flatnonzero(np.diff(flat.astype(int)))+1):
        if flat[run[0]]: keep[run[::3]] = True; keep[run[-1]] = True
    return xm[keep], ym[keep]


# ---------------- Pelanti (2022) Fig. 6: pTg-relax (blue), p-relax (violet), exact 5-eq p-relaxed (black, rho/u/p) ----------------
BLUE = (0.0, 0.0, 100.0)
VIOLET = (69.8, 0.0, 69.8)
BLACK = (12.7, 12.7, 12.7)
# Panel frames on page 15 (page coordinates) and y-value maps calibrated on each panel's tick marks
FRAMES = {"rho": (122.2, 292.1, 66.9, 200.5, lambda y: 10**(3-(y-66.91)/44.533)),
          "u":   (310.5, 480.4, 66.4, 200.5, lambda y: (200.51-y)/134.08*300),
          "p":   (309.5, 479.5, 223.1, 357.2, lambda y: 10**(5+(355.29-y)/41.617)),
          "Y2":  (122.39, 292.32, 225.68, 357.22, lambda y: (357.22-y)/119.61),
          "a2":  (119.9, 289.8, 382.9, 513.9, lambda y: (513.93-y)/119.1),
          "T":   (312.0, 481.9, 379.9, 513.9, lambda y: 300+(513.93-y)/134.07*900)}


def pelanti_curves(pdf, colour):
    with tempfile.TemporaryDirectory() as tmp:
        svg = Path(tmp)/"p15.svg"
        subprocess.run(["pdftocairo", "-svg", "-f", "15", "-l", "15", str(pdf), str(svg)], check=True)
        text = svg.read_text()
    curves = {k: [] for k in ["rho", "u", "p", "Y2", "a2", "TL", "TG"]}
    for m in re.finditer(r"<path ([^>]*)/>", text):
        a = m.group(1)
        st = re.search(r'stroke="rgb\(([^)]*)\)"', a); d = re.search(r' d="([^"]*)"', a)
        if not st or not d: continue
        if tuple(round(float(c.strip(" %")), 1) for c in st.group(1).split(",")) != colour: continue
        # Thin blue lines mark regions where that phase is almost absent: not compared
        if colour == BLUE and float(re.search(r'stroke-width="([^"]*)"', a).group(1)) < 0.5: continue
        dashed = "stroke-dasharray" in a
        # Black also draws frames, ticks and dashed guides: keep only long solid paths (and the exact Y2/alpha2 steps, below)
        if colour == BLACK and dashed: continue
        short = colour == BLACK and d.group(1).count("L") < 50
        T = [float(v) for v in re.search(r'transform="matrix\(([^)]*)\)"', a).group(1).split(",")]
        pts = np.array([(T[0]*float(x)+T[2]*float(y)+T[4], T[1]*float(x)+T[3]*float(y)+T[5])
                        for x, y in re.findall(r"(-?[\d.]+)\s+(-?[\d.]+)", d.group(1))])
        for k, (x0, x1, y0, y1, fy) in FRAMES.items():
            inside = (pts[:, 0] >= x0-0.5) & (pts[:, 0] <= x1+0.5) & (pts[:, 1] >= y0-0.5) & (pts[:, 1] <= y1+0.5)
            if inside.mean() < 0.9: continue
            interior = (pts[:, 0] > x0+3) & (pts[:, 0] < x1-3) & (pts[:, 1] > y0+3) & (pts[:, 1] < y1-3)
            if short and (k not in ["Y2", "a2"] or not interior.any()): continue
            xs = (pts[:, 0]-x0)/(x1-x0)
            # Legend samples: short strokes in a panel corner
            if xs.ptp() < 0.15 and k != "T": continue
            if k == "T" and xs.max() < 0.15: continue
            key = k if k != "T" else ("TG" if dashed else "TL")
            curves[key].append(np.c_[xs, fy(pts[:, 1])])
    # Dashed vapor temperature comes as many short strokes: join them in x order
    for key in ["TL", "TG"]:
        if curves[key]:
            allp = np.vstack(curves[key]); curves[key] = [allp[np.argsort(allp[:, 0])]]
    return curves


# ---------------- NGA2 final plotfile ----------------
def nga2_profiles(case):
    yt.set_log_level(40)
    f = sorted(Path(case).glob("plt.nga2.cell.*"), key=lambda p: int(p.name.split(".")[-1]))[-1]
    ds = yt.load(str(f)); n = int(ds.domain_dimensions[0])
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    d = {k: np.array(cg["boxlib", k])[:, 0, 0] for k in ["VF", "RHOL", "RHOG", "PL", "PG", "TL", "TG", "U", "Yv"]}
    le, re_ = ds.domain_left_edge.to_value(), ds.domain_right_edge.to_value()
    x = le[0]+(np.arange(n)+0.5)*(re_[0]-le[0])/n
    VF = d["VF"]; aV = 1.0-VF
    rho = VF*d["RHOL"]+aV*d["RHOG"]
    q = {"rho": rho, "u": d["U"], "a2": aV, "p": VF*d["PL"]+aV*d["PG"], "Y2": aV*d["RHOG"]*d["Yv"]/rho,
         "TL": np.where(VF > 1.0e-12, d["TL"], np.nan), "TG": np.where(aV > 1.0e-12, d["TG"], np.nan)}
    return x, q, float(ds.current_time)


def mixture_T(a2, TL, TG):
    # Volume-fraction weighted mixture temperature; where one phase is absent (nan) the other's T is used
    T = (1.0-a2)*TL+a2*TG
    return np.where(np.isnan(TL), TG, np.where(np.isnan(TG), TL, T))


def pelanti_T(curves):
    xs = np.linspace(0.0, 1.0, 2001)
    def interp(k):
        if not curves[k]: return np.full_like(xs, np.nan)
        s = np.vstack(curves[k]); s = s[np.argsort(s[:, 0])]
        return np.interp(xs, s[:, 0], s[:, 1], left=np.nan, right=np.nan)
    T = mixture_T(interp("a2"), interp("TL"), interp("TG"))
    return [np.c_[xs, T][~np.isnan(T)]]


refs = {"p": [("p", VIOLET)], "pT": [("p", VIOLET), ("pTg", BLUE)], "pTg": [("pTg", BLUE)]}[REF]
refs = [(lab, pelanti_curves(PDF, col)) for lab, col in refs]
for _, ref in refs: ref["T"] = pelanti_T(ref)
exact = pelanti_curves(PDF, BLACK) if REF == "p" else {}
x, q, t = nga2_profiles(CASE)
q["T"] = mixture_T(q["a2"], q["TL"], q["TG"])

panels = [("rho", r"$\rho\;\left(\mathrm{kg\,m^{-3}}\right)$", "log", (1.0, 1.0e3)),
          ("u",   r"$u\;\left(\mathrm{m\,s^{-1}}\right)$", "linear", (0.0, 340.0)),
          ("Y2",  r"$Y_v$", "linear", (-0.02, 1.05)),
          ("p",   r"$p\;\left(\mathrm{Pa}\right)$", "log", (5.0e4, 2.0e8)),
          ("a2",  r"$\alpha_v$", "linear", (-0.02, 1.05)),
          ("T",   r"$T\;\left(\mathrm{K}\right)$", "linear", (500.0, 1150.0))]
fig, axes = plt.subplots(3, 2, figsize=(9.0, 10.0), sharex=True)
axes = axes.ravel()
for ax, (k, ylabel, scale, ylim) in zip(axes, panels):
    for i, seg in enumerate(exact.get(k, [])):
        # Exact solution as filled circles at uniform x so it stays visible under coinciding curves
        ax.plot(*exact_markers(seg, scale, ylim), "o", color="tab:blue", ms=5, mew=0, zorder=5, label="Exact (5-eq.)" if i == 0 else None)
    for (lab, ref), col in zip(refs, ["tab:red", "tab:blue"] if len(refs) > 1 else ["tab:red"]):
        for i, seg in enumerate(ref.get(k, [])):
            ax.plot(seg[:, 0], seg[:, 1], "--", color=col, lw=2.6, label=(f"Pelanti, ${lab}$" if len(refs) > 1 else "Pelanti") if i == 0 else None)
    ax.plot(x, q[k], "-", color="black", lw=2.2, label="NGA2")
    ax.set_yscale(scale); ax.set_ylim(*ylim); ax.set_xlim(0.0, 1.0)
    ax.set_ylabel(ylabel)
    ax.tick_params(which="both", top=True, right=True, pad=2.5)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
    if scale == "linear": ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
for ax in axes[-2:]:
    ax.set_xlabel(r"$x\;\left(\mathrm{m}\right)$", labelpad=2.5)
axes[0].legend(fontsize=LEGEND_FONTSIZE, frameon=False, loc="lower left")
fig.suptitle(rf"Dodecane shock tube, ${REF}$ relaxation, $t={t*1.0e6:.0f}\,\mu\mathrm{{s}}$", fontsize=TITLE_FONTSIZE)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(OUT+".pdf")
print(f"wrote {OUT}.pdf (t={t*1e6:.1f} us), Pelanti segments:", {lab: {k: len(v) for k, v in ref.items()} for lab, ref in refs})
