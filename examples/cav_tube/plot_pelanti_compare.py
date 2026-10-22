"""Water cavitation tube: NGA2 final-time profiles against Pelanti (2022) Fig. 4, whose curves are extracted
from the vector paths of the paper's page 13 (pdftocairo -svg).
Usage: python3 plot_pelanti_compare.py --p | --pT | --pTg [--case amrviz dir] [--paper pdf]
Reads amrviz/cavtube_relax_<type> by default and writes cav_tube_<type>.pdf"""
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
ap = argparse.ArgumentParser(description="NGA2 water cavitation tube against Pelanti (2022) Fig. 4")
mode = ap.add_mutually_exclusive_group(required=True)
for r in ["p", "pT", "pTg"]:
    mode.add_argument(f"--{r}", dest="ref", action="store_const", const=r, help=f"{r}-relax run")
ap.add_argument("--case", type=Path, default=None, help="amrviz case dir (default: amrviz/cavtube_relax_<type>)")
ap.add_argument("--paper", type=Path, default=Path.home()/"Zotero/storage/RNRA9IMX/Pelanti - 2022 - Arbitrary-rate relaxation techniques for the numer.pdf")
args = ap.parse_args()
REF = args.ref
CASE = args.case if args.case is not None else HERE/"amrviz"/f"cavtube_relax_{REF}"
if not any(CASE.glob("plt.nga2.cell.*")): ap.error(f"no plotfiles in {CASE}")
PDF = args.paper
OUT = str(HERE/f"cav_tube_{REF}")

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


# ---------------- Pelanti (2022) Fig. 4: p-relax (violet), pT (dark green), pTg (blue), exact p-/pT-relaxed (black solid/dashed) ----------------
VIOLET = (69.8, 0.0, 69.8)
DGREEN = (0.0, 60.0, 0.0)
BLUE = (0.0, 0.0, 100.0)
BLACK = (12.7, 12.7, 12.7)
# Panel frames on page 13 (page coordinates) and y-value maps calibrated on each panel's tick marks;
# row 1 (p, u) holds all relaxation levels, row 2 (pX, uX) the exact solutions, row 3 alpha_2 and Y_2
FRAMES = {"p":  (129.35, 290.00, 65.87, 192.63, lambda y: (192.63-y)/24.146*2.0e4),
          "u":  (309.65, 470.30, 65.87, 192.63, lambda y: 2.0-(y-68.87)/15.091*0.5),
          "pX": (129.35, 290.00, 214.09, 340.84, lambda y: (340.84-y)/24.146*2.0e4),
          "uX": (309.65, 470.30, 214.09, 340.84, lambda y: 2.0-(y-217.09)/15.091*0.5),
          "a2": (134.03, 294.68, 365.21, 489.06, lambda y: (489.06-y)/17.697*0.05),
          "Y2": (310.37, 471.02, 364.70, 489.06, lambda y: (489.06-y)/19.137*0.2e-4)}


def pelanti_curves(pdf, colour, dashed_wanted=False):
    with tempfile.TemporaryDirectory() as tmp:
        svg = Path(tmp)/"p13.svg"
        subprocess.run(["pdftocairo", "-svg", "-f", "13", "-l", "13", str(pdf), str(svg)], check=True)
        text = svg.read_text()
    curves = {k: [] for k in FRAMES}
    for m in re.finditer(r"<path ([^>]*)/>", text):
        a = m.group(1)
        st = re.search(r'stroke="rgb\(([^)]*)\)"', a); d = re.search(r' d="([^"]*)"', a)
        if not st or not d: continue
        if tuple(round(float(c.strip(" %")), 1) for c in st.group(1).split(",")) != colour: continue
        if ("stroke-dasharray" in a) != dashed_wanted: continue
        # Black also draws frames and ticks: keep only long paths
        if colour == BLACK and d.group(1).count("L") < 50: continue
        T = [float(v) for v in re.search(r'transform="matrix\(([^)]*)\)"', a).group(1).split(",")]
        pts = np.array([(T[0]*float(x)+T[2]*float(y)+T[4], T[1]*float(x)+T[3]*float(y)+T[5])
                        for x, y in re.findall(r"(-?[\d.]+)\s+(-?[\d.]+)", d.group(1))])
        for k, (x0, x1, y0, y1, fy) in FRAMES.items():
            inside = (pts[:, 0] >= x0-0.5) & (pts[:, 0] <= x1+0.5) & (pts[:, 1] >= y0-0.5) & (pts[:, 1] <= y1+0.5)
            if inside.mean() < 0.9: continue
            xs = (pts[:, 0]-x0)/(x1-x0)
            # Legend samples: short strokes in a panel corner
            if np.ptp(xs) < 0.15: continue
            curves[k].append(np.c_[xs, fy(pts[:, 1])])
    # Dashed curves come as many short strokes: join them in x order; exact solutions live in row 2 for p and u
    for k in curves:
        if len(curves[k]) > 1:
            allp = np.vstack(curves[k]); curves[k] = [allp[np.argsort(allp[:, 0])]]
    for k in ["p", "u"]:
        if not curves[k]: curves[k] = curves[k+"X"]
    return curves


# ---------------- NGA2 final plotfile ----------------
def nga2_profiles(case):
    yt.set_log_level(40)
    f = sorted(Path(case).glob("plt.nga2.cell.*"), key=lambda p: int(p.name.split(".")[-1]))[-1]
    ds = yt.load(str(f)); n = int(ds.domain_dimensions[0])
    cg = ds.covering_grid(level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions)
    d = {k: np.array(cg["boxlib", k])[:, 0, 0] for k in ["VF", "RHOL", "RHOG", "PL", "PG", "U", "Yv"]}
    le, re_ = ds.domain_left_edge.to_value(), ds.domain_right_edge.to_value()
    x = le[0]+(np.arange(n)+0.5)*(re_[0]-le[0])/n
    VF = d["VF"]; aV = 1.0-VF
    rho = VF*d["RHOL"]+aV*d["RHOG"]
    q = {"p": VF*d["PL"]+aV*d["PG"], "u": d["U"], "a2": aV, "Y2": aV*d["RHOG"]*d["Yv"]/rho}
    return x, q, float(ds.current_time)


ref_colour = {"p": VIOLET, "pT": DGREEN, "pTg": BLUE}[REF]
ref = pelanti_curves(PDF, ref_colour)
# Exact solutions: p-relaxed (5-eq., solid black) and pT-relaxed (4-eq., dashed black)
exact = {"p": (pelanti_curves(PDF, BLACK), "Exact 5-eq."), "pT": (pelanti_curves(PDF, BLACK, dashed_wanted=True), "Exact 4-eq.")}.get(REF)
x, q, t = nga2_profiles(CASE)

# Log scale for alpha_v and Y_v: Pelanti's 1e-2/1e-4 levels and NGA2's pure-vapor cavity both stay visible
panels = [("p",  r"$p\;\left(\mathrm{Pa}\right)$", "linear", (0.0, 1.05e5)),
          ("u",  r"$u\;\left(\mathrm{m\,s^{-1}}\right)$", "linear", (-2.1, 2.1)),
          ("a2", r"$\alpha_v$", "log", (1.0e-3, 1.5)),
          ("Y2", r"$Y_v$", "log", (1.0e-7, 1.5))]
fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0), sharex=True)
for ax, (k, ylabel, scale, ylim) in zip(axes.ravel(), panels):
    if exact is not None:
        for i, seg in enumerate(exact[0][k]):
            # Exact solution as filled circles at uniform x so it stays visible under coinciding curves
            ax.plot(*exact_markers(seg, scale, ylim), "o", color="tab:blue", ms=5, mew=0, zorder=5, label=exact[1] if i == 0 else None)
    for i, seg in enumerate(ref[k]):
        ax.plot(seg[:, 0], seg[:, 1], "--", color="tab:red", lw=2.6, label="Pelanti" if i == 0 else None)
    ax.plot(x, q[k], "-", color="black", lw=2.2, label="NGA2")
    # Keep Pelanti's range but show NGA2 values beyond it (log axes keep their positive lower limit)
    ax.set_yscale(scale)
    ax.set_ylim(ylim[0] if scale == "log" else min(ylim[0], np.nanmin(q[k])), max(ylim[1], 1.05*np.nanmax(q[k])))
    ax.set_xlim(0.0, 1.0)
    ax.set_ylabel(ylabel)
    ax.tick_params(which="both", top=True, right=True, pad=2.5)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
for ax in axes[-1]:
    ax.set_xlabel(r"$x\;\left(\mathrm{m}\right)$", labelpad=2.5)
axes[0, 1].legend(fontsize=LEGEND_FONTSIZE, frameon=False, loc="upper left")
fig.suptitle(rf"Water cavitation tube, ${REF}$ relaxation, $t={t*1.0e3:.2f}\,\mathrm{{ms}}$", fontsize=TITLE_FONTSIZE)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT+".pdf")
print(f"wrote {OUT}.pdf (t={t*1e3:.3f} ms), Pelanti segments:", {k: len(v) for k, v in ref.items()},
      "exact:", {k: len(v) for k, v in exact[0].items()} if exact is not None else None)
