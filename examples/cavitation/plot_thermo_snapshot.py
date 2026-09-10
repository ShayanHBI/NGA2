import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter, MaxNLocator
import pyvista as pv
import yt

parser = argparse.ArgumentParser()
parser.add_argument("input", help="case identifier, e.g. wall_pTg")
parser.add_argument("time", type=float, help="simulation time in mus to render (closest available frame is used)")
args = parser.parse_args()

INPUT = args.input
TARGET_TIME = args.time
CASE = f"cavitation_{INPUT}"
AMRVIZ_DIR = Path("amrviz") / CASE

# PL/TL/RHOL are placeholder zeros where VF~0, PG/TG/RHOG where VF~1; those cells are left blank.
VF_EPS = 1.0e-6

# Plotted coordinates in mm, time in mus (raw plotfile/vtp data is in m and s).
M_TO_MM = 1.0e3
S_TO_MUS = 1.0e6

# (name, title, scale)
FIELDS = [
    ("TL", r"$T_l\;\left(\mathrm{K}\right)$", 1.0),
    ("TG", r"$T_g\;\left(\mathrm{K}\right)$", 1.0),
    ("PL", r"$p_l\;\left(\mathrm{MPa}\right)$", 1.0e-6),
    ("PG", r"$p_g\;\left(\mathrm{MPa}\right)$", 1.0e-6),
    ("RHOL", r"$\rho_l\;\left(\mathrm{kg\,m^{-3}}\right)$", 1.0),
    ("RHOG", r"$\rho_g\;\left(\mathrm{kg\,m^{-3}}\right)$", 1.0),
]
LIQUID_FIELDS = {"TL", "PL", "RHOL"}
GAS_FIELDS = {"TG", "PG", "RHOG"}

AUTO_VIEW_LIMITS = True      # zoom to the bubble's extent (from VF) at the rendered time; if False, use VIEW_XLIM/VIEW_YLIM
VIEW_MARGIN = 1.0            # [mm] padding around the bubble's bounding box, auto mode only
VIEW_MIN_HALF_SPAN = 1.0     # [mm] minimum half-span, so a tiny bubble doesn't zoom in absurdly tight
VIEW_SYMMETRIC = True        # square view centered on the domain center, auto mode only
BUBBLE_VF_THRESHOLD = 1.0e-3 # 1-VF below this is treated as noise when locating the bubble

# Fixed view window [mm], used only when AUTO_VIEW_LIMITS is False.
VIEW_XLIM = (-3.0, 3.0)
VIEW_YLIM = (-3.0, 3.0)

CMAP = "jet"
INTERFACE_COLOR = "black"
INTERFACE_LINEWIDTH = 0.3

# Paper-style figure geometry (inches). Each field gets its own panel and
# colorbar since TL/TG/PL/PG live on very different scales.
LEFT_MARGIN_IN   = 0.62
RIGHT_MARGIN_IN  = 0.6
BOTTOM_MARGIN_IN = 0.65
TOP_MARGIN_IN    = 0.75
ROW_GAP_IN       = 0.45
COL_GAP_IN       = 0.55
CBAR_GAP_IN      = 0.12
CBAR_WIDTH_IN    = 0.16
TICK_PAD_PT      = 5.0
TIME_TITLE_PAD_IN = 0.2  # distance from figure top edge down to the "t=..." text

# Panel size is fit to whichever of these two budgets is tighter.
MAX_FIG_WIDTH_IN  = 8.0
MAX_FIG_HEIGHT_IN = 8.5

AXES_LABEL_FONTSIZE = 13.0
NUMBER_FONTSIZE     = 10.5
TITLE_FONTSIZE      = 14.0
SUPTITLE_FONTSIZE   = 15.0

plt.rcParams.update({
    "text.usetex":         True,
    "font.family":         "serif",
    "font.size":           12.0,
    "axes.labelsize":      AXES_LABEL_FONTSIZE,
    "xtick.labelsize":     NUMBER_FONTSIZE,
    "ytick.labelsize":     NUMBER_FONTSIZE,
    "axes.linewidth":      1.3,
    "xtick.major.width":   1.1,
    "ytick.major.width":   1.1,
    "xtick.minor.width":   0.8,
    "ytick.minor.width":   0.8,
    "xtick.major.size":    4.0,
    "ytick.major.size":    4.0,
    "xtick.minor.size":    2.3,
    "ytick.minor.size":    2.3,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}",
})


def tick_formatter(value, _):
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    return rf"${value:.2f}$"


def cbar_tick_formatter(value, _):
    if abs(value) < 1.0e-8:
        return r"$0$"
    # .2g switches to scientific notation above 100, which would overflow into the neighboring panel
    if abs(value) >= 100:
        return rf"${value:.0f}$"
    return rf"${value:.2g}$"


def add_axes_in_inches(fig, left, bottom, width, height):
    fig_w, fig_h = fig.get_size_inches()
    return fig.add_axes([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])


def frame_number(path: Path):
    return int(re.search(r"(\d+)(?:\.vtp)?$", path.name).group(1))


def extract_plic_segments(vtp_path: Path):
    """Reduce the extruded PLIC quads (2D interface segments extruded in z)
    back down to 2D line segments in the x-y plane, in mm."""
    mesh = pv.read(str(vtp_path))
    if mesh.n_cells == 0:
        return np.empty((0, 2, 2))

    pts = mesh.points
    faces = mesh.faces
    segments = []
    idx = 0
    while idx < len(faces):
        n = faces[idx]
        cell_pt_idx = faces[idx + 1: idx + 1 + n]
        idx += 1 + n
        cell_pts = pts[cell_pt_idx]
        zmin = cell_pts[:, 2].min()
        bottom_xy = cell_pts[np.isclose(cell_pts[:, 2], zmin, atol=1.0e-9)][:, :2]
        if len(bottom_xy) >= 2:
            segments.append(bottom_xy[:2])
    segments = np.array(segments) if segments else np.empty((0, 2, 2))
    return segments * M_TO_MM


def closest_frame_index(frames, target_time):
    """Binary search on frame time (monotonic in frame number) so only a
    handful of plotfiles are opened instead of the whole series."""
    cache = {}

    def time_at(i):
        if i not in cache:
            cache[i] = float(yt.load(str(frames[i][0])).current_time.to_value()) * S_TO_MUS
        return cache[i]

    lo, hi = 0, len(frames) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if time_at(mid) < target_time:
            lo = mid + 1
        else:
            hi = mid

    candidates = [lo] + ([lo - 1] if lo > 0 else [])
    idx = min(candidates, key=lambda i: abs(time_at(i) - target_time))
    return idx, time_at(idx)


def compute_view_limits(vf, le, re):
    """Bounding box (with padding) of cells holding gas (1-VF > BUBBLE_VF_THRESHOLD)."""
    ny, nx = vf.shape
    dx = (re[0] - le[0]) / nx
    dy = (re[1] - le[1]) / ny

    mask = vf < 1.0 - BUBBLE_VF_THRESHOLD
    row_idx = np.where(mask.any(axis=1))[0]
    col_idx = np.where(mask.any(axis=0))[0]
    if len(row_idx) == 0 or len(col_idx) == 0:
        # No bubble (not yet cavitated, or fully collapsed) -- show the whole domain
        return (le[0], re[0]), (le[1], re[1])

    x_lo = le[0] + col_idx.min() * dx - VIEW_MARGIN
    x_hi = le[0] + (col_idx.max() + 1) * dx + VIEW_MARGIN
    y_lo = le[1] + row_idx.min() * dy - VIEW_MARGIN
    y_hi = le[1] + (row_idx.max() + 1) * dy + VIEW_MARGIN

    if VIEW_SYMMETRIC:
        xc = 0.5 * (le[0] + re[0])
        yc = 0.5 * (le[1] + re[1])
        half = max(xc - x_lo, x_hi - xc, yc - y_lo, y_hi - yc, VIEW_MIN_HALF_SPAN)
        x_lo, x_hi, y_lo, y_hi = xc - half, xc + half, yc - half, yc + half
    else:
        xc, xh = 0.5 * (x_lo + x_hi), max(0.5 * (x_hi - x_lo), VIEW_MIN_HALF_SPAN)
        yc, yh = 0.5 * (y_lo + y_hi), max(0.5 * (y_hi - y_lo), VIEW_MIN_HALF_SPAN)
        x_lo, x_hi, y_lo, y_hi = xc - xh, xc + xh, yc - yh, yc + yh

    # Clip to the actual domain so padding never requests space outside it.
    x_lo, x_hi = max(x_lo, le[0]), min(x_hi, re[0])
    y_lo, y_hi = max(y_lo, le[1]), min(y_hi, re[1])
    return (x_lo, x_hi), (y_lo, y_hi)


def load_frame(plt_path: Path, vtp_path: Path):
    ds = yt.load(str(plt_path))
    t = float(ds.current_time.to_value()) * S_TO_MUS

    # Sample at the finest AMR level so each pixel matches one leaf cell -- no interpolation.
    max_level = ds.index.max_level
    res = int(ds.domain_dimensions[0]) * 2**max_level

    le_m = ds.domain_left_edge.to_value("code_length")
    re_m = ds.domain_right_edge.to_value("code_length")

    slc = ds.slice("z", 0.0)
    frb = slc.to_frb((re_m[0] - le_m[0], "code_length"), res, height=(re_m[1] - le_m[1], "code_length"))

    vf = np.array(frb["boxlib", "VF"])
    data = {name: np.array(frb["boxlib", name]) * scale for name, _, scale in FIELDS}
    for name in LIQUID_FIELDS:
        data[name] = np.where(vf > VF_EPS, data[name], np.nan)
    for name in GAS_FIELDS:
        data[name] = np.where(vf < 1.0 - VF_EPS, data[name], np.nan)

    segments = extract_plic_segments(vtp_path)
    return t, data, vf, segments, le_m * M_TO_MM, re_m * M_TO_MM


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

frame_idx, _ = closest_frame_index(frames, TARGET_TIME)
t, data, vf, segments, le, re = load_frame(*frames[frame_idx])

if AUTO_VIEW_LIMITS:
    VIEW_XLIM, VIEW_YLIM = compute_view_limits(vf, le, re)

tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
x_ticks = [v for v in tick_locator.tick_values(*VIEW_XLIM) if VIEW_XLIM[0] <= v <= VIEW_XLIM[1]]
y_ticks = [v for v in tick_locator.tick_values(*VIEW_YLIM) if VIEW_YLIM[0] <= v <= VIEW_YLIM[1]]

extent = [le[0], re[0], le[1], re[1]]
data_height_over_width = (VIEW_YLIM[1] - VIEW_YLIM[0]) / (VIEW_XLIM[1] - VIEW_XLIM[0])

N_ROWS = 3

width_budget_in = (
    MAX_FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - COL_GAP_IN - 2 * (CBAR_GAP_IN + CBAR_WIDTH_IN)
) / 2
height_budget_in = (
    MAX_FIG_HEIGHT_IN - BOTTOM_MARGIN_IN - TOP_MARGIN_IN - (N_ROWS - 1) * ROW_GAP_IN
) / N_ROWS

PLOT_WIDTH_IN = min(width_budget_in, height_budget_in / data_height_over_width)
PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width

FIG_WIDTH_IN = LEFT_MARGIN_IN + RIGHT_MARGIN_IN + COL_GAP_IN + 2 * (CBAR_GAP_IN + CBAR_WIDTH_IN + PLOT_WIDTH_IN)
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + N_ROWS * PANEL_HEIGHT_IN + (N_ROWS - 1) * ROW_GAP_IN + TOP_MARGIN_IN

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

row_bottoms = {
    row: BOTTOM_MARGIN_IN + row * (PANEL_HEIGHT_IN + ROW_GAP_IN)
    for row in range(N_ROWS)
}
col_lefts = {
    0: LEFT_MARGIN_IN,
    1: LEFT_MARGIN_IN + PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN + COL_GAP_IN,
}

# (field, title) placed row-major top to bottom: TL,TG / PL,PG / RHOL,RHOG.
LAYOUT = [
    (FIELDS[0], 2, 0), (FIELDS[1], 2, 1),
    (FIELDS[2], 1, 0), (FIELDS[3], 1, 1),
    (FIELDS[4], 0, 0), (FIELDS[5], 0, 1),
]

for (name, label, _), row, col in LAYOUT:
    field_data = data[name]
    if np.isfinite(field_data).any():
        vmin, vmax = np.nanmin(field_data), np.nanmax(field_data)
    else:
        vmin, vmax = 0.0, 1.0

    left = col_lefts[col]
    bottom = row_bottoms[row]

    ax = add_axes_in_inches(fig, left, bottom, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
    cax = add_axes_in_inches(fig, left + PLOT_WIDTH_IN + CBAR_GAP_IN, bottom, CBAR_WIDTH_IN, PANEL_HEIGHT_IN)

    im = ax.imshow(field_data, origin="lower", extent=extent, cmap=CMAP, vmin=vmin, vmax=vmax, interpolation="none")
    ax.add_collection(LineCollection(segments, colors=INTERFACE_COLOR, linewidths=INTERFACE_LINEWIDTH))

    ax.set_xlim(*VIEW_XLIM)
    ax.set_ylim(*VIEW_YLIM)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.tick_params(which="both", top=True, right=True, pad=TICK_PAD_PT)

    if row == 0:
        ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_xlabel(r"$x\;\left(\mathrm{mm}\right)$", labelpad=4.0)
    else:
        ax.tick_params(labelbottom=False)

    if col == 0:
        ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_ylabel(r"$y\;\left(\mathrm{mm}\right)$", labelpad=2.0)
    else:
        ax.tick_params(labelleft=False)

    ax.set_title(label, fontsize=TITLE_FONTSIZE, pad=5.0)

    cbar = fig.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_ticks(np.linspace(vmin, vmax, 4))
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
    cbar.ax.tick_params(pad=3.0)

fig.text(
    0.5, 1.0 - TIME_TITLE_PAD_IN / FIG_HEIGHT_IN,
    rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$",
    ha="center", va="center", fontsize=SUPTITLE_FONTSIZE,
)

out = f"thermo_{INPUT}_t{t:.2f}.pdf"
fig.savefig(out)
print(f"Saved: {out}  (t={t:.4f} mus)")
