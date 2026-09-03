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
parser.add_argument("input", help="case identifier, e.g. pTg")
parser.add_argument("time", type=float, help="simulation time to render (closest available frame is used)")
args = parser.parse_args()

INPUT = args.input
TARGET_TIME = args.time
CASE = f"impact_relax_{INPUT}"
AMRVIZ_DIR = Path("amrviz new") / CASE

FIELD = "cluster_idx"
CBAR_LABEL = r"$\mathrm{cluster\;id}$"

# Case is non-dimensional -- plotted coordinates and time are used as-is.
VIEW_XLIM = (0.0, 3.5)
VIEW_YLIM = (-6, 6)

# Rotate the whole view 90 deg counterclockwise: simulation x runs vertically
# (bottom to top) and simulation y runs horizontally (right to left, since a
# true CCW rotation sends +y to the left). Also flips a tall/narrow view into
# a short/wide one, which is what keeps the panel a reasonable shape here.
ROTATE_CCW = True

CMAP = "jet"
# Most of the domain is unpooled (masked to blank/white), unlike the other
# field plots where color fills the whole domain -- black shows up on both
# the blank background and the colored cluster regions.
INTERFACE_COLOR = "black"
INTERFACE_LINEWIDTH = 0.1

# Paper-style figure geometry (inches).
LEFT_MARGIN_IN   = 0.62
RIGHT_MARGIN_IN  = 0.6
BOTTOM_MARGIN_IN = 0.65
TOP_MARGIN_IN    = 0.45
CBAR_GAP_IN      = 0.13
CBAR_WIDTH_IN    = 0.17
TICK_PAD_PT      = 5.0

# Panel size is fit to whichever of these two budgets is tighter, so a very
# tall/narrow VIEW_XLIM/YLIM (large y_span/x_span, or the reverse after
# rotation) shrinks the panel instead of blowing up the figure size.
MAX_FIG_WIDTH_IN  = 8.0
MAX_FIG_HEIGHT_IN = 8.5

AXES_LABEL_FONTSIZE = 16.0
NUMBER_FONTSIZE     = 13.0
CBAR_LABEL_FONTSIZE = 14.0
TITLE_FONTSIZE      = 15.0

plt.rcParams.update({
    "text.usetex":         True,
    "font.family":         "serif",
    "font.size":           14.0,
    "axes.labelsize":      AXES_LABEL_FONTSIZE,
    "xtick.labelsize":     NUMBER_FONTSIZE,
    "ytick.labelsize":     NUMBER_FONTSIZE,
    "axes.linewidth":      1.5,
    "xtick.major.width":   1.35,
    "ytick.major.width":   1.35,
    "xtick.minor.width":   1.0,
    "ytick.minor.width":   1.0,
    "xtick.major.size":    5.0,
    "ytick.major.size":    5.0,
    "xtick.minor.size":    3.0,
    "ytick.minor.size":    3.0,
    "xtick.direction":     "in",
    "ytick.direction":     "in",
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}",
})


def tick_formatter(value, _):
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    return rf"${value:.2f}$"


def cbar_tick_formatter(value, _):
    return rf"${int(round(value))}$"


def add_axes_in_inches(fig, left, bottom, width, height):
    fig_w, fig_h = fig.get_size_inches()
    return fig.add_axes([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])


def frame_number(path: Path):
    return int(re.search(r"(\d+)(?:\.vtp)?$", path.name).group(1))


def extract_plic_segments(vtp_path: Path):
    """Reduce the extruded PLIC quads (2D interface segments extruded in z)
    back down to 2D line segments in the x-y plane."""
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
    return segments


def closest_frame_index(frames, target_time):
    """Binary search on frame time (monotonic in frame number) so only a
    handful of plotfiles are opened instead of the whole series."""
    cache = {}

    def time_at(i):
        if i not in cache:
            cache[i] = float(yt.load(str(frames[i][0])).current_time.to_value())
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


def load_frame(plt_path: Path, vtp_path: Path):
    ds = yt.load(str(plt_path))
    t = float(ds.current_time.to_value())

    max_level = ds.index.max_level
    res = int(ds.domain_dimensions[0]) * 2**max_level

    le = ds.domain_left_edge.to_value("code_length")
    re = ds.domain_right_edge.to_value("code_length")

    slc = ds.slice("z", 0.0)
    frb = slc.to_frb((re[0] - le[0], "code_length"), res, height=(re[1] - le[1], "code_length"))

    # 0 = unpooled cell (see amrmpcomp_class.f90); mask it out so only cells
    # that were actually part of a relaxation cluster get colored.
    field_data = np.array(frb["boxlib", FIELD])
    field_data = np.where(field_data != 0.0, field_data, np.nan)

    segments = extract_plic_segments(vtp_path)
    return t, field_data, segments, le, re


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

frame_idx, _ = closest_frame_index(frames, TARGET_TIME)
t, field_data, segments, le, re = load_frame(*frames[frame_idx])

tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
x_ticks = [v for v in tick_locator.tick_values(*VIEW_XLIM) if VIEW_XLIM[0] <= v <= VIEW_XLIM[1]]
y_ticks = [v for v in tick_locator.tick_values(*VIEW_YLIM) if VIEW_YLIM[0] <= v <= VIEW_YLIM[1]]

if ROTATE_CCW:
    # (x,y) -> (-y,x); implemented as transpose (row<->col) + set_xlim/invert
    # so the displayed x-axis carries y-values without negating them.
    field_data = field_data.T
    segments = segments[..., ::-1]
    extent = [le[1], re[1], le[0], re[0]]
    plot_xlim, plot_ylim = VIEW_YLIM, VIEW_XLIM
    plot_xticks, plot_yticks = y_ticks, x_ticks
    plot_xlabel, plot_ylabel = r"$y$", r"$x$"
    data_height_over_width = (VIEW_XLIM[1] - VIEW_XLIM[0]) / (VIEW_YLIM[1] - VIEW_YLIM[0])
else:
    extent = [le[0], re[0], le[1], re[1]]
    plot_xlim, plot_ylim = VIEW_XLIM, VIEW_YLIM
    plot_xticks, plot_yticks = x_ticks, y_ticks
    plot_xlabel, plot_ylabel = r"$x$", r"$y$"
    data_height_over_width = (VIEW_YLIM[1] - VIEW_YLIM[0]) / (VIEW_XLIM[1] - VIEW_XLIM[0])

width_budget_in = MAX_FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
height_budget_in = MAX_FIG_HEIGHT_IN - BOTTOM_MARGIN_IN - TOP_MARGIN_IN

PLOT_WIDTH_IN = min(width_budget_in, height_budget_in / data_height_over_width)
PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width

FIG_WIDTH_IN = LEFT_MARGIN_IN + PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN + RIGHT_MARGIN_IN
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + PANEL_HEIGHT_IN + TOP_MARGIN_IN

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
ax = add_axes_in_inches(fig, LEFT_MARGIN_IN, BOTTOM_MARGIN_IN, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
cax = add_axes_in_inches(
    fig, LEFT_MARGIN_IN + PLOT_WIDTH_IN + CBAR_GAP_IN, BOTTOM_MARGIN_IN, CBAR_WIDTH_IN, PANEL_HEIGHT_IN
)

vmin = np.nanmin(field_data) if np.isfinite(field_data).any() else 0.0
vmax = np.nanmax(field_data) if np.isfinite(field_data).any() else 1.0

im = ax.imshow(field_data, origin="lower", extent=extent, cmap=CMAP, vmin=vmin, vmax=vmax, interpolation="none")
ax.add_collection(LineCollection(segments, colors=INTERFACE_COLOR, linewidths=INTERFACE_LINEWIDTH))

ax.set_xlim(*plot_xlim)
ax.set_ylim(*plot_ylim)
if ROTATE_CCW:
    ax.invert_xaxis()
ax.set_aspect("equal", adjustable="box")

ax.set_xticks(plot_xticks)
ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
ax.set_yticks(plot_yticks)
ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
ax.tick_params(which="both", top=True, right=True, pad=TICK_PAD_PT)

ax.set_xlabel(plot_xlabel, labelpad=5.0)
ax.set_ylabel(plot_ylabel, labelpad=5.0)
ax.set_title(rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}$", fontsize=TITLE_FONTSIZE, pad=6.0)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.locator = MaxNLocator(integer=True)
cbar.update_ticks()
cbar.ax.tick_params(pad=3.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=8.0, fontsize=CBAR_LABEL_FONTSIZE)

out = f"cluster_idx_{INPUT}_t{t:.2f}.pdf"
fig.savefig(out)
print(f"Saved: {out}  (t={t:.4f})")
