import argparse
import re
from functools import lru_cache
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.animation import FuncAnimation
import pyvista as pv
import yt

parser = argparse.ArgumentParser()
parser.add_argument("--pdf", action="store_true",
                     help="Only render the first frame (t=0) to a PDF for layout tuning, skip the full video.")
args = parser.parse_args()

# Cases to compare: top row / bottom row.
TOP_INPUT, TOP_LABEL = "pThybrid", "pThybrid"
BOTTOM_INPUT, BOTTOM_LABEL = "pTg", "pTg"

FIELD = "VF"
CBAR_LABEL = r"$\alpha$"
VMIN, VMAX = 0.0, 1.0

# Case is non-dimensional -- plotted coordinates and time are used as-is.
VIEW_XLIM = (0.0, 11.0)
VIEW_YLIM = (-2.5, 2.5)

CMAP = "jet"
INTERFACE_COLOR = "white"
INTERFACE_LINEWIDTH = 0.6

FPS = 12
OUTPUT = f"VF_plic_compare_{TOP_INPUT}_{BOTTOM_INPUT}.mp4"
OUTPUT_PDF = f"VF_plic_compare_{TOP_INPUT}_{BOTTOM_INPUT}_t0.pdf"

# Paper-style figure geometry (inches).
CM_TO_IN = 1.0 / 2.54
FIG_WIDTH_IN     = 9.5 * CM_TO_IN
LEFT_MARGIN_IN   = 0.78
RIGHT_MARGIN_IN  = 0.75
BOTTOM_MARGIN_IN = 0.62
TOP_MARGIN_IN    = 0.40
ROW_GAP_IN       = 0.30
CBAR_GAP_IN      = 0.13
CBAR_WIDTH_IN    = 0.17
TICK_PAD_PT      = 7.0

AXES_LABEL_FONTSIZE = 16.0
NUMBER_FONTSIZE     = 13.0
CBAR_LABEL_FONTSIZE = 14.0
TITLE_FONTSIZE      = 13.0

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
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    return rf"${value:.1f}$"


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


def gather_frames(case_input):
    """Sorted list of (t, plt_path, vtp_path) for a case. Only reads plotfile
    headers (cheap) -- no field/slice data is touched here."""
    amrviz_dir = Path("amrviz") / f"impact_relax_{case_input}"
    plt_files = sorted(amrviz_dir.glob("plt.nga2.cell.*"), key=frame_number)
    vtp_files = {frame_number(p): p for p in amrviz_dir.glob("plic_*.vtp")}
    pairs = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

    frames = []
    for plt_path, vtp_path in pairs:
        ds = yt.load(str(plt_path))
        t = float(ds.current_time.to_value())
        frames.append((t, plt_path, vtp_path))
    frames.sort(key=lambda f: f[0])
    return frames


@lru_cache(maxsize=None)
def load_frame(plt_path_str, vtp_path_str):
    plt_path = Path(plt_path_str)
    vtp_path = Path(vtp_path_str)

    ds = yt.load(str(plt_path))
    t = float(ds.current_time.to_value())

    max_level = ds.index.max_level
    res = int(ds.domain_dimensions[0]) * 2**max_level

    le = ds.domain_left_edge.to_value("code_length")
    re = ds.domain_right_edge.to_value("code_length")

    slc = ds.slice("z", 0.0)
    frb = slc.to_frb((re[0] - le[0], "code_length"), res, height=(re[1] - le[1], "code_length"))
    data = np.array(frb["boxlib", FIELD])

    segments = extract_plic_segments(vtp_path)
    return t, data, segments, le, re


def nearest_index(times, t):
    return int(np.argmin(np.abs(times - t)))


frames_top = gather_frames(TOP_INPUT)
frames_bottom = gather_frames(BOTTOM_INPUT)

times_top = np.array([f[0] for f in frames_top])
times_bottom = np.array([f[0] for f in frames_bottom])

# Drive the animation off whichever case has the finer temporal sampling;
# the other case's nearest-in-time frame is used at each step, and each
# row's title always reports that row's own actual simulation time.
if len(frames_top) >= len(frames_bottom):
    master_times = times_top
    master_is_top = True
else:
    master_times = times_bottom
    master_is_top = False


def paired_frame(i):
    t_master = master_times[i]
    if master_is_top:
        top = frames_top[i]
        bottom = frames_bottom[nearest_index(times_bottom, t_master)]
    else:
        bottom = frames_bottom[i]
        top = frames_top[nearest_index(times_top, t_master)]
    return top, bottom


x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
data_height_over_width = y_span / x_span

tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
x_ticks = tick_locator.tick_values(*VIEW_XLIM)
x_ticks = [v for v in x_ticks if VIEW_XLIM[0] <= v <= VIEW_XLIM[1]]
y_ticks = tick_locator.tick_values(*VIEW_YLIM)
y_ticks = [v for v in y_ticks if VIEW_YLIM[0] <= v <= VIEW_YLIM[1]]

PLOT_WIDTH_IN = FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + 2.0 * PANEL_HEIGHT_IN + ROW_GAP_IN + TOP_MARGIN_IN
CONTENT_WIDTH_IN = PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
LEFT_IN = 0.5 * (FIG_WIDTH_IN - CONTENT_WIDTH_IN)

BOTTOM_ROW_Y_IN = BOTTOM_MARGIN_IN
TOP_ROW_Y_IN = BOTTOM_MARGIN_IN + PANEL_HEIGHT_IN + ROW_GAP_IN

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
ax_top = add_axes_in_inches(fig, LEFT_IN, TOP_ROW_Y_IN, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
ax_bottom = add_axes_in_inches(fig, LEFT_IN, BOTTOM_ROW_Y_IN, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
cax = add_axes_in_inches(
    fig, LEFT_IN + PLOT_WIDTH_IN + CBAR_GAP_IN, BOTTOM_MARGIN_IN,
    CBAR_WIDTH_IN, 2.0 * PANEL_HEIGHT_IN + ROW_GAP_IN,
)


def setup_row(ax, t0, data0, segs0, le, re, label, show_xlabel):
    extent = [le[0], re[0], le[1], re[1]]
    im = ax.imshow(data0, origin="lower", extent=extent, cmap=CMAP, vmin=VMIN, vmax=VMAX, interpolation="none")
    lc = LineCollection(segs0, colors=INTERFACE_COLOR, linewidths=INTERFACE_LINEWIDTH)
    ax.add_collection(lc)

    ax.set_xlim(*VIEW_XLIM)
    ax.set_ylim(*VIEW_YLIM)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.set_yticks(y_ticks)
    ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.tick_params(which="both", top=True, right=True, pad=TICK_PAD_PT)

    if show_xlabel:
        ax.set_xlabel(r"$x$", labelpad=5.0)
    else:
        ax.set_xticklabels([])
    ax.set_ylabel(r"$y$", labelpad=-5.0)

    title = ax.set_title(
        rf"$\mathrm{{{label}}}: t = {f'{t0:.2f}'.rstrip('0').rstrip('.')}$",
        fontsize=TITLE_FONTSIZE, pad=6.0,
    )
    return im, lc, title


(top0, bottom0) = paired_frame(0)
t0_top, data0_top, segs0_top, le_top, re_top = load_frame(str(top0[1]), str(top0[2]))
t0_bottom, data0_bottom, segs0_bottom, le_bottom, re_bottom = load_frame(str(bottom0[1]), str(bottom0[2]))

im_top, lc_top, title_top = setup_row(ax_top, t0_top, data0_top, segs0_top, le_top, re_top, TOP_LABEL, show_xlabel=False)
im_bottom, lc_bottom, title_bottom = setup_row(ax_bottom, t0_bottom, data0_bottom, segs0_bottom, le_bottom, re_bottom, BOTTOM_LABEL, show_xlabel=True)

cbar = fig.colorbar(im_top, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)

if args.pdf:
    fig.savefig(OUTPUT_PDF)
    print(f"Saved: {OUTPUT_PDF}")
    raise SystemExit(0)


def update(i):
    top, bottom = paired_frame(i)
    t_top, data_top, segs_top, _, _ = load_frame(str(top[1]), str(top[2]))
    t_bottom, data_bottom, segs_bottom, _, _ = load_frame(str(bottom[1]), str(bottom[2]))

    im_top.set_data(data_top)
    lc_top.set_segments(segs_top)
    title_top.set_text(rf"$\mathrm{{{TOP_LABEL}}}: t = {f'{t_top:.2f}'.rstrip('0').rstrip('.')}$")

    im_bottom.set_data(data_bottom)
    lc_bottom.set_segments(segs_bottom)
    title_bottom.set_text(rf"$\mathrm{{{BOTTOM_LABEL}}}: t = {f'{t_bottom:.2f}'.rstrip('0').rstrip('.')}$")

    return im_top, lc_top, title_top, im_bottom, lc_bottom, title_bottom


anim = FuncAnimation(fig, update, frames=len(master_times), blit=False)
anim.save(OUTPUT, fps=FPS, dpi=200, writer="ffmpeg")
print(f"Saved: {OUTPUT}  ({len(master_times)} frames)")
