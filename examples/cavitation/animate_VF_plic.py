import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter
from matplotlib.animation import FuncAnimation
import pyvista as pv
import yt

parser = argparse.ArgumentParser()
parser.add_argument("input", help="case identifier, e.g. NASG_relax_pTg")
args = parser.parse_args()

INPUT = args.input
CASE = f"cavitation_{INPUT}"
AMRVIZ_DIR = Path("amrviz") / CASE

FIELD = "VF"
CBAR_LABEL = r"$\alpha$"
VMIN, VMAX = 0.0, 1.0

# View limits and all plotted coordinates are in cm (raw plotfile/vtp data is in meters).
M_TO_CM = 100.0
VIEW_XLIM = (-1.5, 1.5)
VIEW_YLIM = (-1.5, 1.5)

CMAP = "jet"
INTERFACE_COLOR = "white"
INTERFACE_LINEWIDTH = 1.2

FPS = 12
OUTPUT = f"VF_plic_{INPUT}.mp4"

# Paper-style figure geometry (inches). Canvas sized to ~9.5 x 7.0 cm --
# smaller physical size, so fonts/lines are scaled up to stay legible and
# margins are generous enough that the rotated y-label always fits.
CM_TO_IN = 1.0 / 2.54
FIG_WIDTH_IN     = 9.5 * CM_TO_IN
LEFT_MARGIN_IN   = 0.78
RIGHT_MARGIN_IN  = 0.75
BOTTOM_MARGIN_IN = 0.62
TOP_MARGIN_IN    = 0.40
CBAR_GAP_IN      = 0.13
CBAR_WIDTH_IN    = 0.17
TICK_PAD_PT      = 7.0

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
    return segments * M_TO_CM


def load_frame(plt_path: Path, vtp_path: Path):
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


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

t0, data0, segs0, le_m, re_m = load_frame(*frames[0])
le = le_m * M_TO_CM
re = re_m * M_TO_CM

x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
data_height_over_width = y_span / x_span
extent = [le[0], re[0], le[1], re[1]]
x_ticks = [-1, 0.0, 1]
y_ticks = [-1, 0.0, 1]

PLOT_WIDTH_IN = FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + PANEL_HEIGHT_IN + TOP_MARGIN_IN
CONTENT_WIDTH_IN = PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
LEFT_IN = 0.5 * (FIG_WIDTH_IN - CONTENT_WIDTH_IN)

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
ax = add_axes_in_inches(fig, LEFT_IN, BOTTOM_MARGIN_IN, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
cax = add_axes_in_inches(
    fig, LEFT_IN + PLOT_WIDTH_IN + CBAR_GAP_IN, BOTTOM_MARGIN_IN, CBAR_WIDTH_IN, PANEL_HEIGHT_IN
)

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

ax.set_xlabel(r"$x\;\left(\mathrm{cm}\right)$", labelpad=5.0)
ax.set_ylabel(r"$y\;\left(\mathrm{cm}\right)$", labelpad=-5.0)
title = ax.set_title(
    rf"$t = {f'{t0 * 1.0e6:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$", fontsize=TITLE_FONTSIZE, pad=6.0
)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)


def update(i):
    plt_path, vtp_path = frames[i]
    t, data, segments, _, _ = load_frame(plt_path, vtp_path)
    im.set_data(data)
    lc.set_segments(segments)
    title.set_text(rf"$t = {f'{t * 1.0e6:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$")
    return im, lc, title


anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
anim.save(OUTPUT, fps=FPS, dpi=200, writer="ffmpeg")
print(f"Saved: {OUTPUT}  ({len(frames)} frames)")
