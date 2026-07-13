import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.collections import LineCollection
import pyvista as pv
import yt

parser = argparse.ArgumentParser()
parser.add_argument("input", help="case identifier, e.g. wall_NASG_pTg")
args = parser.parse_args()

INPUT = args.input
CASE = f"cavitation_{INPUT}"
PLTFILE = f"amrviz/{CASE}/plt.nga2.cell.000001"

CBAR_LABEL = r"$\left(\mathrm{m\,s^{-1}}\right)$"
CMAP = "jet"

VMIN_UV, VMAX_UV = -2150.0, 2150.0
VMIN_UMAG, VMAX_UMAG = 0.0, 2150.0
UV_TICKS = [-2150, -1000, 0, 1000, 2150]
UMAG_TICKS = [0, 500, 1000, 1500, 2150]

# Paper-style figure geometry (inches), following plot_RHOG_SG_NASG_PTg.py.
FIG_WIDTH_IN     = 6
UV_LEFT_MARGIN_IN  = 0.6
UV_RIGHT_MARGIN_IN = 1.0
PL_LEFT_MARGIN_IN  = 0.75
PL_RIGHT_MARGIN_IN = 0.8
BOTTOM_MARGIN_IN = 0.55
TOP_MARGIN_IN    = 0.3
ROW_GAP_IN       = 0.8
PANEL_GAP_IN     = 0.15
CBAR_GAP_IN      = 0.15
CBAR_WIDTH_IN    = 0.20

AXES_LABEL_FONTSIZE = 17
NUMBER_FONTSIZE     = 15
CBAR_LABEL_FONTSIZE = 16
TITLE_FONTSIZE      = 18

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
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    text = f"{value:.2f}".rstrip("0").rstrip(".")
    return rf"${text}$"


def cbar_tick_formatter(value, _):
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    text = f"{value:.1f}".rstrip("0").rstrip(".")
    return rf"${text}$"


def add_axes_in_inches(fig, left, bottom, width, height):
    fig_w, fig_h = fig.get_size_inches()
    return fig.add_axes([left / fig_w, bottom / fig_h, width / fig_w, height / fig_h])


def style_axes(ax, x_ticks, y_ticks, show_ylabel):
    ax.set_xlim(le[0], re[0])
    ax.set_ylim(le[1], re[1])
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.set_yticks(y_ticks)
    ax.tick_params(which="both", top=True, right=True, pad=2.5)

    if show_ylabel:
        ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_ylabel(r"$y\;\left(\mathrm{cm}\right)$", labelpad=2.5)
    else:
        ax.tick_params(labelleft=False)

    ax.set_xlabel(r"$x\;\left(\mathrm{cm}\right)$", labelpad=2.5)


def style_cbar(cbar, ticks):
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
    cbar.ax.tick_params(pad=2.5)
    cbar.set_label(CBAR_LABEL, rotation=90, labelpad=4.0, fontsize=CBAR_LABEL_FONTSIZE)


ds = yt.load(PLTFILE)

# Raw plotfile coordinates are already meters (yt assumes cm; ignore its unit label).
le_m = ds.domain_left_edge.to_value("code_length")
re_m = ds.domain_right_edge.to_value("code_length")

# Sample at the finest AMR level so each pixel matches one leaf cell -- no interpolation.
max_level = ds.index.max_level
res = int(ds.domain_dimensions[0]) * 2**max_level

slc = ds.slice("z", 0.0)
frb = slc.to_frb((re_m[0] - le_m[0], "code_length"), res, height=(re_m[1] - le_m[1], "code_length"))

data = {field: np.array(frb["boxlib", field]) for field in ("U", "V", "Umag")}

# Everything plotted from here on is in cm -- gives clean integer tick values.
M_TO_CM = 100.0
le = le_m * M_TO_CM
re = re_m * M_TO_CM
extent = [le[0], re[0], le[1], re[1]]

x_span = re[0] - le[0]
y_span = re[1] - le[1]
data_height_over_width = y_span / x_span

x_ticks = [-4, -2, 0, 2, 4]
y_ticks = x_ticks

# Row 1: U and V side by side, sharing one colorbar.
ROW1_PLOT_WIDTH_IN = (FIG_WIDTH_IN - UV_LEFT_MARGIN_IN - UV_RIGHT_MARGIN_IN - PANEL_GAP_IN - CBAR_GAP_IN - CBAR_WIDTH_IN) / 2
ROW1_PANEL_HEIGHT_IN = ROW1_PLOT_WIDTH_IN * data_height_over_width

# Row 2: Umag alone with its own colorbar, same panel width as row 1, centered.
ROW2_PLOT_WIDTH_IN = ROW1_PLOT_WIDTH_IN
ROW2_PANEL_HEIGHT_IN = ROW2_PLOT_WIDTH_IN * data_height_over_width
ROW2_CONTENT_WIDTH_IN = ROW2_PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
ROW2_LEFT_IN = 0.5 * (FIG_WIDTH_IN - ROW2_CONTENT_WIDTH_IN)

FIG_HEIGHT_IN = (BOTTOM_MARGIN_IN + ROW2_PANEL_HEIGHT_IN + ROW_GAP_IN + ROW1_PANEL_HEIGHT_IN + TOP_MARGIN_IN)

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

row2_bottom = BOTTOM_MARGIN_IN
row1_bottom = row2_bottom + ROW2_PANEL_HEIGHT_IN + ROW_GAP_IN

ax_U = add_axes_in_inches(fig, UV_LEFT_MARGIN_IN, row1_bottom, ROW1_PLOT_WIDTH_IN, ROW1_PANEL_HEIGHT_IN)
ax_V = add_axes_in_inches(
    fig,
    UV_LEFT_MARGIN_IN + ROW1_PLOT_WIDTH_IN + PANEL_GAP_IN,
    row1_bottom,
    ROW1_PLOT_WIDTH_IN,
    ROW1_PANEL_HEIGHT_IN,
)
cax_UV = add_axes_in_inches(
    fig,
    UV_LEFT_MARGIN_IN + 2 * ROW1_PLOT_WIDTH_IN + PANEL_GAP_IN + CBAR_GAP_IN,
    row1_bottom,
    CBAR_WIDTH_IN,
    ROW1_PANEL_HEIGHT_IN,
)

ax_Umag = add_axes_in_inches(fig, ROW2_LEFT_IN, row2_bottom, ROW2_PLOT_WIDTH_IN, ROW2_PANEL_HEIGHT_IN)
cax_Umag = add_axes_in_inches(
    fig,
    ROW2_LEFT_IN + ROW2_PLOT_WIDTH_IN + CBAR_GAP_IN,
    row2_bottom,
    CBAR_WIDTH_IN,
    ROW2_PANEL_HEIGHT_IN,
)

im_U = ax_U.imshow(data["U"], origin="lower", extent=extent, cmap=CMAP, vmin=VMIN_UV, vmax=VMAX_UV, interpolation="none")
style_axes(ax_U, x_ticks, y_ticks, show_ylabel=True)
ax_U.set_title(r"$U$", fontsize=TITLE_FONTSIZE, pad=5.0)

im_V = ax_V.imshow(data["V"], origin="lower", extent=extent, cmap=CMAP, vmin=VMIN_UV, vmax=VMAX_UV, interpolation="none")
style_axes(ax_V, x_ticks, y_ticks, show_ylabel=False)
ax_V.set_title(r"$V$", fontsize=TITLE_FONTSIZE, pad=5.0)

cbar_UV = fig.colorbar(im_V, cax=cax_UV, orientation="vertical")
style_cbar(cbar_UV, UV_TICKS)

im_Umag = ax_Umag.imshow(data["Umag"], origin="lower", extent=extent, cmap=CMAP, vmin=VMIN_UMAG, vmax=VMAX_UMAG, interpolation="none")
style_axes(ax_Umag, x_ticks, y_ticks, show_ylabel=True)
ax_Umag.set_title(r"$|\vec{U}|$", fontsize=TITLE_FONTSIZE, pad=5.0)

cbar_Umag = fig.colorbar(im_Umag, cax=cax_Umag, orientation="vertical")
style_cbar(cbar_Umag, UMAG_TICKS)

out = f"UV_Umag0_{INPUT}.pdf"
fig.savefig(out)
print(f"Saved: {out}")


# =============================================================================
# Second figure: 3x3 grid of liquid pressure p_l at 9 different times.
# Same styling/fonts/gaps/FIG_WIDTH_IN as the velocity figure above; panel
# dimensions fall out of the same margin arithmetic, just for a 3x3 grid.
# =============================================================================

if INPUT.startswith("wall_"):
    PL_FRAMES = [1, 12, 23, 34, 45, 56, 66, 76, 86]
elif INPUT.startswith("dirichlet_"):
    PL_FRAMES = [1, 3, 6, 9, 11, 14, 16, 18, 21]
else:
    raise ValueError(f"no PL_FRAMES defined for input {INPUT!r}")
PL_ROW_GAP_IN = PANEL_GAP_IN   # vertical gap between grid rows == horizontal gap between columns

# PL swings from ~1e-4 GPa (ambient, t=0) to O(1) GPa (post-impact) down to
# O(0.1) GPa (cavitation rarefaction spikes at later times); shown here on a
# single shared linear color scale (vmin/vmax computed from all frames below)
# -- most panels will read as close to one flat color since the early/late
# frames differ by orders of magnitude, but every panel uses the exact same
# scale for direct comparison.
PL_SCALE = 1.0e9  # Pa -> GPa
PL_CBAR_LABEL = r"$p_l\;\left(\mathrm{GPa}\right)$"

N_ROWS, N_COLS = 3, 3

PL_VIEW_XLIM = (-5, 5)
PL_VIEW_YLIM = (-5, 5)
PL_X_TICKS = [-3, 0, 3]
PL_Y_TICKS = PL_X_TICKS


def style_axes_grid(ax, x_ticks, y_ticks, xlim, ylim, show_xlabel, show_ylabel):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.set_yticks(y_ticks)
    ax.tick_params(which="both", top=True, right=True, pad=2.5)

    if show_ylabel:
        ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_ylabel(r"$y\;\left(\mathrm{cm}\right)$", labelpad=2.5)
    else:
        ax.tick_params(labelleft=False)

    if show_xlabel:
        ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
        ax.set_xlabel(r"$x\;\left(\mathrm{cm}\right)$", labelpad=2.5)
    else:
        ax.tick_params(labelbottom=False)


def load_pl(frame_number):
    path = f"amrviz/{CASE}/plt.nga2.cell.{frame_number:06d}"
    ds_i = yt.load(path)
    t = float(ds_i.current_time.to_value())
    slc_i = ds_i.slice("z", 0.0)
    frb_i = slc_i.to_frb((re_m[0] - le_m[0], "code_length"), res, height=(re_m[1] - le_m[1], "code_length"))
    pl = np.array(frb_i["boxlib", "PL"]) / PL_SCALE
    return t, pl


PLIC_COLOR = "white"
PLIC_LINEWIDTH = 0.4


def extract_plic_segments(vtp_path):
    """Reduce the extruded PLIC quads (2D interface segments extruded in z)
    back down to 2D line segments in the x-y plane, in cm."""
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


def load_plic(frame_number):
    path = f"amrviz/{CASE}/plic_{frame_number:06d}.vtp"
    return extract_plic_segments(path)


PL_PLOT_WIDTH_IN = (
    FIG_WIDTH_IN - PL_LEFT_MARGIN_IN - PL_RIGHT_MARGIN_IN - (N_COLS - 1) * PANEL_GAP_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
) / N_COLS
PL_PANEL_HEIGHT_IN = PL_PLOT_WIDTH_IN * data_height_over_width
PL_GRID_HEIGHT_IN = N_ROWS * PL_PANEL_HEIGHT_IN + (N_ROWS - 1) * PL_ROW_GAP_IN
PL_FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + PL_GRID_HEIGHT_IN + TOP_MARGIN_IN

fig2 = plt.figure(figsize=(FIG_WIDTH_IN, PL_FIG_HEIGHT_IN))

cax_PL = add_axes_in_inches(
    fig2,
    PL_LEFT_MARGIN_IN + N_COLS * PL_PLOT_WIDTH_IN + (N_COLS - 1) * PANEL_GAP_IN + CBAR_GAP_IN,
    BOTTOM_MARGIN_IN,
    CBAR_WIDTH_IN,
    PL_GRID_HEIGHT_IN,
)

# Load every frame up front so all panels can share one fixed color scale --
# per-panel autoscaling would blow up round-off-level noise (e.g. the near-
# uniform t=0 frame) into a full-colormap checkerboard tracing AMR patch edges.
pl_frames = [load_pl(frame_number) for frame_number in PL_FRAMES]
PL_VMIN = min(pl.min() for _, pl in pl_frames)
PL_VMAX = max(pl.max() for _, pl in pl_frames)

im_PL = None
for panel_idx, frame_number in enumerate(PL_FRAMES):
    visual_row = panel_idx // N_COLS       # 0 = top row, N_ROWS-1 = bottom row
    visual_col = panel_idx % N_COLS        # 0 = left column
    placement_row_from_bottom = (N_ROWS - 1) - visual_row

    left_in = PL_LEFT_MARGIN_IN + visual_col * (PL_PLOT_WIDTH_IN + PANEL_GAP_IN)
    bottom_in = BOTTOM_MARGIN_IN + placement_row_from_bottom * (PL_PANEL_HEIGHT_IN + PL_ROW_GAP_IN)

    ax = add_axes_in_inches(fig2, left_in, bottom_in, PL_PLOT_WIDTH_IN, PL_PANEL_HEIGHT_IN)

    t, pl = pl_frames[panel_idx]
    im_PL = ax.imshow(pl, origin="lower", extent=extent, cmap=CMAP, vmin=PL_VMIN, vmax=PL_VMAX, interpolation="none")

    plic_segments = load_plic(frame_number)
    ax.add_collection(LineCollection(plic_segments, colors=PLIC_COLOR, linewidths=PLIC_LINEWIDTH))

    style_axes_grid(
        ax, PL_X_TICKS, PL_Y_TICKS, PL_VIEW_XLIM, PL_VIEW_YLIM,
        show_xlabel=(visual_row == N_ROWS - 1),
        show_ylabel=(visual_col == 0),
    )

    ax.text(
        0.06, 0.94,
        rf"${f'{t * 1.0e6:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$",
        transform=ax.transAxes, ha="left", va="top", fontsize=NUMBER_FONTSIZE,
        color="black",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.75),
    )

cbar_PL = fig2.colorbar(im_PL, cax=cax_PL, orientation="vertical")
cbar_PL.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar_PL.ax.tick_params(pad=2.5)
cbar_PL.set_label(PL_CBAR_LABEL, rotation=90, labelpad=4.0, fontsize=CBAR_LABEL_FONTSIZE)

out2 = f"PL_{INPUT}.pdf"
fig2.savefig(out2)
print(f"Saved: {out2}")


# =============================================================================
# Third figure: 3x3 grid of vapor volume fraction VF (with PLIC interface) at
# the same 9 frames as the pressure grid above, but zoomed in on the bubble.
# Same layout/margins/panel geometry as the PL grid, just a narrower view.
# =============================================================================

VF_CBAR_LABEL = r"$\alpha$"
VF_VMIN, VF_VMAX = 0.0, 1.0
VF_TICKS = np.linspace(VF_VMIN, VF_VMAX, 5)

VF_VIEW_XLIM = (-1.5, 1.5)
VF_VIEW_YLIM = (-1.5, 1.5)
VF_X_TICKS = [-1, 0, 1]
VF_Y_TICKS = VF_X_TICKS


def load_vf(frame_number):
    path = f"amrviz/{CASE}/plt.nga2.cell.{frame_number:06d}"
    ds_i = yt.load(path)
    t = float(ds_i.current_time.to_value())
    slc_i = ds_i.slice("z", 0.0)
    frb_i = slc_i.to_frb((re_m[0] - le_m[0], "code_length"), res, height=(re_m[1] - le_m[1], "code_length"))
    vf = np.array(frb_i["boxlib", "VF"])
    return t, vf


fig3 = plt.figure(figsize=(FIG_WIDTH_IN, PL_FIG_HEIGHT_IN))

cax_VF = add_axes_in_inches(
    fig3,
    PL_LEFT_MARGIN_IN + N_COLS * PL_PLOT_WIDTH_IN + (N_COLS - 1) * PANEL_GAP_IN + CBAR_GAP_IN,
    BOTTOM_MARGIN_IN,
    CBAR_WIDTH_IN,
    PL_GRID_HEIGHT_IN,
)

im_VF = None
for panel_idx, frame_number in enumerate(PL_FRAMES):
    visual_row = panel_idx // N_COLS       # 0 = top row, N_ROWS-1 = bottom row
    visual_col = panel_idx % N_COLS        # 0 = left column
    placement_row_from_bottom = (N_ROWS - 1) - visual_row

    left_in = PL_LEFT_MARGIN_IN + visual_col * (PL_PLOT_WIDTH_IN + PANEL_GAP_IN)
    bottom_in = BOTTOM_MARGIN_IN + placement_row_from_bottom * (PL_PANEL_HEIGHT_IN + PL_ROW_GAP_IN)

    ax = add_axes_in_inches(fig3, left_in, bottom_in, PL_PLOT_WIDTH_IN, PL_PANEL_HEIGHT_IN)

    t, vf = load_vf(frame_number)
    im_VF = ax.imshow(vf, origin="lower", extent=extent, cmap=CMAP, vmin=VF_VMIN, vmax=VF_VMAX, interpolation="none")

    plic_segments = load_plic(frame_number)
    ax.add_collection(LineCollection(plic_segments, colors=PLIC_COLOR, linewidths=PLIC_LINEWIDTH))

    style_axes_grid(
        ax, VF_X_TICKS, VF_Y_TICKS, VF_VIEW_XLIM, VF_VIEW_YLIM,
        show_xlabel=(visual_row == N_ROWS - 1),
        show_ylabel=(visual_col == 0),
    )

    ax.text(
        0.06, 0.94,
        rf"${f'{t * 1.0e6:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$",
        transform=ax.transAxes, ha="left", va="top", fontsize=NUMBER_FONTSIZE,
        color="black",
        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.75),
    )

cbar_VF = fig3.colorbar(im_VF, cax=cax_VF, orientation="vertical")
cbar_VF.set_ticks(VF_TICKS)
cbar_VF.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar_VF.ax.tick_params(pad=2.5)
cbar_VF.set_label(VF_CBAR_LABEL, rotation=90, labelpad=4.0, fontsize=CBAR_LABEL_FONTSIZE)

out3 = f"VF_plic_{INPUT}.pdf"
fig3.savefig(out3)
print(f"Saved: {out3}")
