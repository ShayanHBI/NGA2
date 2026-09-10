"""
Option 2 from the "how do I avoid hand-picking margins" discussion:
instead of declaring LEFT/RIGHT/BOTTOM/TOP_MARGIN_IN as constants, this
script *measures* how much space the actual rendered text (fonts, tick
labels, axis labels, title, colorbar label) needs, and derives the
margins from that measurement plus one explicit EDGE_PADDING_IN.

Mechanism: matplotlib can only report a text artist's true size after a
render pass. So this does a throwaway "staging" figure/axes styled
exactly like the real panels, forces a draw, and compares
ax.get_tightbbox() (full extent including labels/ticks/title) against
ax.get_window_extent() (just the data box) to get the overflow on each
side. That overflow, plus padding, becomes the real margin.

Trade-off vs. plot_initial.py: margins here are correct-by-construction
for whatever tick values/labels are actually used, but are NOT
guaranteed identical to another figure with different tick content --
see the "hybrid" option (measure once, freeze the result) if you need
that guarantee back.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import yt

parser = argparse.ArgumentParser()
parser.add_argument("input", help="case identifier, e.g. wall_pTg")
args = parser.parse_args()

INPUT = args.input
CASE = f"cavitation_{INPUT}"
PLTFILE = f"amrviz/{CASE}/plt.nga2.cell.000001"

CBAR_LABEL = r"$\left(\mathrm{m\,s^{-1}}\right)$"
CMAP = "jet"

# ---- What you still choose explicitly ("the rest of things") ----------
FIG_WIDTH_IN     = 6
ROW_GAP_IN       = 0.8
PANEL_GAP_IN     = 0.15
CBAR_GAP_IN      = 0.15
CBAR_WIDTH_IN    = 0.20

# The one new free parameter: distance between the outermost rendered
# text and the figure's edge. This replaces LEFT/RIGHT/BOTTOM/TOP_MARGIN_IN.
EDGE_PADDING_IN = 0.06

AXES_LABEL_FONTSIZE = 17#13.0
NUMBER_FONTSIZE     = 15#11.0
CBAR_LABEL_FONTSIZE = 16#12.0
TITLE_FONTSIZE      = 18#14.0

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
    return rf"${value:.2f}$"


def cbar_tick_formatter(value, _):
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    return rf"${value:.1f}$"


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
        ax.set_ylabel(r"$y\;\left(\mathrm{mm}\right)$", labelpad=2.5)
    else:
        ax.tick_params(labelleft=False)

    ax.set_xlabel(r"$x\;\left(\mathrm{mm}\right)$", labelpad=2.5)


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

# Colorbar ranges/ticks follow the actual data range of this run.
VMIN_UV = min(data["U"].min(), data["V"].min())
VMAX_UV = max(data["U"].max(), data["V"].max())
UV_TICKS = np.linspace(VMIN_UV, VMAX_UV, 5)

VMIN_UMAG = data["Umag"].min()
VMAX_UMAG = data["Umag"].max()
UMAG_TICKS = np.linspace(VMIN_UMAG, VMAX_UMAG, 5)

# Everything plotted from here on is in mm -- gives clean integer tick values.
M_TO_MM = 1000.0
le = le_m * M_TO_MM
re = re_m * M_TO_MM
extent = [le[0], re[0], le[1], re[1]]

x_span = re[0] - le[0]
y_span = re[1] - le[1]
data_height_over_width = y_span / x_span

x_ticks = [-10, -5, 0, 5, 10]
y_ticks = x_ticks


# =============================================================================
# Measure the margins instead of hand-picking them.
# =============================================================================

def _overflow_in(artist, renderer, fig):
    """How far an artist's rendered extent (ticks+labels+title) spills
    past its own data box, per side, in inches."""
    full = artist.get_tightbbox(renderer)
    data_box = artist.get_window_extent(renderer)
    dpi = fig.dpi
    return {
        "left":   (data_box.x0 - full.x0) / dpi,
        "right":  (full.x1 - data_box.x1) / dpi,
        "bottom": (data_box.y0 - full.y0) / dpi,
        "top":    (full.y1 - data_box.y1) / dpi,
    }


def measure_required_margins():
    # Throwaway figure: any size works, since text size in inches doesn't
    # depend on how big the axes box happens to be -- only on font size
    # and the actual tick/label/title strings, which we set identically
    # to the real panels.
    staging_fig = plt.figure(figsize=(4, 4))

    ax = staging_fig.add_axes([0.3, 0.3, 0.4, 0.4])
    ax.imshow(np.zeros((2, 2)), origin="lower", extent=extent, cmap=CMAP, vmin=VMIN_UV, vmax=VMAX_UV)
    style_axes(ax, x_ticks, y_ticks, show_ylabel=True)   # worst case: leftmost panel
    ax.set_title(r"$U$", fontsize=TITLE_FONTSIZE, pad=5.0)

    cax = staging_fig.add_axes([0.8, 0.3, 0.05, 0.4])
    sm = plt.cm.ScalarMappable(cmap=CMAP)
    sm.set_array([])
    sm.set_clim(VMIN_UV, VMAX_UV)
    cbar = staging_fig.colorbar(sm, cax=cax, orientation="vertical")
    style_cbar(cbar, UV_TICKS)   # UV_TICKS has more/wider digits than UMAG_TICKS

    staging_fig.canvas.draw()
    renderer = staging_fig.canvas.get_renderer()

    panel = _overflow_in(ax, renderer, staging_fig)
    cbar_ov = _overflow_in(cax, renderer, staging_fig)

    plt.close(staging_fig)

    return {
        "left":   panel["left"] + EDGE_PADDING_IN,
        "bottom": panel["bottom"] + EDGE_PADDING_IN,
        "top":    panel["top"] + EDGE_PADDING_IN,
        "right":  cbar_ov["right"] + EDGE_PADDING_IN,
    }


margins = measure_required_margins()
LEFT_MARGIN_IN   = margins["left"]
RIGHT_MARGIN_IN  = margins["right"]
BOTTOM_MARGIN_IN = margins["bottom"]
TOP_MARGIN_IN    = margins["top"]

print(
    "Measured margins (in): "
    f"left={LEFT_MARGIN_IN:.3f} right={RIGHT_MARGIN_IN:.3f} "
    f"bottom={BOTTOM_MARGIN_IN:.3f} top={TOP_MARGIN_IN:.3f}"
)


# =============================================================================
# From here down this is identical to plot_initial.py -- only the margin
# numbers feeding in are now measured instead of hand-picked.
# =============================================================================

# Row 1: U and V side by side, sharing one colorbar.
ROW1_PLOT_WIDTH_IN = (FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - PANEL_GAP_IN - CBAR_GAP_IN - CBAR_WIDTH_IN) / 2
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

ax_U = add_axes_in_inches(fig, LEFT_MARGIN_IN, row1_bottom, ROW1_PLOT_WIDTH_IN, ROW1_PANEL_HEIGHT_IN)
ax_V = add_axes_in_inches(
    fig,
    LEFT_MARGIN_IN + ROW1_PLOT_WIDTH_IN + PANEL_GAP_IN,
    row1_bottom,
    ROW1_PLOT_WIDTH_IN,
    ROW1_PANEL_HEIGHT_IN,
)
cax_UV = add_axes_in_inches(
    fig,
    LEFT_MARGIN_IN + 2 * ROW1_PLOT_WIDTH_IN + PANEL_GAP_IN + CBAR_GAP_IN,
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

out = f"UV_Umag0_auto_layout_{INPUT}.pdf"
fig.savefig(out)
print(f"Saved: {out}")
