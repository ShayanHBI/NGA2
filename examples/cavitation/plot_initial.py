import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import yt

CASE = "cavitation_NASG_relax_pTg"
PLTFILE = f"amrviz/{CASE}/plt.nga2.cell.000001"

CBAR_LABEL = r"$\left(\mathrm{m\,s^{-1}}\right)$"
CMAP = "jet"

VMIN_UV, VMAX_UV = -50.0, 50.0
VMIN_UMAG, VMAX_UMAG = 35.0, 50.0
UV_TICKS = [-50, -25, 0, 25, 50]
UMAG_TICKS = [35, 40, 45, 50]

# Paper-style figure geometry (inches), following plot_RHOG_SG_NASG_PTg.py.
FIG_WIDTH_IN     = 7.2
LEFT_MARGIN_IN   = 0.68
RIGHT_MARGIN_IN  = 0.80
BOTTOM_MARGIN_IN = 0.50
TOP_MARGIN_IN    = 0.38
ROW_GAP_IN       = 0.85
PANEL_GAP_IN     = 0.40
CBAR_GAP_IN      = 0.15
CBAR_WIDTH_IN    = 0.20

AXES_LABEL_FONTSIZE = 13.0
NUMBER_FONTSIZE     = 11.0
CBAR_LABEL_FONTSIZE = 12.0
TITLE_FONTSIZE      = 14.0

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
        ax.set_ylabel(r"$y\;\left(\mathrm{m}\right)$", labelpad=2.5)
    else:
        ax.tick_params(labelleft=False)

    ax.set_xlabel(r"$x\;\left(\mathrm{m}\right)$", labelpad=2.5)


def style_cbar(cbar, ticks):
    cbar.set_ticks(ticks)
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
    cbar.ax.tick_params(pad=2.5)
    cbar.set_label(CBAR_LABEL, rotation=90, labelpad=4.0, fontsize=CBAR_LABEL_FONTSIZE)


ds = yt.load(PLTFILE)

# Raw plotfile coordinates are already meters (yt assumes cm; ignore its unit label).
le = ds.domain_left_edge.to_value("code_length")
re = ds.domain_right_edge.to_value("code_length")
extent = [le[0], re[0], le[1], re[1]]

# Sample at the finest AMR level so each pixel matches one leaf cell -- no interpolation.
max_level = ds.index.max_level
res = int(ds.domain_dimensions[0]) * 2**max_level

slc = ds.slice("z", 0.0)
frb = slc.to_frb((re[0] - le[0], "code_length"), res, height=(re[1] - le[1], "code_length"))

data = {field: np.array(frb["boxlib", field]) for field in ("U", "V", "Umag")}

x_span = re[0] - le[0]
y_span = re[1] - le[1]
data_height_over_width = y_span / x_span

x_ticks = [le[0], 0.0, re[0]]
y_ticks = [le[1], 0.0, re[1]]

# Row 1: U and V side by side, sharing one colorbar.
ROW1_PLOT_WIDTH_IN = (
    FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - PANEL_GAP_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
) / 2
ROW1_PANEL_HEIGHT_IN = ROW1_PLOT_WIDTH_IN * data_height_over_width

# Row 2: Umag alone with its own colorbar, same panel width as row 1, centered.
ROW2_PLOT_WIDTH_IN = ROW1_PLOT_WIDTH_IN
ROW2_PANEL_HEIGHT_IN = ROW2_PLOT_WIDTH_IN * data_height_over_width
ROW2_CONTENT_WIDTH_IN = ROW2_PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
ROW2_LEFT_IN = 0.5 * (FIG_WIDTH_IN - ROW2_CONTENT_WIDTH_IN)

FIG_HEIGHT_IN = (
    BOTTOM_MARGIN_IN + ROW2_PANEL_HEIGHT_IN + ROW_GAP_IN + ROW1_PANEL_HEIGHT_IN + TOP_MARGIN_IN
)

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

out = "UV_Umag0.pdf"
fig.savefig(out)
print(f"Saved: {out}")
