import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import yt

CASE = "impact_relax_pTg"
PLTFILE = f"amrviz/{CASE}/plt.nga2.cell.000001"
FIELD = "Mach"
CBAR_LABEL = r"$\mathrm{Mach}$"

CMAP = "jet"

# Paper-style figure geometry (inches) -- identical to plot_initial.py in
# examples/cavitation, so panel size, fonts, and line weights match exactly.
FIG_WIDTH_IN     = 7.2
LEFT_MARGIN_IN   = 0.68
RIGHT_MARGIN_IN  = 0.80
BOTTOM_MARGIN_IN = 0.50
TOP_MARGIN_IN    = 0.38
PANEL_GAP_IN     = 0.40
CBAR_GAP_IN      = 0.15
CBAR_WIDTH_IN    = 0.20

AXES_LABEL_FONTSIZE = 13.0
NUMBER_FONTSIZE     = 11.0
CBAR_LABEL_FONTSIZE = 12.0

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
data = np.array(frb["boxlib", FIELD])

vmin = float(data.min())
vmax = float(data.max())

x_span = re[0] - le[0]
y_span = re[1] - le[1]
data_height_over_width = y_span / x_span

x_ticks = [le[0], 0.5 * (le[0] + re[0]), re[0]]
y_ticks = [le[1], 0.5 * (le[1] + re[1]), re[1]]

# Same panel width as the two-panel row in examples/cavitation/plot_initial.py.
PLOT_WIDTH_IN = (
    FIG_WIDTH_IN - LEFT_MARGIN_IN - RIGHT_MARGIN_IN - PANEL_GAP_IN - CBAR_GAP_IN - CBAR_WIDTH_IN
) / 2
PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + PANEL_HEIGHT_IN + TOP_MARGIN_IN

# Single panel, centered on the canvas.
CONTENT_WIDTH_IN = PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
LEFT_IN = 0.5 * (FIG_WIDTH_IN - CONTENT_WIDTH_IN)

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
ax = add_axes_in_inches(fig, LEFT_IN, BOTTOM_MARGIN_IN, PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
cax = add_axes_in_inches(
    fig,
    LEFT_IN + PLOT_WIDTH_IN + CBAR_GAP_IN,
    BOTTOM_MARGIN_IN,
    CBAR_WIDTH_IN,
    PANEL_HEIGHT_IN,
)

im = ax.imshow(data, origin="lower", extent=extent, cmap=CMAP, vmin=vmin, vmax=vmax, interpolation="none")

ax.set_xlim(le[0], re[0])
ax.set_ylim(le[1], re[1])
ax.set_aspect("equal", adjustable="box")

ax.set_xticks(x_ticks)
ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
ax.set_yticks(y_ticks)
ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
ax.tick_params(which="both", top=True, right=True, pad=2.5)

ax.set_xlabel(r"$x$", labelpad=2.5)
ax.set_ylabel(r"$y$", labelpad=2.5)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(vmin, vmax, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=2.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=4.0, fontsize=CBAR_LABEL_FONTSIZE)

out = "Mach0.pdf"
fig.savefig(out)
print(f"Saved: {out}")
