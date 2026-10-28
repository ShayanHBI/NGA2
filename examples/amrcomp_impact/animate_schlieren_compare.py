import re
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.animation import FuncAnimation
import pyvista as pv
import yt

# Cases to compare: top row / bottom row (case identifier -> amrviz/impact_relax_<INPUT>).
TOP_INPUT, TOP_LABEL = "pTg_inst", r"Instantaneous MTC"
BOTTOM_INPUT, BOTTOM_LABEL = "pTg_finite", r"Instantaneous MT, finite C"

CBAR_LABEL = r"$\phi$"
VMIN, VMAX = 0.0, 1.0

# Case is non-dimensional -- plotted coordinates and time are used as-is.
# VIEW_XLIM = (0.0, 11.0)
# VIEW_YLIM = (-5.5, 5.5)
VIEW_XLIM = (2, 7)
VIEW_YLIM = (-2, 2)

CMAP = "gray"
SHOW_INTERFACE = True
INTERFACE_COLOR = "red"
INTERFACE_LINEWIDTH = 0.2

# Set to a simulation time to render a single PDF of the nearest common
# frame instead of the full animation. Set to None for the full animation.
SNAPSHOT_TIME = 7    # e.g. 5.0

FPS = 12
ANIMATION_DPI = 200
OUTPUT = f"schlieren_plic_compare_{TOP_INPUT}_{BOTTOM_INPUT}.mp4"

# --- schlieren options (same definition as animate_field.py) ---------------
# phi = exp(-CONTRAST*[KL_REL*VF + (1-VF)] * |grad(rho)| / NORM), rho the
# mixture density, NORM a high percentile of |grad(rho)| over the gas region
# of each frame. Cells with any liquid are excluded and masked to phi=1.
KL = 120.0
KG = 30.0
KL_REL = KL / KG
CONTRAST = 10.0
GRAD_NORM_PERCENTILE = 99.9
VF_MASK_THRESHOLD = 1.0e-3
# Animation: uniform blur hides the coarse-AMR staircase (each coarse cell is replicated over 2^(Lmax-L) pixels)
SMOOTH_SIGMA_PIX = 5.0
# Snapshot: maximum resolution -- finest-level pixels unsmoothed, level-L pixels blurred by
# SNAPSHOT_SMOOTH_CELLS of their own cell size (just enough to remove the staircase)
SNAPSHOT_SMOOTH_CELLS = 0.5

# Frames of the two cases are paired by simulation time; only times both cases have are animated.
# Pairing tolerance as a fraction of the output interval (outputs land on the first step past each multiple).
TIME_MATCH_FRAC = 0.25

# Paper-style figure geometry (inches).
CM_TO_IN = 1.0 / 2.54
FIG_WIDTH_IN     = 9.5 * CM_TO_IN
LEFT_MARGIN_IN   = 0.78
RIGHT_MARGIN_IN  = 0.75
BOTTOM_MARGIN_IN = 0.62
TOP_MARGIN_IN    = 0.62
ROW_GAP_IN       = 0.30
CBAR_GAP_IN      = 0.13
CBAR_WIDTH_IN    = 0.17
TICK_PAD_PT      = 7.0

AXES_LABEL_FONTSIZE = 16.0
NUMBER_FONTSIZE     = 13.0
CBAR_LABEL_FONTSIZE = 14.0
TITLE_FONTSIZE      = 11.0
TIME_FONTSIZE       = 13.0

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


def header_time(plt_path: Path):
    """Simulation time from the plotfile Header (avoids a yt.load per frame)."""
    lines = (plt_path / "Header").read_text().split("\n")
    nvars = int(lines[1])
    return float(lines[2 + nvars + 1])


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
    """Sorted list of (t, plt_path, vtp_path) for a case."""
    amrviz_dir = Path("amrviz") / f"impact_relax_{case_input}"
    plt_files = sorted(amrviz_dir.glob("plt.nga2.cell.*"), key=frame_number)
    vtp_files = {frame_number(p): p for p in amrviz_dir.glob("plic_*.vtp")}
    frames = [(header_time(p), p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]
    frames.sort(key=lambda f: f[0])
    return frames


def gas_smooth(rho, gas, sigma):
    """Gaussian blur restricted to gas pixels (normalized convolution), so the
    liquid's much larger density never leaks into the average."""
    if sigma <= 0.0:
        # Keep gas pixels as they are; liquid pixels take the nearby gas average so the gradient never sees the liquid
        return np.where(gas, rho, gas_smooth(rho, gas, 1.0))
    weighted = gaussian_filter(np.where(gas, rho, 0.0), sigma=sigma)
    weight = gaussian_filter(gas.astype(rho.dtype), sigma=sigma)
    return weighted / np.maximum(weight, np.finfo(rho.dtype).tiny)


def compute_schlieren(VF, RHOL, RHOG, dx, dy, level=None, max_level=None):
    rho = VF * RHOL + (1.0 - VF) * RHOG
    gas = VF < VF_MASK_THRESHOLD

    if level is None:
        rho_smooth = gas_smooth(rho, gas, SMOOTH_SIGMA_PIX)
    else:
        rho_smooth = np.empty_like(rho)
        for lev in np.unique(level):
            sigma = SNAPSHOT_SMOOTH_CELLS * 2.0**(max_level - lev) if lev < max_level else 0.0
            sel = level == lev
            rho_smooth[sel] = gas_smooth(rho, gas, sigma)[sel]

    drho_dy, drho_dx = np.gradient(rho_smooth, dy, dx)
    grad_mag = np.sqrt(drho_dx**2 + drho_dy**2)

    reference = grad_mag[gas] if gas.any() else grad_mag
    norm = max(np.percentile(reference, GRAD_NORM_PERCENTILE), np.finfo(grad_mag.dtype).tiny)

    K = KL_REL * VF + (1.0 - VF)
    phi = np.exp(-CONTRAST * K * grad_mag / norm)
    return np.where(gas, phi, 1.0)


def load_frame(plt_path: Path, vtp_path: Path):
    ds = yt.load(str(plt_path))
    t = float(ds.current_time.to_value())

    # Cropped to the view at native finest-level resolution
    max_level = ds.index.max_level
    domain_le = ds.domain_left_edge.to_value("code_length")
    dcell = (ds.domain_right_edge.to_value("code_length")[0] - domain_le[0]) / (
        int(ds.domain_dimensions[0]) * 2**max_level
    )
    x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
    y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
    res_x = int(round(x_span / dcell))
    res_y = int(round(y_span / dcell))
    center = [0.5 * (VIEW_XLIM[0] + VIEW_XLIM[1]), 0.5 * (VIEW_YLIM[0] + VIEW_YLIM[1]), 0.0]

    slc = ds.slice("z", 0.0)
    frb = slc.to_frb((x_span, "code_length"), (res_x, res_y), center=center, height=(y_span, "code_length"))
    VF = np.array(frb["boxlib", "VF"])
    RHOL = np.array(frb["boxlib", "RHOL"])
    RHOG = np.array(frb["boxlib", "RHOG"])
    if SNAPSHOT_TIME is not None:
        level = np.array(frb["index", "grid_level"]).astype(int)
        data = compute_schlieren(VF, RHOL, RHOG, dcell, dcell, level=level, max_level=max_level)
    else:
        data = compute_schlieren(VF, RHOL, RHOG, dcell, dcell)

    segments = extract_plic_segments(vtp_path) if SHOW_INTERFACE else np.empty((0, 2, 2))
    return t, data, segments


frames_top = gather_frames(TOP_INPUT)
frames_bottom = gather_frames(BOTTOM_INPUT)

# Pair frames at common simulation times (a running case may not have reached the other's end time yet)
times_bottom = np.array([f[0] for f in frames_bottom])
times_top = np.array([f[0] for f in frames_top])
time_match_tol = TIME_MATCH_FRAC * (np.median(np.diff(times_top)) if len(times_top) > 1 else np.inf)
pairs = []
for top in frames_top:
    j = int(np.argmin(np.abs(times_bottom - top[0])))
    if abs(times_bottom[j] - top[0]) <= time_match_tol:
        pairs.append((top, frames_bottom[j]))
if not pairs:
    raise SystemExit(f"No common frame times between {TOP_INPUT} and {BOTTOM_INPUT}")
print(f"{len(pairs)} common frames, t = {pairs[0][0][0]:.2f} .. {pairs[-1][0][0]:.2f}")

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


def time_text(t):
    return rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}$"


def setup_row(ax, data0, segs0, label, show_xlabel):
    extent = [VIEW_XLIM[0], VIEW_XLIM[1], VIEW_YLIM[0], VIEW_YLIM[1]]
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

    # Horizontal tick numbers and axis label on the bottom row only
    if show_xlabel:
        ax.set_xlabel(r"$x$", labelpad=5.0)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    ax.set_ylabel(r"$y$", labelpad=-5.0)

    ax.set_title(label, fontsize=TITLE_FONTSIZE, pad=5.0)
    return im, lc


start = 0
if SNAPSHOT_TIME is not None:
    start = int(np.argmin([abs(top[0] - SNAPSHOT_TIME) for top, _ in pairs]))
(top0, bottom0) = pairs[start]
t0_top, data0_top, segs0_top = load_frame(top0[1], top0[2])
t0_bottom, data0_bottom, segs0_bottom = load_frame(bottom0[1], bottom0[2])

im_top, lc_top = setup_row(ax_top, data0_top, segs0_top, TOP_LABEL, show_xlabel=False)
im_bottom, lc_bottom = setup_row(ax_bottom, data0_bottom, segs0_bottom, BOTTOM_LABEL, show_xlabel=True)
# Paired frames share the simulation time: one time label above both rows
time_label = fig.text(LEFT_IN / FIG_WIDTH_IN + 0.5 * PLOT_WIDTH_IN / FIG_WIDTH_IN, 1.0 - 0.12 / FIG_HEIGHT_IN, time_text(t0_top),
                      ha="center", va="top", fontsize=TIME_FONTSIZE)

cbar = fig.colorbar(im_top, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 6))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)

if SNAPSHOT_TIME is not None:
    snapshot_output = f"schlieren_plic_compare_{TOP_INPUT}_{BOTTOM_INPUT}_t{t0_top:.2f}.pdf"
    fig.savefig(snapshot_output)
    print(f"Saved snapshot: {snapshot_output}  (t={t0_top:.4f})")
    raise SystemExit(0)


def update(i):
    top, bottom = pairs[i]
    t_top, data_top, segs_top = load_frame(top[1], top[2])
    t_bottom, data_bottom, segs_bottom = load_frame(bottom[1], bottom[2])

    im_top.set_data(data_top)
    lc_top.set_segments(segs_top)

    im_bottom.set_data(data_bottom)
    lc_bottom.set_segments(segs_bottom)
    time_label.set_text(time_text(t_top))

    return im_top, lc_top, im_bottom, lc_bottom, time_label


anim = FuncAnimation(fig, update, frames=len(pairs), blit=False)
anim.save(OUTPUT, fps=FPS, dpi=ANIMATION_DPI, writer="ffmpeg")
print(f"Saved: {OUTPUT}  ({len(pairs)} frames)")
