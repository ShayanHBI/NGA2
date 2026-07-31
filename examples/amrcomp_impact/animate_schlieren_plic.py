import argparse
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

parser = argparse.ArgumentParser()
parser.add_argument("input", help="case identifier, e.g. pTg")
parser.add_argument(
    "--t", type=float, default=None,
    help="Render a single snapshot at (nearest available frame to) this time, instead of the full animation",
)
args = parser.parse_args()
SNAPSHOT = args.t is not None

INPUT = args.input
CASE = f"impact_relax_{INPUT}"
AMRVIZ_DIR = Path("amrviz") / CASE

# Numerical schlieren: phi = exp(-CONTRAST*[Kl_rel*VF + (1-VF)] * |grad(rho)| / NORM),
# with rho = VF*RHOL + (1-VF)*RHOG the mixture density. NORM is not a fixed
# constant (this case is non-dimensional, so a dimensional reference gradient
# would not scale correctly) -- instead it is set per frame to a high
# percentile of |grad(rho)| over the field, so contrast adapts automatically
# to whatever density/length scales and shock strength are present.
#
# Kl=120, Kg=30 (the paper's values) were meant to be used as the absolute
# exponent scale against a *fixed*, dimensional NORM=2.5e5, where typical
# gradients sit far below the reference and only true extremes approach
# exponent~O(1). Against a percentile-based NORM, the reference gradient IS
# the exponent scale by construction, so using Kg=30 directly there means any
# gas pixel near the reference percentile gets exp(-30)~0 -- solid black
# instead of graded gray. Kl_rel = KL/KG keeps their relative liquid:gas
# sensitivity (4:1), while CONTRAST is the actual O(1) exponent scale, tuned
# so the reference percentile itself lands at a visible mid-gray rather than
# fully saturating.
KL = 120.0
KG = 30.0
KL_REL = KL / KG
CONTRAST = 1.0
GRAD_NORM_PERCENTILE = 99.9

# The liquid/gas density jump at the interface is orders of magnitude larger
# than any gas-phase shock gradient. Left in, it would dominate the contrast
# normalization (everything else reads as faded). The PLIC overlay already
# marks the interface unambiguously, so cells with any liquid presence are
# excluded from the normalization percentile and masked to phi=1 (blank) in
# the final image -- schlieren is shown only for the surrounding gas-phase
# wave structure. (KL has no visible effect under this mask; it only matters
# if VF_MASK_THRESHOLD is raised enough to let partially-liquid cells render.)
VF_MASK_THRESHOLD = 1.0e-3

# yt's to_frb pixelizes AMR data by replicating each cell's value across the
# pixels it covers at the finest level present (nearest-neighbor upsampling of
# coarser patches, not interpolation). Differentiating that field directly
# picks up spurious jumps at every AMR patch/cell boundary -- a checkerboard
# artifact of the mesh, not a physical gradient. A small Gaussian blur of rho
# (scale in pixels, smaller than any real shock-front width) removes it.
#
# ds.smoothed_covering_grid was tried as a fix (it interpolates across AMR
# levels instead of replicating, which should avoid the artifact without any
# blur) but produces outright WRONG data for this dataset: the VF field comes
# back with a large spurious blocky region of false liquid presence well
# beyond the actual droplet (verified by comparing VF masks side by side --
# to_frb gives the correct circular drop, smoothed_covering_grid gives a
# corrupted "Pac-Man" shape). That is a correctness bug, not a cosmetic one,
# so to_frb + blur is used despite the softening it costs on real features.
#
# A plain blur, however, does not respect the gas/liquid mask: the kernel
# still averages the drop's (huge) density into gas cells within ~SIGMA
# pixels of the interface, which then show a fake gradient there even though
# VF_MASK_THRESHOLD excludes them from the liquid side. That leakage is what
# turns into a thick halo around the drop. compute_schlieren blurs with a
# normalized convolution restricted to gas-phase values instead (blur the
# gas-masked field and the mask itself, then divide), so the liquid density
# never enters the kernel average near the boundary.
SMOOTH_SIGMA_PIX = 5.0

CBAR_LABEL = r"$\phi$"
VMIN, VMAX = 0.0, 1.0

# Case is non-dimensional -- plotted coordinates and time are used as-is.
VIEW_XLIM = (0.0, 11.0)
VIEW_YLIM = (-7, 7)

CMAP = "gray"
INTERFACE_COLOR = "red"
INTERFACE_LINEWIDTH = 0.25

FPS = 12
OUTPUT = f"schlieren_plic_{INPUT}.mp4"

# load_frame requests the FRB pre-cropped to VIEW_XLIM/VIEW_YLIM (see below),
# so the plotted array's extent matches the axes xlim/ylim exactly -- same as
# plot_initial.py's Mach0.pdf. That lets matplotlib's PDF backend embed the
# image losslessly at its own native pixel size, independent of savefig dpi,
# so the snapshot passes no dpi at all (same as plot_initial.py). The
# animation renders to a raster mp4 via ffmpeg, which has no such lossless
# path -- its dpi is a real quality/file-size/render-time tradeoff.
ANIMATION_DPI = 300

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


def get_time(plt_path: Path) -> float:
    return float(yt.load(str(plt_path)).current_time.to_value())


def nearest_frame_by_time(frames, target_t):
    """Binary search on frame index for the nearest simulation time, relying on
    time increasing monotonically with frame index. Avoids a yt.load (which
    reads the full AMR hierarchy metadata) per frame -- O(log n) loads instead
    of O(n), which matters since --snapshot only needs a single frame."""
    lo, hi = 0, len(frames) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if get_time(frames[mid][0]) < target_t:
            lo = mid + 1
        else:
            hi = mid
    if lo == 0:
        return 0
    t_lo = get_time(frames[lo][0])
    t_prev = get_time(frames[lo - 1][0])
    return lo if abs(t_lo - target_t) < abs(t_prev - target_t) else lo - 1


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


def compute_schlieren(VF, RHOL, RHOG, dx, dy):
    rho = VF * RHOL + (1.0 - VF) * RHOG
    gas = VF < VF_MASK_THRESHOLD

    # Normalized convolution: blur rho and the gas mask separately, then
    # divide, so the kernel average near the interface is only ever informed
    # by gas-phase values (a plain blur of rho would leak the liquid's much
    # larger density into nearby gas cells and fake a gradient there).
    weighted = gaussian_filter(np.where(gas, rho, 0.0), sigma=SMOOTH_SIGMA_PIX)
    weight = gaussian_filter(gas.astype(rho.dtype), sigma=SMOOTH_SIGMA_PIX)
    rho_smooth = weighted / np.maximum(weight, np.finfo(rho.dtype).tiny)

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

    # Request the FRB cropped to VIEW_XLIM/VIEW_YLIM at native finest-level
    # resolution, rather than the full domain (which is what plot_initial.py
    # does for Mach0.pdf). When the axes extent matches xlim/ylim exactly,
    # matplotlib's PDF backend embeds the array losslessly at its own native
    # pixel size regardless of savefig dpi; if the array instead spans a
    # wider region than the view (as it did before, full domain vs. a
    # VIEW_XLIM subset), the PDF backend must resample it to the view, which
    # then goes through the dpi-dependent (and lossy) rasterization path.
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

    phi = compute_schlieren(VF, RHOL, RHOG, dcell, dcell)

    segments = extract_plic_segments(vtp_path)
    return t, phi, segments


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

if SNAPSHOT:
    frames = [frames[nearest_frame_by_time(frames, args.t)]]

t0, data0, segs0 = load_frame(*frames[0])

x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
data_height_over_width = y_span / x_span
extent = [VIEW_XLIM[0], VIEW_XLIM[1], VIEW_YLIM[0], VIEW_YLIM[1]]
tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
x_ticks = tick_locator.tick_values(*VIEW_XLIM)
x_ticks = [t for t in x_ticks if VIEW_XLIM[0] <= t <= VIEW_XLIM[1]]
y_ticks = tick_locator.tick_values(*VIEW_YLIM)
y_ticks = [t for t in y_ticks if VIEW_YLIM[0] <= t <= VIEW_YLIM[1]]

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

im = ax.imshow(
    data0, origin="lower", extent=extent, cmap=CMAP, vmin=VMIN, vmax=VMAX, interpolation="none"
)
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

ax.set_xlabel(r"$x$", labelpad=5.0)
ax.set_ylabel(r"$y$", labelpad=-5.0)
title = ax.set_title(
    rf"$t = {f'{t0:.2f}'.rstrip('0').rstrip('.')}$", fontsize=TITLE_FONTSIZE, pad=6.0
)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(CBAR_LABEL, rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)


def update(i):
    plt_path, vtp_path = frames[i]
    t, data, segments = load_frame(plt_path, vtp_path)
    im.set_data(data)
    lc.set_segments(segments)
    title.set_text(rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}$")
    return im, lc, title


if SNAPSHOT:
    snapshot_output = f"schlieren_plic_{INPUT}_t{t0:.2f}.pdf"
    fig.savefig(snapshot_output)
    print(f"Saved snapshot: {snapshot_output}  (t={t0:.3f})")
else:
    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    anim.save(OUTPUT, fps=FPS, dpi=ANIMATION_DPI, writer="ffmpeg")
    print(f"Saved: {OUTPUT}  ({len(frames)} frames)")
