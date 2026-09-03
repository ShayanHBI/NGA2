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

# ============================================================================
# USER CONFIGURATION -- edit these and just run `python animate_field.py`
# ============================================================================

INPUT = "pTg"          # case identifier -> reads amrviz/impact_relax_<INPUT>

# Field to animate. One of:
#   "VF"        volume fraction                   (whole domain, range fixed 0-1)
#   "Yv"        vapor mass fraction               (whole domain, range fixed 0-1)
#   "RHOG"      gas density                       (gas-only,    blank where liquid)
#   "RHOL"      liquid density                    (liquid-only, blank where gas)
#   "TG"        gas temperature                   (gas-only,    blank where liquid)
#   "TL"        liquid temperature                (liquid-only, blank where gas)
#   "PG"        gas pressure                      (gas-only,    blank where liquid)
#   "PL"        liquid pressure                   (liquid-only, blank where gas)
#   "schlieren" numerical schlieren of mixture density (gas-phase wave structure only)
FIELD = "schlieren"

SHOW_INTERFACE = True   # overlay the PLIC interface contour
INTERFACE_COLOR = "white"
INTERFACE_LINEWIDTH = 0.75
# INTERFACE_LINEWIDTH = 0.025

START_TIME = 0        # skip frames before this simulation time (animation mode only; ignored for a snapshot)

# Set to a simulation time to render a single PDF of the nearest available
# frame instead of the full animation. Set to None for the full animation.
SNAPSHOT_TIME = None    # e.g. 15.0

# Colorbar range. Leave as None to auto-compute: fixed 0-1 for VF/Yv/schlieren,
# or the global min/max over the unmasked region across all frames otherwise.
# Override manually if the auto range is not what you want, e.g. VMAX_OVERRIDE = 30.0
VMIN_OVERRIDE = None
VMAX_OVERRIDE = None

# Case is non-dimensional -- plotted coordinates and time are used as-is.
VIEW_XLIM = (0.0, 11.0)
# VIEW_XLIM = (0.0, 3.5)
VIEW_YLIM = (-5.5, 5.5)

CMAP_OVERRIDE = None    # None -> "jet", except "gray" for schlieren

FPS = 12
ANIMATION_DPI = 200     # bump to ~300 for schlieren (finer gradient detail)

# --- schlieren-only options (ignored unless FIELD == "schlieren") ---------
# phi = exp(-CONTRAST*[KL_REL*VF + (1-VF)] * |grad(rho)| / NORM), rho the
# mixture density. NORM is set per frame to a high percentile of |grad(rho)|
# over the gas region so contrast adapts to whatever density/length scales
# and shock strength are present (this case is non-dimensional, so a fixed
# dimensional reference would not scale correctly).
KL = 120.0
KG = 30.0
KL_REL = KL / KG
CONTRAST = 1.0
GRAD_NORM_PERCENTILE = 99.9
# The liquid/gas density jump at the interface is orders of magnitude larger
# than any gas-phase shock gradient and would dominate the normalization if
# included. Cells with any liquid presence are excluded from the percentile
# and masked to phi=1 (blank) in the final image.
VF_MASK_THRESHOLD = 1.0e-3
# yt's to_frb pixelizes AMR data by replicating each cell's value at the
# finest level present (nearest-neighbor upsampling, not interpolation).
# Differentiating that directly picks up spurious AMR patch-boundary jumps,
# so rho is blurred first (normalized convolution restricted to gas-phase
# values, so the liquid's much larger density never leaks into the kernel
# average near the interface).
SMOOTH_SIGMA_PIX = 5.0

# ============================================================================
# End of user configuration
# ============================================================================

CASE = f"impact_relax_{INPUT}"
AMRVIZ_DIR = Path("amrviz") / CASE

# VF~0/VF~1 cells carry placeholder thermo values (PL/TL/RHOL where VF~0,
# PG/TG/RHOG where VF~1); mask them to NaN so they render blank instead of
# polluting the color scale.
VF_EPS = 1.0e-6

FIELD_INFO = {
    "VF":        dict(label=r"$\alpha$",    mask=None,    fixed_range=(0.0, 1.0)),
    "Yv":        dict(label=r"$Y_v$",       mask=None,    fixed_range=(0.0, 1.0)),
    "RHOG":      dict(label=r"$\rho_g$",    mask="gas",   fixed_range=None),
    "RHOL":      dict(label=r"$\rho_l$",    mask="liquid",fixed_range=None),
    "TG":        dict(label=r"$T_g$",       mask="gas",   fixed_range=None),
    "TL":        dict(label=r"$T_l$",       mask="liquid",fixed_range=None),
    "PG":        dict(label=r"$p_g$",       mask="gas",   fixed_range=None),
    "PL":        dict(label=r"$p_l$",       mask="liquid",fixed_range=None),
    "schlieren": dict(label=r"$\phi$",      mask=None,    fixed_range=(0.0, 1.0)),
}
if FIELD not in FIELD_INFO:
    raise ValueError(f"Unknown FIELD {FIELD!r}; choose one of {sorted(FIELD_INFO)}")
INFO = FIELD_INFO[FIELD]

CMAP = CMAP_OVERRIDE or ("gray" if FIELD == "schlieren" else "jet")
OUTPUT = f"{FIELD}_plic_{INPUT}.mp4"

# Paper-style figure geometry (inches). Canvas sized to ~9.5 x 7.0 cm --
# smaller physical size, so fonts/lines are scaled up to stay legible and
# margins are generous enough that the rotated y-label always fits.
CM_TO_IN         = 1.0 / 2.54
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
    return rf"${value:.2f}$"


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


def field_range(frames, field, mask_side):
    """Global min/max of `field` over the masked (physically meaningful)
    region across the whole series, so the colorbar stays fixed across
    frames instead of rescaling every timestep."""
    vmin, vmax = np.inf, -np.inf
    for plt_path, _ in frames:
        ds = yt.load(str(plt_path))
        ad = ds.all_data()
        vf = np.array(ad["boxlib", "VF"])
        data = np.array(ad["boxlib", field])
        if mask_side == "gas":
            mask = vf < 1.0 - VF_EPS
        elif mask_side == "liquid":
            mask = vf > VF_EPS
        else:
            mask = np.ones_like(vf, dtype=bool)
        if not np.any(mask):
            continue
        vmin = min(vmin, float(data[mask].min()))
        vmax = max(vmax, float(data[mask].max()))
    return vmin, vmax


def compute_schlieren(VF, RHOL, RHOG, dx, dy):
    rho = VF * RHOL + (1.0 - VF) * RHOG
    gas = VF < VF_MASK_THRESHOLD

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

    if FIELD == "schlieren":
        # Cropped to VIEW_XLIM/VIEW_YLIM at native finest-level resolution so
        # the plotted array's extent matches the axes exactly -- matplotlib's
        # PDF backend then embeds it losslessly regardless of savefig dpi.
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
        data = compute_schlieren(VF, RHOL, RHOG, dcell, dcell)
        extent = [VIEW_XLIM[0], VIEW_XLIM[1], VIEW_YLIM[0], VIEW_YLIM[1]]
    else:
        max_level = ds.index.max_level
        res = int(ds.domain_dimensions[0]) * 2**max_level
        le = ds.domain_left_edge.to_value("code_length")
        re = ds.domain_right_edge.to_value("code_length")

        slc = ds.slice("z", 0.0)
        frb = slc.to_frb((re[0] - le[0], "code_length"), res, height=(re[1] - le[1], "code_length"))
        data = np.array(frb["boxlib", FIELD])
        if INFO["mask"] is not None:
            vf = np.array(frb["boxlib", "VF"])
            keep = (vf < 1.0 - VF_EPS) if INFO["mask"] == "gas" else (vf > VF_EPS)
            data = np.where(keep, data, np.nan)
        extent = [le[0], re[0], le[1], re[1]]

    segments = extract_plic_segments(vtp_path) if SHOW_INTERFACE else np.empty((0, 2, 2))
    return t, data, segments, extent


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

if SNAPSHOT_TIME is not None:
    # Binary search directly on the full, unfiltered frame list -- touches
    # only O(log n) plotfiles instead of scanning every timestamp. START_TIME
    # only matters for trimming the animation loop, so it is skipped here.
    frame_idx, _ = closest_frame_index(frames, SNAPSHOT_TIME)
else:
    if START_TIME > 0.0:
        frames = [
            (p, vp) for p, vp in frames
            if float(yt.load(str(p)).current_time.to_value()) >= START_TIME
        ]
    frame_idx = 0

t0, data0, segs0, extent = load_frame(*frames[frame_idx])

if VMIN_OVERRIDE is not None and VMAX_OVERRIDE is not None:
    VMIN, VMAX = VMIN_OVERRIDE, VMAX_OVERRIDE
elif INFO["fixed_range"] is not None:
    VMIN, VMAX = INFO["fixed_range"]
elif SNAPSHOT_TIME is not None:
    # Snapshot mode: range the colorbar off the single displayed frame
    # (already masked by load_frame) rather than scanning the whole series.
    VMIN, VMAX = float(np.nanmin(data0)), float(np.nanmax(data0))
    if VMIN_OVERRIDE is not None:
        VMIN = VMIN_OVERRIDE
    if VMAX_OVERRIDE is not None:
        VMAX = VMAX_OVERRIDE
else:
    VMIN, VMAX = field_range(frames, FIELD, INFO["mask"])
    if VMIN_OVERRIDE is not None:
        VMIN = VMIN_OVERRIDE
    if VMAX_OVERRIDE is not None:
        VMAX = VMAX_OVERRIDE

x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
data_height_over_width = y_span / x_span
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

ax.set_xlabel(r"$x$", labelpad=5.0)
ax.set_ylabel(r"$y$", labelpad=-5.0)
title = ax.set_title(
    rf"$t = {f'{t0:.2f}'.rstrip('0').rstrip('.')}$", fontsize=TITLE_FONTSIZE, pad=6.0
)

cbar = fig.colorbar(im, cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(INFO["label"], rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)


def update(i):
    plt_path, vtp_path = frames[i]
    t, data, segments, _ = load_frame(plt_path, vtp_path)
    im.set_data(data)
    lc.set_segments(segments)
    title.set_text(rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}$")
    return im, lc, title


if SNAPSHOT_TIME is not None:
    snapshot_output = f"{FIELD}_plic_{INPUT}_t{t0:.2f}.pdf"
    fig.savefig(snapshot_output)
    print(f"Saved snapshot: {snapshot_output}  (t={t0:.4f})")
else:
    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    anim.save(OUTPUT, fps=FPS, dpi=ANIMATION_DPI, writer="ffmpeg")
    print(f"Saved: {OUTPUT}  ({len(frames)} frames)")
