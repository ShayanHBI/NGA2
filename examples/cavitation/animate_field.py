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

INPUT = "wall_pTg"     # case identifier -> reads amrviz/cavitation_<INPUT>

# Field to animate. One of:
#   "VF"        liquid volume fraction            (whole domain, range fixed 0-1)
#   "Yv"        vapor mass fraction               (whole domain, range fixed 0-1)
#   "RHOG"      gas density                       (gas-only,    blank where liquid)
#   "RHOL"      liquid density                    (liquid-only, blank where gas)
#   "TG"        gas temperature                   (gas-only,    blank where liquid)
#   "TL"        liquid temperature                (liquid-only, blank where gas)
#   "PG"        gas pressure                      (gas-only,    blank where liquid)
#   "PL"        liquid pressure                   (liquid-only, blank where gas)
#   "Umag"      velocity magnitude                (whole domain)
#   "Mach"      Mach number                       (whole domain)
#   "schlieren" numerical schlieren of mixture density (SCHLIEREN_PHASE only)
FIELD = "VF"

SHOW_INTERFACE = True   # overlay the PLIC interface
INTERFACE_COLOR = "white"
INTERFACE_LINEWIDTH = 1.2

START_TIME = 0.0        # [mus] skip frames before this time (animation mode only)

# [mus] render a single PDF of the nearest frame instead of the animation; None for the animation
SNAPSHOT_TIME = None

# Colorbar range; None -> fixed 0-1 for VF/Yv/schlieren, else global min/max over the unmasked region
VMIN_OVERRIDE = None
VMAX_OVERRIDE = None

# View limits and all plotted coordinates are in mm, time in mus (raw plotfile/vtp data is in m and s).
M_TO_MM = 1.0e3
S_TO_MUS = 1.0e6
VIEW_XLIM = (-3.0, 3.0)
# VIEW_XLIM = (-15.0, 15.0)
VIEW_YLIM = (-3.0, 3.0)

CMAP_OVERRIDE = None    # None -> "jet", except "gray" for schlieren

FPS = 12
ANIMATION_DPI = 200     # bump to ~300 for schlieren (finer gradient detail)

# --- schlieren-only options (ignored unless FIELD == "schlieren") ---------
# phi = exp(-CONTRAST*|grad(rho)|/NORM), NORM a per-frame percentile of |grad(rho)| over SCHLIEREN_PHASE
SCHLIEREN_PHASE = "liquid"   # "liquid" (pressure waves in the water) or "gas"
CONTRAST = 1.0
GRAD_NORM_PERCENTILE = 99.9
VF_MASK_THRESHOLD = 1.0e-3   # cells with more than this of the other phase are blanked
SMOOTH_SIGMA_PIX = 2.0       # blur before differentiating to suppress AMR patch-boundary jumps

# ============================================================================
# End of user configuration
# ============================================================================

CASE = f"cavitation_{INPUT}"
AMRVIZ_DIR = Path("amrviz") / CASE

# VF~0/VF~1 cells carry placeholder thermo values; mask them to NaN so they render blank.
VF_EPS = 1.0e-6

FIELD_INFO = {
    "VF":        dict(label=r"$\alpha$",                                   mask=None,     scale=1.0,    fixed_range=(0.0, 1.0)),
    "Yv":        dict(label=r"$Y_v$",                                      mask=None,     scale=1.0,    fixed_range=(0.0, 1.0)),
    "RHOG":      dict(label=r"$\rho_g\;\left(\mathrm{kg\,m^{-3}}\right)$", mask="gas",    scale=1.0,    fixed_range=None),
    "RHOL":      dict(label=r"$\rho_l\;\left(\mathrm{kg\,m^{-3}}\right)$", mask="liquid", scale=1.0,    fixed_range=None),
    "TG":        dict(label=r"$T_g\;\left(\mathrm{K}\right)$",             mask="gas",    scale=1.0,    fixed_range=None),
    "TL":        dict(label=r"$T_l\;\left(\mathrm{K}\right)$",             mask="liquid", scale=1.0,    fixed_range=None),
    "PG":        dict(label=r"$p_g\;\left(\mathrm{MPa}\right)$",           mask="gas",    scale=1.0e-6, fixed_range=None),
    "PL":        dict(label=r"$p_l\;\left(\mathrm{MPa}\right)$",           mask="liquid", scale=1.0e-6, fixed_range=None),
    "Umag":      dict(label=r"$|\vec{U}|\;\left(\mathrm{m\,s^{-1}}\right)$", mask=None,   scale=1.0,    fixed_range=None),
    "Mach":      dict(label=r"$\mathrm{Mach}$",                            mask=None,     scale=1.0,    fixed_range=None),
    "schlieren": dict(label=r"$\phi$",                                     mask=None,     scale=1.0,    fixed_range=(0.0, 1.0)),
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
    # Fixed-point integers above 100 so MPa/K/kg m^-3 labels stay narrow
    if abs(value) >= 100:
        return rf"${value:.0f}$"
    return rf"${value:.2g}$"


def time_label(t):
    return rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}\;\mu\mathrm{{s}}$"


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
    return segments * M_TO_MM


def frame_time(plt_path: Path):
    return float(yt.load(str(plt_path)).current_time.to_value()) * S_TO_MUS


def closest_frame_index(frames, target_time):
    """Binary search on frame time (monotonic in frame number) so only a
    handful of plotfiles are opened instead of the whole series."""
    cache = {}

    def time_at(i):
        if i not in cache:
            cache[i] = frame_time(frames[i][0])
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


def phase_mask(vf, mask_side):
    if mask_side == "gas":
        return vf < 1.0 - VF_EPS
    if mask_side == "liquid":
        return vf > VF_EPS
    return np.ones_like(vf, dtype=bool)


def finite_range(vmin, vmax):
    # Fall back to 0-1 when the field is fully masked (e.g. gas fields before cavitation)
    if not (np.isfinite(vmin) and np.isfinite(vmax)):
        return 0.0, 1.0
    return vmin, vmax


def field_range(frames, field, mask_side, scale):
    """Global min/max of `field` over the masked (physically meaningful)
    region across the whole series, so the colorbar stays fixed across
    frames instead of rescaling every timestep."""
    vmin, vmax = np.inf, -np.inf
    for plt_path, _ in frames:
        ad = yt.load(str(plt_path)).all_data()
        vf = np.array(ad["boxlib", "VF"])
        data = np.array(ad["boxlib", field]) * scale
        mask = phase_mask(vf, mask_side)
        if not np.any(mask):
            continue
        vmin = min(vmin, float(data[mask].min()))
        vmax = max(vmax, float(data[mask].max()))
    return finite_range(vmin, vmax)


def compute_schlieren(VF, RHOL, RHOG, dx, dy):
    rho = VF * RHOL + (1.0 - VF) * RHOG
    keep = (VF > 1.0 - VF_MASK_THRESHOLD) if SCHLIEREN_PHASE == "liquid" else (VF < VF_MASK_THRESHOLD)

    # Normalized convolution restricted to the kept phase, so the other phase's density never leaks in
    weighted = gaussian_filter(np.where(keep, rho, 0.0), sigma=SMOOTH_SIGMA_PIX)
    weight = gaussian_filter(keep.astype(rho.dtype), sigma=SMOOTH_SIGMA_PIX)
    rho_smooth = weighted / np.maximum(weight, np.finfo(rho.dtype).tiny)

    drho_dy, drho_dx = np.gradient(rho_smooth, dy, dx)
    grad_mag = np.sqrt(drho_dx**2 + drho_dy**2)

    reference = grad_mag[keep] if keep.any() else grad_mag
    norm = max(np.percentile(reference, GRAD_NORM_PERCENTILE), np.finfo(grad_mag.dtype).tiny)

    phi = np.exp(-CONTRAST * grad_mag / norm)
    return np.where(keep, phi, 1.0)


def load_frame(plt_path: Path, vtp_path: Path):
    ds = yt.load(str(plt_path))
    t = float(ds.current_time.to_value()) * S_TO_MUS

    # Crop to the view at native finest-level resolution, snapped outward to cell edges so each pixel is one leaf cell
    ncells = np.array(ds.domain_dimensions[:2], dtype=int) * 2**ds.index.max_level
    dom_lo = ds.domain_left_edge.to_value("code_length")[:2] * M_TO_MM
    dom_hi = ds.domain_right_edge.to_value("code_length")[:2] * M_TO_MM
    dcell = (dom_hi - dom_lo) / ncells
    view_lo = np.array([VIEW_XLIM[0], VIEW_YLIM[0]])
    view_hi = np.array([VIEW_XLIM[1], VIEW_YLIM[1]])
    ilo = np.clip(np.floor((view_lo - dom_lo) / dcell), 0, ncells).astype(int)
    ihi = np.clip(np.ceil((view_hi - dom_lo) / dcell), 0, ncells).astype(int)
    lo = dom_lo + ilo * dcell
    hi = dom_lo + ihi * dcell
    center = [0.5 * (lo[0] + hi[0]) / M_TO_MM, 0.5 * (lo[1] + hi[1]) / M_TO_MM, 0.0]

    slc = ds.slice("z", 0.0)
    frb = slc.to_frb(
        ((hi[0] - lo[0]) / M_TO_MM, "code_length"), (int(ihi[0] - ilo[0]), int(ihi[1] - ilo[1])),
        center=center, height=((hi[1] - lo[1]) / M_TO_MM, "code_length"),
    )
    VF = np.array(frb["boxlib", "VF"])

    if FIELD == "schlieren":
        RHOL = np.array(frb["boxlib", "RHOL"])
        RHOG = np.array(frb["boxlib", "RHOG"])
        data = compute_schlieren(VF, RHOL, RHOG, dcell[0], dcell[1])
    else:
        data = np.array(frb["boxlib", FIELD]) * INFO["scale"]
        if INFO["mask"] is not None:
            data = np.where(phase_mask(VF, INFO["mask"]), data, np.nan)
    extent = [lo[0], hi[0], lo[1], hi[1]]

    segments = extract_plic_segments(vtp_path) if SHOW_INTERFACE else np.empty((0, 2, 2))
    return t, data, segments, extent


plt_files = sorted(AMRVIZ_DIR.glob("plt.nga2.cell.*"), key=frame_number)
vtp_files = {frame_number(p): p for p in AMRVIZ_DIR.glob("plic_*.vtp")}
frames = [(p, vtp_files[frame_number(p)]) for p in plt_files if frame_number(p) in vtp_files]

if SNAPSHOT_TIME is not None:
    # O(log n) plotfiles touched; START_TIME only trims the animation loop
    frame_idx, _ = closest_frame_index(frames, SNAPSHOT_TIME)
else:
    if START_TIME > 0.0:
        frames = [(p, vp) for p, vp in frames if frame_time(p) >= START_TIME]
    frame_idx = 0

t0, data0, segs0, extent = load_frame(*frames[frame_idx])

if VMIN_OVERRIDE is not None and VMAX_OVERRIDE is not None:
    VMIN, VMAX = VMIN_OVERRIDE, VMAX_OVERRIDE
else:
    if INFO["fixed_range"] is not None:
        VMIN, VMAX = INFO["fixed_range"]
    elif SNAPSHOT_TIME is not None:
        # Snapshot mode: range off the single (already masked) frame
        VMIN, VMAX = finite_range(float(np.nanmin(data0)), float(np.nanmax(data0)))
    else:
        VMIN, VMAX = field_range(frames, FIELD, INFO["mask"], INFO["scale"])
    if VMIN_OVERRIDE is not None:
        VMIN = VMIN_OVERRIDE
    if VMAX_OVERRIDE is not None:
        VMAX = VMAX_OVERRIDE

x_span = VIEW_XLIM[1] - VIEW_XLIM[0]
y_span = VIEW_YLIM[1] - VIEW_YLIM[0]
data_height_over_width = y_span / x_span
tick_locator = MaxNLocator(nbins=5, steps=[1, 2, 5, 10])
x_ticks = [v for v in tick_locator.tick_values(*VIEW_XLIM) if VIEW_XLIM[0] <= v <= VIEW_XLIM[1]]
y_ticks = [v for v in tick_locator.tick_values(*VIEW_YLIM) if VIEW_YLIM[0] <= v <= VIEW_YLIM[1]]

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

ax.set_xlabel(r"$x\;\left(\mathrm{mm}\right)$", labelpad=5.0)
ax.set_ylabel(r"$y\;\left(\mathrm{mm}\right)$", labelpad=-5.0)
title = ax.set_title(time_label(t0), fontsize=TITLE_FONTSIZE, pad=6.0)

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
    title.set_text(time_label(t))
    return im, lc, title


if SNAPSHOT_TIME is not None:
    snapshot_output = f"{FIELD}_plic_{INPUT}_t{t0:.2f}.pdf"
    fig.savefig(snapshot_output)
    print(f"Saved snapshot: {snapshot_output}  (t={t0:.4f} mus)")
else:
    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    anim.save(OUTPUT, fps=FPS, dpi=ANIMATION_DPI, writer="ffmpeg")
    print(f"Saved: {OUTPUT}  ({len(frames)} frames)")
