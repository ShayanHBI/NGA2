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

INPUT = "pTg_finite"     # case identifier -> reads amrviz/impact_relax_<INPUT> (<type>_inst or <type>_finite)

# Row-by-row comparison: COMPARE = True stacks COMPARE_CASES top to bottom with a
# shared colorbar (INPUT is then ignored). Frames are paired at common simulation
# times, so a case still running limits the comparison to the times both have.
COMPARE = True
COMPARE_CASES = [
    ("pTg_inst",   r"Instantaneous MTC"),
    ("pTg_finite", r"Instantaneous MT, finite C"),
]

# Field to animate. One of:
#   "VF"        volume fraction                   (whole domain, range fixed 0-1)
#   "Yv"        vapor mass fraction               (whole domain, range fixed 0-1)
#   "RHOG"      gas density                       (gas-only,    blank where liquid)
#   "RHOL"      liquid density                    (liquid-only, blank where gas)
#   "TG"        gas temperature                   (gas-only,    blank where liquid)
#   "TL"        liquid temperature                (liquid-only, blank where gas)
#   "PG"        gas pressure                      (gas-only,    blank where liquid)
#   "PL"        liquid pressure                   (liquid-only, blank where gas)
#   "RHO"       mixture density     VF*RHOL+(1-VF)*RHOG  (whole domain)
#   "P"         mixture pressure    VF*PL  +(1-VF)*PG    (whole domain)
#   "T"         mixture temperature VF*TL  +(1-VF)*TG    (whole domain)
#   "schlieren" numerical schlieren of mixture density (gas-phase wave structure only)
FIELD = "VF"

SHOW_INTERFACE = True   # overlay the interface contour
INTERFACE_COLOR = "black"
# INTERFACE_LINEWIDTH = 0.75
INTERFACE_LINEWIDTH = 0.2

START_TIME = 0        # skip frames before this simulation time (animation mode only; ignored for a snapshot)

# Set to a simulation time to render a single PDF of the nearest available
# frame instead of the full animation. Set to None for the full animation.
SNAPSHOT_TIME = 11    # e.g. 15.0

# Colorbar range. Leave as None to auto-compute: fixed 0-1 for VF/Yv/schlieren,
# or the global min/max over the unmasked region across all frames otherwise.
# Override manually if the auto range is not what you want, e.g. VMAX_OVERRIDE = 30.0
VMIN_OVERRIDE = None
VMAX_OVERRIDE = None

# Case is non-dimensional -- plotted coordinates and time are used as-is.
# VIEW_XLIM = (0, 11)
# VIEW_YLIM = (-5.5, 5.5)

VIEW_XLIM = (0, 2)
VIEW_YLIM = (-2.5, 2.5)

CMAP_OVERRIDE = None    # None -> "jet", except "gray" for schlieren

FPS = 12
ANIMATION_DPI = 200     # minimum; raised automatically so each panel gets at least one pixel per finest cell

# --- schlieren-only options (ignored unless FIELD == "schlieren") ---------
# phi = exp(-CONTRAST*[KL_REL*VF + (1-VF)] * |grad(rho)| / NORM), rho the
# mixture density. NORM is set per frame to a high percentile of |grad(rho)|
# over the gas region so contrast adapts to whatever density/length scales
# and shock strength are present (this case is non-dimensional, so a fixed
# dimensional reference would not scale correctly).
KL = 120.0
KG = 30.0
KL_REL = KL / KG
CONTRAST = 10.0
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
# Maximum resolution: finest-level pixels unsmoothed, level-L pixels blurred by
# SMOOTH_CELLS of their own cell size (just enough to remove the staircase)
SMOOTH_CELLS = 0.5

# ============================================================================
# End of user configuration
# ============================================================================

CASES = COMPARE_CASES if COMPARE else [(INPUT, None)]
NROWS = len(CASES)

# VF~0/VF~1 cells carry placeholder thermo values (PL/TL/RHOL where VF~0,
# PG/TG/RHOG where VF~1); mask them to NaN so they render blank instead of
# polluting the color scale.
VF_EPS = 1.0e-6

# Frames of compared cases are paired by simulation time within this fraction of the
# output interval (outputs land on the first step past each multiple, so times can differ by ~dt).
TIME_MATCH_FRAC = 0.25

# Mixture fields: (liquid, gas) phasic fields weighted by VF and 1-VF
MIXTURE = {"RHO": ("RHOL", "RHOG"), "P": ("PL", "PG"), "T": ("TL", "TG")}

FIELD_INFO = {
    "VF":        dict(label=r"$\alpha$",    mask=None,    fixed_range=(0.0, 1.0)),
    "Yv":        dict(label=r"$Y_v$",       mask=None,    fixed_range=(0.0, 1.0)),
    "RHOG":      dict(label=r"$\rho_g$",    mask="gas",   fixed_range=None),
    "RHOL":      dict(label=r"$\rho_l$",    mask="liquid",fixed_range=None),
    "TG":        dict(label=r"$T_g$",       mask="gas",   fixed_range=None),
    "TL":        dict(label=r"$T_l$",       mask="liquid",fixed_range=None),
    "PG":        dict(label=r"$p_g$",       mask="gas",   fixed_range=None),
    "PL":        dict(label=r"$p_l$",       mask="liquid",fixed_range=None),
    "RHO":       dict(label=r"$\rho$",      mask=None,    fixed_range=None),
    "P":         dict(label=r"$p$",         mask=None,    fixed_range=None),
    "T":         dict(label=r"$T$",         mask=None,    fixed_range=None),
    "schlieren": dict(label=r"$\phi$",      mask=None,    fixed_range=(0.0, 1.0)),
}
if FIELD not in FIELD_INFO:
    raise ValueError(f"Unknown FIELD {FIELD!r}; choose one of {sorted(FIELD_INFO)}")
INFO = FIELD_INFO[FIELD]

CMAP = CMAP_OVERRIDE or ("gray" if FIELD == "schlieren" else "jet")
OUTPUT_TAG = "compare_" + "_".join(c for c, _ in CASES) if COMPARE else INPUT
OUTPUT = f"{FIELD}_plic_{OUTPUT_TAG}.mp4"

# Paper-style figure geometry (inches). Canvas sized to ~9.5 x 7.0 cm --
# smaller physical size, so fonts/lines are scaled up to stay legible and
# margins are generous enough that the rotated y-label always fits.
CM_TO_IN         = 1.0 / 2.54
FIG_WIDTH_IN     = 9.5 * CM_TO_IN
LEFT_MARGIN_IN   = 0.78
RIGHT_MARGIN_IN  = 0.75
BOTTOM_MARGIN_IN = 0.62
TOP_MARGIN_IN    = 0.62 if COMPARE else 0.40   # compare: room for the time label above the row labels
ROW_GAP_IN       = 0.30
CBAR_GAP_IN      = 0.13
CBAR_WIDTH_IN    = 0.17
TICK_PAD_PT      = 7.0

AXES_LABEL_FONTSIZE = 16.0
NUMBER_FONTSIZE     = 13.0
CBAR_LABEL_FONTSIZE = 14.0
TITLE_FONTSIZE      = 15.0
ROW_TITLE_FONTSIZE  = 11.0
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
    return rf"${value:.2f}$"


def time_text(t):
    return rf"$t = {f'{t:.2f}'.rstrip('0').rstrip('.')}$"


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
    if not frames:
        raise SystemExit(f"No frames found in {amrviz_dir}")
    frames.sort(key=lambda f: f[0])
    return frames


def get_field(src, field):
    """Field from a yt data container or FRB; mixture fields are VF-weighted phasic values."""
    if field in MIXTURE:
        vf = np.array(src["boxlib", "VF"])
        liq = np.array(src["boxlib", MIXTURE[field][0]])
        gas = np.array(src["boxlib", MIXTURE[field][1]])
        # The placeholder value of an absent phase gets zero weight
        return np.where(vf > VF_EPS, vf * liq, 0.0) + np.where(vf < 1.0 - VF_EPS, (1.0 - vf) * gas, 0.0)
    return np.array(src["boxlib", field])


def field_range(steps, field, mask_side):
    """Global min/max of `field` over the masked (physically meaningful)
    region across the whole series and all compared cases, so the colorbar
    stays fixed across frames instead of rescaling every timestep."""
    vmin, vmax = np.inf, -np.inf
    for row in steps:
        for _, plt_path, _ in row:
            ds = yt.load(str(plt_path))
            ad = ds.all_data()
            vf = np.array(ad["boxlib", "VF"])
            data = get_field(ad, field)
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


def gas_smooth(rho, gas, sigma):
    """Gaussian blur restricted to gas pixels (normalized convolution), so the
    liquid's much larger density never leaks into the average."""
    if sigma <= 0.0:
        # Keep gas pixels as they are; liquid pixels take the nearby gas average so the gradient never sees the liquid
        return np.where(gas, rho, gas_smooth(rho, gas, 1.0))
    weighted = gaussian_filter(np.where(gas, rho, 0.0), sigma=sigma)
    weight = gaussian_filter(gas.astype(rho.dtype), sigma=sigma)
    return weighted / np.maximum(weight, np.finfo(rho.dtype).tiny)


def compute_schlieren(VF, RHOL, RHOG, dx, dy, level, max_level):
    rho = VF * RHOL + (1.0 - VF) * RHOG
    gas = VF < VF_MASK_THRESHOLD

    rho_smooth = np.empty_like(rho)
    for lev in np.unique(level):
        sigma = SMOOTH_CELLS * 2.0**(max_level - lev) if lev < max_level else 0.0
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
        level = np.array(frb["index", "grid_level"]).astype(int)
        data = compute_schlieren(VF, RHOL, RHOG, dcell, dcell, level=level, max_level=max_level)
        extent = [VIEW_XLIM[0], VIEW_XLIM[1], VIEW_YLIM[0], VIEW_YLIM[1]]
    else:
        max_level = ds.index.max_level
        res = int(ds.domain_dimensions[0]) * 2**max_level
        le = ds.domain_left_edge.to_value("code_length")
        re = ds.domain_right_edge.to_value("code_length")

        slc = ds.slice("z", 0.0)
        frb = slc.to_frb((re[0] - le[0], "code_length"), res, height=(re[1] - le[1], "code_length"))
        data = get_field(frb, FIELD)
        if INFO["mask"] is not None:
            vf = np.array(frb["boxlib", "VF"])
            keep = (vf < 1.0 - VF_EPS) if INFO["mask"] == "gas" else (vf > VF_EPS)
            data = np.where(keep, data, np.nan)
        extent = [le[0], re[0], le[1], re[1]]

    segments = extract_plic_segments(vtp_path) if SHOW_INTERFACE else np.empty((0, 2, 2))
    return t, data, segments, extent


# One step = one frame per case at a common simulation time (the first case drives the pairing)
case_frames = [gather_frames(case) for case, _ in CASES]
first_times = np.array([f[0] for f in case_frames[0]])
time_match_tol = TIME_MATCH_FRAC * (np.median(np.diff(first_times)) if len(first_times) > 1 else np.inf)
steps = []
for first in case_frames[0]:
    row = [first]
    for frames in case_frames[1:]:
        times = np.array([f[0] for f in frames])
        j = int(np.argmin(np.abs(times - first[0])))
        if abs(times[j] - first[0]) > time_match_tol:
            break
        row.append(frames[j])
    if len(row) == NROWS:
        steps.append(row)
if not steps:
    raise SystemExit("No common frame times between " + ", ".join(c for c, _ in CASES))
if COMPARE:
    print(f"{len(steps)} common frames, t = {steps[0][0][0]:.2f} .. {steps[-1][0][0]:.2f}")

if SNAPSHOT_TIME is not None:
    # START_TIME only matters for trimming the animation loop, so it is skipped here.
    steps = [min(steps, key=lambda row: abs(row[0][0] - SNAPSHOT_TIME))]
elif START_TIME > 0.0:
    steps = [row for row in steps if row[0][0] >= START_TIME]

loaded0 = [load_frame(plt_path, vtp_path) for _, plt_path, vtp_path in steps[0]]
t0 = loaded0[0][0]

if VMIN_OVERRIDE is not None and VMAX_OVERRIDE is not None:
    VMIN, VMAX = VMIN_OVERRIDE, VMAX_OVERRIDE
elif INFO["fixed_range"] is not None:
    VMIN, VMAX = INFO["fixed_range"]
elif SNAPSHOT_TIME is not None:
    # Snapshot mode: range the colorbar off the displayed frames
    # (already masked by load_frame) rather than scanning the whole series.
    VMIN = min(float(np.nanmin(d)) for _, d, _, _ in loaded0)
    VMAX = max(float(np.nanmax(d)) for _, d, _, _ in loaded0)
    if VMIN_OVERRIDE is not None:
        VMIN = VMIN_OVERRIDE
    if VMAX_OVERRIDE is not None:
        VMAX = VMAX_OVERRIDE
else:
    VMIN, VMAX = field_range(steps, FIELD, INFO["mask"])
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
STACK_HEIGHT_IN = NROWS * PANEL_HEIGHT_IN + (NROWS - 1) * ROW_GAP_IN
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + STACK_HEIGHT_IN + TOP_MARGIN_IN
CONTENT_WIDTH_IN = PLOT_WIDTH_IN + CBAR_GAP_IN + CBAR_WIDTH_IN
LEFT_IN = 0.5 * (FIG_WIDTH_IN - CONTENT_WIDTH_IN)

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))
# Rows top to bottom
axes = [
    add_axes_in_inches(fig, LEFT_IN, BOTTOM_MARGIN_IN + (NROWS - 1 - r) * (PANEL_HEIGHT_IN + ROW_GAP_IN), PLOT_WIDTH_IN, PANEL_HEIGHT_IN)
    for r in range(NROWS)
]
cax = add_axes_in_inches(
    fig, LEFT_IN + PLOT_WIDTH_IN + CBAR_GAP_IN, BOTTOM_MARGIN_IN, CBAR_WIDTH_IN, STACK_HEIGHT_IN
)

images, lines = [], []
for r, (ax, (_, label), (_, data0, segs0, extent)) in enumerate(zip(axes, CASES, loaded0)):
    im = ax.imshow(data0, origin="lower", extent=extent, cmap=CMAP, vmin=VMIN, vmax=VMAX, interpolation="none")
    lc = LineCollection(segs0, colors=INTERFACE_COLOR, linewidths=INTERFACE_LINEWIDTH)
    ax.add_collection(lc)
    images.append(im)
    lines.append(lc)

    ax.set_xlim(*VIEW_XLIM)
    ax.set_ylim(*VIEW_YLIM)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(x_ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.set_yticks(y_ticks)
    ax.yaxis.set_major_formatter(FuncFormatter(tick_formatter))
    ax.tick_params(which="both", top=True, right=True, pad=TICK_PAD_PT)

    # Horizontal tick numbers and axis label on the bottom row only
    if r == NROWS - 1:
        ax.set_xlabel(r"$x$", labelpad=5.0)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    ax.set_ylabel(r"$y$", labelpad=-5.0)
    if COMPARE:
        ax.set_title(label, fontsize=ROW_TITLE_FONTSIZE, pad=5.0)

if COMPARE:
    # Paired frames share the simulation time: one time label above all rows
    time_label = fig.text((LEFT_IN + 0.5 * PLOT_WIDTH_IN) / FIG_WIDTH_IN, 1.0 - 0.12 / FIG_HEIGHT_IN, time_text(t0),
                          ha="center", va="top", fontsize=TIME_FONTSIZE)
else:
    time_label = axes[0].set_title(time_text(t0), fontsize=TITLE_FONTSIZE, pad=6.0)

cbar = fig.colorbar(images[0], cax=cax, orientation="vertical")
cbar.set_ticks(np.linspace(VMIN, VMAX, 5))
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(pad=3.5)
cbar.set_label(INFO["label"], rotation=90, labelpad=5.0, fontsize=CBAR_LABEL_FONTSIZE)


def update(i):
    loaded = [load_frame(plt_path, vtp_path) for _, plt_path, vtp_path in steps[i]]
    for im, lc, (_, data, segments, _) in zip(images, lines, loaded):
        im.set_data(data)
        lc.set_segments(segments)
    time_label.set_text(time_text(loaded[0][0]))
    return (*images, *lines, time_label)


if SNAPSHOT_TIME is not None:
    snapshot_output = f"{FIELD}_plic_{OUTPUT_TAG}_t{t0:.2f}.pdf"
    fig.savefig(snapshot_output)
    print(f"Saved snapshot: {snapshot_output}  (t={t0:.4f})")
else:
    # Video resolution: at least one pixel per finest cell across each panel
    ds0 = yt.load(str(steps[0][0][1]))
    dcell0 = (ds0.domain_right_edge - ds0.domain_left_edge).to_value("code_length")[0] / (
        int(ds0.domain_dimensions[0]) * 2**ds0.index.max_level
    )
    video_dpi = max(ANIMATION_DPI, int(np.ceil(x_span / dcell0 / PLOT_WIDTH_IN)))
    print(f"Video: {video_dpi} dpi, {int(round(PLOT_WIDTH_IN * video_dpi))} px per panel width for {int(round(x_span / dcell0))} finest cells")
    # Whole, even pixel dimensions: otherwise ffmpeg receives frames of a different size than declared (sheared video)
    fig.set_size_inches(2 * int(np.ceil(FIG_WIDTH_IN * video_dpi / 2)) / video_dpi,
                        2 * int(np.ceil(FIG_HEIGHT_IN * video_dpi / 2)) / video_dpi)
    anim = FuncAnimation(fig, update, frames=len(steps), blit=False)
    anim.save(OUTPUT, fps=FPS, dpi=video_dpi, writer="ffmpeg")
    print(f"Saved: {OUTPUT}  ({len(steps)} frames)")
