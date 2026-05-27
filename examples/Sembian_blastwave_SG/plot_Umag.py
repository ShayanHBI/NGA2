#!/usr/bin/env pvpython

from pathlib import Path
import re
import numpy as np

from paraview.simple import *
from paraview import servermanager

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter, FuncFormatter


# ============================================================
# User settings
# ============================================================

DATASETS = [
    (
        r"\textbf{(a)} $pT$",
        "/Users/shayanhbi/Repositories/NGA2/examples/Sembian_blastwave_NASG/amrviz/Sembian_blastwave_NASG_PT",
    ),
    (
        r"\textbf{(b)} $pTg$",
        "/Users/shayanhbi/Repositories/NGA2/examples/Sembian_blastwave_NASG/amrviz/Sembian_blastwave_NASG_PTg",
    ),
]

FIELD = "Umag"

TARGET_TIME = 4.0007e-4
READER_LEVEL = 100
PLOT_NUMBER = None

X_BOUNDS = (0.14, 0.25)
Y_BOUNDS = (-0.037, 0.037)

VMIN = 0.0
VMAX = 650.0
CBAR_TICKS = [0.0, 162.5, 325.0, 487.5, 650.0]
CMAP = "turbo"

RASTERIZED = True
SAVE_DPI = 600
OUTPUT = "uMag.pdf"


# ============================================================
# Paper layout
# ============================================================
#
# This script preserves the true data aspect ratio:
#
#   physical panel height / physical panel width = Δy / Δx
#
# The figure is tuned to look readable when inserted around 0.75--0.85 textwidth.
# In LaTeX, use for example:
#
#   \includegraphics[width=0.78\textwidth]{uMag}
#
# Do NOT use 0.5\textwidth for this vertical two-panel figure unless you want
# everything to become tiny.
#
# ============================================================

FIG_WIDTH_IN = 5.05

LEFT_MARGIN_IN = 0.62
RIGHT_MARGIN_IN = 0.22
BOTTOM_MARGIN_IN = 0.58
TOP_MARGIN_IN = 0.18

PLOT_WIDTH_IN = 3.72
PANEL_GAP_IN = 0.30

CBAR_GAP_IN = 0.18
CBAR_WIDTH_IN = 0.20


# ============================================================
# Matplotlib / LaTeX style
# ============================================================

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",

    # These are intentionally larger than the previous version because the
    # figure is inserted smaller than \textwidth in the AIAA paper.
    "font.size": 11.0,
    "axes.labelsize": 12.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,

    "axes.linewidth": 1.05,
    "xtick.major.width": 1.05,
    "ytick.major.width": 1.05,
    "xtick.minor.width": 0.85,
    "ytick.minor.width": 0.85,
    "xtick.major.size": 5.0,
    "ytick.major.size": 5.0,
    "xtick.minor.size": 3.0,
    "ytick.minor.size": 3.0,
    "xtick.direction": "in",
    "ytick.direction": "in",

    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}",
})


# ============================================================
# Helper functions
# ============================================================

def plotfile_number(path: Path):
    m = re.search(r"(\d+)$", path.name)
    if m is None:
        return None
    return int(m.group(1))


def find_plotfiles(root):
    root = Path(root)
    plotfiles = sorted(
        [p for p in root.glob("plt.nga2.cell.*") if p.is_dir()],
        key=lambda p: plotfile_number(p) if plotfile_number(p) is not None else -1,
    )

    if not plotfiles:
        raise RuntimeError(f"No plotfiles found in:\n{root}")

    return plotfiles


def read_amrex_header(plotdir: Path):
    header = Path(plotdir) / "Header"
    if not header.exists():
        raise RuntimeError(f"Header file not found: {header}")

    lines = header.read_text(errors="ignore").splitlines()
    if len(lines) < 5:
        raise RuntimeError(f"Header file appears too short: {header}")

    ncomp = int(lines[1].strip())
    names = [lines[2 + i].strip() for i in range(ncomp)]

    ndim_line = 2 + ncomp
    time_line = ndim_line + 1

    try:
        ndim = int(lines[ndim_line].strip())
    except Exception:
        ndim = None

    try:
        time_value = float(lines[time_line].strip())
    except Exception as exc:
        raise RuntimeError(
            f"Could not parse time from {header}. Expected it on line {time_line + 1}."
        ) from exc

    return {
        "header": header,
        "ncomp": ncomp,
        "names": names,
        "ndim": ndim,
        "time": time_value,
    }


def choose_plotfiles(dataset_roots, requested_number=None, target_time=None):
    chosen = []

    for root in dataset_roots:
        plotfiles = find_plotfiles(root)

        if requested_number is not None:
            num_to_path = {plotfile_number(p): p for p in plotfiles}
            if requested_number not in num_to_path:
                raise RuntimeError(
                    f"Requested plot number {requested_number} was not found in {root}."
                )
            chosen.append(num_to_path[requested_number])
            continue

        if target_time is None:
            chosen.append(plotfiles[-1])
            continue

        best_plotfile = None
        best_error = None

        for p in plotfiles:
            try:
                meta = read_amrex_header(p)
                t = meta["time"]
            except Exception:
                continue

            err = abs(t - target_time)
            if best_error is None or err < best_error:
                best_plotfile = p
                best_error = err

        if best_plotfile is None:
            raise RuntimeError(
                f"Could not find any readable AMReX plotfile headers in {root}."
            )

        chosen.append(best_plotfile)

    return chosen


def set_reader_level(reader, requested_level):
    if requested_level is None:
        return

    try:
        reader.UpdatePipelineInformation()
    except Exception:
        pass

    try:
        reader.UpdatePropertyInformation()
    except Exception:
        pass

    preferred_names = [
        "Level",
        "MaxLevel",
        "MaximumLevel",
        "MaximumLevels",
        "MaxLevels",
        "RefinementLevel",
        "MaximumRefinementLevel",
    ]

    candidate_names = []
    props = list(reader.ListProperties())

    for name in preferred_names:
        if name in props:
            candidate_names.append(name)

    for name in props:
        if "level" in name.lower() and name not in candidate_names:
            candidate_names.append(name)

    for prop_name in candidate_names:
        for value in (requested_level, [requested_level]):
            try:
                setattr(reader, prop_name, value)
                return
            except Exception:
                pass


def enable_reader_arrays(reader, array_names):
    if not array_names:
        return

    try:
        reader.UpdatePipelineInformation()
    except Exception:
        pass

    try:
        reader.UpdatePropertyInformation()
    except Exception:
        pass

    candidate_properties = [
        "CellArrayStatus",
        "PointArrayStatus",
        "CellArrays",
        "PointArrays",
    ]

    for prop_name in candidate_properties:
        try:
            reader.GetProperty(prop_name)
            setattr(reader, prop_name, array_names)
        except Exception:
            pass

    try:
        UpdatePipeline(proxy=reader)
    except Exception:
        pass


def iter_leaf_datasets(data_object):
    if data_object is None:
        return

    try:
        is_composite = data_object.IsA("vtkCompositeDataSet")
    except Exception:
        is_composite = False

    if is_composite:
        iterator = data_object.NewIterator()
        iterator.InitTraversal()
        while not iterator.IsDoneWithTraversal():
            child = iterator.GetCurrentDataObject()
            if child is not None:
                yield from iter_leaf_datasets(child)
            iterator.GoToNextItem()
    else:
        yield data_object


def open_amrex_plotfile(plotdir: Path):
    plotdir = Path(plotdir)
    meta = read_amrex_header(plotdir)

    last_error = None
    candidates = [plotdir, meta["header"]]

    for candidate in candidates:
        try:
            src = OpenDataFile(str(candidate))
            set_reader_level(src, READER_LEVEL)
            enable_reader_arrays(src, meta["names"])
            UpdatePipeline(proxy=src)
            return src, meta
        except Exception as e:
            last_error = e

    raise RuntimeError(
        f"Could not open AMReX plotfile:\n{plotdir}\n\nLast error was:\n{last_error}"
    )


def make_surface_proxy_without_interpolation(src):
    filter_attempts = []

    try:
        filter_attempts.append(("ExtractSurface", lambda: ExtractSurface(Input=src)))
    except NameError:
        pass

    try:
        filter_attempts.append(("DataSetSurfaceFilter", lambda: DataSetSurfaceFilter(Input=src)))
    except NameError:
        pass

    last_error = None

    for _, maker in filter_attempts:
        try:
            surf = maker()
            for prop_name in ["PassThroughCellIds", "PassThroughPointIds"]:
                try:
                    setattr(surf, prop_name, 1)
                except Exception:
                    pass
            UpdatePipeline(proxy=surf)
            return surf
        except Exception as e:
            last_error = e

    raise RuntimeError(
        "Could not create a fetchable surface proxy from the AMR dataset. "
        f"Last error: {last_error}"
    )


def cell_intersects_window(bounds, x_bounds, y_bounds):
    x0, x1, y0, y1, _, _ = bounds

    if x_bounds is not None:
        xmin, xmax = x_bounds
        if x1 < xmin or x0 > xmax:
            return False

    if y_bounds is not None:
        ymin, ymax = y_bounds
        if y1 < ymin or y0 > ymax:
            return False

    return True


def extract_cell_rectangles_from_fetched(fetched, field, x_bounds=None, y_bounds=None):
    rectangles = []
    values = []
    areas = []

    for ds in iter_leaf_datasets(fetched):
        if not hasattr(ds, "GetNumberOfCells"):
            continue

        ncell = ds.GetNumberOfCells()
        if ncell == 0:
            continue

        cd = ds.GetCellData()
        if cd is None:
            continue

        arr = cd.GetArray(field)
        if arr is None:
            continue

        ghost = cd.GetArray("vtkGhostType")

        for icell in range(ncell):
            if ghost is not None:
                try:
                    if int(ghost.GetTuple1(icell)) != 0:
                        continue
                except Exception:
                    pass

            cell = ds.GetCell(icell)
            b = [0.0] * 6
            cell.GetBounds(b)
            x0, x1, y0, y1, _, _ = b

            if not cell_intersects_window(b, x_bounds, y_bounds):
                continue

            try:
                val = arr.GetTuple1(icell)
            except Exception:
                val = arr.GetTuple(icell)[0]

            if not np.isfinite(val):
                continue

            area = (x1 - x0) * (y1 - y0)
            if area <= 0.0:
                continue

            rectangles.append([
                (x0, y0),
                (x1, y0),
                (x1, y1),
                (x0, y1),
            ])
            values.append(val)
            areas.append(area)

    if not rectangles:
        raise RuntimeError(
            f"No cell-centered values for field '{field}' were extracted."
        )

    rectangles = np.asarray(rectangles, dtype=float)
    values = np.asarray(values, dtype=float)
    areas = np.asarray(areas, dtype=float)

    # Draw coarse cells first and fine cells last.
    order = np.argsort(areas)[::-1]
    rectangles = rectangles[order]
    values = values[order]

    return rectangles, values


def sample_plotfile_cell_data(plotfile, field, x_bounds=None, y_bounds=None):
    src, meta = open_amrex_plotfile(plotfile)
    surf = make_surface_proxy_without_interpolation(src)
    fetched = servermanager.Fetch(surf)
    rectangles, values = extract_cell_rectangles_from_fetched(
        fetched,
        field=field,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
    )
    return rectangles, values, meta["time"]


def cbar_tick_formatter(value, _):
    if abs(value - round(value)) < 1.0e-8:
        return rf"${int(round(value))}$"
    return rf"${value:.1f}$"


# ============================================================
# Main
# ============================================================

labels = [item[0] for item in DATASETS]
roots = [item[1] for item in DATASETS]

chosen_plotfiles = choose_plotfiles(
    roots,
    requested_number=PLOT_NUMBER,
    target_time=TARGET_TIME,
)

plot_data = []
for label, plotfile in zip(labels, chosen_plotfiles):
    rectangles, values, actual_time = sample_plotfile_cell_data(
        plotfile=plotfile,
        field=FIELD,
        x_bounds=X_BOUNDS,
        y_bounds=Y_BOUNDS,
    )
    plot_data.append((label, plotfile, rectangles, values, actual_time))

norm = Normalize(vmin=VMIN, vmax=VMAX)

# Compute physical figure geometry from the data aspect ratio.
# This is the part that prevents stretching.
x_span = X_BOUNDS[1] - X_BOUNDS[0]
y_span = Y_BOUNDS[1] - Y_BOUNDS[0]
data_height_over_width = y_span / x_span

PANEL_HEIGHT_IN = PLOT_WIDTH_IN * data_height_over_width
CBAR_HEIGHT_IN = 2.0 * PANEL_HEIGHT_IN + PANEL_GAP_IN
FIG_HEIGHT_IN = BOTTOM_MARGIN_IN + CBAR_HEIGHT_IN + TOP_MARGIN_IN

min_width = (
    LEFT_MARGIN_IN
    + PLOT_WIDTH_IN
    + CBAR_GAP_IN
    + CBAR_WIDTH_IN
    + RIGHT_MARGIN_IN
)
if FIG_WIDTH_IN < min_width:
    FIG_WIDTH_IN = min_width

fig = plt.figure(figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN))

def add_axes_in_inches(fig, left, bottom, width, height):
    fig_w, fig_h = fig.get_size_inches()
    return fig.add_axes([
        left / fig_w,
        bottom / fig_h,
        width / fig_w,
        height / fig_h,
    ])

bottom_lower = BOTTOM_MARGIN_IN
bottom_upper = BOTTOM_MARGIN_IN + PANEL_HEIGHT_IN + PANEL_GAP_IN

axes = [
    add_axes_in_inches(fig, LEFT_MARGIN_IN, bottom_upper, PLOT_WIDTH_IN, PANEL_HEIGHT_IN),
    add_axes_in_inches(fig, LEFT_MARGIN_IN, bottom_lower, PLOT_WIDTH_IN, PANEL_HEIGHT_IN),
]

cax = add_axes_in_inches(
    fig,
    LEFT_MARGIN_IN + PLOT_WIDTH_IN + CBAR_GAP_IN,
    bottom_lower,
    CBAR_WIDTH_IN,
    CBAR_HEIGHT_IN,
)

last_collection = None

for i, (ax, (label, plotfile, rectangles, values, actual_time)) in enumerate(zip(axes, plot_data)):
    collection = PolyCollection(
        rectangles,
        array=values,
        cmap=CMAP,
        norm=norm,
        edgecolors="none",
        linewidths=0.0,
        antialiased=False,
        rasterized=RASTERIZED,
    )

    ax.add_collection(collection)
    last_collection = collection

    ax.set_xlim(*X_BOUNDS)
    ax.set_ylim(*Y_BOUNDS)

    # Preserve true data aspect ratio: no stretching.
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(np.arange(0.14, 0.251, 0.02))
    ax.set_yticks([-0.03, 0.00, 0.03])
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(3))

    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(which="major", length=5.0, width=1.05)
    ax.tick_params(which="minor", length=3.0, width=0.85)

    for spine in ax.spines.values():
        spine.set_linewidth(1.05)

    ax.text(
        0.025,
        0.88,
        label,
        transform=ax.transAxes,
        fontsize=13.0,
        ha="left",
        va="center",
    )

    ax.set_ylabel(r"$y\;(\mathrm{m})$", labelpad=3.0)

    if i == 0:
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(r"$x\;(\mathrm{m})$", labelpad=3.0)

cbar = fig.colorbar(last_collection, cax=cax, orientation="vertical")
cbar.set_ticks(CBAR_TICKS)
cbar.ax.yaxis.set_major_formatter(FuncFormatter(cbar_tick_formatter))
cbar.ax.tick_params(direction="in", length=5.0, width=1.05, labelsize=10.5)
cbar.outline.set_linewidth(1.05)
cbar.set_label(
    r"$|\boldsymbol{u}|\;\left(\frac{\mathrm{m}}{\mathrm{s}}\right)$",
    rotation=90,
    labelpad=8.0,
    fontsize=12.5,
)

fig.savefig(OUTPUT, dpi=SAVE_DPI)
print(f"Saved figure to: {OUTPUT}")
