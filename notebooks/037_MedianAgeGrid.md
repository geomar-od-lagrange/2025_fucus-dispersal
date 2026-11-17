---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: py:percent,md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Median particle age: with and without stranding, per regime

A grid of median-age maps: column 1 **with stranding** (each sample
weighted by its un-beached fraction `exp(−A)`, the 024f rate-model member
`member`), column 2 **without stranding** (the unweighted 024 counts), one
row per regime of `rows_csv`. Two columns across the 180 mm page, so each
map is half width; at most two rows, so a figure fits on one page.

Reads only 026's exports — no parquet, no Dask:
`Exports/026/026_medianage_<regime>_<season>_<tag>.geojson`, `tag` the
member or `unweighted`. 026 writes them with the mapped hexes only (share
of the window's weight at or above the 5th percentile, water area > 0), so
the hexes in the file are the hexes drawn; the median age, the colour
classes and the backdrop are 026's. Runs after 026 at `member` and
unweighted for every regime of `rows_csv`.

One figure per run. The rows are a parameter rather than a fixed pair of
figures, so the split across pages (default: surface + Stokes drift and
bottom on one, bottom with bed friction on the other) lives in the job
script, and each figure gets its own executed notebook to cite. Considered
one run drawing every figure from a nested `"a,b;c"` list; rejected, it
adds a second separator to parse and couples the figures' printouts in one
notebook for no saving worth having (the backdrop takes seconds).

```python
import itertools
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm
from shapely.geometry import box
from shapely.ops import unary_union
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (BSH static grids, wet-cell coastline).
data_root = "../data"
# Read root of 026's exports; write root for the figure.
output_root = "../output"
# Regimes drawn as rows, top to bottom; at most two (one page).
rows_csv = "surface_stokes,bottom"
# 024f rate-model member of the "with stranding" column: the 026 weighting
# tag of its exports. The "without stranding" column reads `unweighted`.
member = "step_wc0p1_ts3_tcinf_tf1_tw0_b750"
# Release season: DJF, MAM, JJA, SON or ALL.
season = "ALL"
# Hex radius 026 ran at. Labels the figure file only: the export names carry
# no radius, so this must match the 026 run that wrote them.
hex_radius = 6000

# Map extent (degrees), one for every map of the pipeline: the Baltic with
# the Belt Sea and Kattegat, North Sea cropped (docs/visualisations.md).
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0
```

# Parse parameters

The plain-language labels live here, not in the parameters cell: the
figure names regime and weighting in words, never by run key or member tag
(docs/visualisations.md "Labels in the figure").

```python
SEASONS = ("DJF", "MAM", "JJA", "SON", "ALL")
if season not in SEASONS:
    raise ValueError(f"season {season!r} not one of {SEASONS}")
rows = [r.strip() for r in rows_csv.split(",") if r.strip()]
if not 1 <= len(rows) <= 2:
    raise ValueError(f"rows_csv {rows_csv!r}: 1 or 2 regimes per figure, got {len(rows)}")

REGIME_LABELS = {
    "surface": "surface",
    "surface_stokes": "surface + Stokes drift",
    "bottom": "bottom",
    "bottom_ub49": "bottom, bed friction",
}
# Columns left to right: label and 026 weighting tag.
COLUMNS = {"with stranding": member, "without stranding": "unweighted"}

output_root = Path(output_root)
borders_path = output_root / "SubbasinBorders" / "subbasin_borders.geojson"
data_root = Path(data_root)
export_dir = output_root / "Exports" / "026"

# Print-ready output: one full text-width (180 mm) figure at 300 dpi, so it
# lands in the manuscript at its final size (rationale in
# docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Subbasin border stroke, as on every map of the pipeline (docs/visualisations.md).
basin_border_lw = 0.52

# 026's median-age classes: the default colormap, discretised into eleven
# classes, finer below 100 d; no over class (a median of ages < 180 d cannot
# exceed it).
AGE_BOUNDS_DAYS = [0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 140, 180]
age_cmap = plt.get_cmap().resampled(len(AGE_BOUNDS_DAYS) - 1)
age_norm = BoundaryNorm(AGE_BOUNDS_DAYS, age_cmap.N)

figure_dir = output_root / "Figures" / "037"
figure_dir.mkdir(parents=True, exist_ok=True)
fig_path = figure_dir / f"MedianAgeGrid_{'_'.join(rows)}_{season}_r{hex_radius}m.png"
print(f"rows {rows}, season {season}, columns {COLUMNS}")
```

# Read 026's median-age exports

One GeoJSON per (regime, weighting): hex geometry (EPSG:4326), `weight`
(survival-weighted `surv` or `n_obs`), `percentage` and
`median_age_days`. A missing file raises: the figure would otherwise drop
a panel silently.

```python
gdfs = {}
for regime, (column, tag) in itertools.product(rows, COLUMNS.items()):
    path = export_dir / f"026_medianage_{regime}_{season}_{tag}.geojson"
    gdfs[regime, column] = gpd.read_file(path)
    print(f"read {path.name}: {len(gdfs[regime, column]):,} hexes")
```

# Map backdrop: coastline and subbasin borders

026's backdrop. Coastline is the model's own: the outline of the BSH
HBMnoku wet cells (`coastline.geojson`), merged fine-first as in 024a —
fine polygons, plus the coarse ones outside the fine-grid bbox — drawn minus
a 0.001° strip along that bbox, where the two grids' cells meet in open
water. Subbasin borders: `006`'s `subbasin_borders.geojson` (output_root),
continuous HELCOM level-2 borders ending on the model coastline
(docs/visualisations.md "Subbasin borders").

```python
domain_box = box(extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max)
h0_fine = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
).H0

wet = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
# Fine bbox padded by half a fine cell, so every fine cell lies inside it.
fine_bbox = box(
    float(h0_fine.lon.min()) - 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.min()) - 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
    float(h0_fine.lon.max()) + 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.max()) + 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
)
wet_fine_first = unary_union([
    unary_union(wet.loc[wet["grid"] == "fine", "geometry"].tolist()),
    unary_union(wet.loc[wet["grid"] == "coarse", "geometry"].tolist()).difference(fine_bbox),
])
coast = (
    gpd.GeoSeries([wet_fine_first], crs="EPSG:4326").boundary
    .difference(fine_bbox.exterior.buffer(1e-3))
    .clip(domain_box)
)

basin_borders = gpd.read_file(borders_path).geometry.clip(domain_box)
print(f"subbasin borders: {len(basin_borders)} subbasin pairs (006)")

# Width/height of the extent in visually equal degrees — 1° lon at the
# extent's mid-latitude against 1° lat.
domain_aspect = (
    (extent_lon_max - extent_lon_min)
    * np.cos(np.radians(0.5 * (extent_lat_min + extent_lat_max)))
) / (extent_lat_max - extent_lat_min)
```

# Layout helpers

026's, unchanged: a constrained-layout grid at full page width,
`fit_height` measuring the height at which every map fills its column
width (no blank band between rows), and one horizontal colorbar under the
grid.

```python
def map_grid(nrows, ncols):
    """A constrained-layout panel grid at full page width. The height is only
    a seed (the maps alone, no decorations); ``fit_height`` sets the final
    one once the figure is drawn."""
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, layout="constrained")
    fig.set_size_inches(FIGURE_WIDTH_IN, nrows * FIGURE_WIDTH_IN / ncols / domain_aspect)
    return fig, axes


def fit_height(fig, axes):
    """Set the figure height at which every map fills its column width, with
    no blank band between rows. Per row, the height the map takes at its
    cell's full width against the cell height constrained_layout gave it,
    summed over the rows and added to the figure height; repeated until the
    correction is below half a pixel."""
    for _ in range(5):
        fig.draw_without_rendering()
        width_in, height_in = fig.get_size_inches()
        deficit_in = {}
        for ax in axes.flat:
            cell, drawn = ax.get_position(original=True), ax.get_position()
            full_width_height_in = cell.width * drawn.height / drawn.width * height_in
            deficit_in[ax.get_subplotspec().rowspan.start] = (
                full_width_height_in - cell.height * height_in
            )
        correction_in = sum(deficit_in.values())
        if abs(correction_in) < 0.5 / FIGURE_DPI:
            return
        fig.set_size_inches(width_in, height_in + correction_in)
    raise RuntimeError(f"figure height did not settle, last correction {correction_in:.4f} in")


def shared_colorbar(fig, axes, norm, cmap, label, extend, **kwargs):
    """One horizontal colorbar under the whole grid. pad override: the
    default 0.05 is a fraction of the whole grid's height, a blank band
    between the last row and the bar. panchor=False keeps the maps'
    top anchor (``hex_map``); the default re-anchors them to the bottom."""
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), ax=axes, orientation="horizontal",
        label=label, extend=extend, pad=0.01, panchor=False, **kwargs,
    )


def hex_map(gdf, column, ax, norm, cmap):
    """Draw one panel: hex polygons through ``norm`` and ``cmap``, subbasin
    borders and coastline, on a plain (non-cartopy) lon/lat axis — all
    EPSG:4326, so nothing is reprojected out of register."""
    # Only hexes touching the extent: off-extent patches still count in
    # constrained_layout's bounding box and push the title off the page.
    gdf = gdf.cx[extent_lon_min:extent_lon_max, extent_lat_min:extent_lat_max]
    gdf.plot(
        ax=ax, column=column, norm=norm, cmap=cmap,
        # Default polygon stroke is a black hairline that dominates at
        # Baltic zoom; matching edge to fill makes the hexes read as a
        # continuous field.
        edgecolor="face", linewidth=0.4, zorder=1,
    )
    # Override: neutral mid-grey thin subbasin borders (as 029), so a border
    # reads neither as coastline (black) nor as part of the hex field.
    basin_borders.plot(ax=ax, color="0.5", linewidth=basin_border_lw, zorder=1.5)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent_lon_min, extent_lon_max)
    ax.set_ylim(extent_lat_min, extent_lat_max)
    # Top anchor: spare cell height falls under the map while fit_height
    # measures, so the title above it keeps its place.
    ax.set_aspect(
        1 / np.cos(np.radians(0.5 * (extent_lat_min + extent_lat_max))), anchor="N"
    )
    ax.set_xticks([])
    ax.set_yticks([])
```

# Validation numbers

Per panel: hexes mapped (in the file, and inside the map extent), the
share of the mapped weight outside the extent, and the weight-weighted
median of the hex median ages — the hex age the median unit of
particle-time sits in, 026's printed statistic, recomputed from the export.

```python
def weighted_median(values, weights):
    """Value at which the weights, sorted by value, reach half their total."""
    order = np.argsort(values.to_numpy())
    cum = np.cumsum(weights.to_numpy()[order]) / weights.sum()
    return values.to_numpy()[order][int(np.searchsorted(cum, 0.5))]


summary = []
for (regime, column), gdf in gdfs.items():
    point = gdf.geometry.representative_point()
    inside = point.x.between(extent_lon_min, extent_lon_max) & point.y.between(
        extent_lat_min, extent_lat_max
    )
    summary.append({
        "regime": regime,
        "weighting": column,
        "hexes_mapped": len(gdf),
        "hexes_in_extent": int(inside.sum()),
        "weight_outside_extent_pct": 100 * gdf.loc[~inside, "weight"].sum() / gdf["weight"].sum(),
        "weighted_median_age_d": weighted_median(gdf["median_age_days"], gdf["weight"]),
        "unweighted_median_age_d": gdf["median_age_days"].median(),
    })
summary = pd.DataFrame(summary).set_index(["regime", "weighting"])
print(f"season {season}")
print(summary.round(2).to_string())
```

# Figure

Columns "with stranding" | "without stranding", one row per regime, the
regime written as the row's y label. One colorbar under the grid, every
class edge ticked (the classes are unequal).

```python
fig, axes = map_grid(len(rows), len(COLUMNS))
for (i, regime), (j, column) in itertools.product(enumerate(rows), enumerate(COLUMNS)):
    ax = axes[i, j]
    hex_map(gdfs[regime, column], "median_age_days", ax, age_norm, age_cmap)
    if i == 0:
        ax.set_title(column)
    if j == 0:
        ax.set_ylabel(REGIME_LABELS.get(regime, regime))
shared_colorbar(
    fig, axes, age_norm, age_cmap, "median particle age (days)", "neither",
    ticks=AGE_BOUNDS_DAYS,
)
fit_height(fig, axes)
# Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path} ({fig.get_size_inches()[1] * 25.4:.1f} mm tall)")
plt.show()
```
