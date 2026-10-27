# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,md,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Seeding map: release cells by subbasin, and their land-seeded share
#
# The release set as the dispersal runs see it:
# the HELCOM REDLIST *Fucus vesiculosus* cells `010` seeds from
# (`helcom_fucus_redlist/fucus_release_points.geojson`, baked by `000`), in
# two panels on one figure:
#
# 1. **Release cells by HELCOM subbasin** — each cell filled by the level-2
#    subbasin of its hex under the production assignment (largest overlap
#    plus shoreward snap, the `024a` key's `helcom_subbasin`), the borders
#    between the subbasin polygons drawn.
# 2. **Land-seeded share per cell** — the share of the cell's seeded particles
#    that start on BSH land and never move (the sentinel `release_hex = -1`
#    rows of the `024` counts store), from the trajectory zarrs.
#
# The land-seeded share panel is also saved on its own, one map at the full
# page width, as the paper's supplement figure.
#
# Plus a computed table of release cells (and their land-seeded share) per
# subbasin. Reads the data twin, the `024a` key and the first two time steps
# of one regime-year of trajectory zarrs; no Dask cluster.

# %%
import itertools
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, to_hex
from matplotlib.patches import Patch
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

from hextraj import HexProj

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex key and the trajectory zarrs; write root for the figure.
output_root = "../output"
# Read root of the data twin: the release cells and the HELCOM subbasins.
data_root = "../data"
# Hex radius of the 024a key whose subbasin assignment colours the cells.
hex_radius = 6000

# Trajectory zarrs the land-seeded share is counted from: every release of
# one regime and release year. The BSH land mask is the same for every regime
# and release, so one regime-year (73 releases × particles_per_cell per cell)
# is the sample; more years only shrink the sampling noise of a fixed share.
land_regime = "surface_stokes"
land_release_year = 2019
# Particles 010 seeds per release cell, cell-major (010's particles_per_cell).
particles_per_cell = 100

# Map extent (degrees), as 029/033: the Baltic proper, North Sea cropped.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
data_root = Path(data_root)
release_cells_path = data_root / "helcom_fucus_redlist" / "fucus_release_points.geojson"
subbasins_path = data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp"
zarr_dir = output_root / "Trajectories" / land_regime / str(land_release_year)

# Print-ready output: one full text-width (180 mm) figure at 300 dpi, so the
# panels land in the manuscript at their final size (docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Vertical inches under the two maps for the subbasin legend (17 entries in
# three columns) and the panel titles.
LEGEND_IN = 1.45
# Vertical inches the full-width single land-share map adds to its own
# height: the title above, the horizontal colorbar with ticks and label below.
SINGLE_DECOR_IN = 0.75

figure_dir = output_root / "Figures" / "036"
figure_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Release cells and their subbasin
#
# A cell is assigned to the subbasin of the `024a` hex at its representative
# point, as 029 does for its overlay: the cell (~10 km) is larger than a
# 6 km hex, but its particles share one subbasin in all but a few border
# cells, and the key's assignment is the one every subbasin product
# (026a, 028, 029's afloat table) uses for the particles' `release_hex`. The
# HELCOM polygons are a marine delineation that stops short of the shore, so
# a coastal cell often lies outside every polygon on the map and still gets a
# subbasin through the shoreward snap
# ([docs/hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md));
# `_outside` (-1) is a cell whose hex no rule reaches.

# %%
key_path = output_root / "HexAggregates" / f"HexAgg_key_r{hex_radius}m.parquet"
key = pd.read_parquet(key_path, columns=["hex_id", "helcom_subbasin"])
meta = json.loads(key_path.with_suffix(".json").read_text())
hp = HexProj(**meta["hex_proj"])
subbasin_name = {int(k): v for k, v in meta["subbasin_id_to_name"].items()}
print(f"key: {key_path.name}, subbasin rule {meta['subbasin_rule']!r}, "
      f"snap {meta['subbasin_snap_m']:,.0f} m")

cells = gpd.read_file(release_cells_path)
point = cells.representative_point()
cells["release_hex"] = hp.label(point.x.to_numpy(), point.y.to_numpy()).astype("int64")
cells["subbasin"] = (
    cells["release_hex"]
    .map(key.set_index("hex_id")["helcom_subbasin"])
    .fillna(-1)
    .astype(int)
    .map(subbasin_name)
)
print(f"release cells: {len(cells):,} from {release_cells_path}")

# Subbasin borders as the edges two HELCOM polygons share. The polygons'
# full outlines also trace every island and skerry of the coast they were
# digitised against, which in the Åland and Stockholm archipelagos is a
# solid smear over the very cells the map shows; the shared edges are the
# basin-to-basin borders alone. The one border this drops is the Kattegat's
# open edge to the Skagerrak, which no other polygon shares.
subbasin_outline = gpd.read_file(subbasins_path).to_crs(4326).boundary
basin_borders = gpd.GeoSeries(
    [
        a.intersection(b)
        for a, b in itertools.combinations(subbasin_outline, 2)
        if a.intersects(b)
    ],
    crs="EPSG:4326",
)
print(f"subbasin borders: {len(basin_borders)} shared edges between "
      f"{len(subbasin_outline)} HELCOM polygons")

# %% [markdown]
# # Land-seeded share per cell
#
# `010` seeds `particles_per_cell` particles uniformly inside each cell,
# **cell-major**: trajectory `i` belongs to cell `i // particles_per_cell` in
# the release-cell file's row order. A particle seeded on BSH land never
# moves, and the pipeline's land-seeded test is exactly that — zero lon and lat
# displacement over the first output step (`024`'s `release_hex = -1`). Read
# from the first two time steps only, so a regime-year of zarrs costs seconds.

# %%
zarrs = sorted(zarr_dir.glob("*.zarr"))
if not zarrs:
    raise FileNotFoundError(f"no trajectory zarrs under {zarr_dir}")
n_land = np.zeros(len(cells))
n_seeded = np.zeros(len(cells))
for path in zarrs:
    ds = xr.open_zarr(path)[["lon", "lat"]].isel(obs=slice(0, 2)).compute()
    assert ds.sizes["trajectory"] == particles_per_cell * len(cells), (
        f"{path.name}: {ds.sizes['trajectory']} trajectories, expected "
        f"{particles_per_cell} × {len(cells)} cells"
    )
    on_land = ((ds.lon.diff("obs") == 0) & (ds.lat.diff("obs") == 0)).isel(obs=0)
    cell = (ds.trajectory // particles_per_cell).to_numpy()
    n_land += np.bincount(cell, weights=on_land.to_numpy(), minlength=len(cells))
    n_seeded += np.bincount(cell, minlength=len(cells))
cells["land_share_pct"] = 100 * n_land / n_seeded
print(f"{len(zarrs)} zarrs of {land_regime} {land_release_year}: "
      f"{int(n_seeded.sum()):,} seeded particles, {int(n_land.sum()):,} on land "
      f"({100 * n_land.sum() / n_seeded.sum():.1f} %)")
print(f"cells with no land-seeded particle: {int((cells['land_share_pct'] == 0).sum()):,}; "
      f"≥ 50 %: {int((cells['land_share_pct'] >= 50).sum()):,}; "
      f"all on land: {int((cells['land_share_pct'] == 100).sum()):,}")

# %% [markdown]
# # Release cells per subbasin

# %%
table = (
    cells.groupby("subbasin")
    .agg(cells=("subbasin", "size"), land_share_pct=("land_share_pct", "mean"))
    .sort_values("cells", ascending=False)
)
table.insert(1, "cells_pct", 100 * table["cells"] / table["cells"].sum())
table.loc["all"] = [
    len(cells), 100.0, cells["land_share_pct"].mean()
]
table["cells"] = table["cells"].astype(int)
print(table.round(1).to_string())

# %% [markdown]
# # Map
#
# Left: cells by subbasin. **Categorical colours, one override defended in
# [docs/visualisations.md](../docs/visualisations.md):** 16 subbasins hold
# release cells, and the default cycle has 10 colours, so it would give two
# subbasins one colour and the legend could not tell them apart. The palette is
# the default cycle followed by its lighter tints from `tab20` (the cycle's
# paired extension), grey held back for `_outside`. Subbasins take colours in
# order of their cells' mean latitude, so a dark hue and its light tint go to
# a southern and a northern subbasin and never to neighbours.
#
# Right: the land-seeded share per cell on the default colormap, fixed 0–100 %.
# Both panels draw the subbasin borders and the Natural Earth coastline.

# %%
lon_min, lon_max, lat_min, lat_max = (
    extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max
)
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)

cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# tab20's odd entries are the light tints of the default (tab10) cycle, in
# cycle order; index 7 of both is grey, reserved for `_outside`.
GREY = 7
palette = [c for i, c in enumerate(cycle) if i != GREY] + [
    to_hex(c) for i, c in enumerate(plt.colormaps["tab20"].colors[1::2]) if i != GREY
]
order = (
    cells[cells["subbasin"] != subbasin_name[-1]]
    .assign(lat=point.y)
    .groupby("subbasin")["lat"].mean()
    .sort_values()
    .index
)
assert len(order) <= len(palette), f"{len(order)} subbasins, {len(palette)} colours"
colour = dict(zip(order, palette)) | {subbasin_name[-1]: cycle[GREY]}


def backdrop(ax):
    """Subbasin borders, coastline, extent and aspect on a plain lon/lat
    axis (everything is EPSG:4326, so nothing is reprojected)."""
    # Override: neutral mid-grey hairline borders, so a subbasin border reads
    # neither as coastline (black) nor as one of the category fills.
    basin_borders.plot(ax=ax, color="0.5", linewidth=0.4, zorder=2)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=3)
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (lat_min + lat_max))))
    ax.set_xticks([])
    ax.set_yticks([])


fig, (ax_sub, ax_land) = plt.subplots(ncols=2, layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, 0.8 * FIGURE_WIDTH_IN / 2 / domain_aspect + LEGEND_IN)

cells.plot(ax=ax_sub, color=cells["subbasin"].map(colour), zorder=1)
backdrop(ax_sub)
ax_sub.set_title(f"{len(cells)} release cells by subbasin")

land_norm = Normalize(vmin=0, vmax=100)
# Colorbar right of the map at its full height, as 027/029/033's tall
# Baltic-cropped maps.
cax = ax_land.inset_axes([1.02, 0.0, 0.035, 1.0])
cells.plot(
    ax=ax_land, column="land_share_pct", norm=land_norm, legend=True, cax=cax,
    legend_kwds={"label": "land-seeded share of the cell's particles (%)"}, zorder=1,
)
backdrop(ax_land)
ax_land.set_title("land-seeded share per cell")

counts = cells["subbasin"].value_counts()
fig.legend(
    handles=[
        Patch(facecolor=colour[name], label=f"{name} ({counts[name]})")
        for name in [*order, subbasin_name[-1]] if name in counts
    ],
    loc="outside lower center", ncols=3,
)
fig_path = figure_dir / f"SeedingMap_r{hex_radius}m.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Land-seeded share, single panel
#
# The right-hand panel above on its own for the supplement: one map at the
# full page width, its colorbar under it at the same width (029's
# where-stranded layout), the same fixed 0–100 % default-colormap scale.

# %%
fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, FIGURE_WIDTH_IN / domain_aspect + SINGLE_DECOR_IN)
# Colorbar under the map at its full width, as 029's where-stranded map; the
# map is bottom-anchored so constrained_layout's room for the bar holds
# (centred, the map drifts down and clips the colorbar label).
cax = ax.inset_axes([0.0, -0.045, 1.0, 0.022])
ax.set_anchor("S")
cells.plot(
    ax=ax, column="land_share_pct", norm=land_norm, legend=True, cax=cax,
    legend_kwds={"label": "land-seeded share of the cell's particles (%)",
                 "orientation": "horizontal"},
    zorder=1,
)
backdrop(ax)
ax.set_title("land-seeded share per release cell")
fig_path = figure_dir / f"LandSeededShare_r{hex_radius}m.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Validation prints

# %%
print(f"release cells: {len(cells):,}; subbasins holding cells: {len(order)}; "
      f"_outside cells: {int(counts.get(subbasin_name[-1], 0))}")
print(f"land-seeded share: pooled {100 * n_land.sum() / n_seeded.sum():.2f} %, "
      f"cell mean {cells['land_share_pct'].mean():.2f} %, "
      f"cell median {cells['land_share_pct'].median():.1f} %")
print(f"land sample: {land_regime} {land_release_year}, {len(zarrs)} releases, "
      f"{int(n_seeded.min()):,}–{int(n_seeded.max()):,} particles per cell")
