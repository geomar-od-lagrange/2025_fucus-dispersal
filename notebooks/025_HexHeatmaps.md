---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: min_data (pixi)
    language: python
    name: min_data
---

# Hex heatmaps

Particle-density maps on the hex-aggregated dispersal stores built by
notebook 024. Reads `key.parquet` and `counts/` partitions from the two
pre-built stores (`baltic_r6km_v1` and `de_r4km_v1`); no trajectory
zarrs, no Dask cluster.

Four panels per run:

- **A** — whole-Baltic density (baltic store, all release hexes summed)
- **B** — per-HELCOM-subbasin density (baltic store, one sub-plot per
  release subbasin)
- **C** — German-waters zoom (DE store, all release hexes summed)
- **D** — per-release-quarter density (baltic store, 2×2 grid)

```python
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib.pyplot as plt
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

from helpers import QUARTER_LABELS
```

# Parameters

```python tags=["parameters"]
data_root = "../data"
output_root = "../output"
regime = "surface"
release_year = 2019
baltic_config_name = "baltic_r6km_v1"
de_config_name = "de_r4km_v1"

# Baltic extent
lon_min = 5
lon_max = 32
lat_min = 53
lat_max = 66

# DE zoom extent
de_lon_min = 8
de_lon_max = 15
de_lat_min = 53.2
de_lat_max = 55.5

# Colormap (kept as a named override since matplotlib default colormap
# isn't well-suited to log-density; see plotting convention note below)
cmap = "viridis"

# Panel geometry — kept from current 024 because the LAEA-vs-lonlat
# aspect mismatch would otherwise pad the figure with whitespace.
baltic_panel_height_in = 6
de_panel_height_in = 4.5
```

# Path construction and parquet reads

```python
data_root = Path(data_root)
output_root = Path(output_root)
baltic_store = output_root / "HexAggregates" / baltic_config_name
de_store = output_root / "HexAggregates" / de_config_name

baltic_key = gpd.read_parquet(baltic_store / "key.parquet")
baltic_counts = pd.read_parquet(
    baltic_store / "counts" / f"regime={regime}" / f"release_year={release_year}" / "part.parquet"
)
de_key = gpd.read_parquet(de_store / "key.parquet")
de_counts = pd.read_parquet(
    de_store / "counts" / f"regime={regime}" / f"release_year={release_year}" / "part.parquet"
)

# Read subbasin_id_to_name from the baltic key.parquet file-level metadata.
_baltic_key_meta = json.loads(
    pq.read_table(baltic_store / "key.parquet").schema.metadata[b"hex_aggregate_store"].decode()
)
subbasin_id_to_name = {int(k): v for k, v in _baltic_key_meta["subbasin_id_to_name"].items()}

print(f"baltic_key: {len(baltic_key):,} hexes")
print(f"baltic_counts: {len(baltic_counts):,} rows, sum(n_obs)={int(baltic_counts.n_obs.sum()):,}")
print(f"de_key: {len(de_key):,} hexes")
print(f"de_counts: {len(de_counts):,} rows, sum(n_obs)={int(de_counts.n_obs.sum()):,}")
```

# Rendering helpers

`lonlat_aspect` and `log_density_plot` are copied from the zarr-based
predecessor (notebook renumbered from 024 → 025 during prod-prep). Coast
and subbasin data are loaded once here and reused across all panels.

```python
def lonlat_aspect(extent):
    """Displayed width / height ratio for a lon/lat extent with an
    aspect that keeps 1 deg lon at lat_mean visually equal to 1 deg lat."""
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    lat_mean = 0.5 * (lat_min_ + lat_max_)
    return ((lon_max_ - lon_min_) * np.cos(np.radians(lat_mean))) / (
        lat_max_ - lat_min_
    )


def log_density_plot(gdf, ax, extent, title=None, overlay=None, coast=None):
    """Plot hex density, coastline, and optional subbasin overlay on a
    plain (non-cartopy) lon/lat axis. Everything stays in EPSG:4326 so
    hexes, coastline, and overlay share one coordinate system — no
    drift from mismatched reprojection pipelines."""
    gdf = gdf.copy()
    gdf["log_count"] = np.log10(gdf["count"].where(gdf["count"] > 0))
    gdf.plot(
        ax=ax,
        column="log_count",
        cmap=cmap,
        legend=False,
        missing_kwds={"color": "none"},
        edgecolor="face",
        linewidth=0.4,
        zorder=1,
    )
    if coast is not None:
        coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    if overlay is not None:
        overlay.boundary.plot(
            ax=ax,
            color="magenta",
            linewidth=1.05,
            zorder=3,
        )
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    ax.set_xlim(lon_min_, lon_max_)
    ax.set_ylim(lat_min_, lat_max_)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (lat_min_ + lat_max_))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)
```

```python
# Load Natural Earth 10m coastline via cartopy's shapereader (cached
# locally), clip to each extent once so per-panel plotting is cheap.
_coast_gdf = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
)
coast_baltic = _coast_gdf.clip(box(lon_min, lat_min, lon_max, lat_max))
coast_de = _coast_gdf.clip(box(de_lon_min, de_lat_min, de_lon_max, de_lat_max))

# HELCOM subbasin polygons, smoothed at ~5 km scale so the overlay is
# readable on Baltic-wide panels. Buffer out then in (in lon/lat);
# 0.05 deg ≈ 5.5 km at mid-latitudes.
subbasins = gpd.read_file(
    data_root / "helcom_subbasins/HELCOM_subbasins_2022_level2.shp"
).to_crs(epsg=4326).rename(columns={"level_2": "subbasin"})
subbasins_smoothed = subbasins.assign(
    geometry=subbasins.geometry.buffer(0.05).buffer(-0.05)
)
```

# Panel A — Whole-Baltic density

```python
gdf_baltic = (
    baltic_counts.groupby("target_hex")["n_obs"].sum()
    .rename("count")
    .reset_index()
    .merge(baltic_key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
    .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
)

baltic_extent = [lon_min, lon_max, lat_min, lat_max]
baltic_aspect = lonlat_aspect(baltic_extent)
fig, ax = plt.subplots(
    figsize=(baltic_panel_height_in * baltic_aspect, baltic_panel_height_in),
    layout="constrained",
)
log_density_plot(gdf_baltic, ax, baltic_extent, coast=coast_baltic)
plt.show()
```

# Panel B — Per-HELCOM-subbasin density

Release subbasin comes from the baltic key via hex_id → helcom_subbasin;
subbasin -1 means "outside all HELCOM polygons" and is excluded.

```python
release_subbasin = baltic_key.set_index("hex_id")["helcom_subbasin"]
counts_with_subbasin = baltic_counts.assign(
    release_subbasin=baltic_counts["release_hex"].map(release_subbasin)
)
counts_with_subbasin = counts_with_subbasin[counts_with_subbasin["release_subbasin"] >= 0]

by_subbasin = (
    counts_with_subbasin
    .groupby(["release_subbasin", "target_hex"])["n_obs"].sum()
    .rename("count").reset_index()
)

# Subbasins with at least one count row, ordered by integer ID.
subbasins_ordered = [
    (sid, name) for sid, name in sorted(subbasin_id_to_name.items())
    if sid >= 0 and sid in by_subbasin["release_subbasin"].values
]
print(f"Subbasins with counts: {len(subbasins_ordered)}")

ncols = 4
nrows = int(np.ceil(len(subbasins_ordered) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(
        baltic_panel_height_in * baltic_aspect * ncols,
        baltic_panel_height_in * nrows,
    ),
    layout="constrained",
)
for ax, (sid, sname) in zip(axes.flat, subbasins_ordered):
    sub_counts = by_subbasin[by_subbasin["release_subbasin"] == sid]
    gdf_sub = (
        sub_counts
        .merge(baltic_key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    overlay = subbasins_smoothed[subbasins_smoothed["subbasin"] == sname]
    log_density_plot(
        gdf_sub, ax, baltic_extent, title=sname,
        overlay=overlay, coast=coast_baltic,
    )
for ax in axes.flat[len(subbasins_ordered):]:
    ax.set_visible(False)
plt.show()
```

# Panel C — German-waters zoom

Uses the DE store (`de_r4km_v1`) at 4 km resolution for finer detail in
the German Bight.

```python
gdf_de = (
    de_counts.groupby("target_hex")["n_obs"].sum()
    .rename("count")
    .reset_index()
    .merge(de_key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
    .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
)

de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
de_aspect = lonlat_aspect(de_extent)
fig, ax = plt.subplots(
    figsize=(de_panel_height_in * de_aspect, de_panel_height_in),
    layout="constrained",
)
log_density_plot(gdf_de, ax, de_extent, coast=coast_de)
plt.show()
```

# Panel D — Per-release-quarter density

Quarter is derived from release_doy + release_year; no pre-stored column
needed.

```python
doy_to_quarter = (
    pd.to_datetime(str(release_year)) + pd.to_timedelta(baltic_counts["release_doy"] - 1, unit="D")
).dt.quarter
baltic_counts_q = baltic_counts.assign(quarter=doy_to_quarter)
by_quarter = (
    baltic_counts_q.groupby(["quarter", "target_hex"])["n_obs"].sum()
    .rename("count").reset_index()
)

ncols, nrows = 2, 2
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(
        baltic_panel_height_in * baltic_aspect * ncols,
        baltic_panel_height_in * nrows,
    ),
    layout="constrained",
)
for ax, (q_int, q_label) in zip(axes.flat, QUARTER_LABELS.items()):
    q_counts = by_quarter[by_quarter["quarter"] == q_int]
    if q_counts.empty:
        ax.set_visible(False)
        continue
    gdf_q = (
        q_counts
        .merge(baltic_key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    log_density_plot(gdf_q, ax, baltic_extent, title=q_label, coast=coast_baltic)
plt.show()
```

# Validation prints

```python
print(f"regime={regime}, release_year={release_year}")
print()
print(f"Baltic store ({baltic_config_name}):")
print(f"  rows in counts: {len(baltic_counts):,}")
print(f"  sum(n_obs): {int(baltic_counts.n_obs.sum()):,}")
print(f"  unique target_hex in Panel A: {gdf_baltic['hex_id'].nunique():,}")
print(f"  unique release subbasins in Panel B: {len(subbasins_ordered)}")
print(f"  quarters with data in Panel D: {sorted(by_quarter['quarter'].unique().tolist())}")
print()
print(f"DE store ({de_config_name}):")
print(f"  rows in counts: {len(de_counts):,}")
print(f"  sum(n_obs): {int(de_counts.n_obs.sum()):,}")
print(f"  unique target_hex in Panel C: {gdf_de['hex_id'].nunique():,}")
```
