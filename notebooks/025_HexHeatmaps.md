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

Particle-density maps on the hex-aggregated dispersal store built by
notebook 024. Reads `key.parquet` and `counts/` partitions from one
pre-built store; no trajectory zarrs, no Dask cluster. Panel C is the
same data as Panel A, viewport-clipped to the German Bight — there is
one hex size for the whole basin.

Four panels per run:

- **A** — whole-Baltic density (all release hexes summed)
- **B** — per-HELCOM-subbasin density (one sub-plot per release subbasin)
- **C** — German-waters zoom (Panel A data clipped to `de_extent`)
- **D** — per-release-quarter density (2×2 grid)

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
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (HELCOM polygons).
data_root = "../data"
# Read root of the hex-aggregate stores built by notebook 024.
output_root = "../output"
# Which regime's counts partition to read; one regime per run.
regime = "surface"
# Which release_year partition to read; one year per run.
release_year = 2019
# Hex radius of the aggregate store (built by 024). One radius covers
# the whole BSH domain; Panel C is a viewport clip, not a separate store.
hex_radius = 6000

# Baltic-wide map extent (degrees E / degrees N).
baltic_lon_min = 5
baltic_lon_max = 32
baltic_lat_min = 53
baltic_lat_max = 66

# DE-zoom map extent (degrees E / degrees N).
de_lon_min = 8
de_lon_max = 15
de_lat_min = 53.2
de_lat_max = 55.5

# Colormap (TODO(phase-f): justify in docs/visualisations.md — log
# density is not well-served by the matplotlib default).
cmap = "viridis"

# Per-panel height in inches (panel widths are aspect-derived).
baltic_panel_height_in = 6
de_panel_height_in = 4.5
```

# Path construction and parquet reads

Layout assumption: ``output_root/HexAggregates/<config>/key.parquet``
plus partitions at ``counts/regime=<regime>/release_year=<year>/part.parquet``.

```python
data_root = Path(data_root)
output_root = Path(output_root)
store = output_root / "HexAggregates" / f"r{hex_radius}m"

key = gpd.read_parquet(store / "key.parquet")
counts = pd.read_parquet(
    store / "counts" / f"regime={regime}" / f"release_year={release_year}" / "part.parquet"
)

# Read subbasin_id_to_name from the key.parquet file-level metadata.
_key_meta = json.loads(
    pq.read_table(store / "key.parquet").schema.metadata[b"hex_aggregate_store"].decode()
)
subbasin_id_to_name = {int(k): v for k, v in _key_meta["subbasin_id_to_name"].items()}

print(f"key: {len(key):,} hexes")
print(f"counts: {len(counts):,} rows, sum(n_obs)={int(counts.n_obs.sum()):,}")
```

# Rendering helpers

Coast and subbasin data are loaded once here and reused across all
panels. ``to_hex_gdf`` joins a ``(target_hex, count)`` table to a key
GeoDataFrame; four callsites below justify the local def.

```python
def to_hex_gdf(counts_df, key_df):
    """Build a hex GeoDataFrame from per-target_hex counts + a key
    GeoDataFrame. Sums ``n_obs`` over ``target_hex`` and joins to
    ``key_df`` on ``target_hex`` ↔ ``hex_id``."""
    return (
        counts_df.groupby("target_hex")["n_obs"].sum()
        .rename("count").reset_index()
        .merge(key_df[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def log_density_plot(gdf, ax, extent, title=None, overlay=None, coast=None):
    """Plot hex density, coastline, and optional subbasin overlay on a
    plain (non-cartopy) lon/lat axis. Everything stays in EPSG:4326 so
    hexes, coastline, and overlay share one coordinate system — no
    drift from mismatched reprojection pipelines."""
    gdf = gdf.copy()
    gdf["log_count"] = np.log10(gdf["count"].where(gdf["count"] > 0))
    # TODO(phase-f): justify cmap, edgecolor="face", linewidth in docs/visualisations.md.
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
        # TODO(phase-f): justify color="black", linewidth in docs/visualisations.md.
        coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    if overlay is not None:
        # TODO(phase-f): justify color="magenta", linewidth in docs/visualisations.md.
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
coast_baltic = _coast_gdf.clip(box(baltic_lon_min, baltic_lat_min, baltic_lon_max, baltic_lat_max))
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

baltic_extent = [baltic_lon_min, baltic_lon_max, baltic_lat_min, baltic_lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
# Aspect ratio that keeps 1° lon at lat_mean visually equal to 1° lat.
baltic_aspect = (
    (baltic_lon_max - baltic_lon_min) * np.cos(np.radians(0.5 * (baltic_lat_min + baltic_lat_max)))
) / (baltic_lat_max - baltic_lat_min)
de_aspect = (
    (de_lon_max - de_lon_min) * np.cos(np.radians(0.5 * (de_lat_min + de_lat_max)))
) / (de_lat_max - de_lat_min)
```

# Panel A — Whole-Baltic density

```python
gdf_full = to_hex_gdf(counts, key)

fig, ax = plt.subplots(
    figsize=(baltic_panel_height_in * baltic_aspect, baltic_panel_height_in),
    layout="constrained",
)
log_density_plot(gdf_full, ax, baltic_extent, coast=coast_baltic)
plt.show()
```

# Panel B — Per-HELCOM-subbasin density

Release subbasin comes from the key via hex_id → helcom_subbasin;
subbasin -1 means "outside all HELCOM polygons" and is excluded.

```python
release_subbasin = key.set_index("hex_id")["helcom_subbasin"]
counts_with_subbasin = counts.assign(
    release_subbasin=counts["release_hex"].map(release_subbasin)
)
counts_with_subbasin = counts_with_subbasin[counts_with_subbasin["release_subbasin"] >= 0]

# Iterate every named subbasin in stable id order. Empty panels render
# fine and keep the layout stable across regimes / release_years.
subbasin_pairs = [(sid, name) for sid, name in sorted(subbasin_id_to_name.items()) if sid >= 0]

ncols = 4
nrows = int(np.ceil(len(subbasin_pairs) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(
        baltic_panel_height_in * baltic_aspect * ncols,
        baltic_panel_height_in * nrows,
    ),
    layout="constrained",
)
for ax, (sid, sname) in zip(axes.flat, subbasin_pairs):
    sub_counts = counts_with_subbasin[counts_with_subbasin["release_subbasin"] == sid]
    gdf_sub = to_hex_gdf(sub_counts, key)
    overlay = subbasins_smoothed[subbasins_smoothed["subbasin"] == sname]
    log_density_plot(
        gdf_sub, ax, baltic_extent, title=sname,
        overlay=overlay, coast=coast_baltic,
    )
for ax in axes.flat[len(subbasin_pairs):]:
    ax.set_visible(False)
plt.show()
```

# Panel C — German-waters zoom

Same data as Panel A; viewport clipped to ``de_extent``. Hex size is
fixed by the store, so the zoom is purely visual.

```python
fig, ax = plt.subplots(
    figsize=(de_panel_height_in * de_aspect, de_panel_height_in),
    layout="constrained",
)
log_density_plot(gdf_full, ax, de_extent, coast=coast_de)
plt.show()
```

# Panel D — Per-release-quarter density

Quarter is derived from `release_doy` via the calendar (handles
leap-year boundary asymmetry between Q1/Q2/Q3/Q4 month lengths).

```python
release_dates = pd.to_datetime(
    (release_year * 1000 + counts["release_doy"]).astype(str),
    format="%Y%j",
)
counts_q = counts.assign(quarter=release_dates.dt.quarter)

quarter_labels = {1: "JFM", 2: "AMJ", 3: "JAS", 4: "OND"}
ncols, nrows = 2, 2
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(
        baltic_panel_height_in * baltic_aspect * ncols,
        baltic_panel_height_in * nrows,
    ),
    layout="constrained",
)
for ax, (q_int, q_label) in zip(axes.flat, quarter_labels.items()):
    q_counts = counts_q[counts_q["quarter"] == q_int]
    gdf_q = to_hex_gdf(q_counts, key)
    log_density_plot(gdf_q, ax, baltic_extent, title=q_label, coast=coast_baltic)
plt.show()
```

# Validation prints

```python
print(f"regime={regime}, release_year={release_year}, hex_radius={hex_radius} m")
print()
print(f"  rows in counts: {len(counts):,}")
print(f"  sum(n_obs): {int(counts.n_obs.sum()):,}")
print(f"  unique target_hex in Panel A: {gdf_full['hex_id'].nunique():,}")
print(f"  named release subbasins in Panel B: {len(subbasin_pairs)}")
print(f"  quarters with data in Panel D: {sorted(counts_q['quarter'].unique().tolist())}")
```
