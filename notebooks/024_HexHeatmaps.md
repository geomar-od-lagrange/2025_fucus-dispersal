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

# Hex heatmaps (hextraj)

Particle-density maps on a hexagonal equal-area grid via
[`hextraj`](https://pypi.org/project/hextraj/). One regime per run
(papermill parameter). Scopes: whole Baltic, per HELCOM release
subbasin, German waters, per release quarter (JFM/AMJ/JAS/OND).

Compared to 023 this uses a Lambert Azimuthal Equal-Area hex grid
(equal-area cells, nicer for visual density comparison at Baltic
latitudes) instead of lon/lat rectangles. Counts only — no mean-age
(`hextraj.hex_counts` does unweighted counts).

```python
import dask
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path

from hextraj import HexProj
from hextraj.hex_analysis import hex_counts

from helpers import (
    attach_release_metadata,
    load_trajectories,
    mask_land_seeded,
    relabel_quarter,
    QUARTER_LABELS,
)
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
experiment_type = "surface"

# Baltic bounding box and hex grid.
lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
hex_size_meters = 10_000  # corner-to-centre radius
hex_origin_lon = 18
hex_origin_lat = 59

# German waters zoom.
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
hex_size_meters_de = 4_000

cmap = "viridis"

# Panel sizing probed for ~3" width per panel at standard dpi, Baltic box
# (27x13 deg). GeoDataFrame.plot enforces equal aspect by default.
single_baltic_figsize = (12, 6)
single_de_figsize = (12, 4.5)
subbasin_grid_figsize = (54, 27)
quarter_grid_figsize = (27, 13.5)
```

# Dask cluster

Connect to an external scheduler when ``SCHEDULER_FILE`` is set (written
by the multi-task SLURM job). Otherwise spin up a local cluster on the
current node.

```python
import os
import time
from dask.distributed import Client

scheduler_file = os.environ.get("SCHEDULER_FILE")
if scheduler_file:
    for _ in range(60):
        if os.path.exists(scheduler_file):
            break
        time.sleep(1)
    client = Client(scheduler_file=scheduler_file)
else:
    client = Client(ip="0.0.0.0")
client
```

# Hex projection

Lambert Azimuthal Equal-Area centred on the Baltic; ≈10 km hex radius
covers the basin at visually useful resolution without overwhelming the
aggregator.

```python
hp_baltic = HexProj(
    projection_name="laea",
    lon_origin=hex_origin_lon,
    lat_origin=hex_origin_lat,
    hex_size_meters=hex_size_meters,
)
hp_de = HexProj(
    projection_name="laea",
    lon_origin=(de_lon_min + de_lon_max) / 2,
    lat_origin=(de_lat_min + de_lat_max) / 2,
    hex_size_meters=hex_size_meters_de,
)
```

# Release area and subbasins

```python
base_path = Path(base_path)
subbasins = gpd.read_file(
    base_path / "data" / "HELCOM_subbasins_2022_level2" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(epsg=4326).rename(dict(level_2="subbasin"), axis=1)
subbasins
```

# Load trajectories and attach metadata

```python
trajectory_path = base_path / "output" / "Trajectories" / experiment_type
ds, zarr_files = load_trajectories(trajectory_path)
print(f"{len(zarr_files)} trajectory files for {experiment_type}")
ds, _ = mask_land_seeded(ds)
ds = attach_release_metadata(ds, subbasins)
ds
```

# Hex labelling (lazy)

`hp.label(lon, lat)` is pure numpy; wrap in `apply_ufunc` with
`dask="parallelized"` so the full `(trajectory, obs)` hex-ID array stays
lazy, one chunk at a time. Scope masks and `hex_counts` below reduce
on top of this lazy array.

```python
def label_hexes(ds_, hp):
    return xr.apply_ufunc(
        hp.label, ds_.lon, ds_.lat,
        dask="parallelized",
        output_dtypes=[np.int64],
    )

hex_ids_baltic = label_hexes(ds, hp_baltic)
hex_ids_de = label_hexes(ds, hp_de)
```

# Precompute per-trajectory scope keys

```python
sb_np, quarter_np = dask.compute(ds.subbasin, ds.release_quarter)
sb_np = sb_np.values
quarter_np = quarter_np.values

subbasins_list = sorted({s for s in sb_np if isinstance(s, str)})
quarters = sorted({int(q) for q in quarter_np if not np.isnan(q)})
```

# Count helpers

`counts_for` delegates aggregation and geometry attachment to
`hextraj.hex_counts`, which handles dask-backed `(trajectory, obs)`
inputs lazily without materialising the full array. `counts_by_scope`
adds a scope column (subbasin, release quarter) and replaces N
per-scope scans with a single `groupby(["scope", "hex_id"]).size()`
over the same data — hextraj's `reduce_dims` only groups by
DataArray dims, not by `(trajectory,)`-valued categoricals, so the
scope path stays in the notebook.

```python
def counts_for(hex_ids, hp):
    """Aggregate (trajectory, obs) hex IDs into a GeoDataFrame."""
    return hex_counts(hex_ids, hp=hp).dropna(subset=["geometry"])


def counts_by_scope(hex_ids, scope, hp):
    """Hex counts per scope value, computed in a single groupby pass.

    ``scope`` is a 1-D ``(trajectory,)`` DataArray; xarray broadcasts
    it to ``hex_ids`` shape at dataframe-conversion time. One groupby
    over ``(scope, hex_id)`` replaces N independent value_counts
    scans. Returns ``{scope_value: GeoDataFrame}``.
    """
    # Dask-backed accessors like ``.dt.quarter`` on a NaT-masked time
    # array advertise int64 meta while actually holding float-with-NaN;
    # to_dask_dataframe then crashes on the partition int cast. Repair
    # per-chunk via apply_ufunc, staying lazy so scope keeps its
    # trajectory chunking. (A .compute() here would turn scope into a
    # graph literal and pin the groupby tail to a single worker.)
    sentinel = None
    if np.issubdtype(scope.dtype, np.number):
        sentinel = -1
        scope = xr.apply_ufunc(
            lambda x: np.where(
                np.isnan(np.asarray(x, dtype=float)), -1, x
            ).astype(np.int64),
            scope,
            dask="parallelized",
            output_dtypes=[np.int64],
        )

    frame = xr.Dataset(
        {"hex_id": hex_ids, "scope": scope}
    ).to_dask_dataframe(dim_order=list(hex_ids.dims))
    frame = frame[frame["hex_id"] >= 0]
    if sentinel is not None:
        frame = frame[frame["scope"] != sentinel]
    counts = frame.groupby(["scope", "hex_id"]).size().rename("count").compute()

    unique_ids = np.asarray(counts.index.unique(level="hex_id"), dtype=np.int64)
    geo = hp.to_geodataframe(unique_ids)

    result = {}
    for scope_val, sub in counts.groupby(level="scope"):
        hex_id_values = sub.index.get_level_values("hex_id").to_numpy()
        gdf = gpd.GeoDataFrame(
            {"count": sub.to_numpy(),
             "geometry": geo.geometry.reindex(hex_id_values).values},
            index=pd.Index(hex_id_values, name="hex_id"),
            crs="EPSG:4326",
        ).dropna(subset=["geometry"])
        result[scope_val] = gdf
    return result


def hp_to_cartopy(hp):
    """Cartopy CRS matching a hextraj ``HexProj``."""
    if hp.projection_name == "laea":
        return ccrs.LambertAzimuthalEqualArea(
            central_longitude=hp.lon_origin,
            central_latitude=hp.lat_origin,
        )
    raise ValueError(f"Unsupported HexProj projection: {hp.projection_name}")


def log_density_plot(gdf, ax, extent, title=None, overlay=None):
    # Reproject the hex GDF to the axis projection once so matplotlib
    # can draw native polygons without cartopy re-projecting every edge
    # at draw time. Paired with ``edgecolor="face"`` this removes the
    # anti-aliased seams between adjacent hexes.
    gdf = gdf.to_crs(ax.projection)
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
    ax.coastlines(resolution="50m", color="black", linewidth=0.5, zorder=2)
    if overlay is not None:
        overlay.boundary.plot(
            ax=ax,
            color="black",
            linewidth=0.7,
            transform=ccrs.PlateCarree(),
            zorder=3,
        )
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    if title is not None:
        ax.set_title(title)
```

# Whole Baltic

```python
gdf_baltic = counts_for(hex_ids_baltic, hp_baltic)
fig, ax = plt.subplots(
    figsize=single_baltic_figsize,
    subplot_kw={"projection": hp_to_cartopy(hp_baltic)},
)
log_density_plot(gdf_baltic, ax, [lon_min, lon_max, lat_min, lat_max])
plt.show()
```

# Per HELCOM release subbasin

```python
gdfs_by_basin = counts_by_scope(hex_ids_baltic, ds.subbasin, hp_baltic)

ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=subbasin_grid_figsize,
    subplot_kw={"projection": hp_to_cartopy(hp_baltic)},
)
for ax, basin in zip(axes.flat, subbasins_list):
    gdf = gdfs_by_basin.get(basin)
    if gdf is None or gdf.empty:
        ax.set_visible(False)
        continue
    overlay = subbasins[subbasins["subbasin"] == basin]
    log_density_plot(
        gdf, ax,
        [lon_min, lon_max, lat_min, lat_max],
        title=basin,
        overlay=overlay,
    )
for ax in axes.flat[len(subbasins_list):]:
    ax.set_visible(False)
plt.show()
```

# German waters

Same trajectories, finer hex grid (`hex_size_meters_de`) for the
German-waters zoom.

```python
gdf_de = counts_for(hex_ids_de, hp_de)
fig, ax = plt.subplots(
    figsize=single_de_figsize,
    subplot_kw={"projection": hp_to_cartopy(hp_de)},
)
log_density_plot(gdf_de, ax, [de_lon_min, de_lon_max, de_lat_min, de_lat_max])
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
gdfs_by_quarter = counts_by_scope(hex_ids_baltic, ds.release_quarter, hp_baltic)

ncols, nrows = 2, 2
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=quarter_grid_figsize,
    subplot_kw={"projection": hp_to_cartopy(hp_baltic)},
)
for ax, (q_int, q_label) in zip(axes.flat, QUARTER_LABELS.items()):
    gdf = gdfs_by_quarter.get(q_int)
    if gdf is None or gdf.empty:
        ax.set_visible(False)
        continue
    log_density_plot(gdf, ax, [lon_min, lon_max, lat_min, lat_max], title=q_label)
plt.show()
```
