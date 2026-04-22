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

# Colour scale (log density).
cmap = "YlOrRd"

# Panel sizing probed for ~3" width per panel at standard dpi, Baltic box
# (27x13 deg). GeoDataFrame.plot enforces equal aspect by default.
single_baltic_figsize = (4, 2)
single_de_figsize = (4, 1.5)
subbasin_grid_figsize = (18, 9)
quarter_grid_figsize = (9, 4.5)
```

# Dask cluster

```python
from dask.distributed import Client
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

`hex_counts` materialises the (trajectory, obs) hex-ID array it is
given; we therefore pre-subset with `ds.where(key==v)` so each scope
only computes its own slice. For the full Baltic the scope is the
whole dataset.

```python
def counts_for(hex_ids, hp, mask=None):
    """Distributed-friendly hex counting.

    Instead of materialising the full (trajectory, obs) dask array in one
    shot (OOM on large regimes), flatten to a dask Series and let
    dask.dataframe.value_counts aggregate per partition, then build
    geometries from the small result."""
    if mask is not None:
        hex_ids = hex_ids.where(mask, other=-1)
    flat = hex_ids.data.ravel()
    vc = dd.from_dask_array(flat, columns="hex_id")["hex_id"].value_counts().compute()
    vc = vc[vc.index >= 0]
    gdf = hex_counts(pd.Series(vc.index, dtype=np.int64), hp=hp)
    gdf["count"] = vc.reindex(gdf.index).values
    return gdf.dropna(subset=["geometry"])


def log_density_plot(gdf, ax, extent, title=None):
    gdf = gdf.copy()
    gdf["log_count"] = np.log10(gdf["count"].where(gdf["count"] > 0))
    gdf.plot(
        ax=ax,
        column="log_count",
        cmap=cmap,
        legend=False,
        missing_kwds={"color": "none"},
    )
    ax.set_xlim(*extent[:2])
    ax.set_ylim(*extent[2:])
    if title is not None:
        ax.set_title(title)
```

# Whole Baltic

```python
gdf_baltic = counts_for(hex_ids_baltic, hp_baltic)
fig, ax = plt.subplots(figsize=single_baltic_figsize)
log_density_plot(gdf_baltic, ax, [lon_min, lon_max, lat_min, lat_max])
plt.show()
```

# Per HELCOM release subbasin

```python
ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=subbasin_grid_figsize,
)
for ax, basin in zip(axes.flat, subbasins_list):
    gdf = counts_for(hex_ids_baltic, hp_baltic, mask=(ds.subbasin == basin))
    log_density_plot(gdf, ax, [lon_min, lon_max, lat_min, lat_max], title=basin)
for ax in axes.flat[len(subbasins_list):]:
    ax.set_visible(False)
plt.show()
```

# German waters

Same trajectories, finer hex grid (`hex_size_meters_de`) for the
German-waters zoom.

```python
gdf_de = counts_for(hex_ids_de, hp_de)
fig, ax = plt.subplots(figsize=single_de_figsize)
log_density_plot(gdf_de, ax, [de_lon_min, de_lon_max, de_lat_min, de_lat_max])
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
ncols, nrows = 2, 2
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=quarter_grid_figsize,
)
for ax, (q_int, q_label) in zip(axes.flat, QUARTER_LABELS.items()):
    gdf = counts_for(hex_ids_baltic, hp_baltic, mask=(ds.release_quarter == q_int))
    log_density_plot(gdf, ax, [lon_min, lon_max, lat_min, lat_max], title=q_label)
plt.show()
```
