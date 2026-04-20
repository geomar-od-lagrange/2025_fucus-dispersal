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

# Raw trajectories

Polyline subsets of particle tracks on maps. One regime per run
(papermill parameter). Scopes: per HELCOM release subbasin, German
waters, per release quarter (JFM/AMJ/JAS/OND), per release year.

```python
import dask
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path

from helpers import (
    QUARTER_LABELS,
    attach_release_metadata,
    load_trajectories,
    mask_land_seeded,
)
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
experiment_type = "surface"

n_traj_subset = 1000

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5

panel_size = 4
panel_size_sub = 8
```

# Global RNG

```python
rng = np.random.default_rng()
```

# Dask cluster

```python
from dask.distributed import Client
client = Client(ip="0.0.0.0")
client
```

# Release area

```python
base_path = Path(base_path)
release_area = gpd.read_file(
    base_path / "data" / "Fucus_location_shp" / "REDLIST_SIS_Macrophytes.shp"
)
release_area = release_area.loc[
    release_area.F_vesiculo != 0, ["geometry", "CELLCODE"]
].to_crs(crs=ccrs.Geodetic())
release_area
```

# HELCOM subbasins

```python
subbasins = gpd.read_file(
    base_path / "data" / "HELCOM_subbasins_2022_level2" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(crs=ccrs.Geodetic()).rename(dict(level_2="subbasin"), axis=1)
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

# Precompute per-trajectory scope keys

`release_year` / `release_quarter` are lazy 1-D `(trajectory,)` arrays.
Compute them once so the per-panel loops below don't re-walk the graph.

```python
release_year_np, release_quarter_np, subbasin_np = dask.compute(
    ds.release_year, ds.release_quarter, ds.subbasin,
)
release_year_np = release_year_np.values
release_quarter_np = release_quarter_np.values
subbasin_np = subbasin_np.values
```

# Plot helper

One matplotlib call per panel: `ax.plot(lon, lat)` where each is a
`(obs, trajectory)` 2-D array draws one line per column.

```python
def plot_lines(ds_, sel, ax, n, lw=None):
    if sel is None:
        avail = np.arange(ds_.sizes["trajectory"])
    else:
        avail = np.flatnonzero(sel)
    if avail.size == 0:
        return
    chosen = rng.choice(avail, min(avail.size, n), replace=False)
    ds_plot = ds_.isel(trajectory=chosen).compute()
    lon = ds_plot.lon.to_pandas().T
    lat = ds_plot.lat.to_pandas().T
    ax.plot(lon, lat, lw=lw, transform=ccrs.PlateCarree())
```

# Per HELCOM release subbasin

```python
subbasins_list = sorted({s for s in subbasin_np if isinstance(s, str)})
ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(ncols * panel_size_sub, nrows * panel_size_sub),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, basin in zip(axes.flat, subbasins_list):
    sel = subbasin_np == basin
    plot_lines(ds, sel, ax, n_traj_subset, lw=0.5)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(basin)
for ax in axes.flat[len(subbasins_list):]:
    ax.set_visible(False)
plt.show()
```

# German waters

```python
fig, ax = plt.subplots(
    figsize=(panel_size, panel_size),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
plot_lines(ds, None, ax, n_traj_subset)
ax.set_extent([de_lon_min, de_lon_max, de_lat_min, de_lat_max], crs=ccrs.PlateCarree())
ax.coastlines()
ax.set_title(f"German waters — {experiment_type}")
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
ncols, nrows = 2, 2
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(ncols * panel_size, nrows * panel_size),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, (q_int, q_label) in zip(axes.flat, QUARTER_LABELS.items()):
    sel = release_quarter_np == q_int
    plot_lines(ds, sel, ax, n_traj_subset)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(q_label)
plt.show()
```

# Per release year

```python
years = sorted({int(y) for y in release_year_np if not np.isnan(y)})
ncols = 4
nrows = int(np.ceil(len(years) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(ncols * panel_size, nrows * panel_size),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, y in zip(axes.flat, years):
    sel = release_year_np == y
    plot_lines(ds, sel, ax, n_traj_subset)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(str(y))
for ax in axes.flat[len(years):]:
    ax.set_visible(False)
plt.show()
```
