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
waters, per release season, per release year.

```python
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
experiment_type = "surface"

n_traj_subset = 1000
seed = 0

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
```

# Load regime and filter land-seeded

```python
base_path = Path(base_path)
trajectory_path = base_path / "output" / "Trajectories" / experiment_type
zarr_files = sorted(trajectory_path.glob("**/*.zarr"))
print(f"{len(zarr_files)} trajectory files for {experiment_type}")
ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")

dlon0 = ds.lon.diff("obs").isel(obs=0)
dlat0 = ds.lat.diff("obs").isel(obs=0)
land_seeded = (
    ((dlon0 == 0) & (dlat0 == 0)).drop_vars("obs", errors="ignore").compute()
)
n_land = int(land_seeded.sum())
print(f"Dropping {n_land} land-seeded trajectories of {ds.sizes['trajectory']}")
ds = ds.isel(trajectory=~land_seeded)
```

# Attach release metadata

```python
try:
    ds = ds.rename(dict(CellID="cell_ID"))
except ValueError:
    pass

release_area = gpd.read_file(
    base_path / "data" / "Fucus_location_shp" / "REDLIST_SIS_Macrophytes.shp"
)
release_area = release_area.assign(CELLID=release_area.index.astype(int))
release_area = release_area.loc[
    release_area.F_vesiculo != 0, ["geometry", "CELLCODE", "CELLID"]
].to_crs(crs=ccrs.Geodetic())

subbasins = gpd.read_file(
    base_path / "data" / "HELCOM_subbasins_2022_level2" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(crs=ccrs.Geodetic()).rename(dict(level_2="subbasin"), axis=1)

release_centroid = gpd.GeoDataFrame(
    release_area.CELLCODE, geometry=release_area.geometry.centroid
)
release_to_subbasin = gpd.sjoin_nearest(
    release_centroid, subbasins, how="left"
)[["CELLCODE", "subbasin", "geometry"]]

cell_ids = ds.cell_ID.isel(obs=0).astype(int).compute().values
subbasin_per_traj = release_to_subbasin.loc[cell_ids].subbasin.values
ds = ds.assign(subbasin=("trajectory", subbasin_per_traj))
```

```python
release_time = ds.time.isel(obs=0).compute()
month = release_time.dt.month
season = xr.where(
    month.isin([12, 1, 2]), "DJF",
    xr.where(month.isin([3, 4, 5]), "MAM",
    xr.where(month.isin([6, 7, 8]), "JJA", "SON")),
)
ds = ds.assign(
    release_year=("trajectory", release_time.dt.year.values),
    release_season=("trajectory", season.values),
)
```

# Dask cluster

```python
from dask.distributed import Client
client = Client()
client
```

# Subset sampling

```python
def subset_indices(n_available, n_want, seed):
    rng = np.random.default_rng(seed)
    n = min(n_available, n_want)
    return rng.choice(n_available, n, replace=False)

def plot_lines(ds_sub, ax, n, seed):
    if ds_sub.sizes["trajectory"] == 0:
        return
    idx = subset_indices(ds_sub.sizes["trajectory"], n, seed)
    ds_plot = ds_sub.isel(trajectory=idx).compute()
    lon = ds_plot.lon.values
    lat = ds_plot.lat.values
    for i in range(ds_plot.sizes["trajectory"]):
        ax.plot(lon[i], lat[i], transform=ccrs.PlateCarree())
```

# Per HELCOM release subbasin

```python
subbasins_list = sorted({s for s in ds.subbasin.values if s is not None})
ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, basin in zip(axes.flat, subbasins_list):
    sel = ds.subbasin == basin
    plot_lines(ds.isel(trajectory=sel), ax, n_traj_subset, seed)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(basin)
for ax in axes.flat[len(subbasins_list):]:
    ax.set_visible(False)
plt.show()
```

# German waters

```python
fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
plot_lines(ds, ax, n_traj_subset, seed)
ax.set_extent([de_lon_min, de_lon_max, de_lat_min, de_lat_max], crs=ccrs.PlateCarree())
ax.coastlines()
ax.set_title(f"German waters — {experiment_type}")
plt.show()
```

# Per release season

```python
seasons = ["DJF", "MAM", "JJA", "SON"]
fig, axes = plt.subplots(
    nrows=2, ncols=2,
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, s in zip(axes.flat, seasons):
    sel = ds.release_season == s
    plot_lines(ds.isel(trajectory=sel), ax, n_traj_subset, seed)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(s)
plt.show()
```

# Per release year

```python
years = sorted({int(y) for y in ds.release_year.values})
ncols = 4
nrows = int(np.ceil(len(years) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, y in zip(axes.flat, years):
    sel = ds.release_year == y
    plot_lines(ds.isel(trajectory=sel), ax, n_traj_subset, seed)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(str(y))
for ax in axes.flat[len(years):]:
    ax.set_visible(False)
plt.show()
```
