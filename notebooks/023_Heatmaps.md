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

# Heatmaps

Particle-density + mean-age maps from one lazy xhistogram (obs kept
as a dim) per regime. One regime per run (papermill parameter).
Scopes: whole Baltic, per HELCOM release subbasin, German waters,
per release season, per release year.

```python
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from xhistogram.xarray import histogram as xhist
from pathlib import Path
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
experiment_type = "surface"

output_dt_mins = 60

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
n_bins_baltic = 120

de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
de_bin_width = 0.125
```

# Load regime and filter land-seeded

```python
base_path = Path(base_path)
trajectory_path = base_path / "output" / "Trajectories" / experiment_type
zarr_files = sorted(trajectory_path.glob("**/*.zarr"))
print(f"{len(zarr_files)} trajectory files for {experiment_type}")
ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")

try:
    ds = ds.rename(dict(CellID="cell_ID"))
except ValueError:
    pass

dlon0 = ds.lon.diff("obs").isel(obs=0)
dlat0 = ds.lat.diff("obs").isel(obs=0)
land = (
    ((dlon0 == 0) & (dlat0 == 0)).drop_vars("obs", errors="ignore").compute()
)
print(f"Dropping {int(land.sum())} land-seeded of {ds.sizes['trajectory']}")
ds = ds.isel(trajectory=~land)
```

# Attach release metadata

```python
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
ds = ds.assign(subbasin=("trajectory", release_to_subbasin.loc[cell_ids].subbasin.values))

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

# Loaded data

```python
ds
```

```python
release_area
```

```python
subbasins
```

```python
release_to_subbasin
```

# Dask cluster

```python
from dask.distributed import Client
client = Client()
client
```

# Bins and histogram helpers

```python
lon_bins = np.linspace(lon_min, lon_max, n_bins_baltic)
lat_bins = np.linspace(lat_min, lat_max, n_bins_baltic)
de_lon_bins = np.arange(de_lon_min, de_lon_max + de_bin_width, de_bin_width)
de_lat_bins = np.arange(de_lat_min, de_lat_max + de_bin_width, de_bin_width)

age_hours_per_obs = output_dt_mins / 60
```

```python
def count_hist(ds_, lon_bins, lat_bins):
    return xhist(
        ds_.lon, ds_.lat,
        bins=[lon_bins, lat_bins],
        dim=["trajectory"],
    ).rename(dict(lon_bin="lon", lat_bin="lat"))

def density(h):
    return h.sum("obs")

def mean_age_hours(h):
    return ((h * h.obs).sum("obs") / h.sum("obs")) * age_hours_per_obs
```

```python
def facet_map(da, col, col_wrap=None):
    fg = da.plot(
        x="lon", y="lat",
        col=col, col_wrap=col_wrap,
        subplot_kws=dict(projection=ccrs.PlateCarree()),
        transform=ccrs.PlateCarree(),
    )
    for ax in fg.axs.flat:
        ax.coastlines()
    return fg

def single_map(da, extent):
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))
    da.plot(ax=ax, x="lon", y="lat", transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    return fig
```

# Whole Baltic

```python
h_baltic = count_hist(ds, lon_bins, lat_bins)
single_map(density(h_baltic), [lon_min, lon_max, lat_min, lat_max])
plt.show()
single_map(mean_age_hours(h_baltic), [lon_min, lon_max, lat_min, lat_max])
plt.show()
```

# Per HELCOM release subbasin

```python
h_by_sb = ds.groupby("subbasin").apply(
    lambda d: count_hist(d, lon_bins, lat_bins)
)
facet_map(density(h_by_sb), col="subbasin", col_wrap=4)
plt.show()
facet_map(mean_age_hours(h_by_sb), col="subbasin", col_wrap=4)
plt.show()
```

# German waters

```python
h_de = count_hist(ds, de_lon_bins, de_lat_bins)
single_map(density(h_de), [de_lon_min, de_lon_max, de_lat_min, de_lat_max])
plt.show()
single_map(mean_age_hours(h_de), [de_lon_min, de_lon_max, de_lat_min, de_lat_max])
plt.show()
```

# Per release season

```python
h_by_season = ds.groupby("release_season").apply(
    lambda d: count_hist(d, lon_bins, lat_bins)
)
facet_map(density(h_by_season), col="release_season")
plt.show()
facet_map(mean_age_hours(h_by_season), col="release_season")
plt.show()
```

# Per release year

```python
h_by_year = ds.groupby("release_year").apply(
    lambda d: count_hist(d, lon_bins, lat_bins)
)
facet_map(density(h_by_year), col="release_year", col_wrap=4)
plt.show()
facet_map(mean_age_hours(h_by_year), col="release_year", col_wrap=4)
plt.show()
```
