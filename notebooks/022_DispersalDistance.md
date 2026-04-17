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

# Dispersal distance vs. time

Mean displacement from release point as a function of particle age.
Regimes overlaid via `hue=`. Scopes: global, per HELCOM release
subbasin, per release season, per release year.

```python
import numpy as np
import xarray as xr
import geopandas as gpd
import cartopy.crs as ccrs
from pathlib import Path
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
```

# Load all regimes + land-seed filter

```python
base_path = Path(base_path)
trajectory_root = base_path / "output" / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

```python
def load_regime(regime):
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    try:
        ds = ds.rename(dict(CellID="cell_ID"))
    except ValueError:
        pass
    dlon0 = ds.lon.diff("obs").isel(obs=0)
    dlat0 = ds.lat.diff("obs").isel(obs=0)
    land = ((dlon0 == 0) & (dlat0 == 0)).compute()
    ds = ds.isel(trajectory=~land)
    ds = ds.expand_dims(dim={"regime": [regime]}).squeeze(drop=False)
    return ds
```

# Release metadata

```python
release_area = gpd.read_file(
    base_path / "data" / "Fucus_location_shp" / "REDLIST_SIS_Macrophytes.shp"
)
release_area = release_area.loc[
    release_area.F_vesiculo != 0, ["geometry", "CELLCODE", "CELLID"]
]
subbasins = gpd.read_file(
    base_path / "data" / "HELCOM_subbasins_2022_level2" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(crs=ccrs.Geodetic()).rename(dict(level_2="subbasin"), axis=1)
release_centroid = gpd.GeoDataFrame(
    release_area.CELLCODE, geometry=release_area.geometry.centroid
)
release_to_subbasin = gpd.sjoin_nearest(
    release_centroid, subbasins, how="left"
)[["CELLCODE", "subbasin", "geometry"]]
```

```python
def attach_metadata(ds):
    cell_ids = ds.cell_ID.isel(obs=0).astype(int).compute().values
    subbasin_per_traj = release_to_subbasin.loc[cell_ids].subbasin.values
    release_time = ds.time.isel(obs=0).compute()
    month = release_time.dt.month
    season = xr.where(
        month.isin([12, 1, 2]), "DJF",
        xr.where(month.isin([3, 4, 5]), "MAM",
        xr.where(month.isin([6, 7, 8]), "JJA", "SON")),
    )
    return ds.assign(
        subbasin=("trajectory", subbasin_per_traj),
        release_year=("trajectory", release_time.dt.year.values),
        release_season=("trajectory", season.values),
    )
```

# Distance from release

```python
def distance_km(ds):
    lon0 = ds.lon.isel(obs=0)
    lat0 = ds.lat.isel(obs=0)
    dlat = ds.lat - lat0
    dlon = (ds.lon - lon0) * np.cos(np.deg2rad(lat0))
    return (111.0 * np.sqrt(dlat ** 2 + dlon ** 2)).rename("distance_km")
```

# Load all regimes (lazy)

```python
regime_dsets = {}
for regime in regimes:
    regime_dsets[regime] = attach_metadata(load_regime(regime))
```

# Dask cluster

```python
from dask.distributed import Client
client = Client()
client
```

# Per-regime distance (lazy)

```python
regime_distance = {regime: distance_km(ds) for regime, ds in regime_dsets.items()}
```

# Global

```python
da_global = xr.concat(
    [regime_distance[r].mean("trajectory").expand_dims(regime=[r]) for r in regimes],
    dim="regime",
)
da_global.plot.line(x="obs", hue="regime")
```

# Per HELCOM release subbasin

```python
by_sb = []
for regime in regimes:
    ds = regime_dsets[regime]
    d = regime_distance[regime].assign_coords(subbasin=("trajectory", ds.subbasin.values))
    m = d.groupby("subbasin").mean("trajectory")
    by_sb.append(m.expand_dims(regime=[regime]))
da_sb = xr.concat(by_sb, dim="regime")
da_sb.plot.line(x="obs", hue="regime", col="subbasin", col_wrap=4)
```

# German waters (release cells inside bounding box)

```python
def in_de(ds):
    lon0 = ds.lon.isel(obs=0).compute()
    lat0 = ds.lat.isel(obs=0).compute()
    return ((lon0 >= de_lon_min) & (lon0 <= de_lon_max)
            & (lat0 >= de_lat_min) & (lat0 <= de_lat_max))

da_de = xr.concat(
    [
        regime_distance[r].isel(trajectory=in_de(regime_dsets[r]).values).mean("trajectory").expand_dims(regime=[r])
        for r in regimes
    ],
    dim="regime",
)
da_de.plot.line(x="obs", hue="regime")
```

# Per release season

```python
by_season = []
for regime in regimes:
    ds = regime_dsets[regime]
    d = regime_distance[regime].assign_coords(release_season=("trajectory", ds.release_season.values))
    m = d.groupby("release_season").mean("trajectory")
    by_season.append(m.expand_dims(regime=[regime]))
da_season = xr.concat(by_season, dim="regime")
da_season.plot.line(x="obs", hue="regime", col="release_season")
```

# Per release year

```python
by_year = []
for regime in regimes:
    ds = regime_dsets[regime]
    d = regime_distance[regime].assign_coords(release_year=("trajectory", ds.release_year.values))
    m = d.groupby("release_year").mean("trajectory")
    by_year.append(m.expand_dims(regime=[regime]))
da_year = xr.concat(by_year, dim="regime")
da_year.plot.line(x="obs", hue="regime", col="release_year", col_wrap=4)
```
