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
subbasin, German waters, per release quarter (JFM/AMJ/JAS/OND).

```python
import os
import time

import dask
import numpy as np
import xarray as xr
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from dask.distributed import Client
from pathlib import Path
```

```python
def assign_release_subbasin(ds, subbasins):
    # Lazy per (trajectory,) chunk; STRtree built once, looked up per-chunk.
    # The full-sweep concat reaches 60M+ trajectories — eager would OOM.
    tree = shapely.STRtree(subbasins.geometry.values)
    names = subbasins["subbasin"].to_numpy()

    def _lookup(lon, lat):
        out = np.full(lon.shape, None, dtype=object)
        valid = ~(np.isnan(lon) | np.isnan(lat))
        if valid.any():
            pts = shapely.points(lon[valid], lat[valid])
            out[valid] = names[tree.nearest(pts)]
        return out

    lon0 = ds.lon.isel(obs=0, drop=True)
    lat0 = ds.lat.isel(obs=0, drop=True)
    subbasin = xr.apply_ufunc(
        _lookup, lon0, lat0,
        dask="parallelized", output_dtypes=[object],
    )
    return ds.assign(subbasin=subbasin)
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (HELCOM polygons, Fucus shapefile).
data_root = "../data"
# Read root of trajectory zarrs.
output_root = "../output"

# German-waters bounding box for release-cell membership (degrees).
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5

# Facet-grid panel size + aspect for line plots (xarray FacetGrid kwargs).
facet_line_size = 3.0
facet_line_aspect = 1.0
```

# Dask cluster

Connect to an external scheduler when ``SCHEDULER_FILE`` is set (written
by the multi-task SLURM job). Otherwise spin up a local cluster on the
current node.

```python
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

# List regimes

Layout assumption: ``output_root/Trajectories/<regime>/<release_year>/*.zarr``.
``regimes`` is the list of immediate subdirectories.

```python
data_root = Path(data_root)
output_root = Path(output_root)
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

# Release area

```python
release_area = gpd.read_file(
    data_root / "fucus_redlist_shapefile" / "REDLIST_SIS_Macrophytes.shp"
)
release_area = release_area.loc[
    release_area.F_vesiculo != 0, ["geometry", "CELLCODE"]
].to_crs(crs=ccrs.Geodetic())
release_area
```

# HELCOM subbasins

```python
subbasins = gpd.read_file(
    data_root / "helcom_subbasins" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(crs=ccrs.Geodetic()).rename(dict(level_2="subbasin"), axis=1)
subbasins
```

# Load each regime

```python
regime_dsets = {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    # First-step displacement of zero ⇒ trajectory was seeded on land.
    ds = ds.where(~(
        (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
    ))
    ds = ds.assign(release_quarter=ds.time.isel(obs=0, drop=True).dt.quarter)
    ds = assign_release_subbasin(ds, subbasins)
    regime_dsets[regime] = ds
regime_dsets
```

# Distance from release

Great-circle approximation (111 km per degree lat).

```python
def distance_km(ds):
    lon0 = ds.lon.isel(obs=0, drop=True)
    lat0 = ds.lat.isel(obs=0, drop=True)
    dlat = ds.lat - lat0
    dlon = (ds.lon - lon0) * np.cos(np.deg2rad(lat0))
    return (111.0 * np.sqrt(dlat ** 2 + dlon ** 2)).rename("distance_km")
```

# Compute per-scope means (one shared dask pass)

All scope means walk the same per-regime distance graph. Build them
lazily, then let `dask.compute(*)` evaluate everything in a single pass.

The group-coord (`subbasin`, `release_quarter`) is small —
`(trajectory,)`-shaped — so we materialise it once per regime before the
groupby; that keeps the rest of the graph lazy and avoids forcing a
client-side fancy index into the concat dim.

```python
regime_distance = {r: distance_km(ds) for r, ds in regime_dsets.items()}

da_global_lazy = xr.concat(
    [regime_distance[r].mean("trajectory").expand_dims(regime=[r]) for r in regimes],
    dim="regime",
)

da_de_lazy = xr.concat(
    [
        regime_distance[r].where(
            (regime_dsets[r].lon.isel(obs=0, drop=True) >= de_lon_min)
            & (regime_dsets[r].lon.isel(obs=0, drop=True) <= de_lon_max)
            & (regime_dsets[r].lat.isel(obs=0, drop=True) >= de_lat_min)
            & (regime_dsets[r].lat.isel(obs=0, drop=True) <= de_lat_max)
        ).mean("trajectory").expand_dims(regime=[r])
        for r in regimes
    ],
    dim="regime",
)

# Pre-compute the small (trajectory,)-shaped group coords once per regime.
regime_group_coords = {
    r: dict(
        subbasin=regime_dsets[r].subbasin.compute(),
        release_quarter=regime_dsets[r].release_quarter.compute(),
    )
    for r in regimes
}

da_sb_lazy = xr.concat(
    [
        regime_distance[r]
        .assign_coords(subbasin=regime_group_coords[r]["subbasin"])
        .groupby("subbasin").mean("trajectory")
        .expand_dims(regime=[r])
        for r in regimes
    ],
    dim="regime",
)

da_quarter_lazy = xr.concat(
    [
        regime_distance[r]
        .assign_coords(release_quarter=regime_group_coords[r]["release_quarter"])
        .groupby("release_quarter").mean("trajectory")
        .expand_dims(regime=[r])
        for r in regimes
    ],
    dim="regime",
)

da_global, da_de, da_sb, da_quarter = dask.compute(
    da_global_lazy, da_de_lazy, da_sb_lazy, da_quarter_lazy,
)
```

# Global

```python
fig, ax = plt.subplots(layout="constrained")
da_global.plot.line(x="obs", hue="regime", ax=ax)
```

# Per HELCOM release subbasin

```python
da_sb.plot.line(
    x="obs", hue="regime", col="subbasin", col_wrap=4,
    size=facet_line_size, aspect=facet_line_aspect,
)
```

# German waters (release cells inside bounding box)

```python
fig, ax = plt.subplots(layout="constrained")
da_de.plot.line(x="obs", hue="regime", ax=ax)
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
da_quarter.plot.line(
    x="obs", hue="regime", col="release_quarter", col_wrap=2,
    size=facet_line_size, aspect=facet_line_aspect,
)
```
