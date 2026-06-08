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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
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
def release_subbasin(lon0, lat0, subbasins):
    # Nearest-subbasin lookup on the obs=0 release positions. Lazy per
    # (trajectory,) chunk; STRtree built once, looked up per-chunk. lon0/lat0
    # arrive already sliced to obs=0 at the zarr read (see the load cell), so
    # the full (trajectory, obs) chunks never materialise — the full-sweep
    # concat reaches 60M+ trajectories and reading whole obs-chunks would OOM.
    tree = shapely.STRtree(subbasins.geometry.values)
    names = subbasins["subbasin"].to_numpy()

    def _lookup(lon, lat):
        out = np.full(lon.shape, None, dtype=object)
        valid = ~(np.isnan(lon) | np.isnan(lat))
        if valid.any():
            pts = shapely.points(lon[valid], lat[valid])
            out[valid] = names[tree.nearest(pts)]
        return out

    return xr.apply_ufunc(
        _lookup, lon0, lat0,
        dask="parallelized", output_dtypes=[object],
    )
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
    data_root / "helcom_fucus_redlist" / "REDLIST_SIS_Macrophytes.shp"
)
release_area = release_area.loc[
    release_area.F_vesiculo != 0, ["geometry", "CELLCODE"]
].to_crs(crs=ccrs.Geodetic())
release_area
```

# HELCOM subbasins

```python
subbasins = gpd.read_file(
    data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(crs=ccrs.Geodetic()).rename(dict(level_2="subbasin"), axis=1)
subbasins
```

# Load each regime

```python
regime_dsets = {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    # Cheap release-edge view: push the obs slice into each per-file read so
    # the full (trajectory, obs) chunks never linger as dask task results.
    # Two steps suffice — obs=0 for the release position, obs=1 for the land
    # test. A post-concat ds.isel(obs=0) instead reads and retains every full
    # obs-chunk, exhausting worker memory at the 60M+-trajectory scale.
    edge = xr.concat(
        [
            xr.open_zarr(z)[["lon", "lat", "time"]]
            .isel(obs=slice(0, 2))
            .chunk(obs=2)
            for z in zarr_files
        ],
        dim="trajectory",
    )
    lon0 = edge.lon.isel(obs=0, drop=True)
    lat0 = edge.lat.isel(obs=0, drop=True)
    # First-step displacement of zero ⇒ trajectory was seeded on land.
    on_land = (
        (edge.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (edge.lat.diff("obs").isel(obs=0, drop=True) == 0)
    )
    ds = ds.where(~on_land)
    # Release position as (trajectory,) coords — reused by distance_km and the
    # German-waters filter so neither re-slices obs=0 off the full concat.
    ds = ds.assign_coords(
        release_lon=lon0.where(~on_land),
        release_lat=lat0.where(~on_land),
    )
    ds = ds.assign(
        release_quarter=edge.time.isel(obs=0, drop=True).where(~on_land).dt.quarter,
        subbasin=release_subbasin(lon0.where(~on_land), lat0.where(~on_land), subbasins),
    )
    regime_dsets[regime] = ds
regime_dsets
```

# Distance from release

Great-circle approximation (111 km per degree lat).

```python
def distance_km(ds):
    # release_lon/release_lat are the (trajectory,) obs=0 coords attached at
    # load time, so this never re-reads obs=0 off the full concat.
    dlat = ds.lat - ds.release_lat
    dlon = (ds.lon - ds.release_lon) * np.cos(np.deg2rad(ds.release_lat))
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
            (regime_dsets[r].release_lon >= de_lon_min)
            & (regime_dsets[r].release_lon <= de_lon_max)
            & (regime_dsets[r].release_lat >= de_lat_min)
            & (regime_dsets[r].release_lat <= de_lat_max)
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
