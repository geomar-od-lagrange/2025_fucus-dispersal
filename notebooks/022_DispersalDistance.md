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

# Compute per-scope means (per-batch accumulation)

A `mean("trajectory")` keeping `obs` is a small-output reduction, but one
`dask.compute` over the whole-regime concat builds a single graph that
references every trajectory-chunk at once; the freed chunk buffers accumulate
as unmanaged worker memory and the run stalls at the pause threshold. Instead
we loop over small batches of files: each batch is its own bounded graph that
computes and releases before the next opens, so peak memory plateaus at
`batch_files` chunks regardless of how many files a regime has. A NaN-aware
mean is `sum / count`, and both are additive over any partition of the
trajectory axis, so summing per-batch partials and dividing at the end
reproduces the exact global mean per `obs` per scope.

`batch_files` trades plateau height for cluster parallelism: one file is a
single `(trajectory=10000, obs≈thousands)` chunk (one task, one worker), so a
batch holds `batch_files` such chunks (~1 GiB each for lon+lat). Keep it a
small multiple of the worker count so a pass spreads across the cluster while
staying well under a worker's memory.

`obs` length varies file to file — Parcels trims trailing empty steps — so the
per-batch partials have different `obs` extents. They are outer-aligned on
`obs` with fill `0` before accumulating (`add_obs`): a batch that ends earlier
contributes `0` to both sum and count at another batch's later `obs`, exactly
as those (absent) trajectories would have been NaN-skipped in a single concat.
A plain `+` would inner-join and silently truncate to the shortest batch.

```python
batch_files = 16

quarters = [1, 2, 3, 4]
# Sorted to match the order xarray's groupby would have produced.
all_subbasins = sorted(subbasins["subbasin"].tolist())


def scope_sum_count(dist, mask=None):
    # NaN-aware partials: skipna sum and count of non-NaN. Land-masked
    # trajectories are already NaN and drop out of both, matching .mean().
    d = dist if mask is None else dist.where(mask)
    return d.sum("trajectory"), d.notnull().sum("trajectory")


def add_obs(a, b):
    # Outer-align on the variable-length obs axis, fill 0, then add — so a
    # shorter batch contributes 0 to both sum and count at later obs.
    a, b = xr.align(a, b, join="outer", fill_value=0)
    return a + b


def nanmean(total_sum, total_count):
    # 0 valid trajectories ⇒ NaN, not 0/0.
    return (total_sum / total_count.where(total_count > 0)).rename("distance_km")


da_global, da_de, da_sb, da_quarter = {}, {}, {}, {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))

    acc = None  # dict of running (sum, count) partials, lazily initialised
    for start in range(0, len(zarr_files), batch_files):
        batch = zarr_files[start:start + batch_files]
        ds = xr.concat([xr.open_zarr(z) for z in batch], dim="trajectory")
        # Per-batch release edge (same obs pushdown as the metadata read):
        # obs=0 release position, obs=1 land test.
        edge = xr.concat(
            [
                xr.open_zarr(z)[["lon", "lat", "time"]]
                .isel(obs=slice(0, 2))
                .chunk(obs=2)
                for z in batch
            ],
            dim="trajectory",
        )
        lon0 = edge.lon.isel(obs=0, drop=True)
        lat0 = edge.lat.isel(obs=0, drop=True)
        keep = ~(
            (edge.lon.diff("obs").isel(obs=0, drop=True) == 0)
            & (edge.lat.diff("obs").isel(obs=0, drop=True) == 0)
        )
        ds = ds.assign_coords(
            release_lon=lon0.where(keep), release_lat=lat0.where(keep)
        )
        dist = distance_km(ds).where(keep)
        # Small (trajectory,) group labels — materialise once, reuse per group.
        subbasin = release_subbasin(lon0.where(keep), lat0.where(keep), subbasins).compute()
        quarter = edge.time.isel(obs=0, drop=True).where(keep).dt.quarter.compute()

        de_mask = (
            (ds.release_lon >= de_lon_min) & (ds.release_lon <= de_lon_max)
            & (ds.release_lat >= de_lat_min) & (ds.release_lat <= de_lat_max)
        )
        lazy = {"global": scope_sum_count(dist), "de": scope_sum_count(dist, de_mask)}
        for sb in all_subbasins:
            lazy[("sb", sb)] = scope_sum_count(
                dist, xr.DataArray(subbasin == sb, dims="trajectory")
            )
        for q in quarters:
            lazy[("q", q)] = scope_sum_count(
                dist, xr.DataArray(quarter == q, dims="trajectory")
            )
        # One dask pass per batch; every value is (obs,)-shaped and cheap.
        part = dask.compute(lazy)[0]
        acc = part if acc is None else {
            k: (add_obs(acc[k][0], part[k][0]), add_obs(acc[k][1], part[k][1]))
            for k in acc
        }

    da_global[regime] = nanmean(*acc["global"])
    da_de[regime] = nanmean(*acc["de"])
    da_sb[regime] = xr.concat(
        [nanmean(*acc[("sb", sb)]).expand_dims(subbasin=[sb]) for sb in all_subbasins],
        dim="subbasin",
    )
    da_quarter[regime] = xr.concat(
        [nanmean(*acc[("q", q)]).expand_dims(release_quarter=[q]) for q in quarters],
        dim="release_quarter",
    )

da_global = xr.concat([da_global[r].expand_dims(regime=[r]) for r in regimes], dim="regime")
da_de = xr.concat([da_de[r].expand_dims(regime=[r]) for r in regimes], dim="regime")
da_sb = xr.concat([da_sb[r].expand_dims(regime=[r]) for r in regimes], dim="regime")
da_quarter = xr.concat([da_quarter[r].expand_dims(regime=[r]) for r in regimes], dim="regime")

# Drop subbasins absent from every regime so the facet set matches the data
# (groupby would have omitted them); all four quarters always occur.
da_sb = da_sb.isel(subbasin=np.flatnonzero(da_sb.notnull().any(["regime", "obs"]).values))
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
