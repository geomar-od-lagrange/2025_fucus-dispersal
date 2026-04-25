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

# Raw trajectories

Polyline subsets of particle tracks on maps. All regimes overlaid per
panel, one colour per regime. Scopes: per HELCOM release subbasin,
German waters, per release quarter (JFM/AMJ/JAS/OND).

```python
import os
import time
import warnings

import numpy as np
import xarray as xr
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from dask.distributed import Client
from pathlib import Path

# RuntimeWarnings come from per-trajectory NaN tails (chunks that
# didn't fill up); the polyline plot tolerates them.
warnings.simplefilter("ignore", RuntimeWarning)
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
# Read root of trajectory zarrs (output_root/Trajectories/<regime>/...).
output_root = "../output"

# Number of trajectories sampled per panel (subbasin or quarter).
n_traj_subset = 300

# Baltic-wide map extent (degrees E / degrees N).
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66
# German-waters zoom map extent (degrees E / degrees N).
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5

# Per-panel height in inches (panel widths are aspect-derived).
baltic_panel_height_in = 2
de_panel_height_in = 1.5
```

# Global RNG

```python
seed = np.random.randint(0, 2**31 - 1)
print(f"RNG seed: {seed}")
rng = np.random.default_rng(seed)
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

# Release area

```python
data_root = Path(data_root)
output_root = Path(output_root)

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

# Load all regimes and attach metadata

Layout assumption: ``output_root/Trajectories/<regime>/<release_year>/*.zarr``.
``regimes`` is the list of immediate subdirectories.

```python
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")

# One colour per regime, taken from the matplotlib default cycle by index.
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
regime_color = {r: cycle[i] for i, r in enumerate(regimes)}

regime_dsets = {}
for regime in regimes:
    zarr_files = sorted((trajectory_root / regime).glob("**/*.zarr"))
    print(f"{regime}: {len(zarr_files)} trajectory files")
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    # First-step displacement of zero ⇒ trajectory was seeded on land.
    ds = ds.where(~(
        (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
    ))
    ds = ds.assign(release_quarter=ds.time.isel(obs=0, drop=True).dt.quarter)
    ds = assign_release_subbasin(ds, subbasins)
    regime_dsets[regime] = ds
```

# Map extents and aspects

```python
baltic_extent = [baltic_lon_min, baltic_lon_max, baltic_lat_min, baltic_lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
# Aspect ratio that keeps 1° lon at lat_mean visually equal to 1° lat
# (matches ax.set_aspect(1 / cos(lat_mean))).
baltic_aspect = (
    (baltic_lon_max - baltic_lon_min) * np.cos(np.radians(0.5 * (baltic_lat_min + baltic_lat_max)))
) / (baltic_lat_max - baltic_lat_min)
de_aspect = (
    (de_lon_max - de_lon_min) * np.cos(np.radians(0.5 * (de_lat_min + de_lat_max)))
) / (de_lat_max - de_lat_min)
```

# Per HELCOM release subbasin

One figure per regime, one panel per HELCOM subbasin. Empty panels
(subbasin with no release in this regime) are kept blank to keep the
layout stable across regimes.

```python
ncols = 4
nrows = int(np.ceil(len(subbasins) / ncols))
for regime in regimes:
    ds = regime_dsets[regime]
    # subbasin is a (trajectory,)-shaped coord — small, materialise once
    # per regime so the per-panel filter is a numpy comparison.
    subbasin_np = ds.subbasin.compute().values
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=(baltic_panel_height_in * baltic_aspect * ncols, baltic_panel_height_in * nrows),
        layout="constrained",
        subplot_kw=dict(projection=ccrs.PlateCarree()),
    )
    for ax, basin in zip(axes.flat, subbasins["subbasin"].tolist()):
        avail = np.flatnonzero(subbasin_np == basin)
        idx = rng.choice(avail, size=min(n_traj_subset, avail.size), replace=False)
        ds_plot = ds.isel(trajectory=idx).compute()
        for i in range(ds_plot.sizes["trajectory"]):
            # TODO(phase-f): justify linewidth in docs/visualisations.md.
            ax.plot(
                ds_plot.lon.isel(trajectory=i),
                ds_plot.lat.isel(trajectory=i),
                color=regime_color[regime], linewidth=0.5, alpha=0.3,
                transform=ccrs.PlateCarree(),
            )
        ax.set_extent(baltic_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(basin)
    for ax in axes.flat[len(subbasins):]:
        ax.set_visible(False)
    fig.suptitle(regime)
    plt.show()
```

# German waters

```python
fig, axes = plt.subplots(
    nrows=1, ncols=len(regimes),
    figsize=(de_panel_height_in * de_aspect * len(regimes), de_panel_height_in),
    layout="constrained",
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, regime in zip(axes, regimes):
    ds = regime_dsets[regime]
    idx = rng.choice(ds.sizes["trajectory"], size=min(n_traj_subset, ds.sizes["trajectory"]), replace=False)
    ds_plot = ds.isel(trajectory=idx).compute()
    for i in range(ds_plot.sizes["trajectory"]):
        ax.plot(
            ds_plot.lon.isel(trajectory=i),
            ds_plot.lat.isel(trajectory=i),
            color=regime_color[regime], alpha=0.3,
            transform=ccrs.PlateCarree(),
        )
    ax.set_extent(de_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(regime)
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
quarter_labels = {1: "JFM", 2: "AMJ", 3: "JAS", 4: "OND"}
nrows = len(quarter_labels)
ncols = len(regimes)
regime_quarter_np = {
    regime: regime_dsets[regime].release_quarter.compute().values for regime in regimes
}
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(baltic_panel_height_in * baltic_aspect * ncols, baltic_panel_height_in * nrows),
    layout="constrained",
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for row, (q_int, q_label) in enumerate(quarter_labels.items()):
    for col, regime in enumerate(regimes):
        ax = axes[row, col]
        ds = regime_dsets[regime]
        avail = np.flatnonzero(regime_quarter_np[regime] == q_int)
        idx = rng.choice(avail, size=min(n_traj_subset, avail.size), replace=False)
        ds_plot = ds.isel(trajectory=idx).compute()
        for i in range(ds_plot.sizes["trajectory"]):
            ax.plot(
                ds_plot.lon.isel(trajectory=i),
                ds_plot.lat.isel(trajectory=i),
                color=regime_color[regime], alpha=0.3,
                transform=ccrs.PlateCarree(),
            )
        ax.set_extent(baltic_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        if row == 0:
            ax.set_title(regime)
        if col == 0:
            ax.set_ylabel(q_label)
plt.show()
```
