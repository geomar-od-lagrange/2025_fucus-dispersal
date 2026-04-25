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

Polyline subsets of particle tracks on maps. All regimes overlaid per
panel, one colour per regime. Scopes: per HELCOM release subbasin,
German waters, per release quarter (JFM/AMJ/JAS/OND).

```python
import warnings

import dask
import numpy as np
import xarray as xr
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
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
data_root = "../data"
output_root = "../output"

n_traj_subset = 300

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5

baltic_panel_height_in = 2
de_panel_height_in = 1.5

# Regime colours.
regime_colors = {
    "bottom": "tab:orange",
    "surface": "tab:blue",
    "surface_stokes": "tab:green",
}
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

# Release area

```python
data_root = Path(data_root)
output_root = Path(output_root)

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

# Load all regimes and attach metadata

```python
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")

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

# Precompute per-trajectory scope keys per regime

`release_quarter` / `subbasin` are lazy 1-D `(trajectory,)` arrays.
Compute them once per regime so the per-panel loops below don't re-walk
the graph.

```python
regime_keys = {}
for regime, ds in regime_dsets.items():
    quarter_np, subbasin_np = dask.compute(ds.release_quarter, ds.subbasin)
    regime_keys[regime] = dict(
        quarter=quarter_np.values,
        subbasin=subbasin_np.values,
    )
```

# Plot helpers

```python
def lonlat_aspect(extent):
    """Displayed width / height ratio for a lon/lat ``extent`` with an
    aspect that keeps 1 deg lon at ``lat_mean`` visually equal to 1 deg
    lat. Matches ``ax.set_aspect(1 / cos(lat_mean))`` below."""
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    lat_mean = 0.5 * (lat_min_ + lat_max_)
    return ((lon_max_ - lon_min_) * np.cos(np.radians(lat_mean))) / (
        lat_max_ - lat_min_
    )


baltic_extent = [lon_min, lon_max, lat_min, lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
baltic_aspect = lonlat_aspect(baltic_extent)
de_aspect = lonlat_aspect(de_extent)


def plot_lines(ds_plot, ax, color, lw=None):
    if ds_plot.sizes["trajectory"] == 0:
        return
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for i in range(ds_plot.sizes["trajectory"]):
            ax.plot(
                ds_plot.lon.isel(trajectory=i),
                ds_plot.lat.isel(trajectory=i),
                color=color, linewidth=lw, alpha=0.3,
                transform=ccrs.PlateCarree(),
            )
```

# Per HELCOM release subbasin

Subbasin list comes from the union over regimes (they should agree, but
the union is robust to a regime missing a subbasin entirely).

```python
subbasins_list = sorted({
    s
    for regime in regimes
    for s in regime_keys[regime]["subbasin"]
    if isinstance(s, str)
})

ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
for regime in regimes:
    ds = regime_dsets[regime]
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=(baltic_panel_height_in * baltic_aspect * ncols, baltic_panel_height_in * nrows),
        layout="constrained",
        subplot_kw=dict(projection=ccrs.PlateCarree()),
    )
    for ax, basin in zip(axes.flat, subbasins_list):
        mask = regime_keys[regime]["subbasin"] == basin
        avail = np.flatnonzero(mask)
        if avail.size == 0:
            ax.set_visible(False)
            continue
        idx = rng.choice(avail, size=min(n_traj_subset, avail.size), replace=False)
        ds_plot = ds.isel(trajectory=idx).compute()
        plot_lines(ds_plot, ax, color=regime_colors[regime], lw=0.5)
        ax.set_extent(baltic_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(basin)
    for ax in axes.flat[len(subbasins_list):]:
        ax.set_visible(False)
    fig.suptitle(regime)
    fig.legend(
        handles=[Line2D([], [], color=regime_colors[r], label=r, linewidth=1.5) for r in regimes],
        loc="lower center", ncol=len(regimes),
    )
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
    plot_lines(ds_plot, ax, color=regime_colors[regime])
    ax.set_extent(de_extent, crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(regime)
fig.legend(
    handles=[Line2D([], [], color=regime_colors[r], label=r, linewidth=1.5) for r in regimes],
    loc="lower center", ncol=len(regimes),
)
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
quarter_labels = {1: "JFM", 2: "AMJ", 3: "JAS", 4: "OND"}
nrows = len(quarter_labels)
ncols = len(regimes)
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
        mask = regime_keys[regime]["quarter"] == q_int
        avail = np.flatnonzero(mask)
        if avail.size == 0:
            ax.set_visible(False)
            continue
        idx = rng.choice(avail, size=min(n_traj_subset, avail.size), replace=False)
        ds_plot = ds.isel(trajectory=idx).compute()
        plot_lines(ds_plot, ax, color=regime_colors[regime])
        ax.set_extent(baltic_extent, crs=ccrs.PlateCarree())
        ax.coastlines()
        if row == 0:
            ax.set_title(regime)
        if col == 0:
            ax.set_ylabel(q_label)
fig.legend(
    handles=[Line2D([], [], color=regime_colors[r], label=r, linewidth=1.5) for r in regimes],
    loc="lower center", ncol=len(regimes),
)
plt.show()
```
