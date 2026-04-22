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
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
from pathlib import Path

# Silence two noisy-but-harmless warning classes the structural fixes
# should already avoid, as belt-and-braces for edge cases.
warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in linestrings",
    category=RuntimeWarning,
)
warnings.filterwarnings(
    "ignore",
    message=r"Sending large graph of size",
    category=UserWarning,
)

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

n_traj_subset = 300

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5

# Panel sizing probed for ~3" width per panel at standard dpi, Baltic box
# (27x13 deg) in PlateCarree. Single German panel is narrower in lat
# (2.3 deg), so figsize height is reduced to match.
single_baltic_figsize = (4, 2)
single_de_figsize = (4, 1.5)
subbasin_grid_figsize = (18, 9)   # 4x4, Baltic extent per panel
quarter_grid_figsize = (9, 4.5)   # 2x2, Baltic extent per panel

# Regime colours.
regime_colors = {
    "bottom": "tab:orange",
    "surface": "tab:blue",
    "surface_stokes": "tab:green",
}
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

# Load all regimes and attach metadata

```python
trajectory_root = base_path / "output" / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")

regime_dsets = {}
for regime in regimes:
    ds, zarr_files = load_trajectories(trajectory_root / regime)
    print(f"{regime}: {len(zarr_files)} trajectory files")
    ds, _ = mask_land_seeded(ds)
    ds = attach_release_metadata(ds, subbasins)
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

`sample_subsets` picks up to `n` random trajectories per group and
materialises them in a single `dask.compute` per panel loop — one big
`ds.isel(trajectory=...)` instead of N per-group isels.

`plot_lines` renders a per-trajectory NaN-filtered `LineCollection` so
cartopy's `path_to_shapely` never hands NaN coords to
`shapely.linestrings` (the source of the
``invalid value encountered in linestrings`` warning storm). NaN-only
or <2-valid-point trajectories are dropped.

```python
def sample_subsets(ds_, groups, n):
    """groups: dict name -> bool sel over trajectory (None = all). Returns
    dict name -> locally-materialised Dataset with <= n trajectories."""
    picks = {}
    parts = []
    for name, sel in groups.items():
        avail = (
            np.arange(ds_.sizes["trajectory"])
            if sel is None
            else np.flatnonzero(sel)
        )
        if avail.size == 0:
            picks[name] = np.array([], dtype=int)
            continue
        chosen = rng.choice(avail, min(avail.size, n), replace=False)
        picks[name] = chosen
        parts.append(chosen)
    if not parts:
        empty = ds_.isel(trajectory=[]).compute()
        return {name: empty for name in groups}
    flat = np.concatenate(parts)
    ds_all = ds_.isel(trajectory=flat).compute()
    out = {}
    offset = 0
    for name, chosen in picks.items():
        k = chosen.size
        out[name] = ds_all.isel(trajectory=slice(offset, offset + k))
        offset += k
    return out


def plot_lines(ds_plot, ax, color, lw=None):
    if ds_plot.sizes["trajectory"] == 0:
        return
    lon = ds_plot.lon.values
    lat = ds_plot.lat.values
    segments = []
    for i in range(lon.shape[0]):
        xi, yi = lon[i], lat[i]
        m = ~(np.isnan(xi) | np.isnan(yi))
        if m.sum() < 2:
            continue
        segments.append(np.column_stack([xi[m], yi[m]]))
    if not segments:
        return
    ax.add_collection(
        LineCollection(
            segments, linewidths=lw, colors=color,
            alpha=0.3, transform=ccrs.PlateCarree(),
        )
    )


def regime_legend_handles():
    return [
        Line2D([], [], color=regime_colors[r], label=r, linewidth=1.5)
        for r in regimes
    ]


def sample_per_regime(groups_for):
    """For each regime, run sample_subsets(ds, groups_for(regime), n).
    Returns dict[regime] -> dict[group_name -> materialised ds]."""
    return {
        regime: sample_subsets(regime_dsets[regime], groups_for(regime), n_traj_subset)
        for regime in regimes
    }
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

subsets_per_regime = sample_per_regime(
    lambda regime: {b: regime_keys[regime]["subbasin"] == b for b in subbasins_list}
)

ncols = 4
nrows = int(np.ceil(len(subbasins_list) / ncols))
for regime in regimes:
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols,
        figsize=subbasin_grid_figsize,
        subplot_kw=dict(projection=ccrs.PlateCarree()),
    )
    for ax, basin in zip(axes.flat, subbasins_list):
        plot_lines(
            subsets_per_regime[regime][basin], ax,
            color=regime_colors[regime], lw=0.5,
        )
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.coastlines()
        ax.set_title(basin)
    for ax in axes.flat[len(subbasins_list):]:
        ax.set_visible(False)
    fig.suptitle(regime)
    fig.legend(handles=regime_legend_handles(), loc="lower center", ncol=len(regimes))
    plt.show()
```

# German waters

```python
subsets_per_regime = sample_per_regime(lambda regime: {"all": None})

fig, axes = plt.subplots(
    nrows=1, ncols=len(regimes),
    figsize=(single_de_figsize[0] * len(regimes), single_de_figsize[1]),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for ax, regime in zip(axes, regimes):
    plot_lines(
        subsets_per_regime[regime]["all"], ax,
        color=regime_colors[regime],
    )
    ax.set_extent([de_lon_min, de_lon_max, de_lat_min, de_lat_max], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_title(regime)
fig.legend(handles=regime_legend_handles(), loc="lower center", ncol=len(regimes))
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
subsets_per_regime = sample_per_regime(
    lambda regime: {q_int: regime_keys[regime]["quarter"] == q_int for q_int in QUARTER_LABELS}
)

nrows = len(QUARTER_LABELS)
ncols = len(regimes)
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(single_baltic_figsize[0] * ncols, single_baltic_figsize[1] * nrows),
    subplot_kw=dict(projection=ccrs.PlateCarree()),
)
for row, (q_int, q_label) in enumerate(QUARTER_LABELS.items()):
    for col, regime in enumerate(regimes):
        ax = axes[row, col]
        plot_lines(
            subsets_per_regime[regime][q_int], ax,
            color=regime_colors[regime],
        )
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.coastlines()
        if row == 0:
            ax.set_title(regime)
        if col == 0:
            ax.set_ylabel(q_label)
fig.legend(handles=regime_legend_handles(), loc="lower center", ncol=len(regimes))
plt.show()
```
