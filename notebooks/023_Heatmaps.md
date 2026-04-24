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
per release quarter (JFM/AMJ/JAS/OND).

```python
import dask
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from xhistogram.xarray import histogram as xhist
from pathlib import Path

from helpers import (
    attach_release_metadata,
    load_trajectories,
    mask_land_seeded,
    relabel_quarter,
)
```

# Parameters

```python tags=["parameters"]
base_path = "/gxfs_work/geomar/smomw122/2025_fucus-dispersal"
experiment_type = "surface"

output_dt_mins = 60

lon_min, lon_max = 5, 32
lat_min, lat_max = 53, 66
# Baltic-wide: 200 lon bins across the Baltic.
n_lon_baltic = 200

de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
dlon_de = 0.125
dlat_de = 0.125

baltic_panel_height_in = 2
de_panel_height_in = 1.5
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

# Bins and histogram helpers

```python
# Baltic-wide: exactly n_lon_baltic bins; same degree-resolution in lat.
lon_bins = np.linspace(lon_min, lon_max, n_lon_baltic + 1)
dlon_baltic = (lon_max - lon_min) / n_lon_baltic
n_lat_baltic = int(np.ceil((lat_max - lat_min) / dlon_baltic))
lat_bins = np.linspace(lat_min, lat_max, n_lat_baltic + 1)

de_lon_bins = np.arange(de_lon_min, de_lon_max + dlon_de, dlon_de)
de_lat_bins = np.arange(de_lat_min, de_lat_max + dlat_de, dlat_de)

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
    # Suppress exactly-zero cells so land renders as figure background,
    # not the lowest colormap colour.
    d = h.sum("obs")
    return d.where(d > 0)

def mean_age_hours(h):
    totals = h.sum("obs")
    return (((h * h.obs).sum("obs") / totals) * age_hours_per_obs).where(totals > 0)
```

```python
def lonlat_aspect(extent):
    """Displayed width / height ratio for a lon/lat ``extent`` with an
    aspect that keeps 1 deg lon at ``lat_mean`` visually equal to 1 deg
    lat."""
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    lat_mean = 0.5 * (lat_min_ + lat_max_)
    return ((lon_max_ - lon_min_) * np.cos(np.radians(lat_mean))) / (
        lat_max_ - lat_min_
    )

baltic_extent = [lon_min, lon_max, lat_min, lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
baltic_aspect = lonlat_aspect(baltic_extent)
de_aspect = lonlat_aspect(de_extent)
```

```python
def facet_map(da, col, col_wrap=None):
    fg = da.plot(
        x="lon", y="lat",
        col=col, col_wrap=col_wrap,
        size=baltic_panel_height_in, aspect=baltic_aspect,
        subplot_kws=dict(projection=ccrs.PlateCarree()),
        transform=ccrs.PlateCarree(),
    )
    for ax in fg.axs.flat:
        ax.coastlines()
    return fg

def single_map(da, extent, panel_height_in):
    width = lonlat_aspect(extent) * panel_height_in
    fig, ax = plt.subplots(
        figsize=(width, panel_height_in),
        subplot_kw=dict(projection=ccrs.PlateCarree()),
        layout="constrained",
    )
    da.plot(ax=ax, x="lon", y="lat", transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    return fig
```

# Precompute per-trajectory scope keys

`subbasin`, `release_quarter` are lazy 1-D `(trajectory,)` arrays.
Materialise their unique values once so the per-scope histograms below
can be built as a single lazy graph.

```python
sb_np, quarter_np = dask.compute(ds.subbasin, ds.release_quarter)
sb_np = sb_np.values
quarter_np = quarter_np.values

subbasins_list = sorted({s for s in sb_np if isinstance(s, str)})
quarters = sorted({int(q) for q in quarter_np if not np.isnan(q)})
```

# Compute histograms (one shared dask pass)

Histograms live in bin space (small); compute them all in one
`dask.compute(*)` so the trajectory graph is only walked once. Per-scope
histograms loop over group values with a lazy `.where()` mask — the
subbasin/year/quarter coords stay dask-backed and chunk-aligned with
`ds.lon`/`ds.lat`, so the mask fuses block-wise on the workers instead
of forcing client-side fancy indexing along the concat dim.

```python
def hist_by(key_da, values, name, lon_bins, lat_bins):
    parts = [
        count_hist(ds.where(key_da == v), lon_bins, lat_bins) for v in values
    ]
    return xr.concat(parts, dim=name).assign_coords({name: values})

h_baltic_lazy = count_hist(ds, lon_bins, lat_bins)
h_de_lazy = count_hist(ds, de_lon_bins, de_lat_bins)
h_by_sb_lazy = hist_by(ds.subbasin, subbasins_list, "subbasin", lon_bins, lat_bins)
h_by_quarter_lazy = hist_by(ds.release_quarter, quarters, "release_quarter", lon_bins, lat_bins)

h_baltic, h_de, h_by_sb, h_by_quarter = dask.compute(
    h_baltic_lazy, h_de_lazy, h_by_sb_lazy, h_by_quarter_lazy,
)
h_by_quarter = relabel_quarter(h_by_quarter)
```

# Whole Baltic

```python
single_map(density(h_baltic), baltic_extent, baltic_panel_height_in)
plt.show()
single_map(mean_age_hours(h_baltic), baltic_extent, baltic_panel_height_in)
plt.show()
```

# Per HELCOM release subbasin

```python
facet_map(density(h_by_sb), col="subbasin", col_wrap=4)
plt.show()
facet_map(mean_age_hours(h_by_sb), col="subbasin", col_wrap=4)
plt.show()
```

# German waters

```python
single_map(density(h_de), de_extent, de_panel_height_in)
plt.show()
single_map(mean_age_hours(h_de), de_extent, de_panel_height_in)
plt.show()
```

# Per release quarter (JFM/AMJ/JAS/OND)

```python
facet_map(density(h_by_quarter), col="release_quarter")
plt.show()
facet_map(mean_age_hours(h_by_quarter), col="release_quarter")
plt.show()
```
