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
import os
import time

import dask
import numpy as np
import xarray as xr
import geopandas as gpd
import shapely
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from xhistogram.xarray import histogram as xhist
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
# Which regime's trajectory zarrs to read; one regime per run.
regime = "surface"

# Trajectory output cadence (minutes per obs); used for age conversion.
output_dt_mins = 60

# Baltic-wide map extent (degrees E / degrees N) and bin width (degrees).
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66
dlon_baltic = 0.135
dlat_baltic = 0.135

# German-waters zoom map extent and bin width.
de_lon_min, de_lon_max = 8, 15
de_lat_min, de_lat_max = 53.2, 55.5
dlon_de = 0.125
dlat_de = 0.125

# Per-panel height in inches (panel widths are aspect-derived).
baltic_panel_height_in = 2
de_panel_height_in = 1.5
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

# Load trajectories and attach metadata

Layout assumption: ``output_root/Trajectories/<regime>/<release_year>/*.zarr``.

```python
trajectory_path = output_root / "Trajectories" / regime
zarr_files = sorted(trajectory_path.glob("**/*.zarr"))
print(f"{len(zarr_files)} trajectory files for {regime}")
ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
# First-step displacement of zero ⇒ trajectory was seeded on land.
ds = ds.where(~(
    (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
    & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
))
ds = ds.assign(release_quarter=ds.time.isel(obs=0, drop=True).dt.quarter)
ds = assign_release_subbasin(ds, subbasins)
ds
```

# Bins and extents

```python
lon_bins = np.arange(baltic_lon_min, baltic_lon_max + dlon_baltic, dlon_baltic)
lat_bins = np.arange(baltic_lat_min, baltic_lat_max + dlat_baltic, dlat_baltic)
de_lon_bins = np.arange(de_lon_min, de_lon_max + dlon_de, dlon_de)
de_lat_bins = np.arange(de_lat_min, de_lat_max + dlat_de, dlat_de)

age_hours_per_obs = output_dt_mins / 60

baltic_extent = [baltic_lon_min, baltic_lon_max, baltic_lat_min, baltic_lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]
# Aspect ratio that keeps 1° lon at lat_mean visually equal to 1° lat.
baltic_aspect = (
    (baltic_lon_max - baltic_lon_min) * np.cos(np.radians(0.5 * (baltic_lat_min + baltic_lat_max)))
) / (baltic_lat_max - baltic_lat_min)
de_aspect = (
    (de_lon_max - de_lon_min) * np.cos(np.radians(0.5 * (de_lat_min + de_lat_max)))
) / (de_lat_max - de_lat_min)
```

# Histogram reductions

```python
def density(h):
    # Suppress exactly-zero cells so land renders as figure background,
    # not the lowest colormap colour.
    d = h.sum("obs")
    return d.where(d > 0)

def mean_age_hours(h):
    totals = h.sum("obs")
    return (((h * h.obs).sum("obs") / totals) * age_hours_per_obs).where(totals > 0)
```

# Plot helpers

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
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    aspect = ((lon_max_ - lon_min_) * np.cos(np.radians(0.5 * (lat_min_ + lat_max_)))) / (
        lat_max_ - lat_min_
    )
    fig, ax = plt.subplots(
        figsize=(aspect * panel_height_in, panel_height_in),
        subplot_kw=dict(projection=ccrs.PlateCarree()),
        layout="constrained",
    )
    da.plot(ax=ax, x="lon", y="lat", transform=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    return fig
```

# Compute histograms (one shared dask pass)

Histograms live in bin space (small); compute them all in one
`dask.compute(*)` so the trajectory graph is only walked once. Per-scope
histograms loop over group values with a lazy `.where()` mask — the
subbasin/release_quarter coords stay dask-backed and chunk-aligned with
`ds.lon`/`ds.lat`, so the mask fuses block-wise on the workers instead
of forcing client-side fancy indexing along the concat dim.

Subbasins come from the HELCOM polygons; quarters are literally
`[1, 2, 3, 4]`. No pre-pass over the trajectory data is needed.

```python
def hist_by(key_da, values, name, lon_bins_, lat_bins_):
    parts = [
        xhist(
            ds.lon.where(key_da == v), ds.lat.where(key_da == v),
            bins=[lon_bins_, lat_bins_], dim=["trajectory"],
        ).rename(dict(lon_bin="lon", lat_bin="lat"))
        for v in values
    ]
    return xr.concat(parts, dim=name).assign_coords({name: values})

subbasins_list = subbasins["subbasin"].tolist()
quarters = [1, 2, 3, 4]

h_baltic_lazy = xhist(
    ds.lon, ds.lat, bins=[lon_bins, lat_bins], dim=["trajectory"],
).rename(dict(lon_bin="lon", lat_bin="lat"))
h_de_lazy = xhist(
    ds.lon, ds.lat, bins=[de_lon_bins, de_lat_bins], dim=["trajectory"],
).rename(dict(lon_bin="lon", lat_bin="lat"))
h_by_sb_lazy = hist_by(ds.subbasin, subbasins_list, "subbasin", lon_bins, lat_bins)
h_by_quarter_lazy = hist_by(ds.release_quarter, quarters, "release_quarter", lon_bins, lat_bins)

h_baltic, h_de, h_by_sb, h_by_quarter = dask.compute(
    h_baltic_lazy, h_de_lazy, h_by_sb_lazy, h_by_quarter_lazy,
)
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
