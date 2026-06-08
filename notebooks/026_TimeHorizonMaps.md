---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: py:percent,md,ipynb
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

# Time-horizon dispersal maps

Hex particle-density maps at successive elapsed-time horizons (e.g. 20,
40, 80, 160 days), showing how the dispersal cloud spreads over time.
Reads the hex counts store built by 024a (key) + 024 (counts) — no
trajectory zarrs, no Dask cluster. One regime per run; August/September
releases pooled across all available years (empty month list ⇒ all
releases).

Each horizon selects the matching 10-day `age_bin` of the counts store
(the store already carries elapsed-time bins), so this is 025's
per-quarter panel re-keyed by elapsed time.

```python
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth
```

# Parameters

```python tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Which regime's counts partitions to read; one regime per run.
regime = "surface"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the counts store (must match 024). Horizons must
# be whole multiples of this.
age_bin_days = 10

# Release months to keep (Aug, Sep), pooled across all available years.
# Empty ⇒ all releases.
release_months_csv = "8,9"
# Elapsed-time horizons to map, in days. Each must be a multiple of
# age_bin_days.
time_horizons_days_csv = "20,40,80,160"

# Baltic-wide map extent (degrees E / degrees N).
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66

# Colormap: log-density spans several decades, so a perceptually uniform
# map is load-bearing (justified in docs/visualisations.md, as for 025).
cmap = "viridis"

# Per-panel height in inches (panel widths are aspect-derived).
baltic_panel_height_in = 6
```

# Parse parameters

```python
output_root = Path(output_root)
release_months = [int(x) for x in release_months_csv.split(",") if x]
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would fall in a neighbouring age bin"
    )
```

# Read key + pool counts across years

Layout: flat files under ``output_root/HexAggregates/`` —
``HexAgg_key_r<radius>m.parquet`` and
``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``. Globbing every year
pools across years; the release year is parsed from each filename so the
`release_doy → month` conversion is leap-correct.

```python
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)

# Match exactly four year digits after the regime so `surface` does not also
# glob `surface_stokes_*` files (prefix collision).
_COUNTS_RE = re.compile(rf"HexAgg_counts_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$")
counts_files = sorted(
    store_root.glob(f"HexAgg_counts_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet")
)
if not counts_files:
    raise FileNotFoundError(
        f"no counts partitions for regime {regime!r} at {store_root} — run 024."
    )

parts = []
for f in counts_files:
    year = int(_COUNTS_RE.search(f.name).group(1))
    df = pd.read_parquet(f).reset_index(drop=True)
    df["release_month"] = pd.to_datetime(
        (year * 1000 + df["release_doy"].astype("int32")).astype(str), format="%Y%j"
    ).dt.month
    parts.append(df)
counts = pd.concat(parts, ignore_index=True)
print(f"key: {len(key):,} hexes; pooled {len(counts_files)} year partition(s)")

if release_months:
    counts = counts[counts["release_month"].isin(release_months)]
    print(f"Keeping release months {release_months}: {len(counts):,} rows")
else:
    print(f"All release months: {len(counts):,} rows")
```

# Rendering helpers

Copied from 025 (`to_hex_gdf` joins a per-`target_hex` count table to the
key geometry; `log_density_plot` draws hex polygons + coastline on a
plain EPSG:4326 axis). Subbasin overlay dropped — not used here.

```python
def to_hex_gdf(counts_df, key_df):
    """Build a hex GeoDataFrame from per-target_hex counts + a key
    GeoDataFrame. Sums ``n_obs`` over ``target_hex`` and joins to
    ``key_df`` on ``target_hex`` ↔ ``hex_id``."""
    return (
        counts_df.groupby("target_hex")["n_obs"].sum()
        .rename("count").reset_index()
        .merge(key_df[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def log_density_plot(gdf, ax, extent, title=None, coast=None):
    """Plot hex density + coastline on a plain (non-cartopy) lon/lat axis.
    Everything stays in EPSG:4326 so hexes and coastline share one
    coordinate system. Empty ``gdf`` skips the hex layer so empty panels
    render with stable layout."""
    if not gdf.empty:
        gdf = gdf.copy()
        gdf["log_count"] = np.log10(gdf["count"].where(gdf["count"] > 0))
        # cmap/edgecolor/linewidth justified in docs/visualisations.md (as 025).
        gdf.plot(
            ax=ax, column="log_count", cmap=cmap, legend=False,
            missing_kwds={"color": "none"}, edgecolor="face", linewidth=0.4, zorder=1,
        )
    if coast is not None:
        coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    ax.set_xlim(lon_min_, lon_max_)
    ax.set_ylim(lat_min_, lat_max_)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (lat_min_ + lat_max_))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)
```

```python
# Natural Earth 10m coastline via cartopy's shapereader (cached locally),
# clipped to the Baltic extent once.
_coast_gdf = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
)
coast_baltic = _coast_gdf.clip(box(baltic_lon_min, baltic_lat_min, baltic_lon_max, baltic_lat_max))

baltic_extent = [baltic_lon_min, baltic_lon_max, baltic_lat_min, baltic_lat_max]
# Aspect ratio that keeps 1° lon at lat_mean visually equal to 1° lat.
baltic_aspect = (
    (baltic_lon_max - baltic_lon_min) * np.cos(np.radians(0.5 * (baltic_lat_min + baltic_lat_max)))
) / (baltic_lat_max - baltic_lat_min)
```

# Map per horizon

Each horizon selects `age_bin = horizon // age_bin_days` — the 10-day
occupancy window `[horizon, horizon + age_bin_days)`.

```python
horizon_age_bins = {h: h // age_bin_days for h in time_horizons_days}

ncols = 2
nrows = int(np.ceil(len(time_horizons_days) / ncols))
fig, axes = plt.subplots(
    nrows=nrows, ncols=ncols,
    figsize=(
        baltic_panel_height_in * baltic_aspect * ncols,
        baltic_panel_height_in * nrows,
    ),
    layout="constrained",
)
for ax, h in zip(axes.flat, time_horizons_days):
    counts_h = counts[counts["age_bin"] == horizon_age_bins[h]]
    gdf_h = to_hex_gdf(counts_h, key)
    log_density_plot(gdf_h, ax, baltic_extent, title=f"{h} d", coast=coast_baltic)
for ax in axes.flat[len(time_horizons_days):]:
    ax.set_visible(False)
plt.show()
```

# Validation prints

```python
print(f"regime={regime}, hex_radius={hex_radius} m, age_bin_days={age_bin_days}")
print(f"release months: {release_months or 'all'}")
for h in time_horizons_days:
    counts_h = counts[counts["age_bin"] == horizon_age_bins[h]]
    print(
        f"  {h:>4}d (age_bin {horizon_age_bins[h]}): "
        f"{counts_h['target_hex'].nunique():,} target hexes, "
        f"sum(n_obs)={int(counts_h['n_obs'].sum()):,}"
    )
```
