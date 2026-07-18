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

# Survival-weighted heatmaps: beaching folded into occupancy

Parquet-only consumer of the survival-occupancy store built by
`024e_BuildSurvivalOccupancy` (+ the `024a` key). For each elapsed-time
horizon it draws three panels:

1. **Occupancy** — plain particle-residence density (`occ`), the baseline
   (no beaching), on a log scale.
2. **Survival-weighted** — the same, with beaching progressively removed
   (`surv = Σ exp(−A)`), on the *same* log scale, so the thinning is legible.
3. **Surviving fraction** — `surv / occ` per hex (linear 0–1): where the
   free-drifting cloud is most depleted by stranding (near retentive shores).

One regime per run; `release_month = 0` pools every monthly partition
across years. No Dask, no zarrs — in the `025`/`026`/`029` lineage.

```python
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth
```

# Parameters

```python tags=["parameters"]
output_root = "../output"
regime = "surface_stokes"
hex_radius = 6000
age_bin_days = 10
# 0 = pool all monthly partitions across years; 1..12 = that month only.
release_month = 0

# Elapsed-time horizons to map (days); each a multiple of age_bin_days and
# within the store's occupancy_max_days.
time_horizons_days_csv = "20,50,100"

cmap = "viridis"
panel_height_in = 6
fig_dpi_scale = 2
```

# Parse parameters

```python
output_root = Path(output_root)
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d"
    )
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]

figure_dir = output_root / "Figures" / "030"
figure_dir.mkdir(parents=True, exist_ok=True)
```

# Read key + pool survival-occupancy partitions

```python
store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

month_suffix = f"_m{release_month:02d}" if release_month else ""
month_re = rf"_m{release_month:02d}" if release_month else r"_m\d{2}"
_PART_RE = re.compile(
    rf"HexAgg_survocc_r{hex_radius}m_{regime}_(\d{{4}}){month_re}\.parquet$"
)
survocc_files = [
    f for f in sorted(store_root.glob(f"HexAgg_survocc_r{hex_radius}m_{regime}_*.parquet"))
    if _PART_RE.search(f.name)
]
if not survocc_files:
    raise FileNotFoundError(
        f"no survocc partitions for regime {regime!r}"
        + (f", month {release_month}" if release_month else "")
        + f" at {store_root} — run 024e."
    )
survocc = pd.concat(
    [pd.read_parquet(f).reset_index(drop=True) for f in survocc_files],
    ignore_index=True,
)
print(f"key: {len(key):,} hexes; pooled {len(survocc_files)} (year, month) partition(s)")
print(f"survocc rows: {len(survocc):,}")
```

# Rendering helpers

```python
def hex_gdf(df, value_col):
    """Sum a value column per target_hex, joined to the key geometry."""
    return (
        df.groupby("target_hex")[value_col].sum().rename("value").reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def frac_gdf(df):
    """Surviving fraction surv/occ per target_hex, joined to key geometry."""
    g = df.groupby("target_hex")[["occ", "surv"]].sum()
    g["value"] = g["surv"] / g["occ"]
    return (
        g.reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def log_norm(values):
    """Shared LogNorm floored four decades below the peak (survival weights
    have a long sub-1 tail)."""
    vmax = float(values.max())
    vmin = max(float(values[values > 0].min()), vmax / 1e4)
    return LogNorm(vmin=vmin, vmax=vmax)


def hex_map(gdf, ax, norm=None, vmin=None, vmax=None, title=None):
    if not gdf.empty:
        gdf.plot(
            ax=ax, column="value", cmap=cmap, norm=norm, vmin=vmin, vmax=vmax,
            legend=True, edgecolor="face", linewidth=0.4, zorder=1,
        )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)
```

```python
lon_min, lat_min, lon_max, lat_max = key.total_bounds
extent = [lon_min, lon_max, lat_min, lat_max]
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)
```

# Comparison grid: occupancy vs survival-weighted vs surviving fraction

Rows = horizons; columns = occupancy, survival-weighted, surviving fraction.
Occupancy and survival share one `LogNorm` (across all horizons) so both the
age-thinning and the beaching-removal are directly readable.

```python
horizon_bins = {h: h // age_bin_days for h in time_horizons_days}
occ_gdfs = {h: hex_gdf(survocc[survocc["age_bin"] == b], "occ") for h, b in horizon_bins.items()}
surv_gdfs = {h: hex_gdf(survocc[survocc["age_bin"] == b], "surv") for h, b in horizon_bins.items()}
frac_gdfs = {h: frac_gdf(survocc[survocc["age_bin"] == b]) for h, b in horizon_bins.items()}

shared = pd.concat(
    [g["value"] for g in list(occ_gdfs.values()) + list(surv_gdfs.values()) if not g.empty]
)
dens_norm = log_norm(shared)

nrows = len(time_horizons_days)
fig, axes = plt.subplots(
    nrows=nrows, ncols=3,
    figsize=(panel_height_in * domain_aspect * 3, panel_height_in * nrows),
    layout="constrained", squeeze=False,
)
for i, h in enumerate(time_horizons_days):
    hex_map(occ_gdfs[h], axes[i, 0], norm=dens_norm, title=f"occupancy — {h} d")
    hex_map(surv_gdfs[h], axes[i, 1], norm=dens_norm, title=f"survival-weighted — {h} d")
    # Surviving fraction on a fixed linear 0–1 scale for cross-horizon reading.
    hex_map(frac_gdfs[h], axes[i, 2], vmin=0.0, vmax=1.0, title=f"surviving fraction — {h} d")
fig_path = figure_dir / f"SurvivalHeatmaps_{regime}_r{hex_radius}m{month_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()
```

# Validation / summary

```python
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["drifting_fraction"] = per_bin["surv"] / per_bin["occ"]
print(f"regime={regime}, hex_radius={hex_radius} m, "
      + (f"month={release_month}, " if release_month else "all months, ")
      + f"age_bin_days={age_bin_days}")
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
for h in time_horizons_days:
    b = horizon_bins[h]
    sub = survocc[survocc["age_bin"] == b]
    o, s = sub["occ"].sum(), sub["surv"].sum()
    print(f"  {h:>4} d (bin {b}): drifting fraction {s / max(o, 1):.3f}, "
          f"{sub['target_hex'].nunique():,} hexes")
```
