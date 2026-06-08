# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,md,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # hex0 distance-quantile maps
#
# Per-source-hex crow-flies final-displacement quantiles, drawn on the
# 024a hex grid as one map per quantile level. Reads the distance
# histogram store built by 024b (+ the 024a key for geometry) — no
# trajectory zarrs, no Dask cluster. One regime per run; August/September
# releases pooled across all available years (empty month list ⇒ all
# releases).
#
# Quantiles are derived from the pooled per-hex histogram (cumulative count
# over `distance_bin`). Histograms are additive, so pooling across years is
# summing partitions — unlike pre-computed quantiles, which cannot be
# averaged across years.

# %%
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Which regime's distance partitions to read; one regime per run.
regime = "surface"
# Hex radius of the store (built by 024a/024b). Must match files on disk.
hex_radius = 6000
# Distance histogram bin width (km). Must match the 024b build.
distance_bin_km = 1.0

# Release months to keep (Aug, Sep), pooled across all available years.
# Empty ⇒ all releases.
release_months_csv = "8,9"
# Distance quantile levels to compute per source hex.
quantile_levels_csv = "0.1,0.5,0.9"
# Minimum trajectories per source hex to report a quantile (else dropped —
# avoids unstable quantiles from a handful of particles).
min_traj_per_hex = 30

# Baltic-wide map extent (degrees E / degrees N).
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66

# Per-panel height in inches (panel width is aspect-derived).
baltic_panel_height_in = 6

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
release_months = [int(x) for x in release_months_csv.split(",") if x]
quantile_levels = [float(x) for x in quantile_levels_csv.split(",")]

# %% [markdown]
# # Read key + pool distance histogram across years
#
# Layout: ``output_root/HexAggregates/HexAgg_key_r<radius>m.parquet`` and
# ``HexAgg_distance_r<radius>m_<regime>_<year>.parquet``. The release year
# is parsed from each filename so `release_doy → month` is leap-correct;
# the per-hex histograms are summed across the kept months/years.

# %%
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)

# Match exactly four year digits after the regime so `surface` does not also
# glob `surface_stokes_*` files (prefix collision).
_DIST_RE = re.compile(rf"HexAgg_distance_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$")
dist_files = sorted(
    store_root.glob(f"HexAgg_distance_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet")
)
if not dist_files:
    raise FileNotFoundError(
        f"no distance partitions for regime {regime!r} at {store_root} — run 024b."
    )

parts = []
for f in dist_files:
    year = int(_DIST_RE.search(f.name).group(1))
    df = pd.read_parquet(f).reset_index(drop=True)
    df["release_month"] = pd.to_datetime(
        (year * 1000 + df["release_doy"].astype("int32")).astype(str), format="%Y%j"
    ).dt.month
    parts.append(df)
dist = pd.concat(parts, ignore_index=True)
print(f"key: {len(key):,} hexes; pooled {len(dist_files)} year partition(s)")

if release_months:
    dist = dist[dist["release_month"].isin(release_months)]
    print(f"Keeping release months {release_months}: {len(dist):,} rows")
else:
    print(f"All release months: {len(dist):,} rows")

# Pooled per-hex histogram.
hist = dist.groupby(["release_hex", "distance_bin"])["n_traj"].sum().reset_index()

# %% [markdown]
# # Per-hex quantiles from the histogram
#
# For each source hex: cumulative count over `distance_bin`, then read off
# each quantile at the bin-left-edge (`distance_bin * distance_bin_km` km).
# Hexes with fewer than `min_traj_per_hex` trajectories are dropped. A few
# hundred hexes — plain pandas, no dask.

# %%
def hex_quantiles(hist, levels, bin_km, min_traj):
    rows = []
    for hex_id, g in hist.groupby("release_hex"):
        g = g.sort_values("distance_bin")
        cum = g["n_traj"].cumsum().to_numpy()
        total = int(cum[-1])
        if total < min_traj:
            continue
        edges = g["distance_bin"].to_numpy() * bin_km
        rec = {"release_hex": int(hex_id), "n_traj": total}
        for q in levels:
            # First bin whose cumulative count reaches q*total (side="left"),
            # reported at its left edge. q in [0,1] ⇒ idx < len, so the
            # min() clamp is only defensive.
            idx = int(np.searchsorted(cum, q * total))
            rec[q] = float(edges[min(idx, len(edges) - 1)])
        rows.append(rec)
    return pd.DataFrame(rows)


quant = hex_quantiles(hist, quantile_levels, distance_bin_km, min_traj_per_hex)
print(f"{len(quant):,} source hexes meet min_traj_per_hex={min_traj_per_hex}")

# %% [markdown]
# # Rendering helpers
#
# `to_hex_value_gdf` joins a per-hex value column to the key geometry
# (adapted from 025's `to_hex_gdf`). `hex_value_plot` adapts 025's
# `log_density_plot` — same plain EPSG:4326 axis and layout/registration
# overrides (aspect figsize, hex `edgecolor="face"/linewidth=0.4`, black
# coastline; see docs/visualisations.md) — but **without** `np.log10` /
# `cmap="viridis"`: distance is linear, so the default colormap with a
# matplotlib auto-ranged norm and a km colorbar is appropriate.

# %%
def to_hex_value_gdf(df, value_col, key_df):
    return (
        df[["release_hex", value_col]]
        .merge(key_df[["hex_id", "geometry"]], left_on="release_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def hex_value_plot(gdf, ax, extent, column, title=None, coast=None):
    # No missing_kwds (unlike 025/026): under-sampled hexes are dropped
    # upstream by hex_quantiles, so there are no NaN-valued rows to grey out.
    if not gdf.empty:
        gdf.plot(
            ax=ax, column=column, legend=True,
            legend_kwds={"label": "final displacement (km)"},
            edgecolor="face", linewidth=0.4, zorder=1,
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


# %%
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

# %% [markdown]
# # One map per quantile
#
# Each map colours source hexes by that quantile's final-displacement
# distance, auto-scaled independently (0 → that quantile's max).

# %%
for q in quantile_levels:
    gdf_q = to_hex_value_gdf(
        quant.rename(columns={q: "distance_km"}), "distance_km", key
    )
    fig, ax = plt.subplots(
        figsize=(baltic_panel_height_in * baltic_aspect, baltic_panel_height_in),
        layout="constrained",
    )
    hex_value_plot(
        gdf_q, ax, baltic_extent, "distance_km",
        title=f"quantile {q:g}", coast=coast_baltic,
    )
    plt.show()

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, distance_bin_km={distance_bin_km}")
print(f"release months: {release_months or 'all'}")
print(f"quantile levels: {quantile_levels}")
print(f"min trajectories per hex: {min_traj_per_hex}")
print(f"source hexes meeting the gate: {len(quant):,}")
for q in quantile_levels:
    print(
        f"  quantile {q:g}: "
        f"min={quant[q].min():.1f} "
        f"median={quant[q].median():.1f} "
        f"max={quant[q].max():.1f} km"
    )
