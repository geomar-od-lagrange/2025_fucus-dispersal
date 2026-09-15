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
# 024a hex grid as one panel per quantile level. Reads the distance
# histogram store built by 024b (+ the 024a key for geometry) — no
# trajectory zarrs, no Dask cluster. One regime and one release `season`
# per run, pooled across every available release year.
#
# Quantiles are derived from the pooled per-hex histogram (cumulative count
# over `distance_bin`). Histograms are additive, so pooling across years is
# summing partitions — unlike pre-computed quantiles, which cannot be
# averaged across years.
#
# The figure is written print-ready (full page width, 300 dpi) and the same
# per-hex quantiles are exported as GeoJSON for collaborators.

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

# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year.
season = "ALL"
# Distance quantile levels to compute per source hex.
quantile_levels_csv = "0.1,0.5,0.9"
# Minimum trajectories per source hex to report a quantile (else dropped —
# avoids unstable quantiles from a handful of particles).
min_traj_per_hex = 30

# Map extent (degrees E / degrees N). Cropped to the Baltic proper: the 024a
# key tiles the whole BSH domain including the North Sea, which carries no
# Fucus source hexes and only costs panel height (as 029/030).
baltic_lon_min = 9.0
baltic_lon_max = 30.7
baltic_lat_min = 53.0
baltic_lat_max = 66.0

# %% [markdown]
# # Parse parameters
#
# The season → month mapping lives here rather than in the parameters cell:
# papermill only injects primitives, and the mapping is a fixed convention,
# not a knob.

# %%
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}

output_root = Path(output_root)
if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
quantile_levels = [float(x) for x in quantile_levels_csv.split(",")]

# Print-ready output: every saved figure is one full text-width (180 mm)
# figure at 300 dpi, so panels land in the manuscript at their final size
# (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300

figure_dir = output_root / "Figures" / "027"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "027"
export_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Read key + pool distance histogram across years
#
# Layout: ``output_root/HexAggregates/HexAgg_key_r<radius>m.parquet`` and
# ``HexAgg_distance_r<radius>m_<regime>_<year>.parquet``. The release year
# is parsed from each filename so `release_doy → month` is leap-correct;
# the per-hex histograms are summed across the season's months and all years.

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
    df["release_year"] = year
    df["release_month"] = pd.to_datetime(
        (year * 1000 + df["release_doy"].astype("int32")).astype(str), format="%Y%j"
    ).dt.month
    parts.append(df)
dist = pd.concat(parts, ignore_index=True)
print(f"key: {len(key):,} hexes; pooled {len(dist_files)} year partition(s)")

dist = dist[dist["release_month"].isin(season_months)]
print(f"season {season} (months {season_months}): {len(dist):,} rows")
if dist.empty:
    raise ValueError(f"no releases in season {season} for regime {regime!r}")

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
            rec[f"q{q:g}_km"] = float(edges[min(idx, len(edges) - 1)])
        rows.append(rec)
    return pd.DataFrame(rows)


quant = hex_quantiles(hist, quantile_levels, distance_bin_km, min_traj_per_hex)
print(f"{len(quant):,} source hexes meet min_traj_per_hex={min_traj_per_hex}")
if quant.empty:
    raise ValueError(
        f"no source hex reaches min_traj_per_hex={min_traj_per_hex} in season {season}"
    )

# %% [markdown]
# # Rendering helpers
#
# `hex_value_plot` adapts 025's `log_density_plot` — same plain EPSG:4326
# axis and hex-registration overrides (`edgecolor="face"`, black coastline;
# see docs/visualisations.md) — but **without** the `np.log10`: distance is a
# linear physical quantity, so the default colormap with a matplotlib
# auto-ranged norm and a km colorbar is appropriate.
#
# The colorbar is an inset axes at axes-fraction height 1.0, so it is exactly
# the map height at any figure size — a fixed-aspect map never fills its
# gridspec cell, and a cell-sized colorbar overshoots it.

# %%
def hex_value_plot(gdf, ax, extent, column, title=None, coast=None):
    # No missing_kwds (unlike 025/026): under-sampled hexes are dropped
    # upstream by hex_quantiles, so there are no NaN-valued rows to grey out.
    if not gdf.empty:
        cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
        gdf.plot(
            ax=ax, column=column, legend=True, cax=cax,
            legend_kwds={"label": "final displacement (km)"},
            edgecolor="face", linewidth=hex_seam_lw, zorder=1,
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


def grid_height_in(nrows, ncols, aspect, width_in):
    """Height seed for a grid of fixed-aspect map panels at a fixed figure
    width. The figure width is set by the page, so only the height is free:
    each column gives its map ~72 % of its width (the rest is the inset
    colorbar and its tick labels), and ~10 % is added for titles and padding.
    constrained_layout does the packing; there is no measure-rescale loop."""
    return 1.10 * nrows * (0.72 * width_in / ncols) / aspect


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
# Hex seam stroke, as 029/031: edgecolor="face" means this is not a visible
# outline, it closes the ~1 px anti-aliasing seam between adjacent polygons.
# The seam is a fixed PIXEL artefact, so pin the stroke to ~1 px at the
# savefig dpi rather than to a fixed point width.
hex_seam_lw = 1.1 * 72 / FIGURE_DPI

quant_gdf = quant.merge(
    key[["hex_id", "geometry"]], left_on="release_hex", right_on="hex_id"
).drop(columns="hex_id").pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")

# %% [markdown]
# # One panel per quantile
#
# Each panel colours source hexes by that quantile's final-displacement
# distance, auto-scaled independently (0 → that quantile's max), so the
# quantiles are read as three separate fields rather than one thinning one.

# %%
ncols = len(quantile_levels)
fig, axes = plt.subplots(1, ncols, layout="constrained", squeeze=False)
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(1, ncols, baltic_aspect, FIGURE_WIDTH_IN)
)
for ax, q in zip(axes.flat, quantile_levels):
    hex_value_plot(
        quant_gdf, ax, baltic_extent, f"q{q:g}_km",
        title=f"quantile {q:g}", coast=coast_baltic,
    )
fig_path = figure_dir / f"DistanceQuantiles_{regime}_r{hex_radius}m_{season}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # GeoJSON export
#
# The same per-hex quantiles with the 024a hex geometry, EPSG:4326, one
# feature per source hex and one column per quantile level — the figure's
# values in a form collaborators can open in R or QGIS.

# %%
export_path = export_dir / f"027_{regime}_{season}.geojson"
quant_gdf.to_file(export_path, driver="GeoJSON")
print(f"wrote {export_path} ({len(quant_gdf):,} hexes)")

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, distance_bin_km={distance_bin_km}")
print(f"season: {season} (months {season_months})")
print(f"quantile levels: {quantile_levels}")
print(f"min trajectories per hex: {min_traj_per_hex}")
print(f"source hexes meeting the gate: {len(quant):,}")
for q in quantile_levels:
    col = f"q{q:g}_km"
    print(
        f"  quantile {q:g}: "
        f"min={quant[col].min():.1f} "
        f"median={quant[col].median():.1f} "
        f"max={quant[col].max():.1f} km"
    )
