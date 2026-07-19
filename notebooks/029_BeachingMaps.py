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
# # Beaching maps
#
# Lightweight parquet-only consumer of the beaching store built by
# `024d_BuildBeaching` (+ the `024a` key) — no trajectory zarrs, no Dask.
# One regime per run; release years pooled by globbing the store partitions,
# the release year parsed per filename so `release_doy → month` is
# leap-correct (as `026`). Draws:
#
# 1. **Where-stranded density** — beached weight (expected particles) per
#    `beach_hex` (log scale), pooled over shore types (`trap` is degenerate
#    in 024d, so the wall/flat label is not a model result).
# 2. **Beached fraction per source hex** — of the drifters released in each
#    hex, what fraction strands within the viability window (linear 0–1).
# 3. **Beaching age horizons** — cumulative where-stranded density for
#    strandings at age ≤ T (the beaching analogue of `026`'s horizons).
# 4. **Cumulative beached fraction vs. age** — the stranding time course.

# %%
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

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Which regime's beaching partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024d). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the beaching store (must match 024d). Horizons must
# be whole multiples of this.
age_bin_days = 10
# Release month to analyse: 0 = pool all months (every `_mMM` partition,
# across years); 1..12 = keep just that month. Selects which monthly
# partitions 024d wrote are read.
release_month = 0

# Onshore-Stokes half-saturation (m/s) of the partitions to read — part of
# the store filename, so this selects one member of the w_half sweep.
w_half = 0.05

# Beaching age horizons to map, in days. Each must be a multiple of
# age_bin_days.
time_horizons_days_csv = "10,20,50"

# Colormap: log where-stranded density spans several decades, so a
# perceptually uniform map is load-bearing (as 025/026; docs/visualisations.md).
cmap = "viridis"
# Per-panel height in inches (panel widths are aspect-derived).
panel_height_in = 6
# Figure DPI as a multiple of the matplotlib default (sharpens raster panels;
# the one plotting default overridden here, as 026 — docs/visualisations.md).
fig_dpi_scale = 2

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d"
    )
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]

figure_dir = output_root / "Figures" / "029"
figure_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Read key + pool beaching partitions across years
#
# Layout: flat files under ``output_root/HexAggregates/`` —
# ``HexAgg_key_r<radius>m.parquet`` and the per-(year, month) partitions
# ``HexAgg_beaching_r<radius>m_<regime>_<year>_mMM.parquet`` written by 024d.
# `release_month = 0` pools every month across every year; a nonzero month
# keeps just that month (across years). Matching only `_mMM` files means any
# ad-hoc whole-year build (no suffix) is ignored, so months never double-count.

# %%
store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

# Figure-filename tag: a specific month, or "" when pooling all months.
month_suffix = f"_m{release_month:02d}" if release_month else ""
month_re = rf"_m{release_month:02d}" if release_month else r"_m\d{2}"
wh_suffix = f"_wh{w_half:g}".replace(".", "p")
_PART_RE = re.compile(
    rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}}){month_re}{re.escape(wh_suffix)}\.parquet$"
)
beaching_files = [
    f for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet"))
    if _PART_RE.search(f.name)
]
if not beaching_files:
    raise FileNotFoundError(
        f"no beaching partitions for regime {regime!r}"
        + (f", month {release_month}" if release_month else "")
        + f" at {store_root} — run 024d."
    )

beaching = pd.concat(
    [pd.read_parquet(f).reset_index(drop=True) for f in beaching_files],
    ignore_index=True,
)
print(f"key: {len(key):,} hexes; pooled {len(beaching_files)} (year, month) partition(s)")
print(f"beaching rows: {len(beaching):,}")

# %% [markdown]
# # Rendering helpers
#
# `hex_gdf` sums a value column over a hex column and joins to the key
# geometry (dropping the `-1` sentinel); `hex_map` draws hex polygons +
# coastline on a plain EPSG:4326 axis (as 025/026).

# %%
def hex_gdf(df, hex_col, value_col="weight"):
    grp = df[df[hex_col] >= 0].groupby(hex_col)[value_col].sum()
    return (
        grp.rename("value").reset_index()
        .merge(key[["hex_id", "geometry"]], left_on=hex_col, right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def log_norm(values):
    """Shared LogNorm floored to four decades below the peak. Weighted
    deposition spans many decades — the sparse tail carries fractional
    (< 1) weight — so a floor keeps the colour range readable."""
    vmax = float(values.max())
    vmin = max(float(values[values > 0].min()), vmax / 1e4)
    return LogNorm(vmin=vmin, vmax=vmax)


def hex_map(gdf, ax, norm=None, title=None):
    """Hex choropleth + coastline on a lon/lat axis. `norm=None` → linear
    default scale; a `LogNorm` → shared log scale across panels."""
    if not gdf.empty:
        gdf.plot(
            ax=ax, column="value", cmap=cmap, norm=norm, legend=True,
            edgecolor="face", linewidth=0.4, zorder=1,
        )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


# %%
# The hex key tiles the full BSH domain (North Sea included); the extent
# follows it, and beaching occupies wherever it actually reaches.
lon_min, lat_min, lon_max, lat_max = key.total_bounds
extent = [lon_min, lon_max, lat_min, lat_max]
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)

beached = beaching[beaching["beach_hex"] >= 0]

# %% [markdown]
# # Where-stranded density
#
# Log-scale beached weight (expected particles) per stranding hex, pooled over
# shore types. The store's `shore_type` is *not* split out here: `024d` runs
# with a degenerate `trap` (`trap_flat == trap_wall`), so the `wall`/`flat`
# label does not affect the rate and a split map would invite reading resolved
# coastal morphology into what is only the BSH tidal-flat flag. Split it again
# when a real substrate classification drives `trap`.

# %%
gdf_stranded = hex_gdf(beached, "beach_hex")
strand_norm = log_norm(gdf_stranded["value"])

fig, ax = plt.subplots(
    figsize=(panel_height_in * domain_aspect, panel_height_in),
    layout="constrained",
)
hex_map(gdf_stranded, ax, norm=strand_norm, title="stranded weight")
fig_path = figure_dir / f"WhereStranded_{regime}_r{hex_radius}m{month_suffix}{wh_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Beached fraction per source hex
#
# Of the drifters released in each source hex, the fraction that strands
# within the viability window (linear 0–1). Highlights which source regions
# lose most of their propagules to the coast.

# %%
released = beaching.groupby("release_hex")["weight"].sum().rename("released")
stranded = beached.groupby("release_hex")["weight"].sum().rename("stranded")
frac = (
    pd.concat([released, stranded], axis=1).fillna({"stranded": 0})
    .assign(value=lambda d: d["stranded"] / d["released"])
    .reset_index()
)
frac = frac[frac["release_hex"] >= 0].merge(
    key[["hex_id", "geometry"]], left_on="release_hex", right_on="hex_id"
).pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")

fig, ax = plt.subplots(
    figsize=(panel_height_in * domain_aspect, panel_height_in), layout="constrained"
)
# Linear default scale (fraction in [0, 1]); no norm override needed.
hex_map(frac, ax, norm=None, title="beached fraction per source hex")
fig_path = figure_dir / f"BeachedFraction_{regime}_r{hex_radius}m{month_suffix}{wh_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Beaching age horizons
#
# Cumulative where-stranded density for strandings occurring at age ≤ T.
# `beach_age_bin` bins age at `age_bin_days`, so horizon T keeps
# `beach_age_bin < T // age_bin_days` (strandings in bins fully within T).
# A shared `LogNorm` across horizons makes the fill-in over time legible.

# %%
horizon_bins = {h: h // age_bin_days for h in time_horizons_days}
gdfs = {
    h: hex_gdf(beached[beached["beach_age_bin"] < horizon_bins[h]], "beach_hex")
    for h in time_horizons_days
}
all_h = pd.concat([g["value"] for g in gdfs.values() if not g.empty])
horizon_norm = log_norm(all_h)

ncols = len(time_horizons_days)
fig, axes = plt.subplots(
    1, ncols, figsize=(panel_height_in * domain_aspect * ncols, panel_height_in),
    layout="constrained", squeeze=False,
)
for ax, h in zip(axes.flat, time_horizons_days):
    hex_map(gdfs[h], ax, norm=horizon_norm, title=f"stranded by {h} d")
fig_path = figure_dir / f"BeachingHorizons_{regime}_r{hex_radius}m{month_suffix}{wh_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Cumulative beached fraction vs. age
#
# Share of all released drifters stranded by each elapsed-time horizon — the
# stranding time course. Pooled over shore types, for the reason given at the
# where-stranded map.

# %%
total_released = float(beaching["weight"].sum())
age_bins_present = sorted(beached["beach_age_bin"].unique())
curve = pd.Series(
    [
        beached[beached["beach_age_bin"] <= b]["weight"].sum() / total_released
        for b in age_bins_present
    ],
    index=pd.Index([(b + 1) * age_bin_days for b in age_bins_present], name="age_days"),
)
ax = curve.plot()
ax.set_ylabel("cumulative beached fraction")
plt.show()

# %% [markdown]
# # Validation / summary

# %%
by_type = beaching.groupby("shore_type")["weight"].sum()
n_beached = float(beached["weight"].sum())
print(f"regime={regime}, hex_radius={hex_radius} m, "
      + (f"month={release_month}, " if release_month else "")
      + f"age_bin_days={age_bin_days}")
print(f"  drifters (Σweight): {total_released:,.0f}")
print(f"  beached:           {n_beached:,.0f} "
      f"({100 * n_beached / max(total_released, 1):.1f}%)")
print(f"  shore_type split:  "
      + ", ".join(f"{k}={float(v):,.0f}" for k, v in by_type.items())
      + "  [classification label only; trap is degenerate in 024d]")
print(f"  stranding hexes:   {beached['beach_hex'].nunique():,}")
print(f"  source hexes:      {frac['release_hex'].nunique():,}")
