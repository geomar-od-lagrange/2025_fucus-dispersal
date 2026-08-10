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
# # Survival-weighted heatmaps: beaching folded into occupancy
#
# Parquet-only consumer of the survival-occupancy store built by
# `024e_BuildSurvivalOccupancy` (+ the `024a` key). For each elapsed-time
# horizon it draws three panels:
#
# 1. **Occupancy** — plain particle-residence density (`occ`), the baseline
#    (no beaching), on a log scale.
# 2. **Survival-weighted** — the same, with beaching progressively removed
#    (`surv = Σ exp(−A)`), on the *same* log scale, so the thinning is legible.
# 3. **Surviving fraction** — `surv / occ` per hex (linear 0–1), i.e. the mean
#    un-beached fraction of the particles *occupying* that hex at that age.
#    Read it as "how depleted is the mass here", NOT "how much beaching
#    happens here": `A` is a path integral accumulated before arrival, so a
#    hex reads low because the particles reaching it took wave-exposed
#    near-shore routes. The sink field — where stranding actually occurs — is
#    `029`'s where-stranded map, binned on `beach_hex`. Because each panel is
#    one age bin, every particle in it has the same elapsed time, so the
#    variation is route, not age.
#
# One regime per run; `release_month = 0` pools every monthly partition
# across years. No Dask, no zarrs — in the `025`/`026`/`029` lineage.

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
output_root = "../output"
regime = "surface_stokes"
hex_radius = 6000
age_bin_days = 10
# 0 = pool all monthly partitions across years; 1..12 = that month only.
release_month = 0

# Reference onshore-Stokes forcing (m/s) of the partitions to read — part of
# the store filename. tau(w_tau) == tau0 by construction.
w_tau = 0.05
# Base timescale (h) of the partitions to read — also part of the filename.
tau0_hours = 480.0

# Elapsed-time horizons to map (days); each a multiple of age_bin_days and
# within the store's occupancy_max_days.
time_horizons_days_csv = "20,50,100"

# Map extent (degrees). The 024a key tiles the whole BSH domain including the
# North Sea, which is empty for Fucus and costs ~35% of panel height; cropping
# to the Baltic proper buys ~1.35x px/km for free. Set all four to 0 to fall
# back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0

cmap = "viridis"
panel_height_in = 6
fig_dpi_scale = 3

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
# Hex seam stroke. edgecolor="face" means this is not a visible outline -- it
# closes the ~1 px anti-aliasing seam between adjacent polygons so the grid
# reads as a continuous field. The seam is a fixed PIXEL artifact, so a fixed
# point width makes the resulting hex dilation DPI-invariant (17.5% of hex
# width at every dpi). Pinning it to ~1 px instead lets dilation fall as
# resolution rises: ~8.6% at fig_dpi_scale=3 on the Baltic crop.
hex_seam_lw = 1.1 * 72 / (100 * fig_dpi_scale)

figure_dir = output_root / "Figures" / "030"
figure_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Read key + pool survival-occupancy partitions

# %%
store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

month_suffix = f"_m{release_month:02d}" if release_month else ""
month_re = rf"_m{release_month:02d}" if release_month else r"_m\d{2}"
wh_suffix = f"_t{tau0_hours:g}_wt{w_tau:g}".replace(".", "p")
_PART_RE = re.compile(
    rf"HexAgg_survocc_r{hex_radius}m_{regime}_(\d{{4}}){month_re}{re.escape(wh_suffix)}\.parquet$"
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

# %% [markdown]
# # Rendering helpers

# %%
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
            legend=True, edgecolor="face", linewidth=hex_seam_lw, zorder=1,
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
if any((extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max)):
    lon_min, lat_min, lon_max, lat_max = (
        extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max
    )
else:
    lon_min, lat_min, lon_max, lat_max = key.total_bounds
extent = [lon_min, lon_max, lat_min, lat_max]
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)

# %% [markdown]
# # Comparison grid: occupancy vs survival-weighted
#
# Rows = horizons; columns = occupancy, survival-weighted. Both share one
# `LogNorm` (across all horizons) so the age-thinning and the beaching-removal
# are directly readable against each other.
#
# Surviving fraction is drawn separately below: it is a ratio in [0, 1] on a
# linear scale, so putting it in this grid forces an unrelated quantity to
# share a row with two log-density panels and invites reading it as a third
# density map.

# %%
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
    nrows=nrows, ncols=2,
    figsize=(panel_height_in * domain_aspect * 2, panel_height_in * nrows),
    layout="constrained", squeeze=False,
)
for i, h in enumerate(time_horizons_days):
    hex_map(occ_gdfs[h], axes[i, 0], norm=dens_norm, title=f"occupancy — {h} d")
    hex_map(surv_gdfs[h], axes[i, 1], norm=dens_norm, title=f"survival-weighted — {h} d")
fig_path = figure_dir / f"SurvivalHeatmaps_{regime}_r{hex_radius}m{month_suffix}{wh_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Surviving fraction, one panel per horizon
#
# `surv / occ` per hex: the mean un-beached fraction of the particles
# *occupying* that hex at that age — a history integral, not a local beaching
# rate (see the header). Fixed linear 0–1 across panels so the age progression
# is readable; the summary below prints the realised range, which at a
# rare-event `w_tau` sits well above 0 and makes the fixed scale look flat.

# %%
ncols = len(time_horizons_days)
fig, axes = plt.subplots(
    nrows=1, ncols=ncols,
    figsize=(panel_height_in * domain_aspect * ncols, panel_height_in),
    layout="constrained", squeeze=False,
)
for j, h in enumerate(time_horizons_days):
    hex_map(frac_gdfs[h], axes[0, j], vmin=0.0, vmax=1.0,
            title=f"surviving fraction — {h} d")
fig_path = figure_dir / f"SurvivalFraction_{regime}_r{hex_radius}m{month_suffix}{wh_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
for h in time_horizons_days:
    g = frac_gdfs[h]
    if not g.empty:
        print(f"  {h:>4} d: surviving fraction per hex "
              f"min {g['value'].min():.3f}, median {g['value'].median():.3f}, "
              f"max {g['value'].max():.3f}")
plt.show()

# %% [markdown]
# # Validation / summary

# %%
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
