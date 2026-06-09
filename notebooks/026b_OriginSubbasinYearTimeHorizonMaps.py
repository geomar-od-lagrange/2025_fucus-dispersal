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
# # Time-horizon dispersal maps per origin subbasin and release year
#
# 026a split further by **release year**: for each (origin subbasin, year)
# pair the same four-panel elapsed-time-horizon density map — target-hex
# occupancy at 10/20/50/100 d — restricted to trajectories released from
# that subbasin in that year. This exposes interannual variability of a
# subbasin's dispersal that 026a hides by pooling years.
#
# Within a subbasin the colour scale is **shared across its years** (and
# horizons), so a faint year reads as genuinely lower occupancy rather than
# a rescaling artefact — magnitudes are comparable year to year. The scale
# is independent between subbasins (cross-subbasin totals span orders of
# magnitude). Origin subbasin is the `helcom_subbasin` id of each
# `release_hex` (id→name map from the key's JSON sidecar); `_outside` and
# land-seeded releases carry no origin and are reported, not mapped. Reads
# the 024a key + 024 counts — no trajectory zarrs, no Dask cluster. One
# regime per run; August/September releases kept, split by year.

# %%
import json
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
# Which regime's counts partitions to read; one regime per run.
regime = "surface"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the counts store (must match 024). Horizons must
# be whole multiples of this.
age_bin_days = 10

# Release months to keep (Aug, Sep). Empty ⇒ all releases.
release_months_csv = "8,9"
# Elapsed-time horizons to map, in days. Each must be a multiple of
# age_bin_days.
time_horizons_days_csv = "10,20,50,100"
# Origin subbasins to map, by HELCOM name (comma-separated). Empty ⇒ every
# named subbasin that seeds releases.
origin_subbasins_csv = ""
# Release years to map (comma-separated). Empty ⇒ every year present.
release_years_csv = ""

# Colormap: log-density spans several decades, so a perceptually uniform
# map is load-bearing (justified in docs/visualisations.md, as for 025).
cmap = "viridis"

# Per-panel height in inches (panel widths are aspect-derived).
panel_height_in = 6
# Figure DPI as a multiple of the matplotlib default. >1 sharpens the
# inline/embedded raster panels — the one plotting default we override here
# (rationale in docs/visualisations.md).
fig_dpi_scale = 2

# %% [markdown]
# # Parse parameters

# %%
output_root = Path(output_root)
release_months = [int(x) for x in release_months_csv.split(",") if x]
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
origin_subbasins_filter = [s.strip() for s in origin_subbasins_csv.split(",") if s.strip()]
release_years_filter = [int(y) for y in release_years_csv.split(",") if y]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would fall in a neighbouring age bin"
    )
# Set from the stock default so re-running this cell is idempotent.
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]

# %% [markdown]
# # Read key + subbasin names + pool counts across years
#
# Layout: flat files under ``output_root/HexAggregates/`` —
# ``HexAgg_key_r<radius>m.parquet`` (+ a ``.json`` sidecar carrying the
# subbasin id→name map) and
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``. The release year is
# parsed from each filename and kept as a column so the maps can split on it
# (and the `release_doy → month` conversion is leap-correct).

# %%
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(key_path.with_suffix(".json").read_text())["subbasin_id_to_name"].items()
}

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
    df["release_year"] = year
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

# %% [markdown]
# # Attach origin subbasin; resolve subbasins and years to map
#
# Each trajectory's origin subbasin is the `helcom_subbasin` of its
# `release_hex`. Land-seeded releases (`release_hex` -1) map to NaN; releases
# outside any named polygon carry id -1 (`_outside`). Both lack a named
# origin and are excluded from the maps (counted below).

# %%
counts["origin_subbasin"] = counts["release_hex"].map(key.set_index("hex_id")["helcom_subbasin"])

origin_codes_present = sorted(
    int(c) for c in counts["origin_subbasin"].dropna().unique() if c >= 0
)
name_to_id = {v: k for k, v in subbasin_id_to_name.items()}
unknown = [n for n in origin_subbasins_filter if n not in name_to_id]
if unknown:
    raise KeyError(f"unknown subbasin name(s) {unknown}; known: {sorted(name_to_id)}")
filter_ids = {name_to_id[n] for n in origin_subbasins_filter}
origin_codes = [c for c in origin_codes_present if not filter_ids or c in filter_ids]

years_present = sorted(int(y) for y in counts["release_year"].unique())
missing_years = [y for y in release_years_filter if y not in years_present]
if missing_years:
    print(f"requested years with no partitions (skipped): {missing_years}")
years = [y for y in years_present if not release_years_filter or y in release_years_filter]

n_landseed = int(counts["origin_subbasin"].isna().sum())
n_outside_hex = int(counts.loc[counts["origin_subbasin"] == -1, "release_hex"].nunique())
print(
    f"origin subbasins to map ({len(origin_codes)}): "
    f"{', '.join(subbasin_id_to_name[c] for c in origin_codes)}"
)
print(f"release years to map ({len(years)}): {years}")
print(
    f"excluded origins: {n_outside_hex:,} release hex(es) outside any named subbasin "
    f"(_outside), {n_landseed:,} land-seeded count rows (release_hex -1)"
)

# %% [markdown]
# # Rendering helpers
#
# Identical to 026/026a (`to_hex_gdf` joins a per-`target_hex` count table to
# the key geometry; `log_density_plot` draws hex polygons + coastline on a
# plain EPSG:4326 axis, colouring through a shared `LogNorm`).

# %%
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


def log_density_plot(gdf, ax, extent, norm, title=None, coast=None):
    """Plot hex density + coastline on a plain (non-cartopy) lon/lat axis.
    Everything stays in EPSG:4326 so hexes and coastline share one
    coordinate system. Colour scales the raw ``count`` through a shared
    ``LogNorm``, so each panel's colorbar reads in particle counts (not
    log10) while spanning the several decades of density. Empty ``gdf``
    skips the hex layer so empty panels render with stable layout."""
    if not gdf.empty:
        # cmap + log norm justified in docs/visualisations.md (as 025);
        # legend=True draws a per-panel colorbar keyed to the shared scale.
        gdf.plot(
            ax=ax, column="count", cmap=cmap, norm=norm, legend=True,
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
# The hex key decides the extent — no manual crop. The key tiles the full
# BSH model domain (North Sea included), so this shows the whole covered
# area; the dispersal cloud occupies wherever it actually reaches.
domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = key.total_bounds
domain_extent = [domain_lon_min, domain_lon_max, domain_lat_min, domain_lat_max]

# Natural Earth 10m coastline via cartopy's shapereader (cached locally),
# clipped to the domain extent once.
_coast_gdf = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
)
coast = _coast_gdf.clip(box(domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max))

# Aspect ratio that keeps 1° lon at lat_mean visually equal to 1° lat.
domain_aspect = (
    (domain_lon_max - domain_lon_min) * np.cos(np.radians(0.5 * (domain_lat_min + domain_lat_max)))
) / (domain_lat_max - domain_lat_min)

# %% [markdown]
# # Map per origin subbasin and year
#
# Outer loop over origin subbasins; inner loop over years. Per subbasin the
# four-panel hex GeoDataFrames are built for every (year, horizon) first so
# one `LogNorm` can span the subbasin's whole year set — making years
# magnitude-comparable — then one four-panel figure is drawn per year that
# has data. Target hexes span the whole domain.

# %%
horizon_age_bins = {h: h // age_bin_days for h in time_horizons_days}
ncols = 2
nrows = int(np.ceil(len(time_horizons_days) / ncols))

for code in origin_codes:
    counts_sb = counts[counts["origin_subbasin"] == code]
    # All (year, horizon) panels for this subbasin, so the colour scale can
    # be shared across years (built once, reused for plotting below).
    gdfs = {
        (y, h): to_hex_gdf(
            counts_sb[
                (counts_sb["release_year"] == y)
                & (counts_sb["age_bin"] == horizon_age_bins[h])
            ],
            key,
        )
        for y in years
        for h in time_horizons_days
    }
    nonempty = [g["count"] for g in gdfs.values() if not g.empty]
    if not nonempty:
        print(f"{subbasin_id_to_name[code]}: no occupied target hexes — skipping")
        continue
    all_counts = pd.concat(nonempty)
    norm = LogNorm(vmin=all_counts.min(), vmax=all_counts.max())

    for y in years:
        if all(gdfs[(y, h)].empty for h in time_horizons_days):
            continue
        fig, axes = plt.subplots(
            nrows=nrows, ncols=ncols,
            figsize=(
                panel_height_in * domain_aspect * ncols,
                panel_height_in * nrows,
            ),
            layout="constrained",
        )
        fig.suptitle(f"origin: {subbasin_id_to_name[code]} — {y}")
        for ax, h in zip(axes.flat, time_horizons_days):
            log_density_plot(gdfs[(y, h)], ax, domain_extent, norm, title=f"{h} d", coast=coast)
        for ax in axes.flat[len(time_horizons_days):]:
            ax.set_visible(False)
        plt.show()

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, age_bin_days={age_bin_days}")
print(f"release months: {release_months or 'all'}")
print(f"origin subbasins mapped: {len(origin_codes)}; years: {years}")
for code in origin_codes:
    counts_sb = counts[(counts["origin_subbasin"] == code) & (counts["target_hex"] != -1)]
    in_horizons = counts_sb[counts_sb["age_bin"].isin(horizon_age_bins.values())]
    by_year = in_horizons.groupby("release_year")["n_obs"].sum()
    per_year = ", ".join(f"{y}={int(by_year.get(y, 0)):,}" for y in years)
    print(f"  {subbasin_id_to_name[code]} (sum n_obs over mapped horizons): {per_year}")
