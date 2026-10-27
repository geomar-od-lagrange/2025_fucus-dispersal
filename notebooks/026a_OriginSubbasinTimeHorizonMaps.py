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
# # Time-horizon relative-density maps per origin subbasin
#
# 026 split by **origin** (HELCOM) subbasin: for each subbasin that seeds
# releases, the same elapsed-time-horizon relative-density maps 026
# produces, but restricted to trajectories released from that subbasin.
# Target hexes still range over the whole domain (a subbasin's particles
# disperse anywhere); only the release side is filtered.
#
# Two stores feed this, selected by the `member` parameter, and the value
# read from either is called `weight` from then on:
#
# | `member` | store | `weight` | origin subbasin |
# |---|---|---|---|
# | `""` | 024 counts (`HexAgg_counts_*`) | `n_obs`, every sample counts 1 | key `helcom_subbasin` of `release_hex` |
# | rate-model tag | 024f survocc (`HexAgg_survocc_*_<member>`) | `surv`, each sample weighted by its un-beached fraction `exp(−A)` | the store's `origin_subbasin` |
#
# Both are read over the same age window `[0, max_age_days)`. The
# subbasin id→name map lives in the key's JSON sidecar. Releases outside
# any named subbasin (origin -1) carry no origin and are reported, not
# mapped.
#
# Each origin is normalised over **its own** releases, so a panel answers
# "where does *this* source's particle-time go?", not "how big is this
# source". In the per-origin figures each panel carries its own colorbar, so
# colour is not comparable between them. The overview figures (one per
# `overview_horizons_days_csv` horizon, every origin's dilution map on one
# page) share one colour scale: there a colour is the same share of *that
# origin's* particle-time per km² in every panel.
#
# A horizon `T` is **cumulative**: every `age_bin` whose window starts
# before `T`, i.e. elapsed time in `[0, T)`. The map is therefore the
# particle-time the cloud has accumulated by `T`, not a 10-day snapshot at
# `T`.
#
# Density is relative, as in 025:
#
# - `percentage` — 100 · weight(hex) / Σ weight, normalised **per horizon**,
#   so each panel shows that horizon's own share of particle-time. Linear
#   colour scale.
# - `dilution` — that share per km² of hex water area (`water_area_m2`
#   from the 024a key). Log colour scale.

# %%
import json
import re
import unicodedata
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the data twin (BSH static H0).
data_root = "../data"
# Read root of the hex-aggregate store; write root for figures and exports.
output_root = "../output"
# Which regime's partitions to read; one regime per run.
regime = "surface_stokes"
# Rate-model member tag. "" reads the unweighted 024 counts store (n_obs);
# any other string reads the 024f survocc store for that member and maps its
# survival-weighted `surv`. Opaque — never parsed.
member = ""
# Release season: DJF, MAM, JJA, SON or ALL (every month).
season = "ALL"
# Release years pooled: "first-last" (inclusive) or one year ("2019").
# Partitions of other years stay on disk unopened (2025 is partial).
release_years = "2016-2024"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the store (must match 024 / 024f). Horizons must
# be whole multiples of this.
age_bin_days = 10
# Reduction window (days): only ages in [0, max_age_days) are read. Cuts the
# uncapped counts store to the window 024f was built with, so weighted and
# unweighted maps cover the same ages. A multiple of age_bin_days.
max_age_days = 180
# Elapsed-time horizons to map, in days. Each must be a multiple of
# age_bin_days and ≤ max_age_days.
time_horizons_days_csv = "30,60,90,180"
# Horizons of the overview figures (one figure per horizon, one dilution
# panel per origin on a shared colour scale), in days: same rules as the map
# horizons.
overview_horizons_days_csv = "30,180"
# Map extent of the overview panels, "lon_min,lon_max,lat_min,lat_max" in
# degrees: the Baltic plus the Kattegat/Skagerrak entrance. The run prints
# how much above-floor particle-time falls outside it.
overview_extent_csv = "8,30.7,53,66"

# Origin subbasins to map, by HELCOM name (comma-separated). Empty ⇒ every
# named subbasin that seeds releases.
origin_subbasins_csv = ""

# Isobath to overlay on every panel, in metres of BSH H0 (floor position,
# see docs/h0_semantics.md). Fucus vesiculosus no longer grows below ~3 m.
isobath_m = 3.0

# %% [markdown]
# # Parse parameters
#
# `SEASON_MONTHS` lives here rather than in the parameters cell because
# papermill only injects primitives; `season` selects one entry.

# %%
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}
if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
# `release_years` is "first-last" (inclusive) or a single year.
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"
print(f"release years in the analysis: {analysis_years[0]}–{analysis_years[-1]}")

output_root = Path(output_root)
data_root = Path(data_root)
store_root = output_root / "HexAggregates"
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
overview_horizons_days = [int(x) for x in overview_horizons_days_csv.split(",") if x]
overview_extent = [float(x) for x in overview_extent_csv.split(",")]
assert len(overview_extent) == 4, f"overview_extent_csv {overview_extent_csv!r} needs 4 numbers"
assert max_age_days % age_bin_days == 0, (
    f"max_age_days {max_age_days} d is not a multiple of age_bin_days {age_bin_days} d"
)
for h in [*time_horizons_days, *overview_horizons_days]:
    assert h > 0, f"horizon {h} d holds no age bin"
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would cut a bin in half"
    )
    assert h <= max_age_days, f"horizon {h} d is past max_age_days {max_age_days} d"
n_age_bins = max_age_days // age_bin_days

member_tag = member or "unweighted"
print(f"regime={regime}, season={season}, member={member_tag}, window [0, {max_age_days}) d")

# Print-ready output: every saved figure is one full text-width (180 mm)
# page-width figure at 300 dpi, so panels land in the manuscript at their
# final size (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Vertical inches a panel needs on top of its map for the horizontal
# colorbar, the colorbar label and the panel title.
PANEL_DECORATION_IN = 0.62

origin_subbasins_filter = [s.strip() for s in origin_subbasins_csv.split(",") if s.strip()]

figure_dir = output_root / "Figures" / "026a"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "026a"
export_dir.mkdir(parents=True, exist_ok=True)


def _slug(name):
    """Filename-safe slug of a HELCOM subbasin name. NFKD + ASCII fold
    transliterates accents (Å → A) before collapsing every non-alphanumeric
    run to a single underscore."""
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^0-9A-Za-z]+", "_", ascii_name).strip("_")

# %% [markdown]
# # Read key + pool the store across years
#
# Layout: flat files under ``output_root/HexAggregates/`` — the key
# ``HexAgg_key_r<radius>m.parquet``, the counts partitions
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet`` (one per release
# year) and the survocc partitions
# ``HexAgg_survocc_r<radius>m_<regime>_<year>[_mMM]_<member>.parquet`` (one
# per release month, or one per year when 024f ran with
# `release_month = 0`). Every partition is read and reduced immediately to
# per-`(origin_subbasin, age_bin, target_hex)` sums — a counts partition is O(10^8) rows, the reduction
# O(10^6). The season's release days and the age window are pushed into the
# parquet read as filters; the doy set is built per partition year so the
# month↔doy mapping is leap-correct, and a monthly survocc partition outside
# the season is not opened at all. The release year comes from the filename.
#
# Sentinels: `target_hex = -1` (no in-key position — the drifter left the
# BSH domain) in both stores, plus `release_hex = -1` (land-seeded
# particles) in the counts store; 024f holds water-seeded drifters only.
# Those rows carry no geometry and are excluded from the maps — numerator
# *and* normalisation denominator — and their share is reported.

# %%
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(key_path.with_suffix(".json").read_text())["subbasin_id_to_name"].items()
}
subbasin_of_hex = key.set_index("hex_id")["helcom_subbasin"]

# The glob is loose and the regex selects: it pins exactly four year digits
# after the regime (so `surface` does not also take `surface_stokes_*`) and,
# for a member, the whole tag after the optional month.
if member:
    store_glob = f"HexAgg_survocc_r{hex_radius}m_{regime}_*_{member}.parquet"
    part_re = re.compile(
        rf"HexAgg_survocc_r{hex_radius}m_{regime}_(?P<year>\d{{4}})(?:_m(?P<month>\d{{2}}))?"
        rf"_{re.escape(member)}\.parquet$"
    )
else:
    store_glob = f"HexAgg_counts_r{hex_radius}m_{regime}_*.parquet"
    part_re = re.compile(rf"HexAgg_counts_r{hex_radius}m_{regime}_(?P<year>\d{{4}})\.parquet$")

partitions = []
years_not_read = set()
for f in sorted(store_root.glob(store_glob)):
    m = part_re.search(f.name)
    if m is None:
        continue
    if int(m["year"]) not in analysis_years:
        years_not_read.add(int(m["year"]))
        continue
    month = m.groupdict().get("month")
    partitions.append((f, int(m["year"]), None if month is None else int(month)))
if years_not_read:
    print(f"on disk but outside release_years, not read: {sorted(years_not_read)}")
if not partitions:
    raise FileNotFoundError(
        f"no partitions matching {store_glob!r} for release years {release_years} "
        f"at {store_root} — "
        + ("run 024f." if member else "run 024.")
    )
partition_years = sorted({year for _, year, _ in partitions})
# A whole-year partition next to monthly ones of the same year covers the
# same releases twice.
mixed = sorted(
    year for f, year, month in partitions
    if month is None and sum(y == year for _, y, _ in partitions) > 1
)
if mixed:
    raise ValueError(f"years {mixed} carry both a whole-year and monthly partitions")

parts = []
store_attrs = []
pooled_years = set()
weight_raw = 0.0
weight_sentinel = 0.0
for f, year, month in partitions:
    if month is not None and month not in season_months:
        continue
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
    season_doys = set(days[days.month.isin(season_months)].dayofyear.tolist())
    filters = [("release_doy", "in", season_doys), ("age_bin", "<", n_age_bins)]
    if member:
        df = pd.read_parquet(
            f, columns=["origin_subbasin", "target_hex", "age_bin", "surv"], filters=filters
        )
        store_attrs.append(df.attrs)
        df = df.rename(columns={"surv": "weight"})
        sentinel = df["target_hex"] < 0
    else:
        df = pd.read_parquet(
            f, columns=["release_hex", "target_hex", "age_bin", "n_obs"], filters=filters
        )
        df = df.rename(columns={"n_obs": "weight"})
        sentinel = (df["release_hex"] < 0) | (df["target_hex"] < 0)
    weight_raw += float(df["weight"].sum())
    weight_sentinel += float(df.loc[sentinel, "weight"].sum())
    df = df[~sentinel]
    if not member:
        df["origin_subbasin"] = (
            subbasin_of_hex.reindex(df["release_hex"].to_numpy())
            .fillna(-1).astype("int64").to_numpy()
        )
    parts.append(
        df.groupby(["origin_subbasin", "age_bin", "target_hex"], as_index=False)["weight"].sum()
    )
    if len(df):
        pooled_years.add(year)
    print(f"  {f.name}: {len(df):,} in-season valid rows")
    del df, sentinel
if not parts:
    raise ValueError(f"no {season} partitions among {len(partitions)} on disk")
binned = (
    pd.concat(parts, ignore_index=True)
    .groupby(["origin_subbasin", "age_bin", "target_hex"], as_index=False)["weight"]
    .sum()
)
del parts

# A pooled store must be one reduction window, at least as long as the one
# mapped here (the horizons are ≤ max_age_days, asserted above).
if member:
    # One window, binning and rate-model point across every pooled partition;
    # provenance (builder, git_sha, slurm_job_id) and the partition's own
    # release_year/release_month are not compared.
    pooled_keys = (
        "max_age_days", "age_bin_days", "w_c", "tau_strong_hours",
        "tau_calm_days", "delta", "trap_flat", "trap_wall", "band_m",
    )
    if any("max_age_days" not in a for a in store_attrs):
        raise ValueError(
            "survocc partitions without build attrs (old schema, written before "
            "max_age_days was recorded) — rebuild them with the current 024f reducer"
        )
    settings = {tuple(a[k] for k in pooled_keys) for a in store_attrs}
    assert len(settings) == 1, (
        f"pooled survocc partitions disagree on {pooled_keys}: "
        f"{sorted(settings)} — rebuild 024f"
    )
    store_build = dict(zip(pooled_keys, settings.pop()))
    store_max_age_days = store_build["max_age_days"]
    store_age_bin_days = store_build["age_bin_days"]
    assert store_age_bin_days == age_bin_days, (
        f"store age_bin_days {store_age_bin_days} != {age_bin_days}"
    )
    assert store_max_age_days >= max_age_days, (
        f"store window {store_max_age_days} d < max_age_days {max_age_days} d — "
        f"rebuild 024f or lower max_age_days"
    )
    print(f"024f window: {store_max_age_days} d in {store_age_bin_days} d bins")

# Years are taken after the season filter: a year whose partitions hold no
# release in this season (the partial last release year has no JJA or SON)
# is reported, not pooled.
pooled_years = sorted(pooled_years)
for y in sorted(set(partition_years) - set(pooled_years)):
    print(f"no {season} releases in {y}; dropped")
print(f"key: {len(key):,} hexes; season {season} (months {season_months}), "
      f"release years {pooled_years}")
print(f"{len(binned):,} (origin, age_bin, target) rows")
print(
    f"sentinel (-1) hexes excluded: {weight_sentinel:,.0f} of {weight_raw:,.0f} weight "
    f"({100 * weight_sentinel / weight_raw:.4f} %)"
)

# %% [markdown]
# # Resolve the origin subbasins to map
#
# `origin_subbasin` is -1 for releases whose hex sits outside every named
# HELCOM polygon (or is absent from the key); they carry no named origin, so
# they are reported here and dropped from the maps.

# %%
origin_codes_present = sorted(
    int(c) for c in binned["origin_subbasin"].unique() if c >= 0
)
name_to_id = {v: k for k, v in subbasin_id_to_name.items()}
unknown = [n for n in origin_subbasins_filter if n not in name_to_id]
if unknown:
    raise KeyError(f"unknown subbasin name(s) {unknown}; known: {sorted(name_to_id)}")
filter_ids = {name_to_id[n] for n in origin_subbasins_filter}
origin_codes = [c for c in origin_codes_present if not filter_ids or c in filter_ids]

print(
    f"origin subbasins to map ({len(origin_codes)}): "
    f"{', '.join(subbasin_id_to_name[c] for c in origin_codes)}"
)
print(
    f"excluded origin -1 (no named subbasin): "
    f"{binned.loc[binned['origin_subbasin'] < 0, 'weight'].sum():,.0f} weight"
)

# %% [markdown]
# # Relative density
#
# `to_hex_gdf` sums `weight` per target hex over a selected release set,
# joins the key geometry and derives both relative views. The
# normalisation is over whatever subset is passed in, so each horizon
# panel reads as that horizon's own share of particle-time.
#
# Hexes whose wet polygon has no water area (`water_area_m2` 0 — the hex
# is all land in the BSH coastline, yet a sub-grid trajectory position
# fell inside it) get no `dilution`; they still carry a `percentage`.

# %%
def to_hex_gdf(binned_df, key_df):
    """Per-target-hex relative density as a GeoDataFrame. ``percentage``
    is the share of the passed-in release set's particle-time; ``dilution``
    is that share per km² of hex water area."""
    gdf = (
        binned_df.groupby("target_hex", as_index=False)["weight"].sum()
        .merge(
            key_df[["hex_id", "geometry", "water_area_m2"]],
            left_on="target_hex", right_on="hex_id",
        )
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    gdf["percentage"] = 100 * gdf["weight"] / gdf["weight"].sum()
    gdf["dilution"] = gdf["percentage"] / (
        gdf["water_area_m2"].where(gdf["water_area_m2"] > 0) / 1e6
    )
    return gdf


COLUMN_LABEL = {
    "percentage": "share of particle-time (%)",
    "dilution": "dilution (% per km² water)",
}

# %% [markdown]
# # Map backdrop: coastline and the 3 m isobath
#
# The hex key decides the extent — no manual crop. The key tiles the full
# BSH model domain (North Sea included), so this shows the whole covered
# area; the dispersal cloud occupies wherever it actually reaches.
# Coastline is Natural Earth 10 m via cartopy's shapereader (cached
# locally), clipped to that extent once.
#
# The isobath comes from the BSH static H0 grids, *not* from the key's
# `mean_depth_m` — a hex mean over a 6 km cell cannot resolve a 3 m
# contour. H0 is the floor position (`H0 > 0` is always wet, see
# docs/h0_semantics.md), NaN over land, so `contour` masks land for free.
# Fine grid takes precedence where it covers (German Bight / western
# Baltic); the coarse grid is blanked inside the fine bbox so the two do
# not draw the same contour twice.

# %%
domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = key.total_bounds
domain_extent = [domain_lon_min, domain_lon_max, domain_lat_min, domain_lat_max]

h0_fine = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
).H0
h0_coarse = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_coarse/H0_file_coarse.nc"
).H0
h0_coarse_outside_fine = h0_coarse.where(
    ~(
        (h0_coarse.lon >= h0_fine.lon.min()) & (h0_coarse.lon <= h0_fine.lon.max())
        & (h0_coarse.lat >= h0_fine.lat.min()) & (h0_coarse.lat <= h0_fine.lat.max())
    )
)

coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max))

# Width/height of the extent in visually equal degrees — 1° lon at the
# extent's mid-latitude against 1° lat.
domain_aspect = (
    (domain_lon_max - domain_lon_min)
    * np.cos(np.radians(0.5 * (domain_lat_min + domain_lat_max)))
) / (domain_lat_max - domain_lat_min)


def map_grid(nrows, ncols):
    """A constrained-layout panel grid at full page width. Each map is
    width-limited, so the panel height follows from the page width and the
    extent aspect — no figsize arithmetic against the layout engine and no
    draw-measure-rescale loop. The colorbars sit *under* the maps (see
    ``hex_map``), so the fixed allowance is vertical only."""
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, layout="constrained")
    fig.set_size_inches(
        FIGURE_WIDTH_IN,
        nrows * ((FIGURE_WIDTH_IN / ncols) / domain_aspect + PANEL_DECORATION_IN),
    )
    return fig, axes


def percentage_norm(gdf):
    """Linear colour scale for the ``percentage`` view, saturating at the
    99th percentile. Occupancy is concentrated enough that a few hexes sit
    two decades above the bulk; scaling to the maximum paints every other
    hex the same dark colour and the map says nothing. The colorbar is
    drawn with an over-range arrow so the saturation is explicit."""
    if gdf.empty:
        return Normalize()
    return Normalize(vmin=0, vmax=gdf["percentage"].quantile(0.99))


def dilution_norm(gdf):
    """Log colour scale for the ``dilution`` view — the whole point of the
    view is the several decades the source concentration thins over."""
    finite = gdf["dilution"].where(gdf["dilution"] > 0).dropna()
    if not len(finite):
        return LogNorm()
    # Floor at the 1st percentile: a handful of near-empty hexes sit many
    # decades below the bulk and would otherwise eat most of the colour
    # range. The colorbar arrows mark both saturated ends.
    return LogNorm(vmin=finite.quantile(0.01), vmax=finite.max())


def extend_of(values, norm):
    """Colorbar ``extend`` that marks which ends of ``norm`` the values
    overrun."""
    over = bool(values.max() > norm.vmax)
    under = bool(norm.vmin is not None and values.min() < norm.vmin)
    return {
        (False, False): "neither", (True, False): "max",
        (False, True): "min", (True, True): "both",
    }[(over, under)]


def hex_map(
    gdf, column, ax, norm, title=None, cbar_label=True, colorbar=True,
    isobath=True, extent=domain_extent,
):
    """Draw one relative-density panel: hex polygons through ``norm``,
    coastline, and the H0 isobath, all on a plain (non-cartopy) lon/lat
    axis. Everything is already EPSG:4326, so no reprojection can push
    hexes, coastline and isobath out of register. An empty ``gdf`` skips
    the hex layer and still renders the backdrop, keeping the grid layout
    stable. ``colorbar=False`` leaves the colorbar to the caller (one
    shared bar for a whole figure); ``isobath=False`` drops the isobath;
    ``extent`` is ``[lon_min, lon_max, lat_min, lat_max]``."""
    if not gdf.empty and gdf[column].notna().any():
        colorbar_kwargs = {}
        if colorbar:
            # Colorbar as an inset *under* the map at axes-fraction width
            # 1.0, so it is exactly the map width at any figure size and the
            # panels keep the full page width instead of losing a quarter of
            # it to side-by-side colorbars and their tick labels.
            colorbar_kwargs = {
                "cax": ax.inset_axes([0.0, -0.09, 1.0, 0.035]),
                "legend_kwds": {
                    "label": COLUMN_LABEL[column] if cbar_label else "",
                    "orientation": "horizontal",
                    "extend": extend_of(gdf[column], norm),
                },
            }
        gdf.plot(
            ax=ax, column=column, norm=norm, legend=colorbar, **colorbar_kwargs,
            missing_kwds={"color": "none"},
            # Default polygon stroke is a black hairline that dominates at
            # Baltic zoom; matching edge to fill makes the hexes read as a
            # continuous field.
            edgecolor="face", linewidth=0.4, zorder=1,
        )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    # One styling override for the isobath: magenta at a hairline width is
    # the canonical anti-viridis, so the 3 m line stays legible over every
    # fill colour and clear of the black coastline.
    for h0 in (h0_fine, h0_coarse_outside_fine) if isobath else ():
        ax.contour(
            h0.lon, h0.lat, h0, levels=[isobath_m],
            colors="magenta", linewidths=0.3, zorder=3,
        )
    lon_min, lon_max, lat_min, lat_max = extent
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (lat_min + lat_max))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def export_geojson(gdf, name, **extra_columns):
    """Write one map product as GeoJSON for collaborators (EPSG:4326, hex
    geometry + the value columns)."""
    out = gdf[["hex_id", "geometry", "weight", "percentage", "dilution"]].assign(**extra_columns)
    path = export_dir / f"{name}.geojson"
    out.to_file(path, driver="GeoJSON")
    print(f"wrote {path}")


# %% [markdown]
# # Map per origin subbasin
#
# One figure per (origin subbasin, view) — the 026 horizon grid restricted
# to that subbasin's releases. Horizon `T` keeps every `age_bin` whose
# 10-day window starts before `T` (elapsed time in `[0, T)`), so the panels
# nest. Percentage and dilution go to separate figures rather than one
# eight-panel grid, so each map keeps a readable size at 180 mm page width.

# %%
ncols = 2
nrows = int(np.ceil(len(time_horizons_days) / ncols))

for code in origin_codes:
    name = subbasin_id_to_name[code]
    binned_sb = binned[binned["origin_subbasin"] == code]
    gdfs = {
        h: to_hex_gdf(binned_sb[binned_sb["age_bin"] * age_bin_days < h], key)
        for h in time_horizons_days
    }
    if all(g.empty for g in gdfs.values()):
        print(f"{name}: no occupied target hexes — skipping")
        continue
    for column, norm_of in (("percentage", percentage_norm), ("dilution", dilution_norm)):
        fig, axes = map_grid(nrows, ncols)
        fig.suptitle(f"origin: {name}")
        for j, (ax, h) in enumerate(zip(axes.flat, time_horizons_days)):
            hex_map(
                gdfs[h], column, ax, norm_of(gdfs[h]), title=f"{h} d",
                # The colour column is the same in every panel; naming it
                # once per figure, under the bottom row, keeps the grid
                # uncluttered.
                cbar_label=j >= len(time_horizons_days) - ncols,
            )
        for ax in axes.flat[len(time_horizons_days):]:
            ax.set_visible(False)
        fig_path = (
            figure_dir
            / f"OriginSubbasinTimeHorizonMaps_{regime}_{season}_r{hex_radius}m_{member_tag}"
            f"_{_slug(name)}_{column}.png"
        )
        # Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
        fig.savefig(fig_path, dpi=FIGURE_DPI)
        print(f"wrote {fig_path}")
        plt.show()
    for h in time_horizons_days:
        export_geojson(
            gdfs[h], f"026a_{regime}_{season}_{member_tag}_T{h}d_{_slug(name)}",
            origin_subbasin=name,
        )

# %% [markdown]
# # Overview: every origin on one page
#
# One figure per `overview_horizons_days` horizon: one dilution panel per
# origin subbasin in a 4-column grid, panel title = origin name, panels
# ordered geographically from the entrance (Kattegat) to the north-east.
# Each panel is still normalised over its own origin's releases, but every
# panel draws through **one** log colour scale under a single colorbar for
# the whole figure, so a colour means the same share of an origin's
# particle-time per km² in every panel and a compact cloud reads directly
# against a spread-out one.
#
# The scale spans four decades: top = 99.9th percentile of the dilution
# pooled over all panels (higher values saturate, `extend="max"`), floor =
# top / 10⁴. Hexes below the floor are not drawn, so each panel shows the
# footprint above 10⁻⁴ of the reference. Considered the per-origin maps'
# rule (1st percentile to maximum); rejected here, survival weights
# `exp(−A)` leave vanishing but nonzero values decades down, and a single
# hex with a sliver of water sets the maximum, so the pooled range ran over
# ~28 decades and every cloud rendered in one colour.
#
# Panels are cropped to `overview_extent` (Baltic plus the Kattegat /
# Skagerrak entrance) and carry no isobath: at a quarter of the page width
# the 3 m contour fills the Danish straits.

# %%
# West → east, south → north: entrance, western Baltic, Baltic Proper, gulfs,
# Bothnia. Origins not listed here are appended in code order.
OVERVIEW_ORDER = [
    "Kattegat", "The Sound", "Great Belt", "Kiel Bay", "Bay of Mecklenburg",
    "Arkona Basin", "Bornholm Basin", "Gdansk Basin", "Eastern Gotland Basin",
    "Western Gotland Basin", "Gulf of Riga", "Northern Baltic Proper",
    "Gulf of Finland", "Åland Sea", "Bothnian Sea", "The Quark",
]
overview_codes = [name_to_id[n] for n in OVERVIEW_ORDER if name_to_id[n] in origin_codes]
overview_codes += [c for c in origin_codes if c not in overview_codes]
overview_ncols = 4
overview_nrows = int(np.ceil(len(overview_codes) / overview_ncols))
OVERVIEW_DECADES = 4
ov_lon_min, ov_lon_max, ov_lat_min, ov_lat_max = overview_extent
overview_aspect = (
    (ov_lon_max - ov_lon_min) * np.cos(np.radians(0.5 * (ov_lat_min + ov_lat_max)))
) / (ov_lat_max - ov_lat_min)
# Vertical inches per row for the panel title, and once per figure for the
# shared colorbar with its label.
OVERVIEW_TITLE_IN = 0.25
OVERVIEW_COLORBAR_IN = 0.6

for h in overview_horizons_days:
    gdfs = {
        code: to_hex_gdf(
            binned[(binned["origin_subbasin"] == code) & (binned["age_bin"] * age_bin_days < h)],
            key,
        )
        for code in overview_codes
    }
    pooled = pd.concat(gdfs.values())
    top = pooled["dilution"].where(pooled["dilution"] > 0).quantile(0.999)
    norm = LogNorm(vmin=top / 10**OVERVIEW_DECADES, vmax=top)
    # Below-floor hexes are masked, not drawn at the floor colour.
    drawn = {code: g[g["dilution"] >= norm.vmin] for code, g in gdfs.items()}
    fig, axes = plt.subplots(
        nrows=overview_nrows, ncols=overview_ncols, squeeze=False, layout="constrained"
    )
    # Width-limited maps: the panel height follows from the page width and
    # the extent aspect, as in ``map_grid``.
    fig.set_size_inches(
        FIGURE_WIDTH_IN,
        overview_nrows
        * ((FIGURE_WIDTH_IN / overview_ncols) / overview_aspect + OVERVIEW_TITLE_IN)
        + OVERVIEW_COLORBAR_IN,
    )
    for ax, code in zip(axes.flat, overview_codes):
        hex_map(
            drawn[code], "dilution", ax, norm, colorbar=False, isobath=False,
            extent=overview_extent,
        )
        # Panel titles one step below the default size: at four maps per
        # 180 mm row the longest HELCOM names would run into the neighbours.
        ax.set_title(subbasin_id_to_name[code], fontsize="medium")
    for ax in axes.flat[len(overview_codes):]:
        ax.set_visible(False)
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm), ax=axes,
        orientation="horizontal", label=COLUMN_LABEL["dilution"],
        extend="max" if pooled["dilution"].max() > norm.vmax else "neither",
    )
    fig.suptitle(f"by origin subbasin, ages < {h} d")
    fig_path = (
        figure_dir
        / f"OriginSubbasinOverview_{regime}_{season}_r{hex_radius}m_{member_tag}_T{h}d_dilution.png"
    )
    # Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()
    # What the floor and the crop leave out, per origin: the share of the
    # origin's particle-time below the floor, and of the above-floor part
    # lying outside the extent (hex centroid).
    print(f"  {h} d: shared dilution scale {norm.vmin:.3g}–{norm.vmax:.3g} % per km²")
    for code in overview_codes:
        g = drawn[code]
        centre = g.geometry.representative_point()
        outside = (
            (centre.x < ov_lon_min) | (centre.x > ov_lon_max)
            | (centre.y < ov_lat_min) | (centre.y > ov_lat_max)
        )
        print(
            f"    {subbasin_id_to_name[code]}: {100 - g['percentage'].sum():.2f} % below floor, "
            f"{g.loc[outside, 'percentage'].sum():.3f} % above floor outside the extent"
            + (f" (westmost {centre.x.min():.2f}°E)" if len(g) else "")
        )

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, season={season} {season_months}, hex_radius={hex_radius} m")
print(f"member={member_tag}, age_bin_days={age_bin_days}, max_age_days={max_age_days}")
print(f"release years pooled: {pooled_years}")
print(f"sentinel (-1) share of raw weight: {100 * weight_sentinel / weight_raw:.4f} %")
print(f"origin subbasins mapped: {len(origin_codes)}")
longest = max(time_horizons_days)
for code in origin_codes:
    binned_sb = binned[
        (binned["origin_subbasin"] == code) & (binned["age_bin"] * age_bin_days < longest)
    ]
    gdf_sb = to_hex_gdf(binned_sb, key)
    print(
        f"  {subbasin_id_to_name[code]} (to {longest} d): "
        f"{len(gdf_sb):,} occupied target hexes, "
        f"sum(weight)={binned_sb['weight'].sum():,.0f}, "
        f"max percentage={gdf_sb['percentage'].max():.4f} %"
    )
