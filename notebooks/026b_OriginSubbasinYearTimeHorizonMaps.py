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
# # Time-horizon relative-density maps per origin subbasin and release year
#
# 026a split a second time, by **release year**: for each (origin HELCOM
# subbasin, year) pair the same elapsed-time-horizon relative-density maps,
# restricted to trajectories released from that subbasin in that year. This
# exposes the interannual variability 026a hides by pooling. Target hexes
# still range over the whole domain; only the release side is filtered.
#
# Origin subbasin is the `helcom_subbasin` of each `release_hex` (024a
# attaches it as an integer id; the id→name map lives in the key's JSON
# sidecar). Releases outside any named subbasin (`_outside`, id -1) and
# land-seeded releases (`release_hex` -1) carry no origin subbasin and are
# reported, not mapped.
#
# Each origin is normalised over **its own** releases across **all** years,
# so a year's panels sum to that year's share of the origin's four-year
# particle-time: a faint year reads as a genuinely weaker year rather than
# a rescaling artefact, and the four years' panels add up to 100 % at a
# given horizon. The colour scale is likewise shared across an origin's
# years and horizons, and independent between origins — read each origin's
# figures on their own colorbar.
#
# A horizon `T` is **cumulative**: every `age_bin` whose window starts
# before `T`, i.e. elapsed time in `[0, T)`. The map is therefore the
# particle-time the cloud has accumulated by `T`, not a 10-day snapshot at
# `T`.
#
# Density is relative, as in 025:
#
# - `percentage` — 100 · n_obs(hex) / Σ n_obs, normalised **per horizon**,
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
# Which regime's counts partitions to read; one regime per run.
regime = "surface_stokes"
# Release season: DJF, MAM, JJA, SON or ALL (every month).
season = "ALL"
# Hex radius of the store (built by 024a/024). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the counts store (must match 024). Horizons must
# be whole multiples of this.
age_bin_days = 10
# Elapsed-time horizons to map, in days. Each must be a multiple of
# age_bin_days.
time_horizons_days_csv = "10,20,50,100"

# Origin subbasins to map, by HELCOM name (comma-separated). Empty ⇒ every
# named subbasin that seeds releases.
origin_subbasins_csv = ""
# Release years to map (comma-separated). Empty ⇒ every year present.
release_years_csv = ""

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

output_root = Path(output_root)
data_root = Path(data_root)
store_root = output_root / "HexAggregates"
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would cut a bin in half"
    )

# Print-ready output: every saved figure is one full text-width (180 mm)
# page-width figure at 300 dpi, so panels land in the manuscript at their
# final size (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Vertical inches a panel needs on top of its map for the horizontal
# colorbar, the colorbar label and the panel title.
PANEL_DECORATION_IN = 0.62

origin_subbasins_filter = [s.strip() for s in origin_subbasins_csv.split(",") if s.strip()]
release_years_filter = [int(y) for y in release_years_csv.split(",") if y]

figure_dir = output_root / "Figures" / "026b"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "026b"
export_dir.mkdir(parents=True, exist_ok=True)


def _slug(name):
    """Filename-safe slug of a HELCOM subbasin name. NFKD + ASCII fold
    transliterates accents (Å → A) before collapsing every non-alphanumeric
    run to a single underscore."""
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^0-9A-Za-z]+", "_", ascii_name).strip("_")

# %% [markdown]
# # Read key + pool counts across years
#
# Layout: flat files under ``output_root/HexAggregates/`` —
# ``HexAgg_key_r<radius>m.parquet`` and
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``. Every year
# partition of the regime is read and reduced immediately to
# per-`(release_year, origin_subbasin, age_bin, target_hex)` sums — the store is O(10^8)
# rows per year, the reduction O(10^6). The season month set is pushed into
# the parquet read as a `release_doy` filter, built per partition year so the
# month↔doy mapping is leap-correct.
#
# The `-1` sentinel in either hex column marks land-seeded particles and
# obs that left the BSH domain; those rows carry no geometry and are
# excluded from the maps — numerator *and* normalisation denominator.

# %%
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(key_path.with_suffix(".json").read_text())["subbasin_id_to_name"].items()
}
subbasin_of_hex = key.set_index("hex_id")["helcom_subbasin"]

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
n_obs_total_raw = 0
n_obs_sentinel = 0
for f in counts_files:
    year = int(_COUNTS_RE.search(f.name).group(1))
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
    season_doys = set(days[days.month.isin(season_months)].dayofyear.tolist())
    df = pd.read_parquet(
        f,
        columns=["release_hex", "target_hex", "age_bin", "n_obs"],
        filters=[("release_doy", "in", season_doys)],
    )
    n_obs_total_raw += int(df["n_obs"].sum())
    valid = (df["release_hex"] >= 0) & (df["target_hex"] >= 0)
    n_obs_sentinel += int(df.loc[~valid, "n_obs"].sum())
    df = df[valid]
    df["origin_subbasin"] = subbasin_of_hex.reindex(df["release_hex"].values).to_numpy()
    part = df.groupby(
        ["origin_subbasin", "age_bin", "target_hex"], as_index=False, dropna=False
    )["n_obs"].sum()
    part["release_year"] = year
    parts.append(part)
    print(f"  {f.name}: {len(df):,} in-season valid rows")
    del df
# One partition per year, so the concat needs no second groupby.
counts = pd.concat(parts, ignore_index=True)
del parts
print(f"key: {len(key):,} hexes; pooled {len(counts_files)} year partition(s)")
print(
    f"season {season} (months {season_months}): "
    f"{len(counts):,} (year, origin, age_bin, target) rows"
)
print(
    f"sentinel (-1) hexes excluded: {n_obs_sentinel:,} of {n_obs_total_raw:,} n_obs "
    f"({100 * n_obs_sentinel / n_obs_total_raw:.4f} %)"
)

# %% [markdown]
# # Resolve the origin subbasins to map
#
# `origin_subbasin` is NaN for land-seeded releases (absent from the key)
# and -1 for releases whose hex centroid sits outside every named HELCOM
# polygon; neither carries a named origin, so both are reported here and
# dropped from the maps.

# %%
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

print(f"release years to map ({len(years)}): {years}")
print(
    f"origin subbasins to map ({len(origin_codes)}): "
    f"{', '.join(subbasin_id_to_name[c] for c in origin_codes)}"
)
print(
    f"excluded origins: "
    f"{int(counts.loc[counts['origin_subbasin'] == -1, 'n_obs'].sum()):,} n_obs outside any "
    f"named subbasin (_outside), "
    f"{int(counts.loc[counts['origin_subbasin'].isna(), 'n_obs'].sum()):,} n_obs from "
    f"release hexes absent from the key"
)

# %% [markdown]
# # Relative density
#
# `to_hex_gdf` sums `n_obs` per target hex over a selected release set,
# joins the key geometry and derives both relative views. The
# normalisation is over whatever subset is passed in, so each horizon
# panel reads as that horizon's own share of particle-time.
#
# Hexes whose wet polygon has no water area (`water_area_m2` 0 — the hex
# is all land in the BSH coastline, yet a sub-grid trajectory position
# fell inside it) get no `dilution`; they still carry a `percentage`.

# %%
def to_hex_gdf(counts_df, key_df, total=None):
    """Per-target-hex relative density as a GeoDataFrame. ``percentage``
    is the share of ``total`` particle-time (default: the passed-in release
    set's own sum); ``dilution`` is that share per km² of hex water area.
    A year panel passes its origin's all-year total, so the four years are
    magnitude-comparable instead of each renormalising to 100 %."""
    gdf = (
        counts_df.groupby("target_hex", as_index=False)["n_obs"].sum()
        .merge(
            key_df[["hex_id", "geometry", "water_area_m2"]],
            left_on="target_hex", right_on="hex_id",
        )
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    gdf["percentage"] = 100 * gdf["n_obs"] / (gdf["n_obs"].sum() if total is None else total)
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


def hex_map(gdf, column, ax, norm, title=None, cbar_label=True):
    """Draw one relative-density panel: hex polygons through ``norm``,
    coastline, and the H0 isobath, all on a plain (non-cartopy) lon/lat
    axis. Everything is already EPSG:4326, so no reprojection can push
    hexes, coastline and isobath out of register. An empty ``gdf`` skips
    the hex layer and still renders the backdrop, keeping the grid layout
    stable."""
    if not gdf.empty and gdf[column].notna().any():
        # Colorbar as an inset *under* the map at axes-fraction width 1.0,
        # so it is exactly the map width at any figure size and the panels
        # keep the full page width instead of losing a quarter of it to
        # side-by-side colorbars and their tick labels.
        cax = ax.inset_axes([0.0, -0.09, 1.0, 0.035])
        over = bool(gdf[column].max() > norm.vmax)
        under = bool(norm.vmin is not None and gdf[column].min() < norm.vmin)
        gdf.plot(
            ax=ax, column=column, norm=norm, legend=True, cax=cax,
            legend_kwds={
                "label": COLUMN_LABEL[column] if cbar_label else "",
                "orientation": "horizontal",
                "extend": {
                    (False, False): "neither", (True, False): "max",
                    (False, True): "min", (True, True): "both",
                }[(over, under)],
            },
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
    for h0 in (h0_fine, h0_coarse_outside_fine):
        ax.contour(
            h0.lon, h0.lat, h0, levels=[isobath_m],
            colors="magenta", linewidths=0.3, zorder=3,
        )
    ax.set_xlim(domain_lon_min, domain_lon_max)
    ax.set_ylim(domain_lat_min, domain_lat_max)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (domain_lat_min + domain_lat_max))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def export_geojson(gdf, name, **extra_columns):
    """Write one map product as GeoJSON for collaborators (EPSG:4326, hex
    geometry + the value columns)."""
    out = gdf[["hex_id", "geometry", "n_obs", "percentage", "dilution"]].assign(**extra_columns)
    path = export_dir / f"{name}.geojson"
    out.to_file(path, driver="GeoJSON")
    print(f"wrote {path}")


# %% [markdown]
# # Map per origin subbasin and year
#
# Outer loop over origin subbasins, inner loop over years: one figure per
# (origin, year, view), the 026a horizon grid restricted to one release
# year. Horizon `T` keeps every `age_bin` whose 10-day window starts before
# `T` (elapsed time in `[0, T)`), so the panels nest. Every panel of an
# origin is normalised by, and colour-scaled against, that origin's
# all-year particle-time at the same horizon, so the years are directly
# comparable.
#
# GeoJSON is exported at the longest horizon only: the product's axis here
# is the release year, and 026a already exports the horizon sweep.

# %%
ncols = 2
nrows = int(np.ceil(len(time_horizons_days) / ncols))
longest = max(time_horizons_days)

for code in origin_codes:
    name = subbasin_id_to_name[code]
    counts_sb = counts[counts["origin_subbasin"] == code]
    # The origin's all-year particle-time per horizon is the denominator
    # every one of its year panels is measured against.
    totals = {
        h: counts_sb.loc[counts_sb["age_bin"] * age_bin_days < h, "n_obs"].sum()
        for h in time_horizons_days
    }
    gdfs = {
        (y, h): to_hex_gdf(
            counts_sb[
                (counts_sb["release_year"] == y) & (counts_sb["age_bin"] * age_bin_days < h)
            ],
            key,
            total=totals[h],
        )
        for y in years
        for h in time_horizons_days
    }
    nonempty = [g for g in gdfs.values() if not g.empty]
    if not nonempty:
        print(f"{name}: no occupied target hexes — skipping")
        continue
    pooled = pd.concat(nonempty, ignore_index=True)
    for column, norm_of in (("percentage", percentage_norm), ("dilution", dilution_norm)):
        # One norm per (origin, view), spanning the origin's years and
        # horizons, so colour differences across years are signal.
        norm = norm_of(pooled)
        for y in years:
            if all(gdfs[(y, h)].empty for h in time_horizons_days):
                continue
            fig, axes = map_grid(nrows, ncols)
            fig.suptitle(f"origin: {name} — {y}")
            for j, (ax, h) in enumerate(zip(axes.flat, time_horizons_days)):
                hex_map(
                    gdfs[(y, h)], column, ax, norm, title=f"{h} d",
                    # The colour column is the same in every panel; naming
                    # it once per figure, under the bottom row, keeps the
                    # grid uncluttered.
                    cbar_label=j >= len(time_horizons_days) - ncols,
                )
            for ax in axes.flat[len(time_horizons_days):]:
                ax.set_visible(False)
            fig_path = (
                figure_dir
                / f"OriginSubbasinYearTimeHorizonMaps_{regime}_{season}_r{hex_radius}m_{_slug(name)}_{y}_{column}.png"
            )
            # Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
            fig.savefig(fig_path, dpi=FIGURE_DPI)
            print(f"wrote {fig_path}")
            plt.show()
    for y in years:
        if gdfs[(y, longest)].empty:
            continue
        export_geojson(
            gdfs[(y, longest)],
            f"026b_{regime}_{season}_T{longest}d_{_slug(name)}_{y}",
            origin_subbasin=name, release_year=y,
        )

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, season={season} {season_months}, hex_radius={hex_radius} m")
print(f"age_bin_days={age_bin_days}, release-year partitions pooled: {len(counts_files)}")
print(f"sentinel (-1) share of raw n_obs: {100 * n_obs_sentinel / n_obs_total_raw:.4f} %")
print(f"origin subbasins mapped: {len(origin_codes)}; years: {years}")
for code in origin_codes:
    in_horizon = counts[
        (counts["origin_subbasin"] == code) & (counts["age_bin"] * age_bin_days < longest)
    ]
    by_year = in_horizon.groupby("release_year")["n_obs"].sum()
    total = by_year.sum()
    per_year = ", ".join(
        f"{y}={100 * by_year.get(y, 0) / total:.1f} %" for y in years
    )
    print(f"  {subbasin_id_to_name[code]} (share of its particle-time to {longest} d): {per_year}")
