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
# # Hex relative-density maps
#
# Relative particle-density maps on the hex-aggregated dispersal store
# built by notebooks 024a (key) and 024 (counts). Reads one key file plus
# every release-year counts partition of one regime; no trajectory zarrs,
# no Dask cluster.
#
# Density is reported in **relative** units, not raw particle counts:
#
# - `percentage` — 100 · n_obs(hex) / Σ n_obs over the selected release
#   set (regime, season, all release years pooled): the share of
#   particle-time the cloud spends in each hex. Linear colour scale.
# - `dilution` — `percentage` per km² of hex **water** area
#   (`water_area_m2` from the 024a key): how far the source concentration
#   has thinned, independent of how much of a hex is land. Log colour
#   scale.
#
# Figures:
#
# - **Baltic** — percentage + dilution over the whole basin
# - **German waters** — the same two views, viewport-clipped to
#   `de_extent` (one hex size for the whole basin, so the zoom is purely
#   visual)
# - **Per origin subbasin** — percentage, normalised per panel, so each
#   HELCOM release subbasin shows its own share

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
# Read root of the data twin (HELCOM polygons, BSH static H0).
data_root = "../data"
# Read root of the hex-aggregate stores built by 024a/024, and write root
# for figures and exports.
output_root = "../output"
# Which regime's counts partitions to read; one regime per run.
regime = "surface_stokes"
# Release season: DJF, MAM, JJA, SON or ALL (every month).
season = "ALL"
# Hex radius of the aggregate store (built by 024). One radius covers
# the whole BSH domain; the DE figure is a viewport clip, not a separate
# store.
hex_radius = 6000

# Baltic-wide map extent (degrees E / degrees N).
baltic_lon_min = 5
baltic_lon_max = 32
baltic_lat_min = 53
baltic_lat_max = 66

# DE-zoom map extent (degrees E / degrees N).
de_lon_min = 8
de_lon_max = 15
de_lat_min = 53.2
de_lat_max = 55.5

# Isobath to overlay on every panel, in metres of BSH H0 (floor position,
# see docs/h0_semantics.md). Fucus vesiculosus no longer grows below ~3 m.
isobath_m = 3.0

# Most panels a facet grid may carry before it is split across figures —
# at 180 mm page width more panels shrink each map below readable size.
max_panels_per_figure = 12

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
season_months = SEASON_MONTHS[season]

data_root = Path(data_root)
output_root = Path(output_root)
store_root = output_root / "HexAggregates"

# Print-ready output: every saved figure is one full text-width (180 mm)
# page-width figure at 300 dpi, so panels land in the manuscript at their
# final size. This is the deliberate figsize/dpi override for this
# notebook (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Vertical inches a panel needs on top of its map for the horizontal
# colorbar, the colorbar label and the panel title.
PANEL_DECORATION_IN = 0.62

figure_dir = output_root / "Figures" / "025"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "025"
export_dir.mkdir(parents=True, exist_ok=True)


def _slug(name):
    """Filename-safe slug of a HELCOM subbasin name. NFKD + ASCII fold
    transliterates accents (Å → A) before collapsing every non-alphanumeric
    run to a single underscore."""
    ascii_name = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^0-9A-Za-z]+", "_", ascii_name).strip("_")


# %% [markdown]
# # Read the key
#
# Layout: flat files under ``output_root/HexAggregates/`` —
# ``HexAgg_key_r<radius>m.parquet`` (+ a ``.json`` sidecar carrying the
# subbasin id→name map) and
# ``HexAgg_counts_r<radius>m_<regime>_<year>.parquet``.

# %%
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(key_path.with_suffix(".json").read_text())["subbasin_id_to_name"].items()
}
print(f"key: {len(key):,} hexes")

# %% [markdown]
# # Pool the counts across release years
#
# Every year partition of the regime is read and reduced immediately to
# per-`(origin_subbasin, target_hex)` sums — the store is O(10^8) rows per
# year, the reduction O(10^5). The season month set is pushed into the
# parquet read as a `release_doy` filter, built per partition year so the
# month↔doy mapping is leap-correct.
#
# Origin subbasin is the `helcom_subbasin` of each `release_hex` (024a
# attaches it as an integer id). The `-1` sentinel in either hex column
# marks land-seeded particles and obs that left the BSH domain; those rows
# carry no geometry and are excluded from the maps — from the numerator
# *and* from the normalisation denominator — and accounted for below.

# %%
_COUNTS_RE = re.compile(rf"HexAgg_counts_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$")
counts_files = sorted(
    store_root.glob(f"HexAgg_counts_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet")
)
if not counts_files:
    raise FileNotFoundError(
        f"no counts partitions for regime {regime!r} at {store_root} — run 024."
    )

subbasin_of_hex = key.set_index("hex_id")["helcom_subbasin"]

parts = []
n_obs_total_raw = 0
n_obs_sentinel = 0
for f in counts_files:
    year = int(_COUNTS_RE.search(f.name).group(1))
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
    season_doys = set(days[days.month.isin(season_months)].dayofyear.tolist())
    df = pd.read_parquet(
        f,
        columns=["release_hex", "target_hex", "n_obs"],
        filters=[("release_doy", "in", season_doys)],
    )
    n_obs_total_raw += int(df["n_obs"].sum())
    valid = (df["release_hex"] >= 0) & (df["target_hex"] >= 0)
    n_obs_sentinel += int(df.loc[~valid, "n_obs"].sum())
    df = df[valid]
    df["origin_subbasin"] = subbasin_of_hex.reindex(df["release_hex"].values).to_numpy()
    parts.append(
        df.groupby(["origin_subbasin", "target_hex"], as_index=False, dropna=False)["n_obs"].sum()
    )
    print(f"  {f.name}: {len(df):,} in-season valid rows")
    del df
counts = (
    pd.concat(parts, ignore_index=True)
    .groupby(["origin_subbasin", "target_hex"], as_index=False, dropna=False)["n_obs"]
    .sum()
)
del parts
print(f"season {season} (months {season_months}): {len(counts):,} (origin, target) pairs")
print(
    f"sentinel (-1) hexes excluded: {n_obs_sentinel:,} of {n_obs_total_raw:,} n_obs "
    f"({100 * n_obs_sentinel / n_obs_total_raw:.4f} %)"
)


# %% [markdown]
# # Relative density
#
# `to_hex_gdf` sums `n_obs` per target hex over a selected release set,
# joins the key geometry, and derives both relative views. The
# normalisation is over whatever subset is passed in, so a per-origin panel
# reads as that origin's own share of particle-time.
#
# Hexes whose wet polygon has no water area (`water_area_m2` 0 — the hex
# is all land in the BSH coastline, yet a sub-grid trajectory position
# fell inside it) get no `dilution`; they still carry a `percentage`.

# %%
def to_hex_gdf(counts_df, key_df):
    """Per-target-hex relative density as a GeoDataFrame. ``percentage``
    is the share of the passed-in release set's particle-time; ``dilution``
    is that share per km² of hex water area."""
    gdf = (
        counts_df.groupby("target_hex", as_index=False)["n_obs"].sum()
        .merge(
            key_df[["hex_id", "geometry", "water_area_m2"]],
            left_on="target_hex", right_on="hex_id",
        )
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )
    gdf["percentage"] = 100 * gdf["n_obs"] / gdf["n_obs"].sum()
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
# Coastline is Natural Earth 10 m via cartopy's shapereader (cached
# locally), clipped once per extent.
#
# The isobath comes from the BSH static H0 grids, *not* from the key's
# `mean_depth_m` — a hex mean over a 6 km cell cannot resolve a 3 m
# contour. H0 is the floor position (`H0 > 0` is always wet, see
# docs/h0_semantics.md), NaN over land, so `contour` masks land for free.
# Fine grid takes precedence where it covers (German Bight / western
# Baltic); the coarse grid is blanked inside the fine bbox so the two do
# not draw the same contour twice.

# %%
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

_coast_gdf = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
)
coast_baltic = _coast_gdf.clip(
    box(baltic_lon_min, baltic_lat_min, baltic_lon_max, baltic_lat_max)
)
coast_de = _coast_gdf.clip(box(de_lon_min, de_lat_min, de_lon_max, de_lat_max))

baltic_extent = [baltic_lon_min, baltic_lon_max, baltic_lat_min, baltic_lat_max]
de_extent = [de_lon_min, de_lon_max, de_lat_min, de_lat_max]


def extent_aspect(extent):
    """Width/height of an extent in visually equal degrees — 1° lon at the
    extent's mid-latitude against 1° lat."""
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    return (
        (lon_max_ - lon_min_) * np.cos(np.radians(0.5 * (lat_min_ + lat_max_)))
    ) / (lat_max_ - lat_min_)


def map_grid(nrows, ncols, extent):
    """A constrained-layout panel grid at full page width. Each map is
    width-limited, so the panel height follows from the page width and the
    extent aspect — no figsize arithmetic against the layout engine and no
    draw-measure-rescale loop. The colorbars sit *under* the maps (see
    ``hex_map``), so the fixed allowance is vertical only."""
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, layout="constrained")
    fig.set_size_inches(
        FIGURE_WIDTH_IN,
        nrows * ((FIGURE_WIDTH_IN / ncols) / extent_aspect(extent) + PANEL_DECORATION_IN),
    )
    return fig, axes


def percentage_norm(gdf):
    """Linear colour scale for the ``percentage`` view, saturating at the
    99th percentile. Occupancy is concentrated enough that a few hexes sit
    two decades above the bulk; scaling to the maximum paints every other
    hex the same dark colour and the map says nothing. The colorbar is
    drawn with an over-range arrow so the saturation is explicit."""
    return Normalize(vmax=gdf["percentage"].quantile(0.99)) if not gdf.empty else Normalize()


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


def hex_map(gdf, column, ax, extent, norm, coast, title=None, cbar_label=True):
    """Draw one relative-density panel: hex polygons through ``norm``,
    coastline, and the H0 isobath, all on a plain (non-cartopy) lon/lat
    axis. Everything is already EPSG:4326, so no reprojection can push
    hexes, coastline and isobath out of register. An empty ``gdf`` (a
    subbasin with no releases) skips the hex layer and still renders the
    backdrop, keeping the grid layout stable."""
    if not gdf.empty and gdf[column].notna().any():
        # Colorbar as an inset *under* the map at axes-fraction width 1.0,
        # so it is exactly the map width at any figure size and the panels
        # keep the full page width instead of losing a quarter of it to
        # side-by-side colorbars and their tick labels.
        cax = ax.inset_axes([0.0, -0.09, 1.0, 0.035])
        over = bool(gdf[column].max() > norm.vmax)
        under = bool(gdf[column].min() < norm.vmin)
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
    lon_min_, lon_max_, lat_min_, lat_max_ = extent
    ax.set_xlim(lon_min_, lon_max_)
    ax.set_ylim(lat_min_, lat_max_)
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (lat_min_ + lat_max_))))
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
# # Baltic-wide relative density

# %%
gdf_full = to_hex_gdf(counts, key)

fig, axes = map_grid(1, 2, baltic_extent)
hex_map(
    gdf_full, "percentage", axes[0, 0], baltic_extent,
    percentage_norm(gdf_full), coast_baltic, title="percentage",
)
hex_map(
    gdf_full, "dilution", axes[0, 1], baltic_extent,
    dilution_norm(gdf_full), coast_baltic, title="dilution",
)
fig_path = figure_dir / f"HexRelativeDensity_{regime}_{season}_r{hex_radius}m.png"
# Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

export_geojson(gdf_full, f"025_{regime}_{season}")

# %% [markdown]
# # German waters
#
# Same data and same normalisation as the Baltic figure — only the
# viewport changes, so the percentages still read as shares of the whole
# basin's particle-time.

# %%
fig, axes = map_grid(1, 2, de_extent)
hex_map(
    gdf_full, "percentage", axes[0, 0], de_extent,
    percentage_norm(gdf_full), coast_de, title="percentage",
)
hex_map(
    gdf_full, "dilution", axes[0, 1], de_extent,
    dilution_norm(gdf_full), coast_de, title="dilution",
)
fig_path = figure_dir / f"HexRelativeDensityDE_{regime}_{season}_r{hex_radius}m.png"
# Print-ready: full 180 mm page width at 300 dpi.
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Per origin subbasin
#
# One panel per HELCOM subbasin that seeds releases, each normalised over
# its own releases — the panel answers "where does *this* source's
# particle-time go?", not "how big is this source". Grids are split into
# figures of at most `max_panels_per_figure` panels so each map keeps a
# readable size at 180 mm page width.

# %%
origin_codes = sorted(
    int(c) for c in counts["origin_subbasin"].dropna().unique() if c >= 0
)
print(
    f"origin subbasins with releases ({len(origin_codes)}): "
    f"{', '.join(subbasin_id_to_name[c] for c in origin_codes)}"
)

origin_gdfs = {
    code: to_hex_gdf(counts[counts["origin_subbasin"] == code], key)
    for code in origin_codes
}
ncols = 3
for i in range(0, len(origin_codes), max_panels_per_figure):
    chunk = origin_codes[i : i + max_panels_per_figure]
    nrows = int(np.ceil(len(chunk) / ncols))
    fig, axes = map_grid(nrows, ncols, baltic_extent)
    for j, (ax, code) in enumerate(zip(axes.flat, chunk)):
        gdf_sb = origin_gdfs[code]
        hex_map(
            gdf_sb, "percentage", ax, baltic_extent,
            percentage_norm(gdf_sb), coast_baltic, title=subbasin_id_to_name[code],
            # The colour column is the same in every panel; naming it once
            # per figure, under the bottom row, keeps the grid uncluttered.
            cbar_label=j >= len(chunk) - ncols,
        )
    for ax in axes.flat[len(chunk):]:
        ax.set_visible(False)
    part = i // max_panels_per_figure + 1
    fig_path = (
        figure_dir
        / f"HexRelativeDensityByOrigin_{regime}_{season}_r{hex_radius}m_part{part}.png"
    )
    # Print-ready: full 180 mm page width at 300 dpi.
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

for code in origin_codes:
    export_geojson(
        origin_gdfs[code],
        f"025_{regime}_{season}_{_slug(subbasin_id_to_name[code])}",
        origin_subbasin=subbasin_id_to_name[code],
    )

# %% [markdown]
# # Validation prints

# %%
print(f"regime={regime}, season={season} {season_months}, hex_radius={hex_radius} m")
print(f"release-year partitions pooled: {len(counts_files)}")
print(f"  mapped n_obs: {int(gdf_full['n_obs'].sum()):,}")
print(f"  occupied target hexes: {len(gdf_full):,}")
print(
    f"  sentinel (-1) share of raw n_obs: "
    f"{100 * n_obs_sentinel / n_obs_total_raw:.4f} %"
)
print(
    f"  hexes without water area (no dilution): "
    f"{int(gdf_full['dilution'].isna().sum()):,}"
)
print(f"  max percentage in one hex: {gdf_full['percentage'].max():.4f} %")
print(
    f"  top-10 hexes carry "
    f"{gdf_full['percentage'].nlargest(10).sum():.2f} % of particle-time"
)
print(f"  max dilution: {gdf_full['dilution'].max():.4g} % per km²")
print(f"  origin subbasins mapped: {len(origin_codes)}")
