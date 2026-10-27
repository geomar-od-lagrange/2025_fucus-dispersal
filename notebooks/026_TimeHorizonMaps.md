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

# Time-horizon relative-density maps

Hex relative-density maps at successive elapsed-time horizons (e.g. 30,
60, 90, 180 days), showing how the dispersal cloud spreads and thins
over time, the median particle age per hex over the whole window, and
the water area that holds a given share of the cloud at
every age-bin edge (a dense series for 035). One regime and one release
season per run, pooled over the `release_years` range; no trajectory zarrs,
no Dask cluster.

Two stores feed this, selected by the `member` parameter, and the value
read from either is called `weight` from then on:

| `member` | store | `weight` |
|---|---|---|
| `""` | 024 counts (`HexAgg_counts_*`) | `n_obs`, every sample counts 1 |
| rate-model tag | 024f survocc (`HexAgg_survocc_*_<member>`) | `surv`, each sample weighted by its un-beached fraction `exp(−A)` |

A horizon `T` is **cumulative**: every `age_bin` whose window starts
before `T`, i.e. elapsed time in `[0, T)`. The map is therefore the
particle-time the cloud has accumulated by `T`, not a 10-day snapshot at
`T`.

Density is relative, as in 025:

- `percentage` — 100 · weight(hex) / Σ weight, normalised **per horizon**,
  so each panel shows that horizon's own share of particle-time. Linear
  colour scale, one per figure.
- `dilution` — that share per km² of hex water area (`water_area_m2`
  from the 024a key). The fixed decade scale shared by every dilution
  figure of the pipeline (`dilution_floor` to 100 % per km²).

```python
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, Normalize
from matplotlib.ticker import FuncFormatter
from shapely.geometry import box
from shapely.ops import unary_union
```

# Parameters

```python tags=["parameters"]
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
# Horizons of the area table (Exports/026/026_areas_*.csv), in days: same
# rules as the map horizons. "" = every age-bin edge up to max_age_days
# (10, 20, ..., 180 at the defaults), the dense series 035 plots.
area_horizons_days_csv = ""
# Shares of the weight whose smallest holding water area is reported per
# horizon (comma-separated, each in (0, 1]). 0.9 only: the 0.5 area is set by
# particles stalling in a few coastal hexes (docs/visualisations.md).
area_quantiles_csv = "0.9"

# Isobath to overlay on every panel, in metres of BSH H0 (floor position,
# see docs/h0_semantics.md). Fucus vesiculosus no longer grows below ~3 m.
isobath_m = 3.0

# Floor of the dilution colour scale, % per km², a whole decade. The same in
# every dilution figure (026, 026a, 034): chosen once from the paper figures'
# pooled data so that ≥ 99 % of every panel's weight lies above it
# (docs/visualisations.md). Hexes below it are not drawn.
dilution_floor = 1e-5

# Map extent (degrees), one for every map of the pipeline: the Baltic with
# the Belt Sea and Kattegat, North Sea cropped (docs/visualisations.md). Set
# all four to 0 to fall back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0
```

# Parse parameters

`SEASON_MONTHS` lives here rather than in the parameters cell because
papermill only injects primitives; `season` selects one entry.

```python
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
area_horizons_days = [int(x) for x in area_horizons_days_csv.split(",") if x] or list(
    range(age_bin_days, max_age_days + 1, age_bin_days)
)
assert max_age_days % age_bin_days == 0, (
    f"max_age_days {max_age_days} d is not a multiple of age_bin_days {age_bin_days} d"
)
for h in [*time_horizons_days, *area_horizons_days]:
    assert h > 0, f"horizon {h} d holds no age bin"
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d; "
        f"it would cut a bin in half"
    )
    assert h <= max_age_days, f"horizon {h} d is past max_age_days {max_age_days} d"
n_age_bins = max_age_days // age_bin_days
area_quantiles = [float(x) for x in area_quantiles_csv.split(",") if x]

member_tag = member or "unweighted"
print(f"regime={regime}, season={season}, member={member_tag}, window [0, {max_age_days}) d")

# Print-ready output: every saved figure is one full text-width (180 mm)
# page-width figure at 300 dpi, so panels land in the manuscript at their
# final size (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Vertical inches per row for the panel title, and once per figure for the
# shared colorbar under the grid with its label.
PANEL_TITLE_IN = 0.25
COLORBAR_IN = 0.6
# One page: a taller grid is capped here, its maps turn height-limited and
# the spare width becomes side margin.
MAX_FIGURE_HEIGHT_IN = 230 / 25.4

# Shared dilution scale: one class per decade from `dilution_floor` up to
# 100 % per km² — all of the material in one km² of water.
DILUTION_TOP = 1e2
dilution_exponents = np.arange(round(np.log10(dilution_floor)), round(np.log10(DILUTION_TOP)) + 1)
assert np.isclose(10.0 ** dilution_exponents[0], dilution_floor), (
    f"dilution_floor {dilution_floor} is not a whole decade"
)
# turbo, not the default colormap (docs/visualisations.md): dilution is
# multiplied by a per-site biomass downstream, so it has to read in absolute
# terms over seven decades — one high-contrast colour per decade. Below-floor
# hexes are not drawn, so the colorbar's extend="min" triangle stays empty.
dilution_cmap = plt.get_cmap("turbo").resampled(len(dilution_exponents) - 1).with_extremes(
    under="none"
)
dilution_norm = BoundaryNorm(10.0 ** dilution_exponents, dilution_cmap.N)
decade_format = FuncFormatter(lambda x, _: f"$10^{{{np.log10(x):.0f}}}$")
print(f"dilution scale: {dilution_floor:g}–{DILUTION_TOP:g} % per km², {dilution_cmap.N} decades")

figure_dir = output_root / "Figures" / "026"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "026"
export_dir.mkdir(parents=True, exist_ok=True)
```

# Read key + pool the store across years

Layout: flat files under ``output_root/HexAggregates/`` — the key
``HexAgg_key_r<radius>m.parquet``, the counts partitions
``HexAgg_counts_r<radius>m_<regime>_<year>.parquet`` (one per release
year) and the survocc partitions
``HexAgg_survocc_r<radius>m_<regime>_<year>[_mMM]_<member>.parquet`` (one
per release month, or one per year when 024f ran with
`release_month = 0`). Every partition is read and reduced immediately to
per-`(release_year, age_bin, target_hex)` sums — a counts partition is
O(10^8) rows, the reduction O(10^5). The season's release days and the
age window are pushed into the parquet read as filters; the doy set is
built per partition year so the month↔doy mapping is leap-correct, and a
monthly survocc partition outside the season is not opened at all.

Sentinels: `target_hex = -1` (no in-key position — the drifter left the
BSH domain) in both stores, plus `release_hex = -1` (land-seeded
particles) in the counts store; 024f holds water-seeded drifters only.
Those rows carry no geometry and are excluded from the maps — numerator
*and* normalisation denominator — and their share is reported.

```python
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)

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
weight_raw = 0.0
weight_sentinel = 0.0
for f, year, month in partitions:
    if month is not None and month not in season_months:
        continue
    days = pd.date_range(f"{year}-01-01", f"{year}-12-31")
    season_doys = set(days[days.month.isin(season_months)].dayofyear.tolist())
    filters = [("release_doy", "in", season_doys), ("age_bin", "<", n_age_bins)]
    if member:
        df = pd.read_parquet(f, columns=["target_hex", "age_bin", "surv"], filters=filters)
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
    parts.append(
        df[~sentinel]
        .groupby(["age_bin", "target_hex"], as_index=False)["weight"].sum()
        .assign(release_year=year)
    )
    print(f"  {f.name}: {int((~sentinel).sum()):,} in-season valid rows")
    del df, sentinel
if not parts:
    raise ValueError(f"no {season} partitions among {len(partitions)} on disk")
binned = (
    pd.concat(parts, ignore_index=True)
    .groupby(["release_year", "age_bin", "target_hex"], as_index=False)["weight"]
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
# must not enter the per-year areas.
years = sorted(int(y) for y in binned["release_year"].unique())
for y in sorted(set(partition_years) - set(years)):
    print(f"no {season} releases in {y}; dropped from the per-year areas")
print(f"key: {len(key):,} hexes; season {season} (months {season_months}), years {years}")
print(f"{len(binned):,} (release_year, age_bin, target) rows")
print(
    f"sentinel (-1) hexes excluded: {weight_sentinel:,.0f} of {weight_raw:,.0f} weight "
    f"({100 * weight_sentinel / weight_raw:.4f} %)"
)
```

# Relative density

`to_hex_gdf` sums `weight` per target hex over a selected release set,
joins the key geometry and derives both relative views. The
normalisation is over whatever subset is passed in, so each horizon
panel reads as that horizon's own share of particle-time.

Hexes whose wet polygon has no water area (`water_area_m2` 0 — the hex
is all land in the BSH coastline, yet a sub-grid trajectory position
fell inside it) get no `dilution`; they still carry a `percentage`.

```python
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
```

# Map backdrop: coastline and the 3 m isobath

The `extent_*` parameters decide the extent, shared by every map of the
pipeline. The key tiles the full BSH model domain (North Sea included), so
some weight can lie outside it; each map prints that share.

Coastline is the model's own: the outline of the BSH HBMnoku wet cells
(`coastline.geojson`, the polygons 024a takes the hex water area from),
merged fine-first as in 024a — fine polygons, plus the coarse ones outside
the fine-grid bbox. The outline is drawn minus a 0.001° strip along that
bbox, where the two grids' cells meet in open water.

The isobath comes from the BSH static H0 grids, *not* from the key's
`mean_depth_m` — a hex mean over a 6 km cell cannot resolve a 3 m
contour. H0 is the floor position (`H0 > 0` is always wet, see
docs/h0_semantics.md), NaN over land, so `contour` masks land for free.
Fine grid takes precedence where it covers (German Bight / western
Baltic); the coarse grid is blanked inside the fine bbox so the two do
not draw the same contour twice.

```python
if any((extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max)):
    domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = (
        extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max
    )
else:
    domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max = key.total_bounds

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

wet = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
# Fine bbox padded by half a fine cell, so every fine cell lies inside it.
fine_bbox = box(
    float(h0_fine.lon.min()) - 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.min()) - 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
    float(h0_fine.lon.max()) + 0.5 * abs(float(h0_fine.lon.diff("lon").mean())),
    float(h0_fine.lat.max()) + 0.5 * abs(float(h0_fine.lat.diff("lat").mean())),
)
wet_fine_first = unary_union([
    unary_union(wet.loc[wet["grid"] == "fine", "geometry"].tolist()),
    unary_union(wet.loc[wet["grid"] == "coarse", "geometry"].tolist()).difference(fine_bbox),
])
coast = (
    gpd.GeoSeries([wet_fine_first], crs="EPSG:4326").boundary
    .difference(fine_bbox.exterior.buffer(1e-3))
    .clip(box(domain_lon_min, domain_lat_min, domain_lon_max, domain_lat_max))
)

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
    draw-measure-rescale loop. One colorbar per figure sits *under* the
    grid (``shared_colorbar``), so the fixed allowance is vertical only."""
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=False, layout="constrained")
    fig.set_size_inches(
        FIGURE_WIDTH_IN,
        min(
            nrows * ((FIGURE_WIDTH_IN / ncols) / domain_aspect + PANEL_TITLE_IN) + COLORBAR_IN,
            MAX_FIGURE_HEIGHT_IN,
        ),
    )
    return fig, axes


def shared_colorbar(fig, axes, norm, cmap, label, extend, **kwargs):
    """One horizontal colorbar under the whole grid. pad override: the
    default 0.05 is a fraction of the whole grid's height, a blank band
    between the last row and the bar."""
    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), ax=axes, orientation="horizontal",
        label=label, extend=extend, pad=0.01, **kwargs,
    )


def share_outside_extent(gdf):
    """Share (%) of ``gdf``'s weight on hexes whose representative point lies
    outside the map extent — what the crop cuts off a map."""
    point = gdf.geometry.representative_point()
    inside = point.x.between(domain_lon_min, domain_lon_max) & point.y.between(
        domain_lat_min, domain_lat_max
    )
    return 100 * gdf.loc[~inside, "weight"].sum() / gdf["weight"].sum()


def hex_map(gdf, column, ax, norm, cmap=None, title=None):
    """Draw one relative-density panel: hex polygons through ``norm`` and
    ``cmap`` (None: the default colormap), coastline, and the H0 isobath,
    all on a plain (non-cartopy) lon/lat axis. Everything is already
    EPSG:4326, so no reprojection can push hexes, coastline and isobath out
    of register. An empty ``gdf`` skips the hex layer and still renders the
    backdrop, keeping the grid layout stable. The colorbar is the caller's
    (one per figure, ``shared_colorbar``)."""
    # Only hexes touching the extent: off-extent patches still count in
    # constrained_layout's bounding box and push the title off the page.
    gdf = gdf.cx[domain_lon_min:domain_lon_max, domain_lat_min:domain_lat_max]
    if not gdf.empty and gdf[column].notna().any():
        gdf.plot(
            ax=ax, column=column, norm=norm, cmap=cmap,
            missing_kwds={"color": "none"},
            # Default polygon stroke is a black hairline that dominates at
            # Baltic zoom; matching edge to fill makes the hexes read as a
            # continuous field.
            edgecolor="face", linewidth=0.4, zorder=1,
        )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    # One styling override for the isobath: magenta at a hairline width is
    # the canonical anti-viridis, so the 3 m line stays legible over every
    # fill colour. Drawn last at zorder 10, above every other artist: where
    # it runs along the coast it must sit on the black line, not under it.
    for h0 in (h0_fine, h0_coarse_outside_fine):
        ax.contour(
            h0.lon, h0.lat, h0, levels=[isobath_m],
            colors="magenta", linewidths=0.3, zorder=10,
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
    out = gdf[["hex_id", "geometry", "weight", "percentage", "dilution"]].assign(**extra_columns)
    path = export_dir / f"{name}.geojson"
    out.to_file(path, driver="GeoJSON")
    print(f"wrote {path}")
```

# Map per horizon

Horizon `T` keeps every `age_bin` whose 10-day window starts before `T`
— elapsed time in `[0, T)` — so the panels nest, each adding the
particle-time accumulated since the previous one. One figure per view
(percentage, dilution) rather than one eight-panel grid, so each map
keeps a readable size at 180 mm page width.

Each figure has one colour scale and one colorbar under the grid.
`percentage` is linear from 0 to the 99th percentile pooled over the
panels: occupancy is concentrated enough that a few hexes sit two decades
above the bulk, and scaling to the maximum would paint every other hex one
colour; the over-range arrow marks the saturation. `dilution` is the fixed
decade scale: hexes below `dilution_floor` are not drawn, and the share of
each panel's weight they hold is printed.

```python
gdfs = {
    h: to_hex_gdf(binned[binned["age_bin"] * age_bin_days < h], key)
    for h in time_horizons_days
}
for h in time_horizons_days:
    print(f"{h} d: {share_outside_extent(gdfs[h]):.3f} % of the mapped weight outside the map extent")

for h in time_horizons_days:
    dilution = gdfs[h]["dilution"]
    print(
        f"{h} d: {gdfs[h].loc[dilution < dilution_floor, 'percentage'].sum():.3f} % of the "
        f"weight below the dilution floor (not drawn), "
        f"{gdfs[h].loc[dilution.isna(), 'percentage'].sum():.3f} % on zero-water hexes, "
        f"{gdfs[h].loc[dilution > DILUTION_TOP, 'percentage'].sum():.3f} % above the top"
    )

ncols = 2
nrows = int(np.ceil(len(time_horizons_days) / ncols))
percentage_all = pd.concat([gdfs[h]["percentage"] for h in time_horizons_days])
VIEWS = {
    "percentage": {
        "norm": Normalize(vmin=0, vmax=percentage_all.quantile(0.99)),
        "cmap": None, "extend": "max", "drawn": lambda g: g,
    },
    "dilution": {
        "norm": dilution_norm, "cmap": dilution_cmap, "extend": "min",
        "drawn": lambda g: g[g["dilution"] >= dilution_floor],
        "colorbar_kwargs": {"format": decade_format},
    },
}
for column, view in VIEWS.items():
    fig, axes = map_grid(nrows, ncols)
    for ax, h in zip(axes.flat, time_horizons_days):
        hex_map(
            view["drawn"](gdfs[h]), column, ax, view["norm"], cmap=view["cmap"],
            title=f"{h} d",
        )
    for ax in axes.flat[len(time_horizons_days):]:
        ax.set_visible(False)
    shared_colorbar(
        fig, axes, view["norm"], view["cmap"], COLUMN_LABEL[column], view["extend"],
        **view.get("colorbar_kwargs", {}),
    )
    fig_path = (
        figure_dir
        / f"TimeHorizonMaps_{regime}_{season}_r{hex_radius}m_{member_tag}_{column}.png"
    )
    # Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

for h in time_horizons_days:
    export_geojson(gdfs[h], f"026_{regime}_{season}_{member_tag}_T{h}d")
```

# Median particle age per hex

The same per-`(age_bin, target_hex)` weights read the other way: per hex,
the median age of its particle-time over the whole window
`[0, max_age_days)`, pooled over the release years and the season — the
age at which the hex's weight, accumulated over `age_bin`, reaches half its
total. Ages are taken as uniform within a bin, so the median is
interpolated linearly inside the crossing bin:

    median_age = (b + (W/2 − C_{b−1}) / w_b) · age_bin_days,

with `w_b` the hex's weight in bin `b`, `C_b` its cumulative weight, `W`
the total and `b` the first bin with `C_b ≥ W/2`. The median rather than
the mean: a hex's age distribution need not be symmetric, and the median
reads directly as "half of the particle-time here is younger than this".

Young water sits where the releases are; old water where material arrives
late or lingers. A hex reached by a handful of samples gets an arbitrary
age, so the map floors on **share**: a hex is drawn only if its share of
the window's weight (`percentage`) is at or above the 5th percentile of the
per-hex shares. A dilution floor would let near-empty hexes with a tiny
water area through (dilution = share / area), which in the Elbe and Ems
estuaries are bright isolated old-age dots. Zero-water hexes carry no
dilution and stay out, as on the dilution maps.

Colour: seven 25-day classes of the default colormap, 0 to 175 d, with an
eighth for 175 d up to the window end (`extend="max"`) — fixed, so every
regime and season reads on the same classes. Exported as
`026_medianage_<regime>_<season>_<tag>.geojson`, mapped hexes only.

```python
def median_age_days(weight_by_bin):
    """Median age (d) per row of a (rows × age_bin) weight table: the age at
    which the row's cumulative weight reaches half its total, linear within
    the crossing bin (ages uniform within a bin)."""
    w = weight_by_bin.to_numpy(dtype=float)
    cum = np.cumsum(w, axis=1)
    half = 0.5 * cum[:, -1]
    # First bin whose cumulative weight reaches half the total.
    b = (cum < half[:, None]).sum(axis=1)
    rows = np.arange(len(w))
    below = cum[rows, b] - w[rows, b]
    return pd.Series(
        (b + (half - below) / w[rows, b]) * age_bin_days,
        index=weight_by_bin.index, name="median_age_days",
    )


hex_weight = binned.groupby("target_hex")["weight"].sum()
median_age = median_age_days(
    binned.groupby(["target_hex", "age_bin"])["weight"].sum()
    .unstack("age_bin", fill_value=0)
    .reindex(columns=range(n_age_bins), fill_value=0)
    .loc[hex_weight > 0]
)
gdf_age = to_hex_gdf(binned, key).join(median_age, on="target_hex")
gdf_age = gdf_age[gdf_age["dilution"] > 0]
# 5th percentile, not the 1st that floors the dilution colour scale: the
# isolated estuary hexes sit at the 0.2nd–4.3rd percentile of share (surface,
# ALL, 2016–2024, production member).
age_share_floor = gdf_age["percentage"].quantile(0.05)
gdf_age = gdf_age[gdf_age["percentage"] >= age_share_floor]

# The default colormap, discretised: seven 25 d classes plus the over class.
AGE_BOUNDS_DAYS = np.arange(0, 176, 25)
age_cmap = plt.get_cmap().resampled(len(AGE_BOUNDS_DAYS))
age_norm = BoundaryNorm(AGE_BOUNDS_DAYS, age_cmap.N, extend="max")

fig, axes = map_grid(1, 1)
hex_map(
    gdf_age, "median_age_days", axes[0, 0], age_norm, cmap=age_cmap,
    title=f"median particle age, ages < {max_age_days} d",
)
# The title names the quantity; the bar only needs the unit.
shared_colorbar(fig, axes, age_norm, age_cmap, "days", "max")
print(f"median age: {share_outside_extent(gdf_age):.3f} % of the mapped weight outside the map extent")
fig_path = figure_dir / f"MedianAge_{regime}_{season}_r{hex_radius}m_{member_tag}.png"
# Print-ready: full 180 mm page width at 300 dpi (see the parse cell).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

export_geojson(
    gdf_age, f"026_medianage_{regime}_{season}_{member_tag}",
    median_age_days=gdf_age["median_age_days"],
)

# Weighted median of the hex median ages: the hex age the median unit of
# particle-time sits in.
by_age = gdf_age.sort_values("median_age_days")
cum_weight = by_age["weight"].cumsum() / by_age["weight"].sum()
weighted_median_age = by_age["median_age_days"].iloc[int(np.searchsorted(cum_weight, 0.5))]
basin_median_age = median_age_days(
    binned.groupby("age_bin")["weight"].sum()
    .reindex(range(n_age_bins), fill_value=0)
    .to_frame().T
).iloc[0]
print(
    f"{regime}, {season}, {member_tag}: {len(gdf_age):,} of {len(median_age):,} hexes mapped "
    f"(share ≥ {age_share_floor:.3g} %), holding "
    f"{100 * gdf_age['weight'].sum() / hex_weight.sum():.3f} % of the weight"
)
print(
    f"  hex median age: weight-weighted median {weighted_median_age:.1f} d, "
    f"unweighted median {gdf_age['median_age_days'].median():.1f} d, "
    f"range {gdf_age['median_age_days'].min():.1f}–{gdf_age['median_age_days'].max():.1f} d; "
    f"basin-wide median age {basin_median_age:.1f} d"
)
```

# Water area holding a share of the cloud

For each horizon `T` and share `q`: the smallest water area that holds a
share `q` of the cumulative weight — hexes ranked by density
`weight / water_area_m2` (descending) and accumulated until the share
reaches `q`. The hex that crosses `q` counts only the fraction of its area
needed: a 6 km hex is 93.5 km² of water, so whole hexes would quantise
the small-T areas. Hexes with weight but no water area (all land in the
BSH coastline, yet reached by a sub-grid position) rank first at zero
area. `length_km = sqrt(area_km2)` is the side of the equal-area square —
a length scale, not a distance travelled.

`pooled` ranks the weights summed over every release year; each year is
the same reduction on that year's releases alone, so their min–max is the
interannual range (a pooled value need not lie inside it).

The table runs at `area_horizons_days` (default: every age-bin edge), not
the four map horizons, so 035 can draw it as a time series. Each scope's
per-`(target_hex, age_bin)` weights are accumulated over `age_bin` once;
horizon `T` is then the column of the last bin starting before `T`.

```python
water_area_km2 = key.set_index("hex_id")["water_area_m2"] / 1e6


def area_to_hold(weight_by_hex, quantiles):
    """Water area (km²) of the densest hexes that together hold each share
    ``q`` of ``weight_by_hex`` (a Series indexed by hex id), the hex that
    crosses ``q`` taken fractionally."""
    w = weight_by_hex[weight_by_hex > 0]
    a = water_area_km2.reindex(w.index).to_numpy()
    w = w.to_numpy(dtype=float)
    with np.errstate(divide="ignore"):
        density = w / a  # +inf where a hex has no water area: ranks first
    order = np.argsort(-density, kind="stable")
    share = w[order] / w.sum()
    area = a[order]
    cum_share = np.cumsum(share)
    cum_area = np.cumsum(area)
    out = []
    for q in quantiles:
        i = min(int(np.searchsorted(cum_share, q)), len(share) - 1)
        out.append(cum_area[i] - area[i] + area[i] * (q - (cum_share[i] - share[i])) / share[i])
    return out


def cumulative_by_hex(binned_df):
    """Weight per target hex (rows) accumulated over age bins (columns):
    column ``b`` holds elapsed time in ``[0, (b + 1) · age_bin_days)``."""
    return (
        binned_df.groupby(["target_hex", "age_bin"])["weight"].sum()
        .unstack("age_bin", fill_value=0)
        .reindex(columns=range(n_age_bins), fill_value=0)
        .cumsum(axis=1)
    )


scopes = {"pooled": cumulative_by_hex(binned)}
for y, g in binned.groupby("release_year"):
    scopes[str(y)] = cumulative_by_hex(g)

area_rows = []
for h in area_horizons_days:
    for scope, cumulative in scopes.items():
        weight_by_hex = cumulative[h // age_bin_days - 1]
        for q, a in zip(area_quantiles, area_to_hold(weight_by_hex, area_quantiles)):
            area_rows.append(
                {"horizon_days": h, "q": q, "scope": scope, "area_km2": a, "length_km": np.sqrt(a)}
            )
areas = pd.DataFrame(area_rows)
areas_path = export_dir / f"026_areas_{regime}_{season}_{member_tag}.csv"
areas.to_csv(areas_path, index=False)
print(f"wrote {areas_path}")

year_range = (
    areas[areas["scope"] != "pooled"]
    .groupby(["horizon_days", "q"])["area_km2"].agg(["min", "max"])
    .add_prefix("year_").add_suffix("_km2")
)
summary = (
    areas[areas["scope"] == "pooled"]
    .set_index(["horizon_days", "q"])[["area_km2", "length_km"]]
    .join(year_range)
)
print(f"{regime}, {season}, {member_tag}: pooled over {len(years)} release year(s) {years}")
print(summary.round(1).to_string())
```

# Validation prints

```python
print(f"regime={regime}, season={season} {season_months}, hex_radius={hex_radius} m")
print(f"member={member_tag}, age_bin_days={age_bin_days}, max_age_days={max_age_days}")
print(f"release years pooled: {years}")
print(f"sentinel (-1) share of raw weight: {100 * weight_sentinel / weight_raw:.4f} %")
for h in time_horizons_days:
    gdf_h = gdfs[h]
    print(
        f"  {h:>4} d: {len(gdf_h):,} occupied target hexes, "
        f"sum(weight)={gdf_h['weight'].sum():,.0f}, "
        f"max percentage={gdf_h['percentage'].max():.4f} %, "
        f"top-10 hexes={gdf_h['percentage'].nlargest(10).sum():.2f} %"
    )
```

# Per-target-hex occupancy quantile function

Quantile function of each hex's share of particle-time — one line per
horizon, log y-axis (the share spans decades). A lifting upper tail at
later horizons means a few hexes accumulate an ever larger share: the
coastal/convergence trapping expected from an advection-only run with no
diffusion or beaching kernel. Each horizon is normalised to 100 %, so a
rising upper tail is redistribution into fewer hexes, not more particles.

```python
q = np.linspace(0, 1, 31)
occupancy_quantiles = pd.DataFrame(
    {f"{h} d": gdfs[h]["percentage"].quantile(q).to_numpy() for h in time_horizons_days},
    index=pd.Index(q, name="quantile"),
)
ax = occupancy_quantiles.plot(logy=True)
ax.set_ylabel("share of particle-time per target hex (%)")
plt.show()
```
