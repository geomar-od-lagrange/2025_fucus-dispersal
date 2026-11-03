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

# Beaching maps

Lightweight parquet-only consumer of the beaching store built by
`024e_BuildBeaching` (+ the `024a` key) — no trajectory zarrs, no Dask.
One regime and one release `season` per run; the season's monthly
partitions are pooled across the `release_years` range.

The rate model is a two-state hazard — a background rate in the near-shore
band plus a strong-wave rate switched on where the onshore Stokes forcing
exceeds a threshold (`docs/beaching.md`). Its settings are baked into the
reducer's opaque `member` tag, which selects the partitions read here and
labels every figure; this notebook never parses the tag. Draws:

1. **Where stranded, and by origin** — one figure, two maps side by side:
   the stranded weight per `beach_hex` as a percentage of the released
   water-seeded drifters (log scale), and per release hex the percentage of
   its own drifters stranded by `origin_horizon_days` (linear 0–100), so
   sink and source read off one figure. HELCOM subbasin borders on every
   map; the subbasin codes, their key and the release sites are on `033`'s
   reference map.
2. **Beaching age horizons** — cumulative where-stranded density for
   strandings at ages < T (the beaching analogue of `026`'s horizons).
3. **Afloat by origin subbasin** — per HELCOM origin subbasin, the
   percentage of its released drifters not yet stranded at each horizon,
   pooled and with the per-release-year min–max; table + CSV.
4. **Cumulative beached percentage vs. age** — the stranding time course,
   with the interannual min–max as a band around the pooled curve, exported
   as CSV for `035`'s combined time-series figure. Maps stay pooled; only
   the scalar-per-horizon curve carries the spread.
5. **Displacement at stranding** — stranded share over `disp_bin`, the
   crow-flies displacement from release at the deposit step. Short-travel
   strandings are a real outcome under a threshold rate model, never a
   mask.

Every percentage is of the **released water-seeded drifters** of the
pooled partitions (their summed `weight`: deposits + never-beached residual
per release hex sum to its released drifter count), so the maps, the curve
and the table share one denominator.

After 1. sits the number the source↔sink pair exists for: the **substrate
contrast between source and sink**, the shoreline flat fraction where the
release sites are against the flat fraction where the weight lands, under
this member and under a substrate-inert `contrast_member`.

Figures are print-ready (full page width, 300 dpi); the per-horizon
where-stranded field is exported as GeoJSON and the stranding curve as CSV.

```python
import itertools
import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm, Normalize
from pyproj import Transformer
from shapely.geometry import box
from shapely.ops import unary_union

from hextraj import HexProj
```

# Parameters

```python tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Read root of the data twin: the Fucus release sites 010 seeds from (for the
# substrate contrast), the HELCOM subbasin polygons whose borders every map
# draws, the BSH wet-cell coastline every map draws and the shoreline-class
# table the substrate contrast reads.
data_root = "../data"
# Which regime's beaching partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024e). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the beaching store (must match 024e; asserted
# against the partitions' attrs). Horizons must be whole multiples of this.
age_bin_days = 10
# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year. Selects
# which of the monthly `_mMM` partitions 024e wrote are read and pooled.
season = "ALL"
# Release years pooled: "first-last" (inclusive) or one year ("2019").
# Partitions of other years stay on disk unopened (2025 is partial).
release_years = "2016-2024"

# Rate-model member: the opaque tag the reducer put in the store filename
# (e.g. "step_wc0p1_ts3_tcinf_tf1_tw0"). Selects the partitions and tags
# the figures.
member = "step_wc0p1_ts3_tcinf_tf1_tw0"
# Second member, read only by the source/sink substrate block: a member whose
# trap term is degenerate (trap_flat = trap_wall = 1), so substrate does not
# enter its rate at all and its deposit is where the drift alone puts the
# material. "" skips the comparison.
contrast_member = "step_wc0p1_ts3_tcinf"
# Displacement bin width (km) of the store's `disp_bin` axis (must match
# 024e; asserted against the partitions' attrs).
disp_bin_km = 10.0

# Cumulative beaching age horizons (d) for the horizon maps and the afloat
# table: strandings at ages < T. Each must be a multiple of age_bin_days and
# ≤ the store's max_age_days.
time_horizons_days_csv = "30,60,90,180"
# Age horizon (d) of the by-origin map beside the where-stranded map: per
# release hex, the share of its drifters stranded at ages < T. Same rules as
# the horizons above.
origin_horizon_days = 180

# Map extent (degrees). The 024a key tiles the whole BSH domain including the
# North Sea, which is empty for Fucus and costs ~35% of panel height; cropping
# to the Baltic proper buys ~1.35x px/km for free. Set all four to 0 to fall
# back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0
```

# Parse parameters

The season → month mapping lives here rather than in the parameters cell:
papermill only injects primitives, and the mapping is a fixed convention,
not a knob.

```python
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}

output_root = Path(output_root)
data_root = Path(data_root)
# The file 010 seeds from: the REDLIST cells baked by 000 (F_vesiculo != 0).
release_points_path = data_root / "helcom_fucus_redlist" / "fucus_release_points.geojson"
subbasins_path = data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp"
shoreclass_path = (
    data_root / "shoreclass_bsh_coastline" / "bsh_coastline_k2_flatfraction.parquet"
)
coastline_path = data_root / "bsh_hbmnoku_static" / "coastline.geojson"
h0_fine_path = data_root / "bsh_hbmnoku_static" / "static_file_fine" / "H0_file_fine.nc"
if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
# `release_years` is "first-last" (inclusive) or a single year.
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"
print(f"release years in the analysis: {analysis_years[0]}–{analysis_years[-1]}")
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days + [origin_horizon_days]:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d"
    )
    # The horizon filter is cumulative-to-T (`beach_age_bin < h //
    # age_bin_days`) — strandings in bins fully within T. A horizon below one
    # age bin therefore selects nothing at all, which is a parameter mistake,
    # not an empty map.
    assert h >= age_bin_days, (
        f"horizon {h} d is shorter than one age bin ({age_bin_days} d), so no "
        f"stranding bin falls within it — raise the horizon or lower "
        f"age_bin_days in 024e and here"
    )

# Print-ready output: every saved figure is one full text-width (180 mm)
# figure at 300 dpi, so panels land in the manuscript at their final size
# (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Hex seam stroke. edgecolor="face" means this is not a visible outline -- it
# closes the ~1 px anti-aliasing seam between adjacent polygons so the grid
# reads as a continuous field. The seam is a fixed PIXEL artifact, so pinning
# the stroke to ~1 px at the output dpi keeps the seam closed while letting
# the resulting hex dilation shrink as resolution rises.
hex_seam_lw = 1.1 * 72 / FIGURE_DPI
# Vertical inches a row of maps adds to its own height: the titles above, the
# horizontal colorbar with its ticks and label below.
STRANDED_DECOR_IN = 0.75
# Subbasin border stroke on every map (docs/visualisations.md).
basin_border_lw = 0.52

figure_dir = output_root / "Figures" / "029"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "029"
export_dir.mkdir(parents=True, exist_ok=True)
```

# Read key + pool the season's partitions across years

Layout: flat files under ``output_root/HexAggregates/`` —
``HexAgg_key_r<radius>m.parquet`` and the per-(year, month) partitions
``HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`` written
by 024e. Every partition whose month is in the season is read and pooled;
the release year is parsed from the filename and kept as a column, because
the interannual spread below groups by it.

Each partition carries its build parameters as `df.attrs`. They are taken
per file **before** the concat — pandas drops attrs that differ between the
frames it concatenates — and every pooled partition must agree on the
window, the binning and the rate-model point.

A whole-year partition (no `_mMM` suffix at all, written when 024e itself
ran with `release_month = 0`) carries every month, so it is only read for
`season = "ALL"`.

```python
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
# The key's sidecar carries the HexProj every 024x stage labels with (used to
# put the release cells and the shoreline sub-segments into hex space) and the
# subbasin id → name map the per-subbasin readout below prints.
meta = json.loads(key_path.with_suffix(".json").read_text())
hp = HexProj(**meta["hex_proj"])
subbasin_name = {int(k): v for k, v in meta["subbasin_id_to_name"].items()}

# Keys every pooled partition of one member must share: the reduction window,
# the binning and the rate-model point. Provenance (builder, git_sha,
# slurm_job_id) and the partition's identity (release_year, release_month)
# legitimately differ and are not compared.
POOLED_KEYS = (
    "max_age_days", "age_bin_days", "disp_bin_km", "w_c", "tau_strong_hours",
    "tau_calm_days", "delta", "trap_flat", "trap_wall", "band_m",
)


def pooled_attrs(attrs, store):
    """The build attrs every pooled partition shares, from one `df.attrs`
    per partition."""
    old = sum("max_age_days" not in a for a in attrs)
    if old:
        raise ValueError(
            f"{old} of {len(attrs)} {store} partition(s) carry no build attrs "
            f"(old schema, written before max_age_days was recorded) — rebuild "
            f"them with the current 024e/024f/024g reducers"
        )
    distinct = {tuple(a[k] for k in POOLED_KEYS) for a in attrs}
    if len(distinct) > 1:
        raise ValueError(
            f"pooled {store} partitions disagree on {POOLED_KEYS}: "
            f"{sorted(distinct)} — rebuild them at one setting"
        )
    return dict(zip(POOLED_KEYS, distinct.pop()))


def season_partitions(tag):
    """(path, release_year) of every partition of member `tag` in the season."""
    part_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
        rf"_{re.escape(tag)}\.parquet$"
    )
    files = []
    for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")):
        m = part_re.search(f.name)
        if m is None or int(m.group(1)) not in analysis_years:
            continue
        month = m.group(2)
        if month is None:
            # whole-year partition: carries every month, so only for season ALL
            if season != "ALL":
                continue
        elif int(month) not in season_months:
            continue
        files.append((f, int(m.group(1))))
    if not files:
        raise FileNotFoundError(
            f"no beaching partitions for regime {regime!r}, season {season}, "
            f"release years {release_years}, member {tag!r} at {store_root} — run 024e"
            + (f', or set contrast_member = "" to skip the comparison.'
               if tag == contrast_member else ".")
        )
    return files


beaching_files = season_partitions(member)
parts = []
store_attrs = []
for f, year in beaching_files:
    df = pd.read_parquet(f).reset_index(drop=True)
    store_attrs.append(df.attrs)
    parts.append(df.assign(release_year=year))
beaching = pd.concat(parts, ignore_index=True)
del parts

window = pooled_attrs(store_attrs, f"024e {member}")
assert window["age_bin_days"] == age_bin_days, (
    f"store age_bin_days {window['age_bin_days']} != {age_bin_days}"
)
assert window["disp_bin_km"] == disp_bin_km, (
    f"store disp_bin_km {window['disp_bin_km']} != {disp_bin_km}"
)
max_age_days = window["max_age_days"]
assert max(time_horizons_days) <= max_age_days, (
    f"horizon {max(time_horizons_days)} d is past the store's reduction window "
    f"({max_age_days} d) — rebuild 024e or lower the horizon"
)

pooled_years = sorted(int(y) for y in beaching["release_year"].unique())
# The denominator of every percentage below: the released water-seeded
# drifters of the pooled partitions.
total_released = float(beaching["weight"].sum())
print(f"key: {len(key):,} hexes; pooled {len(beaching_files)} (year, month) partition(s)")
print(f"season {season} (months {season_months}); release years {pooled_years}")
print(f"024e partitions agree on {window}")
print(f"beaching rows: {len(beaching):,}; released drifters (Σweight) "
      f"{total_released:,.0f}")
```

# Rendering helpers

`hex_gdf` sums the stranded weight over a hex column, converts it to a
percentage of the released drifters and joins to the key geometry
(dropping the `-1` sentinel); `hex_map` draws hex polygons, the HELCOM
subbasin borders and the model coastline on a plain EPSG:4326 axis (as
025/026). A map takes its colorbar as an inset axes under it at
axes-fraction width 1.0, so the bar is exactly the map's width — a
fixed-aspect map never fills its gridspec cell, and a cell-sized colorbar
overshoots it. A panel row on one norm passes `cbar=False` and gets one
shared colorbar under the row instead.

`fit_height` sizes a figure from what is drawn rather than from hand
arithmetic: at the fixed page width the maps are width-limited, so any
height beyond their content is a blank band. One draw at a generous seed
height measures the content (titles to colorbar label), and the figure is
cut to it.

```python
def hex_gdf(df, hex_col):
    """Stranded weight per hex (`weight`) and as a percentage of the released
    drifters (`value`, the plotted column), joined to the key geometry."""
    grp = df[df[hex_col] >= 0].groupby(hex_col)["weight"].sum()
    return (
        grp.reset_index()
        .assign(value=lambda d: 100 * d["weight"] / total_released)
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


def extend_of(values, norm):
    """Colorbar `extend` for `values` drawn on `norm`: an arrow wherever values
    fall outside it (the log norm's four-decade floor clips the sparse tail)."""
    below = norm is not None and bool((values < norm.vmin).any())
    above = norm is not None and bool((values > norm.vmax).any())
    return {(False, False): "neither", (True, False): "min",
            (False, True): "max", (True, True): "both"}[below, above]


def hex_map(gdf, ax, norm=None, title=None, label=None, cbar=True):
    """Hex choropleth + subbasin borders + model coastline on a lon/lat axis.
    `norm=None` → linear default scale; a `LogNorm` → shared log scale across
    panels. The colorbar sits under the map at its full width (`cbar=True`)
    or is left to the caller (`False`)."""
    if not gdf.empty:
        legend_kwds = None
        cax = None
        if cbar:
            # Override: a bar ~3.5 mm thick and ~2.5 mm under a half-width map
            # (axes fractions of the map height); the default inset is sized
            # for no map in particular.
            cax = ax.inset_axes([0.0, -0.06, 1.0, 0.035])
            legend_kwds = {
                "orientation": "horizontal",
                "extend": extend_of(gdf["value"], norm),
            } | ({"label": label} if label else {})
            # Bottom-anchored, so constrained_layout's room for the bar under
            # the map holds; centred, the map drifts down and clips the label.
            ax.set_anchor("S")
        gdf.plot(
            ax=ax, column="value", norm=norm, legend=cbar, cax=cax,
            legend_kwds=legend_kwds,
            edgecolor="face", linewidth=hex_seam_lw, zorder=1,
        )
    # Override: neutral mid-grey thin subbasin borders (as 036), under the
    # coastline, so a border reads neither as coastline (black) nor as part of
    # the hex field.
    basin_borders.plot(ax=ax, color="0.5", linewidth=basin_border_lw, zorder=1.5)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def fit_height(fig):
    """Cut the figure height to its drawn content at the fixed page width,
    plus constrained_layout's own pad above and below. Call after drawing,
    before saving; the seed height only has to be generous."""
    fig.draw_without_rendering()
    content_in = fig.get_tightbbox().height
    h_pad = fig.get_layout_engine().get()["h_pad"]
    fig.set_size_inches(FIGURE_WIDTH_IN, content_in + 2 * h_pad)
```

# Model coastline

Every map draws the coastline the particles were advected against: the
outline of the BSH HBMnoku wet cells (`bsh_hbmnoku_static/coastline.geojson`,
the polygons 024a takes the hex water area from), not a generic shoreline
that disagrees with the model grid by up to a coarse cell. The fine nest
takes precedence over the coarse grid inside its bounding box (024a's
merge), and the edges along that box are dropped: there the fine grid's
water simply stops in open sea (the nest seam), which is no coast.

```python
def model_coastline(clip_box):
    """BSH wet-cell outline as lines, fine nest first, nest seam removed,
    clipped to `clip_box`."""
    wet = gpd.read_file(coastline_path)
    h0_fine = xr.open_dataset(h0_fine_path)
    pad_lon = abs(float(h0_fine.lon.diff("lon").mean())) / 2
    pad_lat = abs(float(h0_fine.lat.diff("lat").mean())) / 2
    fine_bbox = box(
        float(h0_fine.lon.min()) - pad_lon, float(h0_fine.lat.min()) - pad_lat,
        float(h0_fine.lon.max()) + pad_lon, float(h0_fine.lat.max()) + pad_lat,
    )
    fine = unary_union(wet.loc[wet["grid"] == "fine", "geometry"].tolist())
    coarse = unary_union(wet.loc[wet["grid"] == "coarse", "geometry"].tolist())
    lines = unary_union([fine.boundary, coarse.difference(fine_bbox).boundary])
    # One fine cell either side of the bbox edge: the seam and nothing else.
    lines = lines.difference(fine_bbox.exterior.buffer(2 * pad_lat))
    return gpd.GeoSeries([lines], crs="EPSG:4326").clip(clip_box)
```

```python
# The hex key tiles the full BSH domain (North Sea included); the extent
# follows it, and beaching occupies wherever it actually reaches.
if any((extent_lon_min, extent_lon_max, extent_lat_min, extent_lat_max)):
    lon_min, lat_min, lon_max, lat_max = (
        extent_lon_min, extent_lat_min, extent_lon_max, extent_lat_max
    )
else:
    lon_min, lat_min, lon_max, lat_max = key.total_bounds
extent = [lon_min, lon_max, lat_min, lat_max]
coast = model_coastline(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)

beached = beaching[beaching["beach_hex"] >= 0]
```

# Release sites

The source side of the substrate contrast below:
`helcom_fucus_redlist/fucus_release_points.geojson`, the REDLIST cells with
`F_vesiculo != 0` baked by `000`. This is the same file `010` opens — it
distributes `particles_per_cell` particles uniformly inside each cell
polygon — so these are the seeding set itself, not a proxy for it. Each
~10 km cell is put into hex space at its representative point.

The key's `fucus_area_m2` is the same source seen through the hex grid (the
per-hex intersection area with those cells); it carries the source weighting
for the substrate contrast, beside the per-cell counts.

```python
release_cells = gpd.read_file(release_points_path)
release_pt = release_cells.representative_point()
release_hex_of_cell = hp.label(
    release_pt.x.to_numpy(), release_pt.y.to_numpy()
).astype("int32")
print(f"release cells: {len(release_cells):,} from {release_points_path}")
print(f"  in {pd.unique(release_hex_of_cell).size:,} hexes; key hexes with "
      f"fucus_area_m2 > 0: {int((key['fucus_area_m2'] > 0).sum()):,}")
```

# Subbasin borders

The HELCOM level-2 subbasins every per-subbasin readout below (and 026a,
028) groups by, drawn on every map as the edges two polygons share — the
same borders as `036`'s seeding map. The polygons' full outlines also trace
every island and skerry of the coast they were digitised against, which
would smear over exactly the coastal hexes the maps are about; the shared
edges are the basin-to-basin borders alone. The Kattegat's open edge to the
Skagerrak, shared with no polygon, is therefore not drawn. The subbasin
codes and their key are on `033`'s shoreline-type map, the paper's
reference map.

```python
subbasin_outline = gpd.read_file(subbasins_path).to_crs(4326).boundary
basin_borders = gpd.GeoSeries(
    [
        a.intersection(b)
        for a, b in itertools.combinations(subbasin_outline, 2)
        if a.intersects(b)
    ],
    crs="EPSG:4326",
)
print(f"subbasin borders: {len(basin_borders)} shared edges between "
      f"{len(subbasin_outline)} HELCOM polygons")
```

# Where stranded, and by origin

**Left:** log-scale stranded share (% of released drifters) per stranding
hex, summed over the store's remaining axes (`disp_bin`, and `shore_type`
where the reducer emits it). Neither is shown on the maps: `shore_type` is a
threshold label on the HELCOM BRISK + CLMS flat fraction `ff` at the
stranding hour (`flat` if `ff >= 50%` else `wall`), a coarser cut than the
rate itself, which scales continuously with `ff` through `trap_wall +
(trap_flat - trap_wall)*ff/100` (production: `trap_flat = 1`, `trap_wall =
0`, so the substrate signal is already inside the plotted weight);
displacement gets its own figure below. The colorbar's `extend` arrow marks
the hexes below the four-decade floor.

**Right:** per release hex, the percentage of its own drifters stranded at
ages below `origin_horizon_days` — Σ `weight` over `beach_hex ≥ 0` and
`beach_age_bin < origin_horizon_days // age_bin_days`, over Σ `weight`, per
`release_hex`, pooled like the left panel — on a linear 0–100 scale. It
shows which source regions lose most of their propagules to the coast. No
release-site overlay on either panel: the right panel's coloured hexes
*are* the release set, and the release sites are on `033`'s reference map.

Two half-width maps side by side, each with its own colorbar under it at
its width, the figure cut to its content (`fit_height`).

```python
gdf_stranded = hex_gdf(beached, "beach_hex")
strand_norm = log_norm(gdf_stranded["value"])
below_floor = gdf_stranded["value"] < strand_norm.vmin
print(f"stranding hexes: {len(gdf_stranded):,}; peak {strand_norm.vmax:.3g} %, "
      f"floor {strand_norm.vmin:.3g} %; {int(below_floor.sum()):,} below the floor "
      f"holding {100 * gdf_stranded['weight'][below_floor].sum() / gdf_stranded['weight'].sum():.2f} % "
      f"of the stranded weight")

released = beaching.groupby("release_hex")["weight"].sum().rename("released")
stranded = (
    beached[beached["beach_age_bin"] < origin_horizon_days // age_bin_days]
    .groupby("release_hex")["weight"].sum().rename("stranded")
)
frac = (
    pd.concat([released, stranded], axis=1).fillna({"stranded": 0})
    .assign(value=lambda d: 100 * d["stranded"] / d["released"])
    .reset_index()
)
frac = frac[frac["release_hex"] >= 0].merge(
    key[["hex_id", "geometry"]], left_on="release_hex", right_on="hex_id"
).pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
print(f"release hexes: {len(frac):,}; stranded by {origin_horizon_days} d per "
      f"release hex: median {frac['value'].median():.1f} %, quartiles "
      f"{frac['value'].quantile(0.25):.1f}–{frac['value'].quantile(0.75):.1f} %; "
      f"pooled {100 * frac['stranded'].sum() / frac['released'].sum():.1f} %")

fig, axes = plt.subplots(1, 2, layout="constrained")
# Seed height only, generous on purpose: fit_height cuts it to the content.
fig.set_size_inches(FIGURE_WIDTH_IN, FIGURE_WIDTH_IN / 2 / domain_aspect + 2 * STRANDED_DECOR_IN)
hex_map(gdf_stranded, axes[0], norm=strand_norm, title="where stranded",
        label="stranded (% of released)")
# Pinned linear 0–100: a share of the hex's own releases, full range on any run.
hex_map(frac, axes[1], norm=Normalize(vmin=0.0, vmax=100.0), title="by release hex",
        label=f"stranded by {origin_horizon_days} d (% of the hex's released)")
fit_height(fig)
fig_path = figure_dir / f"WhereStranded_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Source and sink substrate

The number the source↔sink pair exists for. The shoreline-class table from the data
twin gives every BSH coastline sub-segment a **flat fraction** `ff` — the
share of classified shore that is low and unconsolidated rather than a hard
wall ([docs/beaching.md](../docs/beaching.md)). Aggregated per hex by
`seg_len_m`-weighted mean, exactly as `033` does it: the nest seam dropped
first (the outline of the BSH wet cells closes along the fine nest's edge as
well as along real land, and a shore type in open water is a map artefact),
and an unattributed sub-segment counted 0.5, the convention `024d` applies.

Then three comparisons, all against the domain's own length-weighted mean:

- **Source substrate** — over the hexes carrying release sites, unweighted
  and weighted by the key's `fucus_area_m2`, plus the share of release cells
  whose hex is wall-like (`ff < 0.5`, the threshold the store's `shore_type`
  label uses).
- **Sink substrate** — over the stranding hexes, weighted by stranded
  weight, for `member` and for `contrast_member`.
- **Per subbasin** — release cells against stranded weight, so the transport
  direction is a table and not an impression of the map.

Reading the two members together is what separates the two possible stories:
under `member` the trap factor is `ff/100`, so a flat-coast bias in its
deposit is partly the rate model by construction; `contrast_member`'s trap
is 1 everywhere, so **its** deposit substrate is purely where the drift goes.

```python
FLAT_THRESHOLD = 0.5

sc = pd.read_parquet(shoreclass_path)
sc = sc[~sc["on_seam"]]
sc = sc[np.isfinite(sc["seg_len_m"]) & (sc["seg_len_m"] > 0)]
lon_sc, lat_sc = Transformer.from_crs(3035, 4326, always_xy=True).transform(
    sc["x_3035"].to_numpy("float64"), sc["y_3035"].to_numpy("float64")
)
sc = sc.assign(
    hex_id=hp.label(lon_sc, lat_sc).astype("int32"),
    ff=sc["flat_fraction"].fillna(0.5),
)
grp = sc.assign(_ff_len=sc["ff"] * sc["seg_len_m"]).groupby("hex_id")
shore = pd.DataFrame(
    {
        "coast_len_m": grp["seg_len_m"].sum(),
        "flat_fraction": grp["_ff_len"].sum() / grp["seg_len_m"].sum(),
    }
)
domain_ff = float(np.average(shore["flat_fraction"], weights=shore["coast_len_m"]))
print(f"hexes carrying coastline: {len(shore):,}; "
      f"domain length-weighted mean flat fraction: {domain_ff:.3f}")
```

The contrast member's deposit is the only extra read this block needs: one
pass over its partitions reduced straight to `beach_hex → weight`, the
production member's coming from the rows already in memory.

```python
strand_by_hex = {member: beached.groupby("beach_hex")["weight"].sum()}
if contrast_member:
    parts = []
    contrast_attrs = []
    for f, _ in season_partitions(contrast_member):
        d = pd.read_parquet(f, columns=["beach_hex", "weight"])
        contrast_attrs.append(d.attrs)
        parts.append(d[d["beach_hex"] >= 0].groupby("beach_hex")["weight"].sum())
    # The two members differ in their rate model by design; they must share
    # the window, or the comparison mixes stranding horizons.
    contrast_window = pooled_attrs(contrast_attrs, f"024e {contrast_member}")
    assert contrast_window["max_age_days"] == max_age_days, (
        f"{contrast_member} window {contrast_window['max_age_days']} d != "
        f"{member} window {max_age_days} d"
    )
    strand_by_hex[contrast_member] = pd.concat(parts).groupby(level=0).sum()
for tag, w in strand_by_hex.items():
    print(f"{tag}: {len(w):,} stranding hexes, Σweight {w.sum():,.0f}")
```

```python
source = (
    key.loc[key["fucus_area_m2"] > 0, ["hex_id", "fucus_area_m2"]]
    .join(shore, on="hex_id")
)
source_coastal = source[source["flat_fraction"].notna()]
cell_ff = pd.Series(release_hex_of_cell, name="hex_id").to_frame().join(
    shore["flat_fraction"], on="hex_id"
)["flat_fraction"]

print(f"domain reference: length-weighted mean flat fraction {domain_ff:.3f}")
print(f"\nSOURCE — {len(source):,} hexes with Fucus area, of which "
      f"{len(source_coastal):,} carry coastline "
      f"({100 * source_coastal['fucus_area_m2'].sum() / source['fucus_area_m2'].sum():.1f} %"
      f" of the Fucus area)")
print(f"  mean flat fraction, unweighted:        "
      f"{source_coastal['flat_fraction'].mean():.3f}")
print(f"  mean flat fraction, fucus_area_m2-wtd: "
      f"{np.average(source_coastal['flat_fraction'], weights=source_coastal['fucus_area_m2']):.3f}")
print(f"  mean flat fraction, coast-length-wtd:  "
      f"{np.average(source_coastal['flat_fraction'], weights=source_coastal['coast_len_m']):.3f}"
      f"  (the form the domain reference takes)")
print(f"  release cells on a coastline-carrying hex: {int(cell_ff.notna().sum()):,} "
      f"of {len(cell_ff):,}")
print(f"  of those, on wall-like coast (ff < {FLAT_THRESHOLD}): "
      f"{100 * (cell_ff < FLAT_THRESHOLD).sum() / cell_ff.notna().sum():.1f} %")

for tag, w in strand_by_hex.items():
    d = w.rename("weight").to_frame().join(shore)
    on_coast = d[d["flat_fraction"].notna()]
    print(f"\nSINK — member {tag}")
    print(f"  stranded weight on coastline-carrying hexes: "
          f"{100 * on_coast['weight'].sum() / d['weight'].sum():.1f} % of "
          f"{d['weight'].sum():,.0f}")
    print(f"  mean flat fraction, stranded-weight-wtd: "
          f"{np.average(on_coast['flat_fraction'], weights=on_coast['weight']):.3f}")
    print(f"  weight landing on flat coast (ff >= {FLAT_THRESHOLD}): "
          f"{100 * on_coast.loc[on_coast['flat_fraction'] >= FLAT_THRESHOLD, 'weight'].sum() / on_coast['weight'].sum():.1f} %")
```

# Source and sink per subbasin

Release cells against stranded weight, both labelled by the HELCOM subbasin
of their hex. The HELCOM polygons leave a coastal hex centred on land or
inside an archipelago unassigned, which is a large share of both columns
here, so the percentages are shares **of the assigned rows** and the
unassigned share is printed beside them rather than folded in.

```python
subbasin_of_hex = key.set_index("hex_id")["helcom_subbasin"]
cell_subbasin = pd.Series(release_hex_of_cell).map(subbasin_of_hex)
per_sub = pd.DataFrame(
    {"cells": cell_subbasin[cell_subbasin >= 0].map(subbasin_name).value_counts()}
)
for tag, w in strand_by_hex.items():
    sub = w.index.map(subbasin_of_hex)
    # float64: the store's weight is float32, and a float32 share prints as
    # 49.290001 however it is rounded.
    per_sub[tag] = w[sub >= 0].astype("float64").groupby(
        sub[sub >= 0].map(subbasin_name)
    ).sum()
per_sub = per_sub.fillna(0.0).sort_values("cells", ascending=False)
table = pd.DataFrame({"source cells": per_sub["cells"].astype(int)})
table["source %"] = 100 * per_sub["cells"] / per_sub["cells"].sum()
for tag in strand_by_hex:
    table[f"sink % ({tag})"] = 100 * per_sub[tag] / per_sub[tag].sum()
print(table.round(2).to_string())

unassigned_src = float((cell_subbasin < 0).sum() / len(release_hex_of_cell))
print(f"\nrelease cells HELCOM leaves unassigned: {100 * unassigned_src:.1f} %")
for tag, w in strand_by_hex.items():
    u = w[w.index.map(subbasin_of_hex) < 0].sum() / w.sum()
    print(f"stranded weight HELCOM leaves unassigned ({tag}): {100 * u:.1f} %")
```

# Beaching age horizons

Cumulative where-stranded density for strandings at ages < T, as a
percentage of all released drifters. `beach_age_bin` bins age at
`age_bin_days`, so horizon T keeps `beach_age_bin < T // age_bin_days`
(strandings in bins fully within T). Every panel is on the where-stranded
map's `LogNorm` — the horizons are cumulative, so the full-window deposit is
the largest any panel reaches — with one shared colorbar under the row, so
the fill-in over time reads on the same scale as the map above.

```python
horizon_bins = {h: h // age_bin_days for h in time_horizons_days}
gdfs = {
    h: hex_gdf(beached[beached["beach_age_bin"] < horizon_bins[h]], "beach_hex")
    for h in time_horizons_days
}
all_h = pd.concat([g["value"] for g in gdfs.values() if not g.empty])

ncols = len(time_horizons_days)
fig, axes = plt.subplots(1, ncols, layout="constrained", squeeze=False)
# Each map takes a full column (no per-panel colorbar), plus the titles above
# and the shared colorbar band below.
fig.set_size_inches(
    FIGURE_WIDTH_IN, FIGURE_WIDTH_IN / ncols / domain_aspect + STRANDED_DECOR_IN
)
for ax, h in zip(axes.flat, time_horizons_days):
    hex_map(gdfs[h], ax, norm=strand_norm, title=f"stranded by {h} d", cbar=False)
fig.colorbar(
    ScalarMappable(norm=strand_norm), ax=axes, location="bottom",
    extend=extend_of(all_h, strand_norm), label="stranded (% of released)",
)
fig_path = figure_dir / f"BeachingHorizons_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# GeoJSON export of the per-horizon deposit

One file per age horizon: the cumulative stranded weight per hex
(`stranded_weight`, expected particles) and the same as a percentage of the
released drifters (`pct_released`, the mapped value), with the 024a
geometry, EPSG:4326 — the horizon panels above in a form collaborators can
open in R or QGIS.

```python
for h in time_horizons_days:
    g = (
        gdfs[h]
        .rename(columns={"weight": "stranded_weight", "value": "pct_released"})
        .drop(columns="hex_id")
    )
    export_path = export_dir / f"029_{regime}_{season}_{member}_T{h}d.geojson"
    g.to_file(export_path, driver="GeoJSON")
    print(f"wrote {export_path} ({len(g):,} hexes)")
```


# Afloat by origin subbasin

Per HELCOM origin subbasin (the key's `helcom_subbasin` of `release_hex`),
the percentage of its released drifters **not yet stranded** at each
horizon `T`:

`afloat(T) = 100 · (Σweight − Σweight[beach_hex ≥ 0 & beach_age_bin < T / age_bin_days]) / Σweight`

with Σweight the origin's released drifters (deposits + residual). Pooled
over every release year from the summed weights, plus the min–max of the
same reduction per release year. "Afloat" is everything not stranded: it
includes drifters that left the BSH domain, and no lifetime or sinking is
applied — a biological loss would have to be multiplied in downstream.
Origin `-1` (a release hex outside every named polygon) is kept as its own
row and printed apart from the named subbasins; the `all` row is the
basin-wide value.

```python
origin_of_hex = key.set_index("hex_id")["helcom_subbasin"]
origin = beaching["release_hex"].map(origin_of_hex).fillna(-1).astype(int)

# Released and stranded-by-T weight per (origin, release_year, horizon).
rows = []
for h in time_horizons_days:
    stranded = (beaching["beach_hex"] >= 0) & (beaching["beach_age_bin"] < h // age_bin_days)
    rows.append(
        beaching[["release_year", "weight"]]
        .assign(origin=origin, stranded=beaching["weight"].where(stranded, 0.0))
        .groupby(["origin", "release_year"])[["weight", "stranded"]].sum()
        .assign(horizon_days=h)
    )
afloat_year = pd.concat(rows).reset_index()
afloat_year["afloat"] = afloat_year["weight"] - afloat_year["stranded"]
afloat_year = pd.concat(
    [
        afloat_year,
        afloat_year.groupby(["release_year", "horizon_days"], as_index=False)[
            ["weight", "stranded", "afloat"]
        ].sum().assign(origin="all"),
    ],
    ignore_index=True,
)
afloat_year["afloat_pct"] = 100 * afloat_year["afloat"] / afloat_year["weight"]
pooled = afloat_year.groupby(["origin", "horizon_days"], sort=False)[
    ["weight", "afloat"]
].sum()
spread = afloat_year.groupby(["origin", "horizon_days"], sort=False)["afloat_pct"]
afloat = pd.DataFrame(
    {
        "released_weight": pooled["weight"],
        "afloat_pct": 100 * pooled["afloat"] / pooled["weight"],
        "afloat_pct_min": spread.min(),
        "afloat_pct_max": spread.max(),
        "n_years": spread.size(),
    }
).reset_index()
names = {**subbasin_name, -1: "outside named subbasins", "all": "all origins"}
afloat.insert(1, "origin_name", afloat["origin"].map(names))
afloat = afloat.rename(columns={"origin": "origin_subbasin"})

afloat_path = export_dir / f"afloat_by_origin_{regime}_{season}_{member}.csv"
with afloat_path.open("w") as fh:
    fh.write(
        f"# 029_BeachingMaps: percentage of each origin subbasin's released "
        f"water-seeded drifters not stranded by age horizon_days "
        f"(regime {regime}, season {season}, member {member}, release years "
        f"{pooled_years[0]}-{pooled_years[-1]}, reduction window "
        f"{max_age_days} d)\n"
        "# afloat = released - stranded at ages < horizon_days; it includes "
        "drifters that left the BSH domain, and no lifetime or sinking is "
        "applied\n"
        "# afloat_pct pools the release years (summed weights); "
        "afloat_pct_min/max are the per-release-year range; "
        "origin_subbasin -1 = release hex outside every named HELCOM "
        "subbasin, all = every origin\n"
    )
    afloat.to_csv(fh, index=False)
print(f"wrote {afloat_path} ({len(afloat)} rows)")


def afloat_cells(frame):
    """origin × horizon of 'pooled [min–max]' strings."""
    cells = frame.assign(
        cell=frame.apply(
            lambda r: f"{r.afloat_pct:5.1f} [{r.afloat_pct_min:5.1f}–{r.afloat_pct_max:5.1f}]",
            axis=1,
        )
    )
    return cells.pivot(index="origin_name", columns="horizon_days", values="cell").rename(
        columns=lambda h: f"{h} d"
    )


named = afloat[afloat["origin_subbasin"].map(lambda o: o not in (-1, "all"))]
order = named[named["horizon_days"] == time_horizons_days[0]].sort_values(
    "released_weight", ascending=False
)["origin_name"]
print(f"afloat % of released by origin subbasin, pooled [per-year min–max], "
      f"{season}, {member}:")
print(afloat_cells(named).reindex(order).to_string())
print()
print(afloat_cells(afloat[afloat["origin_subbasin"].isin([-1, "all"])]).to_string())
share_outside = float(
    afloat.loc[
        (afloat["origin_subbasin"] == -1) & (afloat["horizon_days"] == time_horizons_days[0]),
        "released_weight",
    ].sum()
    / total_released
)
print(f"origin -1 holds {100 * share_outside:.1f} % of the released drifters")
```


# Cumulative beached percentage vs. age

Percentage of all released drifters stranded by each elapsed-time horizon —
the stranding time course, summed over the store's other axes. The pooled
curve is drawn with the **interannual min–max** as a shaded band: the
per-year curves are the same reduction restricted to one `release_year`, so
the band is the honest spread behind the single pooled number. Every curve
runs over the full age axis of the window, a bin without strandings
carrying the cumulative value of the bin before, so the per-year band
compares bin for bin.

The curve is exported as
`Exports/029/stranded_curve_<regime>_<season>_<member>.csv` (`age_days` =
the age-bin right edge T, stranding at ages < T; pooled %; per-year min and
max), which `035` draws above its area panels on one age axis.

```python
age_bins = pd.RangeIndex(max_age_days // age_bin_days)


def beached_fraction_curve(df):
    """Cumulative stranded weight as % of released, at each age-bin right edge."""
    b = df[df["beach_hex"] >= 0]
    total = float(df["weight"].sum())
    per_bin = b.groupby("beach_age_bin")["weight"].sum().reindex(age_bins, fill_value=0.0)
    return pd.Series(
        100 * per_bin.cumsum().to_numpy() / max(total, 1.0),
        index=pd.Index((age_bins.to_numpy() + 1) * age_bin_days, name="age (days)"),
        name="stranded (% of released)",
    )


curve = beached_fraction_curve(beaching)
curve_per_year = pd.DataFrame(
    {y: beached_fraction_curve(g) for y, g in beaching.groupby("release_year")}
)

curve_table = pd.DataFrame(
    {
        "stranded_pct": curve,
        "stranded_pct_year_min": curve_per_year.min(axis=1),
        "stranded_pct_year_max": curve_per_year.max(axis=1),
        "n_years": curve_per_year.shape[1],
    }
).rename_axis("age_days")
curve_path = export_dir / f"stranded_curve_{regime}_{season}_{member}.csv"
with curve_path.open("w") as fh:
    fh.write(
        f"# 029_BeachingMaps: cumulative percentage of the released water-seeded "
        f"drifters stranded at ages < age_days (regime {regime}, season {season}, "
        f"member {member}, release years {pooled_years[0]}-{pooled_years[-1]}, "
        f"age bin {age_bin_days} d, reduction window {max_age_days} d)\n"
        "# stranded_pct pools the release years (summed weights); "
        "stranded_pct_year_min/max are the per-release-year range\n"
    )
    curve_table.to_csv(fh)
print(f"wrote {curve_path} ({len(curve_table)} ages)")

fig, ax = plt.subplots(layout="constrained")
# Print-ready width; a line panel needs no map aspect, so half the page width
# in height reads as a normal wide chart.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
curve.plot(ax=ax)
ax.fill_between(
    curve.index, curve_table["stranded_pct_year_min"], curve_table["stranded_pct_year_max"],
    alpha=0.3, label=f"interannual range ({pooled_years[0]}–{pooled_years[-1]})",
)
ax.set_ylabel("stranded (% of released)")
ax.legend()
fig_path = figure_dir / f"BeachedFractionCurve_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Displacement at stranding

Stranded share (% of released drifters) over `disp_bin` — the crow-flies
displacement from the release point at the deposit step, in `disp_bin_km`
bins — pooled over every source hex. The residual (never-beached) rows carry `disp_bin = -1` and drop
out with the `beach_hex = -1` filter. Short-displacement strandings are kept: a
strong-wave event stranding freshly released material near home is a real
outcome of the threshold rate model, not an artefact to mask.

```python
disp_weight = 100 * beached.groupby("disp_bin")["weight"].sum() / total_released
disp_weight.index = (disp_weight.index.to_numpy() + 0.5) * disp_bin_km
disp_weight.index.name = "displacement (km)"
median_disp = float(
    np.interp(0.5, disp_weight.cumsum() / disp_weight.sum(), disp_weight.index.to_numpy())
)

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
# drawstyle="steps-mid": the x axis is a binned quantity, so a step reads as
# the histogram it is; a categorical bar plot would label every bin and the
# axis runs to hundreds of bins (docs/visualisations.md).
disp_weight.plot(ax=ax, drawstyle="steps-mid")
ax.set_ylabel("stranded (% of released)")
fig_path = figure_dir / f"StrandingTravelDistance_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
print(f"weight-weighted median displacement at stranding: {median_disp:.1f} km")
plt.show()
```

# Validation / summary

```python
n_beached = float(beached["weight"].sum())
print(f"regime={regime}, hex_radius={hex_radius} m, season={season} "
      f"(months {season_months}), age_bin_days={age_bin_days}, "
      f"window {max_age_days} d")
print(f"  drifters (Σweight): {total_released:,.0f}")
print(f"  beached:           {n_beached:,.0f} "
      f"({100 * n_beached / max(total_released, 1):.1f}%)")
print(f"  stranding hexes:   {beached['beach_hex'].nunique():,}")
print(f"  source hexes:      {frac['release_hex'].nunique():,}")
print(f"  stranded by {origin_horizon_days} d, per release hex: median "
      f"{frac['value'].median():.1f} %")
print(f"  median displacement at stranding: {median_disp:,.1f} km")
print(f"  member:            {member}")
print("  cumulative stranded % of released by age (pooled, interannual min–max):")
for t, row in curve_table.iterrows():
    print(f"    {t:>5} d: {row['stranded_pct']:6.2f} %  "
          f"[{row['stranded_pct_year_min']:6.2f} .. {row['stranded_pct_year_max']:6.2f}]")
```
