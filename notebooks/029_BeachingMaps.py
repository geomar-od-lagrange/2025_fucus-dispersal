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
# `024e_BuildBeaching` (+ the `024a` key) — no trajectory zarrs, no Dask.
# One regime and one release `season` per run; the season's monthly
# partitions are pooled across every available release year.
#
# The rate model is a two-state hazard — a background rate in the near-shore
# band plus a strong-wave rate switched on where the onshore Stokes forcing
# exceeds a threshold (`docs/beaching.md`). Its settings are baked into the
# reducer's opaque `member` tag, which selects the partitions read here and
# labels every figure; this notebook never parses the tag. Draws:
#
# 1. **Where-stranded density** — beached weight (expected particles) per
#    `beach_hex` (log scale), with the *Fucus* release sites overlaid so the
#    map reads as a source→sink picture rather than a deposition field alone.
# 2. **Beached fraction per source hex** — of the drifters released in each
#    hex, what fraction strands within the viability window (linear 0–1).
# 3. **Beaching age horizons** — cumulative where-stranded density for
#    strandings at age ≤ T (the beaching analogue of `026`'s horizons).
# 4. **Cumulative beached fraction vs. age** — the stranding time course,
#    with the 2016–2019 interannual min–max as a band around the pooled
#    curve. Maps stay pooled; only the scalar-per-horizon curve carries the
#    spread.
# 5. **Travel distance at stranding** — stranded weight over `disp_bin`, the
#    crow-flies displacement from release at the deposit step. Short-travel
#    strandings are a real outcome under a threshold rate model, never a
#    mask.
#
# Between 1. and 2. sits the number the overlay exists for: the **substrate
# contrast between source and sink**, the shoreline flat fraction where the
# release sites are against the flat fraction where the weight lands, under
# this member and under a substrate-inert `contrast_member`.
#
# Figures are print-ready (full page width, 300 dpi) and the per-horizon
# where-stranded field is exported as GeoJSON.

# %%
import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyproj import Transformer
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

from hextraj import HexProj

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
# Read root of the data twin: the Fucus release sites 010 seeds from and the
# shoreline-class table the substrate contrast reads.
data_root = "../data"
# Which regime's beaching partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024e). Must match the files on disk.
hex_radius = 6000
# Age-bin granularity of the beaching store (must match 024e). Horizons must
# be whole multiples of this.
age_bin_days = 10
# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year. Selects
# which of the monthly `_mMM` partitions 024e wrote are read and pooled.
season = "ALL"

# Rate-model member: the opaque tag the reducer put in the store filename
# (e.g. "step_wc0p1_ts3_tcinf_tf1_tw0"). Selects the partitions and tags
# the figures.
member = "step_wc0p1_ts3_tcinf_tf1_tw0"
# Second member, read only by the source/sink substrate block: a member whose
# trap term is degenerate (trap_flat = trap_wall = 1), so substrate does not
# enter its rate at all and its deposit is where the drift alone puts the
# material. "" skips the comparison.
contrast_member = "step_wc0p1_ts3_tcinf"
# Travel-distance bin width (km) of the store's `disp_bin` axis (must match
# the reducer).
disp_bin_km = 10.0

# Overlay the Fucus release sites on the where-stranded map.
show_release_points = True

# Beaching age horizons to map, in days. Each must be a multiple of
# age_bin_days.
time_horizons_days_csv = "10,20,50"

# Map extent (degrees). The 024a key tiles the whole BSH domain including the
# North Sea, which is empty for Fucus and costs ~35% of panel height; cropping
# to the Baltic proper buys ~1.35x px/km for free. Set all four to 0 to fall
# back to the key's full bounds.
extent_lon_min = 9.0
extent_lon_max = 30.7
extent_lat_min = 53.0
extent_lat_max = 66.0

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
data_root = Path(data_root)
# The file 010 seeds from: the REDLIST cells baked by 000 (F_vesiculo != 0).
release_points_path = data_root / "helcom_fucus_redlist" / "fucus_release_points.geojson"
shoreclass_path = (
    data_root / "shoreclass_bsh_coastline" / "bsh_coastline_k2_flatfraction.parquet"
)
if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
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
# Release-site overlay style. AGENTS.md forbids colour/marker overrides; this
# one is defended in docs/visualisations.md. Short version: the field under it
# is a four-decade LogNorm of the default (viridis) map, so the marker has to
# read against both a near-black low end and a near-white/yellow high end —
# magenta is the anti-viridis hue already used for the 3 m isobath in 025/026.
# Hollow, because sources and sinks overlap in the western Baltic and a filled
# marker would hide the very hex value the reader is comparing against. Sized
# to roughly the 10 km REDLIST cell it stands for at 180 mm page width.
release_marker_color = "magenta"
release_marker_size = 2.5
release_marker_lw = 0.4

figure_dir = output_root / "Figures" / "029"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "029"
export_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Read key + pool the season's partitions across years
#
# Layout: flat files under ``output_root/HexAggregates/`` —
# ``HexAgg_key_r<radius>m.parquet`` and the per-(year, month) partitions
# ``HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`` written
# by 024e. Every partition whose month is in the season is read and pooled;
# the release year is parsed from the filename and kept as a column, because
# the interannual spread below groups by it.
#
# A whole-year partition (no `_mMM` suffix at all, written when 024e itself
# ran with `release_month = 0`) carries every month, so it is only read for
# `season = "ALL"`.

# %%
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)
# The key's sidecar carries the HexProj every 024x stage labels with (used to
# put the release cells and the shoreline sub-segments into hex space) and the
# subbasin id → name map the per-subbasin readout below prints.
meta = json.loads(key_path.with_suffix(".json").read_text())
hp = HexProj(**meta["hex_proj"])
subbasin_name = {int(k): v for k, v in meta["subbasin_id_to_name"].items()}

_PART_RE = re.compile(
    rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
    rf"_{re.escape(member)}\.parquet$"
)


beaching_files = []
for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")):
    m = _PART_RE.search(f.name)
    if m is None:
        continue
    month = m.group(2)
    if month is None:
        # whole-year partition: carries every month, so only for season ALL
        if season != "ALL":
            continue
    elif int(month) not in season_months:
        continue
    beaching_files.append((f, int(m.group(1))))
if not beaching_files:
    raise FileNotFoundError(
        f"no beaching partitions for regime {regime!r}, season {season}, "
        f"member {member!r} at {store_root} — run 024e."
    )

beaching = pd.concat(
    [
        pd.read_parquet(f).reset_index(drop=True).assign(release_year=year)
        for f, year in beaching_files
    ],
    ignore_index=True,
)
release_years = sorted(int(y) for y in beaching["release_year"].unique())
print(f"key: {len(key):,} hexes; pooled {len(beaching_files)} (year, month) partition(s)")
print(f"season {season} (months {season_months}); release years {release_years}")
print(f"beaching rows: {len(beaching):,}")


# %% [markdown]
# # Rendering helpers
#
# `hex_gdf` sums a value column over a hex column and joins to the key
# geometry (dropping the `-1` sentinel); `hex_map` draws hex polygons +
# coastline on a plain EPSG:4326 axis (as 025/026), with the colorbar as an
# inset axes at axes-fraction height 1.0 so it is exactly the map height —
# a fixed-aspect map never fills its gridspec cell, and a cell-sized
# colorbar overshoots it.

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


def hex_map(gdf, ax, norm=None, title=None, label=None):
    """Hex choropleth + coastline on a lon/lat axis. `norm=None` → linear
    default scale; a `LogNorm` → shared log scale across panels."""
    if not gdf.empty:
        cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
        gdf.plot(
            ax=ax, column="value", norm=norm, legend=True, cax=cax,
            legend_kwds={"label": label} if label else None,
            edgecolor="face", linewidth=hex_seam_lw, zorder=1,
        )
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
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
# The hex key tiles the full BSH domain (North Sea included); the extent
# follows it, and beaching occupies wherever it actually reaches.
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

beached = beaching[beaching["beach_hex"] >= 0]

# %% [markdown]
# # Release sites
#
# The overlay geometry: `helcom_fucus_redlist/fucus_release_points.geojson`,
# the REDLIST cells with `F_vesiculo != 0` baked by `000`. This is the same
# file `010` opens — it distributes `particles_per_cell` particles uniformly
# inside each cell polygon — so the overlay is the seeding set itself, not a
# proxy for it. The cells are ~10 km squares, so one marker per cell at their
# representative point stands for the cell footprint at page width; the marker
# is drawn, not the polygon, for the reasons in
# [docs/visualisations.md](../docs/visualisations.md).
#
# The key's `fucus_area_m2` is the same source seen through the hex grid (the
# per-hex intersection area with those cells); it carries the source weighting
# for the substrate contrast below, while the overlay stays on the cells.

# %%
release_cells = gpd.read_file(release_points_path)
release_pt = release_cells.representative_point()
release_hex_of_cell = hp.label(
    release_pt.x.to_numpy(), release_pt.y.to_numpy()
).astype("int32")
print(f"release cells: {len(release_cells):,} from {release_points_path}")
print(f"  in {pd.unique(release_hex_of_cell).size:,} hexes; key hexes with "
      f"fucus_area_m2 > 0: {int((key['fucus_area_m2'] > 0).sum()):,}")

# %% [markdown]
# # Where-stranded density
#
# Log-scale beached weight (expected particles) per stranding hex, summed over
# the store's remaining axes (`disp_bin`, and `shore_type` where the reducer
# emits it). Neither is shown on the maps: `shore_type` is a threshold label
# on the HELCOM BRISK + CLMS flat fraction `ff` at the stranding hour (`flat`
# if `ff >= 50%` else `wall`), a coarser cut than the rate itself, which
# scales continuously with `ff` through `trap_wall + (trap_flat -
# trap_wall)*ff/100` (production: `trap_flat = 1`, `trap_wall = 0`, so the
# substrate signal is already inside the plotted weight); travel distance
# gets its own figure below.

# %%
gdf_stranded = hex_gdf(beached, "beach_hex")
strand_norm = log_norm(gdf_stranded["value"])

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, grid_height_in(1, 1, domain_aspect, FIGURE_WIDTH_IN))
hex_map(gdf_stranded, ax, norm=strand_norm,
        title="stranded weight, Fucus release sites overlaid"
        if show_release_points else "stranded weight",
        label="stranded weight (particles)")
if show_release_points:
    # Hollow magenta rings, one per REDLIST release cell; style constants and
    # their rationale are at the top of the notebook.
    ax.plot(
        release_pt.x, release_pt.y, linestyle="none", marker="o",
        markersize=release_marker_size, markerfacecolor="none",
        markeredgecolor=release_marker_color, markeredgewidth=release_marker_lw,
        zorder=3,
    )
fig_path = figure_dir / f"WhereStranded_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Source and sink substrate
#
# The number the overlay exists for. The shoreline-class table from the data
# twin gives every BSH coastline sub-segment a **flat fraction** `ff` — the
# share of classified shore that is low and unconsolidated rather than a hard
# wall ([docs/beaching.md](../docs/beaching.md)). Aggregated per hex by
# `seg_len_m`-weighted mean, exactly as `033` does it: the nest seam dropped
# first (the outline of the BSH wet cells closes along the fine nest's edge as
# well as along real land, and a shore type in open water is a map artefact),
# and an unattributed sub-segment counted 0.5, the convention `024d` applies.
#
# Then three comparisons, all against the domain's own length-weighted mean:
#
# - **Source substrate** — over the hexes carrying release sites, unweighted
#   and weighted by the key's `fucus_area_m2`, plus the share of release cells
#   whose hex is wall-like (`ff < 0.5`, the threshold the store's `shore_type`
#   label uses).
# - **Sink substrate** — over the stranding hexes, weighted by stranded
#   weight, for `member` and for `contrast_member`.
# - **Per subbasin** — release cells against stranded weight, so the transport
#   direction is a table and not an impression of the map.
#
# Reading the two members together is what separates the two possible stories:
# under `member` the trap factor is `ff/100`, so a flat-coast bias in its
# deposit is partly the rate model by construction; `contrast_member`'s trap
# is 1 everywhere, so **its** deposit substrate is purely where the drift goes.

# %%
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

# %% [markdown]
# The contrast member's deposit is the only extra read this block needs: one
# pass over its partitions reduced straight to `beach_hex → weight`, the
# production member's coming from the rows already in memory.

# %%
strand_by_hex = {member: beached.groupby("beach_hex")["weight"].sum()}
if contrast_member:
    contrast_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
        rf"_{re.escape(contrast_member)}\.parquet$"
    )
    parts = []
    for f in sorted(
        store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")
    ):
        m = contrast_re.search(f.name)
        if m is None:
            continue
        month = m.group(2)
        if month is None:
            if season != "ALL":
                continue
        elif int(month) not in season_months:
            continue
        d = pd.read_parquet(f, columns=["beach_hex", "weight"])
        parts.append(d[d["beach_hex"] >= 0].groupby("beach_hex")["weight"].sum())
    if not parts:
        raise FileNotFoundError(
            f"no beaching partitions for contrast_member {contrast_member!r} "
            f"(regime {regime!r}, season {season}) at {store_root} — run 024e, "
            f'or set contrast_member = "" to skip the comparison.'
        )
    strand_by_hex[contrast_member] = pd.concat(parts).groupby(level=0).sum()
for tag, w in strand_by_hex.items():
    print(f"{tag}: {len(w):,} stranding hexes, Σweight {w.sum():,.0f}")

# %%
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

# %% [markdown]
# # Source and sink per subbasin
#
# Release cells against stranded weight, both labelled by the HELCOM subbasin
# of their hex. The HELCOM polygons leave a coastal hex centred on land or
# inside an archipelago unassigned, which is a large share of both columns
# here, so the percentages are shares **of the assigned rows** and the
# unassigned share is printed beside them rather than folded in.

# %%
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

# %% [markdown]
# # Beached fraction per source hex
#
# Of the drifters released in each source hex, the fraction that strands
# within the viability window (linear 0–1). Highlights which source regions
# lose most of their propagules to the coast.
#
# No release-site overlay here: this panel's coloured hexes *are* the release
# set, so the markers would fall on the values they are meant to locate and
# add nothing but ink.

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

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, grid_height_in(1, 1, domain_aspect, FIGURE_WIDTH_IN))
# Linear default scale (fraction in [0, 1]); no norm override needed.
hex_map(frac, ax, norm=None, title="beached fraction per source hex",
        label="beached fraction")
fig_path = figure_dir / f"BeachedFraction_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
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
fig, axes = plt.subplots(1, ncols, layout="constrained", squeeze=False)
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(1, ncols, domain_aspect, FIGURE_WIDTH_IN)
)
for ax, h in zip(axes.flat, time_horizons_days):
    hex_map(gdfs[h], ax, norm=horizon_norm, title=f"stranded by {h} d",
            label="stranded weight (particles)")
fig_path = figure_dir / f"BeachingHorizons_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # GeoJSON export of the per-horizon deposit
#
# One file per age horizon: the cumulative stranded weight per hex with the
# 024a geometry, EPSG:4326 — the horizon panels above in a form
# collaborators can open in R or QGIS.

# %%
for h in time_horizons_days:
    g = gdfs[h].rename(columns={"value": "stranded_weight"}).drop(columns="hex_id")
    export_path = export_dir / f"029_{regime}_{season}_{member}_T{h}d.geojson"
    g.to_file(export_path, driver="GeoJSON")
    print(f"wrote {export_path} ({len(g):,} hexes)")


# %% [markdown]
# # Cumulative beached fraction vs. age
#
# Share of all released drifters stranded by each elapsed-time horizon — the
# stranding time course, summed over the store's other axes. The pooled curve
# is drawn with the 2016–2019 **interannual min–max** as a shaded band: the
# per-year curves are the same reduction restricted to one `release_year`, so
# the band is the honest spread behind the single pooled number.

# %%
def beached_fraction_curve(df):
    """Cumulative stranded / released weight at each age-bin right edge."""
    b = df[df["beach_hex"] >= 0]
    total = float(df["weight"].sum())
    per_bin = b.groupby("beach_age_bin")["weight"].sum().sort_index()
    return pd.Series(
        per_bin.cumsum().to_numpy() / max(total, 1.0),
        index=pd.Index(
            (per_bin.index.to_numpy() + 1) * age_bin_days, name="age (days)"
        ),
        name="cumulative beached fraction",
    )


curve = beached_fraction_curve(beaching)
# Reindexing onto the pooled age axis leaves a gap wherever a year stranded
# nothing in that bin. The curve is cumulative, so the value there is the
# previous bin's, not "missing": forward-fill it, and read the bins before a
# year's first stranding as the 0 they are. Without this the band's min is the
# nanmin over whichever years happen to carry a bin.
per_year = pd.DataFrame(
    {y: beached_fraction_curve(g) for y, g in beaching.groupby("release_year")}
).reindex(curve.index).ffill().fillna(0.0)

fig, ax = plt.subplots(layout="constrained")
# Print-ready width; a line panel needs no map aspect, so half the page width
# in height reads as a normal wide chart.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
curve.plot(ax=ax)
ax.fill_between(
    curve.index, per_year.min(axis=1), per_year.max(axis=1), alpha=0.3,
    label=f"interannual range ({release_years[0]}–{release_years[-1]})",
)
ax.set_ylabel("cumulative beached fraction")
ax.legend()
fig_path = figure_dir / f"BeachedFractionCurve_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Travel distance at stranding
#
# Stranded weight over `disp_bin` — the crow-flies displacement from the
# release point at the deposit step, in `disp_bin_km` bins — pooled over every
# source hex. The residual (never-beached) rows carry `disp_bin = -1` and drop
# out with the `beach_hex = -1` filter. Short-travel strandings are kept: a
# strong-wave event stranding freshly released material near home is a real
# outcome of the threshold rate model, not an artefact to mask.

# %%
disp_weight = beached.groupby("disp_bin")["weight"].sum()
disp_weight.index = (disp_weight.index.to_numpy() + 0.5) * disp_bin_km
disp_weight.index.name = "travel distance (km)"
median_disp = float(
    np.interp(0.5, disp_weight.cumsum() / disp_weight.sum(), disp_weight.index.to_numpy())
)

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
# drawstyle="steps-mid": the x axis is a binned quantity, so a step reads as
# the histogram it is; a categorical bar plot would label every bin and the
# axis runs to hundreds of bins (docs/visualisations.md).
disp_weight.plot(ax=ax, drawstyle="steps-mid")
ax.set_ylabel("stranded weight")
fig_path = figure_dir / f"StrandingTravelDistance_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
print(f"weight-weighted median travel distance at stranding: {median_disp:.1f} km")
plt.show()

# %% [markdown]
# # Validation / summary

# %%
total_released = float(beaching["weight"].sum())
n_beached = float(beached["weight"].sum())
print(f"regime={regime}, hex_radius={hex_radius} m, season={season} "
      f"(months {season_months}), age_bin_days={age_bin_days}")
print(f"  drifters (Σweight): {total_released:,.0f}")
print(f"  beached:           {n_beached:,.0f} "
      f"({100 * n_beached / max(total_released, 1):.1f}%)")
print(f"  stranding hexes:   {beached['beach_hex'].nunique():,}")
print(f"  source hexes:      {frac['release_hex'].nunique():,}")
print(f"  median travel distance at stranding: {median_disp:,.1f} km")
print(f"  member:            {member}")
print("  cumulative beached fraction by age horizon (pooled, interannual min–max):")
for t in curve.index:
    print(f"    {t:>5} d: {100 * curve[t]:6.2f} %  "
          f"[{100 * per_year.loc[t].min():6.2f} .. {100 * per_year.loc[t].max():6.2f}]")
