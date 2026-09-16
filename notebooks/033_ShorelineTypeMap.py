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
# # Shoreline type map, and whether stranding follows substrate
#
# Lightweight parquet-only consumer: the shoreline-class table from the data
# twin (`shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet`), the
# `024a` hex key, and the `024e` beaching store — no trajectory zarrs, no
# sidecar, no Dask.
#
# The substrate axis of the beaching rate model is the **flat fraction** `ff`:
# the share of classified shoreline that is low and unconsolidated (sand,
# shingle, marsh, mudflat) rather than a hard wall (rock, cliff, armoured
# bank) — see [docs/beaching.md](../docs/beaching.md). `024d` carries it per
# particle-hour as the value of the nearest coastal land cell; here it is
# aggregated the way the reader looks at it, per hex.
#
# 1. **Dominant shoreline type per hex** — the `seg_len_m`-weighted mean flat
#    fraction of the coastline sub-segments falling in the hex, plus the total
#    coastline length. A sub-segment with **no class attribution counts 0.5**,
#    the same convention `024d` applies (halfway between the two trap weights,
#    so unattributed coast is pushed neither towards trapping nor away from
#    it); the share of coastline that is unattributed is annotated on the map
#    and printed. Hexes with **no coastline at all carry no value** and render
#    absent — they are not "mixed".
# 2. **Stranded weight against shore type** — per hex the stranded weight from
#    the beaching store, its coastline-length-normalised value, and the
#    Spearman rank correlation against the hex flat fraction, globally and per
#    HELCOM subbasin, with the flat/wall split of stranded weight.
#
# Two rate-model members are shown side by side, and which is which is the
# whole point of the comparison:
#
# - **`member`** (default `step_wc0p1_ts3_tcinf_tf1_tw0`) — the production
#   setting. `trap_flat = 1`, `trap_wall = 0`, so the realised strong-wave
#   factor is `trap = ff/100`: **substrate genuinely drives the trapping
#   term**, and a fully walled hex (`ff = 0`) strands nothing at all.
# - **`contrast_member`** (default `step_wc0p1_ts3_tcinf`) — the degenerate
#   case, `trap_flat = trap_wall = 1`, where **substrate does not enter the
#   rate at all**. Any relationship its panels show is pure geography: where
#   the coastline happens to be flat against where the drift and the waves put
#   the material.
#
# Reading the two together separates "stranding sits on flat coast because the
# model puts it there" from "stranding sits on flat coast because that is where
# the drift goes". Every figure and every printed block names which member it
# is and whether substrate is in its rate.
#
# Figures are print-ready (full page width, 300 dpi) and the per-hex shore-type
# table is exported as GeoJSON.

# %%
import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pyproj import Transformer
from scipy.stats import spearmanr
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth

from hextraj import HexProj

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store (key + beaching partitions).
output_root = "../output"
# Read root of the data twin (the shoreline-class table).
data_root = "../data"
# Shoreline-class table (sub-segment midpoints in EPSG:3035 with a flat
# fraction each). "" = the data twin's copy under `data_root`, exactly as 024d.
shoreclass_path = ""

# Which regime's beaching partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024e). Must match the files on disk.
hex_radius = 6000
# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year.
season = "ALL"

# Production rate-model member: the opaque tag the reducer put in the store
# filename. Substrate drives the trapping term here (trap_flat = 1,
# trap_wall = 0, so trap = ff/100 and pure-wall coast strands nothing).
member = "step_wc0p1_ts3_tcinf_tf1_tw0"
# The contrast: a member whose trap term is degenerate (trap_flat = trap_wall
# = 1), so substrate does not enter the rate at all. "" skips the second pass.
contrast_member = "step_wc0p1_ts3_tcinf"

# Hexes holding less coastline than this are dropped from the correlations and
# the scatter: a length-weighted mean over a sliver of coast is noise, and its
# per-km stranded weight divides by almost nothing.
min_coast_len_m = 1000.0
# Hexes whose coastline is less than this fraction (by length) inside the
# Baltic are dropped from the fits: the key tiles the whole BSH domain, and
# North Sea / Skagerrak coast the Baltic releases never reach would let "no
# Fucus released within reach" masquerade as a substrate signal.
min_baltic_fraction = 0.5
# Subbasins with fewer coastal hexes than this get no per-subbasin Spearman.
min_hexes_per_subbasin = 10
# Flat-fraction bin edges of the binned summary panel (wall → flat).
flat_bin_edges_csv = "0.0,0.2,0.4,0.6,0.8,1.0"
# Subbasins fronting the Swedish east coast, for the flat/wall split readout.
swedish_subbasins_csv = (
    "Western Gotland Basin,Eastern Gotland Basin,"
    "Northern Baltic Proper,Bothnian Sea"
)

# Map extent (degrees), as 029 — the key tiles the whole BSH domain including
# the North Sea, which is empty for Fucus. Set all four to 0 for the key's
# full bounds.
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
shoreclass_path = Path(
    shoreclass_path
    or data_root / "shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet"
)
if not shoreclass_path.exists():
    raise FileNotFoundError(f"shoreline-class table missing: {shoreclass_path}")

if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]

flat_bin_edges = [float(x) for x in flat_bin_edges_csv.split(",") if x]
swedish_subbasins = [s.strip() for s in swedish_subbasins_csv.split(",") if s.strip()]
members = [member] + ([contrast_member] if contrast_member else [])
# Whether each member's rate reacts to substrate at all — stated on every
# figure and in every printed block, because it is what the comparison means.
SUBSTRATE_IN_RATE = "substrate IS in the rate (trap = ff/100)"
SUBSTRATE_INERT = "substrate is NOT in the rate (trap = 1)"
member_note = {member: SUBSTRATE_IN_RATE}
if contrast_member:
    member_note[contrast_member] = SUBSTRATE_INERT

# Print-ready output: every saved figure is one full text-width (180 mm)
# figure at 300 dpi, so panels land in the manuscript at their final size
# (rationale in docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300
# Hex seam stroke, as 029: edgecolor="face" plus a ~1 px stroke at the output
# dpi closes the anti-aliasing seam between adjacent polygons.
hex_seam_lw = 1.1 * 72 / FIGURE_DPI

figure_dir = output_root / "Figures" / "033"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "033"
export_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Key + hex projection
#
# The key's sidecar carries the `HexProj` every `024x` stage labels with, and
# the subbasin id → name map the per-subbasin readouts use.

# %%
store_root = output_root / "HexAggregates"
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
if not key_path.exists() or not meta_path.exists():
    raise FileNotFoundError(
        f"Key file or sidecar missing — run 024a_BuildHexKey.md first.\n"
        f"  expected: {key_path}\n  expected: {meta_path}"
    )
meta = json.loads(meta_path.read_text())
hp = HexProj(**meta["hex_proj"])
subbasin_name = {int(k): v for k, v in meta["subbasin_id_to_name"].items()}
key = gpd.read_parquet(key_path)
print(f"key: {len(key):,} hexes; HexProj {meta['hex_proj']}")

missing = [s for s in swedish_subbasins if s not in set(subbasin_name.values())]
if missing:
    raise ValueError(
        f"subbasin name(s) {missing} not in the key's sidecar; "
        f"available: {sorted(set(subbasin_name.values()))}"
    )

# %% [markdown]
# # Per-hex shore type
#
# The table is sub-segment midpoints in EPSG:3035 (`024d` snaps them onto its
# own 3035 raster without reprojecting). The hex key lives in lon/lat and is
# labelled through `HexProj`, so the midpoints are reprojected 3035 → 4326 —
# the same `pyproj` transform `024d` uses for the reverse direction — and then
# labelled. Per hex: the `seg_len_m`-weighted mean flat fraction, the total
# coastline length, and the length share that carried no class attribution.
#
# `NaN` flat fraction counts as 0.5, as everywhere else in the pipeline
# (docs/beaching.md). The unattributed share is kept as its own column so the
# map can say how much of what it shows is that convention rather than data.

# %%
sc = pd.read_parquet(shoreclass_path)
lon_sc, lat_sc = Transformer.from_crs(3035, 4326, always_xy=True).transform(
    sc["x_3035"].to_numpy("float64"), sc["y_3035"].to_numpy("float64")
)
sc = sc.assign(
    hex_id=hp.label(lon_sc, lat_sc).astype("int32"),
    unattributed=~np.isfinite(sc["flat_fraction"]),
    ff=sc["flat_fraction"].fillna(0.5),
)
sc = sc[np.isfinite(sc["seg_len_m"]) & (sc["seg_len_m"] > 0)]

grp = sc.assign(
    _ff_len=sc["ff"] * sc["seg_len_m"],
    _unatt_len=sc["unattributed"] * sc["seg_len_m"],
    _baltic_len=sc["in_baltic"] * sc["seg_len_m"],
).groupby("hex_id")
shore = pd.DataFrame(
    {
        "coast_len_m": grp["seg_len_m"].sum(),
        "flat_fraction": grp["_ff_len"].sum() / grp["seg_len_m"].sum(),
        "unattributed_fraction": grp["_unatt_len"].sum() / grp["seg_len_m"].sum(),
        "baltic_fraction": grp["_baltic_len"].sum() / grp["seg_len_m"].sum(),
    }
).reset_index()
# Inner join: a hex the key does not carry has no geometry, and a hex no
# sub-segment falls in gets no row at all — absent, not 0.5.
shore = shore.merge(
    key[["hex_id", "geometry", "helcom_subbasin", "water_area_m2"]], on="hex_id"
)
shore["subbasin"] = shore["helcom_subbasin"].map(subbasin_name)

unattributed_share = float(
    sc.loc[sc["unattributed"], "seg_len_m"].sum() / sc["seg_len_m"].sum()
)
print(f"shoreline sub-segments: {len(sc):,}; total length "
      f"{sc['seg_len_m'].sum() / 1e3:,.0f} km")
print(f"unattributed by length: {100 * unattributed_share:.1f}% (counted as ff = 0.5)")
print(f"hexes carrying coastline: {len(shore):,} of {len(key):,} key hexes")

# %% [markdown]
# # Pool the season's beaching partitions across years
#
# Same layout and filename contract as 029: flat per-(year, month) partitions
# `HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet` under
# `output_root/HexAggregates/`. A whole-year partition (no `_mMM`) carries
# every month, so it is only read for `season = "ALL"`.
#
# Each partition is reduced to `(beach_hex, shore_type) → weight` **before**
# pooling, and the never-beached residual (`beach_hex = -1`) is dropped there:
# this notebook has no release-side, age or displacement axis, and a full
# ten-year pool of the raw rows is ~31 M rows per member for a table that
# reduces to a few thousand.

# %%
def load_beaching(tag):
    """Stranded weight per (beach_hex, shore_type), pooled over the season."""
    pat = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
        rf"_{re.escape(tag)}\.parquet$"
    )
    files = []
    for f in sorted(
        store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")
    ):
        m = pat.search(f.name)
        if m is None:
            continue
        month = m.group(2)
        if month is None:
            if season != "ALL":
                continue
        elif int(month) not in season_months:
            continue
        files.append((f, int(m.group(1))))
    if not files:
        raise FileNotFoundError(
            f"no beaching partitions for regime {regime!r}, season {season}, "
            f"member {tag!r} at {store_root} — run 024e."
        )
    parts = []
    for f, _ in files:
        d = pd.read_parquet(f, columns=["beach_hex", "shore_type", "weight"])
        d = d[d["beach_hex"] >= 0]
        parts.append(d.groupby(["beach_hex", "shore_type"], observed=True)["weight"].sum())
    out = (
        pd.concat(parts).groupby(level=[0, 1], observed=True).sum()
        .rename("weight").reset_index()
    )
    print(f"{tag}: pooled {len(files)} partition(s) over release years "
          f"{sorted({y for _, y in files})} → {len(out):,} (hex, shore_type) rows, "
          f"Σweight {out['weight'].sum():,.0f}")
    return out


beaching = {tag: load_beaching(tag) for tag in members}

# %% [markdown]
# # Rendering helpers
#
# 029's hex rendering verbatim — hex polygons + Natural Earth coastline on a
# plain EPSG:4326 axis, colorbar as an inset axes at axes-fraction height 1.0
# so it is exactly the map height — with a `cmap` passthrough for the one
# override this notebook defends (below).

# %%
def hex_map(gdf, ax, norm=None, title=None, label=None, cmap=None):
    """Hex choropleth of `gdf.value` + coastline on a lon/lat axis."""
    if not gdf.empty:
        cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
        gdf.plot(
            ax=ax, column="value", norm=norm, cmap=cmap, legend=True, cax=cax,
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
    width; constrained_layout does the packing, no measure-rescale loop."""
    return 1.10 * nrows * (0.72 * width_in / ncols) / aspect


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
# # Dominant shoreline type
#
# The length-weighted mean flat fraction per hex, on a fixed linear 0–1 scale.
#
# **One styling override, `cmap="coolwarm"` with `Normalize(0, 1)`** (AGENTS.md
# otherwise forbids `cmap=`): this is not a magnitude but a bounded two-ended
# axis — 0 is fully walled shore, 1 fully flat shore — whose **midpoint 0.5 is
# a distinguished value**, being exactly what the pipeline assigns coastline it
# has no class for. The default sequential map has no distinguished midpoint,
# so "half flat" and "unattributed" would land in an arbitrary mid-ramp colour
# and the two ends would not read as opposites. Fixing the norm to [0, 1]
# rather than auto-ranging keeps the white at 0.5 across runs.

# %%
gdf_shore = shore.assign(value=shore["flat_fraction"]).pipe(
    gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326"
)

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(1, 1, domain_aspect, FIGURE_WIDTH_IN)
)
hex_map(
    gdf_shore, ax, norm=Normalize(vmin=0.0, vmax=1.0), cmap="coolwarm",
    title=(
        f"shoreline type per hex "
        f"({100 * unattributed_share:.0f}% of coastline unattributed, counted as 0.5)"
    ),
    label="flat fraction (0 = wall, 1 = flat)",
)
fig_path = figure_dir / f"ShorelineType_r{hex_radius}m.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Coastline length per hex
#
# The denominator of the normalised stranded weight below, and the map that
# says how much of the flat-fraction map is carried by a sliver of coast.

# %%
fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(1, 1, domain_aspect, FIGURE_WIDTH_IN)
)
hex_map(
    shore.assign(value=shore["coast_len_m"] / 1e3).pipe(
        gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326"
    ),
    ax, title="coastline length per hex", label="coastline length (km)",
)
fig_path = figure_dir / f"CoastlineLength_r{hex_radius}m.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Stranded weight against shore type
#
# Per coastal hex: the stranded weight the member deposits there, and the same
# divided by the hex's coastline length. The length-normalised value is the
# fair comparison — a hex with 40 km of coast catches more than one with 4 km
# whatever it is made of — so it carries the headline Spearman; the raw weight
# is reported beside it.
#
# The fit set is the coastal hexes whose coastline is **at least
# `min_baltic_fraction` inside the Baltic** (the table's own `in_baltic` flag)
# and at least `min_coast_len_m` long. The Baltic cut uses `in_baltic` rather
# than the key's `helcom_subbasin >= 0`: the HELCOM polygons leave 859 genuine
# Baltic coastal hexes unassigned (a coastal hex centred on land or inside an
# archipelago falls outside every polygon), and dropping those would throw away
# about half the stranded weight. The per-subbasin panel below necessarily
# keeps only the named subbasins; how much of the fit set that leaves out is
# printed.
#
# Panels: the scatter with its Spearman; the binned split of stranded weight
# and coastline length over flat-fraction bins — **if stranding merely followed
# the coastline's own composition the two series would coincide**, so the gap
# between them is the answer; and the per-subbasin Spearman, full figure width
# because 17 HELCOM names need the room (028's reasoning), so a regional signal
# is not averaged away.

# %%
def shore_stranding(df):
    """Per coastal hex: stranded weight, and weight per km of coastline."""
    w = df.groupby("beach_hex")["weight"].sum().rename("stranded_weight")
    t = shore.merge(w, left_on="hex_id", right_index=True, how="left")
    t["stranded_weight"] = t["stranded_weight"].fillna(0.0)
    t["stranded_per_km"] = t["stranded_weight"] / (t["coast_len_m"] / 1e3)
    return t


def rho(x, y):
    """Spearman rank correlation, NaN where it is not defined."""
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(spearmanr(x, y).statistic)


strand = {tag: shore_stranding(beaching[tag]) for tag in members}
bin_labels = [
    f"{lo:.1f}–{hi:.1f}" for lo, hi in zip(flat_bin_edges[:-1], flat_bin_edges[1:])
]
summary = {}

for tag in members:
    t = strand[tag]
    in_fit = (t["coast_len_m"] >= min_coast_len_m) & (
        t["baltic_fraction"] >= min_baltic_fraction
    )
    fit = t[in_fit]
    dropped_weight_share = float(
        t.loc[~in_fit, "stranded_weight"].sum() / max(t["stranded_weight"].sum(), 1e-9)
    )
    unnamed_weight_share = float(
        fit.loc[fit["helcom_subbasin"] < 0, "stranded_weight"].sum()
        / max(fit["stranded_weight"].sum(), 1e-9)
    )
    rho_per_km = rho(fit["flat_fraction"], fit["stranded_per_km"])
    rho_raw = rho(fit["flat_fraction"], fit["stranded_weight"])

    summed = (
        fit.assign(
            flat_bin=pd.cut(fit["flat_fraction"], flat_bin_edges, include_lowest=True)
        )
        .groupby("flat_bin", observed=False)[["stranded_weight", "coast_len_m"]]
        .sum()
    )
    binned = pd.DataFrame(
        {
            "stranded weight": summed["stranded_weight"] / summed["stranded_weight"].sum(),
            "coastline length": summed["coast_len_m"] / summed["coast_len_m"].sum(),
        }
    )

    per_sub = (
        fit[fit["helcom_subbasin"] >= 0].groupby("subbasin")
        .apply(
            lambda g: pd.Series(
                {
                    "rho": rho(g["flat_fraction"], g["stranded_per_km"]),
                    "n_hex": len(g),
                    "stranded_weight": g["stranded_weight"].sum(),
                    "flat_fraction": np.average(
                        g["flat_fraction"], weights=g["coast_len_m"]
                    ),
                }
            ),
            include_groups=False,
        )
        .dropna(subset=["rho"])
    )
    per_sub = per_sub[per_sub["n_hex"] >= min_hexes_per_subbasin].sort_values("rho")
    n_zero = int((fit["stranded_weight"] == 0).sum())
    summary[tag] = {
        "rho_per_km": rho_per_km, "rho_raw": rho_raw, "binned": binned,
        "per_sub": per_sub, "fit": fit, "n_zero": n_zero,
        "dropped_weight_share": dropped_weight_share,
        "unnamed_weight_share": unnamed_weight_share,
    }

    # Mosaic rather than a 1x3 row: the subbasin panel carries 17 long names on
    # its y axis and needs the full page width, while the scatter and the bars
    # read fine at half width each.
    fig, axd = plt.subplot_mosaic(
        [["scatter", "binned"], ["subbasin", "subbasin"]], layout="constrained"
    )
    fig.set_size_inches(FIGURE_WIDTH_IN, 0.85 * FIGURE_WIDTH_IN)

    ax = axd["scatter"]
    ax.scatter(fit["flat_fraction"], fit["stranded_per_km"], s=2, alpha=0.3)
    ax.set_yscale("log")
    ax.set_xlabel("hex flat fraction")
    ax.set_ylabel("stranded weight per km of coast")
    # Short titles throughout: at 180 mm / 3 panels the default font cannot
    # shrink, so a long string overruns its panel (028 hit the same).
    ax.set_title(f"ρ = {rho_per_km:+.3f}, n = {len(fit):,}")

    ax = axd["binned"]
    binned.plot.bar(ax=ax)
    ax.set_xlabel("hex flat fraction")
    ax.set_ylabel("share of total")
    ax.set_title("stranding vs. coastline")
    ax.set_xticklabels(bin_labels, rotation=45, ha="right")

    ax = axd["subbasin"]
    ax.scatter(per_sub["rho"], np.arange(len(per_sub)))
    ax.set_yticks(np.arange(len(per_sub)))
    ax.set_yticklabels(per_sub.index)
    ax.axvline(0.0, linewidth=0.5)
    ax.set_xlabel(f"Spearman ρ (≥ {min_hexes_per_subbasin} coastal hexes)")

    fig.suptitle(
        f"{regime}, r{hex_radius} m, season {season}\n"
        f"{tag}\n{member_note[tag]}"
    )
    fig_path = (
        figure_dir
        / f"ShoreTypeVsStranding_{regime}_r{hex_radius}m_{season}_{tag}.png"
    )
    # Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

# %% [markdown]
# # Flat/wall split of stranded weight
#
# The store's own `shore_type` column — `flat` where the stranding site's `ff ≥
# 50` at the stranding hour, else `wall` — summed per subbasin of the stranding
# hex, for the Swedish-coast subbasins specifically and for the whole domain.
# The label is a threshold while the trap factor is linear in `ff`, so under
# the production member wall-labelled sites still strand (at `ff/100` of the
# flat rate) and only `ff = 0` is fully non-trapping; under the contrast member
# the label moves no weight at all and the split is the coastline's own
# composition as the drift samples it.

# %%
def flat_wall_split(df):
    """Stranded weight by (subbasin of the stranding hex, shore_type)."""
    per_hex = df.merge(
        key[["hex_id", "helcom_subbasin"]], left_on="beach_hex", right_on="hex_id"
    )
    per_hex["subbasin"] = per_hex["helcom_subbasin"].map(subbasin_name)
    tab = per_hex.pivot_table(
        index="subbasin", columns="shore_type", values="weight", aggfunc="sum",
    ).fillna(0.0)
    for col in ("flat", "wall"):
        if col not in tab:
            tab[col] = 0.0
    tab["wall_share"] = tab["wall"] / (tab["flat"] + tab["wall"])
    return tab


splits = {tag: flat_wall_split(beaching[tag]) for tag in members}

# %% [markdown]
# # GeoJSON export
#
# The per-hex shore-type table with the 024a geometry, EPSG:4326 — the map
# above in a form collaborators can open in R or QGIS, with the production
# member's stranded weight joined on.

# %%
export = (
    strand[member][
        ["hex_id", "flat_fraction", "coast_len_m", "stranded_weight", "geometry"]
    ]
    .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
)
export_path = export_dir / f"033_{regime}_{season}_{member}.geojson"
export.to_file(export_path, driver="GeoJSON")
print(f"wrote {export_path} ({len(export):,} hexes)")

# %% [markdown]
# # Validation / summary

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, season={season} "
      f"(months {season_months})")
print(f"  coastline sub-segments: {len(sc):,}, "
      f"{sc['seg_len_m'].sum() / 1e3:,.0f} km total")
print(f"  unattributed by length: {100 * unattributed_share:.1f}% "
      f"(counted as flat_fraction = 0.5, as 024d)")
print(f"  hexes with coastline:   {len(shore):,} "
      f"(>= {min_coast_len_m:,.0f} m: {int((shore['coast_len_m'] >= min_coast_len_m).sum()):,})")
print(f"  length-weighted mean flat fraction over the domain: "
      f"{np.average(shore['flat_fraction'], weights=shore['coast_len_m']):.3f}")

for tag in members:
    s = summary[tag]
    t = strand[tag]
    tot = float(t["stranded_weight"].sum())
    print(f"\nmember {tag}  ({member_note[tag]})")
    print(f"  stranded weight on coastal hexes: {tot:,.0f} "
          f"over {int((t['stranded_weight'] > 0).sum()):,} hexes")
    print(f"  Spearman(flat_fraction, stranded weight per km coast): "
          f"{s['rho_per_km']:+.3f}")
    print(f"  Spearman(flat_fraction, stranded weight):              "
          f"{s['rho_raw']:+.3f}")
    print(f"  fit set: {len(s['fit']):,} Baltic coastal hexes "
          f"({s['n_zero']:,} of them strand nothing, so they are off the log "
          f"axis of the scatter but in every statistic)")
    print(f"  stranded weight outside the fit set (coastline < "
          f"{100 * min_baltic_fraction:.0f} % Baltic, or < "
          f"{min_coast_len_m:,.0f} m long): "
          f"{100 * s['dropped_weight_share']:.1f} %")
    print(f"  of the fit set, weight on hexes HELCOM leaves unnamed "
          f"(absent from the per-subbasin panel): "
          f"{100 * s['unnamed_weight_share']:.1f} %")
    print("  share of stranded weight / of coastline length, per flat-fraction bin:")
    for lbl, (_, row) in zip(bin_labels, s["binned"].iterrows()):
        print(f"    {lbl}: {100 * row['stranded weight']:5.1f} % weight   "
              f"{100 * row['coastline length']:5.1f} % coast")
    print("  per-subbasin Spearman (flat fraction vs. stranded weight per km):")
    for name, row in s["per_sub"].iterrows():
        mark = " *" if name in swedish_subbasins else ""
        print(f"    {name:<26} rho={row['rho']:+.3f}  n={int(row['n_hex']):4d}  "
              f"mean ff={row['flat_fraction']:.3f}{mark}")

    split = splits[tag]
    dom_flat = float(split["flat"].sum())
    dom_wall = float(split["wall"].sum())
    print(f"  flat/wall split of stranded weight, whole domain: "
          f"flat {100 * dom_flat / (dom_flat + dom_wall):.1f} % / "
          f"wall {100 * dom_wall / (dom_flat + dom_wall):.1f} %")
    swe = split.reindex(swedish_subbasins).dropna(how="all")
    swe_flat = float(swe["flat"].sum())
    swe_wall = float(swe["wall"].sum())
    print(f"  Swedish-coast subbasins ({', '.join(swedish_subbasins)}):")
    for name, row in swe.iterrows():
        print(f"    {name:<26} flat {100 * row['flat'] / (row['flat'] + row['wall']):5.1f} % / "
              f"wall {100 * row['wall'] / (row['flat'] + row['wall']):5.1f} %   "
              f"(Σ weight {row['flat'] + row['wall']:,.0f})")
    print(f"    {'pooled':<26} flat {100 * swe_flat / (swe_flat + swe_wall):5.1f} % / "
          f"wall {100 * swe_wall / (swe_flat + swe_wall):5.1f} %")
    swe_hexes = t[t["subbasin"].isin(swedish_subbasins)]
    lo = swe_hexes[swe_hexes["flat_fraction"] < 0.5]["stranded_weight"].sum()
    hi = swe_hexes["stranded_weight"].sum()
    print(f"    stranded weight on wall-like hexes (hex flat fraction < 0.5): "
          f"{100 * lo / max(hi, 1e-9):.1f} % of the Swedish-coast total")
