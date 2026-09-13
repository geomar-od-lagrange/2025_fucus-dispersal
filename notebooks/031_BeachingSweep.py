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
# # Beaching parameter sweep: reporting a range, not a number
#
# Parquet-only consumer that compares the **members** of the beaching store
# written by `024e_BuildBeaching`. A member is one rate-model setting the
# reducer was run with — thresholds and timescales of the two-state hazard
# (`docs/beaching.md`) — identified by the opaque tag the reducer wrote into
# the store filename. This notebook takes an explicit list of tags plus a
# display label per tag and never parses either.
#
# One member = every `(year, month)` partition carrying that tag; partitions
# are additive over release_doy/month/year exactly as within a single member.
# In the Baltic the beaching scheme can dominate the answer, so the headline
# stranding number is only meaningful as a **range over the members**, next to
# the pattern statistics that discriminate between them: concentration of the
# stranded weight, how many hexes receive any, how long stranding takes, and
# how far the stranded material had travelled.

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
# Travel-distance bin width (km) of the store's `disp_bin` axis (must match
# the reducer).
disp_bin_km = 10.0
# 0 = pool all monthly partitions across years; 1..12 = that month only.
release_month = 0

# Member tags to compare, comma-separated, in the order they should be drawn.
# Each must have been built by 024e; missing members are reported and skipped
# rather than raising.
members_csv = "step_wc0p05_ts3_tcinf,step_wc0p1_ts3_tcinf,step_wc0p15_ts3_tcinf"
# Display labels, comma-separated, one per tag. Empty = use the tags.
labels_csv = ""
# The member treated as the baseline in the narrative summary; "" = none.
baseline_member = ""

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
member_tags = [x.strip() for x in members_csv.split(",") if x.strip()]
member_labels = [x.strip() for x in labels_csv.split(",") if x.strip()] or member_tags
assert len(member_labels) == len(member_tags), (
    f"{len(member_labels)} labels for {len(member_tags)} members"
)
mpl.rcParams["figure.dpi"] = fig_dpi_scale * mpl.rcParamsDefault["figure.dpi"]
# Hex seam stroke. edgecolor="face" means this is not a visible outline -- it
# closes the ~1 px anti-aliasing seam between adjacent polygons so the grid
# reads as a continuous field. The seam is a fixed PIXEL artifact, so a fixed
# point width makes the resulting hex dilation DPI-invariant (17.5% of hex
# width at every dpi). Pinning it to ~1 px instead lets dilation fall as
# resolution rises: ~8.6% at fig_dpi_scale=3 on the Baltic crop.
hex_seam_lw = 1.1 * 72 / (100 * fig_dpi_scale)

figure_dir = output_root / "Figures" / "031"
figure_dir.mkdir(parents=True, exist_ok=True)

store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

month_suffix = f"_m{release_month:02d}" if release_month else ""
month_re = rf"_m{release_month:02d}" if release_month else r"_m\d{2}"


# %% [markdown]
# # Pool each sweep member
#
# One member = every `(year, month)` partition carrying its tag.

# %%
def load_member(tag):
    part_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}}){month_re}"
        rf"_{re.escape(tag)}\.parquet$"
    )
    files = [
        f for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet"))
        if part_re.search(f.name)
    ]
    if not files:
        return None, 0
    return pd.concat(
        [pd.read_parquet(f).reset_index(drop=True) for f in files], ignore_index=True
    ), len(files)


members = {}
for tag, label in zip(member_tags, member_labels):
    df, n = load_member(tag)
    if df is None:
        print(f"{tag}: no partitions found — skipped")
        continue
    members[label] = df
    print(f"{tag} ({label}): {n} partitions, {len(df):,} rows")
if not members:
    raise FileNotFoundError(
        f"no sweep members found for regime {regime!r} at {store_root} — run 024e."
    )

# %% [markdown]
# # Per-member statistics
#
# `beached_fraction` is stranded weight over all released weight (the residual
# rows carry the never-beached remainder, so the denominator is the full
# release pool). `gini` measures how unevenly the stranded weight is spread
# over the coastal hexes that receive any: 0 = uniform, → 1 = concentrated.
# `median_age_days` and `median_travel_km` are weight-weighted medians of the
# store's `beach_age_bin` and `disp_bin` axes — how long the stranded material
# drifted, and how far it got.


# %%
def gini(values):
    """Gini concentration of a non-negative weight vector."""
    v = np.sort(np.asarray(values, dtype="float64"))
    v = v[v > 0]
    if v.size == 0:
        return np.nan
    n = v.size
    return float((2.0 * np.arange(1, n + 1) - n - 1).dot(v) / (n * v.sum()))


def weighted_median_bin(weights_per_bin, bin_width):
    """Weight-weighted median of a binned axis, at bin centres."""
    if weights_per_bin.empty:
        return np.nan
    centres = (weights_per_bin.index.to_numpy() + 0.5) * bin_width
    return float(np.interp(0.5, weights_per_bin.cumsum() / weights_per_bin.sum(), centres))


rows = []
for label, df in members.items():
    beached = df[df["beach_hex"] >= 0]
    per_hex = beached.groupby("beach_hex")["weight"].sum()
    total = float(df["weight"].sum())
    rows.append({
        "member": label,
        "beached_fraction": float(beached["weight"].sum()) / max(total, 1.0),
        "gini": gini(per_hex.to_numpy()),
        "beach_hexes": int(per_hex.size),
        "median_age_days": weighted_median_bin(
            beached.groupby("beach_age_bin")["weight"].sum(), age_bin_days
        ),
        "median_travel_km": weighted_median_bin(
            beached.groupby("disp_bin")["weight"].sum(), disp_bin_km
        ),
    })
stats = pd.DataFrame(rows).set_index("member")
print(stats.to_string(float_format=lambda v: f"{v:,.4f}"))

# %% [markdown]
# # The range across members
#
# One panel per statistic, members on a categorical x axis in the order given.
# Nothing about the member tags is ordinal — they differ in threshold,
# timescale, and edge width at once — so a categorical axis is the honest one:
# a numeric axis would invite reading a slope where there is only a list.

# %%
metrics = ["beached_fraction", "gini", "beach_hexes", "median_age_days", "median_travel_km"]
titles = {
    "beached_fraction": "beached fraction",
    "gini": "Gini of stranded weight",
    "beach_hexes": "stranding hexes",
    "median_age_days": "median age (d)",
    "median_travel_km": "median travel (km)",
}
fig, axes = plt.subplots(2, 3, layout="constrained")
for ax, m in zip(axes.flat, metrics):
    stats[m].plot(ax=ax, marker="o")
    # Panel titles rather than y labels: five stacked panels with long y labels
    # collide across columns at the default figure size (docs/visualisations.md).
    ax.set_title(titles[m])
    ax.set_xlabel("")
    ax.set_xticks(range(len(stats)))
    ax.set_xticklabels(stats.index, rotation=45, ha="right")
for ax in axes.flat[len(metrics):]:
    ax.set_axis_off()
fig_path = figure_dir / f"BeachingSweepStats_{regime}_r{hex_radius}m{month_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Where-stranded maps across members
#
# A shared `LogNorm` across every member, so the change is in the pattern and
# not in the colour scale.

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


def hex_gdf(df):
    grp = df[df["beach_hex"] >= 0].groupby("beach_hex")["weight"].sum()
    return (
        grp.rename("value").reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="beach_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


gdfs = {label: hex_gdf(df) for label, df in members.items()}
# A member can carry no stranded weight at all (a threshold above every
# forcing hour in the partition). Report and drop it rather than feeding an
# empty list to pd.concat.
empty = [label for label, g in gdfs.items() if g.empty]
for label in empty:
    print(f"  {label}: no stranded weight — skipped in the maps")
    del gdfs[label]
if not gdfs:
    raise ValueError(
        "no member has any stranded weight; nothing to map "
        f"({len(empty)} member(s) empty)"
    )
shared = pd.concat([g["value"] for g in gdfs.values()])
vmax = float(shared.max())
norm = LogNorm(vmin=max(float(shared[shared > 0].min()), vmax / 1e4), vmax=vmax)

ncols = len(gdfs)
fig, axes = plt.subplots(
    1, ncols, figsize=(panel_height_in * domain_aspect * ncols, panel_height_in),
    layout="constrained", squeeze=False,
)
for ax, (label, g) in zip(axes[0], gdfs.items()):
    g.plot(ax=ax, column="value", cmap=cmap, norm=norm, legend=True,
           edgecolor="face", linewidth=hex_seam_lw, zorder=1)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(label)
fig_path = figure_dir / f"BeachingSweepMaps_{regime}_r{hex_radius}m{month_suffix}.png"
fig.savefig(fig_path)
print(f"wrote {fig_path}")
plt.show()

# %% [markdown]
# # Summary

# %%
lo, hi = stats["beached_fraction"].min(), stats["beached_fraction"].max()
print(f"regime={regime}, hex_radius={hex_radius} m, "
      + (f"month={release_month}, " if release_month else "all months, ")
      + f"{len(stats)} sweep member(s)")
print(f"  beached fraction range: {100 * lo:.1f}% .. {100 * hi:.1f}% "
      f"(spread {100 * (hi - lo):.1f} points)")
baseline_label = dict(zip(member_tags, member_labels)).get(baseline_member)
if baseline_label in stats.index:
    print(f"  baseline {baseline_member}: "
          f"{100 * stats.loc[baseline_label, 'beached_fraction']:.1f}%")
print(f"  Gini range: {stats['gini'].min():.3f} .. {stats['gini'].max():.3f} "
      "(higher = stranding concentrated on fewer, wave-exposed hexes)")
print(f"  median stranding age: {stats['median_age_days'].min():.0f} .. "
      f"{stats['median_age_days'].max():.0f} d; median travel distance: "
      f"{stats['median_travel_km'].min():.0f} .. {stats['median_travel_km'].max():.0f} km")
