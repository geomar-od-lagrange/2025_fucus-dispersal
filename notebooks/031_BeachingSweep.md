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

# Beaching parameter sweep: reporting a range, not a number

Parquet-only consumer that compares the **members** of the beaching store
written by `024e_BuildBeaching`. A member is one rate-model setting the
reducer was run with — thresholds and timescales of the two-state hazard
(`docs/beaching.md`) — identified by the opaque tag the reducer wrote into
the store filename. This notebook takes an explicit list of tags plus a
display label per tag and never parses either.

One member = every partition carrying that tag whose release month falls in
the run's `season`; partitions are additive over release_doy/month/year
exactly as within a single member. In the Baltic the beaching scheme can
dominate the answer, so the headline stranding number is only meaningful as
a **range over the members**, next to the pattern statistics that
discriminate between them: concentration of the stranded weight, how many
hexes receive any, how long stranding takes, and how far the stranded
material had travelled.

Each statistic additionally carries the **2016–2019 interannual min–max** as
an error bar, so the spread across members is read against the spread the
same member shows across years. Figures are print-ready (full page width,
300 dpi).

```python
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
```

# Parameters

```python tags=["parameters"]
output_root = "../output"
regime = "surface_stokes"
hex_radius = 6000
age_bin_days = 10
# Travel-distance bin width (km) of the store's `disp_bin` axis (must match
# the reducer).
disp_bin_km = 10.0
# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year.
season = "ALL"

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
# Print-ready figure geometry: full page width in mm and the raster dpi.
fig_width_mm = 180
fig_dpi = 300
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
season_months = SEASON_MONTHS[season]
member_tags = [x.strip() for x in members_csv.split(",") if x.strip()]
member_labels = [x.strip() for x in labels_csv.split(",") if x.strip()] or member_tags
assert len(member_labels) == len(member_tags), (
    f"{len(member_labels)} labels for {len(member_tags)} members"
)

# Print-ready raster: figures go into the manuscript at a fixed column width,
# so figure.dpi is pinned to the savefig dpi (deliberate override of the
# AGENTS.md "no dpi/figsize" default — see docs/visualisations.md).
mpl.rcParams["figure.dpi"] = fig_dpi
fig_width_in = fig_width_mm / 25.4
# Hex seam stroke. edgecolor="face" means this is not a visible outline -- it
# closes the ~1 px anti-aliasing seam between adjacent polygons so the grid
# reads as a continuous field. The seam is a fixed PIXEL artifact, so pinning
# the stroke to ~1 px at the output dpi keeps the seam closed while letting
# the resulting hex dilation shrink as resolution rises.
hex_seam_lw = 1.1 * 72 / fig_dpi

figure_dir = output_root / "Figures" / "031"
figure_dir.mkdir(parents=True, exist_ok=True)

store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")
```

# Pool each sweep member

One member = every partition carrying its tag whose release month is in the
season. The release year is parsed from the filename and kept as a column,
because the interannual spread below groups by it. A whole-year partition
(no `_mMM` suffix) carries every month, so it is only read for
`season = "ALL"`.

```python
def load_member(tag):
    part_re = re.compile(
        rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
        rf"_{re.escape(tag)}\.parquet$"
    )
    parts = []
    for f in sorted(store_root.glob(f"HexAgg_beaching_r{hex_radius}m_{regime}_*.parquet")):
        m = part_re.search(f.name)
        if m is None:
            continue
        month = m.group(2)
        if month is None:
            if season != "ALL":
                continue
        elif int(month) not in season_months:
            continue
        parts.append(
            pd.read_parquet(f).reset_index(drop=True).assign(release_year=int(m.group(1)))
        )
    if not parts:
        return None, 0
    return pd.concat(parts, ignore_index=True), len(parts)


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
        f"no sweep members found for regime {regime!r}, season {season} at "
        f"{store_root} — run 024e."
    )
release_years = sorted(
    int(y) for y in pd.concat([df["release_year"] for df in members.values()]).unique()
)
print(f"season {season} (months {season_months}); release years {release_years}")
```

# Per-member statistics

`beached_fraction` is stranded weight over all released weight (the residual
rows carry the never-beached remainder, so the denominator is the full
release pool). `gini` measures how unevenly the stranded weight is spread
over the coastal hexes that receive any: 0 = uniform, → 1 = concentrated.
`median_age_days` and `median_travel_km` are weight-weighted medians of the
store's `beach_age_bin` and `disp_bin` axes — how long the stranded material
drifted, and how far it got.

The same reduction is applied per `release_year`, giving the interannual
min–max drawn as an error bar below.


```python
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


def member_stats(df):
    beached = df[df["beach_hex"] >= 0]
    per_hex = beached.groupby("beach_hex")["weight"].sum()
    total = float(df["weight"].sum())
    return {
        "beached_fraction": float(beached["weight"].sum()) / max(total, 1.0),
        "gini": gini(per_hex.to_numpy()),
        "beach_hexes": int(per_hex.size),
        "median_age_days": weighted_median_bin(
            beached.groupby("beach_age_bin")["weight"].sum(), age_bin_days
        ),
        "median_travel_km": weighted_median_bin(
            beached.groupby("disp_bin")["weight"].sum(), disp_bin_km
        ),
    }


stats = pd.DataFrame(
    [{"member": label, **member_stats(df)} for label, df in members.items()]
).set_index("member")
# Per-(member, year), for the interannual min–max error bars.
stats_year = pd.DataFrame(
    [
        {"member": label, "release_year": y, **member_stats(g)}
        for label, df in members.items()
        for y, g in df.groupby("release_year")
    ]
).set_index(["member", "release_year"])
print(stats.to_string(float_format=lambda v: f"{v:,.4f}"))
```

# The range across members

One panel per statistic, members on a categorical axis in the order given.
Nothing about the member tags is ordinal — they differ in threshold,
timescale, and edge width at once — so a categorical axis is the honest one:
a numeric axis would invite reading a slope where there is only a list. The
categories run down the **y** axis: member labels are long, and horizontal
labels on one shared left column cost a fraction of what rotated labels
under every panel do.

The bar through each marker is the interannual min–max of the same
statistic, so a member-to-member difference smaller than its own bar is not
a difference.

```python
metrics = ["beached_fraction", "gini", "beach_hexes", "median_age_days", "median_travel_km"]
titles = {
    "beached_fraction": "beached fraction",
    "gini": "Gini of stranded weight",
    "beach_hexes": "stranding hexes",
    "median_age_days": "median age (d)",
    "median_travel_km": "median travel (km)",
}
# Members top-to-bottom in the given order.
y = np.arange(len(stats))[::-1]
# Interannual bars take the second default-cycle colour so they read as a
# second series against the pooled markers (no literal colour name).
cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

ncols = 3
nrows = -(-len(metrics) // ncols)
fig, axes = plt.subplots(nrows, ncols, layout="constrained", squeeze=False)
fig.set_size_inches(fig_width_in, 0.22 * nrows * fig_width_in)
for k, m in enumerate(metrics):
    ax = axes.flat[k]
    lo = stats_year[m].groupby("member").min().reindex(stats.index)
    hi = stats_year[m].groupby("member").max().reindex(stats.index)
    # Drawn as an explicit lo-hi bar rather than a symmetric error offset: the
    # pooled value need not lie inside the per-year range (beach_hexes pools
    # as a union of hexes, so it exceeds every single year's count).
    ax.plot(stats[m].to_numpy(), y, marker="o")
    ax.hlines(y, lo.to_numpy(), hi.to_numpy(), color=cycle[1])
    # Panel titles rather than x labels: the quantities share no units, and a
    # title sits clear of the shared member axis (docs/visualisations.md).
    ax.set_title(titles[m])
    ax.set_yticks(y)
    ax.set_yticklabels(stats.index if k % ncols == 0 else [])
for ax in axes.flat[len(metrics):]:
    ax.set_axis_off()
axes.flat[len(metrics)].text(
    0.0, 0.5,
    f"bars: interannual min-max\n{release_years[0]}-{release_years[-1]}\n"
    f"season {season}",
    va="center",
)
fig_path = figure_dir / f"BeachingSweepStats_{regime}_r{hex_radius}m_{season}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=fig_dpi)
print(f"wrote {fig_path}")
plt.show()
```

# Where-stranded maps across members

A shared `LogNorm` across every member, so the change is in the pattern and
not in the colour scale. Each panel's colorbar is an inset axes at
axes-fraction height 1.0, exactly the map height.

```python
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


def grid_height_in(nrows, ncols, aspect, width_in):
    """Height seed for a grid of fixed-aspect map panels at a fixed figure
    width. The figure width is set by the page, so only the height is free:
    each column gives its map ~72 % of its width (the rest is the inset
    colorbar and its tick labels), and ~10 % is added for titles and padding.
    constrained_layout does the packing; there is no measure-rescale loop."""
    return 1.10 * nrows * (0.72 * width_in / ncols) / aspect


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
fig, axes = plt.subplots(1, ncols, layout="constrained", squeeze=False)
fig.set_size_inches(
    fig_width_in, grid_height_in(1, ncols, domain_aspect, fig_width_in)
)
for ax, (label, g) in zip(axes[0], gdfs.items()):
    cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
    g.plot(ax=ax, column="value", cmap=cmap, norm=norm, legend=True, cax=cax,
           legend_kwds={"label": "stranded weight (particles)"},
           edgecolor="face", linewidth=hex_seam_lw, zorder=1)
    coast.plot(ax=ax, color="black", linewidth=0.5, zorder=2)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect(1 / np.cos(np.radians(0.5 * (extent[2] + extent[3]))))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(label)
fig_path = figure_dir / f"BeachingSweepMaps_{regime}_r{hex_radius}m_{season}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=fig_dpi)
print(f"wrote {fig_path}")
plt.show()
```

# Summary

```python
lo, hi = stats["beached_fraction"].min(), stats["beached_fraction"].max()
print(f"regime={regime}, hex_radius={hex_radius} m, season={season} "
      f"(months {season_months}), {len(stats)} sweep member(s)")
print(f"  beached fraction range across members: {100 * lo:.1f}% .. {100 * hi:.1f}% "
      f"(spread {100 * (hi - lo):.1f} points)")
bf_year = stats_year["beached_fraction"]
print(f"  interannual min–max within a member ({release_years[0]}–{release_years[-1]}):")
for label in stats.index:
    print(f"    {label}: {100 * stats.loc[label, 'beached_fraction']:.1f}% "
          f"[{100 * bf_year.loc[label].min():.1f} .. {100 * bf_year.loc[label].max():.1f}]")
baseline_label = dict(zip(member_tags, member_labels)).get(baseline_member)
if baseline_label in stats.index:
    print(f"  baseline {baseline_member}: "
          f"{100 * stats.loc[baseline_label, 'beached_fraction']:.1f}%")
print(f"  Gini range: {stats['gini'].min():.3f} .. {stats['gini'].max():.3f} "
      "(higher = stranding concentrated on fewer, wave-exposed hexes)")
print(f"  median stranding age: {stats['median_age_days'].min():.0f} .. "
      f"{stats['median_age_days'].max():.0f} d; median travel distance: "
      f"{stats['median_travel_km'].min():.0f} .. {stats['median_travel_km'].max():.0f} km")
```
