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
partitions are pooled across every available release year.

The rate model is a two-state hazard — a background rate in the near-shore
band plus a strong-wave rate switched on where the onshore Stokes forcing
exceeds a threshold (`docs/beaching.md`). Its settings are baked into the
reducer's opaque `member` tag, which selects the partitions read here and
labels every figure; this notebook never parses the tag. Draws:

1. **Where-stranded density** — beached weight (expected particles) per
   `beach_hex` (log scale).
2. **Beached fraction per source hex** — of the drifters released in each
   hex, what fraction strands within the viability window (linear 0–1).
3. **Beaching age horizons** — cumulative where-stranded density for
   strandings at age ≤ T (the beaching analogue of `026`'s horizons).
4. **Cumulative beached fraction vs. age** — the stranding time course,
   with the 2016–2019 interannual min–max as a band around the pooled
   curve. Maps stay pooled; only the scalar-per-horizon curve carries the
   spread.
5. **Travel distance at stranding** — stranded weight over `disp_bin`, the
   crow-flies displacement from release at the deposit step. Short-travel
   strandings are a real outcome under a threshold rate model, never a
   mask.

Figures are print-ready (full page width, 300 dpi) and the per-horizon
where-stranded field is exported as GeoJSON.

```python
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from shapely.geometry import box
from cartopy.io.shapereader import natural_earth
```

# Parameters

```python tags=["parameters"]
# Read root of the hex-aggregate store.
output_root = "../output"
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
# (e.g. "step_wc0p1_ts3_tcinf"). Selects the partitions and tags the figures.
member = "step_wc0p1_ts3_tcinf"
# Travel-distance bin width (km) of the store's `disp_bin` axis (must match
# the reducer).
disp_bin_km = 10.0

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

A whole-year partition (no `_mMM` suffix at all, written when 024e itself
ran with `release_month = 0`) carries every month, so it is only read for
`season = "ALL"`.

```python
store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

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
```

# Rendering helpers

`hex_gdf` sums a value column over a hex column and joins to the key
geometry (dropping the `-1` sentinel); `hex_map` draws hex polygons +
coastline on a plain EPSG:4326 axis (as 025/026), with the colorbar as an
inset axes at axes-fraction height 1.0 so it is exactly the map height —
a fixed-aspect map never fills its gridspec cell, and a cell-sized
colorbar overshoots it.

```python
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
coast = gpd.read_file(
    natural_earth(resolution="10m", category="physical", name="coastline")
).clip(box(lon_min, lat_min, lon_max, lat_max))
domain_aspect = (
    (lon_max - lon_min) * np.cos(np.radians(0.5 * (lat_min + lat_max)))
) / (lat_max - lat_min)

beached = beaching[beaching["beach_hex"] >= 0]
```

# Where-stranded density

Log-scale beached weight (expected particles) per stranding hex, summed over
the store's remaining axes (`disp_bin`, and `shore_type` where the reducer
emits it). Neither is shown on the maps: `shore_type` is a threshold label
on the HELCOM BRISK + CLMS flat fraction `ff` at the stranding hour (`flat`
if `ff >= 50%` else `wall`), while `trap_flat`/`trap_wall` are reducer
parameters that only surface in the member tag when they differ from 1;
travel distance gets its own figure below.

```python
gdf_stranded = hex_gdf(beached, "beach_hex")
strand_norm = log_norm(gdf_stranded["value"])

fig, ax = plt.subplots(layout="constrained")
fig.set_size_inches(FIGURE_WIDTH_IN, grid_height_in(1, 1, domain_aspect, FIGURE_WIDTH_IN))
hex_map(gdf_stranded, ax, norm=strand_norm, title="stranded weight",
        label="stranded weight (particles)")
fig_path = figure_dir / f"WhereStranded_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Beached fraction per source hex

Of the drifters released in each source hex, the fraction that strands
within the viability window (linear 0–1). Highlights which source regions
lose most of their propagules to the coast.

```python
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
```

# Beaching age horizons

Cumulative where-stranded density for strandings occurring at age ≤ T.
`beach_age_bin` bins age at `age_bin_days`, so horizon T keeps
`beach_age_bin < T // age_bin_days` (strandings in bins fully within T).
A shared `LogNorm` across horizons makes the fill-in over time legible.

```python
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
```

# GeoJSON export of the per-horizon deposit

One file per age horizon: the cumulative stranded weight per hex with the
024a geometry, EPSG:4326 — the horizon panels above in a form
collaborators can open in R or QGIS.

```python
for h in time_horizons_days:
    g = gdfs[h].rename(columns={"value": "stranded_weight"}).drop(columns="hex_id")
    export_path = export_dir / f"029_{regime}_{season}_{member}_T{h}d.geojson"
    g.to_file(export_path, driver="GeoJSON")
    print(f"wrote {export_path} ({len(g):,} hexes)")
```

# Cumulative beached fraction vs. age

Share of all released drifters stranded by each elapsed-time horizon — the
stranding time course, summed over the store's other axes. The pooled curve
is drawn with the 2016–2019 **interannual min–max** as a shaded band: the
per-year curves are the same reduction restricted to one `release_year`, so
the band is the honest spread behind the single pooled number.

```python
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
```

# Travel distance at stranding

Stranded weight over `disp_bin` — the crow-flies displacement from the
release point at the deposit step, in `disp_bin_km` bins — pooled over every
source hex. The residual (never-beached) rows carry `disp_bin = -1` and drop
out with the `beach_hex = -1` filter. Short-travel strandings are kept: a
strong-wave event stranding freshly released material near home is a real
outcome of the threshold rate model, not an artefact to mask.

```python
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
```

# Validation / summary

```python
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
```
