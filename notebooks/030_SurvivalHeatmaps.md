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

# Survival-weighted heatmaps: beaching folded into occupancy

Parquet-only consumer of the survival-occupancy store built by
`024f_BuildSurvivalOccupancy` (+ the `024a` key). The survival curve behind
it uses the same two-state hazard as the beaching store — background rate in
the near-shore band plus a strong-wave rate above an onshore-Stokes
threshold, its settings baked into the opaque `member` tag that selects the
partitions read here (`docs/beaching.md`). For each elapsed-time horizon `T`
it maps the **snapshot** age bin ending at `T` (ages in
`[T − age_bin_days, T)`), in three panels:

1. **Occupancy** — plain particle-residence density (`occ`), the baseline
   (no beaching), on a log scale.
2. **Survival-weighted** — the same, with beaching progressively removed
   (`surv = Σ exp(−A)`), on the *same* log scale, so the thinning is legible.

   Both are in **% of the released water-seeded drifters**: 024f keeps the
   `target_hex = -1` sentinel rows, so `Σocc` over a bin is exactly the
   released drifters `N` times the bin's hourly samples, and
   `surv / (N · age_bin_days · 24) · 100` is the mean share of released
   drifters afloat in that hex during the bin.
3. **Surviving fraction** — `surv / occ` per hex (linear 0–1), i.e. the mean
   un-beached fraction of the particles *occupying* that hex at that age.
   Read it as "how depleted is the mass here", NOT "how much beaching
   happens here": `A` is a path integral accumulated before arrival, so a
   hex reads low because the particles reaching it took wave-exposed
   near-shore routes. The sink field — where stranding actually occurs — is
   `029`'s where-stranded map, binned on `beach_hex`. Because each panel is
   one age bin, every particle in it has the same elapsed time, so the
   variation is route, not age.

A snapshot, not a cumulative `[0, T)` horizon as in 026/028/029: being
afloat at `T` is an instantaneous state, and a survival map accumulated over
`[0, T)` would be a residence map — which 026 already draws, weighted.

A fourth figure reduces the same store to the domain-wide **drifting
fraction vs. age**, with the interannual min–max as a band around the
pooled curve. Maps stay pooled; only the scalar-per-horizon curve carries
the spread.

One regime and one release `season` per run; the season's monthly
partitions are pooled across the `release_years` range. No Dask, no
zarrs — in the `025`/`026`/`029` lineage. Figures are print-ready (full page
width, 300 dpi).

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
output_root = "../output"
regime = "surface_stokes"
hex_radius = 6000
# Age-bin granularity of the store (must match 024f; asserted against the
# partitions' attrs).
age_bin_days = 10
# Release season: DJF / MAM / JJA / SON, or ALL for the pooled year. Selects
# which of the monthly `_mMM` partitions 024f wrote are read and pooled.
season = "ALL"
# Release years pooled: "first-last" (inclusive) or one year ("2019").
# Partitions of other years stay on disk unopened (2025 is partial).
release_years = "2016-2024"

# Rate-model member: the opaque tag the reducer put in the store filename
# (e.g. "step_wc0p1_ts3_tcinf_tf1_tw0"). Selects the partitions and tags
# the figures.
member = "step_wc0p1_ts3_tcinf_tf1_tw0"

# Elapsed-time horizons to map (days): the snapshot bin [T − age_bin_days, T)
# ending at each. Each a multiple of age_bin_days and ≤ the store's
# max_age_days.
time_horizons_days_csv = "30,60,90,180"

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
# `release_years` is "first-last" (inclusive) or a single year.
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"
print(f"release years in the analysis: {analysis_years[0]}–{analysis_years[-1]}")
time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
# A horizon T names the snapshot bin that *ends* at T — ages in
# [T − age_bin_days, T) — so its index is T/age_bin_days − 1. T must therefore
# land on a bin edge and be at least one bin long; the upper bound is checked
# against the store's own age axis once it is read.
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d, so "
        f"it does not land on an age-bin edge"
    )
    assert h >= age_bin_days, (
        f"horizon {h} d is shorter than one age bin ({age_bin_days} d), so no "
        f"snapshot bin ends at it"
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

figure_dir = output_root / "Figures" / "030"
figure_dir.mkdir(parents=True, exist_ok=True)
```

# Read key + pool the season's survival-occupancy partitions

Every `_mMM` partition whose month is in the season, across every release
year; the year is parsed from the filename and kept as a column, because
the interannual spread below groups by it. A whole-year partition (no
`_mMM` suffix, written when 024f itself ran with `release_month = 0`)
carries every month, so it is only read for `season = "ALL"`.

Each partition carries its build parameters as `df.attrs`. They are taken
per file **before** the concat — pandas drops attrs that differ between the
frames it concatenates — and every pooled partition must agree on the
window, the binning and the rate-model point.

```python
store_root = output_root / "HexAggregates"
key = gpd.read_parquet(store_root / f"HexAgg_key_r{hex_radius}m.parquet")

_PART_RE = re.compile(
    rf"HexAgg_survocc_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
    rf"_{re.escape(member)}\.parquet$"
)
survocc_files = []
for f in sorted(store_root.glob(f"HexAgg_survocc_r{hex_radius}m_{regime}_*.parquet")):
    m = _PART_RE.search(f.name)
    if m is None or int(m.group(1)) not in analysis_years:
        continue
    month = m.group(2)
    if month is None:
        if season != "ALL":
            continue
    elif int(month) not in season_months:
        continue
    survocc_files.append((f, int(m.group(1))))

if not survocc_files:
    raise FileNotFoundError(
        f"no survocc partitions for regime {regime!r}, season {season}, "
        f"release years {release_years}, member {member!r} at {store_root} — run 024f."
    )

# Keys every pooled partition must share: the reduction window, the binning
# and the rate-model point. Provenance (builder, git_sha, slurm_job_id) and
# the partition's identity (release_year, release_month) legitimately differ
# and are not compared.
POOLED_KEYS = (
    "max_age_days", "age_bin_days", "w_c", "tau_strong_hours", "tau_calm_days",
    "delta", "trap_flat", "trap_wall", "band_m",
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


parts = []
store_attrs = []
for f, year in survocc_files:
    df = pd.read_parquet(f).reset_index(drop=True)
    store_attrs.append(df.attrs)
    parts.append(df.assign(release_year=year))
survocc = pd.concat(parts, ignore_index=True)
del parts

window = pooled_attrs(store_attrs, f"024f {member}")
assert window["age_bin_days"] == age_bin_days, (
    f"store age_bin_days {window['age_bin_days']} != {age_bin_days}"
)
max_age_days = window["max_age_days"]
# A horizon past the window selects no bin at all.
assert max(time_horizons_days) <= max_age_days, (
    f"horizon {max(time_horizons_days)} d is past the store's reduction window "
    f"({max_age_days} d) — rebuild 024f or lower the horizon"
)
pooled_years = sorted(int(y) for y in survocc["release_year"].unique())
print(f"key: {len(key):,} hexes; pooled {len(survocc_files)} (year, month) partition(s)")
print(f"season {season} (months {season_months}); release years {pooled_years}")
print(f"024f partitions agree on {window}")
print(f"survocc rows: {len(survocc):,}")
```

# Rendering helpers

The colorbar is an inset axes at axes-fraction height 1.0, so it is exactly
the map height at any figure size — a fixed-aspect map never fills its
gridspec cell, and a cell-sized colorbar overshoots it.

```python
def hex_gdf(df, value_col):
    """A value column per target_hex as % of the released drifters, joined to
    the key geometry. `df` is one snapshot bin, every origin and the -1
    sentinel included, so its Σocc is N × the bin's hourly samples — the
    exact denominator. The key join then drops the sentinel."""
    return (
        df.groupby("target_hex")[value_col].sum()
        .div(df["occ"].sum()).mul(100).rename("value").reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def frac_gdf(df):
    """Surviving fraction surv/occ per target_hex (summed over origins),
    joined to key geometry."""
    g = df.groupby("target_hex")[["occ", "surv"]].sum()
    g["value"] = g["surv"] / g["occ"]
    return (
        g.reset_index()
        .merge(key[["hex_id", "geometry"]], left_on="target_hex", right_on="hex_id")
        .pipe(gpd.GeoDataFrame, geometry="geometry", crs="EPSG:4326")
    )


def log_norm(values):
    """Shared LogNorm floored four decades below the peak (survival weights
    have a long sub-1 tail)."""
    vmax = float(values.max())
    vmin = max(float(values[values > 0].min()), vmax / 1e4)
    return LogNorm(vmin=vmin, vmax=vmax)


def hex_map(gdf, ax, norm=None, vmin=None, vmax=None, title=None, label=None):
    if not gdf.empty:
        cax = ax.inset_axes([1.02, 0.0, 0.035, 1.0])
        gdf.plot(
            ax=ax, column="value", norm=norm, vmin=vmin, vmax=vmax,
            legend=True, cax=cax,
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
```

# Comparison grid: occupancy vs survival-weighted

Rows = horizons; columns = occupancy, survival-weighted, both in % of the
released drifters (mean share present, resp. afloat, in the hex during the
bin). Both share one `LogNorm` (across all horizons, floored four decades
below the peak) so the age-thinning and the beaching-removal are directly
readable against each other.

Surviving fraction is drawn separately below: it is a ratio in [0, 1] on a
linear scale, so putting it in this grid forces an unrelated quantity to
share a row with two log-density panels and invites reading it as a third
density map.

```python
# Snapshot bin ending at the horizon: ages in [h − age_bin_days, h).
horizon_bins = {h: h // age_bin_days - 1 for h in time_horizons_days}
bin_labels = {h: f"{h - age_bin_days}–{h} d" for h in time_horizons_days}
occ_gdfs = {h: hex_gdf(survocc[survocc["age_bin"] == b], "occ") for h, b in horizon_bins.items()}
surv_gdfs = {h: hex_gdf(survocc[survocc["age_bin"] == b], "surv") for h, b in horizon_bins.items()}
frac_gdfs = {h: frac_gdf(survocc[survocc["age_bin"] == b]) for h, b in horizon_bins.items()}

shared = pd.concat(
    [g["value"] for g in list(occ_gdfs.values()) + list(surv_gdfs.values()) if not g.empty]
)
dens_norm = log_norm(shared)

nrows = len(time_horizons_days)
fig, axes = plt.subplots(nrows=nrows, ncols=2, layout="constrained", squeeze=False)
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(nrows, 2, domain_aspect, FIGURE_WIDTH_IN)
)
for i, h in enumerate(time_horizons_days):
    hex_map(occ_gdfs[h], axes[i, 0], norm=dens_norm,
            title=f"occupancy — age {bin_labels[h]}", label="present (% of released)")
    hex_map(surv_gdfs[h], axes[i, 1], norm=dens_norm,
            title=f"survival-weighted — age {bin_labels[h]}",
            label="afloat (% of released)")
fig_path = figure_dir / f"SurvivalHeatmaps_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Surviving fraction, one panel per horizon

`surv / occ` per hex: the mean un-beached fraction of the particles
*occupying* that hex at that age — a history integral, not a local beaching
rate (see the header). Fixed linear 0–1 across panels so the age progression
is readable; the summary below prints the realised range, which for a
rare-event member sits well above 0 and makes the fixed scale look flat.

```python
ncols = len(time_horizons_days)
fig, axes = plt.subplots(nrows=1, ncols=ncols, layout="constrained", squeeze=False)
fig.set_size_inches(
    FIGURE_WIDTH_IN, grid_height_in(1, ncols, domain_aspect, FIGURE_WIDTH_IN)
)
for j, h in enumerate(time_horizons_days):
    # Panel titles carry the age span alone: the quantity is on every
    # colorbar, and "surviving fraction — age 10–20 d" is wider than a panel
    # in a three-column row, so the titles collide.
    hex_map(frac_gdfs[h], axes[0, j], vmin=0.0, vmax=1.0,
            title=f"age {bin_labels[h]}", label="surviving fraction")
fig_path = figure_dir / f"SurvivalFraction_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
for h in time_horizons_days:
    g = frac_gdfs[h]
    if not g.empty:
        print(f"  age {bin_labels[h]}: surviving fraction per hex "
              f"min {g['value'].min():.3f}, median {g['value'].median():.3f}, "
              f"max {g['value'].max():.3f}")
plt.show()
```

# Drifting fraction vs. age

The domain-wide reduction of the same store: `Σ surv / Σ occ` per age bin,
summed over every origin subbasin **and the `target_hex = -1` sentinel** —
so `Σ occ` is exactly the released drifters × the bin's hourly samples and
the curve is the mean surviving fraction of *all* released drifters,
including those that left the BSH domain (they carry the survival they had
when they left; the band never reaches them again). Dropping the sentinel
would make it the surviving fraction of the drifters still in the domain,
a population that shrinks with age. The pooled curve carries the
**interannual min–max** as a shaded band, the same reduction restricted to
one `release_year` at a time.

```python
def drifting_fraction_curve(df):
    g = df.groupby("age_bin")[["occ", "surv"]].sum().sort_index()
    return pd.Series(
        (g["surv"] / g["occ"]).to_numpy(),
        index=pd.Index(g.index.to_numpy() * age_bin_days, name="age (days)"),
        name="drifting fraction",
    )


drift_curve = drifting_fraction_curve(survocc)
drift_per_year = pd.DataFrame(
    {y: drifting_fraction_curve(g) for y, g in survocc.groupby("release_year")}
).reindex(drift_curve.index)

fig, ax = plt.subplots(layout="constrained")
# Print-ready width; a line panel needs no map aspect, so half the page width
# in height reads as a normal wide chart.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.45 * FIGURE_WIDTH_IN)
drift_curve.plot(ax=ax)
ax.fill_between(
    drift_curve.index, drift_per_year.min(axis=1), drift_per_year.max(axis=1),
    alpha=0.3,
    label=f"interannual range ({pooled_years[0]}–{pooled_years[-1]})",
)
ax.set_ylabel("drifting fraction")
ax.legend()
fig_path = figure_dir / f"DriftingFractionCurve_{regime}_r{hex_radius}m_{season}_{member}.png"
# Print-ready: fixed page width at 300 dpi (docs/visualisations.md).
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Validation / summary

```python
# Per age bin, every origin and the -1 sentinel summed: N is the released
# drifter count Σocc / (age_bin_days · 24) (exact for a full bin), and
# `left_domain` the share of the bin's samples at target_hex = -1.
per_bin = survocc.groupby("age_bin")[["occ", "surv"]].sum()
per_bin["N"] = per_bin["occ"] / (age_bin_days * 24)
per_bin["left_domain"] = (
    survocc[survocc["target_hex"] < 0].groupby("age_bin")["occ"].sum()
    .reindex(per_bin.index, fill_value=0.0) / per_bin["occ"]
)
per_bin["drifting_fraction"] = per_bin["surv"] / per_bin["occ"]
print(f"regime={regime}, hex_radius={hex_radius} m, season={season} "
      f"(months {season_months}), age_bin_days={age_bin_days}, "
      f"window {max_age_days} d, member {member}")
print(per_bin.to_string(float_format=lambda v: f"{v:,.3f}"))
for h in time_horizons_days:
    b = horizon_bins[h]
    sub = survocc[survocc["age_bin"] == b]
    o, s = sub["occ"].sum(), sub["surv"].sum()
    # The curve is indexed by each bin's start age, so the bin ending at h
    # reads at h − age_bin_days.
    per_year_bin = drift_per_year.loc[b * age_bin_days]
    print(f"  age {bin_labels[h]} (bin {b}): drifting fraction {s / max(o, 1):.3f} "
          f"[{per_year_bin.min():.3f} .. {per_year_bin.max():.3f}], "
          f"{sub.loc[sub['target_hex'] >= 0, 'target_hex'].nunique():,} hexes")
```
