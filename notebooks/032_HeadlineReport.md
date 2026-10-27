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

# Headline report

The study's headline numbers and figure selection, regenerated from the
stores rather than transcribed. Reads the `024e`/`024g` parquet partitions
(and the filenames of every `024x` store) under `output_root/HexAggregates/`,
the `024a` key, and the GeoJSON/CSV products the `025`–`031` consumers wrote
under
`output_root/Exports/`; writes `output_root/Headline/HEADLINE.md` plus a
numbered, downscaled copy of the selected figures.

Everything that varies with the data is **discovered**: release years and
rate-model member tags come from the partition filenames, figures from a
glob per notebook. The coverage table lists every year on disk; tables A–C
pool only the `release_years` range, D–G read consumer exports that apply
the same parameter.
Adding release years or sweep members changes the report's content, never
this notebook.

Seven tables:

| table | question | source |
|---|---|---|
| A | beached fraction by release season and year, production member | `HexAgg_beaching_*` |
| B | how far that moves across the rate-model sweep | `HexAgg_beaching_*`, every member on disk |
| C | what survival weighting does to subbasin connectivity | `HexAgg_survconn_*` + `Exports/028/*.csv` |
| D | how far propagules travel, per season and regime | `Exports/027/*.geojson` |
| E | how concentrated the cloud is at the long horizon | `Exports/026/*.geojson` + `024a` key |
| F | the water area holding a share of the cloud, per horizon | `Exports/026/026_areas_*.csv` |
| G | how much of each origin's release is still afloat, per horizon | `Exports/029/afloat_by_origin_*.csv` |

The window every number refers to is read from the `024e`/`024g`
partitions' build attrs, and every pooled partition must agree on it.

A consumer export or figure older than the newest store partition it was
reduced from is **stale** and is not used: it lands in the report's Gaps
section, as a missing one does.

No figure is drawn here — the manuscript figures are the `025`–`036` PNGs,
and this stage only selects, numbers and downscales them.

```python
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image
```

# Parameters

```python tags=["parameters"]
# Read/write root of the study's cluster outputs.
output_root = "../output"
# Hex radius of the stores to read (024a/024/024b/024e/024g). Must match disk.
hex_radius = 6000
# Regimes the report covers, in reading order: the first is the primary one
# (every beaching/survival table and most figures), the second the comparison.
regimes_csv = "surface_stokes,bottom"
# The production rate-model member: tables A, C and the beaching/survival
# figures report this one, table B sweeps every member on disk.
member = "step_wc0p1_ts3_tcinf_tf1_tw0"
# Release seasons, in reading order. Contract in docs/visualisations.md.
seasons_csv = "DJF,MAM,JJA,SON,ALL"
# Release years the tables pool: "first-last" (inclusive) or one year
# ("2019"). The coverage table still lists every year on disk (2025 is
# partial); tables A–C read only these years' partitions, D–G the consumers'
# exports (pooled over the consumers' own release_years).
release_years = "2016-2024"
# Elapsed-time horizon (d) for the concentration table E: which 026 GeoJSON
# export (`_T<horizon_days>d`) it reads.
horizon_days = 180
# Cumulative elapsed-time horizons (d) for table C's Σw/Σn ratio; each ≤ the
# 024g window.
connectivity_horizons_csv = "30,60,90,180"
# Horizon (d) of table C's top-target comparison and of the connectivity
# figures: which 028 exports (`_T<top_horizon_days>d`) they read.
top_horizon_days = 60
# Hex mean depth (m) below which a hex counts as within the Fucus viable band.
viable_depth_m = 3.0
# Share q of the cumulative particle-time whose holding area table F reports;
# the 026 areas CSV rows at other q are filtered out.
area_quantile = 0.9
# Horizons (d) table F reports, the 026 map horizons; the CSV's other rows
# (026 exports every age-bin edge) are 035's time series.
area_horizons_csv = "30,60,90,180"
# Number of largest emitters per season in table C's top-target comparison.
n_top_origins = 4
# Release season the map figures are shown for.
figure_season = "SON"
# Release season the year-pooled curve/sweep figures are shown for.
pooled_season = "ALL"
# Width (px) the selected figures are downscaled to for the report.
figure_width_px = 1400
```

# Parse parameters

The season → month mapping is a fixed convention, not a knob, so it lives
here rather than in the parameters cell (papermill injects only primitives).

```python
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}

output_root = Path(output_root)
regimes = [r.strip() for r in regimes_csv.split(",") if r.strip()]
seasons = [s.strip() for s in seasons_csv.split(",") if s.strip()]
unknown = sorted(set(seasons + [figure_season, pooled_season]) - set(SEASON_MONTHS))
if unknown:
    raise ValueError(f"unknown season(s) {unknown}; expected {sorted(SEASON_MONTHS)}")
if not regimes:
    raise ValueError("regimes_csv is empty")
connectivity_horizons = [int(h) for h in connectivity_horizons_csv.split(",") if h.strip()]
if not connectivity_horizons:
    raise ValueError("connectivity_horizons_csv is empty")
area_horizons = [int(h) for h in area_horizons_csv.split(",") if h.strip()]
if not area_horizons:
    raise ValueError("area_horizons_csv is empty")
# `release_years` is "first-last" (inclusive) or a single year.
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"
analysis_span = f"{analysis_years[0]}–{analysis_years[-1]}"
print(f"release years in the analysis: {analysis_span}")

primary_regime = regimes[0]
other_regime = regimes[1] if len(regimes) > 1 else None
# `pooled_season` (ALL) is the union of the others, so a spread "across
# seasons" is taken over the rest.
map_seasons = [s for s in seasons if s != pooled_season]

store_root = output_root / "HexAggregates"
export_root = output_root / "Exports"
figure_root = output_root / "Figures"
headline_dir = output_root / "Headline"
headline_dir.mkdir(parents=True, exist_ok=True)

# Artefacts the report expected and did not find, or found stale. Written
# into HEADLINE.md as a "gaps" section instead of raising: a partially-run
# consumer matrix should still produce a report that says what is missing.
gaps = []
```

```python
def markdown_table(frame):
    """Render a DataFrame of already-formatted cells as a GitHub markdown table."""
    header = "| " + " | ".join(str(c) for c in frame.columns) + " |"
    rule = "|" + "---|" * len(frame.columns)
    body = ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False)]
    return "\n".join([header, rule, *body])


def fmt(value, decimals=1):
    """Format a float for a report cell, with an em dash for missing values."""
    if value is None or not np.isfinite(value):
        return "—"
    return f"{value:.{decimals}f}"
```

# Discover what is on disk

Release years, rate-model members and regimes are parsed out of the
partition filenames. Nothing below hardcodes a year or a member tag, so a
store that grows from four release years to ten, or from one member to a
full sweep, changes the report without changing the notebook. The full
catalogue feeds the coverage table; `analysis_parts`, its rows inside
`release_years`, feeds everything else.

The per-`(regime, year)` stores (`024` counts, `024b` distance, `024c`
connectivity) and the per-`(regime, year, month, member)` ones (`024e`
beaching, `024g` survival connectivity) are catalogued separately, because
regimes differ in how far downstream they were taken.

```python
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)[["hex_id", "mean_depth_m"]]

year_store_pattern = re.compile(
    rf"^HexAgg_(counts|distance|connectivity)_r{hex_radius}m_(.+)_(\d{{4}})\.parquet$"
)
year_parts = pd.DataFrame(
    [
        {
            "kind": m.group(1),
            "regime": m.group(2),
            "year": int(m.group(3)),
            "path": store_root / m.group(0),
        }
        for m in map(year_store_pattern.match, sorted(p.name for p in store_root.glob("*.parquet")))
        if m
    ]
)

member_store_pattern = re.compile(
    rf"^HexAgg_(beaching|survconn|survocc)_r{hex_radius}m_(.+)_(\d{{4}})_m(\d{{2}})_(.+)\.parquet$"
)
member_parts = pd.DataFrame(
    [
        {
            "kind": m.group(1),
            "regime": m.group(2),
            "year": int(m.group(3)),
            "month": int(m.group(4)),
            "member": m.group(5),
            "path": store_root / m.group(0),
        }
        for m in map(
            member_store_pattern.match, sorted(p.name for p in store_root.glob("*.parquet"))
        )
        if m
    ]
)
analysis_parts = member_parts[member_parts.year.isin(analysis_years)]

sweep_members = sorted(
    analysis_parts.loc[analysis_parts.kind == "beaching", "member"].unique().tolist()
)
beaching_regimes = [
    r
    for r in regimes
    if not analysis_parts[(analysis_parts.kind == "beaching") & (analysis_parts.regime == r)].empty
]
if member not in sweep_members:
    raise ValueError(f"member {member!r} has no 024e partitions; on disk: {sweep_members}")
if not beaching_regimes:
    raise ValueError(f"no 024e beaching partitions for any of {regimes}")

print(f"key {key_path.name}: {len(key)} hexes")
print(f"per-(regime, year) partitions: {len(year_parts)}")
print(f"per-(regime, year, month, member) partitions: {len(member_parts)}, "
      f"{len(analysis_parts)} in {analysis_span}")
print(f"024e members on disk ({len(sweep_members)}): {', '.join(sweep_members)}")
print(f"beaching regimes: {', '.join(beaching_regimes)}")

# The 025/026 family runs survival-weighted at the production member where a
# 024f store exists for the regime (only surface_stokes has 024d sidecars)
# and unweighted otherwise — the rule scripts/025-031_ConsumersAll_job.sh
# applies. Its exports and figures carry that tag in their names.
heatmap_tag = {
    r: member
    if not analysis_parts[
        (analysis_parts.kind == "survocc")
        & (analysis_parts.regime == r)
        & (analysis_parts.member == member)
    ].empty
    else "unweighted"
    for r in regimes
}
print(f"025/026 weighting per regime: {heatmap_tag}")
```

# Staleness

Every export and figure this report reads is a reduction of store
partitions, and a store can be rebuilt after the consumer ran (a longer
window, a new key). Such a product still exists under its old name and
would be captioned with the new window, so presence is not enough: a
product is **stale** when its mtime is older than the newest
`release_years` partition of the store it derives from, or than the `024a`
key, and a stale product is a gap like a missing one. The derivation is a
list of `(store kind, regime, member)` sources per product (`member = None`:
every member, for the sweep); the 025/026 family derives from 024f at the
production member where the regime is weighted, from the 024 counts
otherwise. 035 draws 026's area CSVs rather than a store, so its sources
are those CSVs, each with its own store sources: a 035 figure drawn from a
stale CSV is stale however new it is. Considered recording the input
partitions' build ids in every product and comparing those; rejected,
mtime needs no change to the consumers. Copying outputs without preserving
mtimes defeats it.

```python
store_parts = pd.concat([year_parts, member_parts], ignore_index=True)
store_parts = store_parts[store_parts.year.isin(analysis_years)].assign(
    mtime=lambda d: [p.stat().st_mtime for p in d.path]
).reset_index(drop=True)


def stamp(mtime):
    """An mtime as a UTC timestamp for a gap line."""
    return f"{datetime.fromtimestamp(mtime, timezone.utc):%Y-%m-%d %H:%M UTC}"


def heatmap_sources(regime):
    """The store a regime's 025/026 products (and 035's line) reduce."""
    if regime is None:
        return []
    if heatmap_tag[regime] == "unweighted":
        return [("counts", regime, None)]
    return [("survocc", regime, heatmap_tag[regime])]


def newest_input(sources):
    """(mtime, file name) of the newest input among `sources`, the key
    included. A source is a store `(kind, regime, member)` or an upstream
    product `(path, its sources)`; a missing or stale upstream product counts
    as newer than anything (mtime inf)."""
    newest = (key_path.stat().st_mtime, key_path.name)
    for source in sources:
        if isinstance(source[0], Path):
            path, upstream = source
            mtime = path.stat().st_mtime if path.exists() else np.inf
            if mtime < newest_input(upstream)[0]:
                mtime = np.inf
            newest = max(newest, (mtime, path.name))
            continue
        kind, regime, member_tag = source
        sel = store_parts[(store_parts.kind == kind) & (store_parts.regime == regime)]
        if member_tag is not None:
            sel = sel[sel.member == member_tag]
        if not sel.empty:
            row = sel.loc[sel.mtime.idxmax()]
            newest = max(newest, (row.mtime, row.path.name))
    return newest


def usable(path, sources, what):
    """Whether a product exists and is not older than its inputs; records a
    gap otherwise."""
    if not path.exists():
        gaps.append(f"{what}: missing `{path.name}`")
        return False
    mtime, name = newest_input(sources)
    if mtime == np.inf:
        gaps.append(f"{what}: its input `{name}` is missing or stale — rerun both")
        return False
    if path.stat().st_mtime < mtime:
        gaps.append(
            f"{what}: `{path.name}` ({stamp(path.stat().st_mtime)}) is older than its input "
            f"`{name}` ({stamp(mtime)}) — stale, rerun its producer"
        )
        return False
    return True
```


# Coverage

Which release years each store actually holds, per regime. Reported
explicitly because the regimes are not in step: `surface` was not taken
past the aggregates, and the beaching/survival stores exist only where a
`024d` sidecar was built.


```python
def year_span(years):
    """Render a year set as a compact span, marking gaps."""
    years = sorted(set(int(y) for y in years))
    if not years:
        return "—"
    contiguous = years == list(range(years[0], years[-1] + 1))
    span = f"{years[0]}–{years[-1]}" if len(years) > 1 else f"{years[0]}"
    return f"{span} ({len(years)})" if contiguous else f"{', '.join(map(str, years))} ({len(years)})"


coverage_rows = []
for regime in sorted(set(year_parts.regime) | set(member_parts.regime)):
    row = {"regime": regime}
    for kind in ("counts", "distance", "connectivity"):
        sel = year_parts[(year_parts.regime == regime) & (year_parts.kind == kind)]
        row[f"024 {kind}"] = year_span(sel.year)
    for kind, label in (
        ("beaching", f"024e {member}"),
        ("survocc", f"024f {member}"),
        ("survconn", f"024g {member}"),
    ):
        sel = member_parts[
            (member_parts.regime == regime)
            & (member_parts.kind == kind)
            & (member_parts.member == member)
        ]
        row[label] = (
            f"{year_span(sel.year)}, {len(sel)} monthly parts" if len(sel) else "—"
        )
    coverage_rows.append(row)
coverage = pd.DataFrame(coverage_rows)
coverage_md = markdown_table(coverage)
print(coverage.to_string(index=False))
```

# Tables A and B — beached fraction

One pass over every `024e` partition of every member: the beached fraction
of a group is `Σ weight[beach_hex ≥ 0] / Σ weight`, with `beach_hex = -1`
the never-beached residual. Weights are additive across months and years,
so a season pools by summing its months and a pooled multi-year value is
the ratio of the summed weights, not the mean of the yearly ratios.

Each partition's build attrs are read per file (pandas keeps no attrs
through a concat of differing frames). Within a member every partition must
agree on the window, the binning and the rate-model point; across members
the rate model differs by design, but the window and binning must not, or
table B would compare stranding over different elapsed times.

```python
# Keys every pooled partition of one member must share. Provenance (builder,
# git_sha, slurm_job_id) and the partition's identity (release_year,
# release_month) legitimately differ and are not compared.
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


beaching_rows = []
beaching_attrs = {}
for part in analysis_parts[analysis_parts.kind == "beaching"].itertuples(index=False):
    if part.regime not in beaching_regimes:
        continue
    frame = pd.read_parquet(part.path, columns=["beach_hex", "weight"])
    beaching_attrs.setdefault(part.member, []).append(frame.attrs)
    deposited = frame.beach_hex >= 0
    beaching_rows.append(
        {
            "regime": part.regime,
            "member": part.member,
            "year": part.year,
            "month": part.month,
            "beached": float(frame.weight[deposited].sum()),
            "total": float(frame.weight.sum()),
        }
    )
beaching = pd.DataFrame(beaching_rows)
member_windows = {
    m: (a["max_age_days"], a["age_bin_days"])
    for m, a in ((m, pooled_attrs(attrs, f"024e {m}")) for m, attrs in beaching_attrs.items())
}
if len(set(member_windows.values())) != 1:
    raise ValueError(f"024e members disagree on (max_age_days, age_bin_days): {member_windows}")
max_age_days, age_bin_days = member_windows[member]
print(f"read {len(beaching)} beaching partitions; every member reduced over "
      f"{max_age_days} d in {age_bin_days} d bins")
```

```python
def beached_percent(frame, season):
    """Pooled beached percentage and per-release-year series for one season."""
    sub = frame[frame.month.isin(SEASON_MONTHS[season])]
    if sub.empty:
        return np.nan, pd.Series(dtype=float)
    per_year = sub.groupby("year")[["beached", "total"]].sum()
    return (
        100.0 * sub.beached.sum() / sub.total.sum(),
        100.0 * per_year.beached / per_year.total,
    )


table_a_md = {}
table_a_pooled = {}
table_a_spread = {}
for regime in beaching_regimes:
    frame = beaching[(beaching.regime == regime) & (beaching.member == member)]
    years = sorted(frame.year.unique())
    rows = []
    table_a_pooled[regime] = {}
    table_a_spread[regime] = {}
    for season in seasons:
        pooled, per_year = beached_percent(frame, season)
        table_a_pooled[regime][season] = pooled
        table_a_spread[regime][season] = (
            per_year.max() - per_year.min() if len(per_year) else np.nan
        )
        row = {"season": season}
        row.update({str(y): fmt(per_year.get(y, np.nan)) for y in years})
        row["pooled"] = f"**{fmt(pooled)}**"
        row["interannual min–max"] = (
            f"{fmt(per_year.min())} – {fmt(per_year.max())}" if len(per_year) else "—"
        )
        rows.append(row)
    table_a_md[regime] = markdown_table(pd.DataFrame(rows))
    print(f"\ntable A — {regime}\n{table_a_md[regime]}")
```

```python
table_b_md = {}
table_b_pooled = {}
for regime in beaching_regimes:
    rows = []
    table_b_pooled[regime] = {}
    for sweep_member in sweep_members:
        frame = beaching[(beaching.regime == regime) & (beaching.member == sweep_member)]
        if frame.empty:
            continue
        tag = f"**`{sweep_member}`**" if sweep_member == member else f"`{sweep_member}`"
        row = {"member": tag}
        for season in seasons:
            row[season] = fmt(beached_percent(frame, season)[0])
        # Sort key and interannual range both come from the pooled season,
        # whether or not it is one of the reported columns.
        sort_key, per_year = beached_percent(frame, pooled_season)
        table_b_pooled[regime][sweep_member] = sort_key
        row[f"{pooled_season} min–max"] = (
            f"{fmt(per_year.min())} – {fmt(per_year.max())}" if len(per_year) else "—"
        )
        rows.append((sort_key, row))
    table_b = pd.DataFrame([r for _, r in sorted(rows, key=lambda t: -t[0])])
    table_b_md[regime] = markdown_table(table_b)
    print(f"\ntable B — {regime}\n{table_b_md[regime]}")
```

# Table C — survival-weighted vs unweighted connectivity

Two reductions of the same comparison. `Σw_obs / Σn_obs` comes from the
`024g` partitions, where both weights sit on the same rows over the same
population (real drifters only) — `Exports/028` carries one value column
per file, so the ratio cannot be formed there. The top-target comparison
reads `028`'s emission-fraction CSVs, so the table and the connectivity
figures show the same numbers.

```python
survconn_rows = []
survconn_attrs = []
for part in analysis_parts[
    (analysis_parts.kind == "survconn") & (analysis_parts.member == member)
].itertuples(index=False):
    if part.regime not in beaching_regimes:
        continue
    frame = pd.read_parquet(part.path, columns=["age_bin", "n_obs", "w_obs"])
    survconn_attrs.append(frame.attrs)
    grouped = frame.groupby("age_bin", as_index=False)[["n_obs", "w_obs"]].sum()
    grouped["regime"] = part.regime
    grouped["year"] = part.year
    grouped["month"] = part.month
    survconn_rows.append(grouped)
survconn = pd.concat(survconn_rows) if survconn_rows else pd.DataFrame()
if survconn.empty:
    gaps.append(f"no `024g` survconn partitions for member `{member}` — table C's Σw/Σn is absent.")
else:
    survconn_window = pooled_attrs(survconn_attrs, f"024g {member}")
    if (survconn_window["max_age_days"], survconn_window["age_bin_days"]) != (
        max_age_days, age_bin_days
    ):
        raise ValueError(
            f"024g window {survconn_window['max_age_days']} d / "
            f"{survconn_window['age_bin_days']} d bins != 024e {max_age_days} d / "
            f"{age_bin_days} d bins"
        )
    if max(connectivity_horizons) > max_age_days:
        raise ValueError(
            f"connectivity horizon {max(connectivity_horizons)} d is past the "
            f"{max_age_days} d window"
        )
print(f"read {len(survconn_rows)} survconn partitions")
```

```python
table_c_ratio_md = {}
table_c_ratio = {}
for regime in (r for r in beaching_regimes if not survconn.empty):
    frame = survconn[survconn.regime == regime]
    rows = []
    table_c_ratio[regime] = {}
    for season in seasons:
        row = {"season": season}
        for horizon in connectivity_horizons:
            sub = frame[
                frame.month.isin(SEASON_MONTHS[season]) & (frame.age_bin * age_bin_days < horizon)
            ]
            ratio = sub.w_obs.sum() / sub.n_obs.sum()
            table_c_ratio[regime][(season, horizon)] = ratio
            row[f"Σw/Σn at {horizon} d"] = fmt(ratio, 3)
        rows.append(row)
    table_c_ratio_md[regime] = markdown_table(pd.DataFrame(rows))
    print(f"\ntable C ratio — {regime}\n{table_c_ratio_md[regime]}")
```

```python
def top_targets(emission, origin, n=3):
    """The `n` largest targets of one origin, as 'name pct' strings."""
    sub = emission[emission.origin_subbasin == origin].nlargest(n, "emission_fraction")
    return [f"{t.target_subbasin} {100 * t.emission_fraction:.1f}" for t in sub.itertuples()]


table_c_top_md = {}
flips = []
for regime in beaching_regimes:
    export_dir = export_root / "028"
    rows = []
    for season in seasons:
        stem = f"connectivity_{regime}_{season}"
        paths = {
            "n_obs": export_dir / f"{stem}_unweighted_T{top_horizon_days}d_n_obs.csv",
            "unweighted": export_dir / f"{stem}_unweighted_T{top_horizon_days}d_emission_fraction.csv",
            "weighted": export_dir / f"{stem}_{member}_T{top_horizon_days}d_emission_fraction.csv",
        }
        sources = {
            "n_obs": [("connectivity", regime, None)],
            "unweighted": [("connectivity", regime, None)],
            "weighted": [("survconn", regime, member)],
        }
        what = f"table C top targets, {regime} {season}"
        if not all([usable(paths[k], sources[k], what) for k in paths]):
            continue
        counts = pd.read_csv(paths["n_obs"])
        unweighted = pd.read_csv(paths["unweighted"])
        # 028 adds a target "stranded" row to weighted CSVs; it is not a subbasin.
        weighted = pd.read_csv(paths["weighted"]).query("target_subbasin != 'stranded'")
        by_origin = counts.groupby("origin_subbasin")["n_obs"].sum()
        share = 100.0 * by_origin / by_origin.sum()
        for origin in share.nlargest(n_top_origins).index:
            unweighted_top = top_targets(unweighted, origin)
            weighted_top = top_targets(weighted, origin)
            rows.append(
                {
                    "season": season,
                    "origin (share of all)": f"{origin} ({share[origin]:.1f} %)",
                    "unweighted top-3": ", ".join(unweighted_top),
                    "survival-weighted top-3": ", ".join(weighted_top),
                }
            )
            if [t.rsplit(" ", 1)[0] for t in unweighted_top] != [
                t.rsplit(" ", 1)[0] for t in weighted_top
            ]:
                flips.append(f"{regime} {season} {origin}")
    if rows:
        table_c_top_md[regime] = markdown_table(pd.DataFrame(rows))
        print(f"\ntable C top targets — {regime}: {len(rows)} rows, {len(flips)} ranking flips")
```

# Table D — dispersal distance

`027` exports one per-source-hex quantile set per `(regime, season)`; the
table reports the **median over source hexes** of each quantile level, so
one number per (season, regime, level) summarises a whole map. The quantile
levels are read off the GeoJSON columns rather than assumed.

```python
distance_rows = []
for season in seasons:
    for regime in regimes:
        path = export_root / "027" / f"027_{regime}_{season}.geojson"
        if not usable(path, [("distance", regime, None)], "table D"):
            continue
        frame = gpd.read_file(path)
        for column in [c for c in frame.columns if re.fullmatch(r"q[\d.]+_km", c)]:
            distance_rows.append(
                {
                    "season": season,
                    "regime": regime,
                    "level": column.removesuffix("_km"),
                    "median_km": float(frame[column].median()),
                    "source_hexes": len(frame),
                }
            )
# Explicit columns: every 027 export can be missing or stale, and the table
# then renders as an empty row per season rather than raising.
distance = pd.DataFrame(
    distance_rows, columns=["season", "regime", "level", "median_km", "source_hexes"]
)
table_d = (
    distance.assign(
        column=distance.regime + " " + distance.level,
        cell=distance.median_km.round().astype(int).astype(str) + " km",
    )
    .pivot(index="season", columns="column", values="cell")
    .reindex(seasons)
    .reset_index()
    .fillna("\u2014")
)
table_d_md = markdown_table(table_d)
print(distance.to_string(index=False))
```

# Table E — concentration at the long horizon

`026` exports the cumulative relative density per hex at each horizon
(`percentage` = share of that horizon's particle-time, survival-weighted
where the regime has a 024f store — `heatmap_tag` above). Joined to the
`024a` key for `mean_depth_m`, it gives the concentration of the cloud and
the crude share of particle-time inside the Fucus viable band. A hex-mean
depth cannot resolve a 3 m contour, so the last column is a contrast
between regimes, not a rate.

```python
concentration_rows = []
for regime in regimes:
    for season in seasons:
        path = (
            export_root / "026"
            / f"026_{regime}_{season}_{heatmap_tag[regime]}_T{horizon_days}d.geojson"
        )
        if not usable(path, heatmap_sources(regime), "table E"):
            continue
        frame = gpd.read_file(path).drop(columns="geometry").merge(key, on="hex_id", how="left")
        ranked = frame.percentage.sort_values(ascending=False)
        concentration_rows.append(
            {
                "regime": regime,
                "weighting": heatmap_tag[regime],
                "season": season,
                "max_hex_pct": float(ranked.iloc[0]),
                "top10_pct": float(ranked.iloc[:10].sum()),
                "top100_pct": float(ranked.iloc[:100].sum()),
                "occupied_hexes": len(frame),
                "viable_pct": float(frame.percentage[frame.mean_depth_m <= viable_depth_m].sum()),
            }
        )
concentration = pd.DataFrame(
    concentration_rows,
    columns=[
        "regime", "weighting", "season", "max_hex_pct", "top10_pct", "top100_pct",
        "occupied_hexes", "viable_pct",
    ],
)
table_e = pd.DataFrame(
    {
        "regime": concentration.regime,
        "weighting": concentration.weighting,
        "season": concentration.season,
        "max hex %": concentration.max_hex_pct.map(lambda v: fmt(v, 2)),
        "top-10 %": concentration.top10_pct.map(lambda v: fmt(v, 2)),
        "top-100 %": concentration.top100_pct.map(fmt),
        "occupied hexes": concentration.occupied_hexes.astype(str),
        f"share in hexes \u2264 {viable_depth_m:g} m": concentration.viable_pct.map(fmt) + " %",
    }
)
table_e_md = markdown_table(table_e)
print(concentration.to_string(index=False))
```

# Table F — area holding a share of the cloud

`026` writes, per `(regime, season)`, the smallest water area holding a
share `q` of the cumulative particle-time by each horizon, pooled over the
release years (`scope = pooled`) and per release year. The CSV is read
generically: `scope` splits pooled from per-year rows, every float column
other than `q` is a reported value (`area_km2`, `length_km`), and every
remaining column but `q` is a row key — so a key column the export gains
later (an origin subbasin, say) becomes a row of this table without an edit
here. Only the rows at `area_quantile` and at the `area_horizons` (the map
horizons) are reported — 026 exports every age-bin edge for 035's time
series, and further shares; they are filtered, not assumed absent, and a
requested horizon the CSV lacks is a gap. A cell is the pooled value with
the per-year min–max in brackets. The pooled area ranks the years' summed
weights, which spread over more hexes than any one year's, so it can lie
outside the per-year range.

```python
area_rows = []
for regime in regimes:
    for season in seasons:
        path = export_root / "026" / f"026_areas_{regime}_{season}_{heatmap_tag[regime]}.csv"
        if not usable(path, heatmap_sources(regime), "table F"):
            continue
        frame = pd.read_csv(path, dtype={"scope": str})
        frame = frame[np.isclose(frame.q, area_quantile) & frame.horizon_days.isin(area_horizons)]
        absent = sorted(set(area_horizons) - set(frame.horizon_days))
        if absent:
            gaps.append(
                f"table F: `{path.name}` has no q = {area_quantile:g} rows at "
                f"{'/'.join(map(str, absent))} d"
            )
        area_rows.append(frame.assign(regime=regime, weighting=heatmap_tag[regime], season=season))
areas = pd.concat(area_rows, ignore_index=True) if area_rows else pd.DataFrame()

table_f_md = {}
if not areas.empty:
    value_columns = [
        c for c in areas.columns
        if c != "q" and pd.api.types.is_float_dtype(areas[c])
    ]
    run_columns = ["regime", "weighting", "season"]
    key_columns = run_columns + [
        c for c in areas.columns
        if c not in value_columns and c not in ("q", "scope", "horizon_days", *run_columns)
    ]
    pooled_rows = areas[areas.scope == "pooled"].set_index(key_columns + ["horizon_days"])
    year_rows = areas[areas.scope != "pooled"].groupby(key_columns + ["horizon_days"])
    for value in value_columns:
        cells = pd.DataFrame(
            {"pooled": pooled_rows[value], "lo": year_rows[value].min(), "hi": year_rows[value].max()}
        )
        cells["cell"] = [
            f"{fmt(p, 0)} [{fmt(lo, 0)}–{fmt(hi, 0)}]" for p, lo, hi in cells.itertuples(index=False)
        ]
        wide = cells["cell"].unstack("horizon_days")
        wide.columns = [f"{h} d" for h in wide.columns]
        wide = wide.reset_index()
        # Regimes and seasons in reading order, not alphabetically.
        reading_order = {"regime": regimes.index, "season": seasons.index}
        wide = wide.sort_values(
            ["regime", "season"], key=lambda s: s.map(reading_order[s.name]), kind="stable"
        )
        table_f_md[value] = markdown_table(wide)
        print(f"\ntable F — {value}\n{table_f_md[value]}")
```

# Table G — afloat by origin subbasin

`029` writes, per `(regime, season, member)`, the percentage of each HELCOM
origin subbasin's released water-seeded drifters not yet stranded by each
horizon, pooled over the release years and with the per-year min–max.
"Afloat" includes drifters that left the BSH domain, and no lifetime or
sinking is applied (the CSV's header says so too). The table shows the
year-pooled season per origin; the basin-wide `all` row per season follows
it.

```python
afloat_rows = []
for regime in beaching_regimes:
    for season in seasons:
        path = export_root / "029" / f"afloat_by_origin_{regime}_{season}_{member}.csv"
        if not usable(path, [("beaching", regime, member)], "table G"):
            continue
        afloat_rows.append(
            pd.read_csv(path, comment="#", dtype={"origin_subbasin": str}).assign(
                regime=regime, season=season
            )
        )
afloat = pd.concat(afloat_rows, ignore_index=True) if afloat_rows else pd.DataFrame()


def afloat_wide(frame, index):
    """`index` × horizon of 'pooled [min–max]' cells."""
    cells = frame.assign(
        cell=[
            f"{fmt(p)} [{fmt(lo)}–{fmt(hi)}]"
            for p, lo, hi in frame[["afloat_pct", "afloat_pct_min", "afloat_pct_max"]].itertuples(index=False)
        ]
    )
    wide = cells.pivot(index=index, columns="horizon_days", values="cell")
    wide.columns = [f"{h} d" for h in wide.columns]
    return wide


table_g_md = {}
table_g_all_md = {}
afloat_horizons = []
for regime in (r for r in beaching_regimes if not afloat.empty):
    frame = afloat[afloat.regime == regime]
    if frame.empty:
        continue
    afloat_horizons = sorted(frame.horizon_days.unique())
    pooled_frame = frame[frame.season == pooled_season]
    if not pooled_frame.empty:
        first = pooled_frame[pooled_frame.horizon_days == afloat_horizons[0]]
        named = first[~first.origin_subbasin.isin(["-1", "all"])].sort_values(
            "released_weight", ascending=False
        )
        order = (
            named.origin_name.tolist()
            + first[first.origin_subbasin == "-1"].origin_name.tolist()
            + first[first.origin_subbasin == "all"].origin_name.tolist()
        )
        share = (
            100 * pooled_frame.set_index(["origin_name", "horizon_days"]).released_weight
            .xs(afloat_horizons[0], level="horizon_days")
            / first[first.origin_subbasin == "all"].released_weight.sum()
        )
        wide = afloat_wide(pooled_frame, "origin_name").reindex(order)
        wide.insert(0, "share of released %", share.reindex(order).map(fmt))
        table_g_md[regime] = markdown_table(wide.reset_index().rename(columns={"origin_name": "origin"}))
    basin = frame[frame.origin_subbasin == "all"]
    table_g_all_md[regime] = markdown_table(
        afloat_wide(basin, "season").reindex([s for s in seasons if s in set(basin.season)]).reset_index()
    )
    print(f"\ntable G — {regime}\n{table_g_md.get(regime, '')}\n\n{table_g_all_md[regime]}")
```

# Figures

The selection is a **rule**, not a file list: each entry names the
producing notebook and a filename *stem* built from the run's own
parameters (regime, season, hex radius, member or the regime's 025/026
weighting tag, horizon). Resolution takes the exact `<stem>.png` only: every
field of a stem is a run parameter, and a looser `<stem>_*.png` match would
silently resolve one member's figure for another whose tag it extends
(`step_wc0p1_ts3_tcinf` → `step_wc0p1_ts3_tcinf_tf1_tw0`). An entry that
resolves to nothing, or to a stale PNG (its `sources`, see Staleness), is
recorded as a gap rather than raising: a partially-run consumer matrix
should still yield a report. An entry read against another (`pairs_with`)
names that figure's number, or says it is absent.

Captions take their horizons and windows from the parameters, the store
attrs and the 026 exports on disk — nothing is hardcoded.

The resolved PNGs are copied into `Headline/` numbered in reading order and
downscaled to `figure_width_px`; the 180 mm / 300 dpi originals stay under
`Figures/` as the manuscript versions.

```python
figure_fields = {
    "regime": primary_regime,
    "other_regime": other_regime,
    "season": figure_season,
    "pooled": pooled_season,
    "radius": hex_radius,
    "member": member,
    "tag": heatmap_tag[primary_regime],
    "other_tag": heatmap_tag[other_regime] if other_regime else None,
    "horizon": top_horizon_days,
}


def weighting_phrase(tag):
    """How a 025/026 figure stem's tag weights the particle-time."""
    return "unweighted" if tag == "unweighted" else f"survival-weighted (`{tag}`)"


def horizons_on_disk(regime, season):
    """The 026 horizons exported for one run, from the GeoJSON names."""
    if regime is None:
        return "—"
    pattern = re.compile(rf"026_{regime}_{season}_{heatmap_tag[regime]}_T(\d+)d\.geojson")
    matches = (pattern.fullmatch(p.name) for p in (export_root / "026").glob("026_*.geojson"))
    found = sorted(int(m.group(1)) for m in matches if m)
    return "/".join(map(str, found)) if found else "—"


figure_specs = [
    dict(
        # The release set and its subbasins come from the data twin and the
        # 024a key, so the key is the only staleness source.
        notebook="036", stem="SeedingMap_r{radius}m", name="seeding_map",
        sources=[],
        note="The REDLIST release cells by HELCOM subbasin (the 024a key's shoreward rule, "
        "cell counts in the legend), and the share of each cell's seeded particles that start "
        "on BSH land and never move.",
    ),
    dict(
        notebook="025", stem="HexRelativeDensity_{regime}_{pooled}_r{radius}m_{tag}",
        name="density_baltic_{regime}_{pooled}",
        sources=heatmap_sources(primary_regime),
        note=f"Basin-wide relative density, {weighting_phrase(figure_fields['tag'])}, ages < "
        f"{max_age_days} d: share of particle-time (left) and dilution per km² of hex water area "
        "(right). Where the particle-time sits, and how much of it is a coastal-hex artefact.",
    ),
    dict(
        notebook="025", stem="HexRelativeDensityDE_{regime}_{pooled}_r{radius}m_{tag}",
        name="density_german_waters_{regime}_{pooled}",
        sources=heatmap_sources(primary_regime),
        note="The same two fields on the German-waters viewport, normalisation still basin-wide.",
    ),
    dict(
        notebook="026", stem="TimeHorizonMaps_{regime}_{season}_r{radius}m_{tag}_percentage",
        name="horizons_{regime}_{season}",
        sources=heatmap_sources(primary_regime),
        note=f"Cumulative {horizons_on_disk(primary_regime, figure_season)} d horizons, "
        f"{weighting_phrase(figure_fields['tag'])}, each normalised to its own 100 %: how the "
        "cloud fills the basin.",
    ),
    dict(
        notebook="026", stem="TimeHorizonMaps_{other_regime}_{season}_r{radius}m_{other_tag}_percentage",
        name="horizons_{other_regime}_{season}",
        sources=heatmap_sources(other_regime),
        note=f"The same for the comparison regime ({horizons_on_disk(other_regime, figure_season)} d, "
        f"{weighting_phrase(figure_fields['other_tag']) if other_regime else '—'}).",
        pairs_with="horizons_{regime}_{season}",
    ),
    dict(
        notebook="026", stem="MedianAge_{regime}_{pooled}_r{radius}m_{tag}",
        name="median_age_{regime}_{pooled}",
        sources=heatmap_sources(primary_regime),
        note=f"Median particle age per hex over ages < {max_age_days} d, "
        f"{weighting_phrase(figure_fields['tag'])}: the age at which the hex's cumulative "
        "weight over the age bins reaches half its total (linear within the crossing bin), "
        "release years pooled.",
    ),
    dict(
        notebook="026", stem="MedianAge_{other_regime}_{pooled}_r{radius}m_{other_tag}",
        name="median_age_{other_regime}_{pooled}",
        sources=heatmap_sources(other_regime),
        note="The same for the comparison regime "
        f"({weighting_phrase(figure_fields['other_tag']) if other_regime else '—'}).",
        pairs_with="median_age_{regime}_{pooled}",
    ),
    dict(
        notebook="035", stem="AreaTimeSeries_{pooled}_r{radius}m", name="area_time_series_{pooled}",
        # 035 draws 026's area CSVs, so it is only as fresh as they are.
        sources=[
            (
                export_root / "026" / f"026_areas_{r}_{pooled_season}_{heatmap_tag[r]}.csv",
                heatmap_sources(r),
            )
            for r in regimes
        ],
        note="Table F's area (left) and its square root (right) against the cumulative horizon T "
        "(ages in [0, T)) at every age-bin edge, one line per regime at its 025/026 weighting: "
        "pooled line, per-release-year min–max band.",
    ),
    dict(
        notebook="029", stem="WhereStranded_{regime}_r{radius}m_{season}_{member}",
        name="where_stranded_{season}",
        sources=[("beaching", primary_regime, member)],
        note=f"Stranded within {max_age_days} d, % of the released drifters per hex (log, 4 "
        "decades): the stranding coastline of this season.",
    ),
    dict(
        notebook="029", stem="BeachedFractionCurve_{regime}_r{radius}m_{pooled}_{member}",
        name="beached_fraction_curve_{pooled}",
        sources=[("beaching", primary_regime, member)],
        note=f"Cumulative stranded % of the released drifters against age over the {max_age_days} d "
        "window, with the interannual min–max band; the curve behind table A.",
    ),
    dict(
        notebook="031", stem="BeachingSweepStats_{regime}_r{radius}m_{pooled}",
        name="beaching_sweep_stats_{pooled}",
        sources=[("beaching", primary_regime, None)],
        note=f"Stranded % by {max_age_days} d, Gini, stranding-hex count, median stranding age and "
        "median displacement across every member with interannual bars; the visual form of table B.",
    ),
    dict(
        notebook="028", stem="SubbasinConnectivity_emission_{regime}_{season}_unweighted_T{horizon}d",
        name="connectivity_{season}_T{horizon}d_unweighted",
        sources=[("connectivity", primary_regime, None)],
        note=f"Subbasin emission fraction at ages < {top_horizon_days} d, unweighted residence. "
        "Diagonal = self-retention.",
    ),
    dict(
        notebook="028", stem="SubbasinConnectivity_emission_{regime}_{season}_{member}_T{horizon}d",
        name="connectivity_{season}_T{horizon}d_survival_weighted",
        sources=[("survconn", primary_regime, member)],
        note="The same matrix survival-weighted (table C's top-target columns are the same numbers).",
        pairs_with="connectivity_{season}_T{horizon}d_unweighted",
    ),
    dict(
        notebook="027", stem="DistanceQuantiles_{regime}_r{radius}m_{season}",
        name="distance_quantiles_{regime}_{season}",
        sources=[("distance", primary_regime, None)],
        note="Per-source-hex quantiles of final crow-flies displacement; the map behind table D.",
    ),
    dict(
        notebook="030", stem="SurvivalHeatmaps_{regime}_r{radius}m_{season}_{member}",
        name="occupancy_surviving_fraction_{season}",
        sources=[("survocc", primary_regime, member)],
        note=f"Rows = horizons, each the {age_bin_days} d snapshot bin ending there; left column "
        "occupancy (% of released drifters in the hex, log bins), right column the surviving "
        "(un-beached) fraction surv/occ on a fixed 0–1 scale: where propagules still drift after "
        "the coast has taken its share.",
    ),
]
```


```python
# Stale numbering from a previous run would otherwise survive a shrinking
# figure set, since the numbers are positions, not identities.
for stale in headline_dir.glob("[0-9][0-9]_*.png"):
    stale.unlink()

figure_notes = []
# Figure number of each resolved entry, by its formatted name, for the
# `pairs_with` references.
figure_numbers = {}
for spec in figure_specs:
    notebook = spec["notebook"]
    fields = re.findall(r"{(\w+)}", spec["stem"] + spec["name"])
    if any(figure_fields.get(f) in (None, "") for f in fields):
        gaps.append(f"figure `{spec['name']}` ({notebook}): parameters give no value for one of {fields}")
        continue
    stem = spec["stem"].format(**figure_fields)
    name = spec["name"].format(**figure_fields)
    source = figure_root / notebook / f"{stem}.png"
    if not usable(source, spec["sources"], f"figure `{name}` (Figures/{notebook}/)"):
        continue
    note = spec["note"]
    if "pairs_with" in spec:
        partner = spec["pairs_with"].format(**figure_fields)
        note += (
            f" Read against figure {figure_numbers[partner]}."
            if partner in figure_numbers
            else f" Its counterpart `{partner}` is not in this report (see Gaps)."
        )
    figure_numbers[name] = len(figure_notes) + 1
    target = headline_dir / f"{figure_numbers[name]:02d}_{name}.png"
    with Image.open(source) as image:
        if image.width > figure_width_px:
            height = round(image.height * figure_width_px / image.width)
            image.resize((figure_width_px, height), Image.LANCZOS).save(target)
        else:
            shutil.copyfile(source, target)
    figure_notes.append((target.name, notebook, source.name, note))
    print(f"{target.name} <- Figures/{notebook}/{source.name}")
```

# Assemble HEADLINE.md

```python
figures_md = "\n".join(
    f"{i}. `{name}` — {note} _(`Figures/{notebook}/{source}`)_"
    for i, (name, notebook, source, note) in enumerate(figure_notes, start=1)
)

distance_recap = []
if not distance.empty and other_regime is not None:
    per_level = distance[distance.season == pooled_season].pivot(
        index="level", columns="regime", values="median_km"
    )
    if primary_regime in per_level and other_regime in per_level:
        distance_recap = [
            f"`{primary_regime}` against `{other_regime}` in the {pooled_season} pool: "
            + ", ".join(
                f"{level} {per_level.loc[level, primary_regime] / per_level.loc[level, other_regime]:.1f}\u00d7"
                for level in per_level.index
            )
            + ".",
            "",
        ]

concentration_recap = []
if not concentration.empty:
    concentration_recap = [
        "; ".join(
            f"`{r}` {g.occupied_hexes.min()}\u2013{g.occupied_hexes.max()} occupied hexes, top-100 "
            f"{g.top100_pct.min():.1f}\u2013{g.top100_pct.max():.1f} %, \u2264 {viable_depth_m:g} m share "
            f"{g.viable_pct.min():.1f}\u2013{g.viable_pct.max():.1f} %"
            for r, g in concentration.groupby("regime", sort=False)
        )
        + " (range over the reported seasons).",
        "",
    ]

headline_regime = beaching_regimes[0]
production = beaching[(beaching.regime == headline_regime) & (beaching.member == member)]
beaching_years = sorted(int(y) for y in production.year.unique())
sections = [
    "# Fucus dispersal — headline results",
    "",
    f"Regenerated by `notebooks/032_HeadlineReport.py` on "
    f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} from the stores under `{output_root}`. "
    "Every number below is recomputed from the parquet/GeoJSON products; nothing is transcribed.",
    "",
    "## Setup",
    "",
    f"Parcels runs on BSH HBMnoku currents, regimes {', '.join(f'`{r}`' for r in regimes)}; "
    f"aggregation on the {hex_radius / 1000:g} km hex grid (`024a` key, {len(key)} hexes). "
    "Release **season** is the primary temporal axis (DJF/MAM/JJA/SON, ALL = all 12 months), "
    f"release years {analysis_span} pooled. Beaching and survival use the production member "
    f"`{member}` on {', '.join(f'`{r}`' for r in beaching_regimes)}; "
    f"table B sweeps all {len(sweep_members)} members on disk. Every `024e`/`024g` partition "
    f"records a {max_age_days} d reduction window in {age_bin_days} d age bins (build attrs, "
    "asserted equal across the pooled partitions); no lifetime is applied. Beaching release "
    f"years: {year_span(beaching_years)}. The `025`/`026` maps and tables E/F are "
    + ", ".join(f"`{r}` {weighting_phrase(t)}" for r, t in heatmap_tag.items())
    + f", over the same {max_age_days} d window.",
    "",
    "## Coverage",
    "",
    f"Release years per store and regime, as found on disk. Tables A–C read only the "
    f"{analysis_span} partitions (`release_years`); tables D–G read the consumer exports, "
    "pooled over the consumers' own `release_years`.",
    "",
    coverage_md,
    "",
]

sections += [
    "## Table A — beached fraction by release season, production member",
    "",
    f"`Σ weight[beach_hex ≥ 0] / Σ weight` from "
    f"`HexAggregates/HexAgg_beaching_r{hex_radius}m_<regime>_<year>_mMM_{member}.parquet`; "
    f"`beach_hex = -1` is the never-beached residual. Stranded within the {max_age_days} d "
    "window, percent of the released water-seeded drifters.",
    "",
]
for regime, table in table_a_md.items():
    pooled = {k: v for k, v in table_a_pooled[regime].items() if k in map_seasons and np.isfinite(v)}
    spread = {k: v for k, v in table_a_spread[regime].items() if k in map_seasons and np.isfinite(v)}
    sections += [f"**{regime}**", "", table, ""]
    if pooled and spread:
        high, low = max(pooled, key=pooled.get), min(pooled, key=pooled.get)
        widest = max(spread, key=spread.get)
        sections += [
            f"Spread across release seasons: {pooled[high] - pooled[low]:.1f} points "
            f"({high} {pooled[high]:.1f} \u2212 {low} {pooled[low]:.1f}). Widest interannual spread "
            f"within a season: {spread[widest]:.1f} points ({widest}).",
            "",
        ]

sections += [
    "## Table B — sensitivity across the rate-model sweep",
    "",
    f"Same store, every `024e` member on disk, pooled over the release years; the "
    f"`{pooled_season}` column carries the interannual min–max. Percent, sorted by `{pooled_season}`.",
    "",
]
for regime, table in table_b_md.items():
    values = {k: v for k, v in table_b_pooled[regime].items() if np.isfinite(v)}
    sections += [f"**{regime}**", "", table, ""]
    if values:
        sections += [
            f"Across the {len(values)} members on disk the pooled {pooled_season} beached fraction "
            f"spans {min(values.values()):.1f}\u2013{max(values.values()):.1f} % "
            f"(`{min(values, key=values.get)}` to `{max(values, key=values.get)}`); the production "
            f"member `{member}` sits at {values[member]:.1f} %.",
            "",
        ]

sections += [
    "## Table C — survival-weighted vs unweighted connectivity",
    "",
    f"`Σw_obs / Σn_obs` from `HexAggregates/HexAgg_survconn_r{hex_radius}m_*_{member}.parquet`: "
    "the ratio needs both columns on the same rows, which only that store carries. Horizons are "
    "cumulative (ages < T), all rows kept including the `-1` unnamed-subbasin sentinel.",
    "",
]
for regime, table in table_c_ratio_md.items():
    ratios = table_c_ratio[regime]
    sections += [f"**{regime}**", "", table, ""]
    per_horizon = []
    for horizon in connectivity_horizons:
        by_season = {
            k[0]: v for k, v in ratios.items() if k[1] == horizon and k[0] in map_seasons and np.isfinite(v)
        }
        if not by_season:
            continue
        low, high = min(by_season, key=by_season.get), max(by_season, key=by_season.get)
        per_horizon.append(
            f"{horizon} d: {pooled_season} {ratios[(pooled_season, horizon)]:.3f}, "
            f"seasons {by_season[low]:.3f} ({low}) to {by_season[high]:.3f} ({high})"
        )
    if per_horizon:
        sections += ["Surviving share of the unweighted particle-time \u2014 " + "; ".join(per_horizon) + ".", ""]
sections += [
    f"Top-{n_top_origins} emitters per season at the {top_horizon_days} d horizon, from "
    f"`Exports/028/connectivity_<regime>_<season>_{{unweighted, {member}}}_T{top_horizon_days}d_"
    "emission_fraction.csv`; origins ranked by the matching `_n_obs.csv`. Target shares in percent; "
    "survival-weighted shares are fractions of the origin's unweighted particle-time "
    "(the stranded remainder is excluded, so they do not sum to 100).",
    "",
]
for regime, table in table_c_top_md.items():
    sections += [f"**{regime}**", "", table, ""]
sections += [
    f"Target ranking flips between the two weightings: {len(flips)}"
    + (f" ({'; '.join(flips)})." if flips else "."),
    "",
]

sections += [
    "## Table D — dispersal distance",
    "",
    "Median over source hexes of the per-hex quantile of final crow-flies displacement, from "
    "`Exports/027/027_<regime>_<season>.geojson`.",
    "",
    table_d_md,
    "",
    *distance_recap,
    f"## Table E — concentration at the {horizon_days} d horizon",
    "",
    f"From `Exports/026/026_<regime>_<season>_<weighting>_T{horizon_days}d.geojson` (cumulative "
    f"ages < {horizon_days} d; `percentage` = share of that horizon's particle-time, "
    "survival-weighted where the weighting column names a member), joined to "
    f"`HexAggregates/HexAgg_key_r{hex_radius}m.parquet` for `mean_depth_m`. The last column is a "
    f"crude viable-zone proxy — a {hex_radius / 1000:g} km hex mean cannot resolve the real "
    f"{viable_depth_m:g} m band, so read it as a regime contrast, not a rate.",
    "",
    table_e_md,
    "",
    *concentration_recap,
]

sections += [
    "## Table F — area holding a share of the cloud",
    "",
    "From `Exports/026/026_areas_<regime>_<season>_<weighting>.csv`: the smallest water area "
    f"holding a share q = {area_quantile:g} of the cumulative particle-time at ages < each "
    f"horizon ({'/'.join(map(str, area_horizons))} d; hexes ranked by density, the crossing hex "
    "taken fractionally); `length_km` = √area. Pooled over the release years, per-release-year "
    "min–max in brackets. The pooled area ranks the years' summed weights, which spread over "
    "more hexes than any one year's, so it can lie outside the per-year range. It ranks "
    "cumulative particle-time, so it measures concentration, not spread: where drifters "
    "accumulate in few hexes the area can shrink as the horizon grows.",
    "",
]
for value, table in table_f_md.items():
    sections += [f"**{value}**", "", table, ""]

sections += [
    "## Table G — afloat by origin subbasin",
    "",
    f"From `Exports/029/afloat_by_origin_<regime>_<season>_{member}.csv`: percent of each HELCOM "
    "origin subbasin's released water-seeded drifters not stranded at ages < each horizon "
    "(origin = the key's `helcom_subbasin` of the release hex). \"Afloat\" includes drifters "
    "that left the BSH domain, and no lifetime or sinking is applied. Pooled over the release "
    "years (summed weights), per-release-year min–max in brackets; `outside named subbasins` = "
    "origin `-1`, `all origins` = basin-wide.",
    "",
]
for regime in table_g_all_md:
    if regime in table_g_md:
        sections += [f"**{regime}, {pooled_season}**", "", table_g_md[regime], ""]
    sections += [f"**{regime}, all origins by season**", "", table_g_all_md[regime], ""]
if not afloat.empty and afloat_horizons:
    basin = afloat[(afloat.origin_subbasin == "all") & (afloat.season.isin(map_seasons))]
    last = basin[basin.horizon_days == afloat_horizons[-1]]
    if not last.empty:
        sections += [
            f"Basin-wide afloat at {afloat_horizons[-1]} d across release seasons: "
            + "; ".join(
                f"`{r}` {g.afloat_pct.min():.1f}\u2013{g.afloat_pct.max():.1f} %"
                for r, g in last.groupby("regime", sort=False)
            )
            + ".",
            "",
        ]

sections += [
    f"## Figures — {figure_width_px} px wide, downscaled from the 180 mm / 300 dpi manuscript PNGs",
    "",
    f"Maps show the `{figure_season}` release season and regime `{primary_regime}` unless the "
    f"source filename says otherwise; curves and the sweep pool the release year (`{pooled_season}`).",
    "",
    figures_md if figures_md else "_No figure resolved; see the gaps section._",
    "",
]

if gaps:
    sections += [
        "## Gaps",
        "",
        "Artefacts this report expected and did not find on disk, or found older than the "
        "newest store partition (or `024a` key) they derive from — stale, so not used:",
        "",
        *[f"- {g}" for g in gaps],
        "",
    ]

report_path = headline_dir / "HEADLINE.md"
report_path.write_text("\n".join(sections))
print(f"wrote {report_path} ({report_path.stat().st_size} bytes)")
```

# Summary

```python
pooled_all, per_year_all = beached_percent(production, pooled_season)
sweep_pooled = {
    m: beached_percent(
        beaching[(beaching.regime == headline_regime) & (beaching.member == m)], pooled_season
    )[0]
    for m in sweep_members
}
print(f"headline report: {report_path}")
print(f"regimes {regimes}; beaching regimes {beaching_regimes}; member {member}")
print(f"release years in the production beaching store: {beaching_years}")
print(
    f"beached fraction ({pooled_season}, {headline_regime}, {member}): {fmt(pooled_all)} % pooled, "
    f"{fmt(per_year_all.min())}–{fmt(per_year_all.max())} % across release years"
)
print(
    f"sweep across {len(sweep_members)} members ({pooled_season}): "
    f"{fmt(min(sweep_pooled.values()))}–{fmt(max(sweep_pooled.values()))} %"
)
print(f"tables: A {len(table_a_md)}, B {len(table_b_md)}, C {len(table_c_ratio_md)}/{len(table_c_top_md)}, "
      f"D {len(table_d)} rows, E {len(table_e)} rows, F {len(table_f_md)} value(s), "
      f"G {len(table_g_all_md)} regime(s)")
print(f"window: {max_age_days} d in {age_bin_days} d bins")
print(f"figures copied: {len(figure_notes)} of {len(figure_specs)} specified")
print(f"gaps recorded: {len(gaps)}")
for gap in gaps:
    print(f"  - {gap}")
```
