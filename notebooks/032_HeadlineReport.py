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
# # Headline report
#
# The study's headline numbers and figure selection, regenerated from the
# stores rather than transcribed. Reads the `024`/`024b`/`024e`/`024g`
# parquet partitions under `output_root/HexAggregates/`, the `024a` key, and
# the GeoJSON/CSV products the `025`–`031` consumers wrote under
# `output_root/Exports/`; writes `output_root/Headline/HEADLINE.md` plus a
# numbered, downscaled copy of the selected figures.
#
# Everything that varies with the data is **discovered**: release years and
# rate-model member tags come from the partition filenames, figures from a
# glob per notebook. Adding release years or sweep members changes the
# report's content, never this notebook.
#
# Five tables:
#
# | table | question | source |
# |---|---|---|
# | A | beached fraction by release season and year, production member | `HexAgg_beaching_*` |
# | B | how far that moves across the rate-model sweep | `HexAgg_beaching_*`, every member on disk |
# | C | what survival weighting does to subbasin connectivity | `HexAgg_survconn_*` + `Exports/028/*.csv` |
# | D | how far propagules travel, per season and regime | `Exports/027/*.geojson` |
# | E | how concentrated the cloud is at the long horizon | `Exports/026/*.geojson` + `024a` key |
#
# No figure is drawn here — the manuscript figures are the `025`–`031` PNGs,
# and this stage only selects, numbers and downscales them.

# %%
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
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
# Elapsed-time horizon (d) for the concentration table E.
horizon_days = 100
# Elapsed-time horizons (d) for the connectivity table C; the first is also
# the horizon the top-target comparison and the connectivity figures use.
connectivity_horizons_csv = "50,100"
# Age bin width (d) of the counts/connectivity stores. Build↔consumer contract.
age_bin_days = 10
# Hex mean depth (m) below which a hex counts as within the Fucus viable band.
viable_depth_m = 3.0
# Number of largest emitters per season in table C's top-target comparison.
n_top_origins = 4
# Release season the map figures are shown for.
figure_season = "SON"
# Release season the year-pooled curve/sweep figures are shown for.
pooled_season = "ALL"
# Width (px) the selected figures are downscaled to for the report.
figure_width_px = 1400

# %% [markdown]
# # Parse parameters
#
# The season → month mapping is a fixed convention, not a knob, so it lives
# here rather than in the parameters cell (papermill injects only primitives).

# %%
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

# Artefacts the report expected and did not find. Written into HEADLINE.md as
# a "gaps" section instead of raising: a partially-run consumer matrix should
# still produce a report that says what is missing.
gaps = []


# %%
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


# %% [markdown]
# # Discover what is on disk
#
# Release years, rate-model members and regimes are parsed out of the
# partition filenames. Nothing below hardcodes a year or a member tag, so a
# store that grows from four release years to ten, or from one member to a
# full sweep, changes the report without changing the notebook.
#
# The per-`(regime, year)` stores (`024` counts, `024b` distance, `024c`
# connectivity) and the per-`(regime, year, month, member)` ones (`024e`
# beaching, `024g` survival connectivity) are catalogued separately, because
# regimes differ in how far downstream they were taken.

# %%
key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
key = gpd.read_parquet(key_path)[["hex_id", "mean_depth_m"]]

year_store_pattern = re.compile(
    rf"^HexAgg_(counts|distance|connectivity)_r{hex_radius}m_(.+)_(\d{{4}})\.parquet$"
)
year_parts = pd.DataFrame(
    [
        {"kind": m.group(1), "regime": m.group(2), "year": int(m.group(3))}
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

sweep_members = sorted(
    member_parts.loc[member_parts.kind == "beaching", "member"].unique().tolist()
)
beaching_regimes = [
    r
    for r in regimes
    if not member_parts[(member_parts.kind == "beaching") & (member_parts.regime == r)].empty
]
if member not in sweep_members:
    raise ValueError(f"member {member!r} has no 024e partitions; on disk: {sweep_members}")
if not beaching_regimes:
    raise ValueError(f"no 024e beaching partitions for any of {regimes}")

print(f"key {key_path.name}: {len(key)} hexes")
print(f"per-(regime, year) partitions: {len(year_parts)}")
print(f"per-(regime, year, month, member) partitions: {len(member_parts)}")
print(f"024e members on disk ({len(sweep_members)}): {', '.join(sweep_members)}")
print(f"beaching regimes: {', '.join(beaching_regimes)}")

# %% [markdown]
# # Coverage
#
# Which release years each store actually holds, per regime. Reported
# explicitly because the regimes are not in step: `surface` was not taken
# past the aggregates, and the beaching/survival stores exist only where a
# `024d` sidecar was built.


# %%
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
    for kind, label in (("beaching", f"024e {member}"), ("survconn", f"024g {member}")):
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

# %% [markdown]
# # Tables A and B — beached fraction
#
# One pass over every `024e` partition of every member: the beached fraction
# of a group is `Σ weight[beach_hex ≥ 0] / Σ weight`, with `beach_hex = -1`
# the never-beached residual. Weights are additive across months and years,
# so a season pools by summing its months and a pooled multi-year value is
# the ratio of the summed weights, not the mean of the yearly ratios.

# %%
beaching_rows = []
max_beach_age_bin = -1
for part in member_parts[member_parts.kind == "beaching"].itertuples(index=False):
    if part.regime not in beaching_regimes:
        continue
    frame = pd.read_parquet(part.path, columns=["beach_hex", "beach_age_bin", "weight"])
    deposited = frame.beach_hex >= 0
    max_beach_age_bin = max(max_beach_age_bin, int(frame.beach_age_bin.max()))
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
print(f"read {len(beaching)} beaching partitions; max beach_age_bin {max_beach_age_bin}")


# %%
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

# %%
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

# %% [markdown]
# # Table C — survival-weighted vs unweighted connectivity
#
# Two reductions of the same comparison. `Σw_obs / Σn_obs` comes from the
# `024g` partitions, where both weights sit on the same rows over the same
# population (real drifters only) — `Exports/028` carries one value column
# per file, so the ratio cannot be formed there. The top-target comparison
# reads `028`'s emission-fraction CSVs, so the table and the connectivity
# figures show the same numbers.

# %%
survconn_rows = []
for part in member_parts[
    (member_parts.kind == "survconn") & (member_parts.member == member)
].itertuples(index=False):
    if part.regime not in beaching_regimes:
        continue
    frame = pd.read_parquet(part.path, columns=["age_bin", "n_obs", "w_obs"])
    grouped = frame.groupby("age_bin", as_index=False)[["n_obs", "w_obs"]].sum()
    grouped["regime"] = part.regime
    grouped["year"] = part.year
    grouped["month"] = part.month
    survconn_rows.append(grouped)
survconn = pd.concat(survconn_rows) if survconn_rows else pd.DataFrame()
if survconn.empty:
    gaps.append(f"no `024g` survconn partitions for member `{member}` — table C's Σw/Σn is absent.")
print(f"read {len(survconn_rows)} survconn partitions")

# %%
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


# %%
def top_targets(emission, origin, n=3):
    """The `n` largest targets of one origin, as 'name pct' strings."""
    sub = emission[emission.origin_subbasin == origin].nlargest(n, "emission_fraction")
    return [f"{t.target_subbasin} {100 * t.emission_fraction:.1f}" for t in sub.itertuples()]


top_horizon = connectivity_horizons[0]
table_c_top_md = {}
flips = []
for regime in beaching_regimes:
    export_dir = export_root / "028"
    rows = []
    for season in seasons:
        stem = f"connectivity_{regime}_{season}"
        paths = {
            "n_obs": export_dir / f"{stem}_unweighted_T{top_horizon}d_n_obs.csv",
            "unweighted": export_dir / f"{stem}_unweighted_T{top_horizon}d_emission_fraction.csv",
            "weighted": export_dir / f"{stem}_{member}_T{top_horizon}d_emission_fraction.csv",
        }
        missing = [p.name for p in paths.values() if not p.exists()]
        if missing:
            gaps.append(f"table C top targets, {regime} {season}: missing {', '.join(missing)}")
            continue
        counts = pd.read_csv(paths["n_obs"])
        unweighted = pd.read_csv(paths["unweighted"])
        weighted = pd.read_csv(paths["weighted"])
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

# %% [markdown]
# # Table D — dispersal distance
#
# `027` exports one per-source-hex quantile set per `(regime, season)`; the
# table reports the **median over source hexes** of each quantile level, so
# one number per (season, regime, level) summarises a whole map. The quantile
# levels are read off the GeoJSON columns rather than assumed.

# %%
distance_rows = []
for season in seasons:
    for regime in regimes:
        path = export_root / "027" / f"027_{regime}_{season}.geojson"
        if not path.exists():
            gaps.append(f"table D: missing {path.name}")
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
distance = pd.DataFrame(distance_rows)
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

# %% [markdown]
# # Table E — concentration at the long horizon
#
# `026` exports the cumulative relative density per hex at each horizon
# (`percentage` = share of that horizon's particle-time). Joined to the
# `024a` key for `mean_depth_m`, it gives the concentration of the cloud and
# the crude share of particle-time inside the Fucus viable band. A 6 km hex
# mean cannot resolve a 3 m contour, so the last column is a contrast
# between regimes, not a rate.

# %%
concentration_rows = []
for regime in regimes:
    for season in seasons:
        path = export_root / "026" / f"026_{regime}_{season}_T{horizon_days}d.geojson"
        if not path.exists():
            gaps.append(f"table E: missing {path.name}")
            continue
        frame = gpd.read_file(path).drop(columns="geometry").merge(key, on="hex_id", how="left")
        ranked = frame.percentage.sort_values(ascending=False)
        concentration_rows.append(
            {
                "regime": regime,
                "season": season,
                "max_hex_pct": float(ranked.iloc[0]),
                "top10_pct": float(ranked.iloc[:10].sum()),
                "top100_pct": float(ranked.iloc[:100].sum()),
                "occupied_hexes": len(frame),
                "viable_pct": float(frame.percentage[frame.mean_depth_m <= viable_depth_m].sum()),
            }
        )
concentration = pd.DataFrame(concentration_rows)
table_e = pd.DataFrame(
    {
        "regime": concentration.regime,
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

# %% [markdown]
# # Figures
#
# The selection is a **rule**, not a file list: each entry names the
# producing notebook and a filename *stem* built from the run's own
# parameters (regime, season, hex radius, member, horizon). Resolution takes
# the exact `<stem>.png` when it exists and otherwise the first
# `<stem>_*.png` in sorted order, so a consumer that later appends a field —
# a release-year span, say — still resolves. An entry that resolves to
# nothing is recorded as a gap rather than raising: a partially-run consumer
# matrix should still yield a report.
#
# The resolved PNGs are copied into `Headline/` numbered in reading order and
# downscaled to `figure_width_px`; the 180 mm / 300 dpi originals stay under
# `Figures/` as the manuscript versions.

# %%
figure_fields = {
    "regime": primary_regime,
    "other_regime": other_regime,
    "season": figure_season,
    "pooled": pooled_season,
    "radius": hex_radius,
    "member": member,
    "horizon": top_horizon,
}

figure_specs = [
    ("025", "HexRelativeDensity_{regime}_{pooled}_r{radius}m", "density_baltic_{regime}_{pooled}",
     "Basin-wide relative density: share of particle-time (left) and dilution per km² of hex "
     "water area (right). Where the particle-time sits, and how much of it is a coastal-hex artefact."),
    ("025", "HexRelativeDensityDE_{regime}_{pooled}_r{radius}m", "density_german_waters_{regime}_{pooled}",
     "The same two fields on the German-waters viewport, normalisation still basin-wide."),
    ("026", "TimeHorizonMaps_{regime}_{season}_r{radius}m_percentage", "horizons_{regime}_{season}",
     "Cumulative 10/20/50/100 d horizons, each normalised to its own 100 %: how the cloud fills the basin."),
    ("026", "TimeHorizonMaps_{other_regime}_{season}_r{radius}m_percentage", "horizons_{other_regime}_{season}",
     "The same horizons for the comparison regime. Read against the previous figure."),
    ("029", "WhereStranded_{regime}_r{radius}m_{season}_{member}", "where_stranded_{season}",
     "Stranded weight per hex (log, 4 decades): the stranding coastline of this season."),
    ("029", "BeachedFractionCurve_{regime}_r{radius}m_{pooled}_{member}", "beached_fraction_curve_{pooled}",
     "Cumulative beached fraction against age with the interannual min–max band; the curve behind table A."),
    ("031", "BeachingSweepStats_{regime}_r{radius}m_{pooled}", "beaching_sweep_stats_{pooled}",
     "Beached fraction, Gini, stranding-hex count, median stranding age and median travel distance "
     "across every member with interannual bars; the visual form of table B."),
    ("028", "SubbasinConnectivity_emission_{regime}_{season}_unweighted_T{horizon}d",
     "connectivity_{season}_T{horizon}d_unweighted",
     "Subbasin emission fraction at the top horizon, unweighted residence. Diagonal = self-retention."),
    ("028", "SubbasinConnectivity_emission_{regime}_{season}_{member}_T{horizon}d",
     "connectivity_{season}_T{horizon}d_survival_weighted",
     "The same matrix survival-weighted. Compare with the previous figure: the diagonal brightens, "
     "the far off-diagonal fades."),
    ("027", "DistanceQuantiles_{regime}_r{radius}m_{season}", "distance_quantiles_{regime}_{season}",
     "Per-source-hex quantiles of final crow-flies displacement; the map behind table D."),
    ("030", "SurvivalFraction_{regime}_r{radius}m_{season}_{member}", "surviving_fraction_{season}",
     "Surviving (un-beached) fraction per hex at each horizon on a fixed 0–1 scale: where propagules "
     "still drift after the coast has taken its share."),
]


# %%
def pick_figure(directory, stem):
    """Resolve a figure stem to a path: exact `<stem>.png`, else the first `<stem>_*.png`."""
    exact = directory / f"{stem}.png"
    if exact.exists():
        return exact
    extended = sorted(directory.glob(f"{stem}_*.png"))
    return extended[0] if extended else None


# Stale numbering from a previous run would otherwise survive a shrinking
# figure set, since the numbers are positions, not identities.
for stale in headline_dir.glob("[0-9][0-9]_*.png"):
    stale.unlink()

figure_notes = []
for notebook, stem_template, name, note in figure_specs:
    fields = re.findall(r"{(\w+)}", stem_template + name)
    if any(figure_fields.get(f) in (None, "") for f in fields):
        gaps.append(f"figure `{name}` ({notebook}): parameters give no value for one of {fields}")
        continue
    stem = stem_template.format(**figure_fields)
    source = pick_figure(figure_root / notebook, stem)
    if source is None:
        gaps.append(f"figure `{name}`: no `{stem}.png` under Figures/{notebook}/")
        continue
    target = headline_dir / f"{len(figure_notes) + 1:02d}_{name.format(**figure_fields)}.png"
    with Image.open(source) as image:
        if image.width > figure_width_px:
            height = round(image.height * figure_width_px / image.width)
            image.resize((figure_width_px, height), Image.LANCZOS).save(target)
        else:
            shutil.copyfile(source, target)
    figure_notes.append((target.name, notebook, source.name, note))
    print(f"{target.name} <- Figures/{notebook}/{source.name}")

# %% [markdown]
# # Assemble HEADLINE.md

# %%
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
    "release years pooled. Beaching and survival use the production member "
    f"`{member}` on {', '.join(f'`{r}`' for r in beaching_regimes)}; "
    f"table B sweeps all {len(sweep_members)} members on disk. Deposits reach "
    f"`beach_age_bin` {max_beach_age_bin}, i.e. the whole "
    f"{(max_beach_age_bin + 1) * age_bin_days} d window of the "
    f"{'/'.join(str(y) for y in beaching_years)} beaching store.",
    "",
    "## Coverage",
    "",
    "Release years per store and regime, as found on disk.",
    "",
    coverage_md,
    "",
]

sections += [
    "## Table A — beached fraction by release season, production member",
    "",
    f"`Σ weight[beach_hex ≥ 0] / Σ weight` from "
    f"`HexAggregates/HexAgg_beaching_r{hex_radius}m_<regime>_<year>_mMM_{member}.parquet`; "
    "`beach_hex = -1` is the never-beached residual. Percent.",
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
    f"Top-{n_top_origins} emitters per season at the {top_horizon} d horizon, from "
    f"`Exports/028/connectivity_<regime>_<season>_{{unweighted, {member}}}_T{top_horizon}d_"
    "emission_fraction.csv`; origins ranked by the matching `_n_obs.csv`. Target shares in percent.",
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
    f"From `Exports/026/026_<regime>_<season>_T{horizon_days}d.geojson` (cumulative ages "
    f"< {horizon_days} d; `percentage` = share of that horizon's particle-time), joined to "
    f"`HexAggregates/HexAgg_key_r{hex_radius}m.parquet` for `mean_depth_m`. The last column is a "
    "crude viable-zone proxy — a 6 km hex mean cannot resolve the real "
    f"{viable_depth_m:g} m band, so read it as a regime contrast, not a rate.",
    "",
    table_e_md,
    "",
    *concentration_recap,
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
        "Artefacts this report expected and did not find on disk:",
        "",
        *[f"- {g}" for g in gaps],
        "",
    ]

report_path = headline_dir / "HEADLINE.md"
report_path.write_text("\n".join(sections))
print(f"wrote {report_path} ({report_path.stat().st_size} bytes)")

# %% [markdown]
# # Summary

# %%
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
      f"D {len(table_d)} rows, E {len(table_e)} rows")
print(f"figures copied: {len(figure_notes)} of {len(figure_specs)} specified")
print(f"gaps recorded: {len(gaps)}")
for gap in gaps:
    print(f"  - {gap}")
