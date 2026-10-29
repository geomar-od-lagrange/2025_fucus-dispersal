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
# # Subbasin connectivity matrices
#
# HELCOM subbasin→subbasin **residence** connectivity, pooled over release
# years, for one release `season` and one age horizon at a time.
#
# Two stores feed this, selected by the `member` parameter:
#
# | `member` | store | value column | population |
# |---|---|---|---|
# | `""` | `024c` connectivity | `n_obs` | all particles (land seeds as `-1` rows) |
# | rate-model tag | `024g` survconn | `w_obs` | real drifters only, weighted by `exp(−A)` |
#
# The `024g` store also carries `n_obs` over *its* population, so the
# weighted/unweighted comparison in the survival-weighted run is made on an
# identical row set (see
# [hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)).
#
# Three views:
#
# 1. **Raw matrix**, log scale, pooled over the whole reduction window
#    `age < max_age_days` — the magnitude picture. Rows differ by orders of
#    magnitude.
# 2. **Emission fraction** `x / x.sum(axis=1)` — the fraction of an origin's
#    particle-time that sits in each target, one figure per age horizon
#    `T` (cumulative, half-open `age < T` days, the same convention as 026's
#    horizons). Rows sum to 1, so origins are comparable.
# 3. **Interannual range** — the emission fraction recomputed per release
#    year; mean / min / max per (origin, target, horizon).
#
# In a survival-weighted run views 2 and 3 carry one extra target column,
# **stranded**: the particle-time the `exp(−A)` weighting removed from the
# origin's row, rebuilt from the `024e` beaching store of the same member.
# With it each row adds up to the origin's *unweighted* particle-time, so the
# row shows where the released material sits *and* how much of it beached,
# instead of normalising the beaching loss away.
#
# Both stores are cut to the same window, `age < max_age_days`: 024g records
# its own window in each partition's `df.attrs` (asserted to cover this one),
# while the 024c store is uncapped, so without the cut the weighted and the
# unweighted raw matrix would cover different ages.
#
# Connectivity here is **residence** (particle-timesteps in the target
# subbasin), not particle flux. Subbasin id `-1` = unnamed/outside and is
# excluded from every matrix (reported as a dropped fraction). No Dask
# cluster; parquet-only.

# %%
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the hex-aggregate store; also the figure/export root.
output_root = "../output"
# Which regime's connectivity partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024c/024g). Must match files on disk.
hex_radius = 6000
# Age-bin granularity of the store being read (must match 024c/024g; a 024g
# read asserts it against the partitions' attrs).
age_bin_days = 10
# Reduction window (d): both stores are cut to ages < max_age_days, so the
# weighted and unweighted matrices cover the same ages. A 024g read asserts
# its partitions' own window covers it.
max_age_days = 180
# Cumulative age horizons to draw, in days. Each must be a multiple of
# age_bin_days and ≤ max_age_days.
time_horizons_days_csv = "30,60,90,180"
# Release season: one of DJF, MAM, JJA, SON, ALL.
season = "ALL"
# Release years pooled: "first-last" (inclusive) or one year ("2019").
# Partitions of other years stay on disk unopened (2025 is partial).
release_years = "2016-2024"
# Rate-model member tag. "" reads the unweighted 024c connectivity store;
# any other string reads the 024g survconn store for that member and uses its
# survival-weighted `w_obs`. Opaque — never parsed.
member = ""

# %% [markdown]
# # Season months
#
# Derived here, not in the parameters cell (primitives only there). Release
# month is taken from `release_doy` against the partition's year, so the
# conversion is leap-correct.

# %%
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}

# %% [markdown]
# # Parse parameters / layout

# %%
output_root = Path(output_root)
store_root = output_root / "HexAggregates"

if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]
# `release_years` is "first-last" (inclusive) or a single year.
first_year, _, last_year = release_years.partition("-")
analysis_years = range(int(first_year), int(last_year or first_year) + 1)
assert len(analysis_years), f"release_years {release_years!r} is an empty range"
print(f"release years in the analysis: {analysis_years[0]}–{analysis_years[-1]}")

time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d"
    )
    assert h <= max_age_days, f"horizon {h} d is past max_age_days {max_age_days} d"
assert max_age_days % age_bin_days == 0, (
    f"max_age_days {max_age_days} is not a multiple of age_bin_days {age_bin_days}"
)

member_tag = member or "unweighted"
value_col = "w_obs" if member else "n_obs"
print(f"regime={regime}, season={season} (months {season_months}), "
      f"member={member_tag}, value column {value_col!r}")

figure_dir = output_root / "Figures" / "028"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "028"
export_dir.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # Read the store
#
# `024c` writes one partition per `(regime, year)`; `024g` one per
# `(regime, year, month, member)`, or one per `(regime, year, member)` when it
# ran with `release_month = 0`. Both carry `n_obs`; only `024g` carries
# `w_obs`. The release year comes from the filename so the
# `release_doy → month` conversion is leap-correct (same pattern as 026a).
#
# A 024g partition carries its build parameters as `df.attrs`. They are taken
# per file **before** the concat — pandas drops attrs that differ between the
# frames it concatenates, and `release_year`/`release_month`/`git_sha` always
# differ — and every pooled partition must agree on the window, the binning
# and the rate-model point.

# %%
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(
        (store_root / f"HexAgg_key_r{hex_radius}m.json").read_text()
    )["subbasin_id_to_name"].items()
}

if member:
    # The `_mMM` group is optional: 024g run with `release_month = 0` writes
    # one whole-year partition per year instead of twelve monthly ones (as
    # 029/030 handle for their stores). The glob is deliberately loose and the
    # regex is what selects, so each file is matched exactly once.
    pattern = f"HexAgg_survconn_r{hex_radius}m_{regime}_*_{member}.parquet"
    year_re = re.compile(
        rf"HexAgg_survconn_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
        rf"_{re.escape(member)}\.parquet$"
    )
else:
    pattern = f"HexAgg_connectivity_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet"
    year_re = re.compile(
        rf"HexAgg_connectivity_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$"
    )


def find_partitions(pattern, year_re, builder):
    """Partitions in `release_years` and their release years.

    A year must not carry both a whole-year and monthly partitions (the
    optional `_mMM` group, `year_re` group 2) — they cover the same releases,
    and reading both double-counts every row.
    """
    files = []
    years_not_read = set()
    for f in sorted(store_root.glob(pattern)):
        m = year_re.search(f.name)
        if m is None:
            continue
        if int(m.group(1)) in analysis_years:
            files.append(f)
        else:
            years_not_read.add(int(m.group(1)))
    if years_not_read:
        print(f"on disk but outside release_years, not read: {sorted(years_not_read)}")
    if not files:
        raise FileNotFoundError(
            f"no partitions matching {pattern!r} for release years {release_years} "
            f"at {store_root} — run {builder}."
        )
    groups = [year_re.search(f.name).groups() for f in files]
    years = [int(g[0]) for g in groups]
    whole_year = {int(g[0]) for g in groups if year_re.groups > 1 and g[1] is None}
    mixed = sorted(y for y in whole_year if years.count(y) > 1)
    if mixed:
        raise ValueError(
            f"years {mixed} carry both a whole-year and monthly {builder} "
            f"partitions at {store_root}; remove one form before reading"
        )
    return files, years


files, partition_years = find_partitions(
    pattern, year_re, "024g" if member else "024c"
)

# Keys every pooled 024g partition must share: the reduction window, the
# binning and the rate-model point. Provenance (builder, git_sha,
# slurm_job_id) and the partition's identity (release_year, release_month)
# legitimately differ and are not compared.
POOLED_KEYS = (
    "max_age_days", "age_bin_days", "output_dt_hours", "w_c", "tau_strong_hours",
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


def release_month(year, doy):
    """Calendar month of `release_doy` in `year` (leap-correct)."""
    return pd.to_datetime(
        (year * 1000 + doy.astype("int32")).astype(str), format="%Y%j"
    ).dt.month


parts = []
store_attrs = []
for f, year in zip(files, partition_years):
    df = pd.read_parquet(f).reset_index(drop=True)
    store_attrs.append(df.attrs)
    df = df[df["age_bin"] * age_bin_days < max_age_days].copy()
    df["release_year"] = year
    df["release_month"] = release_month(year, df["release_doy"])
    parts.append(df)
conn = pd.concat(parts, ignore_index=True)
print(f"read {len(files)} partition(s), years {sorted(set(partition_years))}; "
      f"{len(conn):,} rows at ages < {max_age_days} d")

if member:
    window = pooled_attrs(store_attrs, "024g survconn")
    assert window["age_bin_days"] == age_bin_days, (
        f"store age_bin_days {window['age_bin_days']} != {age_bin_days}"
    )
    assert window["max_age_days"] >= max_age_days, (
        f"024g window {window['max_age_days']} d < max_age_days {max_age_days} d "
        f"— rebuild 024g or lower max_age_days"
    )
    print(f"024g partitions agree on {window}")

# Release-doy coverage per year, taken before the season filter: it is what
# tells a partial release year from a hole in the build.
doy_span = conn.groupby("release_year")["release_doy"].agg(["min", "max"])

conn = conn[conn["release_month"].isin(season_months)]
if conn.empty:
    raise ValueError(f"no rows left after the {season} season filter")
# Years are taken after the season filter — a year whose partitions hold no
# release in this season must not appear on the interannual axis. A year can
# legitimately miss a season: the last release year is partial (releases stop
# in May, so it has no JJA and no SON). Report every absent year with the
# release-doy range its partitions do hold, so an unexpected gap is visible in
# the log rather than silent, and 032 reads the same coverage off the stores.
years = sorted(int(y) for y in conn["release_year"].unique())
for y in sorted(set(partition_years) - set(years)):
    print(f"no {season} releases in {y} (partitions cover release_doy "
          f"{doy_span.loc[y, 'min']}–{doy_span.loc[y, 'max']}); dropped from "
          f"the interannual axis")
print(f"season {season}: {len(conn):,} rows, years {years}, "
      f"months {sorted(conn['release_month'].unique())}")

# %% [markdown]
# # Stranded deposits (survival-weighted runs only)
#
# `w_obs` counts each hour at the surviving weight `exp(−A)`, so an origin's
# weighted row is short of its unweighted one by `Σ_hours (1 − exp(−A))`.
# `1 − exp(−A(t))` is the weight stranded by age `t`, i.e. the sum of the
# deposits made before `t`, so the shortfall over ages `[0, T)` is
# `Σ_deposits weight · (T − deposit_age)`. The `024e` beaching store of the
# same member holds those deposits per `(release_hex, release_doy,
# beach_age_bin)`; the deposit age is taken at its bin's midpoint, so the
# reconstruction holds to the 10-d bin resolution. Origin is the key's
# `helcom_subbasin` of `release_hex` — the assignment 024g makes — and both
# stores derive from the same 024d sidecars, so the populations match. Same
# release years and season filter as the connectivity rows; the partitions
# must agree with the 024g ones on the window, binning and rate-model point.

# %%
STRANDED = "stranded"

if member:
    beach_files, beach_years = find_partitions(
        f"HexAgg_beaching_r{hex_radius}m_{regime}_*_{member}.parquet",
        re.compile(
            rf"HexAgg_beaching_r{hex_radius}m_{regime}_(\d{{4}})(?:_m(\d{{2}}))?"
            rf"_{re.escape(member)}\.parquet$"
        ),
        "024e",
    )
    key = pd.read_parquet(
        store_root / f"HexAgg_key_r{hex_radius}m.parquet",
        columns=["hex_id", "helcom_subbasin"],
    )
    hex_subbasin = pd.Series(
        key["helcom_subbasin"].fillna(-1).astype(int).to_numpy(),
        index=key["hex_id"].astype(int),
    )

    deposit_parts = []
    beach_attrs = []
    for f, year in zip(beach_files, beach_years):
        df = pd.read_parquet(
            f, columns=["release_hex", "release_doy", "beach_age_bin", "weight"]
        )
        beach_attrs.append(df.attrs)
        # Deposits only (the never-beached residual has beach_age_bin -1).
        df = df[df["beach_age_bin"] >= 0].copy()
        df["release_year"] = year
        df["release_month"] = release_month(year, df["release_doy"])
        deposit_parts.append(df)
    beach_window = pooled_attrs(beach_attrs, "024e beaching")
    if beach_window != window:
        raise ValueError(
            f"024e partitions {beach_window} disagree with 024g {window}"
        )
    deposits = pd.concat(deposit_parts, ignore_index=True)
    deposits = deposits[deposits["release_month"].isin(season_months)]
    deposits["origin_subbasin"] = (
        hex_subbasin.reindex(deposits["release_hex"].to_numpy())
        .fillna(-1).astype(int).to_numpy()
    )
    samples_per_day = 24 / window["output_dt_hours"]
    print(f"read {len(beach_files)} 024e partition(s): {len(deposits):,} "
          f"deposit rows, Σweight {deposits['weight'].sum():,.0f} stranded")

# %% [markdown]
# # Named-subbasin axis and the matrix reduction
#
# The row/column order is fixed once, from the whole season-filtered frame, so
# every per-horizon and per-year matrix is the same shape and directly
# stackable.

# %%
named_ids = sorted(
    set(conn.loc[conn["origin_subbasin"] >= 0, "origin_subbasin"])
    | set(conn.loc[conn["target_subbasin"] >= 0, "target_subbasin"])
)
id_names = [subbasin_id_to_name[i] for i in named_ids]
print(f"{len(named_ids)} named subbasins on each axis")


def matrix(df, horizon, column):
    """Named-subbasin matrix of `column`, cumulative over ages `< horizon` d.

    `horizon=None` pools the whole window read (ages `< max_age_days`).
    """
    sel = df if horizon is None else df[df["age_bin"] * age_bin_days < horizon]
    sel = sel[(sel["origin_subbasin"] >= 0) & (sel["target_subbasin"] >= 0)]
    return (
        sel.pivot_table(
            index="origin_subbasin",
            columns="target_subbasin",
            values=column,
            aggfunc="sum",
        )
        .reindex(index=named_ids, columns=named_ids)
        .fillna(0.0)
    )


def emission_fraction(m):
    """Row-normalise: each origin's share of particle-time per target."""
    return m.div(m.sum(axis=1).replace(0.0, np.nan), axis=0)


def stranded(dep, horizon):
    """Particle-time stranding removed from ages `< horizon` d, per named
    origin, in the store's timestep units: Σ weight · (T − bin midpoint)."""
    sel = dep[dep["beach_age_bin"] * age_bin_days < horizon]
    age_mid = (sel["beach_age_bin"] + 0.5) * age_bin_days
    removed = sel["weight"] * (horizon - age_mid) * samples_per_day
    return (
        removed.groupby(sel["origin_subbasin"]).sum()
        .reindex(named_ids)
        .fillna(0.0)
    )


def emission_matrix(horizon, year=None):
    """Emission fraction at `horizon`, pooled or for one release year.

    A survival-weighted run appends the stranded column before normalising,
    so each row's denominator is the origin's unweighted particle-time.
    """
    rows = conn if year is None else conn[conn["release_year"] == year]
    m = matrix(rows, horizon, value_col)
    if member:
        dep = deposits if year is None else deposits[deposits["release_year"] == year]
        m[STRANDED] = stranded(dep, horizon)
    return emission_fraction(m)


def draw(ax, values, columns, norm=None):
    """imshow one origin × target matrix with subbasin names as tick labels.

    A trailing stranded column is set off from the subbasins by one blank
    (NaN) column.
    """
    values = np.asarray(values, dtype=float)
    xpos = list(range(len(columns)))
    if columns[-1] == STRANDED:
        values = np.insert(values, len(columns) - 1, np.nan, axis=1)
        xpos[-1] += 1
    im = ax.imshow(values, norm=norm)
    # Tick labels set at plot time only (AGENTS.md: don't mutate the frame
    # index for presentation).
    ax.set_xticks(xpos)
    ax.set_xticklabels([subbasin_id_to_name.get(c, c) for c in columns], rotation=90)
    ax.set_yticks(range(len(id_names)))
    ax.set_yticklabels(id_names)
    ax.set_xlabel("target subbasin")
    ax.set_ylabel("origin subbasin")
    return im


# Print-ready: full text width (180 mm) at 300 dpi. This overrides the
# project's "no figsize/dpi" default because these panels go into the
# manuscript at a fixed column width (rationale in docs/visualisations.md).
PAGE_WIDTH_IN = 180 / 25.4


def save(fig, name):
    fig.savefig(figure_dir / name, dpi=300)  # print-ready, see above
    print(f"  wrote {figure_dir / name}")
    plt.show()


# %% [markdown]
# # View 1 — raw residence matrix, log scale
#
# Pooled over the whole window, ages `< max_age_days`. `LogNorm` is not cosmetic: the
# within-subbasin diagonal dominates by orders of magnitude, so a linear scale
# washes out every cross-subbasin entry. Zeros are masked to NaN so empty
# cells render blank rather than as the bottom of the colour scale.

# %%
raw = matrix(conn, None, value_col)
total = conn[value_col].sum()
dropped = 1.0 - raw.to_numpy().sum() / total
print(f"raw matrix {raw.shape}, {dropped:.1%} of {value_col} dropped as -1 "
      f"origin/target")

fig, ax = plt.subplots(layout="constrained")
im = draw(ax, np.where(raw.to_numpy() > 0, raw.to_numpy(), np.nan), named_ids,
          norm=LogNorm())
fig.colorbar(im, ax=ax, label=f"residence ({value_col})")
# Run context on the figure (centred on the page), horizon on the axes: an
# axes title is centred on the axes box, which the long subbasin tick labels
# push right, and the combined string then overruns the 180 mm page.
fig.suptitle(f"{regime} — {season} — {member_tag}")
ax.set_title(f"age < {max_age_days} d")
fig.set_size_inches(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.80)
save(fig, f"SubbasinConnectivity_raw_{regime}_{season}_{member_tag}_T{max_age_days}d.png")

# %% [markdown]
# # View 2 — emission fraction per age horizon
#
# `x / x.sum(axis=1)`: the share of an origin's cumulative particle-time
# (ages `< T` days) that sits in each target. Rows sum to 1, so origins with
# very different release counts are directly comparable. In a
# survival-weighted run the row includes the **stranded** column (last,
# after a blank gap): the denominator is then the origin's unweighted
# particle-time and the stranded cell the share of it lost to beaching.
# One figure per
# horizon, on a `LogNorm` shared across horizons — the off-diagonal shares
# span several decades, and a common scale makes the horizons comparable
# figure to figure.

# %%
emission = {h: emission_matrix(h) for h in time_horizons_days}
targets = emission[time_horizons_days[0]].columns.tolist()

# Shared LogNorm across horizons, floored four decades below the peak: the
# thinnest off-diagonal shares run to 1e-7 and would otherwise stretch the
# scale until every meaningful cell reads the same colour (same treatment as
# 029's stranding maps).
stacked = np.concatenate([e.to_numpy().ravel() for e in emission.values()])
vmax = np.nanmax(stacked)
norm = LogNorm(vmin=vmax / 1e4, vmax=vmax)

# One figure per horizon rather than a panel grid: the axis carries 17 HELCOM
# subbasin names, which are illegible once four panels share a 180 mm page.
for h in time_horizons_days:
    values = emission[h].to_numpy()
    fig, ax = plt.subplots(layout="constrained")
    im = draw(ax, np.where(values > 0, values, np.nan), targets, norm=norm)
    fig.colorbar(im, ax=ax, label="emission fraction")
    fig.suptitle(f"{regime} — {season} — {member_tag}")
    ax.set_title(f"age < {h} d")
    fig.set_size_inches(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.80)
    save(fig, f"SubbasinConnectivity_emission_{regime}_{season}_{member_tag}_T{h}d.png")

# %% [markdown]
# # View 3 — interannual range of the emission fraction
#
# The emission fraction recomputed per release year, reduced to mean / min /
# max per (origin, target, horizon), one figure per horizon. Cells are
# annotated with the mean and the min–max range only when the matrix is small
# enough for the text to fit (≲ 12×12); the HELCOM level-2 axis is larger than
# that, so the annotation is normally off and the range lives in the CSV.

# %%
ANNOTATE_MAX_N = 12

per_year = {
    (h, y): emission_matrix(h, year=y) for h in time_horizons_days for y in years
}

# An origin that releases in no year at all has an all-NaN row in every
# year's matrix; reducing it with nanmin/nanmax is an all-NaN slice (a
# RuntimeWarning for a NaN that carries no information). Mask those rows out
# of the reduction and keep the mask — the CSV below drops them rather than
# writing bare NaN. Normalisation is per row, so a row is either all finite
# or all NaN; the mask is per origin.
spread = {}
for h in time_horizons_days:
    cube = np.stack([per_year[(h, y)].to_numpy() for y in years])
    emitting = np.isfinite(cube).any(axis=(0, 2))
    reduced = {}
    for name, reduce_ in (("mean", np.nanmean), ("min", np.nanmin), ("max", np.nanmax)):
        full = np.full(cube.shape[1:], np.nan)
        full[emitting] = reduce_(cube[:, emitting, :], axis=0)
        reduced[name] = full
    reduced["emitting"] = emitting
    spread[h] = reduced
    if not emitting.all():
        print(f"age < {h} d: {(~emitting).sum()} origin(s) with no release in "
              f"any year, dropped from the interannual export")

for h in time_horizons_days:
    mean = spread[h]["mean"]
    fig, ax = plt.subplots(layout="constrained")
    im = draw(ax, np.where(mean > 0, mean, np.nan), targets, norm=norm)
    if len(named_ids) <= ANNOTATE_MAX_N:
        for i in range(len(named_ids)):
            for j in range(len(named_ids)):
                if mean[i, j] > 0:
                    ax.text(
                        j, i,
                        f"{mean[i, j]:.2g}\n{spread[h]['min'][i, j]:.1g}–"
                        f"{spread[h]['max'][i, j]:.1g}",
                        ha="center", va="center",
                    )
    fig.colorbar(im, ax=ax, label="emission fraction (mean over years)")
    fig.suptitle(f"{regime} — {season} — {member_tag}")
    ax.set_title(f"age < {h} d — {len(years)}-year mean")
    fig.set_size_inches(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.80)
    save(fig, f"SubbasinConnectivity_interannual_{regime}_{season}_{member_tag}_T{h}d.png")

# %% [markdown]
# # CSV exports
#
# Long-form (origin, target, value) rather than wide matrices: the recipients
# join them against their own tables, and the by-year file needs a mean/min/max
# triple per cell that a wide layout cannot hold. In a survival-weighted run
# the two emission-fraction files carry `target_subbasin == "stranded"` rows
# alongside the subbasin targets; the `n_obs` file has none (it is the
# unweighted residence).

# %%
def long_form(values, columns, column):
    names = [subbasin_id_to_name.get(c, c) for c in columns]
    return pd.DataFrame(
        {
            "origin_subbasin": np.repeat(id_names, len(names)),
            "target_subbasin": np.tile(names, len(id_names)),
            column: np.asarray(values).ravel(),
        }
    )


for h in time_horizons_days:
    stem = f"connectivity_{regime}_{season}_{member_tag}_T{h}d"

    n_obs_m = matrix(conn, h, "n_obs")
    long_form(n_obs_m.to_numpy(), named_ids, "n_obs").to_csv(
        export_dir / f"{stem}_n_obs.csv", index=False
    )

    long_form(emission[h].to_numpy(), targets, "emission_fraction").to_csv(
        export_dir / f"{stem}_emission_fraction.csv", index=False
    )

    by_year = long_form(spread[h]["mean"], targets, "emission_fraction_mean")
    by_year["emission_fraction_min"] = spread[h]["min"].ravel()
    by_year["emission_fraction_max"] = spread[h]["max"].ravel()
    by_year["n_years"] = len(years)
    by_year = by_year[np.repeat(spread[h]["emitting"], len(targets))]
    by_year.to_csv(export_dir / f"{stem}_emission_fraction_by_year.csv", index=False)

    print(f"  wrote {stem}_{{n_obs,emission_fraction,emission_fraction_by_year}}.csv")

# %% [markdown]
# # Validation / summary

# %%
print(f"regime={regime}, hex_radius={hex_radius} m, age_bin_days={age_bin_days}, "
      f"max_age_days={max_age_days}, season={season}, member={member_tag}, "
      f"value={value_col}")
print(f"years pooled: {years}")
print(f"matrix: {len(named_ids)}×{len(named_ids)} named subbasins; "
      f"{dropped:.1%} of {value_col} outside them")

for h in time_horizons_days:
    m = matrix(conn, h, value_col)
    e = emission[h]
    # Origins that never release (present on the target axis only) have an
    # all-NaN row; every row that carries residence must normalise to 1.
    emitting = m.sum(axis=1) > 0
    rows = e[emitting].sum(axis=1)
    assert np.allclose(rows.to_numpy(), 1.0), rows
    diag = np.diag(e[named_ids].to_numpy())
    print(
        f"\nage < {h:>4} d: total {value_col} {m.to_numpy().sum():,.0f}, "
        f"mean self-retention {np.nanmean(diag):.3f} "
        f"(min {np.nanmin(diag):.3f}, max {np.nanmax(diag):.3f})"
    )
    top = e[named_ids].stack().sort_values(ascending=False)
    off = top[[i != j for i, j in top.index]].head(3)
    for (o, t), v in off.items():
        print(f"    top off-diagonal: {subbasin_id_to_name[o]} → "
              f"{subbasin_id_to_name[t]}  {v:.3f}")

if member:
    w = conn["w_obs"].sum()
    n = conn["n_obs"].sum()
    print(f"\nsurvival weighting: Σw_obs / Σn_obs = {w / n:.3f} over the "
          f"{season} season, ages < {max_age_days} d")

    # Closure of the stranded column: per named origin, the weighted
    # particle-time over every target (-1 included) plus the stranded
    # particle-time against the unweighted particle-time. Off zero only by
    # the 10-d deposit-age binning.
    print("\nstranded column: closure (Σw_obs + stranded) / Σn_obs − 1 per "
          "origin, all targets; stranded share of the row")
    for h in time_horizons_days:
        sel = conn[(conn["age_bin"] * age_bin_days < h) & (conn["origin_subbasin"] >= 0)]
        per_origin = sel.groupby("origin_subbasin")[["w_obs", "n_obs"]].sum()
        closure = (
            per_origin["w_obs"] + stranded(deposits, h).reindex(per_origin.index)
        ) / per_origin["n_obs"] - 1
        share = emission[h][STRANDED].dropna()
        print(f"  age < {h:>4} d: closure {closure.min():+.4f} … {closure.max():+.4f}; "
              f"stranded share {share.min():.3f} … {share.max():.3f} "
              f"(min {subbasin_id_to_name[share.idxmin()]}, "
              f"max {subbasin_id_to_name[share.idxmax()]})")
