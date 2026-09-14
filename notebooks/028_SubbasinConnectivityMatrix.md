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

# Subbasin connectivity matrices

HELCOM subbasin→subbasin **residence** connectivity, pooled over release
years, for one release `season` and one age horizon at a time.

Two stores feed this, selected by the `member` parameter:

| `member` | store | value column | population |
|---|---|---|---|
| `""` | `024c` connectivity | `n_obs` | all particles (land seeds as `-1` rows) |
| rate-model tag | `024g` survconn | `w_obs` | real drifters only, weighted by `exp(−A)` |

The `024g` store also carries `n_obs` over *its* population, so the
weighted/unweighted comparison in the survival-weighted run is made on an
identical row set (see
[hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)).

Three views:

1. **Raw matrix**, log scale, pooled over every age bin — the magnitude
   picture. Rows differ by orders of magnitude.
2. **Emission fraction** `x / x.sum(axis=1)` — the fraction of an origin's
   particle-time that sits in each target, one figure per age horizon
   `T` (cumulative, half-open `age < T` days, the same convention 026 uses
   for its snapshot bins). Rows sum to 1, so origins are comparable.
3. **Interannual range** — the emission fraction recomputed per release
   year; mean / min / max per (origin, target, horizon).

Connectivity here is **residence** (particle-timesteps in the target
subbasin), not particle flux. Subbasin id `-1` = unnamed/outside and is
excluded from every matrix (reported as a dropped fraction). No Dask
cluster; parquet-only.

```python
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
```

# Parameters

```python tags=["parameters"]
# Read root of the hex-aggregate store; also the figure/export root.
output_root = "../output"
# Which regime's connectivity partitions to read; one regime per run.
regime = "surface_stokes"
# Hex radius of the store (built by 024a/024c/024g). Must match files on disk.
hex_radius = 6000
# Age-bin granularity of the store being read (must match 024/024g).
age_bin_days = 10
# Cumulative age horizons to draw, in days. Each must be a multiple of
# age_bin_days, and (for a survconn member) ≤ the 024g connectivity horizon.
time_horizons_days_csv = "10,20,50,100"
# Release season: one of DJF, MAM, JJA, SON, ALL.
season = "ALL"
# Rate-model member tag. "" reads the unweighted 024c connectivity store;
# any other string reads the 024g survconn store for that member and uses its
# survival-weighted `w_obs`. Opaque — never parsed.
member = ""
```

# Season months

Derived here, not in the parameters cell (primitives only there). Release
month is taken from `release_doy` against the partition's year, so the
conversion is leap-correct.

```python
SEASON_MONTHS = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11],
    "ALL": list(range(1, 13)),
}
```

# Parse parameters / layout

```python
output_root = Path(output_root)
store_root = output_root / "HexAggregates"

if season not in SEASON_MONTHS:
    raise ValueError(f"season {season!r} not one of {sorted(SEASON_MONTHS)}")
season_months = SEASON_MONTHS[season]

time_horizons_days = [int(x) for x in time_horizons_days_csv.split(",") if x]
for h in time_horizons_days:
    assert h % age_bin_days == 0, (
        f"horizon {h} d is not a multiple of age_bin_days {age_bin_days} d"
    )

member_tag = member or "unweighted"
value_col = "w_obs" if member else "n_obs"
print(f"regime={regime}, season={season} (months {season_months}), "
      f"member={member_tag}, value column {value_col!r}")

figure_dir = output_root / "Figures" / "028"
figure_dir.mkdir(parents=True, exist_ok=True)
export_dir = output_root / "Exports" / "028"
export_dir.mkdir(parents=True, exist_ok=True)
```

# Read the store

`024c` writes one partition per `(regime, year)`; `024g` one per
`(regime, year, month, member)`. Both carry `n_obs`; only `024g` carries
`w_obs`. The release year comes from the filename so the
`release_doy → month` conversion is leap-correct (same pattern as 026a).

```python
subbasin_id_to_name = {
    int(k): v
    for k, v in json.loads(
        (store_root / f"HexAgg_key_r{hex_radius}m.json").read_text()
    )["subbasin_id_to_name"].items()
}

if member:
    pattern = (
        f"HexAgg_survconn_r{hex_radius}m_{regime}_"
        f"[0-9][0-9][0-9][0-9]_m[0-9][0-9]_{member}.parquet"
    )
    year_re = re.compile(
        rf"HexAgg_survconn_r{hex_radius}m_{regime}_(\d{{4}})_m\d{{2}}_{re.escape(member)}\.parquet$"
    )
else:
    pattern = f"HexAgg_connectivity_r{hex_radius}m_{regime}_[0-9][0-9][0-9][0-9].parquet"
    year_re = re.compile(
        rf"HexAgg_connectivity_r{hex_radius}m_{regime}_(\d{{4}})\.parquet$"
    )

files = sorted(p for p in store_root.glob(pattern) if year_re.search(p.name))
if not files:
    raise FileNotFoundError(
        f"no partitions matching {pattern!r} at {store_root} — "
        + ("run 024g." if member else "run 024c.")
    )

parts = []
for f in files:
    year = int(year_re.search(f.name).group(1))
    df = pd.read_parquet(f).reset_index(drop=True)
    df["release_year"] = year
    df["release_month"] = pd.to_datetime(
        (year * 1000 + df["release_doy"].astype("int32")).astype(str), format="%Y%j"
    ).dt.month
    parts.append(df)
conn = pd.concat(parts, ignore_index=True)
years = sorted(conn["release_year"].unique())
print(f"read {len(files)} partition(s), years {years}; {len(conn):,} rows")

conn = conn[conn["release_month"].isin(season_months)]
if conn.empty:
    raise ValueError(f"no rows left after the {season} season filter")
print(f"season {season}: {len(conn):,} rows, "
      f"months {sorted(conn['release_month'].unique())}")
```

# Named-subbasin axis and the matrix reduction

The row/column order is fixed once, from the whole season-filtered frame, so
every per-horizon and per-year matrix is the same shape and directly
stackable.

```python
named_ids = sorted(
    set(conn.loc[conn["origin_subbasin"] >= 0, "origin_subbasin"])
    | set(conn.loc[conn["target_subbasin"] >= 0, "target_subbasin"])
)
id_names = [subbasin_id_to_name[i] for i in named_ids]
print(f"{len(named_ids)} named subbasins on each axis")


def matrix(df, horizon, column):
    """Named-subbasin matrix of `column`, cumulative over ages `< horizon` d.

    `horizon=None` pools every age bin in the store.
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


def draw(ax, values, norm=None):
    """imshow one matrix with subbasin names as tick labels."""
    im = ax.imshow(values, norm=norm)
    # Tick labels set at plot time only (AGENTS.md: don't mutate the frame
    # index for presentation).
    ax.set_xticks(range(len(id_names)))
    ax.set_xticklabels(id_names, rotation=90)
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
```

# View 1 — raw residence matrix, log scale

Pooled over every age bin in the store. `LogNorm` is not cosmetic: the
within-subbasin diagonal dominates by orders of magnitude, so a linear scale
washes out every cross-subbasin entry. Zeros are masked to NaN so empty
cells render blank rather than as the bottom of the colour scale.

```python
raw = matrix(conn, None, value_col)
total = conn[value_col].sum()
dropped = 1.0 - raw.to_numpy().sum() / total
print(f"raw matrix {raw.shape}, {dropped:.1%} of {value_col} dropped as -1 "
      f"origin/target")

fig, ax = plt.subplots(layout="constrained")
im = draw(ax, np.where(raw.to_numpy() > 0, raw.to_numpy(), np.nan), norm=LogNorm())
fig.colorbar(im, ax=ax, label=f"residence ({value_col})")
# Run context on the figure (centred on the page), horizon on the axes: an
# axes title is centred on the axes box, which the long subbasin tick labels
# push right, and the combined string then overruns the 180 mm page.
fig.suptitle(f"{regime} — {season} — {member_tag}")
ax.set_title("all ages")
fig.set_size_inches(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.80)
save(fig, f"SubbasinConnectivity_raw_{regime}_{season}_{member_tag}_Tall.png")
```

# View 2 — emission fraction per age horizon

`x / x.sum(axis=1)`: the share of an origin's cumulative particle-time
(ages `< T` days) that sits in each target. Rows sum to 1, so origins with
very different release counts are directly comparable. One figure per
horizon, on a `LogNorm` shared across horizons — the off-diagonal shares
span several decades, and a common scale makes the horizons comparable
figure to figure.

```python
emission = {h: emission_fraction(matrix(conn, h, value_col)) for h in time_horizons_days}

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
    im = draw(ax, np.where(values > 0, values, np.nan), norm=norm)
    fig.colorbar(im, ax=ax, label="emission fraction")
    fig.suptitle(f"{regime} — {season} — {member_tag}")
    ax.set_title(f"age < {h} d")
    fig.set_size_inches(PAGE_WIDTH_IN, PAGE_WIDTH_IN * 0.80)
    save(fig, f"SubbasinConnectivity_emission_{regime}_{season}_{member_tag}_T{h}d.png")
```

# View 3 — interannual range of the emission fraction

The emission fraction recomputed per release year, reduced to mean / min /
max per (origin, target, horizon), one figure per horizon. Cells are
annotated with the mean and the min–max range only when the matrix is small
enough for the text to fit (≲ 12×12); the HELCOM level-2 axis is larger than
that, so the annotation is normally off and the range lives in the CSV.

```python
ANNOTATE_MAX_N = 12

per_year = {
    (h, y): emission_fraction(
        matrix(conn[conn["release_year"] == y], h, value_col)
    )
    for h in time_horizons_days
    for y in years
}

spread = {}
for h in time_horizons_days:
    cube = np.stack([per_year[(h, y)].to_numpy() for y in years])
    spread[h] = {
        "mean": np.nanmean(cube, axis=0),
        "min": np.nanmin(cube, axis=0),
        "max": np.nanmax(cube, axis=0),
    }

for h in time_horizons_days:
    mean = spread[h]["mean"]
    fig, ax = plt.subplots(layout="constrained")
    im = draw(ax, np.where(mean > 0, mean, np.nan), norm=norm)
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
```

# CSV exports

Long-form (origin, target, value) rather than wide matrices: the recipients
join them against their own tables, and the by-year file needs a mean/min/max
triple per cell that a wide layout cannot hold.

```python
def long_form(values, column):
    return pd.DataFrame(
        {
            "origin_subbasin": np.repeat(id_names, len(id_names)),
            "target_subbasin": np.tile(id_names, len(id_names)),
            column: np.asarray(values).ravel(),
        }
    )


for h in time_horizons_days:
    stem = f"connectivity_{regime}_{season}_{member_tag}_T{h}d"

    n_obs_m = matrix(conn, h, "n_obs")
    long_form(n_obs_m.to_numpy(), "n_obs").to_csv(
        export_dir / f"{stem}_n_obs.csv", index=False
    )

    long_form(emission[h].to_numpy(), "emission_fraction").to_csv(
        export_dir / f"{stem}_emission_fraction.csv", index=False
    )

    by_year = long_form(spread[h]["mean"], "emission_fraction_mean")
    by_year["emission_fraction_min"] = spread[h]["min"].ravel()
    by_year["emission_fraction_max"] = spread[h]["max"].ravel()
    by_year["n_years"] = len(years)
    by_year.to_csv(export_dir / f"{stem}_emission_fraction_by_year.csv", index=False)

    print(f"  wrote {stem}_{{n_obs,emission_fraction,emission_fraction_by_year}}.csv")
```

# Validation / summary

```python
print(f"regime={regime}, hex_radius={hex_radius} m, age_bin_days={age_bin_days}, "
      f"season={season}, member={member_tag}, value={value_col}")
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
    diag = np.diag(e.to_numpy())
    print(
        f"\nage < {h:>4} d: total {value_col} {m.to_numpy().sum():,.0f}, "
        f"mean self-retention {np.nanmean(diag):.3f} "
        f"(min {np.nanmin(diag):.3f}, max {np.nanmax(diag):.3f})"
    )
    top = e.stack().sort_values(ascending=False)
    off = top[[i != j for i, j in top.index]].head(3)
    for (o, t), v in off.items():
        print(f"    top off-diagonal: {subbasin_id_to_name[o]} → "
              f"{subbasin_id_to_name[t]}  {v:.3f}")

if member:
    w = conn["w_obs"].sum()
    n = conn["n_obs"].sum()
    print(f"\nsurvival weighting: Σw_obs / Σn_obs = {w / n:.3f} over the "
          f"{season} season")
```
