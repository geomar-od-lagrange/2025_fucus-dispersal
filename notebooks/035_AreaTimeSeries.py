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
# # Stranding and the area holding the cloud vs. elapsed time
#
# Three time series on one age axis, one line per regime. The cumulative
# percentage of released drifters **stranded** by age `T` (029's stranding
# curve, one 024e rate-model member shared by all regimes), and the R1 area of
# plans/done/horizons_and_weighting.md D6: the smallest
# water area holding a share `q` of the particle-time accumulated by horizon
# `T` (cumulative, ages in `[0, T)`), and its square root as a length scale.
# Stranding is cumulative by age just as the area is by horizon, so "stranded
# by T" and "area at horizon T" share the x axis.
#
# Reads only CSV exports — no parquet, no Dask: 026's
# `Exports/026/026_areas_<regime>_<season>_<tag>.csv`, each regime at its own
# weighting tag (a 024f rate-model member, or `unweighted` for the 024
# counts), and 029's `Exports/029/stranded_curve_<regime>_<season>_<member>.csv`
# per regime. Runs after 026 and 029 for every regime it draws.
#
# `pooled` (every release year summed) is the line; the per-release-year
# min–max of the same reduction is the shaded band. A pooled value need not
# lie inside the band: the pooled ranking is over summed weights, not an
# average of yearly areas.

# %%
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of 026's and 029's exports; write root for the figure.
output_root = "../output"
# Hex radius 026 ran at. Labels the figure only: the 026 area CSV names
# carry no radius, so this must match the 026 run that wrote them.
hex_radius = 6000
# Regimes to draw, in drawing order, and each one's 026 weighting tag (a 024f
# rate-model member, or "unweighted"), position by position.
regimes_csv = "surface_stokes,bottom,bottom_ub49"
tags_csv = "step_wc0p1_ts3_tcinf_tf1_tw0_b750,step_wc0p1_ts3_tcinf_tf1_tw0_b750,step_wc0p1_ts3_tcinf_tf1_tw0_b750"
# 024e rate-model member of the stranding panel, shared by every regime of
# regimes_csv: the 029 runs whose stranding curves are drawn.
stranded_member = "step_wc0p1_ts3_tcinf_tf1_tw0_b750"
# Release season: DJF, MAM, JJA, SON or ALL.
season = "ALL"
# Figure layout: "regimes" draws `season`, one line per regime (stranded %,
# area, √area); "seasons" draws the four release seasons as lines, one row
# per regime (stranded %, √area), and ignores `season`.
layout = "regimes"
# Share of the particle-time whose holding area is drawn; must be one of the
# q values 026 exported.
area_quantile = 0.9

# %% [markdown]
# # Parse parameters

# %%
SEASONS = ("DJF", "MAM", "JJA", "SON", "ALL")
if season not in SEASONS:
    raise ValueError(f"season {season!r} not one of {SEASONS}")
LAYOUTS = ("regimes", "seasons")
if layout not in LAYOUTS:
    raise ValueError(f"layout {layout!r} not one of {LAYOUTS}")
# The seasons the run reads: one, or the four release seasons.
read_seasons = [season] if layout == "regimes" else ["DJF", "MAM", "JJA", "SON"]
regimes = [r.strip() for r in regimes_csv.split(",") if r.strip()]
tags = [t.strip() for t in tags_csv.split(",") if t.strip()]
if len(regimes) != len(tags):
    raise ValueError(f"{len(regimes)} regimes {regimes} but {len(tags)} tags {tags}")

output_root = Path(output_root)
export_dir = output_root / "Exports" / "026"
stranded_dir = output_root / "Exports" / "029"

# Print-ready output: one full text-width (180 mm) figure at 300 dpi, so it
# lands in the manuscript at its final size (rationale in
# docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300

figure_dir = output_root / "Figures" / "035"
figure_dir.mkdir(parents=True, exist_ok=True)

print(f"layout {layout}, seasons {read_seasons}, q = {area_quantile:g}, regimes {list(zip(regimes, tags))}; "
      f"stranding member {stranded_member}")

# %% [markdown]
# # Read the 026 area exports
#
# One CSV per regime, columns `horizon_days, q, scope, area_km2, length_km`,
# `scope` = `pooled` or a release year. A missing CSV is reported and that
# regime left out, so a partial consumer run (one regime's 026 done) still
# draws; no CSV at all raises. A CSV without rows at `area_quantile` raises:
# that is a parameter mismatch with 026, not a missing run.

# %%
frames = []
missing = []
for (regime, tag), s in itertools.product(zip(regimes, tags), read_seasons):
    path = export_dir / f"026_areas_{regime}_{s}_{tag}.csv"
    if not path.exists():
        missing.append(path)
        print(f"MISSING: {path} — {regime} ({tag}, {s}) is NOT drawn; run 026 for it")
        continue
    df = pd.read_csv(path, dtype={"scope": str})
    df = df[np.isclose(df["q"], area_quantile)]
    if df.empty:
        raise ValueError(f"{path.name}: no rows at q = {area_quantile:g}")
    frames.append(df.assign(regime=regime, tag=tag, season=s))
    years = sorted(s for s in df["scope"].unique() if s != "pooled")
    print(
        f"{path.name}: horizons {df['horizon_days'].min()}–{df['horizon_days'].max()} d "
        f"({df['horizon_days'].nunique()}), release years {years}"
    )
if not frames:
    raise FileNotFoundError(f"none of the 026 area exports exist: {[str(p) for p in missing]}")
areas = pd.concat(frames, ignore_index=True)
drawn = [r for r in regimes if r in set(areas["regime"])]

# %% [markdown]
# # Read the 029 stranding curves
#
# One CSV per regime, columns `age_days, stranded_pct, stranded_pct_year_min,
# stranded_pct_year_max, n_years`, `age_days` the age-bin right edge `T`
# (stranding at ages < T), `#` header lines carrying the run's provenance. As
# for the area CSVs, a missing curve is reported and that regime left out of
# the stranding panel; no curve at all raises, since the figure's top panel
# would have nothing to show.

# %%
# Keyed (regime, season).
stranded = {}
stranded_missing = []
for regime, s in itertools.product(regimes, read_seasons):
    path = stranded_dir / f"stranded_curve_{regime}_{s}_{stranded_member}.csv"
    if not path.exists():
        stranded_missing.append(path)
        print(f"MISSING: {path} — {regime} ({s}) stranding is NOT drawn; run 029 for it")
        continue
    curve = pd.read_csv(path, comment="#").set_index("age_days")
    stranded[(regime, s)] = curve
    print(f"{path.name}: ages {curve.index.min()}–{curve.index.max()} d ({len(curve)}), "
          f"{int(curve['n_years'].iloc[0])} release years")
if not stranded:
    raise FileNotFoundError(
        f"none of the 029 stranding curves exist: {[str(p) for p in stranded_missing]}"
    )
# %% [markdown]
# # Pooled value and interannual range
#
# Per regime and horizon: the pooled value and the min–max over the release
# years, for both the area and its square root (min/max commute with the
# square root, so the length band is the area band's root).

# %%
pooled = areas[areas["scope"] == "pooled"].set_index(["regime", "season", "horizon_days"])[
    ["area_km2", "length_km"]
]
year_range = (
    areas[areas["scope"] != "pooled"]
    .groupby(["regime", "season", "horizon_days"])[["area_km2", "length_km"]]
    .agg(["min", "max"])
)
year_range.columns = [f"{value}_year_{stat}" for value, stat in year_range.columns]
table = pooled.join(year_range).sort_index()
print(
    f"area holding {100 * area_quantile:g} % of particle-time, seasons {read_seasons}; "
    f"pooled with per-release-year min–max"
)
print(table.round(1).to_string())

# %% [markdown]
# # Figure
#
# Three panels stacked on a shared age axis, from the top: (a) stranded % of
# released, (b) area, (c) `√area`. One colour per regime, by its position in
# `regimes_csv` in the default cycle (so a regime keeps its colour when
# another one's export is missing), shared by the regime's line, its band and
# its legend entry in every panel. One default legend in panel (a), upper
# left, names the regimes for all three panels: with three regimes the
# lines start together at the first horizon, where inline labels collided,
# and (a)'s upper left stays empty (stranding rises from zero). The
# weighting of each regime goes into the caption. Panel letters sit inside
# the panels, top right. No title and no band key: the season and the band
# (per-release-year min–max) go in the caption.

# %%
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
# Plain-language line labels; the run keys stay in filenames and prints.
REGIME_LABELS = {
    "surface": "surface",
    "surface_stokes": "surface + Stokes drift",
    "bottom": "bottom",
    "bottom_ub49": "bottom + bed friction",
}


def regime_color(regime):
    """The regime's default-cycle colour by its position in `regimes_csv`."""
    return colors[regimes.index(regime)]


def series_panel(ax, lines, ylabel):
    """Pooled line + per-year band per series, y from zero; the line carries
    the series label for a legend. `lines` is a list of (label, color, frame,
    value column, min, max)."""
    for label, color, t, value, lo, hi in lines:
        # The band takes the line's colour so each regime reads as one
        # object; alpha keeps the other regime's line visible through it.
        ax.fill_between(t.index, t[lo], t[hi], color=color, alpha=0.3)
        # Markers mark the horizons the exported table actually holds.
        ax.plot(t.index, t[value], marker="o", color=color, label=label)
    ax.set_ylabel(ylabel)
    # From zero, so the regimes compare as ratios by eye; the top leaves a
    # text line of headroom over the highest band (at a third of the page
    # width per panel) for the panel letter.
    ax.set_ylim(0, 1.2 * max(f[hi].max() for _, _, f, _, _, hi in lines))


def area_lines(value, s):
    """The `series_panel` series of one area quantity in season `s`, one per
    drawn regime."""
    out = []
    for regime in drawn:
        tag = areas.loc[areas["regime"] == regime, "tag"].iloc[0]
        weighting = "unweighted" if tag == "unweighted" else "survival-weighted"
        out.append((
            f"{REGIME_LABELS.get(regime, regime)} ({weighting})", regime_color(regime),
            table.loc[(regime, s)], value, f"{value}_year_min", f"{value}_year_max",
        ))
    return out


if layout == "regimes":
    PANELS = [
        (
            [(REGIME_LABELS.get(regime, regime), regime_color(regime), curve,
              "stranded_pct", "stranded_pct_year_min", "stranded_pct_year_max")
             for (regime, _), curve in stranded.items()],
            "stranded (% of released)",
        ),
        (area_lines("area_km2", season), f"area holding {100 * area_quantile:g} % (km²)"),
        (area_lines("length_km", season), f"√area holding {100 * area_quantile:g} % (km)"),
    ]

    fig, axes = plt.subplots(len(PANELS), 1, sharex=True, layout="constrained")
    # Print-ready width; three stacked line panels at a third of the page width
    # each.
    fig.set_size_inches(FIGURE_WIDTH_IN, FIGURE_WIDTH_IN)
    for ax, letter, (lines, ylabel) in zip(axes, "abc", PANELS):
        series_panel(ax, lines, ylabel)
        # Top right; the legend takes (a)'s top left.
        ax.text(0.99, 0.96, f"({letter})", transform=ax.transAxes, ha="right", va="top")
    axes[0].legend(loc="upper left")
    axes[-1].set_xlabel("age / horizon T (d)")
    fig_path = figure_dir / f"AreaTimeSeries_{season}_r{hex_radius}m.png"
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

# %% [markdown]
# # Season figure
#
# `layout = "seasons"`: one row per regime, two columns, (left) stranded %
# against age and (right) `√area` holding `q` against horizon; one line per
# release season, colour by season from the default cycle (DJF, MAM, JJA,
# SON in cycle order). Lines only, no interannual band: four overlapping
# bands per panel hide the lines, and the range is in the "regimes" figure
# and the HEADLINE tables. y from zero, one y range per column so the
# regimes compare by eye. The regime names the row (y label), one legend in
# the top-right panel's empty lower part names the seasons, panel letters
# top right. Same
# 180 mm square page as the "regimes" figure.

# %%
if layout == "seasons":
    season_colors = dict(zip(read_seasons, colors))
    fig, axes = plt.subplots(
        len(regimes), 2, sharex="col", sharey="col", layout="constrained", squeeze=False
    )
    fig.set_size_inches(FIGURE_WIDTH_IN, FIGURE_WIDTH_IN)
    drawn_pairs = set(table.index.droplevel("horizon_days"))
    for row, regime in enumerate(regimes):
        ax_s, ax_l = axes[row]
        for s in read_seasons:
            if (regime, s) in stranded:
                curve = stranded[(regime, s)]
                ax_s.plot(curve.index, curve["stranded_pct"], marker="o",
                          color=season_colors[s], label=s)
            if (regime, s) in drawn_pairs:
                t = table.loc[(regime, s)]
                ax_l.plot(t.index, t["length_km"], marker="o", color=season_colors[s], label=s)
        ax_s.set_ylabel(f"{REGIME_LABELS.get(regime, regime)}\nstranded (% of released)")
        ax_l.set_ylabel(f"√area holding {100 * area_quantile:g} % (km)")
    # From zero with headroom for the panel letters, as in the "regimes"
    # figure; the rows share each column's y range.
    for col in range(2):
        axes[0, col].set_ylim(0, 1.2 * axes[0, col].get_ylim()[1])
    for ax, letter in zip(axes.flat, "abcdefghij"):
        ax.text(0.99, 0.96, f"({letter})", transform=ax.transAxes, ha="right", va="top")
    axes[0, 1].legend(loc="lower right")
    axes[-1, 0].set_xlabel("age (d)")
    axes[-1, 1].set_xlabel("horizon T (d)")
    fig_path = figure_dir / f"AreaTimeSeries_seasons_r{hex_radius}m.png"
    fig.savefig(fig_path, dpi=FIGURE_DPI)
    print(f"wrote {fig_path}")
    plt.show()

# %% [markdown]
# # Validation prints

# %%
print(f"area drawn: {drawn}; stranding drawn: {sorted(set(r for r, _ in stranded))} "
      f"({stranded_member}), seasons {read_seasons}")
if missing or stranded_missing:
    print(f"NOT drawn ({len(missing) + len(stranded_missing)} missing export(s)):")
    for p in missing + stranded_missing:
        print(f"  {p}")
area_pairs = set(table.index.droplevel("horizon_days"))
for regime, s in itertools.product(drawn, read_seasons):
    if (regime, s) not in area_pairs:
        continue
    t = table.loc[(regime, s)]
    outside = int(
        ((t["area_km2"] < t["area_km2_year_min"]) | (t["area_km2"] > t["area_km2_year_max"])).sum()
    )
    print(f"{regime} {s}: {len(t)} horizons, pooled outside the per-year range at {outside}; "
          f"√area {t['length_km'].iloc[0]:.0f} km at {t.index[0]} d, "
          f"{t['length_km'].iloc[-1]:.0f} km at {t.index[-1]} d")
for (regime, s), curve in stranded.items():
    area_ages = set(table.loc[(regime, s)].index) if (regime, s) in area_pairs else set()
    print(f"{regime} {s} stranded by {curve.index.max()} d: "
          f"{curve['stranded_pct'].iloc[-1]:.1f} % "
          f"[{curve['stranded_pct_year_min'].iloc[-1]:.1f}–"
          f"{curve['stranded_pct_year_max'].iloc[-1]:.1f}]; "
          f"area and stranding ages agree: {set(curve.index) == area_ages}")
