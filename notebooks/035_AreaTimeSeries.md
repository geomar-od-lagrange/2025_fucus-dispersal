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

# Area holding the cloud vs. elapsed time

The R1 area of plans/horizons_and_weighting.md D6 as a time series, one
line per regime: the smallest water area holding a share `q` of the
particle-time accumulated by horizon `T` (cumulative, ages in `[0, T)`),
and its square root as a length scale. Reads only 026's area exports
`Exports/026/026_areas_<regime>_<season>_<tag>.csv` — no parquet, no
Dask. Each regime is read at its own weighting tag (a 024f rate-model
member, or `unweighted` for the 024 counts), as 026 wrote it.

`pooled` (every release year summed) is the line; the per-release-year
min–max of the same reduction is the shaded band. A pooled value need not
lie inside the band: the pooled ranking is over summed weights, not an
average of yearly areas.

```python
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
```

# Parameters

```python tags=["parameters"]
# Read root of 026's exports; write root for the figure.
output_root = "../output"
# Hex radius 026 ran at. Labels the figure only: the 026 area CSV names
# carry no radius, so this must match the 026 run that wrote them.
hex_radius = 6000
# Regimes to draw, in legend order, and each one's 026 weighting tag (a 024f
# rate-model member, or "unweighted"), position by position.
regimes_csv = "surface_stokes,bottom"
tags_csv = "step_wc0p1_ts3_tcinf_tf1_tw0,unweighted"
# Release season: DJF, MAM, JJA, SON or ALL.
season = "ALL"
# Share of the particle-time whose holding area is drawn; must be one of the
# q values 026 exported.
area_quantile = 0.9
```

# Parse parameters

```python
SEASONS = ("DJF", "MAM", "JJA", "SON", "ALL")
if season not in SEASONS:
    raise ValueError(f"season {season!r} not one of {SEASONS}")
regimes = [r.strip() for r in regimes_csv.split(",") if r.strip()]
tags = [t.strip() for t in tags_csv.split(",") if t.strip()]
if len(regimes) != len(tags):
    raise ValueError(f"{len(regimes)} regimes {regimes} but {len(tags)} tags {tags}")

output_root = Path(output_root)
export_dir = output_root / "Exports" / "026"

# Print-ready output: one full text-width (180 mm) figure at 300 dpi, so it
# lands in the manuscript at its final size (rationale in
# docs/visualisations.md).
FIGURE_WIDTH_IN = 180 / 25.4
FIGURE_DPI = 300

figure_dir = output_root / "Figures" / "035"
figure_dir.mkdir(parents=True, exist_ok=True)

print(f"season {season}, q = {area_quantile:g}, regimes {list(zip(regimes, tags))}")
```

# Read the 026 area exports

One CSV per regime, columns `horizon_days, q, scope, area_km2, length_km`,
`scope` = `pooled` or a release year. A missing CSV is reported and that
regime left out, so a partial consumer run (one regime's 026 done) still
draws; no CSV at all raises. A CSV without rows at `area_quantile` raises:
that is a parameter mismatch with 026, not a missing run.

```python
frames = []
missing = []
for regime, tag in zip(regimes, tags):
    path = export_dir / f"026_areas_{regime}_{season}_{tag}.csv"
    if not path.exists():
        missing.append(path)
        print(f"MISSING: {path} — {regime} ({tag}) is NOT drawn; run 026 for it")
        continue
    df = pd.read_csv(path, dtype={"scope": str})
    df = df[np.isclose(df["q"], area_quantile)]
    if df.empty:
        raise ValueError(f"{path.name}: no rows at q = {area_quantile:g}")
    frames.append(df.assign(regime=regime, tag=tag))
    years = sorted(s for s in df["scope"].unique() if s != "pooled")
    print(
        f"{path.name}: horizons {df['horizon_days'].min()}–{df['horizon_days'].max()} d "
        f"({df['horizon_days'].nunique()}), release years {years}"
    )
if not frames:
    raise FileNotFoundError(f"none of the 026 area exports exist: {[str(p) for p in missing]}")
areas = pd.concat(frames, ignore_index=True)
drawn = [r for r in regimes if r in set(areas["regime"])]
```

# Pooled value and interannual range

Per regime and horizon: the pooled value and the min–max over the release
years, for both the area and its square root (min/max commute with the
square root, so the length band is the area band's root).

```python
pooled = areas[areas["scope"] == "pooled"].set_index(["regime", "horizon_days"])[
    ["area_km2", "length_km"]
]
year_range = (
    areas[areas["scope"] != "pooled"]
    .groupby(["regime", "horizon_days"])[["area_km2", "length_km"]]
    .agg(["min", "max"])
)
year_range.columns = [f"{value}_year_{stat}" for value, stat in year_range.columns]
table = pooled.join(year_range).sort_index()
print(
    f"area holding {100 * area_quantile:g} % of particle-time, season {season}; "
    f"pooled with per-release-year min–max"
)
print(table.round(1).to_string())
```

# Figure

Two panels, area and `√area`, against the horizon `T`; one colour per
regime, by its position in `regimes_csv` in the default cycle (so a regime
keeps its colour when another one's export is missing), shared by the
regime's line and its band.

```python
PANELS = {
    "area_km2": f"area holding {100 * area_quantile:g} % (km²)",
    "length_km": f"√area holding {100 * area_quantile:g} % (km)",
}
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

fig, axes = plt.subplots(1, len(PANELS), layout="constrained")
# Print-ready width; a line panel at the default 4:3 per-panel aspect.
fig.set_size_inches(FIGURE_WIDTH_IN, 0.75 * FIGURE_WIDTH_IN / len(PANELS))
legend_handles = {}
for ax, (value, ylabel) in zip(axes, PANELS.items()):
    for regime in drawn:
        color = colors[regimes.index(regime)]
        t = table.loc[regime]
        tag = areas.loc[areas["regime"] == regime, "tag"].iloc[0]
        weighting = "unweighted" if tag == "unweighted" else "survival-weighted"
        # The band takes the line's colour so each regime reads as one
        # object; alpha keeps the other regime's line visible through it.
        band = ax.fill_between(
            t.index, t[f"{value}_year_min"], t[f"{value}_year_max"],
            color=color, alpha=0.3,
        )
        # Markers mark the horizons the 026 table actually holds.
        (line,) = ax.plot(t.index, t[value], marker="o", color=color)
        legend_handles[f"{regime} ({weighting})"] = (band, line)
    ax.set_xlabel("horizon T (d), ages in [0, T)")
    ax.set_ylabel(ylabel)
    # From zero, so the two regimes' areas compare as ratios by eye.
    ax.set_ylim(0, None)
# One entry per regime (band + pooled line), below the panels: inside an
# axes it covers the curves at print width.
fig.legend(
    list(legend_handles.values()), list(legend_handles),
    loc="outside lower center", ncols=len(legend_handles),
    title="pooled (line) and per-year min–max (band)",
)
fig.suptitle(f"release season {season}")
fig_path = figure_dir / f"AreaTimeSeries_{season}_r{hex_radius}m.png"
fig.savefig(fig_path, dpi=FIGURE_DPI)
print(f"wrote {fig_path}")
plt.show()
```

# Validation prints

```python
print(f"drawn: {drawn}")
if missing:
    print(f"NOT drawn ({len(missing)} missing export(s)):")
    for p in missing:
        print(f"  {p}")
for regime in drawn:
    t = table.loc[regime]
    outside = int(
        ((t["area_km2"] < t["area_km2_year_min"]) | (t["area_km2"] > t["area_km2_year_max"])).sum()
    )
    print(f"{regime}: {len(t)} horizons, pooled outside the per-year range at {outside}")
```
