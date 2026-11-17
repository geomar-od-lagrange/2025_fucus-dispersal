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

# Emission-fraction grid: with and without stranding, two horizons

One 2×2 figure per regime of 028's subbasin **emission fraction** (the
share of an origin's cumulative particle-time, ages `< T` days, that sits
in each target subbasin; rows sum to 1). Columns: **with stranding** (a
024g rate-model member; the row carries the extra *stranded* target, the
share of the origin's unweighted particle-time lost to beaching) and
**without stranding** (the unweighted 024c connectivity). Rows: the age
horizons of `horizons_csv`.

Reads only 028's CSV exports — no parquet, no Dask:
`Exports/028/connectivity_<regime>_<season>_<tag>_T<h>d_emission_fraction.csv`
for `tag` = `member` and `unweighted`. Runs after 028 for both tags and every
horizon drawn. The matrix rendering (west → east axis order, three-letter
codes, batlow on 1–2–5 classes over 10⁻⁵–1, the stranded column set off by
a gap with whole-percent labels) is 028's, so a cell reads the same here as
in 028's single-matrix figures.

```python
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cmcrameri import cm as cmc
from matplotlib.colors import BoundaryNorm
```

# Parameters

```python tags=["parameters"]
# Read root of 028's exports; write root for the figure.
output_root = "../output"
# Regime to draw; one regime per run (one figure).
regime = "surface_stokes"
# Release season 028 ran at: one of DJF, MAM, JJA, SON, ALL.
season = "ALL"
# Age horizons (d), one figure row each, top to bottom. 028 must have
# exported each of them.
horizons_csv = "30,180"
# Rate-model member of the "with stranding" column (028's 024g tag). The
# "without stranding" column is 028's unweighted run. Opaque — never parsed.
member = "step_wc0p1_ts3_tcinf_tf1_tw0_b750"
```

# Parse parameters / layout

```python
output_root = Path(output_root)
export_dir = output_root / "Exports" / "028"
figure_dir = output_root / "Figures" / "038"
figure_dir.mkdir(parents=True, exist_ok=True)

horizons = [int(x) for x in horizons_csv.split(",") if x]
assert horizons, f"horizons_csv {horizons_csv!r} names no horizon"
# Columns of the grid: (plain-language title, 028 tag).
COLUMNS = [("with stranding", member), ("without stranding", "unweighted")]

# Plain-language figure title; the run keys stay in the filename and prints.
REGIME_LABELS = {
    "surface": "surface",
    "surface_stokes": "surface + Stokes drift",
    "bottom": "bottom",
    "bottom_ub49": "bottom, bed friction",
}
title = REGIME_LABELS.get(regime, regime)
if season != "ALL":
    title = f"{title}, {season} releases"
print(f"regime={regime}, season={season}, horizons {horizons} d, "
      f"columns {[tag for _, tag in COLUMNS]}")
```

# Read 028's emission-fraction exports

Long form (origin, target, fraction), subbasins by full name; the stranded
target only in the member's file. Pivoted to origin × target here.

```python
STRANDED = "stranded"

emission = {}
for _, tag in COLUMNS:
    for h in horizons:
        f = export_dir / f"connectivity_{regime}_{season}_{tag}_T{h}d_emission_fraction.csv"
        if not f.exists():
            raise FileNotFoundError(f"{f} missing — run 028 for {regime}, {season}, {tag}")
        df = pd.read_csv(f)
        emission[(tag, h)] = df.pivot(
            index="origin_subbasin", columns="target_subbasin", values="emission_fraction"
        )
        print(f"read {f.name}: {df['origin_subbasin'].nunique()} origins × "
              f"{df['target_subbasin'].nunique()} targets")
```

# Axis order

West → east along the basin's transport path, as 028 (and 026a); a
subbasin in the exports without a code here is an error, not a silently
unordered row.

```python
# Three-letter subbasin codes, west → east: the axis order and tick labels.
SUBBASIN_CODE = {
    "Kattegat": "KAT", "The Sound": "SOU", "Great Belt": "GRB", "Kiel Bay": "KIE",
    "Bay of Mecklenburg": "MEC", "Arkona Basin": "ARK", "Bornholm Basin": "BOR",
    "Gdansk Basin": "GDA", "Eastern Gotland Basin": "EGB",
    "Western Gotland Basin": "WGB", "Gulf of Riga": "GRI",
    "Northern Baltic Proper": "NBP", "Gulf of Finland": "GFI", "Åland Sea": "ALA",
    "Bothnian Sea": "BOS", "The Quark": "QUA", "Bothnian Bay": "BOB",
}

present = set()
for m in emission.values():
    present |= set(m.index) | (set(m.columns) - {STRANDED})
uncoded = sorted(present - set(SUBBASIN_CODE))
assert not uncoded, f"subbasins without a code / west-east position: {uncoded}"
names = [n for n in SUBBASIN_CODE if n in present]
print(f"{len(names)} subbasins per axis, west → east: "
      f"{[SUBBASIN_CODE[n] for n in names]}")

# Targets per column: the member's file must carry the stranded target, the
# unweighted one must not (028 writes it only in a survival-weighted run).
targets = {
    tag: names + ([STRANDED] if tag != "unweighted" else []) for _, tag in COLUMNS
}
for (tag, h), m in emission.items():
    assert (STRANDED in m.columns) == (tag != "unweighted"), (
        f"{tag} T{h}d: stranded target {'missing' if tag != 'unweighted' else 'unexpected'}"
    )
    emission[(tag, h)] = m.reindex(index=names, columns=targets[tag])
```

# Matrix rendering (028's)

Colour classes, colormap and the per-panel draw as in 028: 1–2–5 classes
per decade over 10⁻⁵–1, shares under 10⁻⁵ in the `extend="min"` colour,
zero cells blank; the stranded column after one blank gap column, each cell
labelled with its whole percentage, the column labelled under the bottom row.

```python
# Colormap override (docs/visualisations.md, 028): batlow, multi-hue with
# monotonic lightness, so neighbouring log classes differ in hue as well as
# lightness and connectivity reads apart from the dilution maps' turbo.
MATRIX_CMAP = cmc.batlow

EMISSION_DECADES = 5
EMISSION_BOUNDS = [
    m * 10.0**e for e in range(-EMISSION_DECADES, 0) for m in (1, 2, 5)
] + [1.0]
norm = BoundaryNorm(EMISSION_BOUNDS, ncolors=MATRIX_CMAP.N, extend="min")
# Every edge labelled, log-axis style: decades as powers of ten, the 2 and 5
# edges by their mantissa.
emission_labels = [
    f"$10^{{{e}}}$" if m == 1 else f"{m}"
    for e in range(-EMISSION_DECADES, 0) for m in (1, 2, 5)
] + ["1"]


def code(c):
    """Tick label of an axis entry: the subbasin code, or the column name."""
    return SUBBASIN_CODE.get(c, c)


def draw(ax, values, columns, label_stranded):
    """imshow one origin × target matrix with subbasin codes as tick labels.

    A trailing stranded column is set off from the subbasins by one blank
    (NaN) column, labelled under the matrix with the target codes when
    `label_stranded` (the bottom row only: a label above the matrix would
    lift that panel's column title out of line with its neighbour's), and
    each of its cells carries
    the share as a whole percentage: the stranded shares crowd the top class
    or two, where the colour alone cannot rank them. The target codes are
    set vertical: at half page width a cell is narrower than a horizontal
    three-letter code.
    """
    values = np.asarray(values, dtype=float)
    xpos = list(range(len(columns)))
    if columns[-1] == STRANDED:
        values = np.insert(values, len(columns) - 1, np.nan, axis=1)
        xpos[-1] += 1
        for i, v in enumerate(values[:, -1]):
            if np.isfinite(v):
                ax.text(xpos[-1], i, f"{100 * v:.0f}", ha="center", va="center")
    im = ax.imshow(values, norm=norm, cmap=MATRIX_CMAP)
    # Tick labels set at plot time only (AGENTS.md: don't mutate the frame
    # index for presentation).
    labels = [code(c) for c in columns]
    if columns[-1] == STRANDED and label_stranded:
        labels[-1] = "str. %"
    ax.set_xticks(xpos, labels, rotation="vertical")
    ax.set_yticks(range(len(names)), [code(n) for n in names])
    return im
```

# The grid

Columns share the origin axis, so only the left column carries origin
codes; rows share the target axis per column, so only the bottom row
carries target codes. Width ratios follow the column counts (17 subbasins,
plus gap and stranded on the left), so every cell is the same square and
the rows line up across columns. Row labels sit left of the origin axis
label, column titles above the top row, one colorbar for all four panels.

```python
# Print-ready: full text width (180 mm) at 300 dpi. This overrides the
# project's "no figsize/dpi" default because the figure goes into the
# manuscript at a fixed column width (rationale in docs/visualisations.md).
PAGE_WIDTH_IN = 180 / 25.4
# Margins (in) around the cell grid: origin codes, axis and row labels at the
# left, the gap between columns; title and column titles above; target codes,
# the stranded label, axis label and colorbar below. The height follows from
# square cells.
MARGIN_LEFT_IN, MARGIN_RIGHT_IN, COLUMN_GAP_IN = 0.95, 0.05, 0.15
MARGIN_TOP_IN, ROW_GAP_IN, MARGIN_BOTTOM_IN = 0.9, 0.1, 1.0

n_cols = [len(targets[tag]) + (tag != "unweighted") for _, tag in COLUMNS]
cell = (PAGE_WIDTH_IN - MARGIN_LEFT_IN - MARGIN_RIGHT_IN - COLUMN_GAP_IN) / sum(n_cols)
fig, axes = plt.subplots(
    len(horizons), len(COLUMNS), layout="constrained", squeeze=False,
    gridspec_kw={"width_ratios": n_cols},
)
fig.set_size_inches(
    PAGE_WIDTH_IN,
    len(horizons) * len(names) * cell + MARGIN_TOP_IN + MARGIN_BOTTOM_IN
    + (len(horizons) - 1) * ROW_GAP_IN,
)
print(f"cell {cell * 25.4:.1f} mm; figure {fig.get_size_inches() * 25.4} mm")

for r, h in enumerate(horizons):
    for c, (col_title, tag) in enumerate(COLUMNS):
        ax = axes[r, c]
        values = emission[(tag, h)].to_numpy()
        im = draw(ax, np.where(values > 0, values, np.nan), targets[tag],
                  label_stranded=r == len(horizons) - 1)
        if r == 0:
            ax.set_title(col_title)
        if r < len(horizons) - 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("target subbasin")
        if c == 0:
            ax.set_ylabel("origin subbasin")
            # Row label left of the origin axis label, placed against it so
            # it follows the tick-label width.
            ax.annotate(
                f"age < {h} d", xy=(0, 0.5), xycoords=ax.yaxis.label,
                xytext=(-6, 0), textcoords="offset points",
                rotation="vertical", ha="right", va="center",
            )
        else:
            ax.tick_params(labelleft=False)

# The "str. %" label lengthens the left column's tick labels; level the
# target-axis labels of the bottom row.
fig.align_xlabels()
# pad 0.01 (default 0.05 is a fraction of the grid height): the default left a
# ~8 mm blank band between the target-axis labels and the bar (as 034).
cb = fig.colorbar(im, ax=axes, location="bottom", label="emission fraction", pad=0.01)
cb.set_ticks(EMISSION_BOUNDS, labels=emission_labels)
fig.suptitle(title)

name = f"EmissionGrid_{regime}_{season}.png"
fig.savefig(figure_dir / name, dpi=300)  # print-ready, see above
print(f"wrote {figure_dir / name}")
plt.show()
```

# Summary

```python
for h in horizons:
    s = emission[(member, h)][STRANDED].dropna()
    print(f"age < {h:>4} d: stranded share {s.min():.3f} … {s.max():.3f} "
          f"(min {s.idxmin()}, max {s.idxmax()})")
    for _, tag in COLUMNS:
        m = emission[(tag, h)]
        rows = m.dropna(how="all").sum(axis=1)
        assert np.allclose(rows, 1.0), (tag, h, rows)
        diag = np.diag(m[names].to_numpy())
        print(f"    {tag}: {len(rows)} emitting origins, self-retention "
              f"{np.nanmin(diag):.3f} … {np.nanmax(diag):.3f}")
```
