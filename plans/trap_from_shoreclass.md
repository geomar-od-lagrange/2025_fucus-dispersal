# Activate `trap` from the shoreclass sub-segment classification

Re-establish `trap` in the beaching rate as a **continuous substrate factor**,
driven by the `flat_fraction` table produced by the sidecar repository
<https://github.com/geomar-od-lagrange/2025_fucus-dispersal_shoreclass>.

The classification itself is not built here. This repository derives none of
it: it consumes one table and turns it into a rate multiplier.

## What the sidecar delivers

One row per sub-segment of the BSH tracer-cell outline — the staircase of cell
faces bounding the model's wet cells, subdivided at k=2 into 65,876
sub-segments of 420–465 m, 29,112 km of coastline of which 18,668 km is
Baltic. Columns this repository reads:

| column | meaning |
|---|---|
| `x_3035`, `y_3035` | sub-segment midpoint, EPSG:3035 |
| `flat_fraction` | length-weighted share of attributed evidence that is *flat* (dissipative, retentive) rather than *wall* (reflective, hard); `NaN` where nothing was attributed |
| `seg_len_m` | sub-segment length; the weight for every length statistic |
| `source` | `brisk`, `clms`, or `none` |
| `in_baltic` | domain flag |

Baltic coverage is 93.3 % attributed (BRISK 79.3 %, CLMS 13.9 %), mean
`flat_fraction` 0.527, and the distribution is strongly bimodal — 41.8 %
near-pure wall (≤0.05), 47.6 % near-pure flat (≥0.95), 2.5 % genuinely mixed
(0.35–0.65).

## Why this unblocks the term

[docs/beaching.md](../docs/beaching.md) switched `trap` off for a stated
reason, not a structural one:

> the only typing available is BSH's `H0 ≤ 0` tidal-flat flag, which is not a
> retentiveness proxy for Baltic shores — the basin is effectively tide-free,
> so the flag fires essentially only in the German Bight.

The sidecar table removes exactly that objection: it types the model's own
coastline from substrate evidence, Baltic-wide, at sub-cell resolution.

## Geometry: no reprojection, no match radius

The sidecar derives its outline from the **same** BSH H0 statics under the
**same** fine-over-coarse merge (`fine ∪ (coarse \ fine_bbox)`) that `024d`
rasterises, and ships midpoints in EPSG:3035 — the CRS `024d`'s raster is
already built in. The two coastlines are the same coastline, so mapping one
onto the other is a snap, not a nearest-feature match with a tolerance.

## Changes to `024d_BuildBeaching` (and `024e`)

`024e_BuildSurvivalOccupancy` carries its own copy of the same raster and rate
model, so every change below lands there too — minus the store schema, which
has no `shore_type` column. The two must stay on identical `trap` settings or
the survival store and the beaching store stop describing the same model; a
cell-for-cell comparison of the two `flat_fraction` planes is the check.


1. **Raster field.** `build_beaching_raster` gains a `flat_fraction` plane:
   snap every sub-segment midpoint to the 500 m raster; per touched cell take
   the `seg_len_m`-weighted mean over *attributed* segments; a second EDT over
   the touched cells propagates the value (and the propagation distance, as a
   diagnostic) to every raster cell.

2. **Drop the nearest-land hop.** `nearest_flat` was indexed through the EDT's
   `return_indices=True` land index. Nearest-*segment* is a better shore proxy
   and is what the new plane already carries, so `return_indices` goes away
   with it.

3. **Unattributed coastline is explicit, never silent.** A raster cell whose
   nearest sub-segment carries `NaN` takes the `ff_unattributed` parameter
   (default 0.5, the neutral midpoint and close to the observed 0.527 mean).
   `024d` prints the share of in-band steps and of beached weight that lands
   on unattributed coast, in the same spirit as the WAM extrapolation
   diagnostics.

4. **`trap` is linear in `flat_fraction`:**

   ```
   trap = trap_wall + (trap_flat - trap_wall) * flat_fraction
   ```

   `trap_flat` is the multiplier at pure-flat shore, `trap_wall` at pure-wall.
   Equal values keep the term inert, so the shipped defaults
   (`1.0`/`1.0`) reproduce the existing sweep bit-for-bit and the term is
   activated by moving them apart.

   Linear rather than binary: the sidecar applies no threshold and stores no
   `trap`, and a binary call is recoverable as `flat_fraction > 0.5` while the
   reverse is not. A threshold would also discard the 2.5 % genuinely mixed
   coastline into one arbitrary side.

   Not hard-zeroing `trap_wall`: an absorbing/reflecting dichotomy is a far
   stronger claim than a rate ratio, and it makes the beached total swing on
   the flat share. Keep it a free parameter with 0 available as one sweep
   member.

5. **`shore_type` in the store becomes real.** Same column, same `flat`/`wall`
   vocabulary, now `flat_fraction > 0.5` — the sidecar's own stated binary
   call, so no new threshold is introduced — plus `unattributed`, and `none`
   for the residual row. It stops being a label that expresses nothing, so
   `029` and [docs/visualisations.md](../docs/visualisations.md) lose the
   standing instruction not to report it.

## Data handoff

The sidecar gitignores its `data/processed/`, so the table is not fetchable
from the GitHub repository. It reaches this study the same way the derived
release-points geojson does: as a **slimmed derived blob in the data twin** at
`data/shoreclass_bsh_coastline/`, carrying only the five columns above, with
[`scripts/obtain/obtain_shoreclass.sh`](../scripts/obtain/obtain_shoreclass.sh)
as the canonical recipe (clone the sidecar, `pixi run pipeline`, slim, copy).

`ATTRIBUTION.md` gains one block: the derived product carries HELCOM BRISK and
Copernicus CLMS terms through.

## Open questions

- **`trap` and `s(w_onshore)` are probably collinear.** In the Baltic the wall
  shores are the exposed Fennoscandian ones and the flat shores are sheltered
  lagoons and bays, so the substrate factor and the wave-forcing factor may be
  measuring the same gradient twice. Measure the correlation between
  `flat_fraction` and `w_onshore` over in-band steps **before** reporting a run
  with both active.
- **Weight values are unconstrained.** There is no observational calibration
  for `trap_flat`/`trap_wall` here any more than for `w_half`. Sweep in `031`
  and report a range.
- **The kernel smooths over one cell face**, so on the 5.5 km coarse grid
  `flat_fraction` is a kilometre-scale quantity. Read the coarse-grid
  substrate signal as regional, not local — the same caveat `band_m` already
  carries.
