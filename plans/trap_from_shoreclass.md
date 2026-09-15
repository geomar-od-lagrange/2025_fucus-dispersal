# Activate `trap` from the shoreclass sub-segment classification

Re-establish `trap` in the beaching rate as a **continuous substrate factor**,
driven by the `flat_fraction` table produced by the sidecar repository
<https://github.com/geomar-od-lagrange/2025_fucus-dispersal_shoreclass>.

The classification is not built here: this repository consumes one table and
turns it into a rate multiplier.

## What the sidecar delivers

One row per sub-segment of the BSH tracer-cell outline — the staircase of cell
faces bounding the wet cells, subdivided at k=2 into 65,876 sub-segments of
420–465 m, 29,112 km of coastline, 18,668 km of it Baltic. The slimmed twin
blob carries six columns:

| column | meaning |
|---|---|
| `x_3035`, `y_3035` | sub-segment midpoint, EPSG:3035 |
| `flat_fraction` | length-weighted share of attributed evidence that is *flat* (dissipative, retentive) rather than *wall* (reflective, hard); `NaN` where nothing was attributed |
| `seg_len_m` | sub-segment length; the weight for every length statistic |
| `source` / `in_baltic` | provenance (`brisk`, `clms`, `none`) and domain flag |

Baltic coverage: 93.3 % attributed (BRISK 79.3 %, CLMS 13.9 %), mean
`flat_fraction` 0.527, strongly bimodal — 41.8 % near-pure wall (≤0.05),
47.6 % near-pure flat (≥0.95), 2.5 % mixed (0.35–0.65).

## Why this unblocks the term

[docs/beaching.md](../docs/beaching.md) switched `trap` off for a stated
reason, not a structural one: the only typing available was BSH's `H0 ≤ 0`
tidal-flat flag, no retentiveness proxy for a tide-free basin. The sidecar
types the model's own coastline from substrate evidence, Baltic-wide — and
derives that outline from the **same** BSH H0 statics under the **same**
fine-over-coarse merge (`fine ∪ (coarse \ fine_bbox)`) that `024d` rasterises,
in EPSG:3035, the CRS of `024d`'s 500 m raster. Mapping one onto the other is
a snap, not a nearest-feature match with a tolerance.

## `024d_BuildBeachingForcing`: the `ff` sidecar variable

`build_beaching_raster` gains a flat-fraction plane, and the sidecar gains one
variable:

1. **Snap.** Round each midpoint `(x_3035, y_3035)` through the raster's own
   affine to `(row, col)` — no reprojection, the table is already in the
   raster's CRS. Should that change, reproject at read time.

2. **Per-cell value.** `seg_len_m`-weighted mean of `flat_fraction` over the
   midpoints landing in the cell. Sub-segments the sidecar left unattributed
   (`NaN`) enter as **0.5** — the neutral midpoint, close to the observed
   0.527 mean — substituted once at build time. Deliberately **no**
   `ff_unattributed` parameter and **no** flag: a parameter invites sweeping a
   quantity carrying no evidence, and a flag propagates a third state through
   the reducer chain for 6.7 % of Baltic coastline.

3. **Propagate.** Cells no midpoint reaches inherit from the nearest cell that
   does; each position then reads **its nearest land cell**, through the same
   distance-EDT `return_indices` index that produces `dist`, so `dist` and
   `ff` always describe the same shore.

4. **Store.** `ff` is `uint8` **percent, 0–100**, on `(trajectory, obs)` like
   every other sidecar variable, and **0 outside the sampling band**,
   mirroring `flat`. Out-of-band steps carry no rate, so it is never read
   there.

5. **Currency.** The sidecar attrs gain the sha256 of the shoreclass parquet
   next to the git sha and SLURM job id, so a sidecar built against a
   superseded table is identifiable without re-deriving it. `024d` prints the
   direct-attribution share of coastal land and the in-band mean of `ff`.

## The reducers: `024e`, `024f`, `024g`

The reducers read `ff` and never touch the parquet. `trap` is linear in it:

```
trap = trap_wall + (trap_flat - trap_wall) * ff / 100
```

`trap_flat` is the multiplier at pure-flat shore, `trap_wall` at pure-wall;
equal values keep the term inert, so the shipped `1.0`/`1.0` reproduces the
existing sweep bit-for-bit and moving them apart activates the term.

Linear rather than binary: the sidecar applies no threshold, a binary call is
recoverable as `ff >= 50` while the reverse is not, and a threshold would
discard the 2.5 % mixed coastline into one arbitrary side. `trap_wall` is not
hard-zeroed — an absorbing/reflecting dichotomy is a far stronger claim than a
rate ratio and makes the beached total swing on the flat share. Keep it free,
with 0 available as one sweep member.

**Member tag.** The member suffix gains `_tf{trap_flat:g}_tw{trap_wall:g}`
**only when either weight differs from 1**, with `.` → `p`. An inert member's
filename is unchanged; only a substrate-sensitive one pays for the name.

**`shore_type` in the beaching store becomes real.** Same column, same
`flat`/`wall` vocabulary, now `ff >= 50` — the sidecar's own stated binary
call, no new threshold — plus `none` for the residual row. There is no
`unattributed` level: that coastline is folded into `ff` at 50 and resolves to
`wall`. The column stops expressing nothing, so `029` and
[docs/visualisations.md](../docs/visualisations.md) lose the standing
instruction not to report it. `024f`/`024g` carry the identical raster and
rate model but no `shore_type` column, and must stay on identical `trap`
settings or the survival and beaching stores describe different models.

## Data handoff

The sidecar gitignores its `data/processed/`, so the table is not fetchable
from its GitHub repository. It reaches this study as the derived
release-points geojson does — a **slimmed blob in the data twin** at
`data/shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet`, with
[`scripts/obtain/obtain_shoreclass.sh`](../scripts/obtain/obtain_shoreclass.sh)
as the canonical recipe (clone the sidecar, `pixi run pipeline`, slim, copy).
`ATTRIBUTION.md` carries the HELCOM BRISK and Copernicus CLMS terms through.
Sidecar stage 1 needs a CLMS service key (EU Login, *My settings → API
Tokens*) and a cut that queues 10–30 min server-side, so the recipe cannot run
unattended.

## Open questions

- **`trap` and `s(w_onshore)` are probably collinear.** The Baltic wall shores
  are the exposed Fennoscandian ones, the flat shores the sheltered lagoons
  and bays, so the two factors may measure one gradient twice. Correlate `ff`
  against `w_on` over in-band steps **before** reporting a run with both
  active.
- **Weight values are unconstrained** — no observational calibration for
  `trap_flat`/`trap_wall` any more than for `w_half`. Sweep in `031`, report a
  range.
- **The kernel smooths over one cell face**, so on the 5.5 km coarse grid
  `flat_fraction` is kilometre-scale. Read the coarse-grid substrate signal as
  regional, not local — the same caveat `band_m` carries.
