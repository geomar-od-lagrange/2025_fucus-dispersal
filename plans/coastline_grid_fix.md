# Coastline grid registration fix and sliver-hex dilution mask

Status: **plan written; implementation next** (update this line as steps land).

## Bug

`notebooks/004_extract_coastline.py` places wet cell `i` at
`lon[0] − dlon/2 + i·dlon` with `dlon = lon[1] − lon[0]` (same for lat) taken
from the H0 coordinates, which are stored rounded to 1e-4°:

| grid | step from first pair | true step | drift at the far edge |
|---|---|---|---|
| coarse lon | 0.0834 | 1/12 = 0.083333 | +0.33 cell (413 cells, 30° E) |
| coarse lat | −0.05 | −1/20 | none |
| fine lon | 0.0139 | 1/72 = 0.013889 | +0.5 cell |
| fine lat | −0.0083 | −1/120 = −0.008333 | −1.5 cell (southern fine rows) |

`data/bsh_hbmnoku_static/coastline{,_always_wet}.geojson` therefore drifts off
the model grid linearly with the cell index. Evidence and a corrected
water-area table: `tmp_outlier_hexes/` (`g_water_area.py`,
`out/water_area_fixed.parquet`, `out/h_area_*.csv`) — 200 of 2,383 partly-land
r6000m hexes change water area by more than 2×; the top dilution outliers
(2301, 1415, 2949, 9468, 1413, 4248, 9131, 16161) are such hexes.

The same index arithmetic sits in `024d`'s `_h0_nearest_sampler`
(`round((lon_q − lon[0]) / (lon[1] − lon[0]))`), which builds the 024d
water/land raster (distance to shore, seaward normal, hex labels) and
`bsh_water_on` for the Stokes extrapolation — same drift, same fix.

`003`'s `lon_f = lon + dlon/2` uses the rounded step too, but per element: a
constant offset of ≤ 3.4e-5° (< 4 m), not a drift. The c_file / 005
wall-field F-coordinates are correct to the 1e-4° storage rounding; left as is
(changing them would mean rerunning 003/005/010 for no measurable gain).

Out of this repo: the shoreclass sidecar
(`../2025_fucus-dispersal_shoreclass/scripts/02_extract_bsh_coastline.py`)
copies 004's arithmetic, so the sub-segment midpoints in
`data/shoreclass_bsh_coastline/` carry the same drift. Reported to the
authors; fixing it is a sidecar change plus a twin blob rebuild.

## Fix

- **004:** cell edges from a linear fit of the stored centres
  (`np.polyfit(index, coord, 1)`: slope = step, intercept = centre 0) — the
  regular grid the rounded coordinates sample. Residual ≤ 5e-5°, i.e. the
  storage rounding, not a drift. Edges = centre ∓ step/2, which are the NE
  F-points the c_files and 005 wall fields carry (T-cell `i` spans
  `F[i−1]..F[i]`). Vectorise the box build (`shapely.box` over
  `np.nonzero(wet)`). Regenerate both geojsons in `./data` and commit inside
  the submodule (no push); stage the gitlink deliberately.
- **024d:** `_h0_nearest_sampler` takes the step the same way.
- **Key:** 024a unchanged in code; rebuilt from the new geojsons.

## Sliver mask (dilution maps only)

`dilution` is NaN where `water_area_m2 < dilution_min_water_fraction ·
area_m2` (parameter, default 0.25, in each producing notebook): a hex with
less than a quarter of its area in water is left blank on dilution maps; it
keeps its `percentage`. Producers: 025, 026, 026a, 026b, 034 (all build
`dilution` in a notebook-local `to_hex_gdf`). The 032 headline embeds 025's
figures; nothing else divides by water area for a map.

**90 %-area metric (026 → 035, 034):** stays unmasked. It ranks hexes by
density but sums *water area*, so a sliver enters at its own (small, now
correct) water area and holds its own weight: the quantity is the smallest
water area holding the share, and a sliver's true area is the right cost.
Masking would either drop real particle-time from the denominator or charge
a hex for area that is land. The noise the mask removes is a display problem
(tiny denominators painting single hexes at the top of a log scale), not a
bias in a sum. Documented in `docs/visualisations.md`.

## Readers of the geojson / key fields

| reads | what | rerun |
|---|---|---|
| 024a | both geojsons → `water_area_m2`, `dist_to_coast_m`, and `helcom_subbasin` (shoreward snap gated on `water_area_m2 < area_m2`) | yes |
| 024c, 024f, 024g | key `helcom_subbasin` embedded in the stores | only if the rebuilt key changes any `helcom_subbasin` (checked in validation) |
| 024, 024b, 024e | key hex ids only (unchanged: the hex set comes from the H0 bbox) | no |
| 024d | own H0 raster (the 024d fix), shoreclass midpoints | yes, with the 024d fix → then 024e/f/g |
| 025, 026, 026a, 026b, 034 | key `water_area_m2` (dilution, 90 % area) + geojson overlay | yes |
| 027, 028, 029, 030, 031, 033, 036 | geojson overlay, key subbasin/geometry | yes (figures) |
| 035 | 026 area CSV, 029 curve CSV | yes, after 026/029 |
| 032 headline | 025 figures, 026/029 exports, key `mean_depth_m` | yes, last |

## Rerun order (sbatch, partition base)

1. 004 (sbatch wrap) → `data/bsh_hbmnoku_static/coastline{,_always_wet}.geojson`
   (done before review; validation below).
2. `sbatch scripts/024a_BuildHexKey_job.sh 6000` (r6000m is the only radius
   on disk).
3. If any `helcom_subbasin` changes: 024c for every (regime, year), and
   024f/024g (with 024e) for every member — via `submit_aggregates.sh` /
   `submit_024efg.sh`.
4. 024d (`scripts/024d_*_job.sh`, ~50 min on 141 nodes) → 024e/f/g
   (`scripts/submit_024efg.sh`) — needed for the 024d fix.
5. `sbatch scripts/025-031_ConsumersAll_job.sh` (default REGIMES), afterok on
   2–4.
6. 033 × {ALL, DJF, MAM, JJA, SON}, 036, then `scripts/headline_report_job.sh`.

## Validation (before the rerun)

- Overlay the new coastline on the 005 wall-field land mask and c_file
  wet/dry at 16.43E 57.09N, 22.8E 59N, 17.3E 62.3N (and one fine-grid
  south-east spot): edges coincide. Measure the max distance between
  polygon vertices and the nearest F-point lattice node.
- Hex water area from the new geojsons (024a's merge) vs
  `tmp_outlier_hexes/out/water_area_fixed.parquet` (`water_h0fix_km2`).
- `helcom_subbasin` and `dist_to_coast_m`: old vs rebuilt key (scratch
  output root), number of changed hexes.
- 024d raster: old vs fixed sampler, share of raster cells whose water/land
  flips.

## Docs

- `docs/visualisations.md`: sliver mask + 90 %-area rationale; coastline
  overlay paragraph.
- `docs/hexbinning_and_connectivity.md`: coastline provenance (grid fit), any
  counts quoted from the key that change after the rebuild.
- `docs/h0_semantics.md` / `docs/2d_field_extraction.md`: check the coastline
  references stay true.
