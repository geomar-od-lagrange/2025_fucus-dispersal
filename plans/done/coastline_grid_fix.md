# Coastline grid registration fix and sliver-hex dilution mask

Implemented and rerun (2026-11-15); the current state is documented in
[hexbinning_and_connectivity.md](../../docs/hexbinning_and_connectivity.md)
(source data, coastline merge) and
[visualisations.md](../../docs/visualisations.md) (sliver mask, 90 %-area
metric, coastline overlay, 033 seam).

Status: **done.** Main-repo commits 9ff8f15 … (this archive); data twin
98022cd (coastline geojsons) and 247bcb3 (shoreclass table), local, not
pushed; sidecar `2025_fucus-dispersal_shoreclass` 6d847d0, 2384256, local,
not pushed.

Rerun job ids: 024a 24930310 → 024d 24930311 → 024efg 24930312 →
consumers 24930313 (afterok 024a + 024efg); 033 ALL/DJF/MAM/JJA/SON
24930314–24930318 (afterok 024efg); 036 24930319 (afterok 024a); headline
24930320 (afterok consumers, 033, 036). Pre-rerun snapshot (Headline/, 026/028/
029/033/034 CSVs, old key): `tmp_coastline_fix/before/`.

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

The shoreclass sidecar
(`../2025_fucus-dispersal_shoreclass/scripts/02_extract_bsh_coastline.py`)
copied 004's arithmetic, so the sub-segment midpoints in
`data/shoreclass_bsh_coastline/` carried the same drift. Fixed there
(`grid.axis_fit`, the same fit; `h0_bbox` on the fitted edge; stage 4 tells
the nest seam from genuine coast by which side of the box edge the water is
on), stages 2/4/5/6/7 rerun, twin blob rebuilt with
`scripts/obtain/obtain_shoreclass.sh`. Old → new: 65,876 → 64,799
sub-segments, 29,112 → 28,649 km; Baltic attributed 94.2 → 97.1 % (sidecar
stats 93.3 → 97.1 %), mean `flat_fraction` 0.527 → 0.537, > 0.5 share
52.2 → 53.3 % of attributed; `on_seam` 472 km (open-water double ring) →
6 km. Per sub-segment vs the nearest old midpoint (median 281 m away):
length-weighted mean |Δff| 0.012, 1.1 % flip sides of 0.5.

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
- **Key:** 024a's `fine_bbox` from the same fitted outer edge
  (`grid_extent`). With the geojson on the true grid, the old
  `min ± half mean step` bbox sat 1e-5° off the fine polygons; the hairline
  gap between fine and clipped coarse water became a self-intersection in
  EPSG:3035 (TopologyException in the water-area intersection).

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

Done: 004 → `data/bsh_hbmnoku_static/coastline{,_always_wet}.geojson`.

1. `sbatch scripts/024a_BuildHexKey_job.sh 6000` (r6000m is the only radius
   on disk; ~23 min).
2. 024c: **not needed** — the rebuilt key changes no `helcom_subbasin`.
   024/024b: hex ids only, unchanged.
3. `scripts/submit_024d.sh --dependency=afterok:<1>` (all surface_stokes
   zarrs, `overwrite True`; ~50 min on 141 nodes last time) then
   `scripts/submit_024efg.sh --dependency=afterok:<3>` — for the 024d fix.
4. `sbatch --dependency=afterok:<1>:<3b> scripts/025-031_ConsumersAll_job.sh`
   (default REGIMES; 120 units, check "N OK, 0 FAIL").
5. afterok on 4: `sbatch scripts/033_ShorelineTypeMap_job.sh surface_stokes 6000
   step_wc0p1_ts3_tcinf_tf1_tw0 <S> step_wc0p1_ts3_tcinf` for S in ALL DJF MAM
   JJA SON; `sbatch scripts/036_SeedingMap_job.sh 6000`; then
   `sbatch scripts/headline_report_job.sh`.

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

## Validation results (2026-11-15, `tmp_coastline_fix/`)

- Vertex → F-lattice distance (cells, max): old fine 0.50 lon / 0.50 lat,
  coarse 0.33 / 0; new fine 0.004 / 0.006, coarse 0.001 / 0.
- Point-in-polygon at every H0 centre vs H0 notnull: old fine 98.61 %
  (3,398 cells differ), new 100 % both grids (also `always_wet` vs H0 > 0).
  vs the 005 surface_stokes land flag: new fine 99.61 % (949), coarse
  99.98 % (33) — the cells 005 closes that H0 calls wet, not registration.
  Figure: `fig/coastline_old_new_vs_walls.png` (16.43E 57.09N, 22.8E 59N,
  17.3E 62.3N, 14.3E 53.9N).
- Scratch key (`out_root/`, `out/key_compare.parquet`): water area vs
  `water_area_fixed.parquet` h0fix max |Δ| 0.035 km² over 12,420 wet hexes
  (the fine-bbox edge, 1e-5°); totals old 1,044,544, new 1,044,258 km².
  159 partly-land hexes change water area by > 2×; slivers (< ¼ water)
  781 → 765. `helcom_subbasin`, `mean_depth_m`, `fucus_area_m2` unchanged.
  `dist_to_coast_m`: median Δ 0.5 km, up to 121 km in the North Sea, where the
  old fine/coarse seam gap was counted as coast; no notebook reads it.
- All eight former outliers are now < 0.1 water fraction → masked.
- 024d 500 m raster: 46,610 of 4.18 M water cells (1.1 %) flip water/land
  between the old and fixed sampler (14,302 fine, 32,308 coarse) — all
  coastal, i.e. in the beaching band.

## Docs

- `docs/visualisations.md`: sliver mask + 90 %-area rationale; coastline
  overlay paragraph.
- `docs/hexbinning_and_connectivity.md`: coastline provenance (grid fit), any
  counts quoted from the key that change after the rebuild.
- `docs/h0_semantics.md` / `docs/2d_field_extraction.md`: check the coastline
  references stay true.

## Rerun outcome (2026-11-15)

| job | id | outcome |
|---|---|---|
| 024a r6000m | 24930310 | COMPLETED 29 min; identical to the validated scratch key |
| 024d | 24930311 | 686 OK, 0 FAIL (32 min, 42 nodes) |
| 024efg | 24930312 | 2825 OK, 0 FAIL (21 min) |
| 025–035 consumers | 24930313 | 120 OK, 0 FAIL (56 min) |
| 033 ALL/DJF/MAM/JJA/SON | 24930314–18 | COMPLETED, no errors |
| 036 | 24930319 | COMPLETED |
| 032 headline | 24930320 | COMPLETED |

- 024d raster: water cells 4,178,297 → 4,176,989; coastal land cells with a
  direct `flat_fraction` 52.7 % → 52.3 %; mean `ff` over the 5 km band
  50.0 % → 50.8 %.
- Headline (pooled 2016–2024, production member, ALL): stranded within 180 d
  67.6 → 70.3 % (seasons 64.5–69.9 → 67.2–73.0 %); afloat at 180 d, all
  origins 32.4 → 29.7 % (seasons 30.1–35.5 → 27.0–32.8 %); inert member
  97.1 → 98.0 %; survival share of unweighted particle-time at 180 d
  0.497 → 0.468. Bottom 90 %-area changes < 1 % (water areas only).
- Former outliers: all masked in the production 026 T90 exports; 0–90 d
  dilution maximum 0.017 → 1.1e-4 (bottom) and 0.017 → 2.6e-4
  (surface_stokes, unweighted); blank hexes 12 → 337 and 21 → 563.
  Figure: `tmp_outlier_hexes/fig/before_after_dilution_T90d.png`.
- Not refreshed: `docs/beaching.md`'s results tables (54.3 % production,
  93.0 % inert, sweep table) predate this rerun and the previous one; they
  need a separate pass from `031`/`032`.
