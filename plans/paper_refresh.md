# Paper refresh: fold the eps = 2 mm/s pipeline into `paper/main.tex`

`paper/main.tex` (physics part) and its 15 PNGs are from the eps = 0 runs
(≈45 `[STALE eps0]` markers). The eps = 2 mm/s rerun 2016–2025, the
`bottom_ub49` regime, the coastline/shoreclass fix and the 750 m beaching
band for every regime finished 2026-11-15 22:01. Production member
`M = step_wc0p1_ts3_tcinf_tf1_tw0_b750`. `O =
/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs`.

## Decisions (authors, 2026-11-16)

1. **Bottom regime: both.** Survival-weighted (member `M`) in the main text,
   as the pipeline runs it; unweighted in a supplement. Remove every claim
   that negatively buoyant material does not strand (bottom: 22.6 % stranded
   within 180 d vs 67.3 % surface_stokes, `bottom_ub49` 8.0 %; HEADLINE.md
   Table A).
2. **bottom_ub49: FACTS plus both figures.** A 034 figure in Results
   "Spread by buoyancy", and ub49 as third regime line in the 035 figure.
3. **Dilution top:** authors pick after a look (see Open for authors; the
   floor turns out to matter more than the top).
4. **2000 m vs 750 m band sensitivity:** not in the paper.
5. **(a) 034 form: both** — T = 90 d dilution map pair (bottom | ub49) and
   the √area-vs-horizon plot.
6. **(b) Supplement (unweighted bottom): all four** — 026 dilution horizons,
   026 median age, 026a origin overview T180d, 028 emission T30d/T180d.
   All four PNGs exist on disk from the 2026-11-15 runs (no rerun unless the
   dilution top or 028 floor changes).
7. **(c) 028 below-floor:** authors asked whether the range should follow
   the new data; measurement and proposal under Open for authors.
8. **(d) P2 combined surface | bottom figure: dropped.** Delete the
   `FIGURE SLOT` comment (main.tex ~912–916); one figure per regime stays.

## Zero-water hex remap (figure_review_r2 G6): not needed

Share of particle-time in zero-water hexes (`water_area_m2 == 0`, 024a key
joined to `Exports/026/026_*_T{30,180}d.geojson`): surface_stokes `M`
1.4e-8 (3–4 hexes), bottom `M` 7e-10 (1 hex), bottom unweighted 0 hexes.
Drop G6. The ¼-hex floor (`dilution_min_water_fraction = 0.25`) blanks far
more: hexes < ¼ water hold 16.0 % of surface_stokes T180d particle-time
(572 of 10 581 hexes; 10.9 % at T30d) and 3.1 % of bottom (363 hexes).
At 0.10 it would be 4.3 % / 0.8 %. Scripts: `tmp_dilution_top/zero_water.py`,
`tmp_dilution_top/sliver_floor.py`. → Open for authors.

## Figure map (paper stem ← source PNG under `O/Figures/`)

| paper stem | source | note |
|---|---|---|
| `methods_shoreline_type` | `033/ShorelineType_r6000m.png` | |
| `r4_where_stranded_ALL` | `029/WhereStranded_surface_stokes_r6000m_ALL_M.png` | member gains `_b750` |
| `r4_r1_stranded_area_ALL` | `035/AreaTimeSeries_ALL_r6000m.png` | rerun with ub49 (step 2) |
| `r1_surface_ALL_dilution` | `026/TimeHorizonMaps_surface_stokes_ALL_r6000m_M_dilution.png` | |
| `r1_bottom_ALL_dilution` | `026/TimeHorizonMaps_bottom_ALL_r6000m_M_dilution.png` | was unweighted |
| `r1_medianage_surface_ALL` | `026/MedianAge_surface_stokes_ALL_r6000m_M.png` | |
| `r1_medianage_bottom_ALL` | `026/MedianAge_bottom_ALL_r6000m_M.png` | was unweighted |
| `r1_friction_map_T90d` (new) | `034/BottomFriction_map_ALL_T90d_r6000m_dilution.png` | |
| `r1_friction_length` (new) | `034/BottomFriction_length_ALL_r6000m.png` | |
| `r2_overview_surface_T30d` | `026a/OriginSubbasinOverview_surface_stokes_ALL_r6000m_M_T30d_dilution.png` | |
| `r2_overview_surface_T180d` | `026a/…_surface_stokes_ALL_r6000m_M_T180d_dilution.png` | |
| `r2_overview_bottom_T180d` | `026a/…_bottom_ALL_r6000m_M_T180d_dilution.png` | was unweighted |
| `r3_emission_surface_ALL_T{30,180}d` | `028/SubbasinConnectivity_emission_surface_stokes_ALL_M_T{30,180}d.png` | |
| `r3_emission_bottom_ALL_T{30,180}d` | `028/…_emission_bottom_ALL_M_T{30,180}d.png` | now with stranded column |
| `r5_occupancy_surviving_ALL` | `030/SurvivalHeatmaps_surface_stokes_r6000m_ALL_M.png` | |
| `s_bottom_unweighted_dilution` (new) | `026/TimeHorizonMaps_bottom_ALL_r6000m_unweighted_dilution.png` | supplement |
| `s_medianage_bottom_unweighted` (new) | `026/MedianAge_bottom_ALL_r6000m_unweighted.png` | supplement |
| `s_overview_bottom_unweighted_T180d` (new) | `026a/OriginSubbasinOverview_bottom_ALL_r6000m_unweighted_T180d_dilution.png` | supplement |
| `s_emission_bottom_unweighted_T{30,180}d` (new) | `028/SubbasinConnectivity_emission_bottom_ALL_unweighted_T{30,180}d.png` | supplement |

Captions: bottom figures drop "unweighted" and gain the stranding-model
clause of the surface captions; colour-scale text follows the chosen top
(026/026a/034) and 028 floor. Each figure's `% Source:` comment gets the new
path and write time.

## Steps

1. **Authors pick** dilution range (floor + top) and 028 floor (Open for authors).
2. **Consumer reruns** (node-local papermill fan-out or `sbatch`):
   - 035 always: `scripts/035_AreaTimeSeries_job.sh` with
     `regimes_csv "surface_stokes,bottom,bottom_ub49"`, `tags_csv "M,M,M"`
     (edit the script; the 026 area CSVs and 029 stranded curves for ub49
     exist). Also set the notebook default `regimes_csv` to the three regimes.
   - If the range changes: set `dilution_floor` / `dilution_top` in 026, 026a, 034 (and the docs
     scale text in docs/visualisations.md), rerun
     `026_TimeHorizonMaps_job.sh {surface_stokes,bottom,bottom_ub49} ALL 6000 {M,""}`,
     `026a_…_job.sh` same sweep, `034_BottomFriction_job.sh`; then 035
     (reads 026 area CSVs; unaffected by the top, but run it last).
   - If the 028 floor changes: edit `EMISSION_BOUNDS` / `emission_labels` in
     028 and the floor in its docs; rerun
     `028_SubbasinConnectivityMatrix_job.sh` for surface_stokes and bottom ×
     {M, unweighted} (+ ub49 for completeness).
   - Then `032` HEADLINE.md if any table source moved.
3. **Copy PNGs** per the figure map into `paper/figures/`; remove stems no
   longer used.
4. **Refresh `[STALE eps0]` FACTS** (agent work, keyword level, source pointers):
   - HEADLINE.md Table A/B (stranded %, sweep): ~1282 block, ~746, ~777.
   - Table C (connectivity): ~1163, ~1168. Table D + 027: ~1038.
   - Table E (top-10 shares): ~970, ~1033. Table F (90 % areas): ~970,
     ~1010, ~1023. Table G (afloat by origin): ~1312, ~1337.
   - 026 ALL executed notebooks (`notebooks_executed/Visualisations/026_*_ALL_r6000m_M.ipynb`):
     ~420 (remaining/released), ~432 (out-of-scale), ~439, ~453, ~1051 (median age), ~981.
   - 026a ALL executed: ~432, ~453, ~1098, ~1103, ~1109.
   - 028 ALL executed + `Exports/028/*_emission_fraction.csv`: ~1187, ~1200,
     ~1217, ~1227, ~1240 (closure, stranded share now also for bottom).
   - 029 / 030 ALL executed + `Exports/029/stranded_curve_*`: ~800, ~831,
     ~844, ~864, ~1359.
   - 031 ALL executed: ~777. 033 ALL executed: ~746, ~884.
   - 024d sidecars / docs/beaching.md: ~676, ~694, ~722; docs/wam_extrapolation.md: ~734.
   - 036 executed (land-seed share): ~108, ~114, ~149.
   - plans/numerical_coastal_trapping.md + HEADLINE.md: ~503 (store description).
   - plans/figure_review_r1.md Q4 diagnostic (rerun): ~1076.
   - Header comment (lines 12–31) and the Results NOTE (~895–903): rewrite
     to "current as of" the new commit; drop "held until" items now decided.
5. **Mark prose `[STALE]`** (humans rewrite; agent adds the marker + pointer):
   - Bottom-does-not-strand claims: ~486 (024g "surface_stokes only"), ~656
     (024d "surface_stokes only … does not strand"), ~929 caption
     ("unweighted"), ~1048 ("bottom is unweighted … does not strand … by
     design"), ~1205 ("Bottom matrices unweighted, no stranded column").
     Source: docs/beaching.md, plans/done/bottom_beaching.md.
   - C-grid-stall mechanism: ~981 ("stalling against the shore"), ~1001
     (90 %-only rationale citing C-grid stalls). Recheck on the new 026
     area CSVs whether the surface 90 % area still falls.
   - §3.4 bottom-friction `[STATUS]` block ~320: friction is implemented and
     run (`bottom_ub49`, U_b = 49 mm/s; docs/bottom_friction.md); replace
     with FACTS + pointer.
   - ~1373 "Supplementary figures: none" → the five supplement figures.
6. **Review agent** over the main.tex diff (facts vs sources), then hand to
   authors; commit/push in `paper/` only when told.

## Open for authors

- **Dilution top (0.1 vs 0.01).** Renders (ALL, member `M`, 026 code
  unchanged, production `Figures/` untouched):
  `tmp_dilution_top/{top0p1,top0p01,top0p001,floor1e-9_top1e-3}/Figures/026/TimeHorizonMaps_{surface_stokes,bottom}_ALL_r6000m_M_dilution.png`.
  **Neither top is in play:** 0 hexes reach the over colour at 0.1 or at
  0.01 (every horizon, both regimes). Max dilution is 1.3–1.9e-4 per km²
  surface, 0.7–1.6e-4 bottom (q99.9 ≤ 6.3e-5); > 1e-4: ≤ 0.06 % of hexes.
  At 0.1 the top six classes are empty, at 0.01 the top four.
  **The floor is the live problem.** Hexes inside the map extent below 1e-7
  (under colour, the same dark turbo end as the lowest class): surface 7 %
  (T30) / 17 % (T180), bottom 32 % / 32 % (unweighted bottom similar);
  they hold 0.2–1.7 % of the particle-time. Below 1e-8: surface 2 / 0.5 %,
  bottom 18 / 15 %; below 1e-9: surface 0.3 / 0.2 %, bottom 10 / 6 %. At
  1e-7 the Bothnian Bay and Bornholm–Gotland interiors read as one dark
  block. Scripts `tmp_dilution_top/{dil_max,dil_low,dil_extent}.py`.
  Options: (i) 1e-7–1e-3 (8 classes, `top0p001`): no over, under share as
  now; (ii) **1e-9–1e-3 (12 classes, `floor1e-9_top1e-3`, agent's
  suggestion):** no over, under ≤ 0.3 % surface / 6–10 % bottom, colour
  range spent on the data, basin interiors resolved; (iii) keep 1e-7–1e-1
  and state the under share in captions. Applies to 026, 026a, 034 (shared
  scale).
- **028 emission floor.** Non-zero cells below 1e-4 (ALL, from
  `Exports/028/*_emission_fraction.csv`): surface_stokes `M` 24 % (T30) /
  31 % (T180); bottom `M` 11 / 20 %; bottom unweighted 13 / 22 %; ub49 `M`
  10 / 13 %. Medians 1e-3–3e-2, q25 3e-5–2e-3. Pooled over these 8
  matrices: 21 % below 1e-4, 14 % below 1e-5, 8 % below 1e-6; the rest is a
  sparse tail to 1e-17, the survival-weighted exp(−A) residue that no
  floor should chase. **Proposal: floor 1e-6, 1-2-5 per decade (18
  classes, top 1).** Under-floor share drops to 7 / 15 % surface, 5 / 8 %
  bottom, ≤ 2 % ub49, and the under colour then marks a real "negligible"
  class instead of a fifth of the matrix. Cost: 18 batlow classes, adjacent
  ones closer in hue. Side-by-side (1e-4 | 1e-6, rough imitation of 028 View
  2, not the production layout): `tmp_028_range/emission_*_floor1e-4_vs_1e-6.png`
  (script `tmp_028_range/compare.py`, numbers `tmp_028_range/distribution.py`).
- **¼-hex water floor blanks 16 % of the surface particle-time at T180d**
  (3 % bottom); 0.10 would blank 4 % / 1 %. Keep 0.25 or lower it?
- **Bottom where-stranded map.** With bottom stranding 22.6 %, 029 writes
  `WhereStranded_bottom_r6000m_ALL_M.png` (and ub49); the paper has the
  surface one only. Add a bottom panel/figure, or FACTS only?
