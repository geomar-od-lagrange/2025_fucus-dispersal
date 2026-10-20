# 180 d horizons, survival weighting, areas and afloat shares

Implements the cross-cutting decisions of
[manuscript_extension.md](manuscript_extension.md) §Phase B (WR,
2026-10-20) plus its items P3, P8 and P9, and folds the pending
[subbasin_rebuild.md](subbasin_rebuild.md) into the same rebuild.

## Status

- Step 1 done (bd23dc4). Verified on 2019-01 `surface_stokes` sidecars at
  180 d: attrs round-trip, 024f Σocc identity exact (also with a partial
  trailing bin at 15 d), 17 origin subbasins present. Deviations: no
  `target_hex = -1` rows occur in practice (no drifter loses its in-key
  position), so the identity holds without them; 030's local
  `occupancy_max_days` renamed and 024d/029 "viability window" comments
  reworded to keep the grep clean; 030's `drifting_fraction_curve` and
  per-bin table now sum over all origins and the sentinel (step 3).
- Step 2 done (cdf7b5c). Verified on production counts (10 years) and
  seven new-schema 024f partitions built into a scratch root (2023/2024
  JJA + 2025-05): 026 bottom JJA and surface_stokes ALL unweighted, 026
  / 025 / 026a / 026b weighted surface_stokes JJA, 026a / 026b
  unweighted. 2025 reported and dropped for JJA; window asserts pass.
  Peak RSS 16.2 GiB (026, ss ALL, 10 years; was 56 GiB at four).
  Deviations: figure stems put the member after `r<radius>m` and before
  the facets, exports after the season; GeoJSON value column `n_obs` →
  `weight`; counts-origin NaN folded into -1; new `area_quantiles_csv`
  parameter in 026; AGENTS.md stage list updated. Finding: the D6 area
  *shrinks* with T (surface_stokes ALL, q = 0.5: 16 084 km² at 30 d →
  4 435 km² at 180 d) — cumulative particle-time concentrates in trap
  hexes, so it measures concentration, not spread (noted in
  visualisations.md). Step 3: 032 table E reads
  `Exports/026/026_<regime>_<season>_<member|unweighted>_T<h>d.geojson`,
  and its 025/026 figure stems need the member tag.
- Step 3 done (754dfdd). Verified on new-schema 024e
  (production + inert), 024f and 024g partitions built into a scratch root
  (surface_stokes 2019 m06–m08 + 2020 m06, 180 d): 028 (JJA/ALL, both
  weightings), 029, 030 (JJA/ALL), 031 (two members), 033, 025/026/026a/026b
  weighted, 026 bottom unweighted, then 032 over all of them. Old-schema
  partition raises the rebuild error; 024f N = 1 363 321 equals 029's
  Σweight; 030 drifting fraction 170–180 d 0.462 ≈ 029 basin afloat at 180 d
  46.0 %; afloat decreases with T for every origin. Deviations: the D8
  four-decade floor is **not** applied to the 025/026 family's dilution maps
  (tried: weighted median hex 5.3 decades below the near-land peak, so half
  the basin hit the floor) — they keep the 1st-percentile floor, documented
  in visualisations.md; 025/026/026a/026b also gained the old-schema error
  and the rate-key comparison (step 2 only compared the window); 029/030/031
  keep `age_bin_days`/`disp_bin_km` parameters and assert them against the
  attrs (the smoke test passes 1 d); 028 raw figure stem `_Tall` →
  `_T<max_age_days>d`; 029 maps/curves/displacement and 031 maps/beached
  panel in % of released, 029 per-source-hex panel 0–100 %; 032 drops its
  `age_bin_days` parameter (read from attrs), resolves figures by exact stem
  only (the `<stem>_*` fallback could pick another member whose tag extends
  this one), Table E/F pick `<member>` or `unweighted` per regime from the
  024f partitions on disk; smoke test's 028 survconn run gains
  `-p max_age_days 1`. Out of scope (024 lock): the 024g comment at
  ~392–394. Folded in after (e9cb195): area metric q = 0.9 only — 026
  `area_quantiles_csv` and 034 `quantiles_csv` default "0.9", 032 table F
  filters the CSV to its `area_quantile` parameter.
- Step 5 stores submitted 2026-10-20 (step 1 reviewed: bit-identical to
  the old stores on the shared window): 024a `24456829`; 024c
  `24456830`–`24456853` (24, afterok 024a); 024efg `24456854` (2825 cells,
  292 tasks, afterok 024a). Consumers wait for steps 2–4. Don't edit
  `notebooks/024e/f/g*.ipynb` or switch branches while 024efg runs: each
  cell reads the notebook from disk when it starts.

## Design decisions

- **D1. One window name, recorded in the store.** `max_float_days` (024e),
  `occupancy_max_days` (024f) and `connectivity_max_days` (024g) become one
  parameter `max_age_days`, production 180, set once in
  `scripts/024efg_ReduceAll_job.sh`. Once no lifetime is applied it is the
  same concept in all three: a reduction window over the sidecar's `obs`
  axis. Each reducer writes `df.attrs` (`max_age_days`, `age_bin_days`,
  the rate-model parameters) into its parquet; every consumer asserts that
  all partitions it pools agree and that `max_age_days` ≥ its largest
  horizon. Reason: today no store records its window, so a partial rerun
  would mix 60 d and 180 d silently. Considered encoding the window in the
  filename; rejected, 032's member regex would absorb it.
- **D2. 024f gains an origin axis and an exact denominator.** Bincount key
  `(release_doy, age_bin, origin_subbasin, target_hex)` with
  `origin_subbasin` via 024g's `to_subidx` on the 024a key, and the
  `target_hex = −1` sentinel rows kept (as 024g does), so
  `Σocc` per age bin = released water-seeded drifters × 240 exactly.
  024f therefore runs after 024a. Considered `release_hex` as the origin
  axis; rejected, it recreates the O(10⁸)-row counts store.
- **D3. Weighted heatmaps via a `member` parameter.** 025, 026, 026a,
  026b take `member` like 028: `""` reads the 024 counts (`n_obs`), a tag
  reads 024f (`surv`); inside the notebook the value is a neutral
  `weight`. Figure and export stems carry the member or `unweighted`.
  Default stays `""` (the smoke test runs 025/026 before 024f). The
  production driver runs `surface_stokes` at the production member and
  `bottom` unweighted.
- **D4. One window for weighted and unweighted.** Consumers reading the
  uncapped counts store cut to `age < max_age_days` (a consumer parameter,
  default 180), so 025's pooled maps and 028's view 1 compare like with
  like.
- **D5. Horizons.** Default `"30,60,90,180"` in 026, 026a, 026b, 028, 029;
  cumulative `[0, T)`. 030 keeps **snapshot** bins `[T − 10, T)` at the
  same T: afloat-at-T is an instantaneous quantity, a cumulative survival
  map would be a residence map (which 026 already is, weighted).
- **D6. P3 lives in 026.** Per-file reduction keeps `release_year`.
  For each T, pooled and per year: rank hexes by density
  (`weight / water_area_m2`, descending; zero-water hexes first at zero
  area), cumulative share, area to reach q ∈ {0.5, 0.9} with the last hex
  taken fractionally (a 6 km hex is 93.5 km², so whole hexes quantise
  small-T areas), `length_km = sqrt(area_km2)`. Pooled from summed weights;
  min–max over years; years lacking the season's releases are dropped and
  reported (as 028 does). Export
  `Exports/026/026_areas_<regime>_<season>_<member|unweighted>.csv`
  (`horizon_days, q, scope, area_km2, length_km`). No figure yet.
- **D7. P8 lives in 029.** `afloat(T) = Σweight − Σweight[beach_hex ≥ 0 &
  beach_age_bin < T/10]` per `release_hex → helcom_subbasin`, as % of the
  origin's Σweight (= its released water-seeded drifters, asserted in
  024e). Pooled + per-year min–max; −1 origin reported separately. "Afloat"
  includes drifters that left the BSH domain — say so in the export
  header. Export `Exports/029/afloat_by_origin_<regime>_<season>_<member>.csv`.
- **D8. P9 units.** 029 maps and curves in % of released water-seeded
  drifters (Σweight of the pooled partitions); GeoJSON gains
  `pct_released`. 030 maps in % of released using the 024f sentinel N:
  `surv / (N × 240) × 100` for a snapshot bin = mean share of released
  drifters afloat in that hex during the bin.
- **D9. 032.** Table E reads T180; table C's top-target horizon becomes
  its own parameter `top_horizon_days = 60` instead of "first entry of a
  list"; new table F (P3 areas) and G (P8 afloat by origin); captions
  derive horizons and window from parameters and store attrs, nothing
  hardcoded.
- **D10. Resources.** 024efg at `--mem-per-cpu=8G` (16 GB per task; measured peak 11.1 GB) (arrays grow ~3× for
  024e, ~1.5× for 024f/g at 180 d).

Deferred: P2's two-regime side-by-side figure (needs a two-regime loader;
after the rebuild), P7 mean age, P5 grouping (open decisions).

## Implementation steps (one commit each, docs in the same commit)

1. **Reducers** — D1, D2, D10: 024e/f/g, their job scripts, 024efg, smoke
   test; `docs/beaching.md` / `docs/survival_occupancy.md` /
   `docs/job_scripts.md` parameter tables. Lifetime prose goes; the
   60 d-measured result tables in `docs/beaching.md` get a one-line
   "measured at 60 d, refreshed after the rebuild" marker until step 5.
2. **Heatmap consumers** — D3, D4, D5, D6: 025, 026, 026a, 026b, their job
   scripts and `025-031_ConsumersAll_job.sh` member plumbing;
   `docs/visualisations.md`.
3. **Other consumers** — D1 asserts, D5, D7, D8, D9: 028, 029, 030, 031,
   033, 032; `docs/visualisations.md`.
4. **Review + smoke test** — separate review agent over steps 1–3;
   `scripts/0-31-smoke-test.sh` green.
5. **Rebuild on SLURM** (one pass, replaces the node-local recipe in
   subbasin_rebuild.md):
   1. `024a` key (6000 m, snap 1.0).
   2. 24 × `024c` afterok 024a.
   3. `024efg` (all 14 + trap members, 180 d) afterok 024a.
   4. `025-031` consumers + `033` afterok 2 and 3.
   5. `032` headline afterok 4.
   Then refresh the result tables in `docs/beaching.md` and the paper
   FACTS from `Headline/HEADLINE.md`; move this plan and
   subbasin_rebuild.md to `plans/done/`.
