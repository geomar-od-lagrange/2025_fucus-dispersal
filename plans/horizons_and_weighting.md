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
