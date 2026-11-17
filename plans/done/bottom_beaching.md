> Implemented — see [docs/beaching.md](../../docs/beaching.md).

# Beaching statistics for bottom particles

Status: **implemented; production build finished 2026-11-15 (below).** Author
decision 2026-11-15: `band_m = 750` is production for every regime (tag
`…_tf1_tw0_b750`), every sweep member rebuilt at 750, the production point at
2000 m kept as the one band sensitivity. Touches
024d–g, the beaching consumers and their job scripts; start only after the
coastline-grid-fix rerun has landed
([coastline_grid_fix.md](coastline_grid_fix.md), landed 2026-11-15), which edited and
reran the same notebooks.

## Design

`bottom` and `bottom_ub49` drifters strand under **the same criterion as
`surface_stokes`** ([beaching.md](../../docs/beaching.md)): the same rate model,
the same thresholded quantity — the surface onshore Stokes `w_on` sampled at
the drifter position — and the same absolute `w_c = 0.10` m/s. A narrow
near-shore band is what makes the surface/bed distinction irrelevant: within
the first raster ring of the coast the water is shallow against the wave
length, so the orbital motion reaches the bed. The Daily et al.
double-counting objection holds a fortiori: the bottom runs carry no Stokes
in their drift at all.

Considered, rejected: scaling `w_on` by the linear-wave bed factor
`F = sech(2kh)` (k from a CMEMS mean period, h = H0) — the authors chose the
same criterion for both regimes; the narrow band stands in for shallowness.
This also drops any new wave-variable download.

## Step 0: band width (first, on existing sidecars)

The production `band_m = 2000` is questioned: under strong onshore Stokes
surface drifters reach the wall anyway (ε rest point 0.01–0.07 cell ≈
10–60 m fine, 50–350 m coarse; [wall_velocity.md](../../docs/wall_velocity.md)),
so the first raster ring may be the right band, for **both** surface and
bottom. Approved test, on the `surface_stokes` sidecars written by the
coastline rerun, no 024d rebuild:

1. **Tag fix first.** `band_m` is in the parquet attrs but not in the member
   tag, so a band member would overwrite production. 024e/f/g append
   `_b<band_m>` (`.` → `p`), absent only at `band_m = 2000` — anchored to the
   value, like the trap suffix at the inert weights, so existing stores keep
   their names. Consumers take the tag as an opaque string, no change there.
2. 024e and 024f at the production point (`…_tf1_tw0`) with `band_m = 750`,
   all release years, beside the existing 2000 m stores: stranded % (031, 032
   table A), afloat-by-age curve (029 `BeachedFractionCurve`), stranding map
   and receiving-hex count (029/031), surviving occupancy (030).
3. Expectation if the argument holds: a **small** change for surface. For
   bottom (after Stage 1 below) a **much larger** reduction is expected,
   since bottom drifters are not pushed onto the wall and spread across the
   band; Validation measures it on one release before the full build.

**Result (2026-11-15, done).** Tag fix `86375ec`; the 2000 m tag is unchanged
and a 2019-01 024e rerun at 2000 m reproduces the production parquet exactly
(351 672 rows, only `git_sha`/`slurm_job_id` differ). `…_tf1_tw0_b750` built
for 024e + 024f, all 113 `(year, month)` cells (job 24930464; rows in
`024efg_ReduceAll_job.sh`). `surface_stokes`, 2016–2024 pooled (37.46 M
drifters), 180 d window; per-year min–max in brackets. Scripts, CSVs and
figures: `tmp_bottom_beaching/` (`compare.py`, `step0_*.csv`, `fig/`).

| quantity | 2000 m | 750 m |
|---|---|---|
| stranded ≤ 30 d | 38.4 % [36.1–40.2] | 35.4 % [33.1–37.3] |
| stranded ≤ 90 d | 57.7 % [53.1–61.8] | 54.4 % [49.6–58.4] |
| stranded ≤ 180 d | 70.3 % [65.2–76.0] | 67.3 % [62.0–73.5] |
| ≤ 180 d DJF / MAM / JJA / SON | 67.2 / 73.0 / 72.8 / 68.0 % | 63.9 / 70.1 / 69.9 / 65.0 % |
| receiving hexes (> 0) | 1 206 | 1 097 |
| hexes holding 90 / 99 % of stranded weight | 488 / 831 | 449 / 757 |
| Gini over receiving hexes; top-10 share | 0.687; 12.8 % | 0.687; 13.4 % |
| survocc `Σsurv / Σocc` (afloat share of particle-time) | 46.8 % | 49.9 % |
| in-band particle-hours (rel. to 2000 m) | 1 | 0.67 |
| strong-wave share of in-band hours (`w_on ≥ 0.10`) | 6.05 % | 7.93 % |
| strong-wave in-band hours (rel. to 2000 m) | 1 | 0.88 |
| drifters ever in band / ever strong-wave in band | 100.0 / 99.2 % | 99.98 / 98.1 % |

The 109 hexes lost at 750 m held 0.6 % of the 2000 m stranded weight;
Spearman rank correlation of per-hex stranded share 0.955; the top-10 hexes
are the same set but one. Per origin subbasin the 180 d drop is 0.1–5.2
points (largest in The Quark 25.1 → 19.9 % and Bothnian Sea 24.6 → 20.6 %,
smallest where stranding is near-total, Kiel Bay / Gdansk Basin ≥ 99.6 %).
The exposure numbers are counted directly from all 686 sidecars
(`exposure_one.py`, 657 of them in 2016–2024), with 024e's gate.

**Read:** as expected, surface stranding moves little — −3.0 points at 180 d,
the same shift in every season and well inside the interannual range; the
band keeps 88 % of strong-wave in-band hours while cutting in-band hours by
a third, so the drifters that strand sit in the first ring already.

**Raster quantisation.** `dist` is the EDT of the 500 m raster
(`raster_dx_m`) in whole cells, and the reducers gate `dist·25 < band_m`, so
`band_m = 500` selects *nothing*; 501–700 m takes the edge-adjacent ring,
701–1000 m adds the diagonals (707 m → 700 stored). The "500 m band" is
therefore `band_m = 750`: water within ~0–500 m of the raster land edge,
which itself sits within ±250 m of the model coast. A finer raster (250 m) is
a full 024d rebuild per regime (~50 min on the ~140-node allocation, 4× the
raster arrays, sidecar size unchanged) and only sharpens the fine nest — the
coarse model coast is 5 km.

**If the band changes:** every production member is rebuilt at the new band
for all regimes, the production tag gains `_b750` (member defaults in
consumers, job scripts and docs follow), and beaching.md states `w_c` as an
**absolute** threshold — today it is introduced as "the p93 of in-band
forcing", a 2 km-band statistic — with its Production/Sensitivity tables
regenerated.

## Pipeline changes

- **024d / 024e / 024f / 024g**: `_ZARR_STEM_RE` accepts
  `bottom_ub\d+(p\d+)?` — today `(surface_stokes|surface|bottom)_dt` rejects
  `Fucus_BSH_…_bottom_ub49_dt60min_…`. 024d otherwise unchanged: it samples
  raw `baltic_highres` Stokes at any drifter position, and `regime` is a
  parameter. Reducers: the `_b<band>` suffix (Step 0) only.
- **Scripts**: `024d_BuildBeachingForcing_job.sh` and
  `024efg_ReduceAll_job.sh` take `regime` as an argument instead of pinning
  `surface_stokes`; `submit_024d.sh` / `submit_024efg.sh` pass it.
  `submit_downstream_chain.sh` stages B/C loop over the chain regimes instead
  of skipping non-`surface_stokes` (B afterany on that regime's 010 jobs).
  `025-031_ConsumersAll_job.sh` gates 029–031 per regime and builds the
  `HexAgg_beaching_…_<regime>_…` member glob per regime.
- **Consumers** (core + weighted, per regime `surface_stokes`, `bottom`,
  `bottom_ub49`):
  - 029, 030, 031: `regime` is a parameter; run per regime.
  - 032: table A (stranded % per year and season) per regime.
  - 035: `stranded_regime` (one) → one stranded curve per regime of
    `regimes_csv`, drawn like its area panels.
  - 025, 026, 026a, 026b (024f member) and 028 (024g member) for the bottom
    regimes: the survival-weighted maps and matrices.
  - 034: a stranded-% row `bottom` vs `bottom_ub49` on the shared releases,
    from the two 024e production stores.
- **Docs**: beaching.md (regimes, per-regime results, absolute `w_c`),
  survival_occupancy.md, bottom_friction.md and job_scripts.md (drop
  "beaching is `surface_stokes`-only"), visualisations.md (035, 034),
  AGENTS.md bullets for 024d–g, 034, 035.

## Cost and disk

- 024d: bottom zarrs 2.8 GB/release vs 1.8 GB (`surface_stokes`), same
  87 200 × 5280 shape; ~1.5–2× the surface per-zarr cost → ~2–3 h for both
  regimes (~1460 zarrs) on the ~140-node allocation; memory stays inside the
  booked 18 GB/task (surface peak 11 GB, dominated by the same lon/lat arrays).
- Disk: ~93 MB/sidecar (2019 `surface_stokes`) → ~70–80 GB per regime,
  ~150 GB for both; check the gxfs quota first.
- 024e/f/g: the surface member list, seconds per zarr per member, ×2 regimes;
  consumers as for `surface_stokes`.

## Validation (`tmp_bottom_beaching/`, release 2019-10-18)

1. 024d (with the regex fix) on `Fucus_BSH_20191018T001500_{bottom,bottom_ub49}_*`
   into a tmp `output_root`; the `surface_stokes` sidecar of the same release
   from production.
2. Per regime: in-band particle-hours at 750 / 2000 m, share of them with
   `w_on ≥ w_c`, stranded % and median stranding age at the production member
   for both bands, receiving hexes.
3. Sanity: the 750 m surface member for this release equals the Step 0
   partition; at 2000 m the tmp reducers reproduce the production
   `surface_stokes` parquet for this release exactly (the tag change must not
   touch numbers).
4. Gate: bring the bottom/surface ratios to the authors before the full
   bottom build.

## Validation result (2026-11-15, job 24930567)

`tmp_bottom_beaching/val_job.sh` + `val_analyse.py`, release 2019-10-18,
tmp `output_root` (`tmp_bottom_beaching/val/`), production rate model, 180 d.

- New 024d on `surface_stokes` reproduces the production sidecar **bit for
  bit** (all five arrays + coords; attrs differ only in `git_sha`,
  `slurm_job_id`, `slurm_step_id`). 024e at 2000 m and at 750 m reproduces
  the production / Step 0 partitions' doy-291 rows exactly (63 882 and
  53 451 rows, `DataFrame.equals`).
- `bottom` and `bottom_ub49` sidecars build (~6 min/zarr wall, MaxRSS
  9.8 GB, 121 / 78 MB on disk); 024e at 750 + 2000, 024f + 024g at 750 run
  for all three regimes.

| regime | band | in-band h / drifter | strong share of in-band | ever strong-wave in band | stranded 30 / 90 / 180 d | median age bin | receiving hexes |
|---|---|---|---|---|---|---|---|
| surface_stokes | 750 | 2 405 | 9.4 % | 99.5 % | 35.7 / 59.0 / 64.1 % | 20–30 d | 734 |
| surface_stokes | 2000 | 3 168 | 8.1 % | 99.7 % | 40.0 / 63.1 / 66.9 % | 10–20 d | 827 |
| bottom | 750 | 460 | 2.2 % | 53.1 % | 7.2 / 21.4 / 28.9 % | 50–60 d | 559 |
| bottom | 2000 | 1 182 | 3.2 % | 77.0 % | 17.3 / 35.5 / 43.0 % | 40–50 d | 729 |
| bottom_ub49 | 750 | 74 | 2.9 % | 16.3 % | 2.5 / 6.0 / 9.7 % | 50–60 d | 428 |
| bottom_ub49 | 2000 | 512 | 3.6 % | 44.0 % | 11.4 / 19.4 / 24.8 % | 40–50 d | 640 |

**Read:** as expected the band cuts bottom stranding far more than surface
(180 d: −2.9 points surface, −14.1 bottom, −15.1 bottom_ub49): bottom
drifters spend 1/5 of the surface in-band time in the first ring, and see
strong waves there a quarter as often. Bottom/surface ratio at 180 d, 750 m:
0.45 (`bottom`), 0.15 (`bottom_ub49`). Gate item 4: bring these to the
authors before the full bottom build.

## Open questions for the authors (answered 2026-11-15)

1. **Band** after Step 0: 2000 m or the first raster ring (`band_m = 750`),
   for all regimes?
2. If the band changes: rebuild the wave- and trap-axis sweep members (031,
   beaching.md tables) at the new band too, or only production?
3. Finer raster (250 m, full 024d rebuild) — only if Step 0 shows the ring
   quantisation matters?

## Production build (launched 2026-11-15 19:48)

Review passed (one data step: stale 2000 m surface sweep stores moved to
`<outputs>/archive_b2000_sweep_20261115/`). Jobs: 024d bottom 24930591,
bottom_ub49 24930592 → 024efg surface_stokes 24930593, bottom 24930594,
bottom_ub49 24930595 → ConsumersAll (160 units) 24930596; 033 ×5
24930597–601; 036 24930602 → headline (3 regimes) 24930603. Log:
`<outputs>/logs/launch_bottom_beaching_20261115.log`. Next: check all OK
(consumer log "160 OK").

Build finished 2026-11-15 22:02: 024d 686/686 OK per bottom regime, 024efg
OK, consumers 160 OK / 0 FAIL, 033 ×5, 036, headline completed. Stranded
within 180 d (production `_b750`, 2016–2024 pooled, ALL seasons):
surface_stokes 67.3 %, bottom 22.6 %.

