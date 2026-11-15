# Beaching statistics for bottom particles

Status: **plan, author decisions in (2026-11-15), not started.** Touches
024d–g, the beaching consumers and their job scripts; start only after the
coastline-grid-fix rerun has landed
([coastline_grid_fix.md](coastline_grid_fix.md)), which is editing and
rerunning the same notebooks.

## Design

`bottom` and `bottom_ub49` drifters strand under **the same criterion as
`surface_stokes`** ([beaching.md](../docs/beaching.md)): the same rate model,
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
10–60 m fine, 50–350 m coarse; [wall_velocity.md](../docs/wall_velocity.md)),
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

## Open questions for the authors

1. **Band** after Step 0: 2000 m or the first raster ring (`band_m = 750`),
   for all regimes?
2. If the band changes: rebuild the wave- and trap-axis sweep members (031,
   beaching.md tables) at the new band too, or only production?
3. Finer raster (250 m, full 024d rebuild) — only if Step 0 shows the ring
   quantisation matters?
