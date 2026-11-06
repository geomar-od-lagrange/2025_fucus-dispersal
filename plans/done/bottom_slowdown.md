> Implemented: see [docs/visualisations.md](../../docs/visualisations.md) §034 and
> [docs/job_scripts.md](../../docs/job_scripts.md) (velocity-scaled variant regimes).
> Result (2019, 12 releases): 0.8x bottom velocity shrinks the 90 % footprint
> length sqrt(area) by 2-3 % at 30-180 d (`Exports/034/034_areas_2019.csv`).

# Bottom slowdown illustration: 0.8× bottom velocity

Replaces the deferred sweep in [analysis_wrapup.md](../analysis_wrapup.md) §7
with one illustrative variant (WR, 2026-10-20): the `bottom` regime at
**0.8×** the BSH bottom-cell velocity, **one release year (2019), one
release per month** (the first `doy = 1 + 5n` in each month: 1, 36, 61,
91, 121, 156, 186, 216, 246, 276, 306, 336). There is no drag model; the
claim is only "how does spread respond to an assumed 20 % slowdown".

## Status

- **Step 1 done** (9d0aaac): `velocity_scale` in 010 via
  `Field.set_scaling_factor` on every nest level's U and V; output regime
  `<regime>_vs<scale>`; `VELOCITY_SCALE` env through the bottom job script and
  `submit_010.sh`; 024's zarr-stem regex accepts any regime word. Checked in
  the production container (Parcels 3.0.6), 2-day / 2-per-cell test: 1.0 is
  bitwise identical to the pre-change notebook; 0.8 gives a first-hour
  displacement ratio of 0.800 (median; 0.8007 fine nest, 0.8000 coarse).
- **Step 2 submitted** (2026-10-20): 010 job **24456827**
  (`010-bottom-vs0.8-2019-n0_7_12_18_24_31_37_43_49_55_61_67`, 12 tasks,
  started at once; a 1x bottom release took ~7 h) and 024 counts
  **24456828** (`afterok:24456827`, writes
  `HexAgg_counts_r6000m_bottom_vs0p8_2019.parquet`).
- **Step 3 done**: `notebooks/034_BottomSlowdown.py`, dry-run green with
  `regime_var = bottom` against symlinked production inputs.
- **Remaining (step 4):** when 24456828 lands, check 12 zarrs under
  `Trajectories/bottom_vs0p8/2019/`, then run 034 —
  `pixi run papermill --cwd notebooks/ notebooks/034_BottomSlowdown.ipynb
  notebooks_executed/Visualisations/034_BottomSlowdown_2019.ipynb -p
  output_root /gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs -k
  python` (seconds, no Dask) — and move this plan to `plans/done/` (docs
  entries are in `docs/visualisations.md` and `docs/job_scripts.md`).

## Design

- **Scaling, not a kernel change.** 010 gets a primitive parameter
  `velocity_scale` (default 1.0) applied with Parcels'
  `Field.set_scaling_factor` to U and V of both nest levels; the kernel
  stays plain `AdvectionRK4`. A uniform factor in a time-varying flow is
  not a slowed-down replay of the 1× trajectories, so it needs runs.
- **Own regime name** `bottom_vs0p8` (output under
  `Trajectories/bottom_vs0p8/2019/`), so every downstream stage keys on it
  without new plumbing. Only the 024 counts store is built for it — no
  024b/c/d (no beaching for bottom).
- **Comparison** in a new notebook `034_BottomSlowdown`: 1× `bottom` vs
  `bottom_vs0p8` on the **same 12 release doys of 2019**, area holding
  50 % / 90 % of particle-time and `sqrt(area)` at 30 / 60 / 90 / 180 d —
  the R1 metric of [horizons_and_weighting.md](horizons_and_weighting.md)
  D6, reimplemented notebook-locally. Plus one map pair at 90 d. Export
  CSV under `Exports/034/`.

## Steps

1. 010: `velocity_scale` parameter; job-script plumbing; fix the stale
   "bottom-slowdown physics live in the kernel" comment in
   `scripts/010_FucusDispersal_bottom_job.sh`.
2. Submit the 12 runs, then the 024 counts partition afterok.
3. 034 notebook, tested on the 1× side while the runs go.
4. Run 034 when the counts land; docs entry (`docs/visualisations.md`,
   `docs/job_scripts.md`); move this plan to `plans/done/`.
