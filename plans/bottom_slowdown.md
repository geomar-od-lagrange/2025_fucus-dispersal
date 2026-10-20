# Bottom slowdown illustration: 0.8× bottom velocity

Replaces the deferred sweep in [analysis_wrapup.md](analysis_wrapup.md) §7
with one illustrative variant (WR, 2026-10-20): the `bottom` regime at
**0.8×** the BSH bottom-cell velocity, **one release year (2019), one
release per month** (the first `doy = 1 + 5n` in each month: 1, 36, 61,
91, 121, 156, 186, 216, 246, 276, 306, 336). There is no drag model; the
claim is only "how does spread respond to an assumed 20 % slowdown".

## Status

- Not started.

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
