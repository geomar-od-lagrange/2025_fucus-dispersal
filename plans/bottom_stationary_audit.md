# Bottom-field stationarity audit

## Why this exists

The 021 TimeStats notebook showed bottom runs dropping more "land-seeded"
trajectories than surface runs (≈39.5% vs ≈34.9%). Deep dive on a single
fine-grid c_file (`2020-01-01T00:15`) found ~9,165 wet cells where the
deepest-non-NaN sigma layer carries `u=v=0.0` exactly, with valid tracers
and a smooth velocity profile ending in a wall value — consistent with an
explicit no-slip BC at the lowest sigma layer. In shallow columns even a
few layers above were bit-exact zero.

`021_TimeStats`'s land-seed heuristic (`dlon0==0 & dlat0==0`) then
registers particles that start in those zero-velocity bottom cells as
"land", even though they are clearly water.

Open question before we change `notebooks/003_prepare_2d_fields.py:load_bottom`:
are these always-immobile bottom cells a persistent structural feature of
the extracted 2D field, or a transient thing at t=0? The audit scans a
whole year of `output/2d_fields` to settle it.

## What the scripts do

`scripts/debug_bottom_stationary.py` walks every
`c_file_{res}_{year}*_{bottom,surface}.nc` in a directory and builds, per
`(lat, lon)`:

- `bottom_ever_moved` / `surface_ever_moved` — `True` if `u!=0 ∨ v!=0` at
  any timestep of any file in the year.
- `bottom_max_speed` / `surface_max_speed` — year-max of
  `sqrt(u²+v²)`, useful for spotting "nearly dead" cells.
- `bottom_dead_but_wet` — `surface_ever_moved ∧ ¬bottom_ever_moved`. The
  headline number. Output written to
  `output/debug/bottom_stationary/bottom_stationary_{res}_{year}.nc`.

`scripts/debug_bottom_stationary_job.sh` is the sbatch wrapper
(`ntasks=1`, `cpus-per-task=4`, `mem-per-cpu=8G`, 2 h). Usage:
`sbatch scripts/debug_bottom_stationary_job.sh [fine|coarse] [year]`,
defaults `fine 2019`.

## Picking up once results are back

1. Read `output/debug/bottom_stationary/bottom_stationary_fine_2019.nc`.
   Print `int(bottom_dead_but_wet.sum())`. Compare to the single-file
   number (5,631 in the smoke test on min_data).
2. Plot `bottom_dead_but_wet` on a map (a quick notebook in
   `notebooks/explore/` is fine). Inspect where the always-immobile
   cells sit — coast-hugging / shallow banks / specific subbasins.
3. Branch on the result:
   - **If the count is large and stable** (≫0, comparable to the
     single-file number): the no-slip wall is a structural feature. Fix
     `load_bottom` in `notebooks/003_prepare_2d_fields.py` to pick the
     deepest layer that is non-NaN **and** not `(0, 0)` — sketch was
     discussed in chat. Re-run 003 for one year to confirm the wet-but-
     dead count drops, then roll out to all years and re-run the
     bottom Parcels jobs and the 020–023 viz notebooks.
   - **If the count is small / sporadic**: stop. The single-file finding
     was a timestep-specific artefact, and `021_TimeStats`'s land-seed
     heuristic is the thing to tighten (use an explicit land mask from
     `land_mask_from_surface` rather than first-step displacement).
4. Either way, also verify: same scan on `coarse` — if the coarse field
   shows the same pattern, the fix needs to be applied identically to
   both.

## Files

- `scripts/debug_bottom_stationary.py` — scanner.
- `scripts/debug_bottom_stationary_job.sh` — sbatch wrapper.
- `notebooks/003_prepare_2d_fields.py:27-44` — `load_bottom`, candidate
  site for the fix.
- `notebooks/021_TimeStats.md` — land-seed heuristic consumer.
