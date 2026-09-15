# Cluster job scripts and notebook output layout

`scripts/0NN_*_job.sh` runs notebook `0NN` on NESH; the rendered copy lands
under `notebooks_executed/Visualisations/`.

```
notebooks/0NN_*.md             ← jupytext source of truth
notebooks/0NN_*.ipynb          ← code-only, in-sync with .md, no outputs
notebooks_executed/
  Visualisations/
    0NN_*_<regime>_<year>.ipynb   ← papermill output, with figures
```

`notebooks/*.ipynb` are stripped before commit (`pixi run jupyter nbconvert
--clear-output --inplace` after any local execution). Cluster runs never
overwrite the source — papermill reads `notebooks/` and writes per-parameter
copies under `notebooks_executed/`, an artifact tree, not source.

## Standalone figure PNGs

The map notebooks write their figures as standalone PNGs under the outputs tree
(outside the repo), one subdir per stage, so they are usable without opening the
executed notebook:

```
output_root/Figures/
  025/  HexRelativeDensity[DE|ByOrigin]_<regime>_<season>_r<radius>m[_part<n>].png
  026/  TimeHorizonMaps_<regime>_<season>_r<radius>m_<percentage|dilution>.png
  026a/ OriginSubbasinTimeHorizonMaps_<regime>_<season>_r<radius>m_<subbasin>_<column>.png
  026b/ OriginSubbasinYearTimeHorizonMaps_<regime>_<season>_r<radius>m_<subbasin>_<year>_<column>.png
  027/  DistanceQuantiles_<regime>_r<radius>m_<season>.png
  028/  connectivity_<regime>_<season>_<member-or-unweighted>_*.png
  029/  <view>_<regime>_r<radius>m_<season>_<member>.png
  030/  SurvivalHeatmaps_<regime>_r<radius>m_<season>_<member>.png
  031/  BeachingSweep{Stats,Maps}_<regime>_r<radius>m_<season>.png
```

`figure_dir` is derived from the `output_root` parameter inside each notebook
(no extra job-script argument); submit once per (regime, season). Every figure
is saved print-ready — one full 180 mm text width at 300 dpi — so the PNG is
the manuscript figure, not a preview (rationale in
[visualisations.md](visualisations.md)).

`025`–`028` additionally write machine-readable exports under
`output_root/Exports/<stage>/`: GeoJSON (hex geometry, EPSG:4326) for the map
products of `025`/`026`/`026a`/`026b`/`027`, long-form CSV for `028`'s
connectivity matrices.

## The Parcels runs (010)

`scripts/010_FucusDispersal_{surface,bottom,surface_stokes}_job.sh` fan out one
single-process papermill run per release (`xargs -P ${SLURM_NTASKS}` over `srun
--ntasks=1 --exact` steps), each writing one trajectory zarr under
`output_root/Trajectories/<regime>/<year>/`. Two positional args:

| arg | default | meaning |
|---|---|---|
| `$1` | `2019` | release year |
| `$2` | `73` | number of releases to run, taken from the front of the schedule |

Releases sit at `doy = 1 + 5n`, `n ∈ [0, n_releases-1]` — 73 per full year, doys
1…361, the same regardless of leap year. Each run integrates `simulation_days=220`
forward from `<doy>T00:15:00`, so the *last* release of year Y needs forcing
through roughly Y+1-08-03. `n_releases` is how a year at the end of the forcing
record runs only the releases whose whole 220-day window is covered: with 2D
fields ending 2025-12-31, the 2025 releases that fit are `n ∈ [0, 28]` (doys
1…141, last window closing 2025-12-27), so `n_releases=29`.

Submit through the wrapper rather than `sbatch` directly:

```bash
scripts/submit_010.sh <regime> <year> [n_releases] [extra sbatch args...]

scripts/submit_010.sh bottom 2022              # full year, --ntasks=73
scripts/submit_010.sh surface_stokes 2025 29   # partial year, --ntasks=29
```

`--ntasks` has to equal the release count and an `#SBATCH` directive cannot read
the job script's own arguments, so the wrapper sets both from one number — a
29-release job books 29 tasks, not 73. It also expands
`scripts/node-blacklist.txt` into `--exclude=` (as `submit_024d.sh` does) and
passes further `sbatch` flags through.

At `--mem-per-cpu=60G` a 73-task job needs ~19–25 nodes, so these are the
largest allocations in the pipeline; submit a couple at a time and let SLURM
queue the rest. `--time=36:00:00` (the `base` cap is 48 h): observed wall times
run 5–19 h, and the fieldset the notebook builds grows with every year of 2D
fields on disk, because `010` globs the whole `2d_fields/` directory rather than
the run window.

### Forcing gaps are silent

`010` globs `2d_fields/`, intersects the fine and coarse stem sets and drops
whatever is missing from either, with a warning. The surviving stems are handed
to `FieldSet.from_netcdf` as an explicit `timestamps` list, so a hole is not a
missing time -- it is a longer interval between two samples, and Parcels
interpolates straight across it. Nothing raises, and
`allow_time_extrapolation=False` does not catch it either, because the hole is
interior to the time axis. Check `2d_fields/` for missing 6-hourly stems across
every year a run window touches *before* submitting, and record what you find
with the results.

Two holes are missing at source (no raw BSH c-file exists, so re-running `003`
cannot fill them) and are crossed by runs in this study:

| gap (model time) | length | crossed by |
|---|---|---|
| 2025-07-24T12:15 - 2025-07-25T06:00 | 18 h | all 2024 releases whose window reaches it, and 2025 releases |
| 2025-10-22T18:15 - 2025-10-23T00:00 | 6 h | 2025 releases |

A third missing stem, 2025-12-31T18, sits past the end of every window that can
be run (the last 2025 release closes 2025-12-27), so it costs nothing beyond
capping the forcing record at 2025-12-31T18:00.

Interpolating over 18 h of a 220-day trajectory is the only option available;
what matters is that it is a known, bounded, documented smoothing of the
velocity field on those dates and not a silent truncation. `2025-07-24/25` is
the one to keep in mind when reading 2024 and 2025 results.

## Multi-task dask layout

`--ntasks=N ≥ 3`, one role per SLURM task: task 0 the dask scheduler plus a
local worker set, task 1 papermill (connecting via `$SCHEDULER_FILE`), tasks
2…N-1 extra worker tasks.

Bootstrap is one file. Task 0 writes `$SCHEDULER_FILE`; the papermill task and
any extra workers poll for it (`Client(scheduler_file=os.environ["SCHEDULER_FILE"])`
in the notebook, `dask worker --scheduler-file ...` in workers). A `trap cleanup
EXIT` removes it on exit, and the hard-coded `sleep 30` after task 0 launches
lets the scheduler come up before papermill starts.

## Per-job constants

Every script sets:

```bash
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export SCHEDULER_FILE=${repo_root}/.scheduler_${SLURM_JOB_ID}.json
```

`output_root` is passed to papermill via `-p output_root ${output_root}` — the
only contract between job script and notebook; there is no `$FUCUS_OUTPUT_ROOT`
env variable. `http_proxy`/`https_proxy` are set for HTTP traffic from the env
(e.g. CMEMS in 002); `no_proxy` excludes intra-cluster TCP so dask node
addresses don't dial the proxy.

## The per-(regime, year) aggregations (024 / 024b / 024c)

`scripts/submit_aggregates.sh <regime> <year> [hex_radius] [sbatch args...]` is
the one call to make per regime-year as the 010 runs land:

```bash
scripts/submit_aggregates.sh bottom 2020
scripts/submit_aggregates.sh surface_stokes 2025 6000 --time=12:00:00
```

It submits `024_BuildHexAggregates` (counts) and `024b_BuildHexDistance`
together — independent Dask jobs over the same trajectory zarrs — and chains
`024c_BuildHexConnectivity` behind the counts job with
`--dependency=afterok:<024 job id>`, because 024c reduces what 024 wrote. It
refuses a (regime, year) with no trajectory zarrs, expands the node blacklist
into `--exclude=` and prints the three job ids. All three job scripts also take
`[regime] [year] [hex_radius]` directly, for a single rerun.

The beaching branch (024d and its reducers) is not part of this: it is one
whole-regime pass over every year at once, submitted separately.

## Horizontal scaling for the GPFS-bound builder (024d)

`scripts/024d_BuildBeachingForcing_job.sh` fans out one independent
single-process papermill run per trajectory zarr (`xargs -P ${SLURM_NTASKS}`, so
`njobs ≠ ntasks`: throttle concurrency without editing the job). No MPI, no
Dask; the bottleneck is streaming trajectory zarrs and hourly Stokes files off
GPFS. Three consequences, all defaults in the script:

- **`--distribution=cyclic`, not `--spread-job`.** Horizontal reach into the
  filesystem beats node locality: packing is a tail-latency disaster, with a
  48-cell month grid running steps min/med/max `11:54 / 17:55 / 54:51` on 8
  nodes against `10:51 / 16:43 / 20:54` on 22 — identical median, 2.6x worse
  tail, and wall time is set by the slowest cell. But `--spread-job` asks for
  one node *per task*, and a few-hundred-node request against ~94 idle nodes
  pends on `(Resources)` forever; `cyclic` round-robins the tasks over however
  many nodes the scheduler grants, which is the spread that matters.
- **No `--constraint`.** Pinning to sapphire (srp) nodes was inherited from
  Dask-backed jobs where IB reliability mattered; here it only shrinks the
  eligible pool and leaves jobs pending indefinitely.
- **`--ntasks` sized from the work list** by `scripts/submit_024d.sh` (one task
  per trajectory zarr, capped at `MAX_TASKS`, default 292), and memory sized
  from measured `sacct MaxRSS` (11.02 GiB peak at `window_days=220`, plus 55%
  headroom: 9G/cpu x 2 = 18 GB).
  Over-requesting either forces the job onto more nodes than it needs and
  pushes it into `(Resources)`.

**Which years run.** No job script names a release year. 024d's work list is
every trajectory zarr under `output_root/Trajectories/<regime>/`; the work list
of 024e/024f/024g and of `024efg_ReduceAll` is every distinct `(release year,
release month)` that has a sidecar under
`output_root/BeachingForcing/<regime>/`. A year therefore joins the grid as
soon as its upstream stage lands, and a partial release year — one whose 010
runs stop short of 73 releases — contributes only the months it has instead of
failing on the empty ones. `YEARS="2020 2021"` in the environment restricts the
set; unset takes everything on disk. `LIST_ONLY=1 bash <job script>` prints the
work list and exits, which is both how the submit wrappers size `--ntasks` and
how to check a selection before submitting.

**Node failures.** A single dying node aborts the whole allocation unless
`--no-kill` is set (job 23834688 lost `nesh-clk542`, "Node unexpectedly
rebooted", at 53/292 done), and the resulting per-step "CANCELLED ... DUE TO
NODE FAILURE" lines name almost every allocated node, so the culprit has to be
read off `sinfo -t down,drain` instead. Culprits go in
`scripts/node-blacklist.txt`, which `scripts/submit_024d.sh` expands into
`sbatch --exclude=` (an `#SBATCH` directive cannot read a file). The exclusion
stays on the allocation: repeating it as a step-level `srun --exclude` made
SLURM serialise the steps. Rebuild the missing stems by listing them
"<year> <stem>" per line and pointing `ONLY_STEMS_FILE` at that file.

**Kernel start-up race.** At high concurrency `jupyter_client` reserves five ZMQ
ports by binding to port 0 and closing them; the kernel re-binds moments later,
and a concurrent kernel on the same host can steal one in that window — the
loser dies with `Address already in use` before executing a cell, and papermill
writes an output notebook with zero executed cells and no exception. Handled by
a bounded retry (3 attempts, 1-10 s jitter) with the exit status propagated
explicitly: without that the loop ends on `sleep`, `bash -c` exits 0, and
`xargs` reports success for a cell that never wrote a partition. A per-cell
`JUPYTER_RUNTIME_DIR` on node-local scratch keeps the connection files off GPFS
but does not fix the race (their names carry UUIDs).

## The sidecar reducers (024e / 024f / 024g)

`scripts/024e_BuildBeaching_job.sh`,
`scripts/024f_BuildSurvivalOccupancy_job.sh` and
`scripts/024g_BuildSurvivalConnectivity_job.sh` fan out the same `(year, month)`
grid, but each cell reads only the compact BeachingForcing sidecar and the key
parquet — seconds per zarr, ~3 min for a whole month — so they are CPU-bound,
not GPFS-bound: small `--ntasks` (12), no `--spread-job`, no `--constraint`.
All three take `[regime] [hex_radius] [stagger_max_s]` positionally and the
rate-model member from the environment (`W_C`, `TAU_STRONG_HOURS`,
`TAU_CALM_DAYS`, `DELTA`, `BAND_M`, plus `MAX_FLOAT_DAYS` for 024e,
`OCCUPANCY_MAX_DAYS` for 024f and `CONNECTIVITY_MAX_DAYS` for 024g — which
must cover the largest age horizon 028 asks for),
so a sweep is a loop of submissions differing only in an export; each
reconstructs the member tag the notebook builds, so an executed notebook and
the parquet it wrote match by eye. The kernel-race retry and the per-cell
`JUPYTER_RUNTIME_DIR` above apply here too.

When *every* member has to be rebuilt at once — after a sidecar rebuild
invalidates all of them — `scripts/024efg_ReduceAll_job.sh` replaces that loop
of one-member submissions with a single 024d-shaped job whose unit of work is
one `(reducer, member, year, month)` cell: 25 members (14 for 024e, 7 for 024f,
4 for 024g) x the `(year, month)` grid — 1200 cells over four release years,
1.2k more per further year — driven by
`xargs -P ${SLURM_NTASKS}` over `srun -n1 -N1 --exact --no-kill` steps with the
same 3-attempt kernel-race retry and `OK`/`FAIL` tally. The members are a plain
table in the script rather than an environment sweep, so the exact set is read
off the source, and the tag each cell produces is rebuilt there exactly as the
notebooks build it (including the `_tf{}_tw{}` suffix a shore-type-sensitive
member carries). `tau_calm_days = 0` is how the notebooks spell an infinite
calm timescale; it tags as `tcinf`. `scripts/submit_024efg.sh` mirrors
`submit_024d.sh` — it expands the node blacklist into `--exclude`, sizes
`--ntasks` from the cell count (capped at `MAX_TASKS`, default 292) and passes
extra `sbatch` args through, so chaining it behind the sidecar build is
`scripts/submit_024efg.sh --dependency=afterok:<024d job id>`.
Reducers are light (`--mem-per-cpu=6G`, `--cpus-per-task=2`), so the 12 h limit
is slack, not an estimate.

## Parquet-only consumers

Single task, no Dask, one submission per run; each runs under
`set -euo pipefail` and passes `--cwd notebooks/` and `-p output_root`.
`season` is one of `DJF`/`MAM`/`JJA`/`SON`/`ALL` and defaults to `ALL`;
`member` is the opaque rate-model tag a 024e/024f/024g reducer put in the
store filename. Both are passed with `-r` (raw string): they are tags, not
Python literals, and `""` is a meaningful `member` that `-p` would parse
away.

| script | positional args |
|---|---|
| `024c_BuildHexConnectivity_job.sh` | `[regime] [year] [hex_radius]` (reads the 024 counts store, writes connectivity) |
| `025_HexHeatmaps_job.sh`, `026*_job.sh` | `[regime] [season] [hex_radius]` (also pass `-p data_root ${repo_root}/data` for the BSH H0 isobath) |
| `027_HexDistanceQuantiles_job.sh` | `[regime] [hex_radius] [season]` |
| `028_SubbasinConnectivityMatrix_job.sh` | `[regime] [hex_radius] [season] [member]` — empty `member` reads the unweighted 024c store, a tag reads 024g |
| `029_BeachingMaps_job.sh`, `030_SurvivalHeatmaps_job.sh` | `[regime] [hex_radius] [member] [season]` |
| `031_BeachingSweep_job.sh` | `[regime] [hex_radius] [members_csv] [labels_csv] [baseline_member] [season]` |

Every default is `regime=surface_stokes`, `hex_radius=6000`. `age_bin_days`,
`time_horizons_days_csv`, `disp_bin_km` and the origin/year filters fall back
to the notebook defaults and must match what the corresponding builder wrote.

`scripts/0-31-smoke-test.sh` runs the whole chain (000–031) on one host
against the demo subset in `data/`, with no SLURM — the integration check for
a fresh clone, not a scientific run.

## Cross-references

- [seeding.md](seeding.md), [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — what the runs produce.
- [../AGENTS.md](../AGENTS.md) — `pixi run` discipline, notebook conventions.
