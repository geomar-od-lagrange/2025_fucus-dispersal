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
  025/  HexRelativeDensity[DE|ByOrigin]_<regime>_<season>_r<radius>m_<member-or-unweighted>[_part<n>].png
  026/  TimeHorizonMaps_<regime>_<season>_r<radius>m_<member-or-unweighted>_<percentage|dilution>.png
        MedianAge_<regime>_<season>_r<radius>m_<member-or-unweighted>.png
  026a/ OriginSubbasinTimeHorizonMaps_<regime>_<season>_r<radius>m_<member-or-unweighted>_<subbasin>_<column>.png
        OriginSubbasinOverview_<regime>_<season>_r<radius>m_<member-or-unweighted>_T<h>d_dilution.png
  026b/ OriginSubbasinYearTimeHorizonMaps_<regime>_<season>_r<radius>m_<member-or-unweighted>_<subbasin>_<year>_<column>.png
  027/  DistanceQuantiles_<regime>_r<radius>m_<season>.png
  028/  SubbasinConnectivity_{raw,emission,interannual}_<regime>_<season>_<member-or-unweighted>_T<h>d.png
  029/  <view>_<regime>_r<radius>m_<season>_<member>.png
  030/  {SurvivalHeatmaps,SurvivalFraction,DriftingFractionCurve}_<regime>_r<radius>m_<season>_<member>.png
  031/  BeachingSweep{Stats,Maps}_<regime>_r<radius>m_<season>.png
  035/  AreaTimeSeries_<season>_r<radius>m.png
  036/  SeedingMap_r<radius>m.png
```

`figure_dir` is derived from the `output_root` parameter inside each notebook
(no extra job-script argument); submit once per (regime, season). Every figure
is saved print-ready — one full 180 mm text width at 300 dpi — so the PNG is
the manuscript figure, not a preview (rationale in
[visualisations.md](visualisations.md)).

`025`–`029` and `033` additionally write machine-readable exports under
`output_root/Exports/<stage>/`: GeoJSON (hex geometry, EPSG:4326) for the map
products of `025`/`026`/`026a`/`026b`/`027`/`029`/`033`, long-form CSV for
`028`'s connectivity matrices, `026`'s median particle age per hex
(`026_medianage_<regime>_<season>_<member-or-unweighted>.geojson`), its area-holding-a-share table
(`026_areas_<regime>_<season>_<member-or-unweighted>.csv`) and `029`'s afloat
share per origin subbasin
(`afloat_by_origin_<regime>_<season>_<member>.csv`).

## The 2D-field preprocessor (003)

`scripts/003_prepare_2d_fields_job.sh` fans out one `srun` step per BSH
c-file (`xargs -P ${SLURM_NTASKS}`), each running
`notebooks/003_prepare_2d_fields.py` inside the parcels container to write
the three variants of that timestamp. Three optional environment variables
pick the work, so a rerun does not mean redoing the decade:

| variable | effect |
|---|---|
| `YEARS="2021 2025"` | restrict the year loop (default 2016 … 2025) |
| `CFILE_LIST=<path>` | process exactly the c-file paths listed there, one per line |
| `MISSING_ONLY=1` | drop every c-file whose three outputs all exist |

```bash
YEARS="2021 2025" MISSING_ONLY=1 sbatch --ntasks=33 --time=02:00:00 \
    scripts/003_prepare_2d_fields_job.sh
```

003 itself already skips a c-file whose outputs are all present, but only
after an `srun` step and a container start-up; `MISSING_ONLY` makes the same
decision up front from a single `ls` of `2d_fields/` (one directory listing
beats three `stat()` calls per candidate on a shared parallel filesystem).
Steps run as `srun -n1 -N1 --exact --no-kill` — without `-N1` a step inherits
the job's node count and dies with "can't run 1 processes on N nodes" the
moment the allocation spans more than one node.

Counting the holes is the same set arithmetic the filter does: a complete
year has `4 timestamps/day x 2 resolutions x 3 variants` files, 8760 in a
common year and 8784 in a leap year, minus whatever the BSH archive is
itself missing.

## The Parcels runs (010)

`scripts/010_FucusDispersal_{surface,bottom,surface_stokes}_job.sh` fan out one
single-process papermill run per release (`xargs -P ${SLURM_NTASKS}` over `srun
--ntasks=1 --exact` steps), each writing one trajectory zarr under
`output_root/Trajectories/<regime>/<year>/`. Two positional args:

| arg | default | meaning |
|---|---|---|
| `$1` | `2019` | release year |
| `$2` | `73` | number of releases to run, taken from the front of the schedule |
| `RELEASE_NS` (env) | unset | run exactly these release indices instead of a prefix |
| `VELOCITY_SCALE` (env, bottom only) | `1.0` | 010's `velocity_scale`: uniform factor on U and V |

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
RELEASE_NS=1,19 scripts/submit_010.sh surface_stokes 2022   # 2 holes, --ntasks=2
```

`--ntasks` has to equal the release count and an `#SBATCH` directive cannot read
the job script's own arguments, so the wrapper sets both from one number — a
29-release job books 29 tasks, not 73. It also expands
`scripts/node-blacklist.txt` into `--exclude=` (as `submit_024d.sh` does) and
passes further `sbatch` flags through.

Submission goes through `robust_sbatch` (`scripts/slurm_submit.sh`), so a
controller timeout that creates the job anyway is adopted rather than queued
twice. The job name is the idempotency key and names the work —
`010-<regime>-<year>`, or `010-<regime>-<year>-n<i>_<j>` for a refill — so
rerunning the same command after a failed submit adopts instead of duplicating.

### Velocity-scaled variant regimes

`VELOCITY_SCALE` runs the bottom 2D fields with Parcels'
`Field.set_scaling_factor` on U and V of both nest levels; the kernel stays
plain `AdvectionRK4`. Off 1.0 the run writes under its own regime,
`<regime>_vs<scale>` with `.` spelled `p` — `0.8` → `bottom_vs0p8`, in the
directory name and the zarr stem — so every downstream stage keys on it like
any other regime. A uniform factor in a time-varying flow is not a slowed
replay of the 1× trajectories, so a variant needs its own runs. The scale goes
into the job name (`010-bottom-vs0.8-<year>…`), so a scaled run neither adopts
nor blocks the 1× one; `submit_010.sh` refuses `VELOCITY_SCALE` for a job
script that does not read it. At 1.0 the trajectories are bitwise identical to
a run without the factor (multiplying by one is exact).

```bash
# 2019, one release per month (doys 1, 36, …, 336), at 0.8x
RELEASE_NS=0,7,12,18,24,31,37,43,49,55,61,67 VELOCITY_SCALE=0.8 \
    scripts/submit_010.sh bottom 2019
```

### Refilling lost releases

`--no-kill` means a release can be lost without failing its job: the job exits
`COMPLETED` and the hole shows up only as a missing zarr. Two causes seen in
practice — a node dying under the job, and a jupyter kernel losing the race for
an ephemeral port when 73 of them start at once on a handful of nodes
(`zmq.error.ZMQError: Address already in use`, then `Kernel died before replying
to kernel_info`). Both are transient and per-step, not per-release-date, so the
fix is to rerun exactly the lost indices.

A hole sits in the middle of the schedule, so `n_releases` — a prefix count —
cannot select it. `RELEASE_NS` takes the indices instead, and `--ntasks` follows
the length of the list:

```bash
RELEASE_NS=1 scripts/submit_010.sh surface_stokes 2021    # doy 6, 2021-01-06
```

Find the holes by diffing the schedule against the zarrs on disk — the stem is
`Fucus_BSH_<YYYYMMDD>T001500_<regime>_dt60min_seed<N>.zarr`, and release `n` is
`doy = 1 + 5n` of the release year. **Check `squeue` first:** a job still
running has steps still writing, and a missing zarr for a release whose step is
alive is work in progress, not a hole. `sacct -j <id> -o JobID,State` settles
it — the step index *is* the release index `n`, so a `FAILED` `<id>.19` is
release `n=19`, and rerunning a release whose step is `RUNNING` would put two
papermills on one output path.

`--mem-per-cpu` is sized from `sacct` peaks over completed runs, not guessed:
11.42 GiB for `bottom`, 10.31 GiB for `surface`, 6.58 GiB for `surface_stokes`,
hence 20G / 20G / 12G at roughly peak x 1.7. This is the single biggest lever on
how fast the set of runs gets through the queue -- at the old 60G a 73-task job
was memory-bound to ~2.4 tasks per node and booked 30 nodes, so a handful of
years saturated the account's share and the rest sat in `(Priority)`.
`--time=36:00:00` (the `base` cap is 48 h): observed wall times run 5-19 h, and
the fieldset the notebook builds grows with every year of 2D fields on disk,
because `010` globs the whole `2d_fields/` directory rather than the run window.

`--no-kill` plus, in the fan-out, `unset SLURM_CPU_BIND*` and `srun -n1 -N1
--exact` with explicit per-step `--cpus-per-task` / `--mem-per-cpu`: without the
unset and the `-N1`, every concurrent step dies at launch with "CPU binding
outside of job step allocation" (the same failure
`024d_BuildBeachingForcing_job.sh` documents).

The batch step's own log (the fan-out's per-release papermill output stays
inline in it, there being no separate step-level `--output`) lands at
`output_root/logs/010/%x_%j.out`, created by `submit_010.sh` before `sbatch`
since `#SBATCH --output` opens its target at job start.

Every job script points `#SBATCH --output` at `output_root/logs/<stage>/`, and
the wrapper that submits it creates that directory first:

| stage | log directory | created by |
|---|---|---|
| 010 | `logs/010/` | `submit_010.sh` |
| 024 / 024b / 024c | `logs/aggregates/` | `submit_aggregates.sh` |
| 024d | `logs/024d/` | `submit_024d.sh` |
| 024e / 024f / 024g | `logs/024efg/` | `submit_024efg.sh` |
| 025-031 | `logs/025-031/` | `submit_downstream_chain.sh` |

The directive lives in the job script rather than in the submitter's argument
list so that a hand `sbatch`, or a wrapper called directly instead of through
`submit_downstream_chain.sh`, still lands its log in the same place instead of
dropping a `slurm-<id>.out` in the repo root.

### Forcing gaps are silent

`010` globs `2d_fields/`, intersects the fine and coarse stem sets and drops
whatever is missing from either, with a warning. The surviving stems are handed
to `FieldSet.from_netcdf` as an explicit `timestamps` list, so a hole is not a
missing time -- it is a longer interval between two samples, and Parcels
interpolates straight across it. Nothing raises, and
`allow_time_extrapolation=False` does not catch it either, because the hole is
interior to the time axis. Check `2d_fields/` for missing 6-hourly stems across
every year a run window touches *before* submitting.

**The 2025 holes.** 2016-2024 are hole-free: every 6-hourly stem is present in
both grids for both regimes. 2025 is not, and the cause is upstream, not the
extraction: 25 raw BSH c-files in the source archive are corrupt (HDF read
error; 1-50 MB against 50-96 MB for their neighbours) and 5 timestamps were
never delivered. Re-running `003` cannot fix either; it needs a BSH
re-delivery. The result, per regime:

| | |
|---|---|
| unusable stems in 2025 | 30, on 29 distinct days |
| distinct holes | 27 |
| isolated 6 h holes | 25 |
| longest hole | 18 h, 2025-07-24T12:15 - 2025-07-25T06:00 |
| second longest | 12 h, 2025-08-24T18:15 - 2025-08-25T06:00 |
| first hole | 2025-01-12T12:15 |
| last usable instant | 2025-12-31T18:00 |

So the holes are scattered single 6 h blocks, not a contiguous outage: the
record is hole-free from 2016-01-01T00:15 to **2025-01-12T12:00**, and from
there on 27 short punctures spread over the rest of 2025.

**Which releases are affected.** A release is hole-free iff its whole 220-day
window closes before the first hole:

| release year | hole-free releases | interpolating releases |
|---|---|---|
| 2016-2023 | all 73 | none |
| 2024 | 32 (n = 0-31, doys 1-156, last window closing 2025-01-10) | 41 (n = 32-72, doys 161-361), crossing 1-17 holes, median 7 |
| 2025 | none | all 29, crossing 17-22 holes, median 20 |

The 2025 runs therefore carry a velocity field that is linearly smoothed across
~20 short windows spread over their 220 days, and the second half of 2024 a
handful. Nothing is truncated -- every submitted release has forcing at both
ends of its window -- but 2024 and 2025 are not forced at the same fidelity as
2016-2023, and the 18 h gap at 2025-07-24/25 is the one worth naming when
reading those two years.

## Multi-task dask layout

`--ntasks=N ≥ 3`, one role per SLURM task: task 0 the dask scheduler plus a
local worker set, task 1 papermill (connecting via `$SCHEDULER_FILE`), tasks
2…N-1 extra worker tasks.

Bootstrap is one file. Task 0 writes `$SCHEDULER_FILE`; the papermill task and
any extra workers poll for it (`Client(scheduler_file=os.environ["SCHEDULER_FILE"])`
in the notebook, `dask worker --scheduler-file ...` in workers). A `trap cleanup
EXIT` removes it on exit, and the batch step polls for the file for up to 300 s
before launching papermill, failing the job if it never appears — the notebook's
`Client(scheduler_file=...)` waits forever otherwise, so a scheduler that dies
at start-up would burn the whole walltime instead of the first minute.

No `--interface` is passed to scheduler or workers: only the srp nodes carry
`ib0`, so pinning it kills the scheduler outright on a clk allocation (`'ib0'
is not a valid network interface`) or, on a mixed one, publishes an address the
papermill task cannot reach (`OSError: [Errno 101] Network is unreachable`).
Dask binds each node's default address instead; these jobs are GPFS-bound and
minutes long, so IB bandwidth never set their wall time.

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
`TAU_CALM_DAYS`, `DELTA`, `TRAP_FLAT`, `TRAP_WALL`, `BAND_M`, and
`MAX_AGE_DAYS`, the reduction window shared by all three, default 180 d,
which must cover the largest age horizon any consumer asks for),
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
calm timescale; it tags as `tcinf`. The reduction window is one variable,
`max_age_days=180`, passed to every cell — the parquets record it in their
attrs, so a partial rerun at another window is detectable downstream.
`scripts/submit_024efg.sh` mirrors
`submit_024d.sh` — it expands the node blacklist into `--exclude`, sizes
`--ntasks` from the cell count (capped at `MAX_TASKS`, default 292) and passes
extra `sbatch` args through, so chaining it behind the sidecar build is
`scripts/submit_024efg.sh --dependency=afterok:<024d job id>`.
Reducers are CPU-light but decompress one sidecar's `(trajectory, obs)`
arrays over the whole window at once, so memory scales with `max_age_days`:
`--mem-per-cpu=8G`, `--cpus-per-task=2` (16 GB per cell) against measured
peaks of 8.6 / 11.1 / 10.9 GB (024e / 024f / 024g, one 2019-01 cell at 180 d).
The 12 h limit is slack, not an estimate.

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
| `025_HexHeatmaps_job.sh`, `026*_job.sh` | `[regime] [season] [hex_radius] [member]` — empty `member` reads the unweighted 024 counts store, a tag reads 024f (also pass `-p data_root ${repo_root}/data` for the BSH H0 isobath) |
| `027_HexDistanceQuantiles_job.sh` | `[regime] [hex_radius] [season]` |
| `028_SubbasinConnectivityMatrix_job.sh` | `[regime] [hex_radius] [season] [member]` — empty `member` reads the unweighted 024c store, a tag reads 024g plus the same member's 024e store (stranded column) |
| `029_BeachingMaps_job.sh`, `030_SurvivalHeatmaps_job.sh` | `[regime] [hex_radius] [member] [season]` |
| `031_BeachingSweep_job.sh` | `[regime] [hex_radius] [members_csv] [labels_csv] [baseline_member] [season]` |
| `035_AreaTimeSeries_job.sh` | `[season] [hex_radius] [member]` — draws `surface_stokes` at 026 tag `member` (default the production member) and `bottom` unweighted, from 026's area CSVs |
| `036_SeedingMap_job.sh` | `[hex_radius]` — no season or member; also reads the first two time steps of one regime-year of trajectory zarrs (notebook defaults `surface_stokes` 2019). Not in the consumer driver: the release set does not change between runs |

Every default is `regime=surface_stokes`, `hex_radius=6000`. `age_bin_days`,
`max_age_days`, `time_horizons_days_csv`, `disp_bin_km` and the origin/year
filters fall back to the notebook defaults. A consumer of the 024e/024f/024g
stores asserts them against the partitions' build attrs (and that every pooled
partition shares one window and rate-model point), so a mismatch fails the
run; the 024/024c counts-derived stores record no attrs, and there the
defaults must match what the builder wrote.
Each ends in `jobinfo || true`: the accounting tool is the last command, so
without that its status, not papermill's, would be the job's — and an
`afterok` dependency would fire on a figure that was never drawn. Same reason
`024c` carries `set -euo pipefail`.

### The whole consumer matrix in one job

`scripts/025-031_ConsumersAll_job.sh` runs every consumer figure of the study
as one `srun` step per `(stage, regime, season, member)`, `xargs -P
${SLURM_NTASKS}` over the same fan-out the 024d/024efg drivers use. Each step
invokes the per-stage script above, so those stay the single source of the
papermill contract and this driver only decides which combinations run:

| stages | coverage | units |
|---|---|---|
| 025, 026, 026a, 026b | surface_stokes production member + bottom unweighted, × 5 seasons | 40 |
| 027 | {surface_stokes, bottom} × 5 seasons | 10 |
| 028 | surface_stokes unweighted + production member, bottom unweighted, × 5 seasons | 15 |
| 029, 030 | surface_stokes, production member, × 5 seasons | 10 |
| 031 | surface_stokes × 5 seasons, `members_csv` = every 024e member on disk | 5 |
| 035 | both regimes in one figure × 5 seasons, second wave | 5 |

The units run as one flat `xargs` pool with no ordering among them, except
035: it reads the 026 area CSVs of both regimes, so the 035 units form a
second wave that starts once every other unit has finished.
`LIST_ONLY=1 bash scripts/025-031_ConsumersAll_job.sh` prints the 85 units;
`STAGES` / `REGIMES` / `SEASONS` restrict the matrix. The 031 member list is
read off `HexAggregates/` inside the job, not at submit time — the reducers
write those stores minutes earlier in the chain below, so any hardcoded list
would be stale by construction.

Per-task shape is `8 cpu × 16G = 128 GB`, sized from the heaviest unit rather
than the median. The counts-store consumers reduce one release-year partition
at a time with the season and the `max_age_days` window pushed into the
parquet read, so their peak is one partition's read, not the pooled store:
measured process-tree peak RSS over ten release years (`surface_stokes`,
`ALL`, unweighted) is 16.2 GiB for 026, 15.5 for 026b and 12.2 for 025; a
single season or a 024f member read stays under 3.2 GiB; 027 and 028 are
0.2 and <0.1 GiB (`sacct`, four years). The shape is generous for the heatmap
family and stays until 029/030/031 are measured at 180 d. One job at that shape beats ~85
one-task allocations queueing independently behind the same dependency for
work that is minutes long; the price is that the cheap two thirds of the
matrix book memory they do not need.

`scripts/0-31-smoke-test.sh` runs the whole chain (000–031) on one host
against the demo subset in `data/`, with no SLURM — the integration check for
a fresh clone, not a scientific run.

## Unattended: the downstream chain

`scripts/submit_downstream_chain.sh <010-job-id> ...` queues everything below
the Parcels runs as one SLURM dependency graph and exits. The 010 runs take
24–48 h, far past any session watching them, so the chain is submitted once
and left alone.

```
A  024 + 024b + 024c per (regime, year)   afterany on that pair's 010 job
B  024d sidecars, surface_stokes          afterany on the surface_stokes 010 jobs
C  024e/024f/024g reducers                afterok on B
D  025-031 consumers, 85 units            afterok on every A job and on C
E  headline regeneration                  afterok on D
```

`afterany` where the stage discovers its own work on disk (a release year
whose 010 run died contributes the zarrs it wrote, and the notebooks pool what
exists); `afterok` where a partial input would silently corrupt the output.
`--kill-on-invalid-dep=yes` on the `afterok` edges turns an unsatisfiable
dependency into a cancelled job rather than one pending forever.

The `(regime, year)` pairs are read off the 010 jobs themselves — `scontrol
show job <id>` still carries the `Command` line, whose script name is the
regime and whose first argument is the release year. Pass
`REGIME_YEARS="bottom:2020 surface_stokes:2020"` when scontrol has forgotten a
completed job; every pair then waits on the whole 010 id set. `ALLOW_MISSING=1`
is what lets `submit_aggregates.sh` queue a year whose trajectory directory
does not exist yet, which is the normal case here.

### Submitting against a sick controller

Every `sbatch` in the chain goes through `robust_sbatch`
(`scripts/slurm_submit.sh`), because a loaded slurmctld fails with `Socket
timed out on send/recv operation` *after* creating the job — a naive retry
then queues the work twice. So each stage carries a job name that is unique to
it and stable across reruns (`<tag>-<regime>-<year>-024b`, `<tag>-024d`, …,
where `<tag>` hashes the 010 id set), and the queue rather than sbatch's exit
status decides what still needs submitting:

| situation | what happens |
| --- | --- |
| sbatch times out, job was created | `squeue --me -n <name>` finds it, id adopted |
| sbatch times out, nothing created | retried, 15→120 s backoff, 6 attempts |
| squeue unreachable too | give up on that stage — "unknown" never resubmits |
| name already queued (rerun) | adopted, never resubmitted |

The chain is therefore idempotent: rerun the identical command after a partial
failure and it adopts what is queued, submits only the rest, and prints an
`adopted`/`submitted` column saying which was which. Stages queued under names
this chain cannot reconstruct — an earlier chain, or a hand submission — are
adopted by id instead: `ADOPT_A="<024/024b/024c ids>" ADOPT_B=<024d id>
ADOPT_C=<024efg id>` skips those stages entirely and hangs D and E off the ids
given, which the chain takes on trust. An aborted run still
prints the stages it did queue plus the command to resume. A stage that has
already *finished* has left the queue and would be submitted again, so resume
while the chain is still pending.

Logs land under `output_root/logs/<stage>/` (`aggregates`, `024d`, `024efg`,
`025-031`, `headline`), the chain's own submission record under
`logs/chain/chain_<timestamp>.log`. Stage E is a hook: if
`scripts/headline_report_job.sh` is absent the chain says so and submits
nothing, and the printed warning carries the one `sbatch` line that wires it
later.

## Cross-references

- [seeding.md](seeding.md), [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — what the runs produce.
- [../AGENTS.md](../AGENTS.md) — `pixi run` discipline, notebook conventions.
