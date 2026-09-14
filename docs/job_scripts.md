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
  026/  TimeHorizonMaps_<regime>_r<radius>m.png
  026a/ OriginSubbasinTimeHorizonMaps_<regime>_r<radius>m_<subbasin>.png
  026b/ OriginSubbasinYearTimeHorizonMaps_<regime>_r<radius>m_<subbasin>_<year>.png
```

`figure_dir` is derived from the `output_root` parameter inside each notebook
(no extra job-script argument); submit once per regime. `savefig` inherits the
notebook's `figure.dpi` (the `fig_dpi_scale` override), so saved panels match
the inline ones. `029`/`030`/`031` write the same way under `Figures/029`,
`Figures/030`, `Figures/031`, tagged with the rate-model member.

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

## Horizontal scaling for the GPFS-bound builder (024d)

`scripts/024d_BuildBeachingForcing_job.sh` fans out one independent
single-process papermill run per `(year, month)` cell (`xargs -P
${SLURM_NTASKS}`, so `njobs ≠ ntasks`: throttle concurrency without editing the
job). No MPI, no Dask; the bottleneck is streaming trajectory zarrs and hourly
Stokes files off GPFS. Three consequences, all defaults in the script:

- **`--spread-job`.** Maximise horizontal reach into the filesystem rather than
  node locality. Packing is a tail-latency disaster: the same 48 cells ran
  steps min/med/max `11:54 / 17:55 / 54:51` on 8 nodes against
  `10:51 / 16:43 / 20:54` on 22 — identical median, 2.6x worse tail, and wall
  time is set by the slowest cell. The flag disables the topology/tree plugin,
  which costs nothing without MPI.
- **No `--constraint`.** Pinning to sapphire (srp) nodes was inherited from
  Dask-backed jobs where IB reliability mattered; here it only shrinks the
  eligible pool and leaves jobs pending indefinitely.
- **`--ntasks` matched to the cell count** (`|YEARS| x 12 = 48`) and memory
  sized from measured `sacct MaxRSS` (024d peaks ~10.6 GB, so 8G/cpu x 2).
  Over-requesting either forces the job onto more nodes than it needs and
  pushes it into `(Resources)`.

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

## The sidecar reducers (024e / 024f)

`scripts/024e_BuildBeaching_job.sh` and
`scripts/024f_BuildSurvivalOccupancy_job.sh` fan out the same `(year, month)`
grid, but each cell reads only the compact BeachingForcing sidecar and the key
parquet — seconds per zarr, ~3 min for a whole month — so they are CPU-bound,
not GPFS-bound: small `--ntasks` (12), no `--spread-job`, no `--constraint`.
Both take `[regime] [hex_radius] [stagger_max_s]` positionally and the
rate-model member from the environment (`W_C`, `TAU_STRONG_HOURS`,
`TAU_CALM_DAYS`, `DELTA`, `BAND_M`, `MAX_FLOAT_DAYS`),
so a sweep is a loop of submissions differing only in an export; each
reconstructs the member tag the notebook builds, so an executed notebook and
the parquet it wrote match by eye. The kernel-race retry and the per-cell
`JUPYTER_RUNTIME_DIR` above apply here too.

## Parquet-only consumers

Single task, no Dask, one submission per run; each passes `--cwd notebooks/`
and `-p output_root`:

| script | positional args |
|---|---|
| `024c_BuildHexConnectivity_job.sh` | `[regime] [year] [hex_radius]` (reads the 024 counts store, writes connectivity; default `regime` is `surface_stokes`) |
| `025_HexHeatmaps_job.sh` | `[regime] [year] [hex_radius]` |
| `026*_job.sh`, `027_HexDistanceQuantiles_job.sh` | `[regime] [hex_radius]` |
| `028_SubbasinConnectivityMatrix_job.sh` | `[regime] [hex_radius]` (default `regime` is `surface_stokes`) |
| `029_BeachingMaps_job.sh`, `030_SurvivalHeatmaps_job.sh` | `[regime] [hex_radius] [member]` |
| `031_BeachingSweep_job.sh` | `[regime] [hex_radius] [members_csv] [labels_csv] [baseline_member]` |

## Cross-references

- [seeding.md](seeding.md), [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — what the runs produce.
- [../AGENTS.md](../AGENTS.md) — `pixi run` discipline, notebook conventions.
