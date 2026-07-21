# Cluster job scripts and notebook output layout

`scripts/0NN_*_job.sh` runs notebook `0NN` on NESH; the rendered copy
lands under `notebooks_executed/Visualisations/`.

```
notebooks/0NN_*.md             ← jupytext source of truth
notebooks/0NN_*.ipynb          ← code-only, in-sync with .md, no outputs
notebooks_executed/
  Visualisations/
    0NN_*_<regime>_<year>.ipynb   ← papermill output, with figures
```

`notebooks/*.ipynb` are stripped before commit
(`pixi run jupyter nbconvert --clear-output --inplace` after any
local execution). Cluster runs never overwrite the source — papermill
reads `notebooks/`, writes per-parameter copies under
`notebooks_executed/`. Treat `notebooks_executed/` as artifact, not
source.

## Standalone figure PNGs

The map notebooks `026`, `026a`, `026b` also write their four-panel
time-horizon maps as standalone PNGs under the outputs tree (outside the
repo), one subdir per stage — so the figures are usable without opening
the executed notebook:

```
output_root/Figures/
  026/  TimeHorizonMaps_<regime>_r<radius>m.png
  026a/ OriginSubbasinTimeHorizonMaps_<regime>_r<radius>m_<subbasin>.png
  026b/ OriginSubbasinYearTimeHorizonMaps_<regime>_r<radius>m_<subbasin>_<year>.png
```

`figure_dir` is derived from the `output_root` parameter inside each
notebook (no extra job-script argument). Submit once per regime
(`surface`, `bottom`, `surface_stokes`) to cover all regimes. `savefig`
inherits the notebook's `figure.dpi` (the `fig_dpi_scale` 2× override),
so saved panels match the inline ones.

## Multi-task dask layout

`--ntasks=N ≥ 3` with one role per SLURM task:

| task    | role |
|---------|------|
| 0       | dask scheduler + one local worker set |
| 1       | papermill (connects via `$SCHEDULER_FILE`) |
| 2…N-1   | extra dask worker tasks |

Bootstrap is one file. Task 0 writes `$SCHEDULER_FILE`; the papermill
task and any extra workers poll for it (via
`Client(scheduler_file=os.environ["SCHEDULER_FILE"])` in the notebook
and `dask worker --scheduler-file ...` in workers). A `trap cleanup
EXIT` removes the file on script exit. The hard-coded `sleep 30`
after task 0 launches gives the scheduler time to come up before
papermill starts.

## Per-job constants

Every script sets:

```bash
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export SCHEDULER_FILE=${repo_root}/.scheduler_${SLURM_JOB_ID}.json
```

`output_root` is passed to papermill via `-p output_root ${output_root}`
— the only contract between job script and notebook. No
`$FUCUS_OUTPUT_ROOT` env variable.

`http_proxy`/`https_proxy` are set for HTTP traffic from the env
(e.g. CMEMS in 002); `no_proxy` excludes intra-cluster TCP so dask
node addresses don't dial the proxy.

## Horizontal scaling for GPFS-bound cells (024d / 024e)

The beaching and survival-occupancy passes fan out one independent
single-process papermill run per `(year, month)` cell. They hold no MPI
communicator and no Dask cluster; the bottleneck is streaming trajectory
zarrs and hourly Stokes files off GPFS. Three consequences, all defaults in
`scripts/024{d,e}_*_job.sh`:

- **`--spread-job`.** Maximise horizontal reach into the filesystem rather
  than node locality. Packing is a tail-latency disaster: the same 48 cells
  ran steps min/med/max `11:54 / 17:55 / 54:51` on 8 nodes against
  `10:51 / 16:43 / 20:54` on 22 — identical median, 2.6x worse tail, and
  since wall time is set by the slowest cell `024e` took 1:13:44 against
  20:54. The flag disables the topology/tree plugin, which costs nothing
  without MPI.
- **No `--constraint`.** Pinning to sapphire (srp) nodes was inherited from
  Dask-backed jobs where IB reliability mattered. For these it only shrinks
  the eligible pool and leaves jobs pending indefinitely.
- **`--ntasks` matched to the cell count** (`|YEARS| x 12 = 48`). More tasks
  than cells inflates the allocation with idle slots and makes the job harder
  to place. Sweeps multiply the grid, so raise it on the command line.

Memory is sized from measured `sacct MaxRSS`, not guessed: 024d peaks
~10.6 GB (8G/cpu x 2), 024e ~20.3 GB (12G/cpu x 2). Over-requesting is what
forces a job onto more nodes than it needs and pushes it into `(Resources)`.

**Kernel start-up race.** At high concurrency `jupyter_client` reserves five
ZMQ ports by binding to port 0 and closing, then the kernel re-binds moments
later; a concurrent kernel on the same host can take one in that window and
the loser dies with `Address already in use` before executing a cell
(papermill then writes an output notebook with zero executed cells and no
exception). The collision domain is the host's ephemeral port space, so it
scales with per-node concurrency. Handled by a bounded retry (3 attempts,
1-10 s jitter) with the exit status propagated explicitly — without that the
loop ends on `sleep`, `bash -c` exits 0, and `xargs` reports success for a
cell that never wrote a partition. A per-cell `JUPYTER_RUNTIME_DIR` on
node-local scratch is also set; it does not fix the port race (connection
filenames carry UUIDs) but keeps those files off GPFS.

## Cross-references

- [seeding.md](seeding.md), [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — what the runs produce.
- [../AGENTS.md](../AGENTS.md) — `pixi run` discipline, notebook conventions.
