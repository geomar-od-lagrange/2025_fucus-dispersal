# Cluster job scripts and notebook output layout

`scripts/0NN_*_job.sh` runs notebook `0NN` on NESH via SLURM; the
executed copy with rendered figures lands under
`notebooks_executed/Visualisations/`.

```
notebooks/0NN_*.md            ← source of truth (jupytext)
notebooks/0NN_*.ipynb         ← code-only, in-sync with .md, no outputs
notebooks_executed/
  Visualisations/
    0NN_*.ipynb               ← papermill output, with rendered figures
    0NN_*_<regime>_<year>.ipynb   (when swept across regime / year)
```

Sources committed in `notebooks/` are stripped of cell outputs (see
[../AGENTS.md](../AGENTS.md) §Conventions › Notebooks). Production
runs go through papermill on the cluster, never overwrite the
sources, and write per-parameter copies under `notebooks_executed/`.

## Multi-task dask layout

Every job script uses the same SLURM shape: `--ntasks=N` ≥ 3 with one
task playing each role.

| SLURM task | role |
|------------|------|
| 0          | dask scheduler + one local worker set |
| 1          | papermill (connects via `$SCHEDULER_FILE`) |
| 2…N-1      | additional dask worker tasks, each pinned to its node |

`SCHEDULER_FILE` is the bootstrap contract: the scheduler writes it,
the papermill task and any extra workers poll for it
(`Client(scheduler_file=os.environ["SCHEDULER_FILE"])` in the notebook;
`dask worker --scheduler-file ...` in workers). A `trap cleanup EXIT`
removes the file and kills child jobs on script exit.

`sleep 30` after launching task 0 is a coarse hand-off to give the
scheduler time to come up; `dask worker --scheduler-file` polls the
file so no explicit wait loop is needed on the worker side.

## Notebook output policy

`notebooks/*.ipynb` files in the repo are code-only:

```sh
pixi run jupytext --sync foo.md                              # default
pixi run jupyter nbconvert --clear-output --inplace foo.ipynb # if executed
```

Local execution for testing is fine; strip outputs before commit.
The cluster path doesn't touch `notebooks/` — papermill reads it,
writes its parameterised copy under `notebooks_executed/`. That keeps
diffs small (no base64 PNGs) and figures fresh (always from the most
recent cluster run).

`notebooks_executed/` is gitignored or thinly tracked depending on
the run; treat the contents as artifacts, not source.

## Per-job constants

Every job script sets two paths and exports `SCHEDULER_FILE`:

```bash
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export SCHEDULER_FILE=${repo_root}/.scheduler_${SLURM_JOB_ID}.json
```

`output_root` is passed to papermill via `-p output_root ${output_root}`
— the single contract between job script and notebook (see notebook
parameters cells). There is no `$FUCUS_OUTPUT_ROOT` env variable;
the papermill `-p` is the only way the notebook learns where to read
from / write to.

`http_proxy` / `https_proxy` are set for HTTP traffic from the env
(e.g. CMEMS downloads inside notebook 002). `no_proxy` defensively
excludes intra-cluster TCP — dask uses raw Tornado streams that
ignore the proxy vars anyway, but anything else that dials a node
address shouldn't touch the proxy.

## Cross-references

- [seeding.md](seeding.md) — what 010 produces, consumed by 020–025.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — store layout 024 writes.
- [../AGENTS.md](../AGENTS.md) — notebook conventions, `pixi run` discipline.
