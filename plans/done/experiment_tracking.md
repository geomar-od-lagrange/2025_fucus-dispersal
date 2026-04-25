> Superseded by `plans/wrapup.md` §1 (filename simplification). Original plan retained as historical record.

# Experiment Tracking: Notebooks and Zarr Stores

## Problem

We run many experiments with near-identical parameters:
- 10 years × ~36 start dates × 3 experiment types × N velocity factors
- Each produces: 1 executed notebook + 1 zarr trajectory store
- Parameters evolve (grid fix, particle count, dt changes)
- Need to distinguish "pre-fix" from "post-fix" runs
- Need to rerun subsets without clobbering good results
- Need to know which parameters produced which output

Current naming: `Fucus_BSH_{date}_{type}_dt{dt}min_vf{vf}_seed{seed}.zarr`
encodes some params in the filename but not all (no grid version, no
particle count, no max_age).

## Option A: Hash-based versioning

Compute a short hash of all parameters. Append to filenames.

```
params = {start_date, experiment_type, velocity_factor, max_age_days,
          calc_dt_mins, particles_per_cell, grid_version, ...}
h = hashlib.sha256(json.dumps(params, sort_keys=True)).hexdigest()[:8]
```

Output: `Fucus_2019-01-01_surface_vf1.0_a3b7c2d1.zarr`

Pros:
- Automatic deduplication — same params always give same hash
- Rerun with different params → different hash → no collision
- Tiny code change (hash in output_filename template)

Cons:
- Hash is opaque — can't tell what changed without a lookup table
- Need a manifest file mapping hash → full params
- Easy to accumulate orphaned outputs

## Option B: Explicit version tag

Add a `--run-tag` or `version` parameter. Bump manually.

```
run_tag = "v2_grid_fix"
```

Output: `Fucus_2019-01-01_surface_vf1.0_v2_grid_fix.zarr`

Pros:
- Human-readable — immediately clear which generation
- Simple grep/ls to find all outputs of a version
- Easy to purge: `rm *_v1_*`

Cons:
- Manual discipline required
- Doesn't capture parameter changes within a version
- Risk of forgetting to bump the tag

## Option C: Directory-based isolation

Each experiment batch gets its own output directory.

```
output/
  Trajectories/
    2025-04-10_grid_fix_50d/
      2019/Fucus_BSH_20190101_surface_*.zarr
    2025-04-09_pre_fix_220d/
      2019/Fucus_BSH_20190101_surface_*.zarr
```

Pros:
- Complete isolation — old and new coexist without collision
- Can diff entire directories
- Natural for "keep old results for comparison"

Cons:
- Path proliferation on shared filesystem
- Downstream notebooks need to know which directory to read
- Doesn't prevent param drift within a directory

## Option D: Metadata in zarr attrs

Write all parameters into the zarr store's `.zattrs` at creation time.
Keep current naming. Query params from the output itself.

```python
ds.attrs["params"] = json.dumps({...})
ds.attrs["grid_version"] = "v2_roll_fpoint"
ds.attrs["created"] = datetime.now().isoformat()
```

Pros:
- Self-describing output — params travel with the data
- No naming convention changes
- Can build a catalog by scanning zarr attrs

Cons:
- Doesn't prevent filename collisions (same name, different params)
- Need tooling to scan/query attrs
- Papermill already captures params in the notebook — duplication?

## Option E: Papermill + notebook as manifest

Papermill already injects parameters into the executed notebook. The
notebook IS the record of what ran. Use the notebook filenames as the
primary index, and derive zarr paths from them.

Add a cell at the end of 010 that writes a small JSON sidecar:

```python
manifest = {
    "params": {all papermill params},
    "zarr_path": str(output_path),
    "notebook": str(output_notebook_path),
    "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"]).strip(),
    "timestamp": datetime.now().isoformat(),
}
json.dump(manifest, open(output_path.with_suffix(".json"), "w"), indent=2)
```

Pros:
- Notebook already exists — just formalize it as the manifest
- Git SHA captures code version (including grid fix)
- JSON sidecar is trivial to scan/aggregate
- Works with any naming convention

Cons:
- One more output file per run
- Git SHA requires clean working tree to be meaningful

## Option F: UUID + nanosecond timestamp + params CSV

Every run gets a UUID and a nanosecond ISO timestamp. Both appear in
the zarr name, the notebook name, and a row in a params CSV. Filenames
retain human-readable metadata (experiment type, start date, etc.) for
quick ls/glob, but the UUID makes each run globally unique.

### Naming

```
010_FucusDispersal_{type}_{start_date}_{timestamp}_{uuid8}.zarr
010_FucusDispersal_{type}_{start_date}_{timestamp}_{uuid8}.ipynb
```

Example:
```
010_FucusDispersal_surface_2019-01-01_20260410T143022.123456789_a3b7c2d1.zarr
010_FucusDispersal_surface_2019-01-01_20260410T143022.123456789_a3b7c2d1.ipynb
010_FucusDispersal_bottom_2019-01-01_20260410T143025.987654321_f9e8d7c6.zarr
```

The `{type}_{start_date}` prefix lets you glob for all surface runs or
all 2019 runs without touching the CSV. The timestamp gives chronological
ordering. The UUID prevents any collision even for identical params
submitted simultaneously.

### Params CSV

One CSV per output directory (or one global), appended per run:

```csv
uuid,timestamp,start_date,experiment_type,velocity_factor,max_age_days,calc_dt_mins,particles_per_cell,git_sha,zarr_path,notebook_path
a3b7c2d1,2026-04-10T14:30:22.123456789,2019-01-01,surface,1.0,220,5,10,9ae0e3d,output/Trajectories/2019/surface_2019-01-01_20260410T143022.123456789_a3b7c2d1.zarr,notebooks_executed/TrajectoryCalc/surface_2019-01-01_20260410T143022.123456789_a3b7c2d1.ipynb
```

### Generation

In the notebook (end cell) or in a wrapper script:

```python
import uuid
from datetime import datetime, timezone

run_uuid = uuid.uuid4().hex[:8]
run_ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%f")  # microsecond

# Build paths
stem = f"010_FucusDispersal_{experiment_type}_{start_date}_{run_ts}_{run_uuid}"
zarr_path = path_trajectories / str(release_date.year) / f"{stem}.zarr"
nb_path = f"notebooks_executed/TrajectoryCalc/{stem}.ipynb"

# Append to CSV
csv_path = path_trajectories / "runs.csv"
row = ",".join(str(x) for x in [
    run_uuid, run_ts, start_date, experiment_type, velocity_factor,
    max_age_days, calc_dt_mins, particles_per_cell, git_sha,
    zarr_path, nb_path,
])
with open(csv_path, "a") as f:
    f.write(row + "\n")
```

### Querying

```python
import pandas as pd
runs = pd.read_csv("output/Trajectories/runs.csv")

# All surface runs for 2019
surface_2019 = runs.query("experiment_type == 'surface' and start_date.str.startswith('2019')")

# Latest run for each (start_date, experiment_type) combo
latest = runs.sort_values("timestamp").groupby(["start_date", "experiment_type"]).last()

# Compare pre-fix vs post-fix
pre = runs[runs.timestamp < "2026-04-10"]
post = runs[runs.timestamp >= "2026-04-10"]
```

### Pros

- Globally unique — no collisions ever, even parallel identical submissions
- Human-scannable filenames — `ls *010_FucusDispersal*surface*2019*` still works
- Full provenance in a flat CSV — trivial to filter, diff, aggregate
- CSV is append-only — concurrent writers are safe (one line per write)
- Chronological ordering from timestamp — easy to find latest
- Old runs preserved — never overwritten, compare at will

### Cons

- Filenames are long (but tabs complete on the type+date prefix)
- CSV needs a header bootstrap (first run writes it, or pre-create)
- CSV can grow large over thousands of runs (but it's just text)
- Requires discipline to always go through the CSV-writing path
- Papermill output notebook name must be set before execution
  (need to generate UUID/timestamp in the job script, pass as params)

### Who generates the UUID, when, and where

The UUID and timestamp must exist BEFORE papermill starts, because
papermill needs the output notebook path at invocation time. Two
things depend on the UUID: the notebook filename and the zarr path.

**The job script generates both UUID and timestamp**, passes them to
papermill as parameters. The notebook uses them to build the zarr
path and the CSV row. This ensures the notebook filename (set by the
job script) and the zarr path (set inside the notebook) share the
same UUID.

```bash
for experiment_type in surface bottom surface_stokes; do
    run_uuid=$(python3 -c "import uuid; print(uuid.uuid4().hex[:8])")
    run_ts=$(date -u +%Y%m%dT%H%M%S.%N)
    stem="010_FucusDispersal_${experiment_type}_${start_date}_${run_ts}_${run_uuid}"

    srun ... papermill \
        notebooks/010_FucusDispersal.ipynb \
        notebooks_executed/TrajectoryCalc/${stem}.ipynb \
        -p start_date ${start_date} \
        -p experiment_type ${experiment_type} \
        -p run_uuid ${run_uuid} \
        -p run_timestamp ${run_ts} \
        ... &
done
```

The notebook receives `run_uuid` and `run_timestamp` as papermill
params. It uses them to:
1. Build `zarr_path = .../010_FucusDispersal_{type}_{date}_{ts}_{uuid}.zarr`
2. Append a row to `runs.csv` with all params + paths

The notebook never generates its own UUID — it always receives one.
This avoids the UUID being different in the notebook name vs the zarr
name, which would happen if the notebook generated it independently.

**Edge case**: if someone runs the notebook interactively (not via
papermill), `run_uuid` and `run_timestamp` should have sensible
defaults in the parameters cell:

```python
# Parameters (overridden by papermill)
run_uuid = ""
run_timestamp = ""
```

Then in the setup cell:
```python
if not run_uuid:
    import uuid
    from datetime import datetime, timezone
    run_uuid = uuid.uuid4().hex[:8]
    run_timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%f")
```

This way interactive runs still get unique IDs, but papermill runs
use the job-script-generated ones that match the notebook filename.

## Recommendation sketch

Combine B + D + E:
1. **Version tag** in the directory path (Option C lite): `output/Trajectories/v2/`
2. **Zarr attrs** with full params (Option D) for self-describing outputs
3. **JSON sidecar** next to each zarr (Option E) for easy catalog building

This gives: human-readable isolation, machine-queryable metadata, and
a paper trail connecting code version to output.

## Open questions

- How many parameter combinations do we actually have? If it's just
  (date × type × vf), the current naming might be sufficient with
  just a version directory.
- Do we need to compare pre-fix and post-fix results side by side?
  If so, directory isolation (Option C) is essential.
- Should the 003 output (2D fields) also be versioned? They're the
  input to 010 and the grid fix changes them.
