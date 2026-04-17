# 010 Notebook Cleanup Plan

TODOs from review, grouped by theme.

## A. Simplify file loading — drop time filtering, keep timestamps

**TODOs**: cells `554c1145`, `e3f2c02c`, `5a762fcf`

The current `years` range, date-string filtering in `get_file_list`, and year-based subfolder logic are unnecessary complexity — just load all available 2D files.

**Why timestamps are essential**: Without explicit `timestamps`, `FieldSet.from_netcdf` falls back to xarray which opens every file to read time coords — minutes to hours for large datasets. With `timestamps`, fieldset construction is < 1 second. But filename-derived timestamps could be wrong (parsing bug, renamed file), and Parcels won't catch this — it only raises on structural mismatches (wrong shape etc.), not on wrong times.

**Action**:
1. Remove `years` variable entirely
2. `get_file_list` becomes a simple sorted glob — no date filtering
3. Construct FieldSet from *all* 2D files for fine and coarse grids
4. Keep `get_timestamp_from_file` and the `timestamps` list
5. **Add plausibility check**: after building the file list and timestamps, open the *first* and *last* file only, read their actual time coordinate, and assert that the first/last derived timestamps match to within 1 minute. Raise with a clear message if not. This catches filename-parsing bugs without paying the cost of reading all files.

## B. Path management

**TODO**: cell `488396bc`

Paths are scattered and assume specific HPC directory layouts (`/gxfs_work/geomar/smomw122/...`, `/gxfs_work/geomar/smomw400/...`). The notebook hardcodes where large data files happen to live.

**Inventory** (what 010 actually needs):

| Path | What | Notes |
|---|---|---|
| `base_path` | Project root | Parameter already |
| `base_path/output/2d_fields/` | 2D velocity fields | Derived from `base_path` |
| `base_path/output/Trajectories/{year}/` | Output zarr | Derived from `base_path` |
| `base_path/data/Fucus_location_shp/` | Release locations geojson | In repo, derived from `base_path` |
| `/gxfs_work/.../bsh_operationalmodel_data/` | H0 bathymetry for analysis plots | Only used in analysis cells |

The two hardcoded `/gxfs_work/` paths (`smomw122` for static files, `smomw400` for raw c_files) are **not needed by 010** once we strip analysis cells (plan E2.7). They're only needed by preprocessing scripts and job files, which already set them as shell variables.

**Action**: Make all input/output directory paths explicit parameters — no derivation from `base_path`. The 010 notebook reads and writes specific directories, and the user must be able to override any of them independently.

Parameters to add (all full directory paths, set by user):
- `path_2d_fields` — where 2D velocity fields live (input)
- `path_release_locations` — where the Fucus geojson lives (input)
- `path_trajectories` — where output zarr files are written (output)

Remove `base_path` / `repo_path` — they only existed to derive the above. Each path is set explicitly in the parameters cell with sensible defaults.

## C. Release logic

### C1. `relative_position_in_cell` (cell `7524c65c`)

**TODO**: Review logic, add docstring.

The function maps `(x_rel, y_rel) ∈ [0,1]²` to a point inside a HELCOM grid cell (quadrilateral). It treats the cell as a parallelogram using two edge vectors from vertex 0. This is a standard bilinear mapping for quads — correct for parallelograms, approximate for general quadrilaterals (misses the cross term). HELCOM cells are regular enough that this is fine.

**Action**: Add a one-line google-style docstring: `"""Map unit-square coordinates to a point inside a quadrilateral cell."""`

### C2. Exact N per cell (cells `079b1447`)

**TODO**: Current code draws `n_total_particles` with uniform random cell selection — cells get approximately `particles_per_cell`, not exactly.

With N=100, random cell selection gives noticeable coverage variability. Since we're running many small local connectivity experiments, low coverage in one region can't be compensated by high coverage elsewhere. We need exact N per cell.

**Action**: Change release logic to exact N per cell with independent random positions:
```python
cell_indices = np.repeat(np.arange(n_release_cells), particles_per_cell)
rand_x = np.random.uniform(size=n_total_particles)
rand_y = np.random.uniform(size=n_total_particles)
```
`np.repeat` assigns cells deterministically (exactly N per cell). `rand_x`/`rand_y` are `n_total_particles` independent draws — no synchronization across cells. Keep `particles_per_cell` parameter name — it now means exactly what it says.

## D. Kernel correctness (cell `4f04e057`)

**TODO**: Verify this is proper RK4.

**Analysis**: The kernel implements standard RK4 for the ODE `dx/dt = u(x,t)`:
- k1 = f(t₀, x₀)
- k2 = f(t₀ + dt/2, x₀ + k1·dt/2)
- k3 = f(t₀ + dt/2, x₀ + k2·dt/2)
- k4 = f(t₀ + dt, x₀ + k3·dt)
- x₁ = x₀ + (k1 + 2k2 + 2k3 + 2k4)/6 · dt

The code matches this exactly. The `particle_dlon +=` / `particle_dlat +=` pattern applies the displacement (scaled by `relative_particle_speed`) to the particle's position via Parcels' displacement convention.

**Verdict**: Correct RK4. The depth argument is hardcoded to `0` which is correct for 2D fields (Parcels still expects the depth dimension in the field lookup). No action needed beyond removing the TODO.

One subtlety: velocities are in deg/s (Parcels convention for C-grid with lon/lat coords), so the RK4 stages correctly add `u·dt` to lon and `v·dt` to lat in degree-space. This is a flat-Earth approximation, acceptable at Baltic scales.

## E. Notebook structure and cleanup

**TODOs**: cells `32c6537d`, `f5bbdace`, `5c8cbe43`, `f781690e`, `d786db87`, `193c5ef9`, `e95be815`, `a668b493`

### E1. Proposed cell order

The notebook should follow a linear flow:

```
# Parameters                    (papermill cell)
# Setup                         (derived params, paths, output path)
# Release locations             (load shapefile, sample positions)
# Fieldset                      (load 2D files, build nested fieldset)
# Kernel + ParticleSet          (define kernel, create pset)
# Output                        (MemoryStore, ParticleFile)
# Execute                       (pset.execute)
# Write to disk                 (memstore → zarr on disk)
# Minimal diagnostics           (n valid trajs, time range — that's it)
```

### E2. Specific actions

1. **Rename `aging_kernel` → `max_age_kernel`** — it deletes particles past max age, not "aging"
2. **Move kernel defs** to after the fieldset section, just before `pset.execute`
3. **Inline kernel list**: pass `[AdvectionRK4_2D_BSH, max_age_kernel]` directly to `pset.execute`, drop `custom_kernel` variable
4. **Remove commented temp/sal code** from particle variables cell — clean dead code
5. **Move `from zarr.storage import MemoryStore`** to the imports cell
6. **Clean up output path logic**: define `output_path` (disk path) early in setup. Create `memory_store = MemoryStore()` just before ParticleFile. After execute, write `memory_store` → `output_path`.
7. **Strip analysis cells**: everything after "Write to disk" should be just: count valid trajectories, print time range. All maps/histograms belong in 020+.
8. **Run `black` on notebooks** (via `black --ipynb notebooks/010*.ipynb` or `nbqa black`)

## F. Chunking strategy

**TODO**: cell `d08a46f2`

Single release time per run → all particles synchronous → obs dimension aligned → free to choose any chunk shape without Parcels rewriting internally ragged chunks.

Sizes (full year, 87,200 particles, hourly output):
- ~700 kB per variable per obs step (all trajectories)
- ~70 kB per variable per trajectory (all obs)
- ~24.4 GB total (4 vars × float64)

**Action**: Add `chunk_traj` and `chunk_obs` as parameters (defaults: 10,000 and 1,000). Last chunks may be smaller — that's fine. Use `None` for MemoryStore during simulation (let Parcels handle in-memory layout), then rechunk to `(chunk_traj, chunk_obs)` on the memstore→disk write.

## Implementation order

1. **B — Explicit path parameters** (sets up the structure everything else depends on)
2. **A — Simplify file loading** (drop date filtering, add timestamp plausibility check)
3. **C2 — Exact N per cell** (change release logic)
4. **C1 — Add docstring** (trivial)
5. **E2.5 — Move zarr import**
6. **E2.6 — Clean up output path** (disk path early, MemoryStore before ParticleFile)
7. **F — Chunking parameters** (add `chunk_traj`, `chunk_obs` to params, rechunk on disk write)
8. **E2.1–E2.4 — Kernel rename, reorder, inline, remove dead code**
9. **E2.7 — Strip analysis cells** (move to 020+)
10. **E2.8 — Run black**
11. **Remove all TODO comments from notebook** (final step)
