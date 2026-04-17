# Rollout Plan

Overall plan for rolling out all changes. Ordered by dependency.

## Phase 1: Cleanup PR (current branch `wr/cleanup`)

What's already done or in progress. Commit and merge.

### Done
- Notebook renames (000, 010, 020, 021) + VoronoiPolygons → explore/ + copy deleted
- Bug fixes: `dlon_to_dkm`, `earth_radius`, `Path(driver=)`, undefined variables, broken paths
- Dead code removal across all notebooks
- `is_papermill` + `TEST_` prefix removed from 010
- `repeated_release` logic removed — single release per run
- Date params simplified to `start_date` (ISO string)
- `base_path` parameterized in all notebooks (papermill-tagged parameters cells)
- All paths derived from `base_path` / `repo_path`
- Job script cleaned up (broken path, dead scaffolding, generic name)
- Job script updated to match new notebook params
- `001_prepare_sigma_files.py` created
- `Fucus_StartLocations.sh` job script created
- Globs in 020/021 fixed (no more `_N2*` hack, no hardcoded `dt60min`)
- `calc_dt_mins` set to 5 (from 15)
- Aging kernel commented out (single release, `endtime` handles lifetime)

### Still to do before merge
- Clean up import comments from HPC debug session (merge numpy.random, drop glob, etc.)
- Remove 010 analysis/vis cells → they move to 020_Diagnostics in Phase 2
- Remove `output/SamplePoints/` tracked files from git

## Phase 2: Restructure visualization notebooks

Split current 020 and 021 into focused notebooks. No Dask — all `.compute()` calls are just triggering lazy xarray evaluation on the default synchronous scheduler.

### 010_FucusDispersal
- Remove all cells after `# Execute Simulation` (analysis + plots)
- Pure function: params in, zarr out

### 020_Diagnostics
- Quick per-run sanity checks (from 010's current analysis section)
- Dead/moving particle classification
- Trajectory heatmap, first/last position maps
- Dead cell count
- Runs per-experiment

### 021_Heatmaps
- Particle count + mean age heatmaps by HELCOM subbasin (4×4 grids)
- From current 020 Maps section
- Aggregates across multiple runs

### 022_HeatmapsSH
- Same as 021 but Schleswig-Holstein regional scale
- Coast segment breakdown (Flensburger Förde, Kieler Bucht, etc.)
- Monthly breakdown
- From current 021 Maps + Regional sections

### 023_Animations
- All animations (from both current 020 and 021)
- Heatmap animations, scatter animations
- Aggregates across runs

### 024_CoastAnalysis
- Near-shore / beaching analysis
- Currently broken in 020 (depends on deleted `gdf_first_step`)
- Will be rewritten to use `nearshore_time` from trajectory output (Phase 4)
- Placeholder until Stokes drift is implemented

### 025_DistanceStats
- Distance quantiles and means by subbasin + coast segment
- Line plots, bar plots, boxplots
- Merged from current 020 and 021 distance sections

## Phase 3: 2D field extraction

Prereq for I/O reduction and Stokes integration. See `plans/2d_field_extraction.md`.

- Create `002_prepare_2d_fields.py`
  - Extract BSH surface layer (sigma=0) as 2D C-grid files
  - Extract BSH bottom layer (deepest valid sigma per cell) as 2D C-grid files
- Simplify 010 kernel: `AdvectionRK4_2D_SIGMABSH` → `AdvectionRK4_2D_BSH` (drop sigma, hardcode depth=0)
- Switch fieldset from 3D nested to 2D nested
- **25x I/O reduction**: ~3 TB/year → ~120 GB/year
- Same notebook for surface and bottom — experiment type is which 2D files are pointed at

## Phase 4: In-memory zarr store

See `plans/notebook_review.md` section E2.

- Parcels writes to `zarr.MemoryStore()` during simulation (already prototyped in HPC session)
- Small chunks during simulation for efficient appending
- After `pset.execute()`, rechunk and write to disk: full obs length, ~100 MB traj chunks
- Memory budget: 87,200 trajectories × 5,280 obs × 48 bytes ≈ 22 GB (fits in 120 GB)

## Phase 5: Parallel job execution

730 independent single-release experiments (10 years × every 5 days × 100 particles/cell).

- Job script generates srun commands in a loop over start dates
- `xargs -P N` for parallelism within a SLURM allocation
- Or SLURM job arrays (one job per start date)
- Each run: ~6 hours at 100 particles/cell with 2D fields
- Aggregate walltime: 730 × 6h. At 100 parallel nodes: ~44 hours

## Phase 6: Stokes drift

See `plans/stokes_drift.md`.

- Download BALTICSEA_MULTIYEAR_WAV_003_015 from CMEMS (~240 GB for 2016–2025)
- Extend `002_prepare_2d_fields.py`:
  - Fill Stokes NaN with zero (sheltered waters = no wave drift, physically correct)
  - Interpolate VSDX/VSDY onto BSH C-grid U/V points
  - Zero at land-adjacent edges
  - Sum with BSH surface: `U_total = U_surface + U_stokes`
- No simulation code changes — same kernel, same notebook, different input files
- Create static coast distance field from H0 bathymetry
- Add `nearshore_time` tracking kernel to 010
- Implement stochastic beaching weights in 024_CoastAnalysis

## Phase 7: Salt/temp fieldset (future)

Salt/temp snippets preserved in 010 for this. Activate when ready:
- Uncomment `salt_temp_variable_ID`, file loading, fieldset creation
- Add S/T as nested fields
- Track as particle variables
- Requires 3D fieldset again (or 2D extraction at the relevant sigma level)
