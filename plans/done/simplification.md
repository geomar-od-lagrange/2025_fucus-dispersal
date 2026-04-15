# Simplification Plan

Treat each notebook as a function of its parameters with explicit side effects.

## 010_FucusDispersal

### Kill `is_papermill` and the `TEST_` prefix

`is_papermill` controls the `TEST_` prefix on output filenames. The viz notebooks then skip `TEST_` files via a `_N2*` glob hack. This is hidden coupling — the consumer silently knows to avoid certain outputs.

A quick dev run is just a run with fewer `particles_per_cell` and shorter `max_age_days`. The parameters already distinguish runs. Remove `is_papermill`, always produce the same filename structure.

Also remove the commented-out `# if not is_papermill: gdf_release_area = gdf_release_area[:25]`.

### Kill `repeated_release` boolean

When `repeated_release = False`, the code builds release dates then immediately overwrites them with `[first_release_date]`. Five parameters become dead weight. A single release is just the case where the date range contains only one release interval. Remove the boolean — let the date range speak for itself.

### Replace 6 date params with 2 ISO strings + interval

Current: `release_year`, `first_release_month`, `first_release_day`, `last_release_month`, `last_release_day`, `repeated_release`, `repeated_release_dt_days` (7 params)

Proposed: `start_date`, `release_period_days`, `release_interval_days` (3 params, already planned in E1)

### Simplify output filename

Current: `Fucus_BSH_{dates}_dt{output_dt_mins}min_d{depth}_s{speed}_N{total}_seed{seed}.zarr` — encodes 6 values, viz notebooks must parse via fragile globs.

Proposed: simpler naming, store run metadata inside the zarr attrs.

Feedback: No let's keep this. Proper would be having a separate csv file per run or so for metadata. But out of scope for now.

## 020 and 021 Visualization

### Fix the glob coupling

Both notebooks use:
```python
zarr_files = sorted(output_path.glob(
    f"20??/Fucus_BSH_20??????-20??????_dt60min_d{release_depth}_s{relative_particle_speed}_N2*"
))
```

Problems:
- `_N2*` is a magic filter that happens to skip TEST_ files
- `dt60min` is hardcoded — breaks if output_dt changes
- Encodes assumptions about 010's filename convention

Once `TEST_` prefix is removed from 010, the `_N2*` hack is unnecessary. Replace with a simpler glob or explicit path list.

### Hardcoded `220` days in 021

`n_obs = int(220 * 24 / dt_in_hours)` — the `220` is `max_age_days` from 010, hardcoded here. Should be derived from the data: `n_obs = len(ds_trajectories.obs)` or read from zarr attrs.

### Extract shared code

020 and 021 share ~100 lines of identical logic:
- Trajectory loading, z-filtering, CellID renaming
- Distance calculation (`dlat_to_dkm`, `dlon_to_dkm`, `latlon_dist_to_km`)
- Histogram functions (`calc_hist_0`, `calc_hist_1`)
- Release area loading + subbasin join

Extract to a shared `fucus_utils.py` module. Each viz notebook imports and calls.

Feedback: No keep these for now.

### Split 020 into focused notebooks

020 currently does 5+ things: heatmaps by subbasin, age maps, animations, coast buffer analysis, distance statistics. Split into:
- Heatmap/agemap visualization
- Coast proximity analysis
- Distance statistics

Each imports the shared loading code.

Feedback: Let's keep together for now.

### Make 021's dependency on 000 explicit

021 reads `output/SamplePoints/sh_baltic_coast_split.geojson` produced by 000. This is invisible from 021's parameters. Either:
- Add it as a parameter: `coast_split_path`
- Or inline the coast construction in 021 (it's ~20 lines)

Feedback: Keep as is for now.

## 000_FucusStartLocations

### Three interleaved concerns

000 does: (a) re-project shapefiles to geojson, (b) construct SH coast segmentation with hardcoded split lines, (c) exploratory Fucus location plotting. These are separate tasks.

The exploratory plotting at the end serves no pipeline purpose — candidate for `explore/`.

## Priority order

1. Kill `is_papermill` + `TEST_` prefix (unblocks glob simplification)
2. Simplify date params (E1, already planned)
3. Fix globs in 020/021 (follows from #1)
4. Extract `fucus_utils.py` (reduces divergence risk)
5. Split 020 (can be done incrementally)
