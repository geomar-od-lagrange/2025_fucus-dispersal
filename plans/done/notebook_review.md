# Notebook Review — Action Plan

Covers the four main notebooks only (000, 010, 020, 021). Explore notebooks are out of scope.

## Bugs to fix

### B1. `Fucus_location_shp` path broken after move to `data/` — 000, 010, 020

`Fucus_location_shp/` was moved to `data/Fucus_location_shp/` in commit `ec84bde`, but notebooks still reference `.../2025_fucus-dispersal/Fucus_location_shp`.

**Action:** Update all references to `.../2025_fucus-dispersal/data/Fucus_location_shp`.

### B2. `dlon_to_dkm` formula inverted — 020, 021

```python
return dlon / r * degrees_to_radians  # WRONG — gives ~2.8e-6 instead of ~109 km
```

Verified numerically: at equator, 1 degree gives 2.8e-6 instead of ~109 km. The deleted copy notebook had the correct formula: `dlon * r * degrees_to_radians`.

**Action:** Change to `return dlon * r * degrees_to_radians` in both 020 and 021.

### B2b. `earth_radius = 6271` should be `6371` — 020, 021

All notebooks use 6271 km. Mean Earth radius is 6371 km — off by 100 km (~1.6% error). Likely a typo propagated from the original code.

**Action:** Change `earth_radius = 6271` to `earth_radius = 6371` in 020 and 021 (and 000 if any distance code survives D9).

### B3. `Path()` receives `driver` kwarg — TypeError — 000

```python
sh_coast_split.to_file(Path(output_path, "sh_baltic_coast_split.geojson", driver="GeoJSON"))
```

**Action:** `sh_coast_split.to_file(Path(output_path, "sh_baltic_coast_split.geojson"), driver="GeoJSON")`

### B4. `ds_traj_0` undefined — NameError — 020

Referenced in scatter animation and distance statistics cells but never assigned.

**Action:** Determine intent (likely `ds_trajectories`), fix or remove the broken cells.

### B5. `sigma_file_fine` / `sigma_file_coarse` undefined — NameError — 000

The sigma cell reads these but they're only defined in 010. The sigma files provide the vertical coordinate (sigma layers) for the Parcels fieldset — this code belongs in 010, not 000.

**Action:** Remove the sigma cell from 000. It's already properly handled in 010.

### B6. Distance filter threshold inconsistency — 020 vs 021

020 uses `> 10`, 021 uses `> 0.01`. With the broken `dlon_to_dkm`, values were ~1e-6 of true km, so both thresholds are meaningless. After fixing B2, recalibrate to a single sensible threshold (e.g. `> 1` km means "particle actually moved").

**Action:** After fixing B2, unify the threshold in both notebooks.

### B7. `ds_scatter_SH_ani` undefined — NameError — 021

Second scatter animation cell references undefined variable.

**Action:** Remove the orphan cell.

### B8. `extent` should be `sh_extent` in 021 animation

**Action:** Replace `extent` with `sh_extent`.

### B9. Hardcoded obs count `5304` in 020

021 already fixed this to `len(ds_SH.obs)`. 020 still hardcodes `5304`.

**Action:** Replace `5304` with `len(ds_trajectories.obs)` in 020.

### B10. `n_obs = 220 * dt_in_hours * 24` inverted — 021

Should be `220 * 24 // dt_in_hours`. Happens to work for dt=1h but breaks otherwise.

**Action:** Fix to `int(220 * 24 / dt_in_hours)`.

### B11. `coast_region_radius` computed in meters, used as degrees — 000

`coast_region_radius_km * 111e3` gives meters, but `.buffer()` on Geodetic CRS expects degrees.

**Action:** Remove or fix. The whole coast-buffer block in 000 feeds dead code (see D3).

### B12. Inconsistent `Administrative_boundaries` path — 000 vs 021

000 uses `_ags_Administrative_boundaries/`; 021 uses `Administrative_boundaries/`. One must be wrong.

**Action:** Check which path exists on HPC and unify.

## Dead code to remove

### D1. First `gdf_release_area` block in 000

Reads shapefile, filters, computes `area_m2`, `n_release_cells`, `n_particles`. Immediately overwritten by the next cell. Results never used.

**Action:** Remove the first block.

### D2. Parcels imports in 020 and 021

Both import `FieldSet`, `JITParticle`, `ScipyParticle`, `ParticleSet`, `AdvectionRK4`, `AdvectionRK4_3D`, `StatusCode`, and `parcels`. None are used — these notebooks only read zarr data.

**Action:** Remove parcels imports from 020 and 021.

### D3. Unused imports in 000

`parcels` (all of it), `dask`, `dask.distributed.Client`, `cmocean`, `xoak` — none used. Copy-pasted from 010.

**Action:** Strip imports to what 000 actually uses.

### D4. `subbasins_simple` in 021

Expensive GIS buffer operation; result never referenced.

**Action:** Remove.

### D5. `last_modeling_date_str` in 010

Computed but never referenced.

**Action:** Remove.

### D6. Ocean clipping + `fucus_release_area.shp` output in 000

The entire section that clips Fucus polygons to ocean and writes `output/ReleaseArea/fucus_release_area.shp` — never read by any downstream notebook.

**Action:** Remove the ocean clipping section and the `output/ReleaseArea` write.

### D7. Duplicate plot cells in 010

Two byte-identical trajectory plot cells.

**Action:** Remove the duplicate.

### D8. Last cell of 021 — orphan references

References `ds_trajectories.total_distance` (never created) and `cell_code` (commented out).

**Action:** Remove the cell.

### D9. `coast_region_radius` block in 000

Computes a coast buffer using the wrong conversion (B11), feeds into `gdf_baltic_unbuffed` which is plotted once and discarded. Dead exploratory code.

**Action:** Remove the block.

## Keep as-is (intentional)

### Salt/temp fieldset snippets in 010

`salt_temp_subfolder_fine/coarse`, `salt_temp_files_fine/coarse`, and the commented-out `salt_temp_fieldset_fine/coarse` construction. These will be activated soon.

**Action:** Keep. Optionally add a `# TODO: activate for salt/temp scenarios` comment.

## Redundancies to address

### R1. Distance helpers duplicated in 020 and 021

`dlat_to_dkm`, `dlon_to_dkm`, `latlon_dist_to_km`, `earth_radius`, `degrees_to_radians`. Bug B2 must be fixed in both.

**Action (minimum):** Fix B2 in both. **Action (optional):** Extract to a shared `helpers.py`.

### R2. `calc_hist_0` / `calc_hist_1` duplicated in 020 and 021

Identical histogram functions.

**Action (minimum):** Keep in sync. **Action (optional):** Extract to shared module.

### R3. Trajectory loading + subbasin join duplicated in 020 and 021

~20 lines of identical setup. 021 has improvements (dynamic obs count, lower threshold).

**Action (minimum):** Port 021's improvements back to 020 during B6/B9 fixes.

## Enhancements to 010_FucusDispersal

### E1. Clean up parameter interface

Replace the 6 time-related params with 4 cleaner ones:

**Current** (6 params):
```python
release_year = 2019
first_release_month = 1
first_release_day = 1
last_release_month = 1
last_release_day = 29
repeated_release = False
repeated_release_dt_days = 7
max_age_days = 50
```

**Proposed** (4 params):
```python
start_date = "2019-01-01"       # first release date (ISO string, papermill-friendly)
release_period_days = 14        # release window length (0 = single release)
release_interval_days = 7       # days between releases within window
max_age_days = 220              # particle lifetime
```

Derived logic:
```python
first_release_date = datetime.fromisoformat(start_date)
end_release_timeframe = first_release_date + timedelta(days=release_period_days)

if release_period_days == 0 or release_interval_days == 0:
    release_dates = [first_release_date]
else:
    release_dates = [
        first_release_date + timedelta(days=d)
        for d in range(0, release_period_days + 1, release_interval_days)
    ]
```

Example: `release_period_days=14, release_interval_days=7` → releases on day 0, 7, 14 = 3 batches.

Eliminates the redundant `repeated_release` bool, the hardcoded `periods=100`, and the year-boundary limitation.

### E2. In-memory zarr store with final-size chunking

Replace the current disk-write-during-simulation with an in-memory zarr store. Chunks are computed upfront from parameters and set at ParticleFile creation — no rechunking needed.

**Current:**
```python
output_chunks = (10000, int(24 * 60 / output_dt_mins) * 40)
output_particle_file = pset.ParticleFile(
    name=output_path,           # writes directly to disk during simulation
    outputdt=output_dt,
    chunks=output_chunks,       # (10000 traj, 960 obs) — arbitrary
)
pset.execute(...)
```

**Proposed:**
```python
import zarr

# Compute final chunk sizes upfront
n_obs_max = int(max_age_days * 24 * 60 / output_dt_mins)
bytes_per_traj = n_obs_max * 48   # 6 vars (lat,lon,z,time,cell_ID,age_sec) * 8 bytes
traj_chunk = max(1, int(100 * 1024**2 / bytes_per_traj))

# Write to memory during simulation
memory_store = zarr.storage.MemoryStore()
output_particle_file = pset.ParticleFile(
    name=memory_store,
    outputdt=output_dt,
    chunks=(traj_chunk, n_obs_max),
)

pset.execute(
    custom_kernel,
    dt=calc_dt,
    endtime=last_modeling_date,
    output_file=output_particle_file,
    verbose_progress=show_verbose_progress,
)

# Single write from RAM to disk — already correctly chunked
zarr.copy_store(memory_store, zarr.storage.DirectoryStore(output_path))
```

For a typical 220-day hourly run: `n_obs_max=5280`, `bytes_per_traj≈253 KB`, `traj_chunk≈395`.

**Memory budget:** 226k particles × 5280 obs × 48 bytes ≈ 55 GB. Fits in the 120 GB allocation.

Parcels' `ParticleFile` accepts `zarr.MemoryStore` as `name` — confirmed.
