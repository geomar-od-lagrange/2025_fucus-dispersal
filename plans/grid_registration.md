# Grid Registration: T-points vs F-points in Parcels

## Problem

BSH velocity files store lon/lat at **T-points** (tracer cell
centers). Parcels with `gridindexingtype="nemo"` and
`interp_method="cgrid_velocity"` expects coordinates at **cell
corners** such that cell `(i, j)` spans `lon[i]..lon[i+1]`,
`lat[j]..lat[j+1]`.

The BSH static data includes a separate `lonlat_file_{fine,coarse}.nc`
with F-point coordinates (NE corners): `F = T + dlon/2` in lon,
`F = T + dlat_orig/2` in lat (where `dlat_orig < 0` for descending
lat, so F is shifted south).

The original bsh-opmod-parcels code (`101_run_experiment_fine.ipynb`)
passes F-point coordinates from the lonlat file. Our
`010_FucusDispersal.ipynb` passes T-point coordinates from the
velocity file. Both may be wrong — the exact relationship between
Parcels' NEMO indexing expectations and the BSH coordinate convention
needs verification.

## Evidence

- Particle trajectories extend into land cells (visible in plots)
- Velocity sampled at known wet cell centers returns near-zero
  (~1e-6 m/s) instead of the expected 0.18 m/s, for ALL three
  coordinate options tested (T-points, F-points, SW-corners)
- Coastline polygons built from H0 don't align with trajectory
  termination points regardless of coordinate choice
- The `spatial_interpolation_UV_c_grid` in `parcels.h` reads
  `U0 = data[yi+1, xi]` (NEMO indexing), which shifts the data
  by one row relative to the coordinate cell

## Impact

All production trajectories may be displaced by half a grid cell
(~2.5 km coarse, ~0.5 km fine). Coastal trapping artifacts are
partly caused by this misregistration.

## Next steps

### Step 1: Controlled velocity test

Create a synthetic velocity field with known values in specific
cells. Use a field where cell (i, j) has `U = i` and `V = j` (or
similar cell-ID encoding). Sample with Parcels and verify the
returned values match. This bypasses all physical interpretation
and directly tests the data-to-coordinate mapping.

Do this for all three coordinate options (T, F, SW) and both
`mesh="flat"` and `mesh="spherical"`.

### Step 2: External interpolation comparison

Implement the C-grid interpolation in Python (matching the
`parcels.h` logic) and compare velocity at test points. This gives
ground truth independent of Parcels' JIT compilation.

### Step 3: Fix and validate

Once the correct coordinate convention is identified:
1. Fix `010_FucusDispersal.ipynb` fieldset construction
2. Fix `004_extract_coastline.py` polygon edges
3. Re-run particle test and verify coastline alignment
4. Re-run production trajectories

## Files

- `plans/grid_registration.md` (this file)
- `tmp_corner_fix/test_grid_options.py` — initial diagnostic
- `tmp_corner_fix/check_grid.py` — coordinate comparison
- `tmp_corner_fix/test_curonian2.py` — particle experiment
- `tmp_corner_fix/plot_curonian.py` — trajectory plotting
- Original fieldset code: `~/Downloads/101_run_experiment_fine.ipynb`
