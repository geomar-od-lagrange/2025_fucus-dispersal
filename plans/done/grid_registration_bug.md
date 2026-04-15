# Grid Registration Bug: BSH T-points vs Parcels F-points

## Discovery

While investigating C-grid corner trapping, we found that the BSH
velocity data is fed to Parcels with wrong coordinate interpretation.

## The problem

BSH velocity files have lon/lat coordinates at **T-points** (tracer
cell centers). The separate `lonlat_file_{fine,coarse}.nc` files
have coordinates at **F-points** (cell corners, = T + dlon/2 in
NEMO convention, NE corner of each T-cell).

Parcels with `gridindexingtype="nemo"` and `interp_method="cgrid_velocity"`
interprets the passed lon/lat as cell corners defining cell extents:
cell `i` spans `lon[i]..lon[i+1]`.

The original bsh-opmod-parcels notebook
(`101_run_experiment_fine.ipynb`) passes F-point coordinates from
the lonlat file. Since F[i] = T[i] + dlon/2 (the NE/eastern corner
of T-cell i), Parcels cell `i` spans F[i]..F[i+1], which
corresponds to **T-cell i+1**, not T-cell i. Everything is shifted
by one cell index.

The 010_FucusDispersal notebook in this repo made it worse by
passing T-point coordinates directly (from the velocity file),
adding another half-cell shift on top.

## Evidence

- Particle trajectories extend into land cells (visible in plots)
- Particles seeded on land get nonzero velocity (interpolated from
  neighboring wet cells due to the shift)
- The stuck-particle locations don't align with the coastline polygon
  (which was built from H0 with correct T-point cell extents)
- Velocity at the last wet T-cell (H0=10.2m) is read by Parcels as
  belonging to the adjacent land cell

## Impact

All production trajectories have been computed with misregistered
grid coordinates. The spatial error is half a grid cell:
- Coarse grid: ~2.5 km in lon, ~2.8 km in lat
- Fine grid: ~0.5 km in both

The coastal trapping we investigated is partly caused by this
misregistration (particles in land cells seeing blended velocity)
and partly by the genuine C-grid corner stagnation (which exists
regardless of registration).

## Fix

The correct coordinates for Parcels NEMO C-grid are the **SW
corners** of each T-cell: `lon_sw = T_lon - dlon/2`,
`lat_sw = T_lat - dlat/2`. These need N+1 values for N cells
(the last point is the eastern/northern boundary).

For the BSH grid: `lon_sw[i] = T[i] - dlon/2` for i=0..N, with
`lon_sw[N] = T[N-1] + dlon/2`. Same for lat.

This needs to be fixed in:
1. `notebooks/010_FucusDispersal.ipynb` — fieldset construction
2. `scripts/004_extract_coastline.py` — coastline polygon edges
3. `notebooks/001_ShowModelCoastline.ipynb` — visualization

## Status

Discovered 2026-04-09. Not yet fixed. Needs verification with a
controlled test (known velocity field, check particle positions
match expected streamlines).
