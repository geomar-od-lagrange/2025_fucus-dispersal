# Next Steps

## Priority 1: Grid registration (blocking)

The BSH velocity data uses T-point (cell center) coordinates but
Parcels NEMO C-grid expects cell corner coordinates. All production
trajectories may be shifted by half a grid cell (~2.5 km coarse,
~0.5 km fine). See `grid_registration.md`.

**Action**: Build a synthetic velocity field with known per-cell
values, verify Parcels samples the correct cell, identify the
right coordinate convention. Then fix `010_FucusDispersal.ipynb`
and rerun.

## Priority 2: Corner rounding (ε-floor)

Once grid registration is correct, implement the ε-floor on land
edges in `parcels.h`. See `corner_rounding.md` and
`corner_theory.md` for the full derivation.

**Action**: Patch `parcels.h`, validate with the scipy prototype
and Parcels test runs at the Polish coast, then integrate into
production.

## Priority 3: Coastline polygons

The `004_extract_coastline.py` script builds cell edges from H0
T-points. After the grid registration fix, verify these edges
match the Parcels cell boundaries. May need adjustment.

**Action**: Re-run coastline extraction, overlay on corrected
trajectories, confirm alignment.

## Priority 4: Production reruns

After fixing grid registration and implementing corner rounding:
1. Rerun `010_FucusDispersal.ipynb` for all release dates
2. Rerun `020_FucusHeatmaps.ipynb` to check for reduced coastal
   artifacts
3. Compare with previous results to assess impact

## Independent: Documentation for Parcels v4

The ε-floor approach should be documented in a form suitable for
upstreaming to Parcels. The `corner_vortex.md` jupytext notebook
contains the full derivation and can serve as the basis for a
technical note or PR description. Not blocking for this project
but important for the community.

## Working artifacts

`tmp_corner_fix/` contains the pixi env (Parcels 3.0.6),
diagnostic scripts, and the executed `corner_vortex.ipynb`. To
resume on another machine: `cd tmp_corner_fix && pixi install`.
