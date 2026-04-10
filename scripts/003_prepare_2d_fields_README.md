# 003_prepare_2d_fields: Grid Registration Hotfix

## What changed

Two fixes applied to `003_prepare_2d_fields.py` to correct how BSH
velocity data maps onto the Parcels NEMO C-grid:

### 1. V-roll: south-face to north-face convention

BSH stores V on the **south face** of each cell (`cellposition=south`
in the NetCDF metadata). Parcels with `gridindexingtype="nemo"` expects
V on the **north face**. The south face of cell j is the north face of
cell j+1, so `V_nemo[j] = V_bsh[j-1]` (a roll by +1 along lat in
descending-lat convention).

Applied via `xarray.DataArray.roll(lat=1, roll_coords=False)` after
land-masking (which uses BSH convention). The boundary row wraps from
the land-masked domain edge, inheriting zero.

### 2. NE F-point coordinates

Parcels NEMO C-grid cell `(yi, xi)` reads data from T-cell
`(yi+1, xi+1)`. The coordinate arrays must be NE F-points
(`lon_T + dlon/2`, `lat_T - dlat/2` where dlat < 0 for descending lat)
so that this index mapping is correct. Previously the output files
stored T-point coordinates, shifting every cell by one index.

### 3. Fine grid crop

Three dead boundaries are removed from the fine grid:
- **North** (lat row 0): V-roll wraps boundary zero here
- **South** (last lat row): land mask appends all-True row for V
- **East** (last lon col): land mask appends all-True column for U

Cropping these (`isel(lat=slice(1,-1), lon=slice(None,-1))`) ensures
particles at the fine grid edges fall through to the coarse grid via
NestedField instead of stalling in a zero-velocity band. The western
boundary has no dead column (no appended mask there).

The coarse grid is not cropped — its dead boundaries are at the overall
domain edge where particles leave anyway.

## Why a hotfix here

The full fix belongs in the bsh-opmod-parcels package, which will be
rewritten from scratch with proper grid handling (see
`grid_reg_bug/PLAN_FOR_BSH_OPMOD_PARCELS.md`). This hotfix applies
the minimum necessary corrections to the 2D field preparation step
so that the fucus-dispersal production runs can proceed with correct
grid registration.

## Impact on downstream

The output 2D files now contain:
- V in NEMO north-face convention (rolled)
- NE F-point lon/lat coordinates (not T-points)
- One fewer lat row for fine grid files

The `010_FucusDispersal` notebook needs no changes — it reads
coordinates from the 2D files via `FieldSet.from_netcdf`, which picks
up the corrected F-point values automatically.

All existing 2D output files must be regenerated before rerunning
particle tracking experiments.

## Validation

See `notebooks/002_ValidateGridFix.md` — 2000 particles across four
regions (Pomeranian Bay, fine N/W/E boundaries) with nested fieldset.
Trajectories cross nest boundaries smoothly, land mask aligns with
coastline, no stagnation bands.
