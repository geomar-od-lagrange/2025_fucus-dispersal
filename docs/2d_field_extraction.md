# 2D field extraction from BSH HBMnoku

`notebooks/003_prepare_2d_fields.py` writes one preprocessed netCDF
per input c-file per variant. The Parcels runs in
`notebooks/010_FucusDispersal.md` consume these 2D fields directly —
sigma-layer selection happens here, not at runtime in the kernel.

| Variant          | U/V content                                | Consumer regime          |
|------------------|--------------------------------------------|--------------------------|
| `surface`        | sigma layer 1                              | `regime="surface"`       |
| `bottom`         | deepest fluid layer per face               | `regime="bottom"`        |
| `surface_stokes` | surface + interpolated CMEMS Stokes drift  | `regime="surface_stokes"`|

Output globs as `c_file_{fine,coarse}_<YYYYMMDDHH>_<i>_<j>_<variant>.nc`
under `output_root/2d_fields/`. 010 picks a variant by changing the
`regime` parameter; the fieldset glob is
`c_file_{fine,coarse}_*_{regime}.nc`.

## Why preprocess to 2D

BSH c-files carry 25/36 sigma layers (fine/coarse) per timestep —
~3 TB/year. Each Parcels run only needs one 2D layer per particle.
Preprocessing once brings I/O down to ~120 GB/year (≈25× reduction)
and lets every regime in the sweep share the same disk-cached field.
Earlier 3D-sigma failure modes (C-grid corner trapping, T-/F-point
grid-registration ambiguity) are archived under
`plans/done/corner_*.md` and `plans/done/grid_registration*.md`; their
fixes were never wired in because 2D preprocessing removed the
conditions that triggered them.

## Bottom: face-by-face deepest fluid

A naive "deepest non-NaN" picker fails on BSH because BSH stores
`0.0` (not NaN) at C-grid U/V faces blocked by a shallower
neighbour's bathymetry — these are no-slip walls, not fluid.
`load_bottom` instead takes, for each face independently, the deepest
layer that is **both non-NaN and non-zero**. U and V are picked
separately because their face bathymetries —
`min(bathy_self, bathy_east)` for U, `min(bathy_self, bathy_south)`
for V — can differ in a single cell. The deepest-fluid index is
computed at `t=0` and applied to every timestep.

## C-grid hygiene (all variants)

After layer selection, two transforms apply uniformly:

1. **Edge land mask** — zero `U` at the eastern face of any cell
   adjacent to land in either neighbour, `V` at the southern face
   analogously.
2. **V-axis roll to NEMO** — BSH stores `V` on the south face;
   `gridindexingtype="nemo"` expects the north face. With descending
   lat, shift `V` by `+1` along lat (`v.roll(lat=1, …)`); the
   wrapped boundary row picks up the land-masked southernmost row.

Output coordinates shift to NE F-points
(`lon_f = lon + dlon/2`, `lat_f = lat - dlat/2`). Fine-grid files crop
three boundary rows/columns after hygiene (lat row 0: V-roll boundary
zero; lat row -1 and lon col -1: edge land mask) so particles fall
through to the coarse grid via `NestedField` rather than stalling.

## Cross-references

- [stokes_drift.md](stokes_drift.md) — `surface_stokes` recipe.
- [h0_semantics.md](h0_semantics.md) — bathymetry sign convention.
- [seeding.md](seeding.md) — how 010 selects the variant per run.
