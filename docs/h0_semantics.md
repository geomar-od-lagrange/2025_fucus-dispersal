# BSH H0 — floor position, not depth

`H0` in the BSH static files is **minus z of the sea floor** (z-up,
`z = 0` at MSL), not water depth. Always-wet cells have `H0 > 0`;
tidal-flat cells have `H0 ≤ 0` (floor at or above MSL; can still be
wet at high water). Same convention on fine and coarse grids.

Restrict to `H0 > 0` before any depth-shaped reduction. Mixing tidal
flats into a regional mean pulls toward zero or flips its sign in
deep subbasins; a "deepest layer" picker that hits `H0 ≤ 0` cells
lands above the air-sea interface. The wet-cell coastline geojsons
at `data/bsh_hbmnoku_static/coastline_*.geojson` (produced by
`notebooks/004_extract_coastline.py`) carry the same constraint as a
vector geometry.

## Cross-references

- [2d_field_extraction.md](2d_field_extraction.md) — face-level
  analogue (non-NaN AND non-zero) for `load_bottom`'s deepest-fluid
  picker.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) —
  `mean_depth_m` per hex over always-wet cells only.
