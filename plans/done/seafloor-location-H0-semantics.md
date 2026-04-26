> Implemented as [../../docs/h0_semantics.md](../../docs/h0_semantics.md). Original plan retained as historical record.

# AGENTS.md

Project-specific facts that don't live in the code itself and that are
easy to get wrong. Agents working on this repo should read this first.

## BSH H0 — floor position, not depth

`H0` in the BSH operational-model static files is **minus z of the sea
floor** on a z-up axis with `z = 0` at mean sea level. It is not the
depth of the water column.

- Always-wet cells: `H0 > 0`, and `H0` coincides with depth.
- Tidal-flat cells that dry at low water: the floor can be above MSL,
  so `H0` can be **negative**.

Anywhere you want "depth", filter to always-wet cells first (`H0 > 0`,
or use the always-wet mask produced by
`scripts/004_extract_coastline.py`). Do not treat `H0` as depth
uncritically; mixing dry flats into a mean over a region will pull the
result toward zero or flip its sign.
