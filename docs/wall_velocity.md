# Wall velocity and land-step reject (010)

The BSH C-grid closes the coast with zero-flux faces. Bilinear C-grid
interpolation drives the wall-normal velocity linearly to zero at such a
face, so onshore forcing pulls drifters exponentially close to it (~1e-12
cell), where float64 rounding freezes them for good. 010 therefore adds a
small **wall velocity ε pointing away from land** on wall faces and drops
any step that would end on land. Diagnosis, prototype and ε sweep:
[numerical_coastal_trapping.md](../plans/numerical_coastal_trapping.md).

| piece | where |
|---|---|
| wall + land fields per nest × regime | `notebooks/005_BuildWallFields.py`, `scripts/005_BuildWallFields_job.sh` → `output_root/wall_fields/walls_{fine,coarse}_<regime>.nc` |
| kernel, parameters | `notebooks/010_FucusDispersal.py`: `wall_velocity_mm_s = 2.0`, `reject_land_steps = True`, `wall_fields = "wall_fields"` |
| per-trajectory outputs | `land_seed` (once), `nrej` (running count of rejected steps, every output step; 0 for land seeds) |
| real-drifter set | `land_seed == 0` in 021–024b, 024d, 036 |

## Wall set and sign (005)

- **Wall face**: velocity 0 or NaN (Parcels fills NaN with 0) at every
  15-min step of every `c_file` 010 loads for the regime (fine ∩ coarse
  stems, 2016–2025). A face that carries flow at any step is not a wall.
- **Wet T-cell**: at least one non-wall face. `land` = not wet.
- **Sign**: wall face between a wet and a non-wet cell → ±1 into the wet
  cell. Thin wall (both sides wet; sills between shallow cells that flood
  from opposite sides) → into the deeper cell by H0, tie → into the cell
  with more open faces, second tie → 0. No wet side → 0.
- Per regime: fine bottom and fine surface_stokes differ (the bottom files
  close some faces the surface files do not); coarse sets agree.
  `bottom_ub49` reads the bottom walls (walls follow `regime`).

The fields share lon/lat with the `c_file`s, so Parcels builds identical
grids and NestedField picks the same nest for `UV`, `UVwall` and `Land`.

## Kernel (010)

`AdvectionRK4WallFriction` is `parcels.AdvectionRK4` with `UV + UVwall`
sampled at every stage (`f·UV + UVwall` with bed friction, `f = 1` at the
default `U_b = 0`; [bottom_friction.md](bottom_friction.md)); `UVwall` is
a NestedField of C-grid VectorFields built like `UV`, and ε (m/s) is the wall
fields' scaling factor. In a wall cell the
wall-normal velocity becomes `u = (1 − ξ)·ε + ξ·U_open`, so under onshore
`U_open` drifters rest at `ξ* = ε / (ε + |U_open|)` instead of 1e-12, stay
resolvable, and leave when the forcing relaxes.

Land seeds (below) are not advected: they start in a cell that never
carries flow, where plain RK4 returns exactly zero, so skipping them changes
no position at ε = 0 and saves ~35 % of the field samples.

With `reject_land_steps`, a step whose end position lies in a land cell
(`Land`, `cgrid_tracer` = constant per cell) is not taken; the particle keeps
its position and `nrej` increments. RK4 stages can still sample land cells;
only the end position is guarded.

`wall_velocity_mm_s = 0` with `reject_land_steps = False` adds exact zeros
and reproduces plain `AdvectionRK4` bit for bit (checked against production
zarrs, plan §Fix).

## ε = 2 mm/s

- Bounds (plan §Fix, "ε bounds"): below ~1 mm/s the int16 velocity storage
  (0.5 mm/s) dominates; above ~2.5–3 mm/s the wall drift over one forcing
  correlation time (~1 d) exceeds ¼ fine cell and bottom rest points move
  past 0.1 cell.
- Constant in m/s on both nests: the same fractional rest point ξ* on the
  0.9 km fine and 5 km coarse cells. ε ∝ Δ would equalise the relaxation
  time but push the coarse rest point ~6× further off.
- Added after the bed friction of `bottom_ub49`, not subject to it: ε is a
  numerical device, not current the drifter feels. The rest point under
  sliding onshore flow becomes `ε / (ε + f·|U_open|)`
  ([bottom_friction.md](bottom_friction.md)).
- Sensitivity: ε = 0 → ε > 0 is a large systematic shift (no pinning, more
  stranding); within ε = 1–3 mm/s results differ beyond the sampling floor,
  non-monotonically, while a 1e-4 relative change of ε or ≤ 1 m release
  jitter stays at the floor (plan §Fix, "Sweep result"). The authors fixed
  ε = 2 mm/s; the paper carries no ε sensitivity.

## Land seeds

A particle is land-seeded when its release position lies in a land cell of
the regime's land set. 010 samples `Land` through the NestedField at
release (same nest and cell as the kernel) and writes `land_seed` once per
trajectory. Considered keeping the zero-first-step-displacement test,
rejected: ε moves land seeds next to a wall face, and a seed in a cell that
is dry at release but floods later is a water seed.

## Considered, rejected

- Post-process stalled drifters as stranded: leaves pinned positions in
  every occupancy store; a methods decision, not a numerical fix.
- Random walk with reflecting walls: adds physics (diffusivity) the study
  does not otherwise carry.
- Zero onshore Stokes in the whole wall cell: moves the pin one cell out.
- Smaller `dt` / AdvectionAnalytical: dt only shrinks land entries (5×
  cost); analytical advection is Scipy-only in Parcels 3.0.6, no NestedField.
- Wall set from a monthly sample: marks intermittently flooded cells as land,
  and drifters entering them while wet then exit only by ε.
- Reject step without the full-period wall set: holds drifters at the edge
  of flooding cells.

Cost: +19 % 010 runtime (prototype); 005 reads ~2.2 TB of `c_file`s per
regime.
