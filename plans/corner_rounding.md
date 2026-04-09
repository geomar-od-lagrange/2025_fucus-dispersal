# Corner Rounding for C-Grid Velocity Fields

## Problem

On a NEMO-style C-grid, velocity components at edges adjacent to land
are zeroed. At concave corners (two land-adjacent edges meeting), both
U and V go to zero. The C-grid interpolation

    U = (1 - xsi) * U0 + xsi * U1
    V = (1 - eta) * V0 + eta * V1

creates a stagnation point at each concave corner. Streamlines converge
asymptotically toward the corner. With finite timesteps, particles
ratchet across streamlines due to truncation error and get permanently
trapped.

Visible everywhere the staircase coastline has concave corners, e.g.
the Curonian Spit.

## Physics

The stream function in a corner cell (U_E = 0, V_S = 0 for a SE
corner, with continuity V_N = U_W) is:

    ψ(x, y) = (1 - x) · y

Streamlines are hyperbolas ψ = K. At the corner, ψ = 0. The ψ = 0
"streamline" is the union of the two land-adjacent edges — once a
particle gets close to ψ = 0, the velocity is near-zero and it's
trapped.

Key insight: in continuous time, a particle stays on its streamline
forever — no trapping. The problem is purely a finite-timestep
artifact. We need a soft wall that prevents ψ from dropping too low.

## Approach: ψ-floor via modified C-grid interpolation

Patch the Parcels C-grid interpolation (`spatial_interpolation_UV_c_grid`
in `parcels.h`) to detect concave corners from the velocity data and
apply a ∇ψ nudge when the particle is in the low-ψ zone.

### Corner detection (automatic, no precomputed flags)

A concave SE corner exists when U1 = 0 and V0 = 0 (both edges toward
the corner are land-masked). Same logic for the other three corners:

| Corner | Condition   | ψ formula     | ∇ψ (cell coords)    |
|--------|-------------|---------------|----------------------|
| SE     | U1=0, V0=0  | (1-xsi)·eta   | (-eta, 1-xsi)       |
| SW     | U0=0, V0=0  | xsi·eta       | (eta, xsi)           |
| NE     | U1=0, V1=0  | (1-xsi)·(1-eta) | (-(1-eta), -(1-xsi)) |
| NW     | U0=0, V1=0  | xsi·(1-eta)   | (1-eta, -xsi)        |

No precomputed corner flags needed — the zero edge values in the
velocity data are sufficient.

### ψ_min value

With r = 0.25 (quarter cell width), ψ_min = r²/4 = 0.015625.

The ψ_min hyperbola's symmetric point (where 1-xsi = eta) is at
(0.875, 0.125). This keeps particles within the inner 87.5% of the
cell, comparable to the quarter-circle arc's exclusion zone.

### Nudge formula

Inside `spatial_interpolation_UV_c_grid`, after computing U and V
(in cell-local flux space), before the Jacobian transform:

    if (is_zero(U1) && is_zero(V0)) {       /* SE corner */
        double psi = (1.0 - xsi) * eta;
        if (psi < PSI_MIN) {
            double deficit = PSI_MIN - psi;
            U += NUDGE_C * deficit * (-eta);        /* ∂ψ/∂xsi */
            V += NUDGE_C * deficit * (1.0 - xsi);   /* ∂ψ/∂eta */
        }
    }
    /* ... repeat for SW, NE, NW */

The nudge is in cell-local flux space. The existing Jacobian transform
converts it to physical (u, v) automatically. The nudge strength is
proportional to how far below ψ_min the particle sits — a soft wall,
not a hard projection.

NUDGE_C controls the wall stiffness. It should be large enough that a
particle at ψ ≈ 0 gets pushed back to ψ_min within a few timesteps.
Order of magnitude: NUDGE_C ≈ 1 / PSI_MIN ≈ 64. To be calibrated in
the test notebook.

### Multiple corners in one cell

A cell can have at most two concave corners (e.g., a peninsula tip).
Each corner's nudge is independent — they act on different parts of
the cell and their ψ formulas don't interfere.

## Implementation plan

### Step 0: scipy prototype (notebook 002)

`notebooks/002_CornerRounding.ipynb`

Before touching Parcels, validate the ψ-floor approach with a pure
scipy/numpy particle integrator. This runs locally on the existing
min_data.

Cells:
1. Load coarse H0 and surface velocity for one timestep
2. Zoom to Curonian Spit (~20–22°E, 54.5–56°N), identify concave
   corners from the velocity data
3. Implement C-grid interpolation in Python (u = (1-xsi)*U0 + xsi*U1,
   v = (1-eta)*V0 + eta*V1) with and without the ψ-floor nudge
4. Plot streamlines in a single corner cell (as in the existing
   `plans/corner_streamlines.png`) — original vs nudged
5. Integrate a handful of trajectories (simple RK4 with scipy or
   manual loop) through the corner region — original vs nudged
6. Show that nudged particles round the corner instead of getting stuck
7. Sweep NUDGE_C values, plot sensitivity — find the range where
   particles escape without oscillating

Data needed (all in min_data/):
- `bsh_operationalmodel_data/static_file_coarse/H0_file_coarse.nc`
- `output_2d/c_file_coarse_2020_surface.nc` (24 timesteps, 15 min)

### Step 1: patch parcels.h

Create `custom_parcels_include/` in the repo with copies of all 4
Parcels 3.0.6 header files from the container image. Patch
`parcels.h`:

1. Add `#define PSI_MIN 0.015625` and `#define NUDGE_C <value from step 0>`
   at the top
2. In `spatial_interpolation_UV_c_grid` (after line ~608 where U and V
   are computed), add the four corner checks with the nudge formula
3. Keep the patch minimal — only touch the interpolation function,
   don't change signatures or data flow

Runtime swap in the notebook/script:

    import parcels.tools.global_statics as _gs
    _gs.get_package_dir = lambda: Path("custom_parcels_include").resolve().parent

This works in the read-only Singularity container because the headers
are read from the bind-mounted repo directory, not the container
filesystem. JIT compilation writes only to the (writable) cache dir.

### Step 2: Parcels validation (notebook 002, continued)

Add cells to `002_CornerRounding.ipynb`:

1. Set up a minimal Parcels FieldSet from the coarse surface velocity
   (one timestep, flat mesh for simplicity)
2. Release ~50 particles near the Curonian Spit concave corners
3. Run with original `parcels.h` — show particles getting stuck
4. Swap to patched `parcels.h` — show particles rounding corners
5. Compare trajectories side by side
6. Verify: particles far from corners are unaffected (diff should be
   zero for particles in open water)

### Step 3: multi-timestep test

Extend the Parcels test to use all 24 available timesteps (6 hours).
Check that:
- The nudge works correctly when flow reverses (tide)
- No particles are ejected to unreasonable positions
- Energy/distance statistics are reasonable

### Step 4: production integration

1. Copy the patched `custom_parcels_include/` into the production
   repo on the HPC
2. Add the `get_package_dir` monkey-patch to `notebooks/010_FucusDispersal.ipynb`
   (in the import/setup section)
3. Re-run a single release date as a test
4. Compare heatmaps (020) with and without the patch
5. If good, re-run all release dates

## Available test data (local)

All in `min_data/`:
- `bsh_operationalmodel_data/static_file_{fine,coarse}/H0_file_{fine,coarse}.nc`
- `output_2d/c_file_{fine,coarse}_2020_surface.nc` — 24 timesteps (6h), 15 min
- `output_2d/c_file_{fine,coarse}_2020_bottom.nc` — same for bottom
- Curonian Spit: coarse grid ~20.5–21°E, 55–56°N

Parcels 3.0.6 available in `tmp_corner_fix/` pixi env and in the
Docker container `quay.io/willirath/parcels-container:2024.10.07-7af7fd0`.

## Sketch

See `plans/corner_rounding_sketch.png`: cell with SE corner rounded by
arc. The blue region is the effective water area, gray is the override
zone. The ψ-floor approach replaces the geometric arc with the
hyperbola ψ = ψ_min, which is a natural streamline shape.

## Open questions

1. **NUDGE_C calibration**: Must be large enough to prevent trapping
   but small enough to avoid oscillations. Step 0 will determine the
   working range.
2. **Interaction with Parcels slip conditions**: The B-grid
   `calculate_slip_conditions_2D` is not used for `cgrid_velocity`,
   so no conflict. But verify.
3. **Cells with two concave corners**: Narrow channels or peninsula
   tips. The nudges are independent and should compose naturally, but
   test explicitly.
4. **ψ_min sensitivity**: The value 0.015625 is a starting point.
   May need tuning per grid resolution (fine grid cells are ~900m,
   coarse ~5km — the physical exclusion zone scales differently).
