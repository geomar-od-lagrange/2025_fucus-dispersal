# Corner Rounding: ε-floor on Land Edges

## Problem

On a NEMO-style C-grid, velocity at land-adjacent edges is zero.
At concave corners (two land edges meeting), both U and V go to
zero, creating a stagnation point. Particles get trapped there
under finite-timestep integration — not by RK4 error (which is
actually repulsive, see `corner_theory.md`) but by float32
precision limits and time-varying velocity fields.

**Prerequisite**: The grid registration bug (`grid_registration.md`)
must be resolved first. The corner fix operates on the interpolation
level and requires correct cell assignment.

## Solution: horizon construction

Replace zero land-edge velocities with a small value ε:

    U₁ = ε,  V₀ = ε

(for a SE corner; analogous for other corners). This is equivalent
to modifying the stream function:

    ψ_new = (1-x) y - ε x (1-y)

which preserves ψ exactly on the free edges (west and north),
creates a smooth horizon (ψ = 0 curve) from (0,0) to (1,1) that
cuts off the corner, and makes ψ < 0 in the shadow zone beyond.
No stagnation point exists. Divergence-free. One parameter.

## Parameters

- **ε**: Land-edge velocity floor. For BSH fine grid (900m cells,
  2 m/s max, 5-min timestep), ε ≥ 0.001 is sufficient (gives
  0.002 m/s at the corner, well above float32 noise). Shadow zone
  is ~3% of cell diagonal. Parasitic land flux is 0.1% of main
  transport. Transit time error is O(ε) for well-behaved
  streamlines.

## Shadow zone repulsion

For numerical overshoots past the horizon (ψ < 0), add a
continuous nudge back toward ψ = 0:

    deficit = max(0, -ψ)
    U += α · deficit · ∂ψ/∂x
    V += α · deficit · ∂ψ/∂y

No if-branch. Uses deficit² for C¹ smoothness at the horizon.

## Implementation

Patch `spatial_interpolation_UV_c_grid` in `parcels.h`. For each
corner type, detect from zero edge values and apply the floor:

    if (is_zero(U1) && is_zero(V0)) {  /* SE corner */
        U1 = EPSILON;
        V0 = EPSILON;
    }

Then compute the standard C-grid interpolation with these modified
values. The shadow zone repulsion adds ~5 more lines.

Runtime swap the patched header via monkey-patching
`get_package_dir`. Works in read-only Singularity containers.

## Implementation steps

1. Resolve grid registration (`grid_registration.md`)
2. Prototype in `notebooks/002_CornerRounding.ipynb`:
   scipy particle integrator on real BSH data near the Polish coast,
   with and without ε-floor
3. Patch `parcels.h`, test in `tmp_corner_fix/` pixi env
4. Validate: 500 particles near coastline, verify no trapping and
   correct coastal trajectories
5. Production: add header swap to `010_FucusDispersal.ipynb`,
   re-run

## Alternatives explored and rejected

- **Circular arc**: Can't be implemented within C-grid bilinear
  interpolation (sub-cell effect)
- **Point vortex** (ln r): Tangential flow, pushes toward east
  land edge
- **Point source** (atan2): Correct direction but asymmetric,
  weak near south edge
- **Dipole/quadrupole**: Wrong angular structure, create new traps
- **Coordinate remapping**: Breaks divergence-free, creates
  attraction on one edge
- **ψ-floor nudge** (∇ψ proportional to deficit): Not
  divergence-free, but viable as shadow zone fallback

See `corner_theory.md` for the full exploration.

## Files

- `plans/corner_rounding.md` (this file)
- `plans/corner_theory.md` — derivations and exploration summary
- `plans/corner_rounding_sketch.png` — cell geometry
- `plans/corner_streamlines.png` — streamline trapping
- `tmp_corner_fix/corner_vortex.md` + `.ipynb` — full theory
  notebook with all plots
