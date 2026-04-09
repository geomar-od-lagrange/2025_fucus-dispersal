# Corner Trapping Theory: From Problem to ε-floor

Summary of the exploration in `tmp_corner_fix/corner_vortex.md`.

## C-grid stream function

For a cell with edge velocities U₀ (west), U₁ (east), V₀ (south),
V₁ (north) and continuity U₁ - U₀ + V₁ - V₀ = 0:

    u(x, y) = (1-x) U₀ + x U₁    (linear in x, constant in y)
    v(x, y) = (1-y) V₀ + y V₁    (linear in y, constant in x)

Stream function: ψ = U₀(1-x)y + U₁xy - V₀x + C

For a concave SE corner (U₁ = 0, V₀ = 0, V₁ = U₀):

    ψ = U₀ (1-x) y

Streamlines are hyperbolas (1-x)y = K, all converging to (1,0).

## Why particles get trapped

NOT because of RK4 error. The RK4 amplification factor for ψ is:

    R(-h) · R(h) = 1 + h⁶/72 + O(h⁸) > 1

So ψ INCREASES each step — RK4 is repulsive at the corner.

The actual mechanisms:
- **Float32 precision**: velocity rounds to zero near the corner
- **Time-varying fields**: velocity updates reshuffle particles
  across streamlines
- **Edge stagnation**: particles landing exactly on a V₀=0 edge
  have zero cross-edge velocity forever

Demonstrated with Parcels: particle 600 in the Curonian Spit test
has lat frozen at 55.3189037106 for all 23 observations — zero
displacement in latitude, oscillating ~4m in longitude with the
tide.

## Potential corrections explored

### Point source (atan2)

    ψ_corr = ε atan2(y, x-1)
    u_corr = ε(x-1)/r², v_corr = εy/r²

Correct direction (pushes away from both land edges). But
asymmetric: strong westward push near east edge, weak northward
push near south edge (v_corr ~ y/r² → 0 as y → 0).

### Point vortex (ln r)

    ψ_corr = ε ln r
    u_corr = εy/r², v_corr = ε(1-x)/r²

Wrong direction: pushes eastward near the east land edge.

### Higher multipoles (dipole, quadrupole)

Steeper concentration at corner, but angular structure creates new
problems. Dipole pushes into the south edge. Quadrupole creates
closed circulation. Only the source (ln z) gives purely radial
outward flow in 2D.

### Coordinate remapping

    x' = x - δ g(x,y),  y' = y + δ g(x,y)

Evaluate velocity at (x', y') instead of (x, y). Edge-preserving
with g = x(1-y). Creates a local minimum in ψ at the corner
(not negative — closed contours, not a true exclusion zone). Also
introduces divergence: positive near south edge (repulsion, good),
negative near east edge (attraction, bad).

## The ε-floor (horizon construction)

The simplest and most correct approach:

    ψ_new = (1-x)y - ε x(1-y)

Velocity: u = 1 - x(1-ε), v = ε + y(1-ε). This is the standard
C-grid interpolation with U₁ = ε, V₀ = ε instead of zero.

Properties:
- ψ matches original exactly on free edges (west, north)
- Horizon (ψ = 0): smooth curve from (0,0) to (1,1)
- Shadow zone (ψ < 0): between horizon and corner
- No stagnation point: u > 0 and v > 0 everywhere
- Divergence-free (continuity ε - 1 + 1 - ε = 0)
- One parameter ε
- Transit time error O(ε) for most streamlines
- Total west→north transport exactly preserved
- Parasitic land flux = ε (small)

On the south edge (y = 0): v = ε > 0. Particles immediately
lifted off the edge. This directly fixes the observed edge
stagnation.

## Shadow zone handling

Continuous nudge for particles past the horizon:

    deficit = max(0, -ψ)
    u += α · deficit · ∂ψ/∂x
    v += α · deficit · ∂ψ/∂y

Pushes along ∇ψ (perpendicular to streamlines, toward the
horizon). Zero correction and zero derivative at ψ = 0 when
using deficit². Not divergence-free in the shadow zone, but
that's virtual land — doesn't matter.

## Choosing ε

Lower bound from float32: ε > 3e-4 (displacement per step must
exceed position noise of ~0.2m).

Practical range: ε = 0.001 to 0.01.

| ε      | Shadow zone | Parasitic flux | Transit error |
|--------|-------------|----------------|---------------|
| 0.001  | 3% diagonal | 0.1%           | 0.1%          |
| 0.01   | 9% diagonal | 1%             | 1%            |

## Flow direction independence

The ε-floor correction is independent of flow direction. Land edges
are always zero regardless of which direction the flow enters/exits.
The correction always pushes away from the land boundary. Verified
analytically for reversed flow (U₀ < 0, V₁ < 0).

## Files

- `tmp_corner_fix/corner_vortex.md` — full jupytext notebook with
  all derivations and plots (20 figures)
- `tmp_corner_fix/corner_vortex.ipynb` — executed version
- `plans/corner_rounding_sketch.png` — cell geometry
- `plans/corner_streamlines.png` — streamline trapping
- `tmp_corner_fix/curonian_test.png` — Parcels 6h displacement map
- `tmp_corner_fix/curonian_test2.png` — 48h trajectory plot
