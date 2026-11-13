# Bottom friction for the `bottom` regime

Add the lab bed-friction model to the `bottom` drifters in 010: a *Fucus*
sliding on the bed moves along the local flow `U` at
`u = max(|U| − U_b, 0)` (Coulomb bed friction + square-law drag; rest for
`|U| ≤ U_b`), `U_b = 0.049 ± 0.005 m/s` from the flume fit. Prototype,
numbers and draft code in `tmp_bottom_friction/` (untracked), runs in
`<outputs>/tmp_bottom_friction/`.

## Start here (state 2026-11-13)

- **Not wired into 010.** The ε = 2 mm/s production rerun is running on
  the current 010; nothing under `notebooks/` or `scripts/` changes until
  it and its refills are done.
- Design chosen (below): friction in the 010 kernel, per RK4 stage, on the
  interpolated velocity; ε added after the friction (authors, fixed).
- Prototype done: JIT kernel, runner, U_b = 0 identity (bit for bit), 3-d
  smoke, dt check, 180-d effect runs (results in §Prototype results).
- Next: authors answer §Open questions → wire into 010 (draft:
  `tmp_bottom_friction/010_wiring.diff`) + `docs/bottom_friction.md` →
  review agent → 010 sweep for the friction regime(s) → 024 chain → 034-type
  comparison.

## Design: friction in the 010 kernel

`tmp_bottom_friction/kernel.py: AdvectionRK4WallFriction` = 010's
`AdvectionRK4Wall` with, at each of the four RK4 stages,

```
(u, v) = UV[stage position]                       # deg/s on the spherical mesh
s      = hypot(u·111120·cos(lat_stage), v·111120)  # |U| in m/s (Parcels' own factors)
f      = 1 − U_b/s  if s > U_b  else 0             # direction-preserving
(u, v) = f·(u, v) + UVwall[stage position]         # ε after friction
```

then the RK4 sum, `Land` reject and land-seed skip unchanged.

- **Why at the interpolated velocity.** The model is pointwise: an object
  at `x` feels the flow at `x`, i.e. the C-grid interpolated `U(x, t)`. The
  ODE is `dx/dt = f(U(x,t))·U(x,t)`; applying `f` at each stage is RK4 of
  that ODE. `f(interp(U)) ≠ interp(f(U))`: transforming face values moves
  drifters where the local flow is sub-threshold and stops them where it
  is not.
- **Units.** `fieldset.UV` returns deg/s; `|U|` is rebuilt in m/s with
  Parcels' converter factors (`1000·1.852·60`, × cos lat for U) at the
  stage latitude. One sample per stage, no `applyConversion=False` resample.
- **Parameters.** `U_b` is a fieldset constant (`add_constant`), so the
  U_b members (0.044 / 0.049 / 0.054) cost CPU only. 010 parameter
  `break_off_speed_mm_s` (0 = off); output regime `<regime>_ub<mm/s>`
  (`bottom_ub49`), composable with `_vs<scale>`.
- **velocity_scale** stays a field scaling factor, applied before the
  friction (it scales the flow the drifter feels).
- **U_b = 0 is the current kernel bit for bit** (`f = 1 − 0/s = 1.0`
  exactly; `f = 0` only where `u = v = 0`), so one kernel serves all
  regimes and the identity is testable.
- **JIT.** `math.sqrt`, `math.cos`, `math.pi` and an `if` are plain C;
  compiles in the production container.
- **Cost.** +5–10 % 010 runtime (3-d smoke: 70–73 s vs 66 s execute,
  21,800 drifters); no new data.

### Alternatives (rejected)

- Field scaling factor (`set_scaling_factor`): linear only; expresses
  `bottom_vs0p8`, cannot express a shift + dead zone.
- Hook on chunk loading (`computeTimeChunk` / `rescale_and_set_minmax`):
  no public hook in Parcels 3.0.6 → subclass or monkey-patch (AGENTS.md);
  plus the face-value problems below.
- Preprocessed bottom c_files with `f(|U_face|)·U_face`: +2.0 TB per U_b
  member (bottom c_files total 2.0 TB; three members 6 TB against ≤ 28 T
  of 29 T). A C-grid face carries one component, so `|U|` at a face needs
  a 4-face stencil for the tangential part (a modelling choice of its
  own). The dead zone zeroes every sub-threshold face, creating zero-flux
  faces in open water — the exact mechanism that pinned drifters at the
  coast (bilinear interpolation drives the normal velocity linearly to 0
  toward such a face).
- Friction on the RK4-averaged velocity instead of per stage: not RK4 of
  the friction ODE; no cheaper.

## Wall velocity with friction (authors' decision: ε after friction)

`u = f(UV)·UV + UVwall`. ε is a numerical device either way. Consequences:

- **Rest is exact rest away from walls** (`f = 0` → `u = 0`). In a wall
  cell a resting drifter still moves at the wall velocity: ≤ ε = 2 mm/s
  away from the wall face, decaying linearly to 0 at the opposite open face,
  so it creeps at most one cell (e-folding Δx/ε ≈ 5 d fine 0.9 km, ≈ 29 d
  coarse 5 km). Where two walls face each other the creep ends mid-cell.
- **The wall fix is still needed.** `f` preserves direction, so with
  sliding flow (`|U| > U_b`, e.g. strong tangential flow along the coast)
  the wall-normal drifter velocity still vanishes linearly at a zero-flux
  face — the pinning mechanism. Off-wall rest point under onshore sliding
  flow: `ξ* = ε / (ε + f·|U_open|)`, further off the wall than without
  friction (f < 1).
- Land-step reject unchanged.

## Physics caveats (accepted; for the paper discussion)

- `U` in the lab model is the flow felt by the object near the bed; the
  forcing is the deepest model layer's mean velocity (layer thickness
  varies; no log-layer correction).
- `U_b` ± 0.005 m/s (lab fit, χ² = 1.0 at 2 dof); one flume, one specimen
  class.
- No hysteresis (one μ for static and sliding contact).
- Kink at `|U| = U_b` (continuous `u`, discontinuous derivative): RK4
  drops to low order on steps that cross it; checked with dt 2.5 vs 5 min
  (§Prototype results).

## How often is the bottom flow below U_b

`tmp_bottom_friction/a_speed_stats.{py,log}`. Eulerian: 2019, the 15th of
each month, 24 steps, wet T-cells (face-averaged speed). Lagrangian: the
ε = 2 mm/s bottom drifters of `<outputs>/tmp_wall_velocity/run_bottom_eps2mm_v3rej.zarr`
(2019-10-18 release, every 4th drifter, water seeds), C-grid interpolated
`|U|` every 6 h over 180 d. "rest" = share of samples with `|U| ≤ U_b`;
"path kept" = Σ max(|U| − U_b, 0) / Σ |U|.

RESULTS_Q4

## Prototype results

`tmp_bottom_friction/` — `kernel.py`, `b_run.py` (release = obs 0 of the
tmp_wall_velocity bottom run, 21,800 drifters, 2019-10-18, production
`wall_fields/walls_*_bottom.nc`, ε = 2 mm/s, reject on), `s_launch.sh`
(base partition, 1 cpu / 15G per job), `c_compare.py`.

RESULTS_PROTO

## Open questions (authors)

1. Does friction replace `bottom` (the paper's bottom regime becomes
   `bottom_ub49`) or run beside it as a sensitivity member?
2. U_b members: only 0.049, or 0.044 / 0.054 too (each a full bottom
   sweep, ~same CPU and disk as `bottom`)?
3. Keep `bottom_vs0p8` as the second bottom-uncertainty member, and does
   friction combine with `velocity_scale` (scale, then friction) in any
   member?
4. Full 2016–2025 sweep or one year (2019, as `bottom_vs0p8`)?

## Test plan

1. U_b = 0 vs 010 kernel verbatim, same fieldset: bit for bit on lon,
   lat, nrej (done, 3 d, 21,800 drifters).
2. Smoke U_b = 0.049, 3 d (done).
3. dt 2.5 vs 5 min, 10 d, with and without friction (done).
4. Effect size, 180 d, one release: median distance from release; area
   holding 90 % of particle-time; overlap of 0–H d occupancy vs `bottom`,
   against `bottom_vs0p8` vs `bottom` and the `bottom` half-vs-half floor
   (done, prototype numbers above; 30 km squares, not the 024a hexes).
5. After wiring: 010 via papermill at `break_off_speed_mm_s = 0` on one
   release = the current production zarr bit for bit; at 49 = the
   prototype zarr bit for bit.
6. Paper numbers: 034-type comparison (`bottom`, `bottom_ub49`,
   `bottom_vs0p8`) from the 024 counts store.
