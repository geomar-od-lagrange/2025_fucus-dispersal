# Bottom friction for the `bottom` regime

Add the lab bed-friction model to the `bottom` drifters in 010: a *Fucus*
sliding on the bed moves along the local flow `U` at
`u = max(|U| − U_b, 0)` (Coulomb bed friction + square-law drag; rest for
`|U| ≤ U_b`), `U_b = 0.049 ± 0.005 m/s` from the flume fit. Prototype,
numbers and draft code in `tmp_bottom_friction/` (untracked), runs in
`<outputs>/tmp_bottom_friction/`.

## Start here (state 2026-11-13, late)

- **Wired** (branch `worktree-agent-addfa4b95e8442caa`, not merged: the main
  checkout serves the ε = 2 rerun's refills). Docs:
  [bottom_friction.md](../docs/bottom_friction.md); §Wiring below is done.
- **Decisions (authors, 2026-11-13):** implementation route approved:
  friction in the 010 kernel, on the interpolated UV at every RK4 stage
  (§Design; no field-level scaling, load hook or preprocessed fields); ε
  added after the friction. Own
  regime `bottom_ub49` beside `bottom` (which stays). One `U_b` = 0.049 m/s,
  no ± members. `velocity_scale` / `bottom_vs0p8` dropped entirely (not a
  yardstick either; bottom-flow uncertainty is discussed separately in the
  paper). Full 2016–2025 sweep, after the ε = 2 rerun.
- Prototype done (§Prototype results); test 5 done (§Test plan).
- Next: merge after the rerun's refills → `submit_010.sh bottom_ub49` for
  2016–2025 → `submit_downstream_chain.sh` on those ids (024/024b/024c,
  025–028 + 034 for `bottom_ub49`) → paper numbers (test 6) → this plan to
  `plans/done/`.

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
- **Parameters.** `U_b` is a fieldset constant (`add_constant`). 010
  parameter `break_off_speed_mm_s` (0 = off); output regime
  `<regime>_ub<mm/s>` = `bottom_ub49`; `regime` still picks the 2D fields
  and walls (`bottom`).
- **U_b = 0 is the current kernel bit for bit** (`f = 1 − 0/s = 1.0`
  exactly; `f = 0` only where `u = v = 0`), so one kernel serves all
  regimes and the identity is testable.
- **JIT.** `math.sqrt`, `math.cos`, `math.pi` and an `if` are plain C;
  compiles in the production container.
- **Cost.** ≤ 10 % 010 runtime (3-d smoke 70–73 s vs 66 s; 180 d 3798 s for both,
  21,800 drifters); no new data.

### Alternatives (rejected)

- Field scaling factor (`set_scaling_factor`): linear only (`u = αU`);
  cannot express a shift + dead zone.
- Hook on chunk loading (`computeTimeChunk` / `rescale_and_set_minmax`):
  no public hook in Parcels 3.0.6 → subclass or monkey-patch (AGENTS.md);
  plus the face-value problems below.
- Preprocessed bottom c_files with `f(|U_face|)·U_face`: +2.0 TB (the
  bottom c_files total 2.0 TB; quota ≤ 28 T of 29 T), and again for any
  other U_b. A C-grid face carries one component, so `|U|` at a face needs
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
- `U_b` = 0.049 ± 0.005 m/s (lab fit, χ² = 1.0 at 2 dof; one flume); the
  run uses the central value only (± 5 mm/s moves the numbers below by
  ~4 points).
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

| sample | n | median \|U\| | rest at U_b 44 / 49 / 54 mm/s | path kept at 44 / 49 / 54 |
|---|---|---|---|---|
| Eulerian, fine wet cells | 35.3 M | 0.103 m/s | 21 / 24 / 27 % | 71 / 69 / 66 % |
| Eulerian, coarse wet cells outside fine | 9.7 M | 0.074 m/s | 35 / 38 / 40 % | 70 / 67 / 65 % |
| Lagrangian, drifters 0–180 d | 10.3 M | 0.035 m/s | 60 / 65 / 69 % | 35 / 31 / 27 % |
| … 0–30 / 30–90 / 90–180 d | | 0.037 / 0.038 / 0.032 | 64 / 61 / 68 % (49) | 31 / 33 / 29 % (49) |

Per drifter, share of time at rest (U_b 49 mm/s), 10 / 50 / 90 %
quantiles: 47 / 63 / 88 %.

- The drifters sit in much slower bottom water than the grid average
  (released on the coast, collect in slow places): the median bottom speed
  they see is *below* U_b. Friction stops them ~2/3 of the time and keeps
  ~1/3 of their along-flow path length. The U_b ± 5 mm/s range moves
  "rest" and "path kept" by ± 4 points.
- Bias of the estimate: positions from the no-friction run; with friction
  drifters stay nearer their release (smoke: median 3-d distance 4.3 →
  0.7 km). Size against the sampling floor: §Prototype results.

## Prototype results

`tmp_bottom_friction/` — `kernel.py`, `b_run.py` (release = obs 0 of the
tmp_wall_velocity bottom run, 21,800 drifters, 2019-10-18, production
`wall_fields/walls_*_bottom.nc`, ε = 2 mm/s, reject on), `s_launch.sh`
(base partition, 1 cpu / 15G per job), `c_compare.py`.

- **Identity** (`c_compare_smoke.log`): friction kernel at U_b = 0 vs
  010's `AdvectionRK4Wall` verbatim, 3 d, 21,800 drifters: lon, lat, nrej
  identical bit for bit.
- **JIT = Python** (`c_compare_mode.log`): friction kernel at 49 mm/s,
  JIT vs Scipy, 109 drifters × 2 d: max separation 0.3 mm.
- **Smoke** (3 d, 14,253 water seeds): median distance from release
  4.33 km (no friction) → 0.72 km (U_b 49); hours with < 1 m displacement
  0 → 43 %; 6 % of drifters never move > 1 m in any hour; rejected land
  steps 307 → 12.
- **dt** (`c_compare_dt.log`; dt 2.5 vs 5 min, position separation after
  10 d): no friction median 7 m / p90 199 m / 5.1 % > 1 km; friction
  1 m / 15 m / 0.8 %. The kink adds no visible dt sensitivity; dt 5 min
  stays.
- **Cost (3 d)**: execute 70–73 s vs 66 s (1 cpu, JIT; job-to-job noise).
- **Effect** (`c_compare_effect.log`; 180-d runs `run_wall` / `run_ub49`,
  14,253 water seeds; 30 km EPSG:3035 squares; overlap = Σ min(p, q) of
  the 0–H d occupancy; floor = `bottom` half vs half, rng 20261103):

  | H (d) | 30 | 60 | 90 | 180 |
  |---|---|---|---|---|
  | median distance from release at H, `bottom` / `bottom_ub49` (km) | 24.8 / 8.4 | 42.6 / 16.3 | 51.1 / 24.1 | 68.3 / 32.9 |
  | area holding 90 % of particle-time 0–H, `bottom` / `bottom_ub49` (10³ km²) | 177 / 167 | 198 / 173 | 205 / 177 | 201 / 178 |
  | overlap `bottom_ub49` vs `bottom`, 0–H | 0.849 | 0.793 | 0.751 | 0.719 |
  | floor, `bottom` half vs half, 0–H | 0.925 | 0.936 | 0.940 | 0.946 |

  The area-90 metric moves little (−6 to −14 %): the spread of the ~500
  release cells along the coast sets it. Drift distance and overlap are
  the metrics that see the friction.
- **Cost (180 d)**: execute 3798 s for both kernels (1 cpu each).

## Wiring (done; as built in docs/bottom_friction.md)

As built, deviating from the draft below: U_b reaches the bottom job script
as a third positional arg pinned by `scripts/010_FucusDispersal_bottom_ub49_job.sh`
(not a `BREAK_OFF_SPEED_MM_S` env var), so `submit_010.sh bottom_ub49` and the
chain's regime discovery work by script name; 034 is `034_BottomFriction`,
pooled over a season of `release_years` on the releases both regimes hold;
024b's stem regex takes any regime word; the consumer matrix runs
`bottom_ub49` (025–028 unweighted + 034) by default, and the chain skips
024d/024efg without a `surface_stokes` job and runs only its own regimes'
consumers.

Draft:

- `notebooks/010_FucusDispersal.py`: as `tmp_bottom_friction/010_wiring.diff`
  — `break_off_speed_mm_s` parameter, `output_regime` suffix `_ub<mm/s>`,
  `add_constant("U_b", …)`, kernel `AdvectionRK4WallFriction` (all
  regimes; U_b = 0 is the current kernel bit for bit), `import math`;
  `velocity_scale` parameter and the "Velocity scale" cell removed. Sync
  `.md` / `.ipynb` (jupytext skill).
- `scripts/010_FucusDispersal_bottom_job.sh`, `scripts/submit_010.sh`:
  `VELOCITY_SCALE` → `BREAK_OFF_SPEED_MM_S` (env, default 0; 49 →
  `Trajectories/bottom_ub49/`), executed-notebook name.
- `notebooks/034_BottomSlowdown.*` repurposed (authors) as the `bottom` vs
  `bottom_ub49` comparison: `regime_var = "bottom_ub49"`, labels and
  markdown (friction, not a uniform slowdown); `release_doys_csv` default
  (now "the set the 0.8× variant was run for") → the doys / years
  `bottom_ub49` is run for (full 2016–2025 sweep: all, or the `year`
  parameter's full set), so its `year` / doy parameters change meaning;
  probably a rename (`034_BottomFriction`); AGENTS.md 034 bullet,
  `docs/visualisations.md` (034 row and §), its job/exec names.
- `bottom_vs0p8` follow-through: comments in `024_BuildHexAggregates.py`, `005_BuildWallFields.py`,
  `scripts/005_BuildWallFields_job.sh`; `docs/job_scripts.md`,
  `docs/wall_velocity.md` (ε vs `velocity_scale`, vs0p8 walls). (The ε = 2
  `bottom_vs0p8` jobs were cancelled; no vs0p8 data from that rerun.)
- `docs/bottom_friction.md` (new), this plan → `plans/done/`.
- Downstream: `bottom_ub49` is one more regime for 024/024b/024c/024d–g and
  the 025–036 consumers; regime discovery is by directory, check the
  notebooks that list regimes explicitly.

## Open

- Paper framing of the friction effect: deferred (authors) until the full
  2016–2025 `bottom_ub49` results exist.

## Test plan

1. U_b = 0 vs 010 kernel verbatim, same fieldset: bit for bit on lon,
   lat, nrej (done, 3 d, 21,800 drifters).
2. Smoke U_b = 0.049, 3 d (done).
3. dt 2.5 vs 5 min, 10 d, with and without friction (done).
4. Effect size, 180 d, one release: median distance from release; area
   holding 90 % of particle-time; overlap of 0–H d occupancy with `bottom`
   against the `bottom` half-vs-half floor (done, prototype numbers above;
   30 km squares, not the 024a hexes).
5. After wiring (done; `<outputs>/tmp_friction_integration/`,
   `compare_AB.log`): 010 via papermill, 2019-10-18 release (seed 32036115,
   87,200 drifters, 3 d): `break_off_speed_mm_s = 0` vs the pre-friction 010
   notebook, lon/lat/nrej/time/land_seed identical bit for bit (6.28 M
   values); at 49 vs the prototype `smoke_ub49` and the first 72 h of
   `run_ub49` (21,800 mapped trajectories), lon/lat/nrej/land_seed identical.
   End to end: `submit_010.sh bottom_ub49 2019` (RELEASE_NS=58,
   SIMULATION_DAYS=5, tmp OUTPUT_ROOT) → `Trajectories/bottom_ub49/2019/`;
   024/024b/024c (bottom and bottom_ub49), 027/028 (bottom_ub49) and 034
   (0–10 d bin) on that root all run clean (`logs/D_24918793.out`). Zarr
   size, 3 d: 53.0 MB friction vs 68.9 MB frictionless (0.77); prototype
   180-d runs mid-way: 0.72.
6. Paper numbers: the same two headline metrics (area holding 90 % of
   particle-time per horizon, occupancy overlap vs the half-vs-half floor)
   for `bottom` vs `bottom_ub49` from the 024 counts store, all years.
