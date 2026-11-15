# Bed friction for bottom drifters (010, `bottom_ub49`)

A *Fucus* sliding on the bed moves along the local bottom current `U` at
`max(|U| − U_b, 0)`: Coulomb bed friction plus square-law drag, at rest for
`|U| ≤ U_b`. The break-off speed `U_b = 0.049 ± 0.005 m/s` is the lab flume
fit; the study runs the central value only, as its own regime `bottom_ub49`
beside the frictionless `bottom`. Design, prototype numbers and the authors'
decisions: [bottom_friction.md](../plans/done/bottom_friction.md).

| piece | where |
|---|---|
| parameter | `notebooks/010_FucusDispersal.py`: `break_off_speed_mm_s` (0 = off), fieldset constant `U_b` (m/s) |
| kernel | `AdvectionRK4WallFriction` in 010, every regime |
| output regime | `<regime>_ub<U_b>` with `.` spelled `p`: 49 → `Trajectories/bottom_ub49/`; forcing and walls still follow `regime` (`bottom`) |
| job scripts | `scripts/010_FucusDispersal_bottom_ub49_job.sh` pins `U_b = 49` and runs `scripts/010_FucusDispersal_bottom_job.sh` (third arg); `scripts/submit_010.sh bottom_ub49 <year> [n]` |
| comparison | `notebooks/034_BottomFriction.py`: `bottom` vs `bottom_ub49` from the 024 counts store ([visualisations.md](visualisations.md) §034) |

## Kernel

At each of the four RK4 stages:

```
(u, v) = UV[stage position]                        # deg/s on the spherical mesh
s      = hypot(u·111120·cos(lat_stage), v·111120)   # |U| in m/s (Parcels' factors)
f      = 1 − U_b/s  if s > U_b  else 0              # direction-preserving
(u, v) = f·(u, v) + UVwall[stage position]          # wall velocity ε after friction
```

then the RK4 sum, the land-step reject and the land-seed skip as in
[wall_velocity.md](wall_velocity.md).

- **Interpolated velocity, every stage.** The model is pointwise: an object
  at `x` feels the C-grid interpolated `U(x, t)`, so the step is RK4 of
  `dx/dt = f(U)·U`. `f(interp(U)) ≠ interp(f(U))`.
- **Units.** `fieldset.UV` returns deg/s; `|U|` is rebuilt in m/s from
  Parcels' converter factors (1852·60 = 111 120 m per degree, × cos lat for
  U) at the stage latitude, one sample per stage. `f` is a scalar and scales
  the deg/s components directly.
- **ε after friction.** The wall velocity is a numerical device, not current
  the drifter feels. Rest is exact rest away from walls; in a wall cell a
  resting drifter still creeps at ≤ ε, at most one cell. Under sliding flow
  the wall fix is still needed: `f` preserves direction, so the wall-normal
  velocity still vanishes at a zero-flux face (rest point
  `ξ* = ε / (ε + f·|U_open|)`).
- **`U_b = 0` is the frictionless kernel bit for bit** (`f = 1 − 0/s = 1.0`
  exactly; `f = 0` only where `u = v = 0`), so `bottom` and `surface_stokes`
  run the same kernel.
- **Cost.** +5–10 % 010 runtime (prototype, JIT); no new input data.

## Job script per regime

`bottom_ub49` has its own job script, a thin wrapper that pins `U_b` and runs
the bottom one: `submit_010.sh` maps the regime argument to
`scripts/010_FucusDispersal_<regime>_job.sh` and names the job after it, and
`submit_downstream_chain.sh` reads the regime off the script name of each 010
job, so one script per output regime keeps all three right with no
regime-name parsing; the run logic stays in one file.

## Paired seeds

Each `bottom_ub49` release reuses the RNG seed of the `bottom` release with
the same start time, read off its zarr name
`Trajectories/bottom/<year>/Fucus_BSH_<YYYYMMDDTHHMMSS>_bottom_dt60min_seed<N>.zarr`
under the same `output_root` (`paired_seed` in the bottom job script). Both
regimes then start from identical positions, so they differ by the friction
alone (at U_b = 0 the run is the `bottom` zarr bit for bit) and per-release
differences need no sampling floor for the release jitter. A release with no
or several matching `bottom` zarrs is not run: the job prints
`FAIL <stem> bottom_ub49: …` and exits non-zero once the other releases are
done; it never draws a fresh seed. Hence `bottom_ub49` is submitted only once
the `bottom` zarrs of its years are on disk (a refilled `bottom` release is
paired with its refill's seed). Considered: fresh seeds as for every other
regime, rejected: friction would then be compared across different release
positions in every cell.

Downstream, `bottom_ub49` is one more regime for every stage: 024/024b/024c,
the beaching branch (024d–g, 029–031, same criterion as every regime,
[beaching.md](beaching.md)) and the survival-weighted 025–028 consumers.
`scripts/025-031_ConsumersAll_job.sh` runs it by default and adds 034, which
compares `bottom` and `bottom_ub49` in area and in stranded %.

## Considered, rejected

- Uniform velocity factor (`set_scaling_factor`): linear only (`u = αU`);
  cannot express a threshold with a dead zone.
- Friction on preprocessed bottom fields (`f(|U_face|)·U_face`): +2 TB per
  `U_b`, `|U|` at a C-grid face needs a stencil of its own, and zeroing
  sub-threshold faces creates zero-flux faces in open water — the coastal
  pinning mechanism.
- Hook on chunk loading: no public hook in Parcels 3.0.6.
- Friction on the RK4-averaged velocity: not RK4 of the friction ODE.

## Caveats

`U` is the deepest model layer's mean velocity (no log-layer correction);
one flume fit, central `U_b` only; no static/sliding hysteresis; the kink at
`|U| = U_b` drops RK4 to low order on steps that cross it (dt 2.5 vs 5 min
showed no added sensitivity, plan §Prototype results).
