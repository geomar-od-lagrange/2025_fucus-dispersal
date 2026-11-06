# Numerical coastal trapping

Drifters in 010 pile up against the model's zero-flux coast faces and freeze
there. Findings below (2026-10-27/28, all from the 2019 sample, scratch in
`tmp_stall_audit/`, `tmp_stokes_faces/`, `tmp_pinned_time/`,
`tmp_land_entry/`, `tmp_pin_sites/`, `tmp_bottom_check/`; all untracked); the fix (a wall velocity applied in the
010 kernel) is planned at the end. Found during the figure review:
[figure_review_r2.md](figure_review_r2.md).

## Start here (state 2026-11-03)

- **Diagnosis done** (sections below, newest last). Pinning at zero-flux
  coast faces is the dominant fate with Stokes (52 % of particle-hours over
  180 d, 82 % of drifters pinned at day 180) and moderate without (surface
  21 %, bottom 7 % of particle-hours). Which coast stretch pins is physical
  (≥ 30 km); the exact pin cell is a grid artefact. 0.9 % of drifters also
  enter land cells (RK4 stage into a zero-velocity land cell; nest-seam
  mismatch). Bottom: no grid error, pins follow the bottom return flow.
- **Chosen direction (authors, 2026-11-03):** fix in the 010 kernel with a
  small wall velocity ε pointing away from land — no new 2D fields. Not
  started. Beaching stays with the rate model (`trap_wall` = hard shore
  substrate, not a C-grid wall); post-processing stalls as stranded is the
  fallback, not the plan.
- **Next steps:** (1) `tmp_` prototype: static wall field per nest
  (+ε on permanent zero faces) + RK4 kernel sampling `UV + UVwall`; one
  release, ε = 0 / 1 / 5 mm/s; metrics in §Fix. (2) Authors pick ε.
  (3) Wire into 010 + docs, review agent. (4) Full rerun (§Fix: cost,
  quota) and downstream chain; then the paper figures and numbers refresh
  (styling from the R2 round is final).
- **Held until after the fix:** dilution colour top (1e-1 vs 1e-2), the ¼-hex
  dilution area floor, the zero-water hex remap (G6) — see
  [figure_review_r2.md](figure_review_r2.md).

## Finding: coastal dilution hot spots and stalled drifters (2026-10-27)

Diagnostics in `tmp_figreview_r2/dilution_coast/` (untracked; 2019 trajectory
sample = every third zarr, 1.42 M drifters per regime).

- Every hex above 1e-3 per km² is a sliver of model water (0.03–0.45 km² of a
  93.5 km² hex); its weight is ordinary, the 1/area factor is ×10^2.5–10^3.4.
  The top 1 % of hexes (skerry coasts: Swedish east coast, Archipelago Sea,
  Gulfs of Finland and Riga) hold 9.5 % (30 d) – 38.5 % (180 d) of the
  surface weight.
- Underneath is real retention: 88.6 % of surface drifters never move ≥ 10 m
  again by 180 d (median stop day 76; bottom 25.5 %). Stalled drifters keep
  survival weight (mean exp(−A) 0.44, same as moving ones; `trap_wall = 0`
  gives walled shores zero stranding rate). The first diagnostic called this
  the C-grid corner trap; that is unverified — its "on F-points" test ran
  only inside the top hexes, as a half-cell offset from the H0 coarse
  coordinates. Superseded by the stall audit below.
- Land seeds are excluded correctly (both regimes count water-seeded
  drifters: surface N = 37,296,956, bottom 37,277,160).
- Options, by soundness of the result: (1) floor the dilution area at ¼ hex
  (parquet consumers only, minutes; display fix only); (2) treat
  stopped-for-good drifters as stranded (024d moved-flag + 024e/f/g
  reducers, ~20 core-hours, plus the 024 counts rebuild for bottom; a
  methods decision); (3) fix the corner trap in 010 and rerun everything.
## Stall audit (2026-10-27, `tmp_stall_audit/`, verified with the 010 fieldset)

- Stopped for good (no hourly step ≥ 10 m for the last ≥ 7 d before
  day 180; water-seeded drifters; 2019 sample, 4 releases × 87,200
  drifters): surface_stokes 83.4 % (median stop day 75; 44.8 % of
  particle-hours after the stop), surface 37.8 %, bottom 17.9 %.
- Mechanism: C-grid interpolation drives the wall-normal velocity linearly
  to zero at a zero-flux coast face, so onshore flow pulls drifters
  exponentially close to the wall (~1e-12 cell). They stall for good where
  the along-wall component is zero too: concave coast corners (59.8 % of
  stopped surface_stokes drifters) or cells open on one side (32.9 %).
  The 2D surface (+ Stokes) field is convergent, so corners are genuinely
  attracting — not the divergence-free saddle of corner_theory.md, not a
  Parcels bug, not the edge mask (a no-op on coarse), not registration
  (99.9–100 % match), not float32 (positions are float64). Production
  reruns bit-identically (Parcels 3.0.6 container).
- Minor: nest-boundary coast mismatch puts ~1 % into fine land cells;
  int16 velocity storage zeroes |u| < 0.5 mm/s and 003 then zeroes Stokes
  there (0.5 % coarse / 2 % fine of open-water face steps).
- Reading: beaching under onshore convergence that never reverses. Fix
  candidates: post-process stalled drifters as stranded (no rerun of 010),
  or a beaching kernel in 010. Superseded by §Start here.

## Stokes face mask audit (2026-10-27, `tmp_stokes_faces/`, verified)

- 003 (`notebooks/003_prepare_2d_fields.py` ~L305–309) zeroes Stokes where
  `u_surf` is NaN or exactly 0, sums, then re-applies the land edge mask to
  the total. On land faces the Stokes mask is redundant (the edge mask
  already zeroes them); its only independent effect is dropping Stokes on
  open-water faces where the int16 current rounds to 0 (0.38 % coarse /
  2.45 % fine of face steps, 0.3 % / 1.6 % of Stokes |flux|). Re-running
  003's functions on 24 stems of 2019 reproduces production to 2.2e-16;
  no wall face is nonzero.
- A zero wall face stops flux *through* the wall, not flux *toward* it:
  the wall-normal velocity is linear between the open face (onshore Stokes
  ~60 % of the time) and 0 at the wall (u/d constant from d = 0.9 to 1e-12,
  sampled with Parcels). Drifters approach exponentially (0.2–3.3 e-folds/d
  at the sites, +104 to +1490 net over 450 d) and are then frozen by
  float64 rounding: at ~1e-12 cell a 5-min step is below half an ulp unless
  the open-face velocity exceeds ~0.2 m/s.
- Trapping condition per cell: no face carries outflow. Share of time it
  holds (current+Stokes / current only): 2-wall cells 14 % / 9 % coarse,
  3-wall cells 45 % / 43 % coarse. Stokes roughly doubles it at corners.
- Land entry: 0.93 % of water-seeded surface_stokes drifters (2019) enter a
  land cell, 1,175 at nest switches, 957 within one grid by diagonal steps
  (suspected RK4 sub-step leak at convex corners, not traced).
- Options unchanged: post-process stalls as stranded (recommended, no 010
  rerun); random walk with reflecting walls (rerun); zero onshore Stokes in
  the whole wall cell (moves the problem one cell, unphysical). Cheap 003
  cleanup: mask Stokes on the land mask, not `u_surf == 0`; fix the
  "no-slip walls" wording in docs/stokes_drift.md.

## Pinned time (2026-10-27, `tmp_pinned_time/`, headline numbers re-derived)

Hourly classification of every water-seeded drifter (2019, 4 releases,
~226 k per regime, 0–180 d): pinned = hourly step < 1 m and within 1e-3
cell of a permanent zero-flux face of a wet cell (spot checks: 1e-13–1e-12).

| | bottom | surface | surface_stokes |
|---|---|---|---|
| pinned share of particle-hours, 0–30 / 0–90 / 0–180 d | 0.4 / 3.0 / 7.1 % | 2.6 / 11.1 / 20.7 % | 12.9 / 33.3 / 52.0 % |
| drifters pinned at day 30 / 90 / 180 | 1.0 / 7.5 / 16.4 % | 5.8 / 23.2 / 39.4 % | 22.0 / 61.3 / 82.5 % |
| final pin length, median (pinned at end) | 47 d | 59 d | 97 d |
| final pin ≤ 7 d / ≥ 30 d | 17 / 61 % | 18 / 66 % | 6 / 82 % |

- Not a late artefact: with Stokes half of all particle-hours are pinned and
  the median final pin starts at day 82. Earlier pin episodes mostly release
  within hours (median 4 h with Stokes); ~6 % last > 7 d.
- 10 m step threshold changes the shares by ≤ 4 points. Release spread
  (surface_stokes, 0–180 d): 38 % (April) to 65 % (October).
- Stalled in land cells ≤ 0.7 % of particle-hours; stalled in open water
  ≤ 3 %.
- Spatial signal (`tmp_pinned_time/e_spatial.py`, ~1 km boxes of the
  day-180 position): pins collect at few persistent sites. surface_stokes:
  186,577 pinned drifters at 1,174 sites; 64 sites hold 50 %, 345 hold
  90 %; the top 8 sites (2–2.6 % each) are hit by all 4 releases (e.g.
  22.17E 60.45N, 24.58E 58.35N, 29.83E 60.20N, 21.42E 62.30N). surface and
  bottom: 19 sites hold 50 %. Median release→pin distance 125 km
  (surface_stokes), 141 km (surface), 28 km (bottom): drifters travel for
  weeks before a downstream attractor catches them; pinning is absorbing,
  so the pinned fraction can only grow with age.

## Within-grid land entries (2026-10-27, `tmp_land_entry/`, verified)

- 957 of 228,822 water-seeded surface_stokes drifters (0.42 %, 2019) enter a
  land cell without a nest switch; all stay bit-frozen there. Plus 1,175 at
  nest switches (0.51 %). Particle-hours in land within 180 d: 0.29 % /
  0.42 %.
- Mechanism (all 957 entry steps rerun bit-identically in the 3.0.6
  container, every RK4 stage traced): AdvectionRK4 samples stage positions
  wherever they land, without checking that they stay in the cell or in
  water. At a convex corner of a one-cell island (144) or spit tip (813) a
  stage crosses into the land cell L, where all four faces are zero and the
  C-grid returns u = v = 0 exactly. That truncates the RK4 average to a
  point already inside L (91 %: patterns PLPL, PPPL). The other 9 % pick up
  a neighbouring water cell's wall-normal velocity and cross the start
  cell's closed face. Not RK4 overshoot: k·dt ≤ 0.26 at a 5-min step, and the
  one-cell RK4 wall distance has no real root.
- Generic: a synthetic 1-km C-grid island in uniform 0.3 m/s flow at 225°
  traps 20/4/0 of 4000 drifters at dt = 300/60/12 s, same stage pattern.
- Stores: counted as real afloat drifters (024, 024f); 40 % sit at
  flat_fraction 0 and never strand under trap weights 1/0.
- Fix options: post-process as stranded (free with the stall rule); a
  land-guard kernel with sub-stepping in 010 (rerun); smaller dt (5×
  cost, only reduces); AdvectionAnalytical (Scipy-only in 3.0.6, no
  NestedField).

## Pin sites: grid or physics? (2026-10-28, `tmp_pin_sites/`, key numbers checked in the logs)

- Verdict: mixed. Which coast stretch pins is physical; the exact pin cell is
  a grid artefact.
- Physical: 94 % of surface_stokes pins sit on walls with land to the E/N
  (52 % of walls), where the 2019 mean Stokes is onshore (+1.9 / +1.1 cm/s
  coarse). Bottom pins sit on W/S walls. Top sites are bay heads and
  downwind archipelago coasts.
- Grid: pins form only in corner / dead-end cells (99.9 %; straight walls
  0.1 %); staircase corners of diagonal coasts hold 24 %, one-cell dead
  ends 26 %. 30 % of top-90 % pins sit on a wall standing on ≥25 % real
  water (GSHHS full).
- Resolution test (the BSH coarse grid is land under the fine nest, so a
  coarse-only run is impossible): the fine surface_stokes field coarsened
  6×6 conserving face transport, two releases × ~13.8 k drifters, 180 d.
  Same drifter's pin site on the two grids: median 13–31 km apart.
  Distribution overlap cross-grid vs same-grid: 2 km 0.00 vs 0.52–0.66;
  6 km 0.15 vs 0.63–0.66; 15 km 0.32 vs 0.68; 30 km 0.63–0.65 vs 0.72–0.75.
- Consequence: per-cell / per-hex hot spots (incl. dilution peaks) are not
  robust; coast segments ≥ 30 km or subbasins are.

## Fix: wall velocity in the 010 kernel (proposed, not started)

- Idea: replace the zero wall face by a small velocity ε pointing away from
  land. In the wall cell `u = (1 − xsi)·ε + xsi·U_E`, so drifters converge on
  `xsi* = ε / (ε + |U_E|)` (≈ 0.05 cell at ε = 1 mm/s, 2 cm/s onshore)
  instead of 1e-12, stay resolvable in float64, and leave once the onshore
  forcing relaxes. A stage that lands in a land cell picks up outward
  velocity instead of exactly 0 (expected to remove most land entries;
  to be tested).
- Implementation without new 2D fields: a static wall field per nest and
  regime (+ε on permanent zero faces, 0 elsewhere; `cgrid_velocity`,
  NestedField) and an RK4 kernel sampling `UV + UVwall` at every stage. 003
  and the 6 TB of fields stay as they are; only 010 and downstream rerun.
  Bottom uses the same coastline walls (checked, see Open questions;
  authors agreed 2026-11-03).
- Beaching stays with the rate model (docs/beaching.md): in band,
  `trap = trap_wall + (trap_flat − trap_wall)·ff/100` on the shore
  substrate. `trap_wall` there means a hard *shore type* (rock, cliff,
  armoured bank; `flat_fraction` 0), not a C-grid wall face — every model
  coast is a C-grid wall. With pinning gone, drifters off rocky shores stay
  afloat and move on, as intended.
- Test first (`tmp_` prototype, one release, ε = 0 / 1 / 5 mm/s): pinned
  fraction and duration, coastal residence, land entries, open-water
  dispersal bias, robustness of the ≥ 30 km beaching pattern.
- Cost of the full rerun: 010 ~1–1.5 d wall clock (24 regime-years,
  ~15–20 k core-hours), downstream chain ~3 h. Quota (`workquota`
  2026-11-03): 25.1 T of 29 T soft; a full trajectory set is ~3–3.4 T
  (≈100 GB/yr surface_stokes, ≈200 GB/yr bottom).

## Open questions

- Bottom pins favour walls with land to the W (coarse 52 %) — checked
  2026-11-03 (`tmp_bottom_check/`), not a grid error: (a) on faces where
  only layer 1 is fluid, bottom and surface 2D files are bit-identical (652 /
  643 coarse, ~7.8 k fine faces per component), and an independent
  re-derivation of the deepest-fluid value from raw BSH matches the bottom
  file on 100 % of nonzero faces (every one-cell shift: 2–5 %); (b) the
  permanent zero faces of the bottom files are exactly the coastline walls
  (3,574 coarse, same as surface — blocked deep faces take a shallower
  layer); (c) 73 snapshots of 2019, mean / median flow toward the wall,
  onshore share: surface E +1.9 / +1.3 cm/s, 0.70; W −1.4 / −1.1, 0.31.
  Bottom W +0.12 / +0.23, 0.61; E −0.36 / −0.26, 0.36. The bottom flow runs
  opposite to the surface drift at the coast (return flow), and pins follow
  it. The bottom wall field for the kernel fix = the coastline walls.
- 003 drops Stokes on open faces where the int16 current rounds to 0
  (0.3 % coarse / 1.6 % fine of Stokes |flux|); the kernel fix leaves this
  in place unless 003 reruns.
