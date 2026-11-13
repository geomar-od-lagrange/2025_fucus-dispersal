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
  with the R3 styling (`wr/figures-r3`, `plans/figure_review_r3.md`; merge it
  before the post-rerun figure pass — PNGs it renders now, and the R3 push
  to `paper/`, show pre-fix numbers). R3 adds a 029 stranded-curve export
  read by 035: run 026 and 029 before 035.
- **Decisions (authors, 2026-11-13):** ε = 2 mm/s, v3 walls (all steps of
  the forcing period + thin-wall rule) + reject-step kernel. No ε
  sensitivity in the paper. Go for step 3 and launch the full rerun this
  weekend (NESH base has ~2,400 idle cores spread over mix nodes; use many
  small-task jobs). Disk: keep ≥ 1 T headroom under the 29 T soft quota
  (i.e. ≤ 28 T used); delete `tmp_*` outputs, and delete old ε = 0
  trajectories before the rerun only where needed to stay under that,
  otherwise after the new ones exist.
- **Step 3 (wiring into the pipeline), in progress 2026-11-13:** design in
  §Fix → "Step 3: production integration". Regimes: surface_stokes, bottom,
  bottom_vs0p8 (plain `surface` is retired).
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
- Prototype design (step 1, `tmp_wall_velocity/`, in progress 2026-11-03):
  - Wall faces per nest and regime = faces exactly 0 at every step of the
    12 monthly 2019 files (`c_file_<nest>_2019MM1500_*_<regime>.nc`, as
    `tmp_pin_sites/a_grid.py`). Wet T-cell = at least one non-wall face.
    A wall face between one wet and one non-wet cell gets ±ε pointing into
    the wet cell; wall faces with both or neither side wet stay 0 (counted).
  - Stored as a time-constant NetCDF per nest with the c_file layout (same
    lon/lat order, `uvel`/`vvel`), loaded like 010's fieldset, a second
    NestedField pair `UVwall` (`cgrid_velocity`, same grids → same nest
    choice as `UV`). Kernel: `AdvectionRK4` with every stage sampling
    `UV + UVwall`.
  - Runs (production container, 010 settings, 180 d): release 2019-10-18
    (highest surface_stokes pin share), every 4th particle; regimes
    surface_stokes and bottom (ε = 1 mm/s is the size of the mean bottom
    flow at walls) × ε = 0 / 1 / 5 mm/s. ε = 0 must match the production
    zarr bit for bit on the subset.
  - Status 2026-11-03 09:00: wall fields, kernel, sign check (40/40) and
    ε = 0 bit-identity vs production (21,800 traj × 5 d, max |Δ| = 0) done.
    Wall = 0-or-NaN at every step (Parcels fills NaN with 0); fine bottom
    differs from fine surface_stokes by ~156 faces, so walls are per regime.
    Overhead +19 % runtime. ε = 5 mm/s smoke: formerly pinned drifters sit
    0.12–0.38 cell (≈ 100–350 m) off the wall and keep moving. 6 runs (180 d)
    write to `<outputs>/tmp_wall_velocity/run_<regime>_eps<mm>mm.zarr`,
    ~1 h each; metrics done (results below). Running (2026-11-03 evening): ε = 1/2/3/5 mm/s with v3 walls + reject, both regimes, plus 024d→024e→024f survival weighting on a tmp output root (`tmp_wall_velocity/t_sweep_compare.*`). Then authors pick ε.
  - Metrics: stalled share of particle-hours (step < 1 m / 10 m) at
    30/90/180 d and drifters stalled at day 180; pinned (≤ 1e-3 cell of a
    wall); wall-distance distribution in wall cells; coastal residence
    (share in cells with ≥ 1 wall); land entries (within grid / at nest
    switch); crow-flies distance quantiles; overlap of coastal-contact and
    day-180 position distributions at 2 / 6 / 15 / 30 km between ε runs;
    runtime.
- Prototype results (2026-11-03, `tmp_wall_velocity/e_metrics.{py,log,csv}`,
  figs `fig/e_overview.png`, `fig/e_trapped.png`; water-seeded = release
  cell has a non-wall face: 14,230 surface_stokes / 14,240 bottom):

  | surface_stokes, ε = 0 / 1 / 5 mm/s | 0 | 1 | 5 |
  |---|---|---|---|
  | pinned particle-hours 0–180 d / drifters at d180 | 65 / 95 % | 0 / 0 | 0 / 0 |
  | cell-trapped (≥ 7 d in one wall cell + 4-nbrs) drifters d30/90/180 | 50/93/99 % | 32/86/77 % | 24/76/49 % |
  | median final trap | 131 d | 41 d | 45 d |
  | trap episodes ≥ 7 d ending before d180; of those reaching 20 km in 30 d | 56; 74 % | 82; 71 % | 88; 77 % |
  | particle-hours in wall cells 0–180 d | 88 % | 82 % | 77 % |
  | land: drifters entering (within grid / nest switch); particle-hours | 181/128; 1.35 % | 573/0; 0.58 % | 479/0; 0.12 % |
  | distance from release d180, median | 165 km | 310 km | 351 km |
  | overlap of d180 positions vs ε = 0 at 2/30 km (noise floor .89/.93) | – | .10/.44 | .01/.28 |

  Bottom: pinned 24 % → 0 at d180; cell-trapped at d180 43/25/16 % (the
  metric also counts slow bottom drifters that never pinned); d180 median
  distance 59/68/71 km; d180 overlap vs ε = 0 at 30 km .78 / .65 (floor .91).
  - Reading: ε ends pinning and frozen land residence, and trap episodes
    now end. Under persistent onshore forcing surface_stokes drifters sit at
    the off-wall fixed point (0.02–0.09 cell, 17–31 m/h jitter, ~1 cell)
    — the corner fixed point is a stable node while both components are
    onshore. Dispersal still depends strongly on ε (ε = 1 vs 5 d180 overlap
    .06–.52, floor ~.9): ε is not a negligible numerical knob.
  - Holes: (a) wall faces with wet cells on both sides get no ε (fine: 3
    surface_stokes / 20 bottom); 2 bottom drifters sit there ~20 d at ε = 5.
    (b) Within-grid land entries rise at ε > 0 (untraced); exit from a land
    cell takes ~cell size / ε (median 322 h at ε = 1).
  - (a) explained (`tmp_wall_velocity/f_thin_walls.log`): all 23 are sills
    between tidal-flat / shallow cells (H0 ≤ ~1 m; Wadden Sea, Little Belt)
    that flood from opposite sides; bottom files are NaN there all year.
    Fix (agreed 2026-11-03): ε toward the deeper cell (larger H0); tie →
    cell with more open faces. Done (`wall_*_v2.nc`, `a4_v2_*.log`): 3 / 20
    faces assigned, sign check OK, the 2 stuck bottom drifters leave within
    1 h (`g_stuck_eval.log`).
  - (b) traced (`h_entries.log`, `l_stages.log`, `m_run_walls.log`, 436
    entry hours rerun bit-identically, all stages traced):
    - Most counted entries are clones: drifters collapse onto
      bit-identical positions at the fixed points (within-grid events
      193 / 915 / 710 → 104 / 216 / 116 distinct at ε = 0 / 1 / 5). Distinct
      surface_stokes positions at d180: 15,231 / 11,091 / 8,895 of 21,800 —
      the effective sample shrinks with ε (deterministic advection, no
      diffusion).
    - Main cause (58 / 68 / 84 % of distinct entries): "intermittent land" —
      the 12-day wall sample marks 153 fine cells (303 faces; Wadden Sea,
      Danish straits) as land that do flood in the run window. Drifters flow
      in while flooded and, once dry, exit only by ε (median 129 h at ε = 1).
    - Rest: the ε = 0 convex-corner RK4 stage leak (PLPL / PPPL), not caused
      by the wall field (0.4–1.8 % of the step); exit 1–2 h at ε > 0.
  - Fix candidates: (1) wall set from every step of the forcing period run
    (`wall_*_surface_stokes_v3.nc`, run window only: fine walls −303, thin
    walls 0); (2) reject a step that ends in a land cell (cell-constant
    `Land` field; holds drifters at flooding cells if used without (1)).
    Test of v3 / v3 + reject / v1 + reject (d40–60, ε = 1 / 5) in
    `<outputs>/tmp_wall_velocity/fixtest/`, eval `q_fixtest_eval.py`.
  - Fix test result (`q_fixtest_eval.log`; 14,233 water-seeded drifters,
    d40–60; "genuine land" = v3 land set):

    | ε = 1 / 5 mm/s | base (v1) | v3 walls | v3 + reject | v1 + reject |
    |---|---|---|---|---|
    | genuine-land entries | 17 / 2 | 17 / 2 | 0 / 0 | 0 / 0 |
    | particle-hours in intermittent cells | 29.7k / 5.4k | 86.7k / 11.3k | 86.6k / 11.3k | 0 / 0 |
    | drifters with a rejected step | – | – | 22 / 11 | 234 / 142 |
    | rejections next to intermittent cells | – | – | 0 % | 95 / 97 % |
    | trajectories identical to base at d20 | 100 % | 94 / 95 % | 94 / 95 % | 98 / 99 % |

    v3 alone leaves the corner leaks; reject alone (v1) holds drifters at
    flooding cells. v3 + reject: no land entries, rejections rare (max
    ~15 h per drifter). Cost: drifters now sit in tidal-flat cells while
    they are dry (no ε there, faces not walls) — ×3 particle-hours there,
    ~1.3 % of the total at ε = 1. Recommended: wall set from every step of
    the full forcing period per regime + reject step on its land set.
- ε bounds (2026-11-03, `tmp_wall_velocity/r_onshore_timescales.{py,log}`;
  2019 hourly, 3,000 wall cells per nest, wall-normal velocity at the open
  face opposite the wall = the U in xsi* = ε/(ε+U)):
  - Cells: fine 0.65–0.91 × 0.93 km, coarse 3.9–5.4 × 5.6 km. ε moves
    86 m/d per mm/s at the wall; relaxation across a cell by ε alone
    Δ/ε = 10 d fine / 58 d coarse at 1 mm/s, 2 d / 12 d at 5 mm/s.
  - Time scales: wall-normal ACF e-folding (median) Stokes 24–26 h,
    current 15–17 h, current + Stokes 19–21 h, bottom 10–11 h. Onshore spells
    (current + Stokes) median 5–6 h, q90 1 d, q99 4.5 d; Stokes alone median
    8–14 h, q99 7–8 d.
  - |U| median current + Stokes 5.2 cm/s fine / 8.4 cm/s coarse; median rest
    point xsi* during onshore hours 0.011–0.014 / 0.05–0.07 / 0.18–0.22 cell
    at ε = 1 / 5 / 20 mm/s (bottom: 0.03 / 0.15 / 0.41).
  - Bounds: lower ~1 mm/s (int16 storage resolves 0.5 mm/s); upper from
    ε·T_forcing ≲ ¼ fine cell with T ≈ 1 d → ~2.5 mm/s, and xsi* ≲ 0.1 for
    bottom → ~3 mm/s. Constant ε in m/s gives the same fractional rest
    point on both nests; ε ∝ Δ would equalise the relaxation time but push
    the coarse rest point ~6× further off — keep constant ε.
- Sweep result (2026-11-03 night, `tmp_wall_velocity/t_sweep_compare.{py,log,csv}`,
  `e_metrics_sweep.*`, `fig/t_sweep.png`; v3 walls + reject, 180 d; 024d
  harness on tmp roots reproduces the production sidecar bit for bit on the
  eps0 subset; survival weighting surface_stokes only — no production bottom
  sidecars; members prod `tf1_tw0`, inert `tf1_tw1`):

  | surface_stokes, ε (mm/s) | 0 | 1 | 2 | 3 | 5 |
  |---|---|---|---|---|---|
  | pinned drifters d180 / land entries | 96 % / 186 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
  | cell-trapped d180 / median final trap | 99 % / 131 d | 77 % / 42 d | 85 % / 42 d | 73 % / 38 d | 51 % / 45 d |
  | d180 distance median | 165 km | 313 | 319 | 347 | 350 |
  | afloat weight (prod) d30 / d90 / d180 | 68 / 54 / 51 % | 64 / 38 / 31 | 64 / 41 / 35 | 63 / 35 / 31 | 61 / 32 / 26 |
  | stranded weight (prod) | 49 % | 69 | 65 | 69 | 74 |

  Bottom: pinned 24 % → 0, land entries 89 → 0, cell-trapped d180 43 → 25 /
  23 / 19 / 16 %. v3 + reject vs v1 at equal ε: dispersal unchanged (d180
  overlap .945, floor .946) — the land fix does not move the result.
  - ε = 0 → ε > 0 is a large systematic shift: unpinned drifters keep
    meeting stranding shores (afloat d180 51 → ~31 %).
  - Within ε = 1–3, overlaps between ε pairs stay below the same-run
    half-vs-half floor (hex .74–.77 vs .95; 30 km .81–.83 vs .96; subbasin
    .93 vs .99), weighted or not; non-monotonic in ε. Survival weighting does
    not make them converge.
  - Systematic, not intrinsic (perturbation test, `y_perturb_compare.*`,
    `fig/y_perturb.png`): ε = 2 vs 2.0002 mm/s and vs ≤ 1 m release jitter
    sit at the half-vs-half floor (surface_stokes unweighted occupancy
    hex .94/.96 vs floor .96; afloat d180 34.9 % ± 0.06 pt vs 26–35 %
    across ε). Same-drifter median separation from ε2 at d180: < 0.01 km
    for the perturbations, 56–80 km for ε1/3/5. A divergent minority exists
    (pertEps: 23 % of drifters > 1 km apart at d180, q90 120 km), and 18 %
    collapse bit-identically onto ε2's coastal fixed points — the bulk
    contracts. Bottom: same picture, smaller ε gaps (hex .87–.89 vs .93);
    subbasin level indistinguishable. So ε choice moves results beyond
    noise, non-monotonically (ε2 has the most afloat weight); cause
    untested.
  - 010 integration note: 024d's real-drifter test (zero first-step
    displacement) fails for land seeds once ε acts on them without reject;
    flag land seeds from the land set at release instead.
- Step 3: production integration (design, 2026-11-13):
  - New input stage `notebooks/005_BuildWallFields.py` +
    `scripts/005_BuildWallFields_job.sh`, one run per regime (Dask over every
    c_file of the regime, fine ∩ coarse stems as 010 loads them, all 24
    steps per file): `<output_root>/wall_fields/walls_<nest>_<regime>.nc`
    with `uvel`/`vvel` (unit ±1 on signed wall faces, v3 rules incl. thin
    walls) and `land` (1 on non-wet T-cells). Ported from
    `tmp_wall_velocity/a_wall_field.py --version v3`.
  - 010: parameters `wall_velocity_mm_s = 2.0`, `reject_land_steps = True`,
    `wall_fields = "wall_fields"` (relative to `output_root`). One kernel:
    RK4 with `UV + UVwall` at every stage (ε = the wall fields' scaling
    factor) and the Land reject guarded by the `reject_land_steps` fieldset
    constant. ε = 0 with `reject_land_steps = False` is the old kernel, bit
    for bit. Walls read by `regime` (bottom_vs0p8 uses the bottom walls);
    ε is not multiplied by `velocity_scale` (a numerical device, not
    physics; 034 then compares the physical slowdown alone).
  - 010 writes `nrej` (cumulative rejected steps) and `land_seed`
    (`to_write="once"`, Land at the release position, sampled through the
    Parcels NestedField so nest choice and cell search are the kernel's).
    Every consumer of the zero-first-step land test (021, 022, 023, 024,
    024b, 024d, 036) reads `land_seed` instead: once ε acts on land seeds
    the old test misclassifies them, and with v3 walls a seed in a cell
    that is dry at release is a water seed.
  - Status 2026-11-13 evening: code committed (005, 010, downstream
    land_seed, docs/wall_velocity.md). Walls built
    (`<outputs>/wall_fields/`, 14,582 stems 2016-01-01..2025-12-31, ~9 min
    per regime on 48 one-core tasks): fine surface_stokes 9,664 signed
    faces / 122,564 wet cells, fine bottom 9,662 / 122,550 (differ by 52
    faces, 14 land cells), coarse 4,044 / 34,053 in both; no thin walls.
    005 on the prototype run window reproduces the prototype v3 files bit
    for bit; the full set opens 68 (ss) / 79 (bottom) fine cells the
    window set called land (`tmp_wall_integration/e_wallset_diff.log`).
    Validation (`tmp_wall_integration/`, logs there): ε = 0 without reject
    reproduces production bit for bit (release 2019-10-18, all 87,200
    trajectories × 72 h, surface_stokes and bottom; `b_compare_prod_*.log`).
    `land_seed` ⊂ the old zero-first-step set; the old test additionally
    flagged 173 (ss) / 205 (bottom) water seeds at rest in step 1. Job
    memory set to 16G per task for both regimes (pinning no longer makes
    surface_stokes output compress better than bottom).
  - Review fixes (2026-11-13 night): land seeds skipped in the kernel (ε = 0
    bit identity re-checked, `b_compare_prod_skip.log`), forcing-span
    assert, `surface` job retired; memory 12G surface_stokes / 15G bottom
    (fit: peak ≈ 1.8 GiB + 3.5 × zarr GB; ε = 2 grows the 180-d prototype
    zarr 1.40× ss / 1.05× bottom). Job-script e2e (surface_stokes 2019
    release 58, 10 d) + 024d: no land entries, no pinned drifters, 024d
    real set = `land_seed == 0` = numpy land-set lookup
    (`f_check_024d.log`); before the skip 6 land seeds within ~ε·dt of a
    wet face stepped out of land.
- Launch recipe (step 4; assumes `<outputs>/Trajectories/`,
  `BeachingForcing/` and the derived HexAggregates stores are moved aside,
  the 024a key `HexAgg_key_r6000m.*` stays; walls in `<outputs>/wall_fields/`):

  ```bash
  cd /gxfs_work/geomar/smomw122/2025_fucus-dispersal
  ids=()
  for r in surface_stokes bottom; do
      for y in $(seq 2016 2024); do
          ids+=($(scripts/submit_010.sh $r $y | tail -1 | cut -d' ' -f1))
      done
      ids+=($(scripts/submit_010.sh $r 2025 29 | tail -1 | cut -d' ' -f1))
  done
  vs=$(RELEASE_NS=0,7,12,18,24,31,37,43,49,55,61,67 VELOCITY_SCALE=0.8 \
       scripts/submit_010.sh bottom 2019 | tail -1 | cut -d' ' -f1)
  # 024/024b/024c per (regime, year) -> 024d -> 024e/f/g -> 025-031 + 035
  # (035 in the second wave, after 026 and 029) -> 032 headline
  scripts/submit_downstream_chain.sh "${ids[@]}"
  # bottom_vs0p8 is not discoverable by the chain (its job script is bottom's)
  ALLOW_MISSING=1 scripts/submit_aggregates.sh bottom_vs0p8 2019 6000 \
      --dependency=afterany:${vs}
  # then by hand once their inputs exist: 034 (bottom + bottom_vs0p8 2019
  # counts), sbatch scripts/036_SeedingMap_job.sh 6000 (surface_stokes 2019
  # zarrs), 033 after 024e (scripts/033_ShorelineTypeMap_job.sh, its
  # header example)
  ```
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
