# Figure review R2: author remarks on the R1 draft

Remarks (2026-10-27) on `paper/main.tex` at `3871373` (19 figures). Previous
round: [figure_review_r1.md](figure_review_r1.md).

## Figure map (draft `3871373`)

| Fig | label | notebook |
|---|---|---|
| 1 | `fig:r4_where_stranded` | 029 |
| 2 | `fig:r4_beached_curve` | 029 |
| 3, 5 | `fig:r1_{bottom,surface}_dilution` | 026 |
| 4 | `fig:r1_area_timeseries` | 035 |
| 6, 7 | `fig:r1_medianage_{surface,bottom}` | 026 |
| 8, 9 | `fig:r1_slowdown_{length,map}` | 034 |
| 10–12 | `fig:r2_overview_*` | 026a |
| 13–16 | `fig:r3_emission_*` | 028 |
| 17 | `fig:r5_drifting_curve` | 030 |
| 18 | `fig:r5_occupancy_surviving` | 030 |
| 19 | `fig:supp_land_seeded_share` | 036 |

## Decisions

- **Dilution is a transfer function.** The authors will multiply it by a
  biomass production per release site, so it stays absolute (% of released
  particle-time per km² of water) on one fixed scale for every dilution
  figure, topped at its physical maximum 100 % per km² (all material in one
  km²). Discrete, one class per decade, high-contrast colormap.
- Fig 9 stays a stacked 1× / 0.8× pair (a full-width pair would be ~420 mm
  tall), switched to dilution on the shared scale.
- Connectivity axes use 3-letter subbasin codes; Fig 1 labels each subbasin
  with its code as the key.
- Fig 17 is exactly 100 % minus Fig 2 (45.7 % + 54.3 % at 180 d) → dropped.
- Fig 7's remark "same as 8" read as "same as 6" (Fig 8 is the dropped
  ratio plot).

## Work items

**Cross-cutting**

- **G1 — model coastline everywhere.** Replace Natural Earth with the BSH
  HBMnoku wet-cell coastline `data/bsh_hbmnoku_static/coastline.geojson`
  (what 024a uses for wet area) in every map notebook.
- **G2 — 3 m isobath genuinely on top.** It is still drawn under or
  outlined by the black coastline (zorder 3 vs 2 does not read: the wider
  black line flanks the hairline). Draw the isobath last at a zorder above
  every other artist, and verify on zoomed crops of each paper figure.
- **G3 — one shared colourbar at the bottom of every multipanel figure;**
  the same range for the same quantity across figures.
- **G4 — dilution scale.** Fixed `BoundaryNorm` decades from a floor F to
  100 % per km² (F chosen once from the pooled data so ≥ 99 % of weight is
  above it, same F on every figure), `extend="min"`, sub-floor hexes not
  drawn, discrete high-contrast colormap (`turbo`, one colour per decade;
  documented override). Applies to Fig 3, 5, 9, 10–12.
- **G5 — captions without session lingo** ("production beaching member" →
  "survival-weighted with the stranding model", etc.).

**Per figure**

- **Fig 1 (029):** subbasin borders 1.3 × thicker (0.4 → 0.52); 3-letter code
  per subbasin as label (the key for Fig 13–16).
- **Fig 4 (035):** two panels stacked, shared x-axis, no title clash with
  the y-label; inline labels in the line colour instead of a legend box;
  shading explained in the caption.
- **Fig 6, 7 (026):** discrete 7-class colormap (0, 25, …, 175 d,
  `extend="max"`), same limits on both; colourbar label "days"; max width.
- **Fig 8 (034):** dropped; its ratios stay in FACTS.
- **Fig 13–16 (028):** west → east order (as the 026a overview), 3-letter
  codes on both axes, bigger cells filling the width, discrete log colormap
  with more structure (half-decade classes over the four decades).
- **Fig 17 (030):** dropped; FACTS note it is the complement of Fig 2.
- **Fig 18 (030):** transposed to 2 rows (occupancy, surviving fraction) ×
  4 columns (horizons), full width, one colourbar per row at the bottom.
- **Fig 19 (036):** dropped. In its place the coastline classification map
  (033, flat fraction along the model coast), placed in Methods §Coastline
  Classification.

**Checks**

- **G6 — subbasin coverage.** Count weight with origin or target `-1`
  inside the map extent, and where it is; the Gdansk Basin zero-water share
  (3–4.5 % of its surface particle-time in hexes with zero water area) is
  part of this.

## Subbasin codes (west → east)

KAT Kattegat · SOU The Sound · GRB Great Belt · KIE Kiel Bay · MEC Bay of
Mecklenburg · ARK Arkona Basin · BOR Bornholm Basin · GDA Gdansk Basin ·
EGB Eastern Gotland Basin · WGB Western Gotland Basin · GRI Gulf of Riga ·
NBP Northern Baltic Proper · GFI Gulf of Finland · ALA Åland Sea · BOS
Bothnian Sea · QUA The Quark · BOB Bothnian Bay

## Status

- 2026-10-27: plan written.
- 2026-10-27: **P1 done** (026, 026a, 034, 035; G1–G4, Fig 3–7, 9–12) —
  `bd9a5e4` (026, 034), `3f918b6` (035), `fbf4cb3` (026a), `20088dd` (docs).
  Dilution floor **F = 10⁻⁵ % per km²** (`dilution_floor`; 10⁻⁴ would cut up
  to 13 % of a panel). Weight below F, worst panel per figure: Fig 5 surface
  0.53 % (180 d), Fig 3 bottom 0.34 % (180 d), Fig 10 overview surface T30
  0.47 % (Kattegat), Fig 11 surface T180 0.80 % (Arkona Basin), Fig 12
  bottom T180 0.22 % (Eastern Gotland Basin), Fig 9 pair 0.26 % (1×);
  nothing above 10² anywhere. Coastline: BSH wet-cell outline, fine-first as
  in 024a, minus a 0.001° strip on the fine-bbox seam (otherwise a line runs
  across the Kattegat at 56.45 °N). Isobath at zorder 10, verified on 2×
  crops (magenta over black where they coincide). 034's colourbar stays an
  inset under the lower map (a figure-wide bar overruns the height-limited
  pair). Reruns: 026 ALL/DJF/MAM/JJA/SON × both regimes, 026a ALL × both,
  034 2019, 035 ALL. Notable: at the fixed scale most open water falls in
  the two lowest decades (10⁻⁵–10⁻³), and the darkest turbo class sits close
  to the black coastline; captions for Fig 4 must carry the season and the
  band (per-release-year min–max).
- 2026-10-27: **P2 done** (028, 029, 033, 036; G1–G3, Fig 1, 13–16, new
  Fig 19) — `61da86c`, `a23994a` (028), `05e741a` (029), `b15e018` (033),
  `13992f4`, `9e64c31` (036), `b9d5b86` (docs; it also swept in two
  025/027 "model coastline" doc lines another agent had left unstaged).
  Coastline: same fine-first wet-cell outline, seam strip one fine cell
  (0.0083°) wide. G2 n/a: none of 028/029/033/036 draws the isobath.
  Fig 1: borders 0.52, codes at `fontsize=6` at the representative point;
  KAT, MEC, GDA, ALA nudged off hexes/coast; SOU, GRB, KIE (straits full of
  hexes) on a white alpha-0.8 backing — leader lines from adjacent land were
  tried and crossed more hexes. Fig 13–16: west→east codes, cells fill the
  width (180 × 200 mm with stranded column, 180 × 210 mm without), half-decade
  `BoundaryNorm` edges 1e-4, 3.16e-4, …, 0.316, 1 (eight classes, peak
  rounded up to its half-decade edge, four-decade floor, `extend="min"`) in
  every run; the stranded column sits in the top one or two classes, so each
  cell prints its share in % and the column is labelled `stranded (%)` on a
  top axis. Fig 19 replacement: `Figures/033/ShorelineType_r6000m.png`
  (180 × 231 mm, colourbar below, one-line title). Reruns: 028 ALL/DJF/MAM/
  JJA/SON × (surface_stokes production member, bottom unweighted), 029
  ALL + four seasons, 033 ALL, 036. The 028 raw view (not in the paper) spans
  ~36 decades for the weighted run (tiny `w_obs`); it keeps the full span.
- 2026-10-27: **P1b done** (026, 026a, 034). Decision: dilution is the
  released-normalised transfer function — multiply by the per-site released
  biomass. `dilution` = Σ weight in hex over `[0, T)` / (`N_released` · `T`)
  / water km², `N_released · T` = age-bin-0 samples (sentinels included; per
  origin in 026a) × `T / age_bin_days`; 90 %-area and median age keep the
  cloud-present normalisation. Scale `dilution_floor` 1e-7 to `dilution_top`
  (parameter, default 0.1), two classes per decade, `extend="both"`, every
  hex with weight drawn (under/over colours). Commits `c8fb573`, `343695c`
  (notebooks), `aabcda5`, `5dd8973` (docs). Remaining/released: surface 026
  0.78/0.70/0.66/0.57 at 30/60/90/180 d; 026a surface per origin 0.32–0.99
  (30 d) and 0.10–0.92 (180 d; Gdansk 0.10, Kiel Bay 0.11); bottom and 034
  1.000. Weight below 1e-7 exceeds 1 % in 026 surface 180 d (1.82 %), 026a
  surface T30 Kiel Bay (1.05 %) and 026a surface T180 (up to 4.09 %,
  Gdansk; mean 1.49 %); above 0.1: 026 bottom 180 d 0.52 %, 026a bottom
  T180 up to 3.6 %, 026a surface up to 0.97 %. With sub-floor hexes drawn,
  the 026a surface T180 panels fill most of the basin in the under colour.
- 2026-10-27: **P3 done** (025, 026b, 027, 030, 031; G1/G2, Fig 18, G6) —
  `288c675` (030), `25fa0ac` (025), `771ff9f` (027), `6186945` (026b, 031),
  `925dd82` (docs; two 025/027 doc lines landed in `b9d5b86`), `0f69ed1`
  (032 Fig 18 caption). Coastline: fine-first wet-cell outline unioned on a
  1e-3° precision grid (closes the seam slivers without a strip). Isobath
  (025/026b) at zorder 10, drawn last, checked on 2× crops. Fig 18:
  2 × 4, one `aspect=40` colourbar under each row, 180 × 149 mm
  (2125 × 1757 px). Reruns: 030 ALL + 4 seasons; 025, 027, 031 ALL + 4
  seasons; 026b ALL × both regimes; 035 ALL (stale after the 026 rerun); 032
  last — HEADLINE has no Gaps. 020–023 left on cartopy `coastlines()` (not
  on the paper path, heavy Dask reruns).
- **G6 findings (read-only, `tmp_figreview_r2/g6*.py`).** (a) Subbasin −1,
  2016–2024 ALL: 024g surface member origin −1 0.12 % n_obs / 0.16 % w_obs,
  target −1 2.03 % / 1.82 %. At hex level (024f) target −1 inside the
  extent 1.45 % occ / 1.34 % surv, of which 1.38 / 1.29 % in the Skagerrak
  (W of 12 °E, N of 57.5 °N: 265 of the 306 wet key −1 hexes in the extent,
  top hexes 11.3 °E 58.9 °N, 11.5 °E 58.3 °N, 11.6 °E 58.1 °N, 0.1–0.17 %
  each); outside the Skagerrak only 0.05–0.07 %, mostly one hex at 23.0 °E
  63.7 °N (Quark coast, 10 km² water, 0.04–0.05 %). Bottom (024 counts,
  land seeds excluded — 024c's −1/−1 row is the 35 % land-seed sentinel):
  origin −1 0.12 %, target −1 in extent 0.41 % (0.39 % Skagerrak). Origin
  −1 = four release hexes at 11.0–11.14 °E, 58.8–58.9 °N (Koster area, N of
  the HELCOM Kattegat polygon). (b) Gdansk zero-water share is one hex,
  13128 (18.25–18.43 °E, 54.68–54.78 °N, west shore of Puck Bay): 1.68 %
  of Gdansk occ, 4.5 % of surv, 0.56 % bottom; its particle-time peaks at
  age 0–20 d and fades by 100 d. Every coarse H0 cell in it is land; the
  nearest wet cell (Puck Bay, H0 2 m) lies 0.0015° east. So the key is
  right (no model water); the particles sit on model land within half a
  cell of the wet edge, where the interpolated velocity is still non-zero.
  Not lagoons (Vistula/Curonian are absent from the coarse grid). Other
  origins: ≤ 0.46 % (NBP), overall 0.15 % surface / 0.13 % bottom, 28
  zero-water hexes with weight. **Recommendation:** not a dilated wet mask
  — a half-cell dilation gives 217 in-extent zero-water hexes water but
  inflates coastal hexes' water area by a median 1.38× (p90 5.4×), biasing
  coastal dilution low. Instead add a key column mapping each zero-water hex
  to its nearest hex with water (024a, seconds) and remap `target_hex`
  through it in the dilution consumers (025/026/026a/026b/034; reruns of
  minutes, 026b ~20 min each); or leave as is and state the 0.15 % / Gdansk
  4.5 % in the caption.

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
  coordinates. An independent audit on actual trajectories with the 010
  fieldset is running (`tmp_stall_audit/`).
- Land seeds are excluded correctly (both regimes count water-seeded
  drifters: surface N = 37,296,956, bottom 37,277,160).
- Options, by soundness of the result: (1) floor the dilution area at ¼ hex
  (parquet consumers only, minutes; display fix only); (2) treat
  stopped-for-good drifters as stranded (024d moved-flag + 024e/f/g
  reducers, ~20 core-hours, plus the 024 counts rebuild for bottom; a
  methods decision); (3) fix the corner trap in 010 and rerun everything.
- Status: awaiting the authors' decision; dilution colour range on hold.

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
  or a beaching kernel in 010. Awaiting the authors' decision.

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
