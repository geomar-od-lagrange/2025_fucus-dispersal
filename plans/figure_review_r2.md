# Figure review R2: author remarks on the R1 draft

Remarks (2026-10-27) on `paper/main.tex` at `3871373` (19 figures). Previous
round: [figure_review_r1.md](done/figure_review_r1.md).

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

## Coastal hot spots → numerical trapping

The coastal dilution hot spots led to an investigation of drifters that pin
against coast walls; findings and the fix live in
[numerical_coastal_trapping.md](numerical_coastal_trapping.md). The dilution
colour range and the ¼-hex area floor wait on that fix.

## Paper pass (2026-11-03)

- `paper` `2be9796` (local): R2 figures copied (executed-notebook code
  matches the committed notebooks), Fig 8 / 17 / 19 dropped, shoreline
  classification map added in Methods as the new Fig 1 (so the code key is
  now Fig 2; numbering above is the old draft's). Captions without session
  lingo; FACTS refreshed; `% NOTE` at Results: values change after the
  trapping fix, styling final. Builds (21 pp); only the three pre-existing
  `TODO_*` citations are undefined. Four image heights reduced to fit.
- Follow-up: internal names inside the PNG titles (028 member tag, 034
  `bottom_vs0p8`) scrubbed in the notebooks and recopied (in progress).
- **For the authors:** (1) dilution top 1e-1 vs 1e-2 (`% TODO(authors)` at
  the first dilution figure); (2) the biomass reading: every dilution figure
  is normalised over the pooled release set (per origin subbasin in the
  overviews), not per release site — "× released biomass per site" needs
  per-site dilution, which no figure shows; wording to decide.
