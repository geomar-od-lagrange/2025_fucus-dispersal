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
