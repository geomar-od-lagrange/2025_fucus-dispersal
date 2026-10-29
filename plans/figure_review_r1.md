# Figure review R1: author remarks on the 2026-10-20 draft

The authors' remarks (2026-10-27) on the figures of `paper/main.tex` at
`baa633c`, turned into work items. Figure numbers are in the order the figures
appear in that draft; the remarks are quoted in short form.

## Figure map (draft `baa633c`)

| Fig | paper file | notebook |
|---|---|---|
| 1 | `data_seeding_map` | 036 |
| 2 | `r4_where_stranded_ALL` | 029 |
| 3 | `r4_beached_fraction_curve_ALL` | 029 |
| 4 | `r4_shoretype_vs_stranding_ALL` | 033 |
| 5, 6 | `r1_bottom_ALL_{percentage,dilution}` | 026 |
| 7 | `r1_area_timeseries_ALL` | 035 (areas from 026) |
| 8, 9 | `r1_surface_ALL_{percentage,dilution}` | 026 |
| 10, 11 | `r1_meanage_{surface,bottom}_ALL` | 026 |
| 12, 13 | `r1_bottom_slowdown_{length,map90}` | 034 |
| 14–29 | `r2_<origin>_{surface,bottom}` (16 origins) | 026a |
| 30–33 | `r3_emission_{surface,bottom}_ALL_T{30,180}d` | 028 |
| 34 | `r5_drifting_fraction_curve_ALL` | 030 |
| 35 | `r5_afloat_maps_ALL` | 030 |

## Work items

**Removals.** Fig 5 and Fig 8 (percentage maps): the dilution style (Fig 6,
Fig 9) is the reference for every dispersal map. Fig 14–29 (per-origin 2×2
pairs) give way to W5.

- **W1 — seeding into the stranding map (Fig 1 + Fig 2).** Fig 2 carries the
  release cells (circles) plus unannotated HELCOM subbasin outlines; no
  colour coding by subbasin. Colourbar at the bottom, full page width.
  Caption must state what a colourbar value means in plain terms (Q1).
  Fig 1 as a separate figure goes away (fate of its land-share panel: Q0).
- **W2 — Fig 4 (shore type vs. stranding).** Answer Q2, then either remove it
  or rebuild it around a substantial explainer.
- **W3 — Fig 7 (area time series).** Area holding 99 % instead of 90 %, if
  that removes the rise-then-fall of the surface curve (Q3).
- **W4 — Fig 10, Fig 11 (age maps).** Median age instead of mean age, from the
  `age_bin` axis. Fix the highlighted points in the German Wadden Sea (Q4).
- **W5 — per-origin maps (replaces Fig 14–29).** One multipanel figure with
  one log-dilution map per origin subbasin (16 panels), surface regime,
  180 d horizon. Same for the bottom regime at 180 d, and optionally the
  surface regime at 30 d.
- **W6 — Fig 12 (bottom slowdown length).** Plot the 0.8×/1× ratio, not two
  absolute lines.
- **W7 — Fig 13 (bottom slowdown maps).** Stack the two panels vertically so
  each is wider.
- **W8 — Fig 30–33 (emission matrices).** Caption/FACTS: how beached
  particles count in the survival-weighted connectivity (Q5).
- **W9 — Fig 35 (afloat maps).** A stronger, possibly discrete, colormap.
  This overrides the AGENTS.md default-colormap rule; record the rationale
  in `docs/visualisations.md`.

## Open questions

- **Q0** — Does the land-seeded share panel of Fig 1 survive (supplement), or
  is it dropped?
- **Q1** — What a value of e.g. 1e-3 on the Fig 2 colourbar means physically.
- **Q2** — Does stranding occur at near-zero flat fraction under the
  production member (`trap_flat = 1`, `trap_wall = 0`), and why?
- **Q3** — Is the 90 % surface area curve non-monotone because of survival
  weighting, and does 99 % fix it?
- **Q4** — What are the highlighted points in the German Wadden Sea in
  Fig 10: release-site markers drawn over the map, or real data?
- **Q5** — How do stranded particles enter the 024g weighted connectivity?

## Findings

Diagnostics in `tmp_figreview_dispersal/` (untracked, replicates 026 exactly;
90 % areas match `Exports/026/026_areas_*.csv`).

- **Q3:** the decline is not a survival-weighting artefact — the unweighted
  surface curve falls too (90 %: 199k → 97k km² from T = 20 d to 180 d). The
  surface cloud concentrates with age: per 10-d bin the 90 % area drops from
  194k km² (0–10 d) to 12k km² (170–180 d), with the 200 densest hexes going
  from 26 % to 84 % of the weight, mostly partly-land coastal hexes (stalling
  against the shore under onshore Stokes). Bottom does not concentrate and is
  monotone at every threshold. At 99 % the pooled surface curve still peaks
  (T = 50 d, 378k km²) and falls 6 % to 355k km² by 180 d, and every
  single year is non-monotone. So 99 % flattens the decline but does not
  remove it; a curve that cannot fall needs a different metric (e.g. the
  cumulative footprint of hexes reached, untested).
- **Q4:** the Wadden Sea dots are real hexes, not markers — three hexes in
  the Elbe and Ems estuaries with mean age ≥ 170 d and almost no weight (1,
  49, 213 of 9.2e10), passing 026's map floor because that floor is on
  dilution (1st percentile, `026:564-565`) and their water area is tiny.
  Fix: floor on share, not dilution. Median age from the 10-d `age_bin`
  CDF is ~5 lines next to `026:556`; surface median is 10–15 d younger than
  the mean across the Baltic Proper (median − mean p50 = −8.7 d), bottom
  within ±5 d; spatial pattern unchanged (Spearman ≥ 0.987).

Beaching diagnostics in `tmp_figreview_beaching/` (untracked).

- **Q1:** Fig 2 value per hex = `100 · stranded weight / total released`
  (`029:321`, `:341-348`), deposits at any age < 180 d. Production member,
  2016–2024, season ALL: 37,296,956 water-seeded drifters released, 54.3 %
  stranded, 1,221 hexes receive weight. 1e-3 % = 1 in 100,000 ≈ 373
  drifters, i.e. ~10 per million released strand on that 6 km stretch
  within 180 d. Peak hex ≈ 0.98 %; hexes ≥ 1e-3 % hold 99.9 % of stranded
  weight. 110 hexes sit below the four-decade floor without an `extend`
  arrow — fix in W1.
- **Q2:** at tc = inf the rate is `ff/100 · r(w_on) / 3 h` (`024e:229-250`),
  so ff = 0 strands nothing (verified on 2019 sidecars). The x ≈ 0 column
  in Fig 4 comes from (a) small-but-nonzero ff integrating over 180 d and
  (b) a definition mismatch: the rate uses the nearest 500 m land cell's ff
  (`024d:301-323`), the deposit is binned at the particle's water hex
  (`024e:321`), and 033's hex ff is the length-weighted mean of sub-segments
  in that hex (`033:264-286`). Hexes with hex ff < 0.05 receive 1.4 % of
  stranded weight; 83 % of stranded weight has ff ≥ 50 at deposit. With the
  rate ∝ ff, the plot mostly hands back the input assumption → **remove**
  (supplement only with the explainer and the `tcinf` substrate-inert
  member alongside).
- **Q5:** 024g weights particle-hours by `exp(−A)` at the particle's current
  hex (`024g:288-303`); stranded mass vanishes from particle-time, has no
  sink column, and 028 row-normalises (`028:334-336`), so Fig 30–31 show
  shares of *surviving* particle-time and per-origin beaching loss
  (surviving share 0.10–0.92 at 180 d) is normalised away. No double
  counting. Fig 32–33 (bottom) are unweighted. A "stranded" sink column is
  available from 024e's `beach_hex` mapped to subbasins.
- **Q6:** Fig 35 is viridis with one shared LogNorm (vmax 2.02 %, four
  decades), no `extend`; 31–63 % of hexes per panel fall below the floor
  and render as uniform dark slabs. Diverging is not meaningful for
  positive shares. Proposal: discrete `BoundaryNorm` on log-spaced bounds
  (1e-4 … 3 %) with `extend="both"`; right column as the surviving
  fraction surv/occ on discrete bins 0–1 (`030:399-411` already computes it).

## Decisions (authors, 2026-10-27)

- Q0: Fig 1's land-seeded share panel moves to a supplement figure.
- Fig 4: removed (no supplement).
- Fig 7: stays at 90 %; the caption/FACTS explain the decline as coastal
  accumulation of the surviving surface cloud (Q3 numbers).
- Fig 30–31: gain a "stranded" target column.

## Implementation

One agent per notebook group (disjoint files); each commits its own paths,
reruns its notebook(s) on NESH for the paper configuration (release years
2016–2024, season ALL, production member for surface), and reports output
PNG paths. Every saved figure stays print-ready (180 mm, 300 dpi,
`od-cartopy-maps` for map grids). Then one paper pass, then review.

- **A — 029 + 036 (W1, Q0).** 029 where-stranded map: HELCOM subbasin
  outlines (unannotated, no fill), release cells as now, one horizontal
  colourbar at the bottom spanning the full width, `extend="min"` for the
  110 below-floor hexes. 036: a single-panel land-seeded share figure for
  the supplement; the subbasin-coloured panel is no longer needed in the
  paper (036 may keep producing it).
- **B — 026 + 026a (W4, W5).** 026: median age per hex (linear
  interpolation within the crossing 10-d bin of the per-hex cumulative
  weight) replaces mean age in the age maps and exports; the age-map mask
  floors on share, not dilution. 026a: new parameter
  `overview_horizons_days_csv` (default `"30,180"`); per horizon one
  multipanel figure with one log-dilution map per origin subbasin (all
  named origins, 4 × 4 for 16), shared colour scale across panels. The
  per-origin 2 × 2 figures stay as collaborator deliverables.
- **C — 034 + 030 (W6, W7, W9).** 034: length panel shows the ratio
  0.8× / 1× against horizon (single line, 1 marked); map pair stacked 2 × 1.
  030: discrete `BoundaryNorm` on log-spaced bounds (1e-4 … 3 %),
  `extend="both"`, a sequential map with wide lightness range; the right
  column becomes the surviving fraction surv/occ on discrete bins 0–1.
  Colormap override rationale in `docs/visualisations.md`.
- **D — 028 (W8).** Surface (weighted) emission matrices gain a
  `stranded` column so each row adds up to the origin's *unweighted*
  particle-time over [0, T). Stranded particle-time per origin =
  Σ over 024e deposits with `beach_age_bin` < T/10 d of
  `weight · (T − age_mid) · 24 h`, origin from the key's `helcom_subbasin`
  of `release_hex` (same assignment as 024g). This equals
  Σ_hours (1 − exp(−A)), the particle-time removed by the weighting, at
  10-d bin resolution. Unweighted (bottom) matrices unchanged.
- **E — paper.** After A–D: copy figures, drop Fig 1 (→ supplement), Fig 4,
  Fig 5, Fig 8, Fig 14–29; add the 026a overview figures (surface 180 d,
  surface 30 d, bottom 180 d); rewrite affected captions and FACTS (Fig 2
  value translation from Q1, Fig 7 decline from Q3, Fig 10/11 median,
  Fig 30–31 stranded column). Commit inside `paper/` locally; push only
  when told.

## Status

- 2026-10-27: Q1–Q6 answered, decisions taken; A–D launched.
- D done (`a5d3296`, `8ff6804`): per-origin closure −4.8 %…+0.05 % at
  30 d, −0.82 %…+0.01 % at 180 d (028 ALL print); stranded share of row
  0.08 (Bothnian Sea) – 0.90 (Great Belt) at 180 d. 032 Table C drops the
  `stranded` rows (`bc924f7`).
- B: 026 done (`8f0eec5`, median age, share floor at the 5th percentile —
  the 1st kept two estuary hexes); 026a overview code in (`e2d6ef4`),
  production runs checked; overview top at the pooled 99th percentile (`43cd079`;
  the 99.9th left basin interiors in the bottom decade).
- C done (`c5c79ec`, `2d41525`): 034 ratio 0.97–0.98 (length) across
  30–180 d; map pair stacked at 220 mm height on a Baltic crop (9–30.7°E,
  53–66°N), each map ~1.6× the old area. 030 on `magma_r` discrete log bins,
  right column = surviving fraction; Fig 35 caption must say so.
- A done (`69ad70d`, `c287fa0`, docs `a472962`, `114966c`): 029 full-width
  map with shared-edge subbasin borders, bottom colourbar, `extend="min"`
  (figure is 180 × 231 mm, LaTeX scales it to 0.85 textheight); 036
  `LandSeededShare_r6000m.png` for the supplement.
- 026a overview v1 unreadable (scale 1e-26…1e2, North Sea half the panel,
  isobath blobs); redone (`e2d6ef4`, `43cd079`): 99th-percentile top, floor
  top/1e4 with sub-floor hexes masked, Baltic crop, no isobath, geographic
  panel order.
- E done (paper `5360545`, `331fc64`; 35 → 19 figures incl. supplement); 032
  Table C drops `stranded` rows (`bc924f7`).
- Review (2026-10-27) — must fix: seasonal 028/030/026 outputs are still
  the 10-20 versions, so 032 Table C mixes denominators and its Gaps entry
  (`SurvivalFraction_*`, spec at `032:1066-1072`) is a stale consumer;
  paper Table C FACTS stale; row-sum captions must say ≈ with the printed
  per-origin closure (−4.8 %…+0.05 % at 30 d, −0.82 %…+0.01 % at 180 d).
  Also: delete `MeanAge_*`/`026_meanage_*` outputs, rename
  `fig:r5_surviving`, surface overview order T30 → T180, 026a colourbar gap.
  Fixed in pass F (next bullet).
- Fix pass F done: 028/030/026 seasonal reruns (DJF–SON, surface_stokes
  member + bottom unweighted for 026), `MeanAge_*`/`026_meanage_*` deleted;
  032 spec → `SurvivalHeatmaps_*` (`fb16a03`), rerun with 0 gaps, Table C on
  the new denominator throughout; 026a colourbar `pad=0.01` and zero-water
  share printed apart from below-floor (`9789a28`), ALL reruns; 034 comment
  (`34da18c`); manuscript_extension P7 line (`2678249`); paper `c86f721`
  (bump `54c7b51`): Table C FACTS, ≈ row-sum captions with per-origin
  closure, `fig:r5_occupancy_surviving`, overviews T30 → T180 regenerated.
  New: Gdansk Basin puts 3.1 % (30 d) / 4.5 % (180 d) of its particle-time
  on zero-water hexes (026a surface ALL print) — check before quoting.
- W10 done: one extent 9–30.7°E, 53–66°N (`extent_*` params) in 025, 026,
  026a (per-origin + overview, `overview_extent_csv` gone), 026b, 027, 033
  (`79ca1d9`); 3 m isobath is the top layer on every map that carries it,
  back on the 026a overview (`737982a`); only hexes touching the extent are
  drawn (off-extent patches pushed the 026 median-age title off the page,
  `8a0ef03`); grids capped at 230 mm, 025 per-origin grid 4 columns; docs
  `6036658`, `bf29d2e`. Reruns: 026 ALL+seasons × both regimes, 026a ALL,
  025/027/033 ALL, 034, 035 (stale after 026), 032 with no Gaps. 026b
  code-only (~15–22 min per run). Weight outside the extent (max): 026
  surface 0.49 %, bottom 0.004 %; 026a overview surface T30 3.3 %, T180
  6.4 %, bottom T180 0.05 % (all Kattegat origin); 034 0 %. Paper `3871373`
  (bump `97899d2`). Open: 025's German-waters zoom (surface) clips its title
  and colourbar label; 020/023 stay at 5–32°E.

## Addendum (authors, 2026-10-27): one map extent everywhere

- **W10** — every map uses the narrow extent of the stranding map (029):
  9–30.7°E, 53–66°N, via the same four `extent_lon_min/max`,
  `extent_lat_min/max` parameters 029/030/034/036 carry. Affected: 025, 026
  (now the key's total bounds), 026a (per-origin maps and the overview,
  replacing `overview_extent_csv`), 026b, 027, 033. Each map notebook prints
  the share of the mapped weight outside the extent, for the captions' FACTS.
- The purple lines the authors asked about are the 3 m isobath of BSH H0
  (magenta contour, fine and coarse grid) in 026/026a/034; along steep
  coasts it hugs the coastline and reads as a second coastline.
