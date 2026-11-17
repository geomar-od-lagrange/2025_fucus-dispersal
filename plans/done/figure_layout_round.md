> Implemented 2026-11-17 (G open) — docs/visualisations.md "Labels in the figure", §033, §034, §037, §038.

# Figure layout round (author remarks on the PDF of paper 8f418d2)

Authors 2026-11-17. Figure numbers refer to that PDF (1 shoreline, 2 curves,
3 seasons, 4–9 dilution, 10–15 median age, 16–17 friction, 18–26
overviews, 27–38 matrices, 39–41 where stranded, 42–44 afloat).

## Decisions (authors)

- **A** Regime and weighting written inside every map / matrix figure, in
  plain words ("surface + Stokes drift, with stranding", "bottom, bed
  friction, without stranding"): suptitle or panel title, no run keys.
- **B** Median age: columns with | without stranding, rows = regime,
  half-width maps, at most two rows per page → figure A (surface, bottom),
  figure B (bottom, bed friction). Same colour classes.
- **C** Emission matrices: one 2×2 per regime (columns with | without
  stranding, rows 30 d | 180 d); three figures. Report if labels do not read
  at half width rather than shrink fonts.
- **D** Friction length ratio (Fig 17): no ratio-1 line.
- **E** Fig 2 regimes unchanged.
- **F** (low) Full-width single-map figures (Fig 1): colorbar as an inset in
  the empty lower right.
- **G** (low) Subbasin borders look broken: find the cause; fix if simple
  and consistent, else report options.

## How

- **A** in the notebooks that draw each family: 026 (dilution 2×2:
  suptitle), 026a (overview suptitle gains the regime), 029, 030
  (suptitle). One plain-language label per (regime, weighting):
  `REGIME_LABELS` + "with stranding" (024f member) / "without stranding"
  (unweighted), as the existing 028 title context does.
- **B, C** need two weightings (and for C two horizons) of one regime in
  one figure; 026 and 028 run once per (regime, member). New lightweight
  consumers reading their exports, by AGENTS.md's tests (different inputs
  per figure, no shared state with the per-run notebooks):
  `037_MedianAgeGrid` (026's `026_medianage_*` GeoJSONs) and
  `038_EmissionGrid` (028's `*_emission_fraction.csv`). Considered a
  multi-member mode in 026 / 028; rejected, a per-run notebook would then
  read other runs' outputs.
- **D** 034: drop the `axhline`.
- **F** 033 (and any other single-map full-width figure): `inset_axes`
  colorbar lower right.
- **G** investigate the HELCOM shared-edge construction.

## Steps

1. Plan (this file). 2. A + D code, reruns (026, 026a, 029, 030 ALL and
seasons; 034). 3. 037 / 038, run. 4. Paper: replace the median-age and
matrix figures, captions, FACTS. 5. F, G. 6. Move to `plans/done/`.

## Outcome (2026-11-17)

- A: 7007fdf (labels in 026 / 026a suptitles, 029 panel titles — a 029
  suptitle left a blank band under its tight-bbox height cut — 030
  suptitle, e03e3a2 height allowance), ce8fb28 (035 "bottom, bed friction").
  D: 7007fdf. B: 037 (f55e470). C: 038 (f9352e4); labels read at half width
  with default fonts; two-digit stranded % fill a 4.2 mm cell, black "0"
  on the darkest classes reads poorly (inherited from 028). F: 4d2c211.
  Docs: c95e7af, 4d2c211.
- Reruns: critical path (season ALL, paper units) job 24946478, 18 / 18 OK,
  plus 029 bottom / ub49, 030 S/B/F, 033, 035 in the session; other seasons
  job 24946597 (85 units) in the background.
- Paper: 31 figures, 36 pages.
- G (borders, not fixed): the borders are the pairwise intersections of the
  HELCOM polygon outlines. They end where the HELCOM polygons meet the
  HELCOM coastline, which is not the BSH model coastline, so in
  archipelagos and coarse-grid coasts a border stops in model water or runs
  over model land; pairs whose outlines do not share vertices exactly give
  point / short-segment intersections (gaps). Options: (1) snap the HELCOM
  polygons to a common topology (shapely `snap` / `set_precision`) before
  intersecting — closes vertex gaps, not the coast mismatch; (2) clip the
  borders to BSH water and extend their open ends to the BSH coastline
  along their last direction; (3) derive borders from the 024a hex
  subbasin assignment (edges between hexes of different subbasins) — on
  the model grid by construction, but hex-stepped. Each changes every map
  notebook (copied block); authors to pick.
