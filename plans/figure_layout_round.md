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
