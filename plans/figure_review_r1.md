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

## Status

- 2026-10-27: plan written; Q1–Q6 answered. Awaiting author decisions on
  Q0, Fig 4, Fig 7 metric, Fig 30 sink column before implementation.
