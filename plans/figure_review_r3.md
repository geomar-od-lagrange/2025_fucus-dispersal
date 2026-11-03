# Figure review R3: author remarks on the R2 draft

Remarks (2026-11-03) on the R2 paper pass (`paper` `2be9796`, 17 figures).
Previous round: [figure_review_r2.md](figure_review_r2.md). Values still
refresh after the coastal-trapping rerun
([numerical_coastal_trapping.md](numerical_coastal_trapping.md)). This round
changes styling only and runs on the current stores. All consumers here are
parquet-only, so runs stay small (the node is shared with the trapping
prototype).

## Figure map (R2 draft)

| Fig | label | notebook | output stem |
|---|---|---|---|
| 1 | `fig:shoreline_type` | 033 | `methods_shoreline_type` |
| 2 | `fig:r4_where_stranded` | 029 | `r4_where_stranded_ALL` |
| 3 | `fig:r4_beached_curve` | 029 | `r4_beached_fraction_curve_ALL` |
| 4, 6 | `fig:r1_{bottom,surface}_dilution` | 026 | `r1_*_ALL_dilution` |
| 5 | `fig:r1_area_timeseries` | 035 | `r1_area_timeseries_ALL` |
| 7, 8 | `fig:r1_medianage_{surface,bottom}` | 026 | `r1_medianage_*_ALL` |
| 9 | `fig:r1_slowdown_map` | 034 | `r1_bottom_slowdown_map90` |
| 10–12 | `fig:r2_overview_*` | 026a | `r2_overview_*` |
| 13–16 | `fig:r3_emission_*` | 028 | `r3_emission_*` |
| 17 | `fig:r5_occupancy_surviving` | 030 | `r5_occupancy_surviving_ALL` |

## Decisions

- **Borders on every map, codes on Fig 1 only.** Fig 1 is the reference map:
  it carries the subbasin codes, a code → name key and the release sites.
  Every other map shows only the borders.
- **Fig 1 shows the value the model uses.** Unattributed coast enters 024d
  as `flat_fraction` 0.5 (docs/beaching.md), so a sequential map coloured by
  the effective per-hex value is accurate. The unattributed share stays in
  the title or caption. A diverging map implies a meaningful midpoint, and
  0.5 is not one.
- **Fig 2 right panel = stranded fraction by origin.** For each release hex,
  the share of its drifters stranded within 180 d. It is computed from the
  024e table: Σ `weight` with `beach_hex ≥ 0` and `beach_age_bin < 18`,
  divided by Σ `weight`, per `release_hex`, pooled like the left panel.
  Confirmed by the authors.
- **Connectivity keeps its own colormap** (authors: distinct from the
  dilution turbo, but with more structure). Emission fraction and dilution
  are different quantities, so they should not share colours.
- **Fig 3 + 5 merge in 035.** 029 exports the stranded curve as CSV and 035
  draws three stacked panels on one shared x-axis. Stranding is cumulative by
  age, so "stranded by age T" and "area at horizon T" share the x-axis.

## Work items

**Cross-cutting**

- **H1 — subbasin borders on all maps:** 025, 026, 026a, 026b, 027, 029, 030,
  033, 034, 036. Use the shared-edge block from 029 (`basin_borders`, grey
  0.5, lw 0.52), copied into each notebook (notebook-local utilities). Draw
  the borders under the isobath and coastline.
- **H2 — codes only on Fig 1:** move `SUBBASIN_CODE` labels (`LABEL_AT`,
  `BOXED`) from 029 to 033; 029 keeps borders only.
- **H3 — full width.** Every map figure is 180 mm wide with no dead margins.
  In `main.tex`, use `width=\textwidth` and drop the `height=` caps that
  shrink maps. Check that each figure fits its page in the compiled PDF. If
  a single map overflows, report it rather than shrink it.
- **H4 — multi-panel map grids fill the width** (026a overview and
  per-subbasin grids, 030, 034 pair, Fig 2 pair). Size the figure from the
  measured axes extent rather than hand arithmetic. The od-cartopy-maps skill
  covers this method, which carries over to these non-cartopy lon/lat axes.

**Per figure**

- **Fig 1 (033):**
  - default sequential colormap on `Normalize(0, 1)`;
  - subbasin borders plus the codes as in 029;
  - release sites (the 029 marker style, documented override);
  - an on-map key "KAT Kattegat …" in the empty land west of the Bothnian
    Sea, so the figure stands alone when extracted. Long names come from the
    key sidecar `subbasin_id_to_name`.
- **Fig 2 (029):** one full-width figure with two panels side by side, each
  with its own colourbar below.
  - Left: the current stranded map, without codes or release markers.
  - Right: stranded % by release hex at 180 d, linear 0–100 %, default
    colormap.
- **Fig 3 (029 + 035):**
  - 029 writes `Exports/029/stranded_curve_{regime}_{season}_{member}.csv`
    with age, pooled %, per-year min and max.
  - 035 draws three stacked panels on a shared x-axis: (a) stranded %, (b)
    area holding 90 %, (c) √area.
  - Panel letters go in the panels. The x label reads "age / horizon T (d)".
  - 035 now runs after 026 and 029; update the pipeline line in AGENTS.md.
- **Fig 7, 8 (026 median age):** non-uniform classes with more resolution
  below 100 d: 0, 5, 10, 15, 20, 30, 40, 60, 80, 100, 140, 180 d,
  `extend="max"`, default colormap resampled. The same bounds go on every
  "days" colourbar.
- **Fig 9 (034):** the 1× / 0.8× maps side by side at full width, with one
  shared colourbar below.
- **Fig 10–12 (026a):** the grid fills the width (H4).
- **Fig 13–16 (028), and the other 028 matrices:**
  - classes at 1–2–5 per decade (12 over 10⁻⁴–1, `extend="min"`);
  - a multi-hue sequential colormap with monotonic lightness (`cmc.batlow`
    from conda-forge `cmcrameri`; fallback `cubehelix`);
  - a documented override.
- **Fig 17 (030):** repeat the "age …" column titles on the second row, so
  the lower row can be cut out on its own.

**Paper pass (after the notebooks)**

- Copy the new PNGs.
- Merge the Fig 3 and Fig 5 floats into one at the Fig 3 position, keeping
  both labels as panel references. Update the `\ref`s.
- Make the caption facts match the new panels (Fig 1 key, release sites;
  Fig 2 panels; merged Fig 3).
- Build, and check page fit (H3).
- Commit inside `paper/` only when asked; push only when told.

## Status

- 2026-11-03: plan written (branch `wr/figures-r3`, worktree
  `.claude/worktrees/figures-r3`).
