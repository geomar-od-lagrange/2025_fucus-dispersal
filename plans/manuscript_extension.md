# Manuscript extension: from the outline notes to the physics part

Intent: extend `paper/main.tex` (the physics part of a larger, joint
manuscript) along the requests in
[`paper/notes/dispersal_f_vesi.md`](../paper/notes/dispersal_f_vesi.md),
keeping to the AGENTS.md §Manuscript rule — humans write the prose, agents
keep `[FACTS]` blocks and figure slots true to the code. Anything the notes
assume that this text does not found (biology, biomass production, carbon)
is the bio part's job and is not chased here.

`notes/` stays as it is: it is the request record. Items below carry the
note's figure/question number so annotations can point at either.

Annotate inline with `> **WR:** …` under any item; open choices are marked
**Decide:**.

## Status

- Phase A done as `paper/` commit `c20edb2` (comment-only: FACTS blocks
  current to `6d854cb`, `[STALE]` flags on four prose sentences). Not yet
  pushed to Overleaf; the superproject gitlink moves only after that push.
- Phases B and C decided yes (WR remarks); Phase D open.
- **Blocking every subbasin-keyed item (P2, P5, P6):** the subbasin rebuild
  in [subbasin_rebuild.md](subbasin_rebuild.md). The 024c/024g stores and
  026a/026b/028/029 figures on NESH are stale until it runs.
- Beaching and survival exist for `surface_stokes` only (the only regime
  with 024d sidecars). That matches the biology — negatively buoyant
  material does not strand — and is stated once in Results, not treated
  as a gap.

## Phase A — bring `main.tex` to the code's current state

Pure fact maintenance, no decisions. One commit inside `paper/`.

- A1. Header: replace "current as of `f8a7996` on `wr/beaching-trapterm`"
  with the current commit on `wr/analysis-wrapup`.
- A2. Release span 2016–2019 → 2016–2025 (2025 partial, 29 releases; 686
  releases per regime = 9 × 73 + 29). [plans/run_extension_2020_2025.md]
- A3. Beaching `[FACTS]`: drop the `tau0` / `w_half` / `s(w)` rate model and
  the "trap ships inert" note; replace with the step model the prose
  equations already show — production member
  `step_wc0p1_ts3_tcinf_tf1_tw0` (`w_c` = 0.10 m/s, `tau_strong` = 3 h,
  `tau_calm` = ∞, `trap_flat`/`trap_wall` = 1/0), `band_m` 2000,
  `max_float_days` 60. [docs/beaching.md]
- A4. Coastline-classification `[FACTS]`: the numbers the prose asks for
  ("GIVE A FEW NUMBERS ON FLAT, STEEP, MIXED FRACTION") from the shoreclass
  table and 033. [docs/beaching.md, notebooks/033]
- A5. Hex/heatmap `[FACTS]`: seasons (DJF/MAM/JJA/SON/ALL) as the temporal
  axis, relative density (% of particle-time, dilution per km²), the
  shoreward subbasin rule. [docs/visualisations.md,
  docs/hexbinning_and_connectivity.md]

## Phase B — research questions (human decision)

Mapping between the two lists; the merged list and its wording are the
authors' call.

| main.tex RQ | notes Q | Results section (proposed) |
|---|---|---|
| — | Q1 surface vs. bottom | R1 Spread by buoyancy |
<!-- Need to quantify this by e.g. a weighted area of the dispersal after a sequence of durations or by median travel distance or similar. -->
| RQ1 where does it end up; RQ4 subbasin differences | Q2 sources/sinks | R2 Where material goes, by origin subbasin |
| RQ4 | Q5 connectivity | R3 Connectivity between subbasins |
<!-- Just as a connectivity table. Resolving different time horizons. (1, 2, 3, 6, 12 months) -->
| RQ3 where it beaches; RQ2 beaching vs. open-sea distribution | Q4 beaching | R4 Beaching |
| — | Q6 not beached | R5 What remains afloat |
<!-- Remaining afloat as a total number per subbasin could be calculated. And then, the heatmaps and connectivity tables only count what managed to float to where it's counted. -->
| — | Q3 seasons | cross-cutting axis in R1–R5, not its own section |

**Decide:** adopt R1–R5 as the Results skeleton?

<!-- YES. -->

## Phase C — Results skeleton in `main.tex`

One `\subsection` per agreed R-item, each with a commented figure slot
(producing notebook + `Figures/` filename stem) and a `[FACTS]` block whose
numbers come from `032`'s `Headline/HEADLINE.md`, never transcribed by hand.
Figures enter `paper/figures/` only once final (Overleaf size); until then
the slot is a comment.

<!-- Yes. -->

## Phase D — figures and analyses, per request

Each item: what the note asks, what exists, what is proposed.

### P1. Seeding map by subbasin (notes Fig 1)

- Exists: release sites overlaid in 029; Fucus area and subbasin per hex in
  the 024a key.
- Proposed: one map of the 872 REDLIST cells coloured by HELCOM subbasin,
  marking the ~35 % land-seeded share. Small new consumer, or a panel of
  033.
- **Decide:** ours, or the bio side's map?

### P2. Surface vs. bottom spread over time (notes Q1, Fig 2)

- Exists: 026 relative density at 10/20/50/100 d, one figure per regime.
- Proposed: both regimes side by side on a shared scale, one row per
  horizon.
- The note asks for 200 d. The viability window is 60 d
  (`max_float_days`), so horizons past 60 d are physics-only. Proposed:
  keep 100 d marked "beyond viability", drop 200 d (run length is 220 d, so
  200 d would also sit at the edge of the data).
- **Decide:** horizon set.

### P3. Area covered, table (notes Q1)

- Exists: nothing.
- "Hexes touched by any particle" grows with particle count. Proposed: area
  of the smallest hex set holding 50 % / 90 % of particle-time, per regime ×
  season × horizon, with interannual min–max. Cheap reduction over the 024
  counts store + `water_area_m2` from the key; lands in 026 (export) and 032
  (table).

### P4. Paths per origin subbasin (notes Q2, Fig 3)

- Exists: 020 raw lines per release subbasin; 026a density per origin
  subbasin.
- Proposed: 026a as the figure (distribution, not sampling); at most one
  illustrative spaghetti panel from 020. `surface_stokes` and `bottom` only,
  as the note says.

### P5. Beaching by origin → destination over time (notes Q4, Fig 5)

- Exists: 024e store carries `release_hex`, `beach_hex`, `beach_age_bin`;
  029 maps stranding; 028 does emission fraction per age horizon on the
  connectivity stores.
- The note's pie-per-cell (Nadal et al. 2024) does not read at 17 × 17.
  Proposed: an origin-subbasin → stranding-subbasin matrix of beached weight
  per age horizon (028's view 2 pattern on the 024e store), on a grouped
  axis of ~6–8 regions for the main text; full 17 in supplement + CSV.
- **Decide:** grouping of the 17 HELCOM level-2 subbasins (proposal to
  follow as a table if adopted).
- Substrate map ("stones vs. sand"): 033, supplement.

### P6. Connectivity by season (notes Q5, Fig 6)

- Exists: 028 × `season` × `member`, unweighted (024c) and survival-weighted
  (024g).
- Proposed: DJF vs. JJA (or SON) on the grouped axis from P5, weighted.
  "Stacked distributions in cells" (note's Geo-thesis S4 reference) only
  on the grouped axis.

### P7. Monthly differences (notes Q3, Fig 4)

- Exists: seasons everywhere; stores are monthly partitions, so a monthly
  view is a parameter, not a rebuild.
- Proposed: seasons in the main text with interannual spread; monthly maps
  as supplement. Notes Fig 4 shows **mean age** per cell, not density —
  computable per hex from the counts `age_bin` axis. The post-reproduction
  biomass loss argues for foregrounding SON.
- **Decide:** mean-age maps wanted?

### P8. What does not beach (notes Q6, Fig 7)

- Exists: 030 drifting-fraction curve; never-beached residual weight kept in
  the 024e store.
- "Not beached" is not "sequestered". Proposed defensible physics product:
  share still afloat at 60 d, and where it is at that point — by subbasin
  and by water depth below it (H0). Bottom regime adds where sinking
  material accumulates. The carbon reading stays in the bio part.

### P9. Units: "% particles everywhere"

- 025–028 already use % of particle-time. Beaching (029) and survival (030)
  report expected particles. Proposed: express both as % of released
  water-seeded particles, so every figure scales to biomass the same way.

## Order of work

1. Phase A (no blockers).
2. P3, P7 mean age, P9 units — no subbasin dependence.
3. Subbasin rebuild, then P2/P4 re-render, P5, P6, P8.
4. Phase C skeleton once Phase B is decided; numbers via 032.
