# Paper: figures grouped by family, parallel across regimes

Authors 2026-11-17. Builds on [done/paper_bottom_variants.md](done/paper_bottom_variants.md)
(paper 9593dc9, main cb4feaf). All figures stay in the main text at their
current size (180 mm; dilution 2×2, median age one map, where-stranded 1×2,
friction pair 1×2, overview 4×4, emission one matrix, occupancy 2×4,
curves as the 035 figure).

## Order

Within a family: S (`surface_stokes`), B (`bottom`), F (`bottom_ub49`);
survival-weighted (+s, `M = step_wc0p1_ts3_tcinf_tf1_tw0_b750`) before
unweighted (−s). Each family sits together in the section it fits; prose
order is not aligned yet.

| # | family | figures | new |
|---|---|---|---|
| 0 | reference: shoreline type (033) | 1 | — |
| 1 | curves (035): stranded %, 90 % area, √area; S/B/F +s | 1 | — |
| 1s | seasons (035 `layout = "seasons"`) | 1 | new |
| 2 | dilution 2×2 (026): S, B, F +s; S, B, F −s | 6 | S−s |
| 3 | friction (034) −s: map pair, length ratio | 2 | +s pair dropped from the paper |
| 4 | median age (026): as 2 | 6 | S−s |
| 5 | origin overview (026a): S, B, F × 30 / 180 d +s; B, F 180 d −s | 8 | B, F 30 d +s |
| 6 | emission (028): (S, B, F) × (+s, −s) × (30, 180 d) | 12 | S−s 30 / 180 d |
| 7 | where stranded (029): S, B, F in Results "Beaching" | 3 | S moves out of Methods |
| 8 | afloat (030, occupancy / surviving 2×4): S, B, F | 3 | B, F |

43 figures.

## Season figure (1s): 035 option, not a new notebook

Same inputs (026's area CSVs, 029's stranded curves), same metric, the same
band/line panel code; only the grouping differs (lines = seasons, rows =
regimes). A new notebook would copy the CSV reading and `series_panel`
(~100 lines) for one figure; a `layout` parameter (`"regimes"` = one season,
one line per regime, the current figure; `"seasons"` = every season, one row
per regime, panels stranded % and √area) keeps one reader. Considered a
second figure always drawn in the ALL run; rejected, a run's output would
then depend on whether another season's CSVs exist. The seasons layout is
drawn once (the job script passes `season=ALL`); the legend names the
seasons.

## Steps

1. Plan (this file).
2. 035 `layout` parameter + job-script argument; docs. Rerun 026
   `surface_stokes` unweighted × 5 seasons (job script) — the 026 area CSVs
   the seasons figure reads exist for +s already. Run 035 seasons.
3. Verify the reused PNGs are current: 026a B/F T30d +s (job 24940971),
   028 S−s (job 24940130, no isobath), 030 B/F (2026-11-15, no isobath,
   no dilution scale).
4. Archive stale outputs (below), after checking no `% Source:` or doc
   points at them.
5. Paper (uncommitted): regroup, add figures and FACTS, drop friction +s,
   move S where-stranded to Results, fix the overview caption, header;
   copy byte-identically; build.
6. Move this plan to `plans/done/`.

## Status (2026-11-17, scope cut by the authors)

Done with existing PNGs (paper uncommitted): friction +s pair dropped;
026a B, F 30 d +s, 028 S−s 30 / 180 d, 030 B, F added; where-stranded S, B,
F together in Results "Beaching"; families regrouped; fig 23 caption fixed.
40 figures + 3 `\missingfigure` placeholders (todonotes).

Deferred items, resumed 2026-11-17:
- S−s dilution / median age: 026 `surface_stokes` unweighted × 5 seasons,
  job 24946249; S−s 180 d overview: 026a, job 24946255 (review round 4).
- Season figure 1s: 035 `layout = "seasons"` (dc02ee6), run in the session;
  folded into the paper.
- Archive: 3 930 files, 4.49 GB moved to
  `<output_root>/archive_pre_b750_20261117/` (relative paths kept,
  `MANIFEST.txt` lists source -> archive): Figures and Exports of 025, 026,
  026a, 026b, 028, 029, 030, 033 whose member tag lacks `_b750` (2000 m band
  `…_tf1_tw0`, pre-750 inert `…_ts3_tcinf`; 2 989 files), the dropped 0.8×
  slowdown (`BottomSlowdown_*`, `034_areas_2019.csv`), and 941 superseded
  executed notebooks (`notebooks_executed/`, stale tags, `034_BottomSlowdown`,
  `034_prelim`, malformed `_6000_r6000m` names). Checked first: no paper
  `% Source:` points at them; docs/beaching.md now points its 2000 m band
  sources at the archive.
