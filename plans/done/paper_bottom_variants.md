> Implemented 2026-11-17 — 034 weighting in [docs/visualisations.md](../../docs/visualisations.md)
> §034; the paper diff sits uncommitted in `paper/` for review.

# Paper: every bottom figure in four variants; weighted 034

Authors 2026-11-17. Builds on [done/paper_refresh.md](paper_refresh.md)
(paper c1334a0, main e31e383).

## Decisions (authors)

1. Each bottom figure type in the main text in four variants: regime
   `bottom` / `bottom_ub49` (bed friction, U_b = 49 mm/s) × survival-weighted
   (production member `M = step_wc0p1_ts3_tcinf_tf1_tw0_b750`, "with
   beaching") / unweighted ("without beaching"). Types: 026 dilution
   horizons ALL, 026 median age ALL, 026a origin overview T180d, 028
   emission T30d and T180d. The appendix figures move into Results; the
   appendix goes.
2. 034 (bottom vs bottom_ub49, T90 map pair + length ratio): a
   survival-weighted version beside the unweighted one, both in the main
   text.
3. bottom_ub49 where-stranded map (029) in Results "Beaching", after the
   bottom one.

## Order within a figure type

Weighted first (the pipeline's default, the main-text reading), then
unweighted; inside each, `bottom` before `bottom_ub49`. The reader then sees
friction's effect twice side by side (consecutive figures differ in one
knob) and the beaching effect as the second pair. Considered regime-major
order (bottom w/u, ub49 w/u); rejected, it splits the friction comparison
the 034 figures make across the gap.

## Steps

1. **PNGs on disk** — check each ALL PNG of the four variants exists and
   comes from job 24940971 (025/026/026a, isobath-free) or 24940130 (028,
   which never drew an isobath: a matrix, not a map).
2. **034 weighting** — parameter `weighting` (`"unweighted"` = 024 counts
   `n_obs`; a 024f member tag = survocc `surv`, the released denominator
   from `occ` at age bin 0, as 026). The stranded-% table stays on the 024e
   `member`. Output names carry the weighting tag
   (`BottomFriction_{length,map}_..._<weighting>`; areas CSV too). Job script
   takes it as an argument; the 025-031 driver lists 034 once per weighting.
   docs/visualisations.md §034 and docs/job_scripts.md updated.
3. **Rerun** 034 × 5 seasons × {unweighted, M} (driver, WORK_LIST), verify OK.
4. **Paper** (uncommitted): copy PNGs byte-identically to new stems
   (`r1_*_bottom_{weighted,unweighted}`, `r1_*_ub49_*`, …), group per type,
   drop the appendix, add the ub49 where-stranded map, FACTS from the
   executed notebooks (remaining/released, out-of-scale and no-dilution
   shares, 026a over-top share in the caption, 028 closure / stranded share
   / under-floor, weighted 034 area ratios), weighting-convention bullet,
   `[AGENT DRAFT]` captions, header. Build.
5. Move this plan to `plans/done/` with an outcome line.

## Outcome (2026-11-17)

- Code bca225c (034 `weighting`, job script, driver, docs). Rerun job
  24945894: 034 × 5 seasons × {unweighted, M}, 10 / 10 OK. Stale untagged
  034 outputs removed.
- All 20 four-variant PNGs (026, 026a, 028) found on disk from jobs
  24940971 / 24940130; copied byte-identically with the 034 pair and the
  ub49 where-stranded map (36 paper PNGs, all identical to source).
- Weighted 034 ratios equal the unweighted within 0.002 (length 0.865 /
  0.842 / 0.833 / 0.827 at 30 / 60 / 90 / 180 d); the weighted areas
  reproduce HEADLINE Table F.
- Paper: 13 figures added (ub49 and unweighted variants, weighted 034 pair
  and length, ub49 where-stranded), the appendix removed; 40 pages.
