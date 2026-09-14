# Beaching review response

**Implemented.** Current state in [../../docs/beaching.md](../../docs/beaching.md)
(rate model, production setting, sensitivity),
[../../docs/visualisations.md](../../docs/visualisations.md) (028-031 plot
rationale) and
[../../docs/hexbinning_and_connectivity.md](../../docs/hexbinning_and_connectivity.md)
(the connectivity schema). This note records what the PR review asked for and
why each call was made.

Acting on the PR #9 review plus the defects the branch review turned up.
Five decisions, then a mechanical sweep.

## 1. `storm` → `strong wave`

`w_c = 0.10` m/s is the p93 of in-band onshore Stokes — a windy afternoon
with the right fetch, not a meteorological storm. Beaching is observed in
ordinary conditions, and direction and fetch set the onshore forcing as much
as wind strength does. Naming the upper state "storm" imports a weather claim
the threshold does not carry, and it propagates into the physics narrative
("storm climatology", "storm season", "the exposure limit").

Name the state for the forcing that defines it:

| was | is |
|---|---|
| `tau_storm_hours` | `tau_strong_hours` |
| `TAU_STORM_HOURS` | `TAU_STRONG_HOURS` |
| "storm term / rate / hours" | "strong-wave term / rate / hours" |
| "storm season", "storm climatology" | "strong-wave season", "strong-wave climatology" |
| "storm-selective member" | "wave-selective member" |

Member tag keeps `_ts…` — `ts` abbreviates `tau_strong` as naturally as it
did the old name, so tags built before and after read the same. Existing
stores on disk are not preserved, so no compatibility constraint applies.

Out of scope: `docs/stokes_drift.md:7` ("comparable to the Eulerian current
in storm conditions") is a statement about wave physics under actual storms,
not about the rate model's state.

## 2. Delete the saturating rate form

`rate_form = "saturating"` exists only to reproduce pre-sidecar stores, and
no produced data needs preserving. Remove `rate_form`, `tau0_hours`, `w_tau`
and every branch on them from `024e`, `024f`, both job scripts, the smoke
test, and the docs — including the sensitivity table's saturating row.
The two-state hazard becomes the only rate model, so `beaching_exponent`
loses its `else` and the member tag loses its `sat_…` form.

## 3. Say what "`shore_type` is a diagnostic label, never a result" means

Spell out the consequence: the store records `wall`/`flat` at each stranding
site, but `trap_flat == trap_wall` means the label cannot move any weight, so
a difference between the two classes in the output measures the coastline's
own composition and nothing about substrate. It is there to be joined against
a real substrate classification later, not to be reported.

Split the surrounding paragraph, which currently runs the Daily et al.
double-counting rebuttal, the degeneracy argument, and the `shore_type`
statement together.

## 4. Anchor the Gini

State what the number is over and how to read it: `031` computes it over
hexes that receive weight (zeros dropped), so it is concentration *among
receiving hexes*, not over the coast — which is why `beach_hexes` is reported
beside it. 0 = every receiving hex takes an equal share, 1 = one hex takes
everything. It rises with `w_c` (0.57 → 0.70) because a higher threshold
keeps only the wave-exposed shores.

## 5. Connectivity runs `surface_stokes` and `bottom` only

`024c`/`028` and their job scripts default to `regime = "surface"`. Change
the defaults to `surface_stokes` and trim the plan's regime list.

## Defects from the branch review

- **`024d_BuildBeachingForcing.py:446-451`** — the donor BFS dilates with
  `np.roll`, which wraps, so an edge cell can take a donor from the opposite
  edge: a land-blind jump across the domain, exactly what propagating through
  water exists to prevent. It does not bite on the 9–30 °E / 53–66 °N grid
  (every mirrored edge index is inland), but that is geometry, not a
  guarantee. Shift with `-1`-filled slices instead, and soften
  `docs/wam_extrapolation.md`'s "cannot … by construction".
- **`029`/`030`/`031`** match only `_mMM` partitions, so a whole-year
  reducer build (`release_month = 0`) is unreadable. Make the month group
  optional.
- **`plans/done/beaching_sidecar.md`** links notebooks as `../notebooks/…`;
  from `plans/done/` that needs `../../`.
- **`plans/done/beaching.md`** pointer names `024d_BuildBeaching` + `029`;
  the pair is `024d_BuildBeachingForcing` → `024e_BuildBeaching` → `029`.
- **`docs/visualisations.md`** — states `029`/`031` reuse 025's
  `linewidth=0.4` where both compute `hex_seam_lw = 1.1·72/(100·fig_dpi_scale)`
  ≈ 0.26; the `031` section sits after `Cross-references`, which should close
  the file; there is no `028` entry although it overrides `figure.dpi`.
- **`docs/hexbinning_and_connectivity.md`** never mentions `024c`/`028`.
  Add the connectivity section, then move
  `plans/subbasin_connectivity.md` to `plans/done/` with a pointer.
- **`scripts/024e_*_job.sh:22`, `scripts/024f_*_job.sh:23`** — usage comment
  omits the `stagger_max_s` third argument both scripts accept.
