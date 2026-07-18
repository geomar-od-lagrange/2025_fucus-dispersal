# Analysis wrapup: dispersal deliverables for the study

Post-first-results plan for the `wr/biomass-scaling-and-connectivity`
branch — the remaining scientific analyses and collaborator deliverables
built on the completed trajectory runs. All heavy aggregation is already
on NESH: `HexAggregates/` holds `counts` + `distance` + `connectivity`
parquet for surface / bottom / surface_stokes × 2016–2019 at `r6000m`,
plus the `024a` key + sidecar. So most items below are **parquet-only**
consumers (the `025`/`026`/`028` family) — no Dask cluster, runnable on a
compute node via `pixi run`.

## Status and sequencing

- **Done:** beaching (§5) — implemented as `024d`/`029`, see
  [../docs/beaching.md](../docs/beaching.md) (design note archived at
  [done/beaching.md](done/beaching.md)).
- **Next (parquet-only, unblocked):** finish connectivity (§3), then the
  relative-density + 3 m maps (§1) and their per-subbasin split (§2). The
  connectivity CSV falls out of §3 into §6.
- **Blocked / deferred:** biomass maps (§4) wait on the Estonian
  production figures (JV); the bottom-velocity sweep (§7) needs a kernel
  change plus fresh runs.

## 1. Depth cutoff and relative-density maps

**Overlay the 3 m isobath on the hex density maps.** Below ~3 m *Fucus
vesiculosus* no longer grows and degrades, so the growth-relevant target
zone is the shallow shelf. The `024a` key already carries `mean_depth_m`
per hex (mean over wet BSH `H0 > 0` cells); draw the `H0 = 3 m` contour
as an overlay on the `025`/`026` heatmaps so the reader sees which
occupancy sits inside the viable band.

**Rescale the maps from absolute counts to relative density.** Today
`025`/`026` plot raw particle-`n_obs` on a `LogNorm`. Switch to relative
units: a **percentage** view on a linear scale (share of released /
arriving particle-time) and a **dilution** view on a log scale (how many
orders of magnitude the source concentration has thinned). This makes
maps comparable across sources and is the substrate for the biomass
scaling in §4.

## 2. Per-subregion percentage and dilution maps

**Emit the percentage/dilution maps split by origin subbasin.** Same
rescaling as §1 but faceted per HELCOM origin subbasin, extending the
existing `026a` (per-origin) / `026b` (per-origin × year) pattern. Keeps
each source's dilution field separate so it can later be weighted by
**regionally varying biomass production** (§4) — a national or per-basin
production estimate multiplies that basin's dilution map.

## 3. Subbasin → subbasin connectivity  *(store + POC matrix already landed)*

**Store and proof-of-concept viz exist; finish the analysis views.**
Commit `b4802f9` added the `024c` residence-connectivity store
(`(origin_subbasin, target_subbasin, release_doy, age_bin) → n_obs`) and
the `028` matrix heatmap (linear + log). Still open in
`plans/subbasin_connectivity.md`:

- Row-normalised **emission-fraction** matrix (`n_obs / row.sum()` — what
  fraction of an origin's particle-time reaches each target; raw rows
  differ by orders of magnitude).
- Age-horizon matrices (cumulative `n_obs` over `age_bin ≤ T`, the
  connectivity analogue of `026`'s horizons).
- Docs: add a Connectivity section to
  `docs/hexbinning_and_connectivity.md` + a `028` entry in
  `docs/visualisations.md`; move the plan to `plans/done/`.
- CSV export of the matrix (see §6).

## 4. Biomass maps

**Turn relative-dilution maps into biomass (kg / m²) maps.** Depends on
the §1/§2 relative maps and a production estimate:

- Typical **biomass production** figures from Estonian *Fucus* sites (via
  JV) as the per-source input rate.
- An **open-water biomass** estimate for drifting *Fucus* (standing
  biomass of the drifting fraction).
- Combine with the per-subregion dilution maps → **biomass (kg / m²)**
  target maps. The `024a` key already carries `fucus_area_m2` per hex
  (REDLIST intersection), usable to tie production to source area.

## 5. Beaching and remobilisation — ✅ implemented, see [../docs/beaching.md](../docs/beaching.md)

**Implemented (weighted scheme).** `024d_BuildBeaching` (heavy pass over the
zarrs + raw `baltic_highres` Stokes → `HexAgg_beaching_*.parquet`) and
`029_BeachingMaps` (parquet-only consumer: where-stranded wall/flat maps,
beached-fraction-per-source-hex, age-horizon maps) both landed; the design
note is archived at [done/beaching.md](done/beaching.md). Beaching uses the
**weighted / fractional** deposition scheme (each particle deposits
fractional stranding weight along its coastal path; deterministic, noise-free
at the high-age tail, composes with a Fucus lifetime `L(t)`) — the retired
stochastic first-stranding draw agreed to <0.2 % on totals. A POC run
(surface_stokes, Aug 2019, defaults) strands ~57 % of drifter weight within
the 60-day viability window, wall-dominated with flat strandings correctly
concentrated in the western Baltic / German Bight.

Open follow-ups (not blocking): a **parameter sweep** of `τ0` / band width /
`trap(flat,wall)` reported as a range (totals are highly parameter-sensitive
in the Baltic); a **Fucus lifetime `L(t)`** replacing the step-function
viability cutoff; and recomputing the free-drifting occupancy/connectivity
weighted by the surviving (un-beached) fraction.

## 6. Data deliverables to send out

**Package the headline products for collaborators.**

- **Connectivity matrices → CSV** (from the `024c` / `028` store; §3).
- **Dilution maps → R-compatible** (e.g. gridded CSV or GeoTIFF the
  recipients can read in R; §1/§2).
- **Biomass map → R-compatible** (same export path; §4).

## 7. Bottom-velocity sensitivity sweep (deferred)

**Deferred — it needs a kernel change plus fresh runs, not just a plot.**
Two costs make this a poor fit for a quick turnaround: (1) it reintroduces
a velocity scale into the `010` kernel (removed in pre-prod cleanup — both
regimes now run plain `parcels.AdvectionRK4`, so the job-script "physics
in the kernel" comment is aspirational), and (2) it re-runs the bottom
regime — a full year is 73 releases × 220-day sims (`--ntasks=73`, up to
24 h on a `base` node), not a quick job.

To make the point cheaply without the full suite: submit a **reduced**
one-year bottom run — a handful of release days at, say, 0.5× bottom
velocity — and compare its spreading against the existing 1× bottom
store. The claim (bottom slowdown barely changes an already-small spread)
doesn't need all 73 releases × 4 years to land.
