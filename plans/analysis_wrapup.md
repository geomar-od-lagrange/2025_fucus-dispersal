# Analysis wrapup: dispersal deliverables for the study

Post-first-results plan for the `wr/analysis-wrapup`
branch — the remaining scientific analyses and collaborator deliverables
built on the completed trajectory runs. All heavy aggregation is already
on NESH: `HexAggregates/` holds `counts` + `distance` + `connectivity`
parquet for surface / bottom / surface_stokes × 2016–2019 at `r6000m`,
plus the `024a` key + sidecar. So most items below are **parquet-only**
consumers (the `025`/`026`/`028` family) — no Dask cluster, runnable on a
compute node via `pixi run`.

## Status and sequencing

- **Implemented:** beaching (§5), the relative-density rescaling and 3 m
  isobath overlay (§1), the per-origin-subbasin split (§2), connectivity
  including survival-weighted `024g` (§3), the `season` roll-out across
  every consumer (§5a) and the GeoJSON/CSV exports (§6). Docs:
  [../docs/visualisations.md](../docs/visualisations.md) (season contract,
  relative density, print-ready figures, per-notebook rationale),
  [../docs/hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)
  (024c/024g stores and the 028 views),
  [../docs/beaching.md](../docs/beaching.md),
  [../docs/survival_occupancy.md](../docs/survival_occupancy.md),
  [../docs/job_scripts.md](../docs/job_scripts.md) (how each is submitted).
  Production beaching member `step_wc0p1_ts3_tcinf`; `w_c` is a working
  choice, results are reported as a range.
- **Open:** §4 (biomass maps) is blocked on the Estonian production
  figures (JV); §7 (bottom-velocity sweep) is deferred — it needs a kernel
  change plus fresh runs, not just a plot. Nothing else in this plan is
  outstanding, which is why the file stays in `plans/`.

## 0. Decisions (2026-09-14)

- **Seasons.** Four meteorological seasons by release month — DJF / MAM /
  JJA / SON — plus `ALL` (pooled year). Single consumer parameter
  `season` (str, one of those five) replaces `release_months_csv` /
  `release_month` in every consumer (025, 026, 026a, 026b, 027, 028, 029,
  030, 031). The month set per season lives in one code cell per notebook
  (a dict literal); figure and export filenames carry the season tag.
  Interannual spread (2016–2019) is shown as the error bar / range on
  seasonal panels where the figure allows.
- **Relative density.** `percentage` = 100 · n_obs(hex) / Σ_hex n_obs over
  the selected release set (regime, season, years, age horizon) — share
  of particle-time. `dilution` view = the same share divided by hex water
  area (`water_area_m2` from the 024a key), log10 scale. Absolute counts
  are dropped from the headline maps.
- **Regimes.** Headline products cover `surface_stokes` and `bottom`
  only; `surface` (no Stokes) is not produced further.
- **Exports (§6).** GeoJSON per hex map (hex geometry + value columns)
  and CSV for connectivity matrices, under
  `output_root/Exports/<notebook>/`. No GeoTIFF.
- **Survival-weighted connectivity.** New reducer
  `024g_BuildSurvivalConnectivity` over the 024d sidecar (hex per hour →
  HELCOM subbasin via key; weight `exp(−A)` at that hour), store
  `HexAgg_survconn_r<radius>m_<regime>_<year>_m<MM>_<member>.parquet`
  with columns `(origin_subbasin, target_subbasin, release_doy, age_bin,
  w_obs)`. 028 takes a `member` parameter: empty string → unweighted
  024c store, else the 024g store. Only `surface_stokes` has sidecars.
- **3 m isobath overlay.** Drawn from the BSH static H0 grid
  (`data/bsh_hbmnoku_static/`), not from the hex mean depth.

## 1. Depth cutoff and relative-density maps — ✅ implemented, see [../docs/visualisations.md](../docs/visualisations.md)

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

**As built.** `percentage` is linear and `dilution` log; dilution divides
by `water_area_m2` rather than full hex area so a mostly-land hex is not
credited with a whole hex of volume. The isobath is contoured from the BSH
static H0 grids, not from the key's `mean_depth_m` — a 6 km hex mean
cannot resolve a 3 m contour. One judgement call worth recording: the
linear percentage scale **saturates at the 99th percentile** with an
over-range colorbar arrow. Occupancy is concentrated enough that a handful
of hexes sit ~two decades above the bulk, and scaling to the maximum
paints every other hex the same dark colour; the cost is that the top
percentile is not readable off the map and must be read from the printed
validation numbers or the GeoJSON export. The dilution floor at the 1st
percentile is the same trade at the other end.

## 2. Per-subregion percentage and dilution maps — ✅ implemented, see [../docs/visualisations.md](../docs/visualisations.md)

**Emit the percentage/dilution maps split by origin subbasin.** Same
rescaling as §1 but faceted per HELCOM origin subbasin, extending the
existing `026a` (per-origin) / `026b` (per-origin × year) pattern, plus
the per-origin percentage grid in `025`. Each origin is normalised over
its own releases (cross-origin totals differ by orders of magnitude), and
`026b` normalises each year against its origin's all-year total so the
years stay magnitude-comparable. Keeps
each source's dilution field separate so it can later be weighted by
**regionally varying biomass production** (§4) — a national or per-basin
production estimate multiplies that basin's dilution map.

## 3. Subbasin → subbasin connectivity — ✅ implemented, see [../docs/hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)

**All five follow-ups landed.**
Commit `b4802f9` added the `024c` residence-connectivity store
(`(origin_subbasin, target_subbasin, release_doy, age_bin) → n_obs`) and
the `028` matrix heatmap (linear + log). Remaining:

- Row-normalised **emission-fraction** matrix — `x / x.sum(axis=1)`, one
  figure per cumulative age horizon.
- Age-horizon matrices — half-open and cumulative, ages `< T` days.
- Survival-weighted connectivity — `024g_BuildSurvivalConnectivity` plus
  the `028` `member` parameter (`""` → unweighted 024c, a rate-model tag →
  024g `w_obs`).
- Interannual mean/min/max of the emission fraction over the four release
  years.
- Docs and CSV export (§6), both in place.

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

**Implemented (weighted scheme).** The heavy pass is now
`024d_BuildBeachingForcing`, which caches the rate-model-free ingredients
(onshore Stokes, distance to land, hex id, shore flag, displacement) per real
drifter and hour; `024e_BuildBeaching` reduces it to the stranding store and
`024f_BuildSurvivalOccupancy` to the survival-weighted occupancy store, each
in seconds per zarr per rate-model member. `029_BeachingMaps`,
`030_SurvivalHeatmaps` and `031_BeachingSweep` are the parquet-only consumers.
Beaching uses the **weighted / fractional** deposition scheme (each particle
deposits fractional stranding weight along its coastal path; deterministic,
noise-free at the high-age tail, composes with a Fucus lifetime `L(t)`) — the
retired stochastic first-stranding draw agreed to <0.2 % on totals. The rate is
a two-state hazard, `1/τ = 1/τ_calm + trap·r(w_on)/τ_strong` in band, with
`trap` wired but degenerate (`trap_flat = trap_wall = 1.0`), so onshore wave
forcing is the only term modulating it — see
[beaching.md](../docs/beaching.md).

Remaining follow-ups (not blocking): a **Fucus lifetime `L(t)`** replacing the
step-function viability cutoff; **survival-weighted connectivity** — the
occupancy side is done in `024f`, and the connectivity side is now done in
`024g` (§3). `w_c = 0.10` is a working choice with no observational
constraint; revisit if the paper needs a defended threshold.
`trap(flat, wall)` is *not* a sweep axis while degenerate.

## 5a. Seasonal view as the default product — ✅ implemented, see [../docs/visualisations.md](../docs/visualisations.md)

**Make release season the primary axis of every headline map; pool the
year only as a summary.** The beaching sweep
([../docs/beaching.md](../docs/beaching.md), "Production setting") showed
that release month moves the beached fraction by ~25 points (64 % for
Feb–Apr releases vs 89 % for Sep–Nov at `w_c = 0.10`) while release year
moves it by ~5. Autumn releases meet the strong-wave season inside their
viability window; spring releases do not. Pooling twelve monthly releases
averages two different regimes into a number that describes neither.

The same asymmetry is expected, untested, in the plain dispersal products:
occupancy, distance quantiles, connectivity. Every `024x` store is
partitioned by `release_doy`, so this is a consumer-side change only.

- **Seasons: locked, see §0.** Four meteorological seasons (DJF/MAM/JJA/
  SON) plus `ALL`.
- **Consumers.** 025/026/026a/026b/027/028/029/030/031 gain the `season`
  parameter (§0) and facet by it where the panel count allows; `ALL`
  stays available as the pooled member. 029/030 already partition per
  month, so they only need the facet.
- **Interannual spread** becomes the error bar on each seasonal panel
  (four years), not a separate product.
- **Docs.** `visualisations.md` gets the season contract (which months
  form which season, where the parameter lives); `beaching.md` already
  carries the monthly table.

**As built.** Every consumer takes one `season` parameter and pools all
four release years; the month set is a code-cell dict (papermill injects
only primitives) and the season tag goes in every figure and export
filename. The month is derived from `release_doy` against the partition's
own year, so it is leap-correct and DJF pairs a year's December with its
own January and February. Interannual spread is drawn as an explicit
min–max range (band in 029/030, lo–hi bar in 031, CSV column in 028)
rather than a symmetric `yerr`: a pooled value need not lie inside it
(`beach_hexes` pools as a union of hexes and exceeds every single year).

## 6. Data deliverables to send out — ✅ implemented for §1/§2/§3

**Package the headline products for collaborators.** Everything lands
under `output_root/Exports/<notebook>/`.

- **Connectivity matrices → CSV** — long-form
  `connectivity_<regime>_<season>_<member-or-unweighted>_T<horizon>d_` +
  `{n_obs, emission_fraction, emission_fraction_by_year}.csv` from `028`.
- **Percentage / dilution maps → GeoJSON** — EPSG:4326 hex geometry with
  `hex_id, n_obs, percentage, dilution` (plus `origin_subbasin` /
  `release_year` where faceted) from `025`/`026`/`026a`/`026b`, and the
  distance quantiles from `027`. GeoJSON rather than GeoTIFF because the
  hex tessellation is vector; R reads it through `sf::st_read`.
- **Biomass map → R-compatible** (same export path) — waits on §4.

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
