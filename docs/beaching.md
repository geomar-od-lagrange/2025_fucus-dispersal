# Beaching diagnostic

A post-simulation estimate of where and when drifting *Fucus* propagules
strand on the coast. Each real drifter carries unit surviving (free-drifting)
weight; inside a near-shore band a **fraction** of that weight strands at the
current coastal hex each step and leaves the drifting pool, the remainder
drifting on. This **weighted / fractional** deposition is the deterministic
expectation of a per-step first-stranding process, with no Monte-Carlo noise
— which matters because releases carry only ~100 particles per seeding cell,
too few for a random-strand rule to resolve the high-age tail where the
free-drifting/beached split is decided. It also composes multiplicatively
with a Fucus **lifetime** `L(t)`: survival is `exp(−∫dt/τ)·L(t)`, today's
`max_float_days` being a step-function `L`. Run as a cheap re-runnable pass
over the trajectory zarrs — the same "aggregate once, consume cheaply" split
the [hex store](hexbinning_and_connectivity.md) uses — so the trapping
parameters can be tuned without touching the physics runs. Design rationale
and literature basis: [../plans/done/beaching.md](../plans/done/beaching.md).

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Build | [`024d_BuildBeaching`](../notebooks/024d_BuildBeaching.py) | trajectory zarrs + raw `baltic_highres` Stokes + `024a` key | `HexAgg_beaching_*.parquet` |
| Consume | [`029_BeachingMaps`](../notebooks/029_BeachingMaps.py) | beaching parquet + key | PNGs under `Figures/029/` |

`024d` is a single-process numpy notebook (no Dask): each zarr fits in
memory and the bottleneck is Stokes I/O, so it processes a partition's zarrs
sequentially (~70 s each). Production partitions per `(regime, year, month)`
— the job fans out the whole `year × month` grid with `xargs`, throttled to
`--ntasks` (so `njobs ≠ ntasks`: request fewer tasks under load without
editing the job). Each task writes its own `_mMM` file; `029` pools the
monthly partitions. `release_month = 0` still builds a single whole-year
file for ad-hoc use, but the pooled store is the monthly partitions.

## The three-ingredient rate model

For each real drifter, positions are walked up to a viability window
(`max_float_days`, capping the biologically live fraction of each 220-day
trajectory). Inside a near-shore band the beaching hazard is a Δt-invariant
e-folding rate (Onink et al. 2021):

```
τ = τ0 / (trap(shore_type) · g(w_onshore))
```

with `trap ≡ 1` as configured, so in practice `τ = τ0 / g(w_onshore)` and
**onshore wave forcing is the only term that modulates the rate** (see the
shore-type bullet for why the factor is kept but left degenerate).

The surviving weight decays by `exp(−Δt/τ)` each in-band step; the weight
deposited at step `h` is the telescoping survival difference
`dep = exp(−A_before) − exp(−A_after)` (with `A` the cumulative
`Σ Δt/τ`). This equals the probability that a first-stranding process
strands at `h`, so the ensemble sum reproduces the stranding field as a
smooth expectation. The surviving weight at the window's end is the
never-beached residual. The three factors:

- **Distance to shore** gates the band. Built as a rasterised
  distance-to-coast field (EPSG:3035, `raster_dx_m`) from the **BSH H0
  land-sea mask** — finite `H0` = water, NaN = land, `H0 ≤ 0` = tidal flat
  ([h0_semantics.md](h0_semantics.md)), fine grid over coarse. This is the
  mask the particles were advected on; the coastline geojson polygons miss a
  large fraction of genuine water positions (their `shapely.contains` drops
  ~36 % of release points) and are not used. `scipy.ndimage` EDT gives the
  distance and the nearest-land index; `∇distance` gives the seaward unit
  normal `n_out`, rotated from the projected grid frame into geographic
  east/north so it dots consistently with the geographic Stokes components
  (LAEA meridian convergence reaches ~17° in the eastern Baltic).
- **Shore type** sets `trap` — **currently degenerate and contributing
  nothing** (`trap_flat = trap_wall = 1.0`), so the rate is uniform along the
  coast for a given wave forcing. The classification still runs: the nearest
  land cell fronted by a tidal-flat cell reads as `flat`, else `wall`, and the
  label is carried into the store. It is kept wired, not deleted, because the
  per-step classification lookup is the expensive part and a real substrate
  dataset should be able to drive `trap` without rebuilding it — set the two
  weights apart to activate. It is kept *inert* because the only typing
  available is BSH's `H0 ≤ 0` tidal-flat flag, which is not a retentiveness
  proxy for Baltic shores: the basin is effectively tide-free, so the flag
  fires essentially only in the German Bight — which is outside the wave grid
  and therefore has zero forcing regardless. Two classes would read as
  resolved coastal morphology while resolving none. Treat `shore_type` in the
  store and in `029`'s summary as a diagnostic label, never as a model result;
  `029` deliberately does not map the split.
- **Onshore wave forcing** sets `g`, and with `trap` degenerate it is the
  *only* term that modulates the rate. The onshore Stokes component
  `w_onshore = max(0, −(VSDX, VSDY) · n_out)`, sampled from the **raw
  `baltic_highres` field** (CMEMS `BALTICSEA_MULTIYEAR_WAV_003_015`, FMI-WAM,
  1 nmi ≈ 1.6 × 1.9 km, hourly) by nearest hour and cell — the cross-shore
  transport the `surface_stokes` runs zeroed at blocked faces
  ([stokes_drift.md](stokes_drift.md)), sampled *before* that mask and
  *without* the N=5 spread (which bridges thin land barriers).
  `g(w) = w / (w + w_half)` is a saturating ramp.

  **The WAM field is extrapolated to full BSH coverage.** WAM's water mask is
  not BSH's: ~20 % (coarse) / ~34 % (fine) of near-shore BSH water cells
  inside WAM's bbox sit on WAM *static land*, and WAM's bbox starts at
  9.01 °E, excluding the German Bight. Left alone these return `w_onshore=0`,
  which for beaching means *rate zero* — a coastline that cannot strand at any
  `τ0`/`w_half`. That is a structural bias, not a parameter choice, and it
  falls hardest on sheltered fjords and archipelago, i.e. prime Fucus habitat.
  So static-land samples read the **nearest static-water cell** instead.

  The distinction that makes this safe is that **WAM NaN means land *or*
  ice** — the wet-cell count varies hour to hour, with ~4.4 % of the grid
  seasonally ice-blanked in the Bothnian Bay and Gulf of Finland. Only
  *static* land is filled, identified by an `ever_wet` mask (finite in any
  hour of a monthly sample spanning the seasonal cycle). **Ice keeps
  `w_onshore = 0`**, which is the physical answer: ice suppresses waves, so
  no-beaching-under-ice falls out for free, consistent with the currents side
  where Fucus simply rides the upper-cell velocity. Fill distance is recorded
  per run and printed — mostly trivial (median ≈ 1 cell, the two coastlines
  disagreeing by one) with a thin tail reaching tens of km into lagoons WAM
  does not represent.

**Land-seeded particles are dropped** (zero first-step displacement), as
`024`/`024b` — ~35 % never enter the water and would strand instantly.
`surface_stokes` is the baseline (beaching driver = the runs' own wave
field); `surface`/`bottom` are sensitivity variants.

## Store schema

`HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_wh<w_half>.parquet` (one per
`(regime, year, month, w_half)`) — a grouped weight table, additive across
`release_doy`/month/year like the other `024x` stores, so `029` pools the
monthly partitions by summing. A single drifter contributes fractional
weight to *every* coastal hex/age-bin it strands in, plus a residual row:

| column | meaning |
|--------|---------|
| `release_hex` | release hex of the drifter (`024a` key space) |
| `release_doy` | release day-of-year of the originating zarr |
| `beach_hex` | hex where the weight stranded; `-1` = never-beached residual |
| `beach_age_bin` | `floor(deposit_age_days / age_bin_days)`; `-1` for residual |
| `shore_type` | `wall` / `flat` label at the stranding site (`none` for residual) — diagnostic only while `trap` is degenerate |
| `weight` | summed stranded weight (expected particles) in the group |

The never-beached residual is kept (`beach_hex = -1`, …), and deposits +
residual per source hex sum to its released drifter count, so the beached
fraction is `sum(weight | beach_hex ≥ 0) / sum(weight)` per source hex.

## Parameters and their defaults

`max_float_days = 60`, `band_m = 2000`, `tau0_hours = 24`,
`w_half = 0.05` m/s, `trap_flat = trap_wall = 1.0` (degenerate),
`age_bin_days = 10`, `raster_dx_m = 500`. The scheme is deterministic (no
RNG). Totals are highly parameter-sensitive — in the Baltic the beaching
scheme can dominate the answer (Siht et al. 2024) — so these warrant a sweep
reported as a range, not a single number. The live knobs are `tau0_hours` and
`w_half` (jointly the rate scale) plus `band_m` (how much trajectory time is
exposed to any rate at all); `trap_*` is not a sweep axis while degenerate.

`w_half` is part of the store filename, so sweep members coexist and
[`031_BeachingSweep`](../notebooks/031_BeachingSweep.py) pools them into a
range. Three things shape how to sweep:

- **`τ0` and `w_half` are not independent.** `τ = τ0 + τ0·w_half/w`: `τ0` is
  an additive floor, the *product* `τ0·w_half` is the weak-wave coefficient.
  Where `w ≪ w_half` only the product matters, so totals trade off along a
  ridge and cannot identify the two separately. Spatial *selectivity* (031's
  Gini of stranded weight) is what discriminates: large `w_half` concentrates
  stranding on wave-exposed coast, small `w_half` saturates `g → 1` and gives
  a pure near-shore-residence map.
- **`τ0` is free.** `A = (1/τ0)·Δt·Σ g(w)`, so `τ0` only scales the exponent —
  caching the `g`-integral would make the whole `τ0` axis a re-reduction with
  no zarr or Stokes I/O. Not yet implemented; `024d` runs one `τ0` per pass.
- **`band_m` is quantised by the mask.** On the 5 km coarse grid (≈70 % of
  release sites) the coastline is a 5 km staircase, so 1–4 km bands all select
  a sub-cell strip of the same first cell ring; only ~5 km reaches a second
  ring. It resolves genuinely only inside the 0.9 km fine nest (western
  Baltic). Sweep it coarsely and report fine-nest and coarse regions apart.

## Limitations

Real beaching is a surf/swash process below the BSH grid and coastline;
where Fucus sits on the coarse (~5 km) grid the band is finer than the
physics, so read the spatial *pattern*, not absolute rates or exact
locations (stranding ≠ source proximity; López et al. 2017). Stranding is
terminal by design (deposited weight leaves the drifting pool and is not
refloated) — beach-cast wrack in the tide-free Baltic is wind/water-level
remobilised (Hammann & Zimmer 2014); a reversible resuspension timescale
bounds the effect and is deferred.

## Cross-references

- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the
  `024x` store pattern + key schema this reuses.
- [h0_semantics.md](h0_semantics.md) — the `H0 ≤ 0` tidal-flat rule behind
  the wall/flat split and the land-sea mask.
- [stokes_drift.md](stokes_drift.md) — the wave field that drives beaching
  and the blocked-face mask this reverses at the coast.
- [visualisations.md](visualisations.md) — `029`'s plot rationale.
