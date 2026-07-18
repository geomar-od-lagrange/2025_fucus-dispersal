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
- **Shore type** sets `trap`. The nearest land cell fronted by a tidal-flat
  cell reads as `flat` (dissipative, retentive: `trap_flat`); else `wall`
  (reflective: `trap_wall`). The Baltic is wall-dominated; flats concentrate
  in the German Bight.
- **Onshore wave forcing** sets `g`. The onshore Stokes component
  `w_onshore = max(0, −(VSDX, VSDY) · n_out)`, sampled from the **raw
  `baltic_highres` field** by nearest hour and cell — the cross-shore
  transport the `surface_stokes` runs zeroed at blocked faces
  ([stokes_drift.md](stokes_drift.md)), sampled *before* that mask and
  *without* the N=5 spread (which bridges thin land barriers). Positions
  outside the Baltic wave grid (German Bight < 9 °E) get `w_onshore = 0`.
  `g(w) = w / (w + w_half)` is a saturating ramp.

**Land-seeded particles are dropped** (zero first-step displacement), as
`024`/`024b` — ~35 % never enter the water and would strand instantly.
`surface_stokes` is the baseline (beaching driver = the runs' own wave
field); `surface`/`bottom` are sensitivity variants.

## Store schema

`HexAgg_beaching_r<radius>m_<regime>_<year>_mMM.parquet` (one per
`(regime, year, month)`) — a grouped weight table, additive across
`release_doy`/month/year like the other `024x` stores, so `029` pools the
monthly partitions by summing. A single drifter contributes fractional
weight to *every* coastal hex/age-bin it strands in, plus a residual row:

| column | meaning |
|--------|---------|
| `release_hex` | release hex of the drifter (`024a` key space) |
| `release_doy` | release day-of-year of the originating zarr |
| `beach_hex` | hex where the weight stranded; `-1` = never-beached residual |
| `beach_age_bin` | `floor(deposit_age_days / age_bin_days)`; `-1` for residual |
| `shore_type` | `wall` / `flat` at the stranding site; `none` for residual |
| `weight` | summed stranded weight (expected particles) in the group |

The never-beached residual is kept (`beach_hex = -1`, …), and deposits +
residual per source hex sum to its released drifter count, so the beached
fraction is `sum(weight | beach_hex ≥ 0) / sum(weight)` per source hex.

## Parameters and their defaults

`max_float_days = 60`, `band_m = 2000`, `tau0_hours = 24`,
`trap_flat = 2.0`, `trap_wall = 1.0`, `w_half = 0.05` m/s,
`age_bin_days = 10`, `raster_dx_m = 500`. The scheme is deterministic (no
RNG). Totals are highly parameter-sensitive — in the Baltic the beaching
scheme can dominate the answer (Siht et al. 2024) — so these warrant a sweep
reported as a range, not a single number.

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
