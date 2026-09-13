# Beaching diagnostic

A post-simulation estimate of where and when drifting *Fucus* propagules strand
on the coast. Each real drifter carries unit surviving (free-drifting) weight;
inside a near-shore band a **fraction** of that weight strands at the current
coastal hex each step and leaves the drifting pool, the rest drifting on. This
**weighted / fractional** deposition is the deterministic expectation of a
per-step first-stranding process: releases carry only ~100 particles per seeding
cell, too few for a random-strand rule to resolve the high-age tail where the
free-drifting/beached split is decided. It composes multiplicatively with a
Fucus **lifetime** `L(t)` — survival is `exp(−∫dt/τ)·L(t)`, today's
`max_float_days` being a step-function `L`. It runs after the physics as a
re-runnable reduction, the "aggregate once, consume cheaply" split the
[hex store](hexbinning_and_connectivity.md) uses, so the rate parameters sweep
without touching the runs.

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Sidecar | [`024d_BuildBeachingForcing`](../notebooks/024d_BuildBeachingForcing.py) | trajectory zarrs + raw `baltic_highres` Stokes + BSH H0 statics + `024a` key | one forcing zarr per trajectory zarr |
| Reduce | [`024e_BuildBeaching`](../notebooks/024e_BuildBeaching.py) | sidecar + key | `HexAgg_beaching_*.parquet` |
| Consume | [`029_BeachingMaps`](../notebooks/029_BeachingMaps.py) | beaching parquet + key | PNGs under `Figures/029/` |
| Sweep | [`031_BeachingSweep`](../notebooks/031_BeachingSweep.py) | several members + key | PNGs under `Figures/031/` |

[`024f_BuildSurvivalOccupancy`](../notebooks/024f_BuildSurvivalOccupancy.py) and
`030` read the same sidecar with the same rate model, recording where the
*surviving* weight is rather than where it leaves
([survival_occupancy.md](survival_occupancy.md)).

## The forcing sidecar

`output_root/BeachingForcing/<regime>/<year>/<trajectory zarr stem>.zarr`,
dims `(trajectory, obs)` with `trajectory` the water-seeded ("real") drifters
only and `obs` hours `0 … window_days·24`. Coords: `trajectory` (int32 index
into the source zarr) and `release_hex` (int32, `024a` key space).

| variable | encoding | semantics | why |
|---|---|---|---|
| `w_on` | float16, Zstd | onshore Stokes (m/s); exactly 0 outside the sampled band | exact for the zeros; max abs error on the beached fraction 2e-8, per-hex ≤ 4e-5 rel. |
| `dist` | uint8, 25 m steps | distance to BSH land; 255 = beyond `band_max_m` or no position | 25 m ≪ `raster_dx_m` = 500 m; nulling beyond the band halves the array |
| `hex` | int32 | `024a` hex id of the position, everywhere; -1 = no position | survival occupancy needs the hex in open water too; 190 k ids do not fit int16; runs compress ~200× |
| `flat` | bool | nearest land is fronted by a tidal-flat (`H0 ≤ 0`) cell | the seam for a substrate map; costs nothing |
| `disp` | uint8, km | crow-flies displacement from release; 255 = saturated | slowly varying, so it compresses like `hex` |

Attrs record provenance and the sampling contract the reducers assert against
(`band_max_m`, `window_days`, `raster_dx_m`, `stokes_fill_max_cells`,
`hex_radius`, `release_time`/`release_doy`/`regime`/`release_year`, source
trajectory and obs counts, `builder`, `git_sha`); chunks are `(10000, nobs)`,
every array written with `numcodecs.Zstd(level=5)`.

All five are **parameter-free with respect to the rate model** — `τ`, the band
width, the viability window, the trap weights and the binning are reductions
over them — so a member costs ~2.5 s per zarr against the ~90 s a re-read of the
trajectories and the raw wave field costs. ~66 MB per zarr at
`window_days = 120` (~19 GB for `surface_stokes` 2016–2019, against 98 GB of
trajectories per regime-year), ~94 s per zarr to build, single-process numpy, no
Dask. Considered and rejected: a sparse (COO) layout — 73–81 % of particle-hours
are in band, dense wins; storing `s(w)` or `Σ s` instead of `w_on` — freezes the
functional form, which is what the sidecar exists to free; zarr's default
Blosc-zstd-bitshuffle — 2–3× larger than raw Zstd on every array here.

`band_max_m = 5000` is the widest band the coarse grid resolves (81 % of
particle-hours, 73 % within 2 km); beyond it `w_on` is exactly 0 and `dist`
reads 255. **Land-seeded particles are dropped** (zero first-step displacement),
as `024`/`024b`. `w_on` comes from the WAM field geodesically extrapolated to
full BSH coverage, without which 78.7 % of in-band samples would read zero
forcing ([wam_extrapolation.md](wam_extrapolation.md)).

## The rate model

Inside the near-shore band the beaching hazard is a Δt-invariant e-folding
rate with two states — a forcing-independent background and a storm term:

```
1/τ = 1/τ_calm + trap · r(w_on) / τ_storm          in band, else 0
r(w) = clip((w − w_c + δ/2) / δ, 0, 1)             step at w_c, optional
                                                   linear edge of width δ
```

Nothing in the data constrains the shape of `r` between calm and storm, so a
shape parameter (Hill exponent, half-saturation) is pure nuisance; the step is
the limiting case of any such form and its one knob means what it says. The
surviving weight decays by `exp(−a)` each in-band step (`a = Δt/τ`), and the
weight deposited at step `h` is the telescoping survival difference
`dep = exp(−A_before) − exp(−A_after)`, `A` the running `Σ a` — the probability
that a first-stranding process strands at `h`, so the ensemble sum is the
stranding field as a smooth expectation. Weight still surviving at the end of
the window is the never-beached residual.

| parameter | meaning | how it is chosen |
|---|---|---|
| `w_c` | onshore Stokes above which stranding is on (m/s) | against the measured in-band forcing (p75 0.057, p90 0.088, p99 0.15 m/s); **the** sweep axis |
| `tau_storm_hours` | e-folding time while `w_on ≥ w_c` | hours; below ~6 h everything in band during a Baltic storm strands and the value stops mattering |
| `tau_calm_days` | background in-band rate regardless of forcing; 0 = off | ∞ or O(1 yr); a year loses 15 % over 60 d, so it is a real axis |
| `delta` | linear edge width around `w_c` | 0 (hard step) by default; `≈0.02` checks against nearest-hour / nearest-cell sampling noise |
| `trap_flat`, `trap_wall` | shore-type factor | deliberately degenerate (both 1.0); the seam for a substrate map |
| `band_m` | near-shore band width (m), ≤ the sidecar's `band_max_m` | 2000; quantised by the mask — on the 5 km coarse grid 1–4 km all select the same first cell ring, so sweep it coarsely and report fine nest and coarse regions apart |
| `max_float_days` | viability window (d); slices the sidecar's `obs` axis | 60; a step-function `L(t)` |
| `age_bin_days`, `disp_bin_km` | store axis granularity | 10 d, 10 km |

Once `τ_storm ≪` storm duration the model degenerates gracefully into an
**exposure model**: the beached fraction is the share of particles ever in band
during a storm hour inside the viability window, and the stranding map is
coastal residence × storm climatology — so the per-month / per-year splits, not
the pooled view, are the primary products. Storm stranding of freshly released
material is signal, not artefact: Fucus releases are coastal, so a particle can
be exposed from `t = 0`. `disp_bin` records the crow-flies travel distance at
stranding precisely so that can be *shown* rather than filtered out — a
diagnostic axis, never a mask. The superseded saturating ramp
`s(w) = 2w/(w + w_tau)`, `τ = τ0/s`, survives in `024e`/`024f` only as
`rate_form = "saturating"`, to reproduce the pre-sidecar stores.

### Provenance

The Δt-invariant e-folding hazard `p = 1 − exp(−Δt/τ)` is standard (Onink et al.
2021; also Hernandez et al. 2024, Daily et al. 2021, Siht et al. 2025). Making
the rate depend on **onshore wave forcing** is this study's extension, with no
precedent in that set — all of them effectively use a constant rate — so the
scheme must not be presented as following Onink et al. Daily et al. object that
a stochastic beaching parameterisation is *meant* to stand in for unresolved
near-shore processes, so feeding a resolved variable back in risks
double-counting it; that does not apply here, because the `surface_stokes` runs
zero the cross-shore Stokes transport at blocked faces
([stokes_drift.md](stokes_drift.md)) — the onshore push is resolved in the
forcing but removed from the drift. That is an argument, not a validation: there
is no observational constraint on the rate here, hence the sweep range. `trap`
is degenerate because BSH's `H0 ≤ 0` tidal-flat flag
([h0_semantics.md](h0_semantics.md)) is no retentiveness proxy for the tide-free
Baltic, because moving `trap_flat` from 2.0 to 1.0 shifted the beached total by
0.9 points, and because Daily et al. and Onink et al. both report terrain
variation mattering little. `shore_type` is a diagnostic label, never a result.

## Store schema

`HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`, one per
`(regime, year, month, member)` — a grouped weight table, additive across
`release_doy`/month/year like the other `024x` stores, so `029` pools the
monthly partitions by summing.

| column | meaning |
|--------|---------|
| `release_hex` | release hex of the drifter (`024a` key space) |
| `release_doy` | release day-of-year of the originating zarr |
| `beach_hex` | hex where the weight stranded; `-1` = never-beached residual |
| `beach_age_bin` | `floor(deposit_age_days / age_bin_days)`; `-1` for residual |
| `shore_type` | `wall` / `flat` at the stranding site (`none` for residual); diagnostic only while `trap` is degenerate |
| `disp_bin` | `floor(disp_km / disp_bin_km)`; `-1` for residual. The sidecar saturates `disp` at 255 km, so the top bin pools everything beyond |
| `weight` | summed stranded weight (expected particles) in the group |

Deposits + residual per source hex sum to that hex's released drifter count, so
the beached fraction is `sum(weight | beach_hex ≥ 0) / sum(weight)`. The
`member` tag names the rate-model point and consumers take it as an opaque
string: `step_wc<w_c>_ts<tau_storm_hours>_tc<tau_calm_days|inf>[_d<delta>]`, or
`sat_t<tau0_hours>_wt<w_tau>` for the superseded ramp, `.` → `p` (e.g.
`step_wc0p1_ts3_tcinf`). `disp_bin_km` is a build↔consumer contract *not* stored
in the parquet — pass `029`/`031` the value `024e` built with, exactly as
`distance_bin_km` works in the distance store
([hexbinning_and_connectivity.md](hexbinning_and_connectivity.md)).

## Production setting

The reducer defaults: `w_c = 0.10` m/s, `τ_storm = 3 h`, `τ_calm = ∞`,
`δ = 0`, `trap ≡ 1`, `band_m = 2000`, `max_float_days = 60` — member tag
`step_wc0p1_ts3_tcinf`. `w_c = 0.10` is the p93 of onshore Stokes over
in-band hours, so 3.9 % of in-band hours count as storm hours. Pooled over
surface_stokes 2016–2019 (16.58 M real drifters):

| quantity | value |
|---|---|
| beached within 60 d | 75.9 % |
| drifters ever in band during a storm hour | 83.7 % |
| stranding hexes (weight > 1) | 1 736 |
| Gini of stranded weight over hexes | 0.62 |
| median stranding age / travel distance | 10–20 d / 40–50 km |

The two exposure numbers being close is the exposure limit at work: at
`τ_storm = 3 h` almost everything in band during a storm strands, so `w_c`
alone decides the total. The setting is a working choice, not a
calibration — there is no observational constraint on `w_c` in this study,
so results are reported as the range below, and `w_c` is the number to argue
about.

**Seasonality is the dominant signal**, interannual variability is not
(beached %, releases pooled over years / months):

| release month | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `w_c = 0.10` | 65 | 64 | 66 | 64 | 72 | 72 | 74 | 84 | 89 | 89 | 89 | 85 |
| `w_c = 0.15` | 39 | 31 | 26 | 22 | 25 | 23 | 23 | 40 | 51 | 52 | 54 | 54 |

| release year | 2016 | 2017 | 2018 | 2019 |
|---|---|---|---|---|
| `w_c = 0.10` | 74 | 79 | 74 | 77 |

Autumn releases meet the storm season inside their viability window; the
per-month partitions are the primary products, the pooled year a summary.

## Sensitivity

Sweep over the 48 partitions per member (`031`, `Figures/031/`); `τ_storm`
in hours, `τ_calm` in days, `ever` = share of drifters ever exposed to a
storm hour in band, Spearman = per-hex stranded weight against the
production member:

| member | beached % | ever % | Gini | Spearman |
|---|---|---|---|---|
| `w_c = 0.05` | 94.0 | 95.6 | 0.57 | 0.93 |
| `w_c = 0.075` | 88.6 | 92.5 | 0.58 | 0.98 |
| **`w_c = 0.10`** | **75.9** | **83.7** | **0.62** | 1 |
| `w_c = 0.15` | 36.7 | 47.3 | 0.70 | 0.89 |
| `w_c = 0.10`, `τ_calm = 365` | 78.1 | 83.7 | 0.69 | — |
| `w_c = 0.10`, `τ_storm = 1` | 81.6 | 83.7 | 0.59 | 0.99 |
| `w_c = 0.10`, `τ_storm = 12` | 60.2 | 83.7 | 0.66 | 0.98 |
| `w_c = 0.10`, `δ = 0.02` | 77.2 | 88.0 | 0.61 | 0.9997 |
| saturating `τ0 = 480 h, w_tau = 0.05` | 42.8 | — | 0.76 | — |

- **`w_c` is the axis.** It moves the total from 37 % to 94 % across
  p72–p99 of the forcing and is the only knob that reorders the map:
  adjacent values share 84 of the top-100 stranding hexes, `0.10` vs `0.15`
  share 59. Per-source beached fractions are more stable (Spearman ≥ 0.92
  between any two members).
- **`τ_storm` scales the total, not the pattern** (Spearman ≥ 0.98 from
  1 h to 12 h); below ~6 h it is irrelevant.
- **`τ_calm = 1 yr`** adds 2–7 points and ~180 low-weight hexes by draining
  calm-water residents; it changes where nothing else strands.
- **`δ`** is noise-level: keep the hard step.
- **The saturating form** (kept in 024e for reproduction) gives a lower total
  at a much older, farther stranding profile — its rate is highest at weak
  forcing, so it strands slowly everywhere instead of fast where waves hit.

**The beached fraction is not a reportable number**; the pattern and the
seasonal contrast are.

## Limitations

Real beaching is a surf/swash process below the BSH grid and coastline; where
Fucus sits on the coarse (~5 km) grid the band is finer than the physics, so
read the spatial *pattern*, not absolute rates or exact locations (stranding ≠
source proximity; López et al. 2017). Totals are highly parameter-sensitive — in
the Baltic the beaching scheme can dominate the answer (Siht et al. 2025) — so
the beached fraction is reported as a sweep range, not a value. Stranding is
terminal by design; beach-cast wrack in the tide-free Baltic is wind- and
water-level remobilised (Hammann & Zimmer 2014), and a reversible resuspension
timescale bounding that effect is deferred.

## Cross-references

- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the `024x`
  store pattern and key schema this reuses.
- [survival_occupancy.md](survival_occupancy.md) — same sidecar and rate model,
  aggregated as surviving occupancy.
- [h0_semantics.md](h0_semantics.md) — the land-sea mask behind the raster;
  [stokes_drift.md](stokes_drift.md) and
  [wam_extrapolation.md](wam_extrapolation.md) — the wave field and its fill.
- [visualisations.md](visualisations.md) — `029`/`031` plot rationale;
  [job_scripts.md](job_scripts.md) — how builder and reducers are scheduled.
