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
| Sidecar | [`024d_BuildBeachingForcing`](../notebooks/024d_BuildBeachingForcing.py) | trajectory zarrs + raw `baltic_highres` Stokes + BSH H0 statics + shoreline-class table + `024a` key | one forcing zarr per trajectory zarr |
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
| `ff` | uint8, percent | flat fraction of the nearest coastal land cell, in band (0 outside) | the substrate axis of the rate model; one byte, compresses like `dist` |
| `disp` | uint8, km | crow-flies displacement from release; 255 = saturated | slowly varying, so it compresses like `hex` |

Attrs record provenance and the sampling contract the reducers assert against
(`band_max_m`, `window_days`, `raster_dx_m`, `stokes_fill_max_cells`,
`hex_radius`, `shoreclass_sha256`, `release_time`/`release_doy`/`regime`/`release_year`,
source trajectory and obs counts, `builder`, `git_sha`); chunks are `(10000, nobs)`,
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
rate with two states — a forcing-independent background and a strong-wave
term. `w_c = 0.10` m/s is the p93 of in-band onshore Stokes — a windy day
with the right fetch, not a meteorological storm — and direction and fetch
set the onshore forcing as much as wind strength does, hence "strong wave"
rather than "storm":

```
1/τ = 1/τ_calm + trap · r(w_on) / τ_strong    in band, else 0
r(w) = clip((w − w_c + δ/2) / δ, 0, 1)        step at w_c, optional
                                               linear edge of width δ
```

Nothing in the data constrains the shape of `r` between calm and strong-wave,
so a shape parameter (Hill exponent, half-saturation) is pure nuisance; the
step is the limiting case of any such form and its one knob means what it
says. The surviving weight decays by `exp(−a)` each in-band step (`a = Δt/τ`), and the
weight deposited at step `h` is the telescoping survival difference
`dep = exp(−A_before) − exp(−A_after)`, `A` the running `Σ a` — the probability
that a first-stranding process strands at `h`, so the ensemble sum is the
stranding field as a smooth expectation. Weight still surviving at the end of
the window is the never-beached residual.

| parameter | meaning | how it is chosen |
|---|---|---|
| `w_c` | onshore Stokes above which stranding is on (m/s) | against the measured in-band forcing (p75 0.057, p90 0.088, p99 0.15 m/s); **the** sweep axis |
| `tau_strong_hours` | e-folding time while `w_on ≥ w_c` | hours; below ~6 h everything in band during a strong-wave event strands and the value stops mattering |
| `tau_calm_days` | background in-band rate regardless of forcing; 0 = off | ∞ or O(1 yr); a year loses 15 % over 60 d, so it is a real axis |
| `delta` | linear edge width around `w_c` | 0 (hard step) by default; `≈0.02` checks against nearest-hour / nearest-cell sampling noise |
| `trap_flat`, `trap_wall` | strong-wave rate factor on a fully flat / fully walled cell; the realised factor is `trap_wall + (trap_flat − trap_wall)·ff/100` | **1.0 / 0.0** at the production point: the rate is proportional to the substrate's flat fraction and pure-wall shore retains nothing. `1.0 / 1.0` keeps the term inert (the wave field alone sets the rate); every member with a weight ≠ 1 names both in its tag |
| `band_m` | near-shore band width (m), ≤ the sidecar's `band_max_m` | 2000; quantised by the mask — on the 5 km coarse grid 1–4 km all select the same first cell ring, so sweep it coarsely and report fine nest and coarse regions apart |
| `max_float_days` | viability window (d); slices the sidecar's `obs` axis | 60; a step-function `L(t)` |
| `age_bin_days`, `disp_bin_km` | store axis granularity | 10 d, 10 km |

Once `τ_strong` ≪ strong-wave event duration the model degenerates gracefully
into a **substrate-weighted exposure model**: the beached fraction is the share
of particles ever in band during a strong-wave hour inside the viability
window, each such encounter counted at its `trap` factor, and the stranding map
is coastal residence × strong-wave climatology × flat fraction. At the inert
weights the trap factor drops out and exposure alone sets the total; at the
production weights it is the binding term (see Production setting).
Strong-wave stranding of freshly released material is signal, not artefact:
Fucus releases are coastal, so a particle can be exposed from `t = 0`.
`disp_bin` records the crow-flies travel distance at stranding precisely so
that can be *shown* rather than filtered out — a diagnostic axis, never a
mask.

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
is no observational constraint on the rate here, hence the sweep range.

### Shore type

`ff` is the **flat fraction** of a coastal cell: the share of classified
shoreline inside it that is low and unconsolidated (sand, shingle, marsh,
mudflat) rather than a hard wall (rock, cliff, armoured bank). It comes from
`data/shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet` —
sub-segment midpoints in EPSG:3035, one `flat_fraction` and one `seg_len_m`
each (attribution in `ATTRIBUTION.md`). `024d` snaps the midpoints onto its
own 500 m EPSG:3035 raster (the table's CRS, so no reprojection), takes the
`seg_len_m`-weighted mean per cell, and lets cells the table does not reach
inherit from the nearest that does. A sub-segment with **no class
attribution counts 0.5** — halfway between the two trap weights, so
unattributed coast is neither pushed towards trapping nor away from it, and
no separate "unknown" flag has to be carried through the reducers. Per
position, `ff` is the value of the nearest land cell, as percent.

The rate model scales the strong-wave term by
`trap = trap_wall + (trap_flat − trap_wall)·ff/100`, linear in the fraction,
so a half-flat cell traps halfway between the two weights rather than being
forced into a class. **The production point is `trap_flat = 1`,
`trap_wall = 0`**, i.e. `trap = ff/100`: retention is proportional to how
much of the local shoreline is low and unconsolidated, and a stretch
classified as pure wall retains nothing at all. Substrate is therefore a
first-order term in the production physics, not a switched-off one — it
decides *whether* a coastal hex can strand, while the wave field decides
*when*. A member with a weight other than 1 names both in its tag
(`…_tf<flat>_tw<wall>`, `.` → `p`, so the production tag is
`step_wc0p1_ts3_tcinf_tf1_tw0`) and writes its own store files; the suffix is
absent only from the substrate-inert `1.0/1.0` member.

**Substrate and wave forcing are not the same gradient twice.** The worry
with two multiplicative terms is collinearity — exposed shore being the wall
shore, so `ff` and `w_on` measure one thing. Over 1.25 G in-band
particle-hours (5 sidecars, 2016–2019) they correlate weakly and *positively*:
Pearson 0.14, Spearman 0.13, with mean `ff` 61 % during strong-wave hours
against 38 % otherwise. Strong-wave exposure in this domain falls
preferentially on the flatter, more retentive coast, so the trap term
sharpens the wave signal rather than cancelling it.

`shore_type` in the store is a threshold label on the same fraction —
`flat` where `ff ≥ 50` at the stranding hour, else `wall` — and stays
diagnostic: it is a coarser cut than the rate, which varies continuously with
`ff`, so a `wall`-labelled hex at `ff = 40` still strands at 0.4× the
pure-flat rate. At the production weights 86.0 % of stranded weight lands on
`flat` hexes against 41.8 % at the inert weights, which is the substrate term
showing up in the deposition field rather than in the coastline's own
composition.

## Store schema

`HexAgg_beaching_r<radius>m_<regime>_<year>_mMM_<member>.parquet`, one per
`(regime, year, month, member)` — a grouped weight table, additive across
`release_doy`/month/year like the other `024x` stores, so `029` pools by
summing the monthly partitions whose month falls in its `season` parameter
(DJF/MAM/JJA/SON, or ALL for the whole year).

| column | meaning |
|--------|---------|
| `release_hex` | release hex of the drifter (`024a` key space) |
| `release_doy` | release day-of-year of the originating zarr |
| `beach_hex` | hex where the weight stranded; `-1` = never-beached residual |
| `beach_age_bin` | `floor(deposit_age_days / age_bin_days)`; `-1` for residual |
| `shore_type` | `flat` where the stranding site's `ff ≥ 50`, else `wall` (`none` for residual); diagnostic |
| `disp_bin` | `floor(disp_km / disp_bin_km)`; `-1` for residual. The sidecar saturates `disp` at 255 km, so the top bin pools everything beyond |
| `weight` | summed stranded weight (expected particles) in the group |

Deposits + residual per source hex sum to that hex's released drifter count, so
the beached fraction is `sum(weight | beach_hex ≥ 0) / sum(weight)`. The
`member` tag names the rate-model point and consumers take it as an opaque
string:
`step_wc<w_c>_ts<tau_strong_hours>_tc<tau_calm_days|inf>[_d<delta>][_tf<trap_flat>_tw<trap_wall>]`,
`.` → `p` (production: `step_wc0p1_ts3_tcinf_tf1_tw0`). `disp_bin_km` is a build↔consumer contract *not* stored
in the parquet — pass `029`/`031` the value `024e` built with, exactly as
`distance_bin_km` works in the distance store
([hexbinning_and_connectivity.md](hexbinning_and_connectivity.md)).

## Production setting

The reducer defaults: `w_c = 0.10` m/s, `τ_strong = 3 h`, `τ_calm = ∞`,
`δ = 0`, **`trap_flat = 1`, `trap_wall = 0`**, `band_m = 2000`,
`max_float_days = 60` — member tag `step_wc0p1_ts3_tcinf_tf1_tw0`.
`w_c = 0.10` is the p93 of onshore Stokes over in-band hours, so 3.9 % of
in-band hours are strong-wave hours; `trap = ff/100` then decides how much of
each such hour a given shore can convert into stranding. Pooled over
`surface_stokes` release years 2016–2025 (38.94 M real drifters), with the
substrate-inert member (`trap_flat = trap_wall = 1`, tag
`step_wc0p1_ts3_tcinf`) beside it for scale:

| quantity | production `…_tf1_tw0` | substrate inert |
|---|---|---|
| beached within 60 d | **41.6 %** | 76.9 % |
| hexes receiving stranded weight | 1 207 (1 179 with weight > 1) | 1 824 (1 781) |
| Gini of stranded weight over hexes | 0.68 | 0.63 |
| median stranding age (weight-weighted) | 7.0 d | 6.4 d |
| median travel distance at stranding | 46.6 km | 40.9 km |
| stranded weight on `flat` hexes (`ff ≥ 50`) | 86.0 % | 41.8 % |

The Gini is computed over hexes that *receive* weight — zero-weight hexes are
dropped, so it measures concentration among receiving hexes, not over the
whole coast; the hex count above it is what anchors that scale. 0 = every
receiving hex takes an equal share, 1 = one hex takes everything.

83.7 % of drifters are ever in band during a strong-wave hour (a property of
`w_c` and `band_m` alone, measured on the 2016–2019 sidecars; the trap weights
cannot change it). At the inert weights the beached total nearly saturates
that exposure share — 76.9 % against 83.7 %, the exposure limit at work, with
`w_c` alone setting the total. At the production point the two decouple:
41.6 % against the same 83.7 %, because half of the exposed coastline is wall
enough to retain little or nothing. **Substrate, not wave exposure, is what
caps the production beached fraction.**

Neither `w_c` nor the trap weights are calibrated — there is no observational
constraint on either in this study — so both are reported as ranges below.

Release season still orders the beached fraction, but it no longer dominates
it (beached %, releases pooled over years / months):

| release month | Jan | Feb | Mar | Apr | May | Jun | Jul | Aug | Sep | Oct | Nov | Dec |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| production `…_tf1_tw0` | 41 | 40 | 39 | 38 | 39 | 39 | 39 | 41 | 46 | 45 | 46 | 46 |
| substrate inert | 71 | 71 | 72 | 66 | 69 | 71 | 75 | 85 | 90 | 89 | 88 | 81 |

| release year | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 |
|---|---|---|---|---|---|---|---|---|---|---|
| production `…_tf1_tw0` | 43 | 42 | 41 | 42 | 42 | 45 | 39 | 40 | 39 | 46 |
| substrate inert | 74 | 79 | 74 | 77 | 85 | 79 | 72 | 76 | 75 | 80 |

Cumulative beached fraction by 60 d per release season, pooled 2016–2025,
with the per-release-year min–max in brackets:

| season | production `…_tf1_tw0` | substrate inert |
|---|---|---|
| DJF | 42 % [35–51] | 74 % [65–89] |
| MAM | 39 % [33–49] | 69 % [62–80] |
| JJA | 40 % [34–45] | 77 % [66–81] |
| SON | 46 % [39–53] | 89 % [83–93] |
| ALL | 42 % [39–46] | 77 % [72–85] |

This is the second substantive consequence of the trap term. At the inert
weights the seasonal spread is 20 points (MAM 69 % → SON 89 %) against ~5
points of interannual spread, which is why the seasonal split is the primary
product. At the production point the spread is 7 points (39 % → 46 %) and the
interannual min–max *within* a season is as wide or wider (SON 39–53). Autumn
still meets the strong-wave season inside the viability window, but the extra
exposure falls largely on shore the trap term does not let strand, so the
**total** is no longer a strongly seasonal quantity. The seasonal split stays
the primary product for *where* material strands — that is a map question, not
a scalar one — but the pooled-year beached fraction is no longer an average of
two distinct regimes.

## Sensitivity

Two independent axes, each swept over the 113 partitions per member (`031`,
`Figures/031/`), `surface_stokes`, release years 2016–2025. Every member on the
wave axis was built at the **inert** trap weights, so their Spearman column is
against the inert `w_c = 0.10` member, not against production; no wave member
has been rebuilt at the production trap point, and mixing the two axes in one
table would compare members that differ in two knobs at once.

### Wave axis, at `trap_flat = trap_wall = 1`

`τ_strong` in hours, `τ_calm` in days, Spearman = per-hex stranded weight
against the inert `w_c = 0.10` member over the union of receiving hexes. Gini
is over hexes that receive weight (see Production setting above); it rises with
`w_c` because a higher threshold keeps only the wave-exposed shores:

| member | beached % | Gini | hexes | Spearman |
|---|---|---|---|---|
| `w_c = 0.05` | 94.5 | 0.58 | 1 828 | 0.94 |
| `w_c = 0.075` | 89.4 | 0.60 | 1 829 | 0.98 |
| **`w_c = 0.10`** (reference) | **76.9** | **0.63** | **1 824** | 1 |
| `w_c = 0.15` | 37.4 | 0.71 | 1 780 | 0.92 |
| `w_c = 0.10`, `τ_calm = 365` | 79.1 | 0.70 | 2 331 | 0.99 |
| `w_c = 0.10`, `τ_strong = 1` | 82.5 | 0.61 | 1 819 | 0.99 |
| `w_c = 0.10`, `τ_strong = 12` | 61.1 | 0.67 | 1 826 | 0.98 |
| `w_c = 0.10`, `δ = 0.02` | 78.2 | 0.63 | 1 826 | 0.9998 |

- **`w_c` is the axis on this side.** It moves the total from 37 % to 94 %
  across p72–p99 of the forcing and is the only wave knob that reorders the
  map: adjacent values share 84–88 of the top-100 stranding hexes, `0.10` vs
  `0.15` share 55.
- **`τ_strong` scales the total, not the pattern** (Spearman ≥ 0.98 from
  1 h to 12 h); below ~6 h it is irrelevant.
- **`τ_calm = 1 yr`** adds ~2 points and ~500 low-weight hexes (2 331 against
  1 824) by draining calm-water residents; it changes where nothing else
  strands.
- **`δ`** is noise-level: keep the hard step.

### Trap (substrate) axis, at `w_c = 0.10`, `τ_strong = 3 h`

`trap_wall` swept at `trap_flat = 1` — pure-flat shore keeps the full
strong-wave rate, pure-wall shore traps at a fraction of it. Spearman and the
top-100 overlap are per `beach_hex` stranded weight against the **production**
member, over the union of receiving hexes (n = 1 824):

| member | beached, ALL | beached, SON | hexes | Spearman | top-100 shared |
|---|---|---|---|---|---|
| `step_wc0p1_ts3_tcinf` (inert) | 0.769 | 0.891 | 1 824 | 0.32 | 40 |
| `…_tf1_tw0p5` | 0.725 | 0.851 | 1 824 | 0.37 | 45 |
| `…_tf1_tw0p25` | 0.668 | 0.793 | 1 824 | 0.45 | 56 |
| **`…_tf1_tw0`** (production) | **0.416** | **0.456** | **1 207** | 1 | 100 |

Suppressing the wall rate to a quarter rescales the total (77 % → 67 %) and
leaves the map almost alone — against the inert member those two intermediates
score Spearman 0.996 and 0.981. The last step, to `trap_wall = 0`, is where the
model changes character rather than degree: 617 of the 1 824 receiving hexes
stop receiving anything, the total falls 46 % relative, the stranded weight
moves from 42 % to 86 % `flat`, and the per-hex ranking against the inert
member drops to Spearman 0.32 with 40 of its top-100 hexes shared. **`trap_wall`
is not a smooth knob near zero**: `tw0` is a different model, not a
perturbation of the inert baseline, and the intermediates
(`…_tf1_tw0p25`, `…_tf1_tw0p5`) are what bracket the space between them.

**The beached fraction is not a reportable number**; the pattern and the
member range are.

## Limitations

**`trap_wall = 0` is an end-member assumption, not a measured one.** The
production physics says a shoreline classified as entirely rock, cliff or
armoured bank *never* retains a Fucus propagule. Nothing in the HELCOM BRISK /
CLMS classification, or in any observation available here, says that: it says
which substrate is there, not what its retention rate is. The choice is
defensible as one end of the plausible range — hard reflective shore in a
tide-free basin gives drifting thalli little to lodge in — but it is the
strongest statement the parameterisation can make, and it costs 46 % of the
beached total and 617 of 1 824 stranding hexes relative to the substrate-inert
member (see Trap sensitivity above, where the per-hex Spearman between the two
is 0.32). The intermediate members `…_tf1_tw0p25` and `…_tf1_tw0p5` are built
for every release year and bracket the cases where wall shore retains at a
reduced but non-zero rate: they put the beached fraction at 67 % and 72 %
against production's 42 % and the inert member's 77 %, while their stranding
map stays within Spearman 0.98 of the inert one. Anything that depends on the
*absolute* stranded amount, or on stranding at a specific wall-dominated
stretch, must be read across that set of four rather than off production
alone.

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
