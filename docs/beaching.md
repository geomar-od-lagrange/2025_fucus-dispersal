# Beaching diagnostic

A post-simulation estimate of where and when drifting *Fucus* propagules strand
on the coast. Each real drifter carries unit surviving (free-drifting) weight;
inside a near-shore band a **fraction** of that weight strands at the current
coastal hex each step and leaves the drifting pool, the rest drifting on. This
**weighted / fractional** deposition is the deterministic expectation of a
per-step first-stranding process: releases carry only ~100 particles per seeding
cell, too few for a random-strand rule to resolve the high-age tail where the
free-drifting/beached split is decided. **No lifetime is applied**:
`max_age_days` is a reduction window over the sidecar's `obs` axis, and the
store stays resolved per `beach_age_bin`, so a biological survival `L(t)` can
be multiplied in downstream without a rebuild. It runs after the physics as a
re-runnable reduction, the "aggregate once, consume cheaply" split the
[hex store](hexbinning_and_connectivity.md) uses, so the rate parameters sweep
without touching the runs.

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Sidecar | [`024d_BuildBeachingForcing`](../notebooks/024d_BuildBeachingForcing.py) | trajectory zarrs + raw `baltic_highres` Stokes + BSH H0 statics + shoreline-class table + `024a` key | one forcing zarr per trajectory zarr |
| Reduce | [`024e_BuildBeaching`](../notebooks/024e_BuildBeaching.py) | sidecar + key | `HexAgg_beaching_*.parquet` |
| Consume | [`029_BeachingMaps`](../notebooks/029_BeachingMaps.py) | beaching parquet + key | PNGs under `Figures/029/`; per-horizon GeoJSON and `afloat_by_origin_*.csv` under `Exports/029/` |
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
width, the reduction window, the trap weights and the binning are reductions
over them — so a member costs ~2.5 s per zarr against the ~90 s a re-read of the
trajectories and the raw wave field costs. ~110 MB per zarr at
`window_days = 220` (~74 GB for `surface_stokes` 2016–2025), ~94 s per zarr to build, single-process numpy, no
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
| `tau_calm_days` | background in-band rate regardless of forcing; 0 = off | ∞ or O(1 yr); a year drains 39 % over the 180 d window, so it is a real axis |
| `delta` | linear edge width around `w_c` | 0 (hard step) by default; `≈0.02` checks against nearest-hour / nearest-cell sampling noise |
| `trap_flat`, `trap_wall` | strong-wave rate factor on a fully flat / fully walled cell; the realised factor is `trap_wall + (trap_flat − trap_wall)·ff/100` | **1.0 / 0.0** at the production point: the rate is proportional to the substrate's flat fraction and pure-wall shore retains nothing. `1.0 / 1.0` keeps the term inert (the wave field alone sets the rate); every member with a weight ≠ 1 names both in its tag |
| `band_m` | near-shore band width (m), ≤ the sidecar's `band_max_m` | 2000; quantised by the mask — on the 5 km coarse grid 1–4 km all select the same first cell ring, so sweep it coarsely and report fine nest and coarse regions apart |
| `max_age_days` | reduction window (d): the half-open `[0, max_age_days)` hours of the sidecar's `obs` axis; one value for 024e/024f/024g | 180, ≤ the sidecar's `window_days` (220); not a lifetime — weight still drifting at its end is the residual |
| `age_bin_days`, `disp_bin_km` | store axis granularity | 10 d, 10 km |

Once `τ_strong` ≪ strong-wave event duration the model degenerates gracefully
into a **substrate-weighted exposure model**: the beached fraction is the share
of particles ever in band during a strong-wave hour inside the reduction
window, each such encounter counted at its `trap` factor, and the stranding map
is coastal residence × strong-wave climatology × flat fraction. At the inert
weights the trap factor drops out and exposure alone sets the total; at the
production weights it is the binding term (see Production setting).
Strong-wave stranding of freshly released material is signal, not artefact:
Fucus releases are coastal, so a particle can be exposed from `t = 0`.
`disp_bin` records the crow-flies displacement at stranding precisely so
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
`.` → `p` (production: `step_wc0p1_ts3_tcinf_tf1_tw0`).

Each partition carries its build parameters as `df.attrs`, which pandas writes
into the parquet metadata and `pd.read_parquet` restores: `builder`, `regime`,
`release_year`, `release_month`, `hex_radius`, `member`, `max_age_days`,
`age_bin_days`, `disp_bin_km`, `output_dt_hours`, the rate-model parameters
(`w_c`, `tau_strong_hours`, `tau_calm_days`, `delta`, `trap_flat`,
`trap_wall`, `band_m`), `git_sha` and `slurm_job_id` (empty outside SLURM).
024f and 024g write the same set minus `disp_bin_km`. Every consumer
(`028`–`033`, and `025`/`026`/`026a`/`026b` for a 024f member) reads the attrs
per file **before** concatenating — pandas drops attrs that differ between
the frames it concatenates — and asserts that all pooled partitions agree on
`max_age_days`, `age_bin_days`, `disp_bin_km` (024e) and the rate-model keys,
that `age_bin_days`/`disp_bin_km` equal its own parameters, and that the
window covers its largest horizon. `builder`, `git_sha`, `slurm_job_id`,
`release_year` and `release_month` legitimately differ and are never
compared; a partition with no attrs (written before the window was recorded)
raises, asking for a rebuild. Across members (`029`'s contrast member, `031`,
`032`, `033`) only the window and binning must agree. Considered encoding the
window in the filename; rejected because consumers' member regexes would
absorb it.

## Production setting

The reducer defaults: `w_c = 0.10` m/s, `τ_strong = 3 h`, `τ_calm = ∞`,
`δ = 0`, **`trap_flat = 1`, `trap_wall = 0`**, `band_m = 2000`,
`max_age_days = 180` — member tag `step_wc0p1_ts3_tcinf_tf1_tw0`.
`w_c = 0.10` is the p93 of onshore Stokes over in-band hours, so 3.9 % of
in-band hours are strong-wave hours; `trap = ff/100` then decides how much of
each such hour a given shore can convert into stranding. Pooled over
`surface_stokes` release years 2016–2024 (37.30 M released water-seeded
drifters, Σweight) over the 180 d window, with the substrate-inert member
(`trap_flat = trap_wall = 1`, tag `step_wc0p1_ts3_tcinf`) beside it for
scale — the first five rows from `031` (season `ALL`), the last from `033`
(season `ALL`):

| quantity | production `…_tf1_tw0` | substrate inert |
|---|---|---|
| stranded within 180 d | **54.3 %** | 93.0 % |
| hexes receiving stranded weight | 1 221 | 1 841 |
| Gini of stranded weight over hexes | 0.695 | 0.639 |
| median stranding age (weight-weighted, 10 d bins) | 14.9 d | 11.4 d |
| median displacement at stranding | 71.1 km | 51.5 km |
| stranded weight with `shore_type = flat` (stranding-site `ff ≥ 50`) | 82.8 % | 39.2 % |

The Gini is computed over hexes that *receive* weight — zero-weight hexes are
dropped, so it measures concentration among receiving hexes, not over the
whole coast; the hex count above it is what anchors that scale. 0 = every
receiving hex takes an equal share, 1 = one hex takes everything.

83.7 % of drifters are ever in band during a strong-wave hour within their
first 60 d (a property of `w_c` and `band_m` alone, measured on the 2016–2019
sidecars; the trap weights cannot change it). Over the 180 d window the inert
member strands 93.0 %, so that 60 d exposure share is no ceiling here. At the
same `w_c` the production total is 54.3 %: the trap axis below takes the total
from 93.0 % (`trap_wall = 1`) to 54.3 % (`trap_wall = 0`) with the wave
forcing unchanged. **Substrate, not wave exposure, is what caps the
production beached fraction.**

Neither `w_c` nor the trap weights are calibrated — there is no observational
constraint on either in this study — so both are reported as ranges below.

Stranded within 180 d per release season, pooled over 2016–2024 with the
per-release-year min–max in brackets (production: `032` table A; inert: the
`031` run of each season):

| season | production `…_tf1_tw0` | substrate inert |
|---|---|---|
| DJF | 50.9 % [45.2–60.0] | 88.8 % [83.9–94.7] |
| MAM | 57.3 % [48.6–68.4] | 93.0 % [90.0–95.7] |
| JJA | 55.6 % [49.1–63.1] | 95.9 % [93.4–97.5] |
| SON | 53.1 % [43.5–64.5] | 94.4 % [90.3–97.4] |
| ALL | 54.3 % [49.3–59.1] | 93.0 % [90.1–95.8] |

Over 180 d the release season barely orders the total: the seasonal spread is
6.4 points at production (DJF → MAM) and 7.1 points at the inert weights
(DJF → JJA), while the interannual min–max within one season reaches 21.0
points (production SON). The seasonal split stays the primary product for
*where* material strands — a map question, not a scalar one. The production
fraction per release year and season is `032`'s table A
(`Headline/HEADLINE.md`); the cumulative curve against age, with its
interannual band, is `029`'s `BeachedFractionCurve_*` (production `ALL`:
32.4 % by 30 d, 41.4 % by 60 d, 46.4 % by 90 d, 54.3 % by 180 d). No consumer
tabulates the fraction per release month.

## Sensitivity

Two independent axes, each swept over the 113 partitions per member (`031`,
`Figures/031/`), `surface_stokes`, release years 2016–2024, 180 d window,
season `ALL`. Every member on the wave axis was built at the **inert** trap
weights, so it compares against the inert `w_c = 0.10` member, not against
production; no wave member has been rebuilt at the production trap point, and
mixing the two axes in one table would compare members that differ in two
knobs at once. `031` reports totals and per-member summaries, not per-hex rank
correlations between members; the pattern is read off its
`BeachingSweepMaps_*` panels.

### Wave axis, at `trap_flat = trap_wall = 1`

`τ_strong` in hours, `τ_calm` in days, median age and displacement
weight-weighted over the stranded weight. Gini is over hexes that receive
weight (see Production setting above); it rises with `w_c` because a higher
threshold keeps only the wave-exposed shores:

| member | stranded % | Gini | hexes | median age (d) | median displacement (km) |
|---|---|---|---|---|---|
| `w_c = 0.05` | 99.4 | 0.587 | 1 840 | 5.0 | 25.2 |
| `w_c = 0.075` | 98.3 | 0.606 | 1 844 | 5.0 | 35.8 |
| **`w_c = 0.10`** (reference) | **93.0** | **0.639** | **1 841** | **11.4** | **51.5** |
| `w_c = 0.15` | 57.8 | 0.720 | 1 804 | 30.6 | 79.3 |
| `w_c = 0.05`, `τ_calm = 365` | 99.6 | 0.680 | 2 386 | 5.0 | 25.0 |
| `w_c = 0.075`, `τ_calm = 365` | 98.8 | 0.694 | 2 397 | 5.0 | 35.2 |
| `w_c = 0.10`, `τ_calm = 365` | 95.3 | 0.718 | 2 397 | 11.5 | 50.1 |
| `w_c = 0.15`, `τ_calm = 365` | 71.7 | 0.764 | 2 399 | 36.2 | 75.8 |
| `w_c = 0.10`, `τ_strong = 1` | 95.4 | 0.615 | 1 832 | 7.7 | 43.0 |
| `w_c = 0.10`, `τ_strong = 12` | 84.0 | 0.690 | 1 842 | 21.8 | 68.4 |
| `w_c = 0.10`, `δ = 0.02` | 93.7 | 0.637 | 1 844 | 10.9 | 50.3 |

- **`w_c` is the axis on this side.** Across p72–p99 of the forcing it moves
  the total from 57.8 % to 99.4 % and the median stranding age from 31 d to
  5 d.
- **`τ_strong` scales the total** (95.4 % at 1 h, 93.0 % at 3 h, 84.0 % at
  12 h); below ~6 h everything in band during a strong-wave event strands.
- **`τ_calm = 1 yr`** adds 0.2–13.9 points (most at `w_c = 0.15`, where the
  waves strand least) and ~550 receiving hexes (2 386–2 399 against
  1 804–1 844) by draining calm-water residents.
- **`δ`** is noise-level (+0.7 points): keep the hard step.

### Trap (substrate) axis, at `w_c = 0.10`, `τ_strong = 3 h`

`trap_wall` swept at `trap_flat = 1` — pure-flat shore keeps the full
strong-wave rate, pure-wall shore traps at a fraction of it (stranded %, `032`
table B; hexes, `031` `ALL`):

| member | stranded, ALL | stranded, SON | hexes |
|---|---|---|---|
| `step_wc0p1_ts3_tcinf` (inert) | 93.0 | 94.4 | 1 841 |
| `…_tf1_tw0p5` | 90.6 | 92.0 | 1 842 |
| `…_tf1_tw0p25` | 86.9 | 88.2 | 1 842 |
| **`…_tf1_tw0`** (production) | **54.3** | **53.1** | **1 221** |

Suppressing the wall rate to a quarter rescales the total (93.0 % → 86.9 %)
and keeps the receiving set. The last step, to `trap_wall = 0`, is where the
model changes character rather than degree: the receiving set shrinks from
1 841 to 1 221 hexes, the total falls to 54.3 %, and the stranded weight
moves from 39.2 % to 82.8 % `flat` (`033`). **`trap_wall` is not a smooth
knob near zero**: `tw0` is a different model, not a perturbation of the inert
baseline, and the intermediates (`…_tf1_tw0p25`, `…_tf1_tw0p5`) are what
bracket the space between them.

**The beached fraction is not a reportable number**; the pattern and the
member range (54.3–99.6 % over the 14 members, `ALL`) are.

## Limitations

**`trap_wall = 0` is an end-member assumption, not a measured one.** The
production physics says a shoreline classified as entirely rock, cliff or
armoured bank *never* retains a Fucus propagule. Nothing in the HELCOM BRISK /
CLMS classification, or in any observation available here, says that: it says
which substrate is there, not what its retention rate is. The choice is
defensible as one end of the plausible range — hard reflective shore in a
tide-free basin gives drifting thalli little to lodge in — but it is the
strongest statement the parameterisation can make: relative to the
substrate-inert member it takes the beached total from 93.0 % to 54.3 % and
the stranding hexes from 1 841 to 1 221 (see Trap sensitivity above). The
intermediate members `…_tf1_tw0p25` and `…_tf1_tw0p5` are built for every
release year and bracket the cases where wall shore retains at a reduced but
non-zero rate: they put the beached fraction at 86.9 % and 90.6 % against
production's 54.3 % and the inert member's 93.0 %, on the inert member's
receiving set (1 842 hexes). Anything that depends on the
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
