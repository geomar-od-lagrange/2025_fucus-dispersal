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
without touching the runs. All three regimes — `surface_stokes`, `bottom`,
`bottom_ub49` — strand under the **same criterion**: the same rate model, the
same thresholded quantity (the surface onshore Stokes `w_on` at the drifter
position) and the same absolute `w_c`, inside the same 750 m band (see
[The band](#the-band)).

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Sidecar | [`024d_BuildBeachingForcing`](../notebooks/024d_BuildBeachingForcing.py) | trajectory zarrs + raw `baltic_highres` Stokes + BSH H0 statics + shoreline-class table + `024a` key | one forcing zarr per trajectory zarr |
| Reduce | [`024e_BuildBeaching`](../notebooks/024e_BuildBeaching.py) | sidecar + key | `HexAgg_beaching_*.parquet` |
| Consume | [`029_BeachingMaps`](../notebooks/029_BeachingMaps.py) | beaching parquet + key | PNGs under `Figures/029/`; per-horizon GeoJSON, `afloat_by_origin_*.csv` and `stranded_curve_*.csv` (read by `035`) under `Exports/029/` |
| Sweep | [`031_BeachingSweep`](../notebooks/031_BeachingSweep.py) | several members + key | PNGs under `Figures/031/` |

[`024f_BuildSurvivalOccupancy`](../notebooks/024f_BuildSurvivalOccupancy.py) and
`030` read the same sidecar with the same rate model, recording where the
*surviving* weight is rather than where it leaves
([survival_occupancy.md](survival_occupancy.md)).

## The forcing sidecar

`output_root/BeachingForcing/<regime>/<year>/<trajectory zarr stem>.zarr`, one
per trajectory zarr of every regime,
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
`window_days = 220` (~74 GB for `surface_stokes` 2016–2025, about the same
per bottom regime), ~94 s per zarr to build, single-process numpy, no
Dask. The sidecar is regime-agnostic: it samples the raw wave field at
whatever positions the trajectory zarr holds. Considered and rejected: a sparse (COO) layout — 73–81 % of particle-hours
are in band, dense wins; storing `s(w)` or `Σ s` instead of `w_on` — freezes the
functional form, which is what the sidecar exists to free; zarr's default
Blosc-zstd-bitshuffle — 2–3× larger than raw Zstd on every array here.

`band_max_m = 5000` is the widest band the coarse grid resolves (81 % of
particle-hours, 73 % within 2 km); beyond it `w_on` is exactly 0 and `dist`
reads 255. **Land-seeded particles are dropped** (010's `land_seed`: release
position in a land cell of the 005 land set), as `024`/`024b`. `w_on` comes from the WAM field geodesically extrapolated to
full BSH coverage, without which 78.7 % of in-band samples would read zero
forcing ([wam_extrapolation.md](wam_extrapolation.md)).

## The rate model

Inside the near-shore band the beaching hazard is a Δt-invariant e-folding
rate with two states — a forcing-independent background and a strong-wave
term. `w_c = 0.10` m/s is an **absolute** threshold on the onshore Stokes,
the same for every regime and band — a windy day with the right fetch, not a
meteorological storm (8 % of in-band `surface_stokes` hours at 750 m reach
it); direction and fetch set the onshore forcing as much as wind strength
does, hence "strong wave" rather than "storm":

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
| `w_c` | onshore Stokes above which stranding is on (m/s) | absolute, 0.10; **the** sweep axis (0.05–0.15) |
| `tau_strong_hours` | e-folding time while `w_on ≥ w_c` | hours; below ~6 h everything in band during a strong-wave event strands and the value stops mattering |
| `tau_calm_days` | background in-band rate regardless of forcing; 0 = off | ∞ or O(1 yr); a year would drain 39 % of a drifter held in band for the whole 180 d window, but in the 750 m band it moves the pooled total ≤ 1.1 points (see Sensitivity) |
| `delta` | linear edge width around `w_c` | 0 (hard step) by default; `≈0.02` checks against nearest-hour / nearest-cell sampling noise |
| `trap_flat`, `trap_wall` | strong-wave rate factor on a fully flat / fully walled cell; the realised factor is `trap_wall + (trap_flat − trap_wall)·ff/100` | **1.0 / 0.0** at the production point: the rate is proportional to the substrate's flat fraction and pure-wall shore retains nothing. `1.0 / 1.0` keeps the term inert (the wave field alone sets the rate); every member with a weight ≠ 1 names both in its tag |
| `band_m` | near-shore band width (m), ≤ the sidecar's `band_max_m` | **750** = the first 500 m raster ring (see [The band](#the-band)); 2000 is the one sensitivity member; any value but 2000 tags `_b<band_m>` |
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

### The band

`band_m = 750` keeps the drifters within the first ring of 500 m raster
cells of the coast. `dist` is the EDT of the 500 m raster (`raster_dx_m`) in
whole cells and the reducers gate `dist·25 < band_m`, so 501–700 m selects
the edge-adjacent ring, 701–1000 m adds its diagonals (707 m), and 500 m
selects nothing; "the 500 m band" is therefore `band_m = 750`, water within
~0–500 m of the raster land edge, itself within ±250 m of the model coast.

Why the first ring, for every regime: within it the water is shallow against
the wave length, so the orbital motion reaches the bed and the surface/bed
distinction stops mattering — the same `w_on` criterion is defensible for
`bottom` drifters, which a 2 km band (several cells out on the fine nest)
would not be. For `surface_stokes` the narrower band costs little: under
strong onshore forcing surface drifters reach the wall anyway (ε rest point
0.01–0.07 cell, [wall_velocity.md](wall_velocity.md)). Measured on `surface_stokes`, release years 2016–2024 pooled (37.46 M
drifters, 180 d, production trap point; in-band hours and the lost-hex share
from the [Step 0 analysis](../plans/done/bottom_beaching.md)):

| quantity (`surface_stokes`) | 2000 m | 750 m |
|---|---|---|
| stranded ≤ 30 / 90 / 180 d | 38.4 / 57.7 / 70.3 % | 35.4 / 54.4 / 67.3 % |
| ≤ 180 d DJF / MAM / JJA / SON | 67.2 / 73.0 / 72.8 / 68.0 % | 63.9 / 70.1 / 69.9 / 65.0 % |
| receiving hexes | 1 206 | 1 097 |
| in-band particle-hours (rel.) | 1 | 0.67 |
| strong-wave in-band hours (rel.); share of in-band | 1; 6.05 % | 0.88; 7.93 % |

Sources: the 30 / 90 / 180 d row is `Exports/029/stranded_curve_surface_stokes_ALL_*`
(the 2000 m file has no `_b750` suffix); seasons are `HEADLINE.md` table A
(750 m) and table B (2000 m member `step_wc0p1_ts3_tcinf_tf1_tw0`); hexes
count the features of `029_surface_stokes_ALL_*_T180d.geojson`.

The band cuts in-band hours by a third but keeps 88 % of the strong-wave
in-band hours, so surface stranding moves by −3 points at 180 d, the same
shift in every season and well inside the interannual range; the 109 hexes
lost held 0.6 % of the 2000 m stranded weight (Spearman 0.955 of per-hex
share). The bed regimes are more band-sensitive: pooled ALL at 180 d,
`bottom` 37.2 → 22.6 % and `bottom_ub49` 21.9 → 8.0 % from 2000 m to 750 m
(table B): the 750–2000 m ring holds a larger share of their strong-wave
exposure than of the surface drifters'. A finer raster (250 m)
would only sharpen the fine nest — the coarse model coast is 5 km — and is a
full 024d rebuild; not done. The 2000 m member stays built for every reducer
and regime as the band sensitivity.

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
forcing but removed from the drift — and the bottom runs carry no Stokes in
their drift at all. That is an argument, not a validation: there
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
`step_wc0p1_ts3_tcinf_tf1_tw0_b750`) and writes its own store files; the suffix is
absent only from the substrate-inert `1.0/1.0` member.

**Substrate and wave forcing are not the same gradient twice.** The worry
with two multiplicative terms is collinearity — exposed shore being the wall
shore, so `ff` and `w_on` measure one thing. Over 1.25 G in-band
particle-hours (2 km band, 5 sidecars, 2016–2019) they correlate weakly and *positively*:
Pearson 0.14, Spearman 0.13, with mean `ff` 61 % during strong-wave hours
against 38 % otherwise. Strong-wave exposure in this domain falls
preferentially on the flatter, more retentive coast, so the trap term
sharpens the wave signal rather than cancelling it.

`shore_type` in the store is a threshold label on the same fraction —
`flat` where `ff ≥ 50` at the stranding hour, else `wall` — and stays
diagnostic: it is a coarser cut than the rate, which varies continuously with
`ff`, so a `wall`-labelled hex at `ff = 40` still strands at 0.4× the
pure-flat rate. The `flat` share of stranded weight at the production
weights against the inert ones (`033`) is the substrate term showing up in
the deposition field rather than in the coastline's own composition.

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
`step_wc<w_c>_ts<tau_strong_hours>_tc<tau_calm_days|inf>[_d<delta>][_tf<trap_flat>_tw<trap_wall>][_b<band_m>]`,
`.` → `p` (production: `step_wc0p1_ts3_tcinf_tf1_tw0_b750`). `_b<band_m>` appears
only off the 2000 m band — anchored to that value so the 2000 m sensitivity
member keeps its name — so the two bands never share a store file.

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
`δ = 0`, **`trap_flat = 1`, `trap_wall = 0`**, **`band_m = 750`**,
`max_age_days = 180` — member tag `step_wc0p1_ts3_tcinf_tf1_tw0_b750`, the
same for `surface_stokes`, `bottom` and `bottom_ub49`. `trap = ff/100`
decides how much of each strong-wave hour a given shore can convert into
stranding.

### Results

Production member, release years 2016–2024 pooled (9 years, 180 d window,
water-seeded drifters), all three regimes. Each row's source is in the last
column; `HEADLINE.md` is `<outputs>/Headline/HEADLINE.md` ([`032`](../notebooks/032_HeadlineReport.py)),
`029` means `<outputs>/Exports/029/`.

| quantity | `surface_stokes` | `bottom` | `bottom_ub49` | source |
|---|---|---|---|---|
| stranded ≤ 30 / 90 / 180 d (%) | 35.4 / 54.4 / 67.3 | 6.7 / 14.6 / 22.6 | 2.7 / 5.2 / 8.0 | `029/stranded_curve_<regime>_ALL_*_b750.csv` |
| per-year min–max at 180 d (%) | 62.0–73.5 | 19.6–25.4 | 7.0–9.3 | same; HEADLINE table A |
| ≤ 180 d, DJF / MAM / JJA / SON releases (%) | 63.9 / 70.1 / 69.9 / 65.0 | 21.6 / 17.5 / 24.1 / 27.6 | 7.9 / 6.0 / 8.1 / 10.0 | HEADLINE table A |
| spread across seasons (points) | 6.2 (DJF → MAM) | 10.0 (MAM → SON) | 4.0 (MAM → SON) | HEADLINE table A |
| receiving hexes (of 43 601) | 1 097 | 876 | 689 | features of `029/029_<regime>_ALL_*_b750_T180d.geojson` |
| Gini over receiving hexes; top-10 share (%) | 0.687; 13.4 | 0.749; 13.0 | 0.739; 16.3 | same geojsons, `stranded_weight` |
| particle-time afloat ≤ 180 d, `Σw/Σn` | 0.499 | 0.863 | 0.950 | HEADLINE table C |
| not stranded at 180 d, basin-wide (%) [year range] | 32.7 [26.5–38.0] | 77.4 [74.6–80.4] | 92.0 [90.7–93.0] | HEADLINE table G (`029/afloat_by_origin_*`) |
| not stranded at 180 d by origin subbasin (%), lowest – highest | 0.2 (Gdansk Basin) – 80.1 (The Quark) | 28.0 (Great Belt) – 99.9 (Åland Sea) | 64.1 (Great Belt) – 100.0 (Åland Sea) | HEADLINE table G |

"Not stranded" includes drifters that left the BSH domain and applies no
lifetime. The release season orders the surface total only weakly (6.2
points) against the interannual range within one season (up to 23.1 points,
SON); the bed regimes order the seasons more clearly, SON highest and MAM
lowest. The cumulative curve with its interannual band is panel (a) of `035`;
median stranding age and displacement per member are `031`'s stats panels
(`Figures/031/`), the `flat` share is `033`.

The Gini is computed over hexes that *receive* weight — zero-weight hexes are
dropped, so it measures concentration among receiving hexes, not over the
whole coast; the hex count beside it anchors that scale. 0 = every
receiving hex takes an equal share, 1 = one hex takes everything.

Neither `w_c` nor the trap weights are calibrated — there is no observational
constraint on either in this study — so both are reported as ranges (below).
No consumer tabulates the fraction per release month.

## Sensitivity

Every member is built on the production band (750 m) for every regime, by
`024efg_ReduceAll_job.sh` ([job_scripts.md](job_scripts.md)); `031` sweeps
the members of one regime (`Figures/031/`) and `032` table B tabulates the
pooled and per-season beached fraction of **every** `024e` member on disk
(15 per regime: the sweep below plus the 2000 m member), sorted by the `ALL`
column. Three independent axes, plus the band:

| axis | members (`step_…`) | at |
|---|---|---|
| wave | `w_c` 0.05 / 0.075 / 0.10 / 0.15, each with `τ_calm` ∞ and 365 d; `τ_strong` 1 / 12 h; `δ = 0.02` | inert trap weights (`trap_flat = trap_wall = 1`) |
| trap | `trap_wall` 0 (production) / 0.25 / 0.5 / 1 (inert) | `w_c = 0.10`, `τ_strong = 3 h`, `trap_flat = 1` |
| band | `band_m` 2000 (tag without `_b`) | the production point |

Every wave member sits at the **inert** trap weights, so it compares against
the inert `w_c = 0.10` member, not production; mixing the axes in one table
would compare members that differ in two knobs at once. 024f carries the
trap axis plus `w_c` 0.075 / 0.15 and `τ_calm` 365 d; 024g the trap axis. The
2000 m production member is the only band sensitivity; the full 2000 m sweep
would add a store per member and regime for no new axis.

Stranded within 180 d, pooled ALL seasons and release years 2016–2024 (%;
`HEADLINE.md` table B, `_b750` unless noted; the full table adds the
per-season columns and the interannual range):

| member (`step_…`) | `surface_stokes` | `bottom` | `bottom_ub49` |
|---|---|---|---|
| **production** `wc0p1_ts3_tcinf_tf1_tw0_b750` | **67.3** | **22.6** | **8.0** |
| wave: `wc0p05_ts3_tcinf` / `wc0p075_ts3_tcinf` | 99.6 / 99.0 | 56.5 / 46.8 | 16.8 / 13.5 |
| wave: `wc0p1_ts3_tcinf` (inert reference) | 96.2 | 34.7 | 9.8 |
| wave: `wc0p15_ts3_tcinf` | 70.7 | 12.1 | 3.5 |
| `τ_strong` 1 h / 12 h: `wc0p1_ts1_tcinf` / `wc0p1_ts12_tcinf` | 97.6 / 90.3 | 40.9 / 21.4 | 11.9 / 5.9 |
| `τ_calm` 365 d: `wc0p1_ts3_tc365` | 96.6 | 35.8 | 10.0 |
| `δ = 0.02`: `wc0p1_ts3_tcinf_d0p02` | 96.6 | 35.5 | 10.0 |
| trap: `…_tf1_tw0p5` / `…_tf1_tw0p25` | 94.6 / 91.9 | 31.9 / 29.0 | 9.3 / 8.9 |
| band 2000 m: `wc0p1_ts3_tcinf_tf1_tw0` (no `_b`) | 70.3 | 37.2 | 21.9 |

Across all 15 members the pooled ALL value spans 67.3–99.7 % on
`surface_stokes`, 12.1–56.7 % on `bottom` and 3.5–21.9 % on `bottom_ub49`.
**`w_c` is the axis**: 0.05 → 0.15 moves the inert total by 29 (`surface_stokes`),
44 (`bottom`) and 13 (`bottom_ub49`) points. `τ_strong` scales the total and
matters less below ~6 h; `τ_calm` = 1 yr and `δ` move it by ≤ 1.1 points at
the 750 m band, so the hard step and no calm background stay. **`trap_wall` is
not a smooth knob near zero** on the surface regime: 1 → 0.5 → 0.25 costs 1.6
and 2.7 points, 0.25 → 0 costs 24.6. The bed regimes, whose inert totals are
already low, step more evenly (`bottom` 2.8 / 2.9 / 6.4; `bottom_ub49` 0.5 /
0.4 / 0.9). **The beached fraction is not a reportable number**; the pattern and the
member range are.

## Limitations

**`trap_wall = 0` is an end-member assumption, not a measured one.** The
production physics says a shoreline classified as entirely rock, cliff or
armoured bank *never* retains a Fucus propagule. Nothing in the HELCOM BRISK /
CLMS classification, or in any observation available here, says that: it says
which substrate is there, not what its retention rate is. The choice is
defensible as one end of the plausible range — hard reflective shore in a
tide-free basin gives drifting thalli little to lodge in — but it is the
strongest statement the parameterisation can make. The intermediate members
`…_tf1_tw0p25_b750` and `…_tf1_tw0p5_b750` bracket the cases where wall shore
retains at a reduced but non-zero rate. Anything that depends on the
*absolute* stranded amount, or on stranding at a specific wall-dominated
stretch, must be read across the trap axis rather than off production alone.

**Bottom drifters use the surface criterion.** `bottom` and `bottom_ub49`
strand on the surface onshore Stokes at their position, which assumes the
wave orbital motion reaches the bed inside the 750 m band. Considered:
scaling `w_on` by the linear-wave bed factor `sech(2kh)` (k from a mean wave
period, h = H0); rejected by the authors in favour of one criterion for all
regimes, the narrow band standing in for shallowness.

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
