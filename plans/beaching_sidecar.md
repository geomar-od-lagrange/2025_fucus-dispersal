# Beaching forcing sidecar: make the rate model a free re-reduction

Cache the per-(particle, hour) ingredients of the beaching rate once, so
that every choice of rate model — `τ0`, `w_tau`, the functional form of
`s(w)`, band width, viability window, trap weights — is a seconds-per-zarr
numpy pass with no trajectory or Stokes I/O.

## Problem

The pre-sidecar beaching and survival-occupancy builders (now
[`024e_BuildBeaching`](../notebooks/024e_BuildBeaching.py) and
[`024f_BuildSurvivalOccupancy`](../notebooks/024f_BuildSurvivalOccupancy.py))
each re-read the trajectory zarrs and re-sample the raw WAM Stokes field for
every parameter setting. Measured on the production run (surface_stokes,
2016–2019, 292 zarrs, 16.6 M real drifters):

| | per zarr | per sweep member |
|---|---|---|
| today (Stokes sampling + reduction) | 87 s | ~7 h CPU |
| sidecar reducer | 2.5 s | ~13 min CPU |

The consequence today is that the sensitivity sweep is stale relative to
production (a `w_tau` sweep at `τ0 = 48 h` against a production point at
`τ0 = 480 h`) and that the functional form of `s(w)` cannot be revisited
without another 48-partition job. Both `024d` and `024e` compute the same
survival curve `S(t) = exp(−A)` from the same ingredients; only the final
aggregation differs (where weight leaves vs. where surviving weight is).

## What the rate model is a function of

For a real drifter at hour `h`, the rate needs only:

| ingredient | source | needed for |
|---|---|---|
| `w_on` — onshore Stokes (m/s) | raw WAM field, hourly, geodesically filled | `s(w)` |
| `dist` — distance to BSH land (m) | 500 m EPSG:3035 raster | band gate, `band_m` sweeps |
| `hex` — `024a` hex id of the position | raster lookup | deposit / occupancy target |
| `flat` — nearest land is tidal-flat-fronted | raster lookup | `trap` (degenerate today) |
| `disp` — crow-flies displacement from release (km) | positions | travel distance at stranding (diagnostic axis) |

All five are parameter-free with respect to the rate model. Everything
else — `s(w)`, `τ0`, `w_tau`, `band_m`, `max_float_days`, `trap_*`,
`age_bin_days`, the deposit telescoping, the groupbys — is a reduction over
these arrays.

## Sidecar store

One zarr per trajectory zarr, real (water-seeded) particles only:

```
output_root/BeachingForcing/<regime>/<year>/<trajectory zarr stem>.zarr
  dims:   trajectory (real only), obs (0 .. window_days*24)
  coords: trajectory   int32   index into the trajectory zarr's `trajectory` dim
          release_hex  int32   per particle, 024a key space
  vars:   w_on   float16  onshore Stokes, 0 outside the sampled band
          dist   uint8    distance to land in 25 m steps, 255 = beyond band_max_m
          hex    int32    hex id everywhere (-1 = NaN position / outside key)
          flat   bool
          disp   uint8    displacement from release in km, 255 = saturated
  attrs:  release_time, release_doy, regime, band_max_m, window_days,
          raster_dx_m, stokes_fill_max_cells, builder git sha
  chunks: (10000, nobs)   one obs chunk per trajectory block
  codec:  numcodecs.Zstd(level=5), no Blosc
```

Encodings were measured on one production zarr (56,810 real drifters ×
1441 h; scripts kept under `tmp_beaching_sidecar/` until this lands):

| variable | encoding | MB / zarr (60 d) | why |
|---|---|---|---|
| `w_on` | float16, Zstd | 21 | exact zeros, no custom codec; max abs error in beached fraction 2e-8, per-hex ≤ 4e-5 rel. above one unit of weight |
| `dist` | uint8 / 25 m, 255 beyond 5 km | 6.5 | 25 m ≪ raster dx; nulling beyond the band halves it |
| `hex` | int32, everywhere | 5.4 | survival occupancy needs the target hex in open water too, so no nulling; key has 190 k ids so int16 does not fit; runs compress 200× |
| `flat` | bool | 0.3 | |
| `disp` | uint8 km | ~3 (est.) | slowly varying → compresses like `hex` |

Rejected: uint8 log-quantised `w_on` (11 MB, 1.2e-3 per-hex error — fine,
but a custom codec for 10 MB/zarr); Blosc-zstd-bitshuffle (2–3× larger than
raw Zstd on every array here); sparse (COO) layout (73–81 % of particle-hours
are in band, dense wins); storing `s(w)` or `Σ s` instead of `w_on` (freezes
`w_tau` and the functional form).

Nulling `dist` and `w_on` beyond a distance threshold buys little because Fucus releases are
coastal and particles hug the shore: 73 % of particle-hours lie within 2 km,
81 % within 5 km. Sample at `band_max_m = 5000` (the widest band the coarse
grid can resolve, see `beaching.md`) and accept ~19 % dead cells.

**Window.** `window_days = 120` so the same sidecar serves the survival
occupancy horizons (`occupancy_max_days = 120`). Beaching's viability window
(`max_float_days = 60`) is a reducer parameter that slices `obs`.

**Size.** ~30 MB per zarr at 60 d, ~60 MB at 120 d → ~18 GB for
surface_stokes 2016–2019. Trajectories are 98 GB per regime-year.

## Notebooks

Renumber so the dependency order reads top to bottom (no produced data
needs preserving):

| stage | notebook | reads | writes | cost / zarr |
|---|---|---|---|---|
| build | `024d_BuildBeachingForcing` | trajectory zarr, raw Stokes, H0 statics, 024a key | sidecar zarr | ~50 s (17 s raster lookups + 32 s Stokes at 5 km) |
| reduce | `024e_BuildBeaching` | sidecar, key | `HexAgg_beaching_*.parquet` | 2.5 s |
| reduce | `024f_BuildSurvivalOccupancy` | sidecar, key | `HexAgg_survocc_*.parquet` | ~2 s |
| consume | `029`, `030`, `031` | parquet | PNGs | unchanged |

`024d` keeps the raster build, the WAM geodesic fill (`OnshoreStokes`), the
land-seeded filter, and the Stokes loop; it drops the rate model and the
groupby entirely. `024e` is today's `deposit_one_zarr` from the `s(w)` line
onward, reading arrays from the sidecar instead of computing them. `024f`
is today's `024e` bincount pass on `exp(−A)`. Both reducers take the same
rate-model parameters and must share one implementation of `A(t)`; since
notebooks own their utilities and this is ~10 lines (`s`, `a`, `cumsum`),
keep it inline in both and let the validation step below hold them equal.

### The rate model becomes a two-state hazard

The saturating ramp `s = 2w/(w + w_tau)` cannot express "drifts freely for a
year under ordinary coastal conditions, strands within hours in a storm":
it is steepest at weak forcing and flat at strong forcing, so `τ` between
the median and the strongest measured onshore Stokes differs by only ~2.7×
whatever `τ0` is (`τ0 = 480 h` gives 30 d at p50 and 11 d at the maximum).
Replace it with

```
1/τ = 1/τ_calm + trap · r(w_on) / τ_storm          in band, else 0
r(w) = clip((w − w_c + δ/2) / δ, 0, 1)             step at w_c, optional
                                                   linear edge of width δ
```

| parameter | meaning | how to choose |
|---|---|---|
| `w_c` | onshore Stokes above which stranding is on | against the measured in-band forcing (p75 0.057, p90 0.088, p99 0.15 m/s); **the** sweep axis |
| `τ_storm` | e-folding time while `w_on ≥ w_c` | hours; below ~6 h everything in band during a Baltic storm (tens of hours) strands and the value stops mattering |
| `τ_calm` | background rate in band regardless of forcing | ∞ or O(1 yr); a year loses 15 % over 60 d and 45 % over 220 d, so it is a real parameter — sweep {∞, 1 yr} |
| `δ` | edge width | 0 by default; check once whether `δ ≈ 0.02` changes results (guards against nearest-hour / nearest-cell sampling noise in `w_on`); if not, keep the step |
| `trap` | shore type factor | degenerate at 1 as today, seam for a substrate map |

Nothing in the data constrains the shape of `r` between calm and storm, so
a shape parameter (Hill exponent, half-saturation) is pure nuisance; the
step is the n → ∞ limit of any such form and its one knob means what it
says. Once `τ_storm ≪` storm duration the model degenerates gracefully into
an **exposure model**: the beached fraction is the share of particles that
were ever in band during a storm hour inside the viability window, and the
stranding map is coastal residence × storm climatology. Expect strong
seasonal and interannual signal; the pooled-year view stops being the
default and the per-month / per-year splits become the primary products.

The saturating form stays in the reducer only until the validation below
has reproduced the existing stores from the sidecar, then it is deleted.
The store filename tag becomes a generic `member` string
(`_wc<w_c>_ts<τ_storm>_tc<τ_calm>[_d<δ>]`) and `031` takes an explicit
list of member tags plus a display label per member, instead of deriving
members from a `w_tau` list at a fixed `τ0`.

### Travel distance at stranding

`disp` in the sidecar gives every deposit a crow-flies travel distance, so
the beaching reducer can add a `disp_bin` axis (km) to the store or emit a
per-source distribution of "how far had the stranded weight travelled". A
storm stranding freshly released material at home is a real outcome and is
*not* to be filtered out — the "confounded by release geometry" paragraph in
`beaching.md` conflates two things: calm-water stranding at release under a
nearly forcing-independent rate (an artefact of the saturating form at low
`w_tau`), and storm stranding at release (signal). With a threshold form only
the second survives. `disp` is kept as a diagnostic axis, never as a mask.

## Job scripts

- `024d_BuildBeachingForcing_job.sh`: today's 024d fan-out (`year × month`
  cells via `xargs`, 3-attempt kernel retry, `--spread-job`), no rate
  parameters, one sidecar per zarr. Runs once per `(regime, year)`.
- `024e_BuildBeaching_job.sh`, `024f_BuildSurvivalOccupancy_job.sh`: take a
  member specification, fan out `member × year × month`; each cell is
  seconds, so `--ntasks` small and no sapphire pin needed. Consider one task
  per `(member, year)` looping months inside the notebook (`release_month = 0`
  already does a whole year) — 73 zarrs × 2.5 s ≈ 3 min.
- Add `029`, `030`, `031` wrappers for parity with the other parquet-only
  consumers (out of this plan's scope but noted in the PR wrap-up).

## Validation

1. Build the sidecar for one production zarr and reduce with
   `(saturating, τ0 = 480 h, w_tau = 0.05, band 2 km, 60 d)`; the per-zarr
   line must match the executed production notebook
   (`56,810 drifters, 31,206 beached (54.9 %)` for
   `Fucus_BSH_20190804T001500_surface_stokes_dt60min_seed294209566`).
2. Reduce the full 2019 surface_stokes set and compare the pooled beached
   fraction and per-hex weights against the existing `t480_wt0p05` store
   (expect float16-level differences only).
3. `024f` against the existing `survocc` store, same tolerance.
4. Then delete `HexAggregates/HexAgg_beaching_*`, `HexAgg_survocc_*`,
   `HexAggregates/superseded_*`, and rebuild from the sidecar.

## Docs and follow-through

- `docs/beaching.md`: pipeline table gains the sidecar stage; the "`τ0` is
  free … not yet implemented" paragraph goes; rate-model section gains the
  `s_form` table; the sweep section is rewritten around whichever members
  are actually built. Fold the sidecar schema in here rather than a new doc
  (it is the beaching store's substrate, as the key is the hex stores').
- `docs/survival_occupancy.md`: reads the sidecar, no Stokes machinery.
- `docs/job_scripts.md`: 024d/024e/024f scheduling; the GPFS-bound notes
  now apply to the builder only.
- `AGENTS.md` pipeline stages: add 024c–024f, 028–031.
- Grep for `024d`/`024e` in notebooks 029–031, scripts, docs, plans and
  update every reference in the same pass.

## Open decisions

- Production values of `w_c`, `τ_storm`, `τ_calm`: decided from the sweep on
  the realised exposure statistics (share of in-band hours above `w_c`,
  stranding-age and travel-distance distributions, seasonal split), not on
  the beached total.
- Whether `δ > 0` is needed (one comparison, see above).
- Whether to keep `flat`. It costs nothing; kept as the seam for a substrate
  classification, as today.
