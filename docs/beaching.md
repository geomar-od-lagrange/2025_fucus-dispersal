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

## The rate model

For each real drifter, positions are walked up to a viability window
(`max_float_days`, capping the biologically live fraction of each 220-day
trajectory). Inside a near-shore band the beaching hazard is a Δt-invariant
e-folding rate:

```
τ = τ0 / (trap(shore_type) · s(w_onshore))
```

with `trap ≡ 1` as configured, so in practice `τ = τ0 / s(w_onshore)` and
**onshore wave forcing is the only term that modulates the rate**.

### What is borrowed and what is new

Being precise about provenance, because the two halves have very different
standing:

| element | status |
|---|---|
| Δt-invariant e-folding hazard, `p = 1 − exp(−Δt/τ)` | Standard. Onink et al. 2021 (also Hernandez et al. 2024, Daily et al. 2021, Siht et al. 2025). |
| `trap(shore_type)` — terrain-dependent rate | Has precedent (Onink, Siht). **Deliberately switched off here** — see below. |
| `s(w_onshore)` — forcing-dependent rate | **This study's extension. No precedent in the cited set**, all of which effectively use a constant rate (`s ≡ 1`). |

So the scheme **inverts** the literature's structure: the cited work varies
the rate by terrain and holds it constant in forcing; this varies it by
forcing and holds it constant in terrain. That is a deliberate choice, not an
inherited one, and it must not be presented as following Onink et al.

Two of the cited papers argue *against* forcing modulation. Onink et al.
2021 keep the scheme "simplest possible" on principle. Daily et al. 2021
object that a stochastic beaching parameterisation is *intended* to represent
the unresolved near-shore processes, so feeding a resolved variable back in
risks double-counting them.

**The answer to the double-counting objection is specific to this setup, and
is the reason the extension is defensible here.** The `surface_stokes`
trajectories do not contain the cross-shore Stokes transport: it is zeroed at
blocked faces by the advection scheme's land mask
([stokes_drift.md](stokes_drift.md)). The onshore push that would physically
drive propagules ashore is therefore *resolved in the forcing data but removed
from the trajectories*. Modulating the beaching rate by `w_onshore` restores
a transport that was deliberately suppressed — it is not a resolved variable
layered on top of a stochastic term already covering it. Daily's objection
applies to schemes where the near-shore transport is present in the drift;
here it is not.

That is an argument, not a validation. There is no observational constraint
on `s` in this study, which is why the headline number is reported as a sweep
range over `w_half` ([`031`](../notebooks/031_BeachingSweep.py)) rather than
as a value.

`s` is named for "saturation"; it was `g` earlier, which collides with
gravitational acceleration and read as a generic helper rather than the
model's one novel term.

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
  weights apart to activate. Treat `shore_type` in the store and in `029`'s
  summary as a diagnostic label, never as a model result; `029` deliberately
  does not map the split.

  Switching off a term the literature *does* support needs a defence, and
  there are three. **Validity**: the only typing available is BSH's `H0 ≤ 0`
  tidal-flat flag, which is not a retentiveness proxy for Baltic shores — the
  basin is effectively tide-free, so the flag fires essentially only in the
  German Bight, which is outside the wave grid and has zero forcing anyway. A
  two-class split would read as resolved coastal morphology while resolving
  none. **Measured effect size**: switching `trap_flat` from 2.0 to 1.0 moved
  the beached total by 0.9 points (57.3 % → 56.4 %, Aug 2019, Euclidean fill),
  because `flat` covers 8,301 of 4,178,297 water cells and takes 6.5 % of
  beached weight. **Literature effect size**: Daily et al. 2021 and Onink et
  al. 2021 both report terrain variation having minimal impact on outcomes.
  So the term is being dropped where it is least load-bearing and least
  supportable, not to simplify the model.
- **Onshore wave forcing** sets `g`, and with `trap` degenerate it is the
  *only* term that modulates the rate. The onshore Stokes component
  `w_onshore = max(0, −(VSDX, VSDY) · n_out)`, sampled from the **raw
  `baltic_highres` field** (CMEMS `BALTICSEA_MULTIYEAR_WAV_003_015`, FMI-WAM,
  1 nmi ≈ 1.6 × 1.9 km, hourly) by nearest hour and cell — the cross-shore
  transport the `surface_stokes` runs zeroed at blocked faces
  ([stokes_drift.md](stokes_drift.md)), sampled *before* that mask and
  *without* the N=5 spread (which bridges thin land barriers).
  `g(w) = w / (w + w_half)` is a saturating ramp.

  **The WAM field is extrapolated to full BSH coverage**, by geodesic
  (through-water) propagation. WAM's water mask is not BSH's, and its bbox
  starts at 9.01 °E, excluding the German Bight. Left alone, positions on WAM
  land return `w_onshore = 0`, which for beaching means *rate zero* — a
  coastline that cannot strand at any `τ0`/`w_half`. That is a structural
  bias, not a parameter choice, and it is large: **78.7 % of in-band samples
  need filling**, because WAM masks its own coastal cells and the 2 km band
  hugs the shore.

  The fill is a breadth-first expansion of the **donor index** over the static
  masks — one 4-neighbour dilation per round, front advancing ~1.6 km, BSH
  land blocking it. Propagating indices rather than values keeps it a one-off
  precompute; per-hour sampling stays a single gather. Three details carry
  weight:

  - **Geodesic, not Euclidean.** A nearest-wet lookup by straight-line
    distance can draw from across a headland or outside a fjord mouth,
    importing an open-water wave climate into sheltered water that WAM
    excluded *because* it is not open water. Geodesic donors instead share
    the particle's water body. **Measured effect: none.** Rebuilding the full
    sweep on the geodesic fill reproduced the Euclidean numbers to within
    0.15 points on every member (e.g. 60.1 % → 60.2 % at `w_half = 0.4`), and
    fill distance barely moved (2.74 → 2.67 km). So Euclidean nearest-wet
    rarely crossed land in a way that mattered. Geodesic is kept because it
    cannot do so *by construction* and because it distinguishes unreachable
    water (capped and counted) from reachable, which the Euclidean version
    could not — not because it changes the answer.
  - **4-neighbour, not 8.** A 3×3 dilation squeezes between diagonally
    touching land cells — the same thin-barrier bridging the `surface_stokes`
    N=5 Stokes spread is faulted for.
  - **Capped and counted.** Cells not reached within `stokes_fill_max_cells`
    rounds keep `w_onshore = 0` and are reported (currently 0.09 % of
    samples), rather than silently importing from an absurd distance as the
    Euclidean version did.

  The mask that makes this safe is that **WAM NaN means land *or* ice** — the
  wet-cell count varies hour to hour, ~4.4 % of the grid seasonally
  ice-blanked in the Bothnian Bay and Gulf of Finland. Only *static* land is
  filled, identified by an `ever_wet` mask (finite in any hour of a monthly
  sample spanning the seasonal cycle). **Ice keeps `w_onshore = 0`**, the
  physical answer: ice suppresses waves, so no-beaching-under-ice falls out
  for free, consistent with the currents side where Fucus rides the upper-cell
  velocity.

  The BSH water mask on the WAM grid uses a **3×3 footprint** rule, not centre
  sampling: WAM cells are ~1.6 km and coastal ones are part land, so a centre
  landing on BSH land would exclude cells particles legitimately occupy, which
  would then never receive a donor and be silently forced to zero.

  Known bias: filled values are read 1–3 km offshore, where Stokes drift is
  stronger than at the shoreline. Zero was a large negative bias; this is a
  smaller positive one. A depth taper is *not* the obvious remedy — Baltic
  wind seas are short-period (Tp 3–6 s), so depth-limited breaking is confined
  to a strip far narrower than the 500 m raster or the 1.6 km WAM cell.

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

- **`τ0` and `w_half` are not independent, and here they are nearly
  degenerate.** `τ = τ0 + τ0·w_half/w`: `τ0` is an additive floor, the
  *product* `τ0·w_half` is the weak-wave coefficient. Where `w ≪ w_half` only
  the product matters. The measured forcing puts us there: over in-band steps
  **only 57.5 % have `w_onshore > 0` at all** (elsewhere waves are offshore,
  alongshore, or ice-suppressed), and among those the distribution is
  p10/p50/p90 = **0.0037 / 0.0256 / 0.0884 m/s** — a median below every sweep
  member. So for `w_half ≳ 0.2` the `w_half` sweep is effectively a sweep of
  the rate scale `τ0·w_half`, not an independent axis; reporting it as
  "sensitivity to `w_half`" would overstate what varies. `024d` prints this
  distribution and the realized `τ` quantiles every run so the regime is
  visible rather than assumed.
- **Spatial selectivity is what discriminates**, since totals trade along the
  ridge: 031's Gini of stranded weight rises 0.668 → 0.772 across the sweep
  while the stranding support stays flat (~2,200 hexes). Large `w_half`
  concentrates stranding on wave-exposed coast; small `w_half` saturates
  `s → 1` and gives a pure near-shore-residence map.
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
