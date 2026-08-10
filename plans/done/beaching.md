# Beaching (post-simulation) — design note

**Implemented.** See [../../docs/beaching.md](../../docs/beaching.md) for
the current state (`024d_BuildBeaching` + `029_BeachingMaps`). This note is
the original design rationale and literature basis. **The production scheme
is the weighted / fractional variant** — the deterministic expectation of
the first-stranding rule sketched below, chosen for noise-free high-age
tails and clean composition with a lifetime distribution; the stochastic
per-particle draw described here was retired after a code-correctness
comparison (the two agree to <0.2 % on totals, within MC noise per hex).

## Purpose

Add a beaching diagnostic to the *Fucus* dispersal study: estimate where
and when drifting propagules strand on the coast, and re-aggregate
occupancy/connectivity with beached particles removed from the free-
drifting pool. Beaching is modelled as a **first-stranding endpoint** —
once a particle strands it leaves the drifting pool and is not refloated.
This is framed as a *dispersal endpoint*, not a claim of permanent
physical attachment: beach-cast wrack in the tide-free Baltic is
demonstrably remobilised by wind and water-level changes
(Hammann & Zimmer 2014), but a fragment that strands then desiccates,
degrades, or is reworked has largely stopped contributing to onward
dispersal, so treating first stranding as terminal is defensible on
viability grounds. A reversible-resuspension sensitivity case bounds the
effect (see caveats).

## Why post-simulation, not a Parcels kernel

The existing runs are the expensive part (73 releases/yr × 220-day sims ×
3 regimes × 4 years, all on NESH under
`2025_fucus_dispersal_outputs/Trajectories/<regime>/<year>/`). Beaching
is a cheap function of positions we already have, and its parameters
(trapping thresholds, wave-forcing weighting) will need tuning. Baking it into a
kernel would force a full re-run per parameter choice. As a
post-processing pass over the trajectory zarrs it is re-runnable in
minutes and keeps the physics runs untouched — the same "aggregate once,
consume cheaply" split the [hex store](../../docs/hexbinning_and_connectivity.md)
already uses.

## The three-ingredient beaching model

For each trajectory, walk its `(obs)` positions in time. At every step
inside a near-shore band, compute an instantaneous beaching propensity
from three factors; strand the particle at the first step that trips the
rule; drop all later obs. The three ingredients and their data readiness:

**Viability gate (applied first).** Cap each trajectory's contributing
age at a *Fucus*-realistic float-and-viability window before scoring
beaching. Floating *F. vesiculosus* stays buoyant and physiologically
viable for weeks to a few months, buoyancy and viability declining through
the season as gametes are shed and vesicles flood (Rothäusler et al. 2019);
a raw 220-day trajectory vastly outlives that window, so its late positions
are biologically dead and would inflate long-range strandings. The clock
is a single `max_float_days` parameter (seasonally varying if wanted)
applied on the `age`/`obs` axis already present in the zarrs.

### 1. Distance to shore — *ready*

Score against a **per-position** distance to the fine BSH coastline, not
the per-hex `dist_to_coast_m` in the [`024a` key](../../notebooks/024a_BuildHexKey.md).
Beaching is a metres-scale swash-zone process; a per-hex distance
quantises it to the hex radius and ties the near-shore band width to the
aggregation grid rather than to a physical length. A rasterised
distance-to-coast field (from the fine coastline) sampled at each obs
keeps the band a named physical parameter. The per-hex value stays useful
as a coarse cross-check only.

### 2. Shore classification: wall vs. flat — *data in hand*

Separable from the two BSH coastlines already shipped in the data twin
(`data/bsh_hbmnoku_static/`, produced by
[`004_extract_coastline.py`](../../notebooks/004_extract_coastline.py)):

- `coastline_always_wet.geojson` — cells with `H0 > 0` only.
- `coastline.geojson` — includes tidal-flat cells (`H0 ≤ 0`, see
  [h0_semantics.md](../../docs/h0_semantics.md)).

A shore fronted by tidal flats (in the second but not the first) reads as
**flat**; a steep always-wet edge reads as **wall**. The weight
`trap(shore_type)` should capture **retention/reflectivity at given wave
forcing** — does what the waves deliver stay put (flat: dissipative,
accumulates wrack) or wash back out (wall: reflective) — **not** how
exposed the shore is. Exposure/fetch is already carried by the
onshore-Stokes magnitude (ingredient 3); letting `trap` re-encode it would
confound the two (see caveats). This yields a per-coastal-hex (or
per-segment) `shore_type ∈ {flat, wall}` weight with no new data.

### 3. Onshore wave forcing — *cross-shore Stokes drift*

The shoreward driver is the onshore component of the wave **Stokes
drift**, not a wind proxy. It is the same CMEMS wave field the
`surface_stokes` runs already use (`VSDX/VSDY`, downloaded by
[`002_download_stokes.py`](../../notebooks/002_download_stokes.py)) — and
specifically the *cross-shore transport the runs deliberately suppressed
at the coast*: [`003`](../../notebooks/003_prepare_2d_fields.py) zeros the
Stokes contribution on blocked BSH faces so it can't push particles
through no-slip walls (see [stokes_drift.md](../../docs/stokes_drift.md)),
which is exactly the onshore push that strands material. Beaching puts it
back:

```
w_onshore = max(0, -(VSDX, VSDY) · n_out)
```

with `n_out` the outward shore normal (from coastline segment geometry or
`∇(dist_to_coast)`). Offshore Stokes (`w_onshore = 0`) suppresses beaching.

**Why Stokes, not 10 m wind.** Wave Stokes is the physically correct
cross-shore transport, and it is **fetch-aware**: the Baltic wave model
grows waves along open fetch and returns ~0/NaN in fetch-starved enclosed
water, so a sheltered lee shore correctly sees little onshore transport. A
raw 10 m-wind proxy is fetch-blind and would over-force sheltered shores.

**Sample the *unmasked* Stokes.** Two implementation points: (a) read the
Stokes *before* the blocked-face mask — the 2D field the runs wrote is
zero at exactly the coastal faces beaching needs — i.e. sample the raw
CMEMS `VSDX/VSDY`; (b) use the raw 2 km wave-model field honouring its
land/NaN topology, **not** the N=5 rolling-mean spread, whose known
land-bridging (open-ocean Stokes leaking across the Curonian Spit;
[stokes_drift.md](../../docs/stokes_drift.md) "open concern") would
manufacture onshore transport where the real fetch is zero.

## Beaching rule (rate-based, Δt-invariant)

Express the propensity as a **beaching rate** `1/τ`, not a per-step
constant, so the outcome is invariant to the zarr output interval `Δt`
(a raw per-step probability would make total beaching depend on the
sampling cadence). Inside the near-shore band the per-step beaching
probability is

```
p_step = 1 − exp(−Δt / τ),    τ = τ0 / (trap(shore_type) · g(w_onshore))
```

with `trap(flat) > trap(wall)` and `g` monotonic in onshore Stokes (a
saturating ramp), so a flat shore under strong onshore waves gives a short
`τ` (fast beaching) and a wall under offshore Stokes a long one. Strand at
the first step whose draw against `p_step` succeeds; drop later obs. This
is the plastics community's e-folding scheme (Onink et al. 2021), which
anchors `τ0` at O(days) and the band at ~one coastal grid cell. One RNG
seed per run, printed, per the repo convention.

## Optional fallback: ERA5 wind via the CDS ARCO store

The Stokes driver (ingredient 3) is primary. ERA5 10 m wind is kept only
as an **optional fallback / cross-check** — a fetch-blind sensitivity
comparison, or filling a gap where wave Stokes is unavailable. It reads
lazily from the Copernicus Climate Data Store analysis-ready
cloud-optimised (ARCO) Zarr store. The store is **Zarr v2 (consolidated)**
and reads with the project's existing `zarr 2.18.7` — no separate
environment. Geo-chunked surface store:

```
sfc: https://arco.datastores.ecmwf.int/cadl-arco-geo-002/arco/reanalysis_era5_single_levels/sfc/geoChunked.zarr
```

Access needs the CDS API key (`~/.cdsapirc`) as a Bearer token, plus
`trust_env` so the async HTTP client uses the NESH `https_proxy`:

```python
import pathlib, xarray as xr
key = next(l.split(":", 1)[1].strip()
           for l in (pathlib.Path.home() / ".cdsapirc").read_text().splitlines()
           if l.strip().startswith("key:"))
ds = xr.open_zarr(
    sfc_geo_url,
    consolidated=True,
    storage_options={
        "headers": {"Authorization": f"Bearer {key}"},
        "client_kwargs": {"trust_env": True},
    },
)
```

`u10`/`v10` sit on coords `time` / `latitude` (descending) / `longitude`
(−180…180). (DestinE Earth Data Hub offers a fuller ERA5 Zarr mirror, but
needs a separate account and is Zarr v3.)

## Scientific basis and caveats

The design follows established Lagrangian beaching practice, mostly from
the marine-litter community, with *Fucus*-specific adaptations:

- **Form is literature-consistent.** A coastal-buffer + stochastic
  e-folding rate is the dominant scheme (Onink et al. 2021); permanent
  removal on a distance-to-shore crossing is the LOCATE/Parcels approach
  (Hernandez et al. 2024); terrain-dependent trapping (the wall/flat
  weight) has direct precedent (Daily et al. 2021). Totals are highly
  parameter-sensitive — in the Baltic the beaching scheme, not the
  hydrodynamics, can dominate the answer (Siht et al. 2024) — so `τ0`, the
  band width, and `trap(flat/wall)` warrant a small sensitivity sweep,
  reported as a range.
- **Driver is consistent for `surface_stokes`.** Onshore Stokes is not a
  foreign windage term: it is the same wave field the `surface_stokes` runs
  advected with, restricted to the near-shore cross-shore component they
  masked. (A raw 10 m-wind/windage driver would be inconsistent — floating
  *Fucus* has real windage ~2–4 %, Wagner et al. 2022 / Wang et al. 2026,
  but the drift never felt it.) For the pure `surface`/`bottom` regimes the
  Stokes driver adds a wave field the trajectories never felt, so those are
  sensitivity variants, not the baseline.
- **Exposure is single-sourced by design.** Onshore-Stokes magnitude and
  the `H0`-slope `wall/flat` split both correlate with wave exposure, but
  enter the rate with opposing signs (exposed = strong push, low retention;
  sheltered = weak push, high retention) and partly cancel in the product.
  To keep them separable, Stokes carries exposure/fetch and `trap` carries
  retention/reflectivity only. The `H0 ≤ 0` split is itself a slope proxy
  (imperfect — dissipative open beaches are gentle yet exposed), so
  validate the *pattern*, don't over-read `trap`.
- **Sub-grid resolution is the main validity threat.** Real beaching is a
  surf/swash-zone process below the BSH grid and coastline; particles pile
  up at the wet–dry boundary (the near-shore-slowdown artefact) and
  beaching scored on top inherits it. Mitigate with the fine coastline +
  raster distance, and validate the spatial *pattern* against independent
  wrack observations rather than absolute rates or locations (stranding ≠
  source proximity; López et al. 2017).
- **Simplifications genuinely safe here.** Baltic tides are negligible (an
  advantage over ocean-coast studies); wind-driven water-level/seiche
  stranding, the real Baltic wrack driver (Hammann & Zimmer 2014), is
  partly captured through the onshore-Stokes term.

## Pipeline placement

Mirror the store-then-consume split:

- **`024d_BuildBeaching`** — heavy per-`(regime, release_year)` pass: read
  the trajectory zarrs, the `024a` key, and the **raw (unmasked) CMEMS
  Stokes** already on disk (`002` output under `output_root/stokes/…`),
  apply the beaching rule, write a beaching store. A Dask-cluster notebook
  like `024`/`024b`. No new download — the wave field is already staged.
- **`029_BeachingMaps`** — lightweight parquet-only consumer: beaching
  density / where-stranded maps and a beached-vs-drifting fraction, in the
  `025`–`028` lineage.

Downstream, occupancy/connectivity can be recomputed against beached-
truncated trajectories to show how stranding reshapes the reachable set.

## Beaching store schema (draft)

Per stranded trajectory, one row (`-1`/NaN sentinels as in the counts
store):

| column | meaning |
|--------|---------|
| `release_hex` | release hex of the trajectory (`024a` key space) |
| `release_doy` | release day-of-year of the originating zarr |
| `beach_hex` | hex where it stranded; `-1` if never beached |
| `beach_age_days` | elapsed time from release to stranding; NaN if never |
| `shore_type` | `flat` / `wall` at the stranding site |

Additive across `release_doy` and year like the other `024x` stores, so
`029` pools the same way `026`/`028` do.

## Open questions / decisions for the human

- **Float-and-viability window** — `max_float_days`: a single value vs.
  seasonally varying; anchor from Rothäusler et al. 2019 (weeks to a few
  months). Sets how much of each 220-day trajectory counts.
- **Rate parameters** — `τ0` (O(days)), the near-shore band width
  (~one coastal cell), and `trap(flat)` / `trap(wall)` as a
  reflectivity/retention contrast. Sweep and report a range.
- **Stokes sampling** — sample the raw 2 km CMEMS Stokes honouring its
  land/NaN topology (lagoon-NaN → ~0 onshore transport); do **not** use the
  N=5 spread field, which bridges thin land barriers.
- **Terminal vs. reversible** — keep first-stranding as the headline;
  test one Onink-style resuspension timescale to bound how much permanence
  inflates coastal retention.
- **Coast normal** — from coastline segment geometry vs. `∇(dist_to_coast)`.
- **Baseline regime** — `surface_stokes` (beaching driver = the runs' own
  wave field); `surface`/`bottom` are sensitivity variants (Stokes was not
  in their drift).

## References

DOIs for the lit database; all verified.

1. Onink et al. (2021), "Global simulations of marine plastic transport
   show plastic trapping in coastal zones," *Environ. Res. Lett.* 16(6):064053 —
   `10.1088/1748-9326/abecbd`. Canonical e-folding beaching+resuspension scheme.
2. Hernandez et al. (2024), "LOCATE v1.0: numerical modelling of floating
   marine debris dispersion in coastal regions using Parcels," *Geosci.
   Model Dev.* 17:2221–2245 — `10.5194/gmd-17-2221-2024`. Parcels
   distance-to-shore beaching with terminal removal.
3. Daily et al. (2021), "Incorporating terrain specific beaching within a
   Lagrangian transport plastics model for Lake Erie," *Microplast.
   Nanoplast.* 1:19 — `10.1186/s43591-021-00019-7`. Terrain/shore-type
   trapping precedent (wall/flat).
4. Wagner et al. (2022), "How Winds and Ocean Currents Influence the Drift
   of Floating Objects," *J. Phys. Oceanogr.* 52(5):907–916 —
   `10.1175/JPO-D-20-0275.1`. Derives the ~2–4 %-of-wind windage rule.
5. Wang et al. (2026), "Quantifying Drag on Floating Sargassum Patches
   Under Combined Air and Water Forcing," *JGR Oceans* 131 —
   `10.1029/2025JC023529`. Measured macroalgal windage 1.1–2.4 %.
6. Rothäusler et al. (2015), "Abundance and dispersal trajectories of
   floating *Fucus vesiculosus* in the Northern Baltic Sea," *Limnol.
   Oceanogr.* 60(6):2173–2184 — `10.1002/lno.10195`. Closest prior art.
7. Rothäusler et al. (2019), "It takes two to stay afloat: … long-term
   floating dispersal of the bladderwrack *Fucus vesiculosus*," *Eur. J.
   Phycol.* 55(2) — `10.1080/09670262.2019.1694706`. Finite,
   morphology-dependent floating duration (viability clock).
8. López et al. (2017), "The variable routes of rafting: stranding
   dynamics of floating bull kelp *Durvillaea antarctica* …," *J. Phycol.*
   53(1) — `10.1111/jpy.12479`. Stranding governed by shore morphology;
   stranding ≠ source proximity.
9. Hammann & Zimmer (2014), "Wind-Driven Dynamics of Beach-Cast Wrack in a
   Tide-Free System," *Open J. Mar. Sci.* 4:68–79 —
   `10.4236/ojms.2014.42009`. Western-Baltic wrack is wind-remobilised
   (challenge to no-remobilisation).
10. Siht et al. (2024), "Modeling the pathways of microplastics in the
    Gulf of Finland, Baltic Sea – sensitivity of parametrizations," *Ocean
    Dyn.* 75 — `10.1007/s10236-024-01649-0`. Beaching parametrisation
    dominates fate in this basin.

## Cross-references

- [hexbinning_and_connectivity.md](../../docs/hexbinning_and_connectivity.md) — the `024x` store pattern and key schema this reuses.
- [h0_semantics.md](../../docs/h0_semantics.md) — the `H0 ≤ 0` tidal-flat rule behind the wall/flat split.
- [stokes_drift.md](../../docs/stokes_drift.md) — the wave Stokes field that drives beaching, and the blocked-face mask this diagnostic reverses at the coast.
- [../notebooks/004_extract_coastline.py](../../notebooks/004_extract_coastline.py) — produces the two coastlines the shore classification reads.
