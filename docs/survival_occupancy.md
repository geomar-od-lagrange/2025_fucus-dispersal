# Survival-weighted occupancy

The free-drifting particle density with **beaching progressively removed**.
Standard occupancy ([024 counts](hexbinning_and_connectivity.md)) weights every
`(trajectory, obs)` sample equally; here each sample is weighted by the
particle's surviving (un-beached) fraction at that age,

```
S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
```

with the same two-state near-shore rate `τ` as the
[beaching diagnostic](beaching.md) — including its degenerate `trap`, so
onshore wave forcing is the only term that modulates the rate. `A` grows only
inside the near-shore band, so open-water residence is undiluted while weight
lingering in a wave-exposed band decays fast. This is the occupancy analogue of
024e's fractional stranding: **024e records where the weight leaves
(`beach_hex`); 024f records where the still-drifting weight is
(`target_hex`)**. It composes with a future Fucus lifetime `L(t)` — survival
becomes `exp(−A)·L`.

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Build | [`024f_BuildSurvivalOccupancy`](../notebooks/024f_BuildSurvivalOccupancy.py) | the `024d` beaching-forcing sidecar + `024a` key | `HexAgg_survocc_*.parquet` |
| Consume | [`030_SurvivalHeatmaps`](../notebooks/030_SurvivalHeatmaps.py) | survocc parquet + key | PNGs under `Figures/030/` |

`024f` reads **only** the sidecar and the key — no raster, no raw Stokes, no
trajectory zarrs; the per-(particle, hour) ingredients of the rate are cached
by [`024d_BuildBeachingForcing`](../notebooks/024d_BuildBeachingForcing.py)
(schema and rate model: [beaching.md](beaching.md)), so a rate-model member is
a seconds-per-zarr numpy pass. It aggregates with `np.bincount` into a
contiguous hex index (occupancy is *dense* — every obs, every hex — unlike
024e's sparse in-band deposits) and runs to `occupancy_max_days` (default
120 d, ≥ the largest 030 horizon and ≤ the sidecar's `window_days`), **not**
024e's `max_float_days` viability cutoff. Partitioned per
`(regime, year, month)`; `030` pools the monthly partitions selected by its `season` parameter.

## Store schema

`HexAgg_survocc_r<radius>m_<regime>_<year>_mMM_<member>.parquet`, additive
across `release_doy`/month/year. The `_<member>` tag is the same rate-model
string `024e` builds (e.g. `step_wc0p1_ts3_tcinf`), and `030` takes it as an
opaque `member` parameter.

| column | meaning |
|--------|---------|
| `release_doy` | release day-of-year of the originating zarr |
| `age_bin` | `floor(age_days / age_bin_days)` |
| `target_hex` | occupied hex (`024a` key space) |
| `occ` | plain occupancy (samples, = 1 each) |
| `surv` | survival-weighted occupancy (`Σ exp(−A)`) |

Both weights share the same window and hexing, so the surviving fraction
`surv/occ` is a self-consistent per-hex, per-age comparison. `release_hex` is
dropped (unlike 024 counts) to keep the bincount aggregation cheap —
per-origin survival maps would need it back and a different aggregation.

## Reading the maps

The surviving-fraction panel is the mean un-beached fraction of the particles
*occupying* a hex, not a local beaching rate — `A` accumulates along the path
before arrival, so a hex reads low when its supply routed through wave-exposed
near-shore water. The sink field is [beaching.md](beaching.md)'s `beach_hex`.
Conditioning on age bin means every particle shares an elapsed time, so the
contrast is route, not age.

At each horizon `T`, `030` selects the snapshot bin *ending* at `T` —
`age_bin = T/age_bin_days − 1`, i.e. ages in `[T − age_bin_days, T)` —
and draws occupancy vs. survival-weighted (shared `LogNorm`) vs. surviving
fraction (linear 0–1). The drifting fraction falls with age, most steeply where
near-shore residence coincides with onshore waves — with `trap` degenerate the
depletion pattern is a *wave-exposure* field, not a substrate one. Absolute
rates are parameter-sensitive ([beaching.md](beaching.md)); read the *pattern*.

## Cross-references

- [beaching.md](beaching.md) — the shared sidecar, rate model and the endpoint
  (stranding) counterpart store.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the plain
  occupancy (024 counts) this re-weights, and the key schema.
- [visualisations.md](visualisations.md) — `030`'s plot rationale.
