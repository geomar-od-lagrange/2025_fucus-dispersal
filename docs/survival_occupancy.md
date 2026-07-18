# Survival-weighted occupancy

The free-drifting particle density with **beaching progressively removed**.
Standard occupancy ([024 counts](hexbinning_and_connectivity.md)) weights
every `(trajectory, obs)` sample equally; here each sample is weighted by the
particle's surviving (un-beached) fraction at that age,

```
S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
```

with the same near-shore beaching rate `τ = τ0/(trap·g(w_onshore))` as the
[beaching diagnostic](beaching.md). `A` grows only inside the near-shore
band, so open-water residence is undiluted while mass lingering near a
retentive shore decays fast. This is the occupancy analogue of 024d's
fractional stranding: **024d records where the weight leaves (`beach_hex`);
024e records where the still-drifting weight is (`target_hex`)**. It composes
with a future Fucus lifetime `L(t)` — survival becomes `exp(−A)·L`.

## Pipeline

| Stage | File | Reads | Writes |
|-------|------|-------|--------|
| Build | [`024e_BuildSurvivalOccupancy`](../notebooks/024e_BuildSurvivalOccupancy.py) | trajectory zarrs + raw `baltic_highres` Stokes + `024a` key | `HexAgg_survocc_*.parquet` |
| Consume | [`030_SurvivalHeatmaps`](../notebooks/030_SurvivalHeatmaps.py) | survocc parquet + key | PNGs under `Figures/030/` |

`024e` reuses 024d's raster + Stokes machinery but aggregates with
`np.bincount` (occupancy is *dense* — every obs, every hex — unlike 024d's
sparse in-band deposits) into a contiguous hex index. It runs to
`occupancy_max_days` (default 120 d, ≥ the largest 030 horizon), **not**
024d's `max_float_days` viability cutoff. Partitioned per
`(regime, year, month)` and fanned out over the full year × month grid
(sapphire-pinned job); `030` pools the monthly partitions.

## Store schema

`HexAgg_survocc_r<radius>m_<regime>_<year>_mMM.parquet`, additive across
`release_doy`/month/year:

| column | meaning |
|--------|---------|
| `release_doy` | release day-of-year of the originating zarr |
| `age_bin` | `floor(age_days / age_bin_days)` |
| `target_hex` | occupied hex (`024a` key space) |
| `occ` | plain occupancy (samples, = 1 each) |
| `surv` | survival-weighted occupancy (`Σ exp(−A)`) |

Both weights share the same window and hexing, so the surviving fraction
`surv/occ` is a self-consistent per-hex, per-age comparison. `release_hex` is
dropped (unlike 024 counts) to keep the bincount aggregation cheap — per-origin
survival maps would need it back and a different aggregation.

## Reading the maps

At each horizon `T`, `030` selects the snapshot bin `age_bin = T/age_bin_days`
and draws occupancy vs. survival-weighted (shared `LogNorm`) vs. surviving
fraction (linear 0–1). The drifting fraction falls with age — e.g. ~0.86 at
0–10 d to ~0.40 by 110–120 d for the default parameters — most steeply near
retentive (flat) shores. Absolute rates are parameter-sensitive (see
[beaching.md](beaching.md)); read the *pattern*.

## Cross-references

- [beaching.md](beaching.md) — the shared rate model and the endpoint
  (stranding) counterpart store.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the
  plain occupancy (024 counts) this re-weights, and the key schema.
- [visualisations.md](visualisations.md) — `030`'s plot rationale.
