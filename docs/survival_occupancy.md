# Survival-weighted occupancy

The free-drifting particle density with **beaching progressively removed**.
Standard occupancy ([024 counts](hexbinning_and_connectivity.md)) weights every
`(trajectory, obs)` sample equally; here each sample is weighted by the
particle's surviving (un-beached) fraction at that age,

```
S(t) = exp(−A(t)),   A(t) = cumsum(Δt/τ) over in-band steps
```

with the same two-state near-shore rate `τ` as the
[beaching diagnostic](beaching.md) — including its `trap` shore-type factor,
which at the production weights (`trap_flat = 1`, `trap_wall = 0`) is
`ff/100`, so the rate is modulated by onshore wave forcing **and** by how flat
the local shoreline is; weight sitting in a wave-exposed band off a rock coast
loses nothing. `A` grows only inside the near-shore band, so open-water
residence is undiluted while weight lingering in a wave-exposed band off flat
shore decays fast. This is the occupancy analogue of
024e's fractional stranding: **024e records where the weight leaves
(`beach_hex`); 024f records where the still-drifting weight is
(`target_hex`)**. **No lifetime is applied**: `max_age_days` is a reduction
window, and the store stays resolved per `age_bin`, so a biological survival
`L(t)` can be multiplied in downstream as `exp(−A)·L`.

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
024e's sparse in-band deposits) over the same reduction window as 024e and
024g, `max_age_days` (180 d, ≥ the largest 030 horizon and ≤ the sidecar's
`window_days` of 220 d). The `024a` key supplies both the dense hex index and
the origin axis (`helcom_subbasin` of the release hex). Partitioned per
`(regime, year, month)`; `030` pools the monthly partitions selected by its
`season` parameter.

## Store schema

`HexAgg_survocc_r<radius>m_<regime>_<year>_mMM_<member>.parquet`, additive
across `release_doy`/month/year. The `_<member>` tag is the same rate-model
string `024e` builds (production: `step_wc0p1_ts3_tcinf_tf1_tw0`), and `030`
takes it as an opaque `member` parameter. `024f` must be run at the same rate
parameters as `024e` and `024g`, or the survival and stranding stores describe
different models.

| column | meaning |
|--------|---------|
| `release_doy` | release day-of-year of the originating zarr |
| `age_bin` | `floor(age_hours / (age_bin_days · 24))`, `0 … ⌈max_age_days / age_bin_days⌉ − 1` |
| `origin_subbasin` | HELCOM subbasin id of the drifter's release hex (key `helcom_subbasin`); `-1` = release hex outside every named polygon or not in the key — the same assignment as 024g |
| `target_hex` | occupied hex (`024a` key space); `-1` = no in-key position (the drifter left the BSH domain or has no position) |
| `occ` | plain occupancy (samples, = 1 each) |
| `surv` | survival-weighted occupancy (`Σ exp(−A)`) |

Both weights share the same window and hexing, so the surviving fraction
`surv/occ` is a self-consistent per-hex, per-age comparison.

**Exact denominator.** The `target_hex = -1` sentinel rows are kept, so every
released drifter contributes `age_bin_days · 24` samples to every full age
bin, and

```
Σ occ over (origin_subbasin, target_hex) at fixed (release_doy, age_bin)
  = released water-seeded drifters of that release_doy × age_bin_days · 24
```

holds exactly; `024f` asserts it before writing. A consumer therefore reads
the release-pool size `N` off the store itself, and `surv / (N · age_bin_days
· 24)` per hex is the mean share of released drifters afloat there during the
bin. Per-hex maps drop the sentinel by joining on the key.

**Origin, not release hex.** Considered keying by `release_hex`; rejected, it
recreates the O(10⁸)-row 024 counts store. The subbasin axis has 18 slots
(17 named + `-1`), so the dense bincount stays ~15 M cells per sidecar.

Each partition carries its build parameters as `df.attrs` (restored by
`pd.read_parquet`): `builder`, `regime`, `release_year`, `release_month`,
`hex_radius`, `member`, `max_age_days`, `age_bin_days`, `output_dt_hours`,
`w_c`, `tau_strong_hours`, `tau_calm_days`, `delta`, `trap_flat`,
`trap_wall`, `band_m`, `git_sha`, `slurm_job_id` — the same set as
[beaching.md](beaching.md) § Store schema minus `disp_bin_km`.

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
near-shore residence coincides with onshore waves **over flat shore** — at the
production `trap_flat = 1`, `trap_wall = 0` the depletion pattern is a
wave-exposure field multiplied by the substrate's flat fraction, so material
routed past exposed rock coast survives where the substrate-inert member
(`step_wc0p1_ts3_tcinf`) removes it. Absolute rates are parameter-sensitive ([beaching.md](beaching.md)); read the *pattern*.

## Cross-references

- [beaching.md](beaching.md) — the shared sidecar, rate model and the endpoint
  (stranding) counterpart store.
- [hexbinning_and_connectivity.md](hexbinning_and_connectivity.md) — the plain
  occupancy (024 counts) this re-weights, and the key schema.
- [visualisations.md](visualisations.md) — `030`'s plot rationale.
