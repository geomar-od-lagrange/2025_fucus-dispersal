# Subbasin→subbasin connectivity store

## Purpose

A HELCOM subbasin-resolved connectivity aggregate: for each
`(origin_subbasin, target_subbasin)` pair, how much particle-time
released from the origin's start hexes lands in the target's hexes.
Sits beside the per-hex counts/distance stores as a third
`(regime, year)` partition, keeping the same params so later notebooks
pool it (all-year, age horizons) the same way they pool counts.

## Semantics: residence, re-aggregated from 024 counts

Connectivity here is **residence** — particle-timesteps, not particle
flux. It is a pure re-partition of the 024 counts store
(`(release_hex, release_doy, age_bin, target_hex) → n_obs`):

```python
origin_sb = release_hex → key.helcom_subbasin   (NaN/-1 → -1)
target_sb = target_hex  → key.helcom_subbasin   (NaN/-1 → -1)
conn = counts.groupby([origin_sb, target_sb, release_doy, age_bin])["n_obs"].sum()
```

`n_obs` counts `(trajectory, obs)` pairs, so summed over a subbasin's
hexes it is the time particles from the origin spend in the target. It
is **fully additive** across `age_bin`, `release_doy`, and year — the
later all-year / age-horizon aggregates are just `groupby().sum()`,
exactly as for counts.

Considered, rejected: **particle flux** (count each trajectory once per
target subbasin, by first-arrival age). That is the standard ecological
connectivity matrix and answers "how many particles reached B", but it
cannot be derived from the counts store — `n_obs` has collapsed
trajectory identity, so deduping would require a fresh Dask
aggregation off the zarrs. Deferred; residence is the agreed metric and
is consistent with the 025/026 density semantics.

## No double counting — by construction

Each hex has exactly one centroid subbasin (the key's
`helcom_subbasin`, assigned by centroid — boundary hexes go to one
side, never split by polygon area). So every counts row maps to exactly
one `(origin_sb, target_sb)`: the groupby re-partitions `n_obs` with no
overlap. The invariant is checkable in one line and is the validation
gate:

```python
assert conn["n_obs"].sum() == counts["n_obs"].sum()
```

## Unnamed / outside handling — single `-1` sentinel

The key carries two unnamed states: land-seed / NaN positions
(`hex == -1`, absent from the key → lookup yields NaN) and
in-key-but-outside-every-polygon (`helcom_subbasin == -1`). At subbasin
granularity these collapse to one category, so both fold to `-1` on
each side. `origin_subbasin`/`target_subbasin ∈ {named ids ≥ 0} ∪ {-1}`.
Kept on disk (queryable), filtered with `>= 0` at read time when only
named subbasins are wanted — same convention as the counts store's `-1`.

## Layout & schema

New build `notebooks/024c_BuildHexConnectivity.{py,md,ipynb}`. It is a
024x-family store producer keyed by `(regime, release_year)` with the
same params as 024 (`hex_radius`, `age_bin_days`). Unlike 024/024b it
reads the **counts parquet**, not the zarrs — so no Dask cluster, a
lightweight parquet-only notebook like 025–027.

```
output_root/HexAggregates/
  HexAgg_connectivity_r<radius>m_<regime>_<year>.parquet
```

| column            | meaning                                                   |
|-------------------|-----------------------------------------------------------|
| `origin_subbasin` | HELCOM subbasin id of the release hex; `-1` = unnamed     |
| `target_subbasin` | HELCOM subbasin id of the position hex; `-1` = unnamed    |
| `release_doy`     | release day-of-year of the originating zarr (kept axis)   |
| `age_bin`         | `floor(age / age_bin_days)` (kept axis)                   |
| `n_obs`           | summed `(trajectory, obs)` particle-timesteps             |

Subbasin id→name lives in the existing key sidecar
(`HexAgg_key_r<radius>m.json`, `subbasin_id_to_name`); the connectivity
file carries no metadata — filename encodes radius/regime/year, names
come from the key, same as counts.

## Params (parameters cell)

```python
output_root = "../output"
regime = "surface"
release_year = 2019
hex_radius = 6000
age_bin_days = 10   # must match the counts file being read
```

`age_bin_days` is the build↔counts contract (counts already bins at this
granularity); 024c does not re-bin, it groups the existing `age_bin`.
`output_dt_mins` is not needed — `age_bin` is already materialised in
counts.

## Build outline

1. Resolve layout; require the counts partition + key + sidecar to
   exist (hard prerequisite, like 024 requires the key).
2. Read counts parquet + key (`hex_id`, `helcom_subbasin`).
3. Map `release_hex`/`target_hex` → subbasin via the key, fill NaN with
   `-1`, cast to int.
4. `groupby([origin_subbasin, target_subbasin, release_doy, age_bin])
   ["n_obs"].sum()`.
5. Write parquet.
6. Validation cell: assert `n_obs` sum preserved; print matrix shape,
   #named subbasins present, fraction of `n_obs` in `-1` origin/target,
   diagonal (within-subbasin) share.

## Job script

`scripts/024c_*_job.sh` papermill-swept across `(regime, year)`, mirroring
024/024b. Cheap (parquet-only), so a single small step — no Dask cluster
bootstrap.

## Visualisation: 028 connectivity matrix (POC)

`notebooks/028_SubbasinConnectivityMatrix.{py,md,ipynb}` — parquet-only
consumer of the 024c store, no Dask, in the 025–027 lineage. A
proof-of-concept: print the matrix and draw the heatmap; richer views
(row-normalised emission fractions, age-horizon matrices, chord
diagrams) come later.

### Scopes

Two release-month windows in one run, mirroring 026's pooling:

- **all-year** — every release month.
- **aug_sep** — months 8, 9.

Both pooled across all available years (glob `*_<regime>_*.parquet`,
parse year from filename for leap-correct `release_doy → month`, exactly
as 026a). Encoded as a derived `scopes = {"all_year": [], "aug_sep":
[8, 9]}` cell, not in the parameters cell — the parameters cell stays
primitive (`output_root`, `regime`, `hex_radius`, `age_bin_days`).

### Per scope

1. Filter to the scope's release months; **pool over `release_doy` and
   `age_bin`** (POC ignores elapsed time — sums all ages) →
   `(origin_subbasin, target_subbasin) → n_obs`.
2. Restrict to named subbasins (`origin >= 0 & target >= 0`); report the
   `-1` (unnamed origin/target) `n_obs` fraction dropped, as 026a does.
3. Square matrix over the **union** of named ids present in the pooled
   data, ordered by id (stable). Pivot to a pandas DataFrame with
   subbasin **names** as index (origin, rows) and columns (target).
   Origins are a subset — only seeding subbasins have non-zero rows; a
   square frame keeps source/sink symmetry readable.
4. **Print** the DataFrame (the matrix).
5. **Heatmap** twice — linear (default norm) and log (`LogNorm`, zeros
   masked to NaN so they render blank; the diagonal retention term
   dominates by orders of magnitude, which is exactly why both scales
   are shown). `imshow` on the matrix; subbasin names as tick labels set
   **at plot time** (don't mutate coords). Default colormap (imshow's
   `viridis`) — no `cmap=`/`figsize=` override; `norm=LogNorm(...)` is
   the one non-default, required for the log panel and noted inline.
6. Save PNGs under `output_root/Figures/028`
   (`SubbasinConnectivityMatrix_<regime>_r<radius>m_<scope>_<lin|log>.png`).

### Params (parameters cell)

```python
output_root = "../output"
regime = "surface"
hex_radius = 6000
age_bin_days = 10   # store contract; pooled away here but kept explicit
```

`age_bin_days` is unused arithmetically (we sum all ages) but kept so the
notebook reads the store under the same contract as its siblings; drop
it only if a no-op param reads as noise.

## Later (out of scope)

- Age-horizon matrices: cumulative `n_obs` over `age_bin ≤
  T/age_bin_days` (the connectivity analogue of 026's horizons).
- Row-normalised matrix: `n_obs / row.sum()` = fraction of an origin's
  particle-time reaching each target (the comparable connectivity view;
  raw `n_obs` rows differ by orders of magnitude across origins).
- Flux (deduped-trajectory) connectivity from the zarrs, if residence
  proves insufficient (see "Semantics" above).

## Docs follow-up

On implementation, add a "Connectivity" section to
[../docs/hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)
(schema table + the residence-vs-flux note + the sum-preserved
invariant) and a 028 entry to
[../docs/visualisations.md](../docs/visualisations.md); move this plan to
`plans/done/`.
