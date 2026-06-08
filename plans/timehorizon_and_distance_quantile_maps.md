# Time-horizon dispersal maps + hex0 distance-quantile maps

Source intent: [analysis.md](analysis.md) ("Dispersal maps" and
"Distance statistics" sections). Two deliverables:

1. **Time-horizon dispersal maps** — where the dispersal cloud sits at
   successive elapsed times (20, 40, 80, 160 d), for August/September
   releases pooled across all available years, per regime.
2. **hex0 distance-quantile maps** — for every source hex (hex0), the
   crow-flies final displacement from release for all its trajectories
   (Aug/Sep, pooled across years), reduced to per-hex0 quantiles
   (0.1, 0.5, 0.9) and drawn as one map per quantile level, per regime.

## Architecture: aggregate once, plot cheap

The heavy trajectory-zarr reads live in **024-style aggregation**
notebooks; every map notebook is a lightweight parquet-only consumer of a
hex store (no Dask cluster). This mirrors the existing
[024a](../notebooks/024a_BuildHexKey.md) (static key) +
[024](../notebooks/024_BuildHexAggregates.md) (counts) + [025](../notebooks/025_HexHeatmaps.md)
(plot) split.

| Quantity | Aggregation (heavy, dask, per regime/year) | Map (lightweight, parquet-only) |
|----------|---------------------------------------------|---------------------------------|
| Density  | 024 → counts `(release_hex, age_bin, target_hex, release_doy)→n_obs` | 025 (sum all ages) · **026 (select `age_bin`=horizon)** |
| Distance | **024b → distance histogram `(release_hex, release_doy, distance_bin)→n_traj`** | **027 (derive quantiles)** |

Both 026 and 027 read the static key from 024a for hex geometry. They
pool across years by globbing every `*_{regime}_*.parquet` partition.

Why a histogram store for distance (not pre-baked quantiles): histograms
are **additive**, so pooling Aug/Sep across years is summing partitions
then deriving quantiles from the pooled histogram. Per-year quantiles
cannot be correctly averaged across years. This is the same "store raw
counts, derive downstream" choice 024 makes for density, and follows
AGENTS.md's "trust derivations; let downstream raise."

Per the notebook-local-utility rule, the map notebooks **copy** 025's
`to_hex_gdf` / `log_density_plot` idioms inline. This reaches three
callsites (025/026/027); if they start to drift, revisit extraction —
but the repo has no shared notebook module today, so copy is the default.

---

## Notebook 024b — build hex distance histogram

**New notebook `notebooks/024b_BuildHexDistance.md` (recommended).**
Structural sibling of 024: same zarr glob, `parse_zarr_stem`, HexProj
from the 024a sidecar, Dask cluster, per-`(regime, release_year)` run.
024 aggregates per-obs occupancy; 024b aggregates a per-trajectory
scalar (final displacement) into a distance histogram. Separate notebook
(not a 024 extension) keeps one concern per aggregation — occupancy vs.
displacement are different reductions over different axes.

### Parameters cell

```python tags=["parameters"]
output_root = "../output"            # str: zarr read + store write root
regime = "surface"                   # str (papermill-swept)
release_year = 2019                  # int (papermill-swept)
hex_radius = 6000                    # int: must match a 024a key
distance_bin_km = 1.0                # float: histogram bin width
```

### Cell-by-cell outline

1. **md** — title + description (mirror 024's intro, "occupancy" →
   "final-displacement histogram").
2. **code (imports)** — copy 024's import block
   ([024:30-44](../notebooks/024_BuildHexAggregates.md)); add `numpy`.
3. **code (`parse_zarr_stem`)** — copy verbatim from
   [024:46-63](../notebooks/024_BuildHexAggregates.md).
4. **code (`release_hex_id`)** — NaN-safe `hp.label` wrapper (the idiom
   of 023's `release_subbasin`): returns -1 for NaN release positions,
   labels valid ones, lazy per chunk via `apply_ufunc`.
5. **md + code (parameters cell)** — as above.
6. **code (derived layout + key/HexProj)** — copy
   [024:87-108](../notebooks/024_BuildHexAggregates.md), swapping the
   counts path for `HexAgg_distance_r{radius}m_{regime}_{year}.parquet`.
   `hp = HexProj(**meta["hex_proj"])` from the sidecar.
7. **code (Dask cluster)** — copy verbatim from
   [024:112-123](../notebooks/024_BuildHexAggregates.md).
8. **code (zarr glob + parse)** — copy
   [024:128-143](../notebooks/024_BuildHexAggregates.md).
9. **code (`zarr_to_distance_frame`)** — **new core logic.** Per zarr:
   `on_land` mask, release `lon0`/`lat0` (obs=0), per-trajectory final
   displacement `distance_km(ds).ffill("obs").isel(obs=-1)` (equirect.
   111 km metric, [docs/distance_calculation.md](../docs/distance_calculation.md)),
   `release_hex = release_hex_id(lon0.where(~on_land), …)`,
   `distance_bin = (final_dist // distance_bin_km)`. Build a
   `(trajectory,)` `Dataset`, `.to_dask_dataframe()`, drop NaN
   `distance_bin` / `release_hex == -1`, tag `release_doy`.
10. **code (concat → groupby → size)** — `dd.concat([...]).groupby(
    ["release_hex","release_doy","distance_bin"]).size().rename("n_traj")
    .reset_index().compute()`, mirroring
    [024:165-173](../notebooks/024_BuildHexAggregates.md).
11. **code (write parquet)** — `to_parquet`, mirror
    [024:177-178](../notebooks/024_BuildHexAggregates.md).
12. **code (validation)** — rows, `sum(n_traj)`, release_doy span,
    distance_bin span, `#` release hexes. Computed dynamically.

### Do-not-improvise

- **`distance_bin` has no upper cap** — `floor(distance/bin)`, bins
  emerge from the data exactly like 024's `age_bin`
  ([024:160](../notebooks/024_BuildHexAggregates.md)).
- **Final displacement = `ffill("obs").isel(obs=-1)`** — last valid
  position's crow-flies distance; exited/NaN-tail trajectories carry
  their last in-domain distance forward.
- **`distance_bin_km` is the store↔consumer contract** — 027 must read
  with the same value (same pattern as `hex_radius` / `output_dt_mins`).

---

## Notebook 026 — time-horizon dispersal maps

**Rewrite `notebooks/026_TimeHorizonMaps.md` as a counts-store consumer.**
Reads 024a key + pooled 024 counts partitions; no zarrs, no cluster. ≈
025's per-quarter panel (Panel D) but faceted by elapsed-time horizon and
filtered to Aug/Sep releases.

### Parameters cell

```python tags=["parameters"]
data_root = "../data"                # str: HELCOM polygons (coast context)
output_root = "../output"            # str: hex store root
regime = "surface"                   # str (papermill-swept)
hex_radius = 6000                    # int: must match the store
age_bin_days = 10                    # int: must match 024
release_months_csv = "8,9"           # str ("" ⇒ all releases)
time_horizons_days_csv = "20,40,80,160"  # str
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66
cmap = "viridis"                     # log-density, justified like 025
baltic_panel_height_in = 6           # int
```

### Cell-by-cell outline

1. **md** — title + description (hex density at successive horizons).
2. **code (imports)** — copy 025's import block
   ([025:31-41](../notebooks/025_HexHeatmaps.md)) + `re` for the
   partition-year parse.
3. **code (parse params + Paths)** — parse the two `*_csv` lists; cast
   Paths.
4. **md + code (parameters cell)** — as above.
5. **code (read key + pool counts)** — read the key
   ([025:90-101](../notebooks/025_HexHeatmaps.md)); glob every
   `HexAgg_counts_r{radius}m_{regime}_*.parquet`, parse `release_year`
   from each filename, derive `release_month` from
   `(year, release_doy)` via `pd.to_datetime(year*1000+doy, "%Y%j").dt.month`
   ([025:289-293](../notebooks/025_HexHeatmaps.md) extended to month),
   filter to `release_months`, concat. Pools across years.
6. **code (`to_hex_gdf` + `log_density_plot`)** — copy from
   [025:128-183](../notebooks/025_HexHeatmaps.md). Copy the Natural Earth
   coast clip + extent/aspect block
   ([025:186-213](../notebooks/025_HexHeatmaps.md), Baltic only; drop the
   subbasin overlay — not used here).
7. **code (horizon → age_bin + panel grid)** — for each horizon,
   `age_bin = horizon_days // age_bin_days`; select
   `counts[counts.age_bin == age_bin]`, `to_hex_gdf`, `log_density_plot`
   into one subplot per horizon (grid like
   [025:296-309](../notebooks/025_HexHeatmaps.md), titled `"{h} d"`).
8. **md + code (validation)** — regime, months, horizons→age_bins,
   per-horizon `sum(n_obs)`. Computed dynamically.

### Do-not-improvise

- **Horizon → `age_bin = horizon_days // age_bin_days`.** Each map is the
  10-day occupancy window `[horizon, horizon+age_bin_days)`. Horizons must
  be multiples of `age_bin_days` (20/40/80/160 with bins of 10 are) — else
  the floor lands the map in a neighbouring window; assert and print.
- **Month filter needs the year** for leap-correct doy→month — parse it
  from each partition filename, never assume one year.

---

## Notebook 027 — hex0 distance-quantile maps

**Rewrite `notebooks/027_HexDistanceQuantiles.md` as a distance-store
consumer.** Reads 024a key + pooled 024b distance partitions; no zarrs,
no cluster.

### Parameters cell

```python tags=["parameters"]
data_root = "../data"                # str: coast context
output_root = "../output"            # str: hex store root
regime = "surface"                   # str (papermill-swept)
hex_radius = 6000                    # int: must match the store
distance_bin_km = 1.0                # float: must match 024b
release_months_csv = "8,9"           # str ("" ⇒ all releases)
quantile_levels_csv = "0.1,0.5,0.9"  # str
min_traj_per_hex = 30                # int
baltic_lon_min, baltic_lon_max = 5, 32
baltic_lat_min, baltic_lat_max = 53, 66
baltic_panel_height_in = 6           # int
```

### Cell-by-cell outline

1. **md** — title + description.
2. **code (imports)** — 025's block ([025:31-41](../notebooks/025_HexHeatmaps.md))
   + `re`.
3. **code (parse params + Paths)** — parse `quantile_levels`,
   `release_months`; cast Paths.
4. **md + code (parameters cell)** — as above.
5. **code (read key + pool distance store)** — read key; glob every
   `HexAgg_distance_r{radius}m_{regime}_*.parquet`, parse year, derive
   `release_month`, filter to `release_months`, concat, then
   `groupby(["release_hex","distance_bin"]).n_traj.sum()` — the pooled
   per-hex histogram.
6. **code (per-hex quantiles from the histogram)** — **new core logic.**
   Per `release_hex`, sort by `distance_bin`, `cumsum` `n_traj`, and read
   off each quantile via `np.searchsorted(cum, q*total)` → bin-left-edge
   `distance_bin * distance_bin_km` km. Apply `min_traj_per_hex` (drop
   hexes whose `total < min`). Pure pandas groupby over a few hundred
   hexes — instant, no dask.
7. **code (`to_hex_value_gdf` + `hex_value_plot` + coast)** — adapt
   025's `to_hex_gdf` (merge value series onto key geometry) and
   `log_density_plot` ([025:128-183](../notebooks/025_HexHeatmaps.md)),
   renamed and **without** `np.log10`/`cmap="viridis"`: distance is
   linear → default colormap, matplotlib auto-ranged, colorbar in km.
8. **code (one map per quantile)** — loop `quantile_levels`, one
   Baltic-extent map each.
9. **md + code (validation)** — regime, months, quantile levels, `#`
   source hexes meeting the gate, min/median/max per band.

### Do-not-improvise

- **Quantiles derived from the pooled histogram** (cumulative count over
  `distance_bin`), never averaged from per-year quantiles. Bin-resolution
  approximate — fine at 1 km.
- **`distance_bin_km` must match 024b** to convert bins back to km.
- **Plotting:** default colormap + matplotlib auto-range + km colorbar
  (`legend_kwds={"label": ...}` — geopandas derives no unit label). Carry
  025's layout/registration overrides (aspect figsize, hex
  `edgecolor="face"/linewidth=0.4`, black coastline) per
  [docs/visualisations.md](../docs/visualisations.md); do **not** carry
  `cmap="viridis"` or `np.log10`.

---

## Resolved decisions

- Time horizons 20/40/80/160 d; quantiles 0.1/0.5/0.9; distance = final
  displacement; `min_traj_per_hex` 30; Aug/Sep pooled across all years;
  bottom seeds the same Fucus cells as surface (hex0 sets identical →
  regimes comparable); figure sizing follows 025.
- 026 reads the 024 counts store; 027 reads a new 024b distance
  histogram store; both are lightweight consumers.
- 160 d fits the 220 d runs.

## Open questions / decisions for the human

- **`distance_bin_km`:** 1 km assumed — fine-grained enough for
  quantiles; bump if the store size matters.
- **Horizon temporal window:** tied to `age_bin_days`=10, so each 026 map
  is a 10-day occupancy window. If tighter snapshots are wanted, rebuild
  024 at a smaller `age_bin_days`.
- **Helper extraction:** `to_hex_gdf`/`log_density_plot` now copied in
  025/026/027. Copy stands as the AGENTS.md default; flag if they drift.

## Out of scope

The deferred work in [analysis.md](analysis.md) — the everywhere→
everywhere bottom transfer matrix `T` and the bottom-Ekman comparison —
is not covered here.

## Implementation / verification notes

- Build with the [jupytext skill](../.agents/skills/jupytext/SKILL.md):
  new notebooks are `py:percent,md,ipynb` with the `.py` authoritative;
  commit all three forms, the `.ipynb` code-only.
- 024b uses the multi-task Dask job-script pattern of 023/024; 026/027 use
  the lightweight parquet-only pattern of 025.
- Test-run via papermill `--cwd notebooks/`, sweeping `regime` (and
  `release_year` for 024b). Build 024b before 027; build 024 before 026.
- Lighter model for mechanical authoring; review after.
- When implemented: docs (extend [hexbinning_and_connectivity.md](../docs/hexbinning_and_connectivity.md)
  for 024b; the 026/027 entries already live in
  [docs/visualisations.md](../docs/visualisations.md)), move this plan to
  `plans/done/`.
