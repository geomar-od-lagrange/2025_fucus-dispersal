# Hex-aggregated dispersal store

`notebooks/024_BuildHexAggregates.md` produces one store per
`(hex_radius)` build:

```
output_root/HexAggregates/r<hex_radius>m/
  key.parquet
  counts/regime=<regime>/release_year=<year>/part.parquet
```

Compact, query-friendly sister of the multi-TB raw zarrs; substrate
behind every map in `notebooks/025_HexHeatmaps.md`. Raw zarrs stay
the source of truth for per-trajectory diagnostics.

## Counts schema

| column        | type   | meaning                                              |
|---------------|--------|------------------------------------------------------|
| `release_hex` | int32  | hex containing the particle's release point          |
| `release_doy` | int16  | release day-of-year (one per partition)              |
| `age_bin`     | int8   | floor(particle age / `age_bin_days`) — default 10 d  |
| `target_hex`  | int32  | hex containing the particle position at this obs     |
| `n_obs`       | int32  | number of `(trajectory, obs)` pairs in this bin      |

`n_obs` is the residence proxy. For a "first-arrival" alternative,
dedupe `(trajectory_id, target_hex)` at query time against the raw
zarrs. Subbasin is **not** a grouping dim — recover via the key file.

## Key file (`key.parquet`, geoparquet, one row per hex)

| column            | type             | meaning                                        |
|-------------------|------------------|------------------------------------------------|
| `hex_id`          | int32 (PK)       | matches `release_hex` / `target_hex`           |
| `geometry`        | geoparquet       | hex polygon in EPSG:4326                       |
| `area_m2`         | float32          | full hex area (m², EPSG:3035)                  |
| `water_area_m2`   | float32          | wet intersection (m², EPSG:3035)               |
| `mean_depth_m`    | float32          | mean over `H0 > 0` cells inside the hex        |
| `dist_to_coast_m` | float32          | centroid → always-wet coast (m, EPSG:3035)    |
| `fucus_area_m2`   | float32          | intersection with REDLIST `F_vesiculo != 0`    |
| `helcom_subbasin` | int8 (nullable)  | HELCOM level-2 ID by centroid; -1 = outside    |

Schema is additive — new attributes mean rebuilding `key.parquet` with
a bumped version field; counts partitions don't move.

## Unified hex grid

Release and target hexes share one `hex_id` space — clicking a hex as
source ("where do particles go?") is the same operation as clicking
it as sink ("where do they come from?"), with `release_hex` and
`target_hex` swapped on a self-join. Considered separate
release/target grids; rejected — id translation kills self-joins.

```python
# Source: where do particles from H go?
counts.query("release_hex == H").groupby("target_hex")["n_obs"].sum()

# Sink: where do particles arriving at H come from?
counts.query("target_hex == H").groupby("release_hex")["n_obs"].sum()
```

## Equal-area projection

Hexes live in Lambert Azimuthal Equal-Area centred on the coarse-grid
BSH bounding box midpoint (computed from H0 at build time). LAEA
keeps hex areas comparable basin-wide; Mercator was rejected (4× area
variation), Geodetic was rejected (hexes wouldn't be hexagons). The
HexProj configuration travels as parquet file-level metadata on
`key.parquet` and every counts partition — without it `hex_id` values
are opaque, with it geometries are reproducible via `hextraj`.

`hex_radius` is the production parameter; multiple radii come from
papermill sweeps of 024.

## Domain coverage

The key file is built over the full BSH model domain (fine + coarse,
including land), not the wet polygon. A trajectory cannot escape the
BSH lat/lon extent, so every `hex_id` `hp.label` can ever assign is
pre-populated. Wet vs dry is captured by `water_area_m2` (zero for
land) and `mean_depth_m` (NaN for land); the key never filters them
out. Counts can therefore assert a hard key-completeness invariant
(every `release_hex`/`target_hex` is in the key) as a one-line
set-difference check.

## Source data and provenance

Key file derives from:

- BSH H0 statics (`data/bsh_hbmnoku_static/static_file_{fine,coarse}/H0_file_*.nc`)
- BSH wet-cell coastline geojsons under `data/bsh_hbmnoku_static/`
  (produced by `notebooks/004_extract_coastline.py`)
- HELCOM subbasins shapefile (`HELCOM_subbasins_2022_level2.shp`)
- REDLIST_SIS_Macrophytes shapefile (`F_vesiculo != 0`)

`mean_depth_m` filters to `H0 > 0` (see
[h0_semantics.md](h0_semantics.md)). Fine grid takes priority over
coarse for the depth join — fine cells contribute where covered;
coarse fills the rest, avoiding double-counting in the nest overlap.

Counts partitions carry self-identifying scope metadata (`regime`,
`release_year`, `release_doy`, `source_zarr`) so a lost partition
rebuilds from its metadata alone.

## Cross-references

- [seeding.md](seeding.md) — release set the `release_hex` derives from.
- [h0_semantics.md](h0_semantics.md) — `mean_depth_m` filter.
- [2d_field_extraction.md](2d_field_extraction.md) — coastline geojsons.
- [visualisations.md](visualisations.md) — notebook 025.
