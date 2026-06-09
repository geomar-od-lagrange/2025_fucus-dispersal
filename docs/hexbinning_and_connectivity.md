# Hex-aggregated dispersal store

The key file is shared by every counts/distance file at the same
`hex_radius`. Two aggregations write per-`(regime, year)` partitions:
024 occupancy counts and 024b distance histograms.

```
output_root/HexAggregates/
  HexAgg_key_r<radius>m.parquet         (notebook 024a — once per radius)
  HexAgg_key_r<radius>m.json            ( ditto, sidecar metadata     )
  HexAgg_counts_r<radius>m_<regime>_<year>.parquet
                                        (notebook 024  — once per (regime, year))
  HexAgg_distance_r<radius>m_<regime>_<year>.parquet
                                        (notebook 024b — once per (regime, year))
```

Flat layout deliberately: parallel 024/024b jobs (multiple regimes/years)
write disjoint filenames, so they can't race. Compact, query-friendly
sister of the multi-TB raw zarrs. The counts store is the substrate behind
`notebooks/025_HexHeatmaps.md` (density) and `notebooks/026_TimeHorizonMaps.md`
(density at a selected `age_bin`); the distance store is the substrate
behind `notebooks/027_HexDistanceQuantiles.md`. Raw zarrs stay the source
of truth for per-trajectory diagnostics.

## Counts schema

| column        | meaning                                                          |
|---------------|------------------------------------------------------------------|
| `release_hex` | hex containing the particle's release point; `-1` if land-seeded |
| `release_doy` | release day-of-year of the originating zarr                      |
| `age_bin`     | floor(particle age / `age_bin_days`) — default 10 d, no upper cap|
| `target_hex`  | hex containing the particle position at this obs; `-1` if NaN    |
| `n_obs`       | number of `(trajectory, obs)` pairs in this bin                  |

Dtypes are whatever pandas/dask pick by default (typically `int64`) —
parquet's RLE/dictionary encoding flattens the size cost, and counts
files are tiny next to the trajectory zarrs upstream.

`n_obs` is the residence proxy. For a "first-arrival" alternative,
dedupe `(trajectory_id, target_hex)` at query time against the raw
zarrs. Subbasin is **not** a grouping dim — recover via the key file.

`-1` is `hextraj.INVALID_HEX_ID`, surfaced when lon/lat is NaN — i.e.
land-seeded trajectories (zero first-step displacement) and any
out-of-domain positions. Preserved on disk so the land-seeded fraction
is queryable (`counts[counts.release_hex == -1].n_obs.sum()`); filter
with `release_hex >= 0` / `target_hex >= 0` at query time when only
valid hexes are wanted.

## Distance schema

One row per `(release_hex, release_doy, distance_bin)`: the distribution
of per-trajectory crow-flies **final displacement** (release point → last
valid position, equirectangular 111 km/deg metric) binned per source hex.

| column         | meaning                                                         |
|----------------|-----------------------------------------------------------------|
| `release_hex`  | hex containing the particle's release point (≥ 0; land dropped) |
| `release_doy`  | release day-of-year of the originating zarr                     |
| `distance_bin` | floor(final displacement / `distance_bin_km`) — default 1 km, no cap |
| `n_traj`       | number of trajectories in this bin                              |

Unlike counts, land-seeded particles are **dropped** here (a displacement
from a land seed is meaningless), so `release_hex` carries no `-1`.

A histogram, not pre-computed quantiles, because histograms are
**additive**: pooling Aug/Sep across years is summing partitions, then
deriving quantiles from the pooled cumulative count (027). Per-year
quantiles cannot be correctly averaged across years. `distance_bin_km` is
the build↔consumer contract (like `hex_radius` / `age_bin_days`); it is
not stored, so 027 must read with the same value.

## Key file (geoparquet, one row per hex)

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

Schema is additive — new attributes mean rebuilding the key with a
bumped version field; counts files don't move.

## Sidecar metadata

`HexAgg_key_r<radius>m.json` next to the parquet carries:

- `hex_proj` — `projection_name`, `lon_origin`, `lat_origin`,
  `hex_size_meters`. `HexProj(**meta["hex_proj"])` rematerialises
  geometry from any `hex_id`.
- `area_crs` — `"EPSG:3035"`, the projection used for all area columns.
- `subbasin_id_to_name` — int → name lookup for the `helcom_subbasin`
  column. Keys are JSON strings (cast to int on load).

Sidecar instead of parquet schema metadata so the layer that writes
the parquet (`GeoDataFrame.to_parquet`, `dask.dataframe.to_parquet`)
doesn't need to grow custom-metadata plumbing. Counts files carry no
metadata: the filename encodes `hex_radius`, `regime`, `release_year`,
and the projection comes from the matching key file.

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
projection lives in the key sidecar — without it `hex_id` values are
opaque, with it geometries are reproducible via `hextraj`.

`hex_radius` is the production parameter; multiple radii come from
papermill sweeps of 024a (key) and 024 (counts).

## Domain coverage

Key hexes = `hp.rectangle_of_hexes(BSH coarse bbox + 1 coarse cell)`.
Sampling H0 grid points directly under-counts because trajectories
use bilinear-interpolated currents and can land in hexes between
points (coastal pockets, coarse-grid periphery). `mean_depth_m` is
attached only where ≥1 wet H0 cell falls inside (fine priority via
`combine_first`), NaN elsewhere — so `hexes with H0` ⊆ `key`.

Out-of-key counts are warned but kept: particles can drift past the
bbox+margin (Stokes, numerical excursion). The build does not abort.

## Coastline merge (fine-precedes-coarse)

`water_area_m2` and `dist_to_coast_m` come from per-hex intersections
with `fine_polys ∪ (coarse_polys \ fine_bbox)`, where `fine_bbox` is
the H0-fine bbox padded by half a fine cell. Coarse is clipped to
coarse-only territory before union, making the two sides spatially
disjoint — fine has sole authority inside its bbox, coarse fills the
rest. A plain `unary_union` would wet pixels fine resolves as land.

## Source data and provenance

Key file derives from:

- BSH H0 statics (`data/bsh_hbmnoku_static/static_file_{fine,coarse}/H0_file_*.nc`)
- BSH wet-cell coastline geojsons under `data/bsh_hbmnoku_static/`
  (produced by `notebooks/004_extract_coastline.py`)
- HELCOM subbasins shapefile (`HELCOM_subbasins_2022_level2.shp`)
- REDLIST_SIS_Macrophytes shapefile (`F_vesiculo != 0`)

`mean_depth_m` filters to `H0 > 0` (see [h0_semantics.md](h0_semantics.md));
fine grid takes priority via `combine_first`, with coarse filling
where fine has no coverage.

Counts files self-identify via filename. Per-zarr provenance lives in
the `release_doy` column, not the metadata: a counts file aggregates
whatever zarrs were on disk at build time, so NESH job failures cause
coverage gaps (visible as missing `release_doy` values) rather than
aborting the build. Same `(release_time, regime)` zarrs from
independent reseeded reruns are summed additively by the groupby.

## Cross-references

- [seeding.md](seeding.md) — release set the `release_hex` derives from.
- [h0_semantics.md](h0_semantics.md) — `mean_depth_m` filter.
- [2d_field_extraction.md](2d_field_extraction.md) — coastline geojsons.
- [visualisations.md](visualisations.md) — notebooks 025–027.
