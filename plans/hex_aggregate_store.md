# Hex-aggregated dispersal store

## Purpose

Distil the multi-TB raw trajectory zarrs into a compact, query-friendly
aggregate so targeted visualisations (per source, per sink, per season,
per age slice) are interactive — no re-scan of the zarrs, no per-plot
aggregation. Designed to support an `ostrea.geomar.de`-style source↔sink
explorer as the downstream use case: click any hex as origin, see the
dispersal kernel; click any hex as destination, see the source kernel.

Raw trajectories stay as the source of truth for per-trajectory diagnostics
(beaching, first-passage times, etc.). This store is strictly for
density/reach aggregates on a fixed spatial grid.

## Design decisions

### Grid: unified hex resolution for source and target

Release and target hexes share **the same grid and the same `hex_id`
space**. A release hex and a target hex are just two roles the same hex
plays in a query. This is what makes ostrea-style "click to flip
source/sink" queries a symmetric self-join on one table.

Hex projection: Lambert Azimuthal Equal-Area centred on the Baltic
(matches 024). Stored as parquet-level metadata on every file so
geometries can be regenerated from `hex_id` alone via `hextraj`.

The key covers the **full BSH model domain** (fine + coarse grids,
including land cells), not just the wet polygon. A trajectory can't
escape the BSH lat/lon extent, so every hex `hp.label` can assign is
pre-populated in the key. Wet vs. dry is captured by `water_area_m2`
(zero for land hexes) and `mean_depth_m` (NaN for land hexes). This
lets the counts store assert hard key-completeness: every `release_hex`
and `target_hex` is guaranteed in the key, and any violation is a data
integrity bug rather than a "just buffer the key" workaround.

**Open: what resolution.** 10 km works for Baltic-wide maps but is
coarse for German-waters detail; 4 km resolves the German Bight but
multiplies Baltic-wide storage. A single resolution has to serve both
maps. Decide after the test decomposition below.

### Counts store: schema

One row per combination of:

| dim             | cardinality          |
|-----------------|----------------------|
| `regime`        | 3 (surface / surface_stokes / bottom) |
| `release_year`  | 9                    |
| `release_doy`   | 73 (day-of-year, ignoring leap-year 74th) |
| `release_hex`   | ~number of release hexes (depends on resolution) |
| `age_bin`       | 22 (10-day bins up to 220 days) |
| `target_hex`    | ~number of visited hexes (depends on resolution) |

Values per row: `n_obs` (hit count, residence proxy). Subbasin is not
a grouping dim; a release hex maps uniquely to a HELCOM subbasin,
recovered via the key-file join.

Partition files by `(regime, release_year)` → 27 parquet files,
pushdown-friendly over the remaining dims.

Size expectation at 5 km (central estimate): ~20 GB parquet, ~80–100 GB
as a polars / pyarrow table.

### Key file: geoparquet, one row per hex

Same hex grid as the counts store. Built once; joined into queries
that need geometry or metadata.

Units are SI throughout — all areas in m², distances in m. Float32
carries the range without ambiguity about prefixes.

The HexProj configuration (projection name, origin lon/lat, hex size in
metres) is attached as parquet file-level metadata on both the key file
and every counts partition. Without it, `hex_id` values are opaque
integers; with it, geometries are reproducible from `hex_id` alone via
`hextraj`. The key file's `geometry` column is a convenience — the
metadata is the authoritative grid definition.

Columns:

- `hex_id` *(int32, PK)* — matches `release_hex` / `target_hex`
- `geometry` *(geoparquet)* — hexagon polygon in EPSG:4326
- `area_m2` *(float32)* — full hex area
- `water_area_m2` *(float32)* — water-filled intersection with the BSH land mask
- `mean_depth_m` *(float32)* — mean depth over always-wet cells inside the hex
- `dist_to_coast_m` *(float32)* — centroid distance to nearest BSH coast
- `fucus_area_m2` *(float32)* — intersection with REDLIST_SIS_Macrophytes `F_vesiculo != 0`
- `helcom_subbasin` *(int8, nullable)* — HELCOM level-2 ID by centroid

Centroids are not stored; derive them from `geometry` at query time via
`gdf.geometry.centroid`.

The list is expected to grow (beaching propensity, coastal-habitat
indicators, country-EEZ, etc.). Treat the schema as additive —
rebuild the key file, bump a version in its metadata, counts store
stays untouched.

Dataset provenance fields (schema version, build-time git SHA, build
timestamp) are planned for when this logic graduates into `hextraj`
proper, but are not implemented in the current build — the custom
metadata only carries `hex_proj`, `subbasin_id_to_name`, `area_crs`
on the key file and the per-partition scope fields (`regime`,
`release_year`, `release_doy`, `source_zarr`) on counts.

Source data: BSH H0 (depth, land mask), HELCOM subbasins shapefile,
REDLIST_SIS_Macrophytes shapefile, and a canonical BSH-model coastline
geojson produced by `scripts/004_extract_coastline.py` (drafts already
under `min_data/data/BSH_model_coastline/`). Finalising that geojson
is a prerequisite of this work — the coastline drives both
`dist_to_coast_m` and the wet-region decomposition used to pick the
resolution. All paths covered by `plans/portable_data_paths.md`.

`mean_depth_m` is computed over always-wet cells only. See
`plans/seafloor-location-H0-semantics.md` for H0 semantics — tidal-flat
cells can have negative H0 and must not be mixed into a depth average.

## Open questions

### Hex resolution — pre-build decision

Decompose the BSH wet region into hexes at candidate radii (4, 6,
10 km). For each:

- Render the tessellation over the whole Baltic at 12×12 inch figsize,
  and over German waters at 4×4 and 6×6 inch figsizes. Visual
  inspection picks which resolutions give usable detail in both scopes.
- Compute expected data volume assuming (a) all wet hexes are potential
  target hexes, (b) release hexes are the subset overlapping Fucus
  release cells. This replaces the earlier heuristic upper bound on
  the target-hex count with a direct count.

Then decide.

### Size-estimate constants — pin before committing

The wet-hex decomposition above pins the visited-target-hex count.
The remaining soft constant — the average number of distinct target
hexes touched per `(release_hex, release_doy, age_bin)` group — still
needs one real aggregation pass before committing to the full 9-year ×
3-regime build. Run one regime × one year at the chosen resolution,
aggregate, read that constant off the result.

### Fucus encoding

`fucus_area_m2` keeps everything derivable (`has_fucus = area > 0`,
fraction of hex, etc.) and matches the SI-units convention. If the
downstream viz turns out to only ever use a boolean, collapse then —
don't prematurely discretise.

## Out of scope (for now)

- Trajectory-level diagnostics (beaching, first-passage, age at
  stranding) — stays on the raw zarrs.
- Reach probability / per-trajectory distinct-count metrics — can live
  in a companion store later, recoverable from the raw zarrs.
- German-waters-only fine grid as a separate store — only needed if the
  unified-resolution decision turns out to be unworkable.
- Calendar-month / season at observation as a stored dim — derivable
  from `release_doy + age_bin_mid`.
