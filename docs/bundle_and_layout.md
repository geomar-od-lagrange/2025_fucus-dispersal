# Data bundle, path layout, and output locations

The pipeline separates inputs, outputs, and the heavy on-NESH BSH store
into three explicit roots passed as papermill parameters. There are no
`FUCUS_OUTPUT_ROOT` or similar environment variables: the papermill
parameter is the single contract, visible in each notebook's parameters
cell.

## The three roots

| Parameter | Default (relative to `notebooks/`) | Holds |
|---|---|---|
| `data_root` | `"../data"` | Redistributed inputs from the data twin submodule. |
| `output_root` | `"../output"` | Heavy outputs, produced by the pipeline. |
| `bsh_root` | `"../data/bsh_hbmnoku_static"` | BSH HBMnoku store — the demo default points at the bundled static subset (statics + coastlines); NESH runs override to the full multi-year store. |

Each notebook's parameters cell holds only the roots it consumes. Every
job script passes them on the papermill command line, resolved to
absolute NESH paths.

## The data twin submodule

`./data/` is a git submodule pointing at
<https://git.geomar.de/od-lagrange/2025_fucus_dispersal_data>, with
top-level dirs named `<source>_<dataset>` and derived blobs co-located
with their source:

```
data/
  helcom_subbasins_2022/      # HELCOM level-2 polygons
  helcom_fucus_redlist/       # REDLIST_SIS_Macrophytes.* + the derived
                              # fucus_release_points.geojson (baked by 000)
  bsh_hbmnoku_static/         # BSH static grids (lonlat, H0) + wet-cell
                              # coastline geojsons (derived via scripts/004)
  bsh_hbmnoku_demo/           # BSH HBMnoku demo subset: full-day c/h/t/z
  cmems_stokes_sample/        # 1-day CMEMS Stokes drift sample
  ATTRIBUTION.md              # per-dataset attribution blocks
  README.md
```

Binary blobs (shapefile parts, NetCDF) live in git-LFS. Plain text
(GeoJSON, metadata XML, attribution) stays in regular git.

## Clone and fetch

```
git clone --recurse-submodules https://github.com/geomar-od-lagrange/2025_fucus-dispersal.git
cd 2025_fucus-dispersal
git -C data lfs pull
```

`git-lfs` must be installed locally first (e.g. `pixi global install
git-lfs`), then run `git lfs install` once per user.

After a plain clone:

```
git submodule update --init data
git -C data lfs pull
```

## Rebuild inputs from upstream

`scripts/fetch_data.sh` reconstructs `./data/` from public upstream
sources. The recipe chain:

```
scripts/obtain/download_helcom_subbasins.sh     # maps.helcom.fi MADS GP service
scripts/obtain/download_fucus_shapefile.sh      # same mechanism (MADS)
scripts/obtain/download_stokes_sample.sh        # via copernicusmarine.subset
scripts/obtain/download_bsh_hbmnoku_demo.sh     # curl from BSH OpenData (public, no auth)
```

After the obtain chain, `fetch_data.sh` invokes notebook 000 via
jupytext to bake the Fucus release-points geojson alongside the
shapefile in `helcom_fucus_redlist/`. The same chain is what the twin's
CI runs to keep its LFS blobs in sync with the recipe.

All four obtain steps work from any networked machine — no NESH mount
required. The CMEMS step additionally needs `copernicusmarine`
credentials; the others are unauthenticated.

## NESH output layout

Cluster outputs live **outside** the repo tree, at
`<work>/2025_fucus_dispersal_outputs/`. Every NESH job script sets
`output_root` to this absolute path and passes it to papermill.

```
output_root/
  2d_fields/                        # from scripts/003
  sigma/                            # from scripts/001
  stokes/                           # from scripts/002 (full years)
  Trajectories/<regime>/<year>/*.zarr   # from notebook 010
  HexAggregates/<config_name>/          # from notebook 024
    key.parquet                         # one row per hex in the BSH-model domain
    counts/regime=<regime>/release_year=<year>/part.parquet
```

Two aggregate stores are produced: `baltic_r6km_v1` (basin reference,
6 km hex radius, Baltic origin) and `de_r4km_v1` (German-waters zoom,
4 km hex radius, DE origin). Notebook 025 renders heatmap panels from
both.

## Pipeline stages

Notebook numbering is linear, with 10-step gaps for inserts.

| Stage | File | Purpose |
|---|---|---|
| 000 | `notebooks/000_FucusStartLocations.md` | One-shot bake: Fucus shapefile → `data/helcom_fucus_redlist/fucus_release_points.geojson`. Driven by `scripts/fetch_data.sh`. |
| 001 | `scripts/001_prepare_sigma_files.py` | BSH sigma file cleanup. |
| 002 | `scripts/002_download_stokes.py` | Full-year Stokes drift download. |
| 003 | `scripts/003_prepare_2d_fields.py` | BSH → 2D surface + bottom current / T / S fields. |
| 004 | `scripts/004_extract_coastline.py` | Extract BSH wet-cell coastline from H0 (the bundled result lives in the twin's `bsh_hbmnoku_static/`). |
| 010 | `notebooks/010_FucusDispersal.md` | Parcels runs, swept over `(year, doy, regime, velocity_factor)` via papermill. |
| 020 | `notebooks/020_RawTrajectories.md` | Raw trajectory lines. |
| 021 | `notebooks/021_TimeStats.md` | Per-trajectory time stats. |
| 022 | `notebooks/022_DispersalDistance.md` | Dispersal distance vs. time. |
| 023 | `notebooks/023_Heatmaps.md` | Density + mean-age heatmaps on a lon/lat grid. |
| 024 | `notebooks/024_BuildHexAggregates.md` | Build the hex-aggregated dispersal stores (parquet key + counts partitions). Iterates regimes via discovery over `output_root/Trajectories/*/`. |
| 025 | `notebooks/025_HexHeatmaps.md` | Render hex heatmaps from the aggregate stores (no Dask; geopandas + parquet reads). |

## Cross-references

- [../plans/done/bundle_and_layout.md](../plans/done/bundle_and_layout.md)
  — the migration plan this layout came out of.
- [../ATTRIBUTION.md](../ATTRIBUTION.md) — per-dataset attribution
  requirements (repo-root copy; the twin carries a verbatim copy).
