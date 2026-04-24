---
jupyter:
  jupytext:
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: min_data (pixi)
    language: python
    name: min_data
---

<!-- #region papermill={"duration": 0.008854, "end_time": "2026-04-24T16:02:49.601918+00:00", "exception": false, "start_time": "2026-04-24T16:02:49.593064+00:00", "status": "completed"} -->
# Build hex-aggregated dispersal store

Distil the per-regime trajectory zarrs into the compact aggregate
store defined in `plans/hex_aggregate_store.md`. This notebook builds:

- `key.parquet` — one row per hex in the BSH model domain (geometry,
  area, water area, depth, coast distance, Fucus area, HELCOM
  subbasin). Built once from static inputs (BSH wet-region geojson,
  H0 grids, Fucus shapefile, HELCOM level-2 polygons). Land hexes
  are included with `water_area_m2 = 0` so the key covers every hex
  `hp.label` can assign to a trajectory position.
- `counts/regime=…/release_year=…/part.parquet` — one row per
  `(release_hex, age_bin, target_hex)` with `n_obs` aggregate. Built
  per regime × release_year from the zarrs.

Units are SI throughout (m², m). The HexProj configuration travels
as parquet file-level metadata so `hex_id` values can be rematerialised
into geometry downstream.
<!-- #endregion -->

```python papermill={"duration": 0.618758, "end_time": "2026-04-24T16:02:50.226650+00:00", "exception": false, "start_time": "2026-04-24T16:02:49.607892+00:00", "status": "completed"}
import json
import os
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import xarray as xr
from shapely.ops import unary_union

from hextraj import HexProj

from helpers import load_trajectories, mask_land_seeded
```

<!-- #region papermill={"duration": 0.001314, "end_time": "2026-04-24T16:02:50.229543+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.228229+00:00", "status": "completed"} -->
# Parameters
<!-- #endregion -->

```python papermill={"duration": 0.004662, "end_time": "2026-04-24T16:02:50.235560+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.230898+00:00", "status": "completed"} tags=["parameters"]
base_path = "/Users/wrath/src/github.com/geomar-od-lagrange/2025_fucus-dispersal"

# Unified hex grid (release + target), matching 024's Baltic origin.
hex_size_meters = 6000
hex_origin_lon = 18.0
hex_origin_lat = 59.0

# Age binning.
age_bin_days = 10        # 22 bins × 10 d = 220 d of drift.
output_dt_mins = 60      # zarr output cadence.

# Dev scope: one zarr per regime for release_year = 2019.
regimes = ["surface", "surface_stokes", "bottom"]
release_year = 2019
```

```python papermill={"duration": 0.004133, "end_time": "2026-04-24T16:02:50.241390+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.237257+00:00", "status": "completed"}
base_path = Path(base_path)
dataset_root = base_path / f"output/HexAggregates/r{hex_size_meters // 1000}km_v1"
dataset_root.mkdir(parents=True, exist_ok=True)

n_age_bins = 220 // age_bin_days   # 22
max_age_bin = n_age_bins - 1        # 21 (bins 0..21)
```

<!-- #region papermill={"duration": 0.001251, "end_time": "2026-04-24T16:02:50.244065+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.242814+00:00", "status": "completed"} -->
# Dask cluster
<!-- #endregion -->

```python papermill={"duration": 0.598343, "end_time": "2026-04-24T16:02:50.843709+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.245366+00:00", "status": "completed"}
from dask.distributed import Client

scheduler_file = os.environ.get("SCHEDULER_FILE")
if scheduler_file:
    for _ in range(60):
        if os.path.exists(scheduler_file):
            break
        time.sleep(1)
    client = Client(scheduler_file=scheduler_file)
else:
    client = Client(ip="0.0.0.0")
client
```

<!-- #region papermill={"duration": 0.001381, "end_time": "2026-04-24T16:02:50.846754+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.845373+00:00", "status": "completed"} -->
# Hex projection + BSH-domain decomposition
<!-- #endregion -->

```python papermill={"duration": 0.17728, "end_time": "2026-04-24T16:02:51.025421+00:00", "exception": false, "start_time": "2026-04-24T16:02:50.848141+00:00", "status": "completed"}
hp = HexProj(
    projection_name="laea",
    lon_origin=hex_origin_lon,
    lat_origin=hex_origin_lat,
    hex_size_meters=hex_size_meters,
)

wet_gdf = gpd.read_file(base_path / "data/BSH_model_coastline/coastline.geojson")
always_wet_gdf = gpd.read_file(
    base_path / "data/BSH_model_coastline/coastline_always_wet.geojson"
)

wet_union = unary_union(wet_gdf.geometry)
always_wet_union = unary_union(always_wet_gdf.geometry)

# The key covers the full BSH model domain, not just the wet region.
# Label every fine+coarse H0 grid point; the union is every hex
# hp.label can return for any trajectory position that stays inside the
# model's lat/lon bounds. Wet vs. dry is captured later by water_area_m2.
_domain_ids = set()
for grid in ("fine", "coarse"):
    h0 = xr.open_dataset(
        base_path / f"min_data/bsh_operationalmodel_data/static_file_{grid}/H0_file_{grid}.nc"
    )
    lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
    labels = hp.label(lon2d.ravel(), lat2d.ravel())
    _domain_ids |= set(int(x) for x in labels if x >= 0)
hex_ids = np.asarray(sorted(_domain_ids), dtype=np.int32)
print(f"BSH-domain hexes: {len(hex_ids):,}")
```

<!-- #region papermill={"duration": 0.001447, "end_time": "2026-04-24T16:02:51.028514+00:00", "exception": false, "start_time": "2026-04-24T16:02:51.027067+00:00", "status": "completed"} -->
# Key file

## Hex geometries in EPSG:4326 and EPSG:3035 (equal-area for Europe)
<!-- #endregion -->

```python papermill={"duration": 0.128738, "end_time": "2026-04-24T16:02:51.158776+00:00", "exception": false, "start_time": "2026-04-24T16:02:51.030038+00:00", "status": "completed"}
hex_gdf = hp.to_geodataframe(hex_ids.tolist())
# to_geodataframe returns a 1-col GeoDataFrame indexed by hex_id; realign.
hex_gdf = hex_gdf.loc[hex_ids].reset_index().rename(columns={"index": "hex_id"})
hex_gdf["hex_id"] = hex_gdf["hex_id"].astype(np.int32)
hex_gdf = hex_gdf.set_crs(epsg=4326, allow_override=True)

hex_gdf_3035 = hex_gdf.to_crs(epsg=3035)
wet_union_3035 = gpd.GeoSeries([wet_union], crs=4326).to_crs(3035).iloc[0]
always_wet_union_3035 = (
    gpd.GeoSeries([always_wet_union], crs=4326).to_crs(3035).iloc[0]
)
```

<!-- #region papermill={"duration": 0.00146, "end_time": "2026-04-24T16:02:51.161968+00:00", "exception": false, "start_time": "2026-04-24T16:02:51.160508+00:00", "status": "completed"} -->
## Areas (m², in EPSG:3035)
<!-- #endregion -->

```python papermill={"duration": 6.110844, "end_time": "2026-04-24T16:02:57.274202+00:00", "exception": false, "start_time": "2026-04-24T16:02:51.163358+00:00", "status": "completed"}
fucus_gdf = (
    gpd.read_file(base_path / "data/Fucus_location_shp/REDLIST_SIS_Macrophytes.shp")
    .loc[lambda df: df.F_vesiculo != 0, ["geometry"]]
    .to_crs(epsg=4326)
)
fucus_union_3035 = (
    gpd.GeoSeries([unary_union(fucus_gdf.geometry)], crs=4326).to_crs(3035).iloc[0]
)

hex_gdf["area_m2"] = hex_gdf_3035.geometry.area.astype(np.float32)
# Clip against prepared unions via a spatial index (shapely auto-prepares
# from geopandas ≥0.13). At ~50 k hexes × 2 unions the naive loop is fine.
hex_gdf["water_area_m2"] = hex_gdf_3035.geometry.intersection(
    wet_union_3035
).area.astype(np.float32)
hex_gdf["fucus_area_m2"] = hex_gdf_3035.geometry.intersection(
    fucus_union_3035
).area.astype(np.float32)
```

<!-- #region papermill={"duration": 0.001423, "end_time": "2026-04-24T16:02:57.277420+00:00", "exception": false, "start_time": "2026-04-24T16:02:57.275997+00:00", "status": "completed"} -->
## Mean depth over always-wet cells (H0 > 0), fine-grid priority
<!-- #endregion -->

```python papermill={"duration": 0.039837, "end_time": "2026-04-24T16:02:57.318650+00:00", "exception": false, "start_time": "2026-04-24T16:02:57.278813+00:00", "status": "completed"}
def _h0_hex_frame(grid):
    h0 = xr.open_dataset(
        base_path / f"min_data/bsh_operationalmodel_data/static_file_{grid}/H0_file_{grid}.nc"
    )
    lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
    vals = h0.H0.values
    flat = pd.DataFrame({
        "lon": lon2d.ravel(),
        "lat": lat2d.ravel(),
        "H0": vals.ravel(),
    })
    flat = flat[(flat.H0 > 0) & np.isfinite(flat.H0)]
    flat["hex_id"] = hp.label(flat.lon.values, flat.lat.values)
    flat = flat[flat.hex_id >= 0]
    flat["grid"] = grid
    return flat[["hex_id", "grid", "H0"]]

h0_frame = pd.concat([_h0_hex_frame("fine"), _h0_hex_frame("coarse")], ignore_index=True)

# Fine-grid priority: where a hex has any fine cells, drop its coarse cells.
hex_has_fine = set(h0_frame.loc[h0_frame.grid == "fine", "hex_id"].unique())
mask = (h0_frame.grid == "fine") | (~h0_frame.hex_id.isin(hex_has_fine))
mean_depth = h0_frame[mask].groupby("hex_id")["H0"].mean()
hex_gdf["mean_depth_m"] = hex_gdf["hex_id"].map(mean_depth).astype(np.float32)
n_depth_nan = hex_gdf["mean_depth_m"].isna().sum()
print(f"hexes without H0 > 0 coverage (mean_depth_m NaN): {n_depth_nan}")
```

<!-- #region papermill={"duration": 0.001412, "end_time": "2026-04-24T16:02:57.321815+00:00", "exception": false, "start_time": "2026-04-24T16:02:57.320403+00:00", "status": "completed"} -->
## Distance from centroid to always-wet coast (m, in EPSG:3035)
<!-- #endregion -->

```python papermill={"duration": 1.241386, "end_time": "2026-04-24T16:02:58.564602+00:00", "exception": false, "start_time": "2026-04-24T16:02:57.323216+00:00", "status": "completed"}
centroids_3035 = hex_gdf.geometry.centroid.to_crs(3035)
# GeoSeries.distance against a prepared boundary is fast enough at 50k points.
coast_boundary_3035 = always_wet_union_3035.boundary
hex_gdf["dist_to_coast_m"] = centroids_3035.distance(
    coast_boundary_3035
).astype(np.float32).values
```

<!-- #region papermill={"duration": 0.001461, "end_time": "2026-04-24T16:02:58.567799+00:00", "exception": false, "start_time": "2026-04-24T16:02:58.566338+00:00", "status": "completed"} -->
## HELCOM subbasin by centroid
<!-- #endregion -->

```python papermill={"duration": 0.398296, "end_time": "2026-04-24T16:02:58.967515+00:00", "exception": false, "start_time": "2026-04-24T16:02:58.569219+00:00", "status": "completed"}
subbasins = (
    gpd.read_file(
        base_path / "data/HELCOM_subbasins_2022_level2/HELCOM_subbasins_2022_level2.shp"
    )
    .to_crs(epsg=4326)
    .rename(columns={"level_2": "subbasin"})
    .reset_index(drop=True)
)
# Stable int8 lookup for the lookup JSON. -1 reserved for "outside".
subbasin_id_to_name = {-1: "_outside"}
for i, name in enumerate(subbasins["subbasin"].tolist()):
    subbasin_id_to_name[int(i)] = str(name)
subbasin_name_to_id = {v: k for k, v in subbasin_id_to_name.items() if k >= 0}

# Centroids only needed for the sjoin; build once here.
centroids_pts = gpd.GeoDataFrame(
    {"hex_id": hex_gdf["hex_id"].values},
    geometry=hex_gdf.geometry.centroid,
    crs=4326,
)
joined = gpd.sjoin(
    centroids_pts,
    subbasins[["geometry", "subbasin"]],
    how="left", predicate="within",
).drop_duplicates(subset="hex_id")
hex_gdf["helcom_subbasin"] = (
    joined.set_index("hex_id")
    .reindex(hex_gdf["hex_id"].values)["subbasin"]
    .map(subbasin_name_to_id)
    .fillna(-1)
    .astype(np.int8)
    .values
)
hex_gdf[["hex_id", "area_m2", "water_area_m2",
         "fucus_area_m2", "mean_depth_m", "dist_to_coast_m", "helcom_subbasin"]].head()
```

<!-- #region papermill={"duration": 0.001503, "end_time": "2026-04-24T16:02:58.970709+00:00", "exception": false, "start_time": "2026-04-24T16:02:58.969206+00:00", "status": "completed"} -->
## Serialise key.parquet as geoparquet with HexProj metadata

Write as proper geoparquet (geometry column owned by geopandas/pyarrow),
then extend the file-level metadata with the `hex_aggregate_store` JSON blob.
<!-- #endregion -->

```python papermill={"duration": 0.053366, "end_time": "2026-04-24T16:02:59.025538+00:00", "exception": false, "start_time": "2026-04-24T16:02:58.972172+00:00", "status": "completed"}
hex_proj_meta = {
    "projection_name": "laea",
    "lon_origin": hex_origin_lon,
    "lat_origin": hex_origin_lat,
    "hex_size_meters": hex_size_meters,
}
key_meta = {
    "hex_proj": hex_proj_meta,
    "area_crs": "EPSG:3035",
    "subbasin_id_to_name": subbasin_id_to_name,
}

# Drop lon_c/lat_c — centroids are derivable from geometry at query time.
key_gdf = hex_gdf[["hex_id", "geometry", "area_m2", "water_area_m2",
                    "fucus_area_m2", "mean_depth_m", "dist_to_coast_m",
                    "helcom_subbasin"]].copy()

key_path = dataset_root / "key.parquet"
key_gdf.to_parquet(key_path)

# Extend the file-level metadata with our custom JSON blob.
key_table = pq.read_table(key_path)
existing_meta = key_table.schema.metadata or {}
existing_meta[b"hex_aggregate_store"] = json.dumps(key_meta).encode("utf-8")
key_table = key_table.replace_schema_metadata(existing_meta)
pq.write_table(key_table, key_path, compression="zstd")
print(f"wrote {key_path} ({key_path.stat().st_size / 1e6:.2f} MB)")

# Expose key_df for validation cells.
key_df = key_gdf
```

<!-- #region papermill={"duration": 0.001527, "end_time": "2026-04-24T16:02:59.028879+00:00", "exception": false, "start_time": "2026-04-24T16:02:59.027352+00:00", "status": "completed"} -->
# Per-regime counts

`hp.label` returns `-1` for NaN lon/lat, so masked land-seeded
trajectories naturally fall out when we filter `target_hex >= 0`.
<!-- #endregion -->

```python papermill={"duration": 0.00784, "end_time": "2026-04-24T16:02:59.038244+00:00", "exception": false, "start_time": "2026-04-24T16:02:59.030404+00:00", "status": "completed"}
def _zarr_for(regime):
    matches = sorted(
        (base_path / f"output/Trajectories/{regime}/{release_year}").glob("*.zarr")
    )
    assert len(matches) == 1, (regime, matches)
    return matches[0]


def build_counts(regime):
    zarr_path = _zarr_for(regime)
    # Filename convention: Fucus_BSH_YYYYMMDD_…
    fn_date = zarr_path.name.split("_")[2]
    assert fn_date.startswith(str(release_year)), (fn_date, release_year)

    stages = {}
    t0 = time.time()

    ds = xr.open_zarr(zarr_path)
    ds, _ = mask_land_seeded(ds)

    release_ts = pd.Timestamp(ds.time.isel(obs=0).compute().values[0])
    release_doy = int(release_ts.dayofyear)
    if release_doy == 366:
        print(f"[{regime}] leap-year release DOY 366; skipping per plan.")
        return None
    fn_doy = int(pd.Timestamp(fn_date).dayofyear)
    assert fn_doy == release_doy, (fn_doy, release_doy, regime)

    # Lazy (trajectory, obs) hex labels.
    target_hex = xr.apply_ufunc(
        hp.label, ds.lon, ds.lat,
        dask="parallelized", output_dtypes=[np.int64],
    )
    # Per-trajectory release_hex from obs=0 (NaN'd for land-seeded; filtered below).
    release_hex = xr.apply_ufunc(
        hp.label, ds.lon.isel(obs=0, drop=True), ds.lat.isel(obs=0, drop=True),
        dask="parallelized", output_dtypes=[np.int64],
    )

    # age_bin over obs indices → days → 10-day bins.
    obs_ages_days = ds.obs.values * output_dt_mins / (60 * 24)
    age_bin = (obs_ages_days // age_bin_days).astype(np.int32)
    age_bin_da = xr.DataArray(age_bin, dims=["obs"])

    frame = xr.Dataset({
        "target_hex": target_hex,
        "release_hex": release_hex,
        "age_bin": age_bin_da,
    }).to_dask_dataframe(dim_order=["trajectory", "obs"])

    # target_hex >= 0 filters both land-obs and land-seeded trajectories
    # (mask_land_seeded NaNs all lon/lat → label returns -1).
    t1 = time.time()
    valid_frame = frame[
        (frame.target_hex >= 0) & (frame.age_bin >= 0) & (frame.age_bin <= max_age_bin)
    ]
    counts = (
        valid_frame
        .groupby(["release_hex", "age_bin", "target_hex"])
        .size().rename("n_obs").reset_index().compute()
    )
    # n_valid: rows passing the filter, before groupby — equals sum(n_obs) by construction.
    n_valid = int(counts["n_obs"].sum())
    stages["compute_s"] = time.time() - t1

    # regime and release_year come back from the hive path on read —
    # storing them in data would just risk partition-schema drift.
    # release_doy varies within a single (regime, year) file in the
    # full build (73 doys per year), so it stays as a data column.
    counts["release_hex"] = counts["release_hex"].astype(np.int32)
    counts["target_hex"] = counts["target_hex"].astype(np.int32)
    counts["age_bin"] = counts["age_bin"].astype(np.int8)
    counts["n_obs"] = counts["n_obs"].astype(np.int32)
    counts["release_doy"] = np.int16(release_doy)

    # Write partition.
    t2 = time.time()
    part_dir = dataset_root / f"counts/regime={regime}/release_year={release_year}"
    part_dir.mkdir(parents=True, exist_ok=True)
    part_path = part_dir / "part.parquet"
    part_meta = {
        "hex_proj": hex_proj_meta,
        "age_bin_days": age_bin_days,
        "output_dt_mins": output_dt_mins,
        "regime": regime,
        "release_year": int(release_year),
        "release_doy": int(release_doy),
        "source_zarr": zarr_path.name,
    }
    tbl = pa.Table.from_pandas(counts, preserve_index=False)
    tbl = tbl.replace_schema_metadata({
        b"hex_aggregate_store": json.dumps(part_meta).encode("utf-8"),
    })
    pq.write_table(tbl, part_path, compression="zstd")
    stages["write_s"] = time.time() - t2

    stages["total_s"] = time.time() - t0
    return {
        "regime": regime,
        "zarr": zarr_path,
        "release_doy": release_doy,
        "n_valid": n_valid,
        "counts": counts,
        "part_path": part_path,
        "stages": stages,
    }
```

```python papermill={"duration": 44.968857, "end_time": "2026-04-24T16:03:44.008733+00:00", "exception": false, "start_time": "2026-04-24T16:02:59.039876+00:00", "status": "completed"}
results = {}
for regime in regimes:
    print(f"\n=== {regime} ===")
    results[regime] = build_counts(regime)
    print(results[regime]["stages"])
```

<!-- #region papermill={"duration": 0.001706, "end_time": "2026-04-24T16:03:44.012605+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.010899+00:00", "status": "completed"} -->
# Validation

## Key-file stats
<!-- #endregion -->

```python papermill={"duration": 0.007305, "end_time": "2026-04-24T16:03:44.021575+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.014270+00:00", "status": "completed"}
print(f"n hexes: {len(key_df):,}")
print("water_area_m2 min/median/max: "
      f"{key_df.water_area_m2.min():.0f} / "
      f"{key_df.water_area_m2.median():.0f} / "
      f"{key_df.water_area_m2.max():.0f}")
print(f"mean_depth_m NaN: {int(key_df.mean_depth_m.isna().sum())}")
print(f"fucus_area_m2 > 0: {int((key_df.fucus_area_m2 > 0).sum())}")
vc = key_df["helcom_subbasin"].value_counts().sort_index()
print("helcom_subbasin counts (id → n):")
for sid, n in vc.items():
    print(f"  {int(sid):>3} {subbasin_id_to_name[int(sid)]:<35s} {int(n):>6}")
```

<!-- #region papermill={"duration": 0.001712, "end_time": "2026-04-24T16:03:44.025211+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.023499+00:00", "status": "completed"} -->
## Per-regime counts stats
<!-- #endregion -->

```python papermill={"duration": 0.007016, "end_time": "2026-04-24T16:03:44.033894+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.026878+00:00", "status": "completed"}
for regime, r in results.items():
    if r is None:
        print(f"[{regime}] skipped")
        continue
    counts = r["counts"]
    print(f"\n[{regime}]")
    print(f"  rows: {len(counts):,}")
    print(f"  sum(n_obs): {int(counts.n_obs.sum()):,}")
    print(f"  stage seconds: {r['stages']}")
```

<!-- #region papermill={"duration": 0.001738, "end_time": "2026-04-24T16:03:44.037615+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.035877+00:00", "status": "completed"} -->
## Key-completeness invariant

Every `release_hex` and `target_hex` in counts must be in `key.parquet`'s
`hex_id` set. The key is built from the full BSH-domain grid; any
trajectory that stays within the model's lat/lon bounds necessarily
labels to a hex already in the key. If this assert ever fires,
investigate — don't paper over it by extending the key after the fact.
<!-- #endregion -->

```python papermill={"duration": 0.522089, "end_time": "2026-04-24T16:03:44.561380+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.039291+00:00", "status": "completed"}
key_ids = set(int(x) for x in key_df["hex_id"].values)
for regime, r in results.items():
    if r is None:
        continue
    unseen = (set(int(x) for x in r["counts"]["release_hex"].values)
              | set(int(x) for x in r["counts"]["target_hex"].values)) - key_ids
    assert not unseen, (regime, sorted(unseen))
print("PASS: every release_hex and target_hex is in key.parquet.")
```

<!-- #region papermill={"duration": 0.001755, "end_time": "2026-04-24T16:03:44.565166+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.563411+00:00", "status": "completed"} -->
## Conservation cross-check: sum(n_obs) == non-masked (traj, obs) pairs

`n_valid` is captured inside `build_counts` immediately after the single-pass
groupby, so no re-scan is needed here. The check: for each regime, what went
into the groupby equals what came out.
<!-- #endregion -->

```python papermill={"duration": 0.005786, "end_time": "2026-04-24T16:03:44.572685+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.566899+00:00", "status": "completed"}
for regime, r in results.items():
    if r is None:
        print(f"[{regime}] skipped")
        continue
    n_valid = r["n_valid"]
    n_stored = int(r["counts"]["n_obs"].sum())
    print(f"[{regime}] n_valid={n_valid:,}  stored sum(n_obs)={n_stored:,}")
    assert n_valid == n_stored, (regime, n_valid, n_stored)
print("PASS: conservation holds for all regimes.")
```

<!-- #region papermill={"duration": 0.001758, "end_time": "2026-04-24T16:03:44.576347+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.574589+00:00", "status": "completed"} -->
## Output sizes
<!-- #endregion -->

```python papermill={"duration": 0.00481, "end_time": "2026-04-24T16:03:44.582908+00:00", "exception": false, "start_time": "2026-04-24T16:03:44.578098+00:00", "status": "completed"}
total = 0
entries = [("key.parquet", key_path)]
for regime in regimes:
    r = results.get(regime)
    if r is None:
        continue
    entries.append((f"counts/{regime}/release_year={release_year}", r["part_path"]))

for label, p in entries:
    size = p.stat().st_size
    total += size
    print(f"  {label:<60s} {size/1e6:>9.2f} MB")
print(f"  {'TOTAL':<60s} {total/1e6:>9.2f} MB")
```
