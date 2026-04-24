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

# Build hex-aggregated dispersal store

Distil the per-regime trajectory zarrs into the compact aggregate
store defined in `plans/hex_aggregate_store.md`. This notebook builds
one store per entry in `hex_configs` (see the cell below the parameters
cell). Each store contains:

- `key.parquet` — one row per hex in the store's domain (geometry,
  area, water area, depth, coast distance, Fucus area, HELCOM
  subbasin). Built once from static inputs (BSH wet-region geojson,
  H0 grids, Fucus shapefile, HELCOM level-2 polygons). Land hexes
  are included with `water_area_m2 = 0` so the key covers every hex
  `hp.label` can assign to a trajectory position within the domain.
- `counts/regime=…/release_year=…/part.parquet` — one row per
  `(release_hex, age_bin, target_hex)` with `n_obs` aggregate. Built
  per regime × release_year from the zarrs. Only (release_hex,
  target_hex) pairs whose centroids fall within the store's domain
  are written; this keeps specialised stores (e.g. German-waters zoom)
  small without violating the key-completeness invariant.

Units are SI throughout (m², m). The HexProj configuration travels
as parquet file-level metadata so `hex_id` values can be rematerialised
into geometry downstream.

```python
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

from helpers import load_trajectories, mask_land_seeded, parse_zarr_stem
```

# Parameters

```python tags=["parameters"]
data_root = "../data"
output_root = "../output"
bsh_root = "../data/bsh_minimal"

# Age binning.
age_bin_days = 10        # 22 bins × 10 d = 220 d of drift.
output_dt_mins = 60      # zarr output cadence.

# Dev scope: one release year to build.
release_year = 2019
```

```python
# Aggregate stores to build. Each config spawns one key+counts set with its
# own HexProj. baltic_r6km_v1 is the basin-wide reference; de_r4km_v1 is the
# German-waters zoom (matches notebook 025's DE panel extent).
hex_configs = [
    {
        "name": "baltic_r6km_v1",
        "hex_size_meters": 6000,
        "hex_origin_lon": 18.0,
        "hex_origin_lat": 59.0,
        "lon_bounds": None,          # full BSH domain
        "lat_bounds": None,
    },
    {
        "name": "de_r4km_v1",
        "hex_size_meters": 4000,
        "hex_origin_lon": 11.5,
        "hex_origin_lat": 54.35,
        "lon_bounds": (7.0, 16.0),   # DE zoom + ~1° margin
        "lat_bounds": (52.5, 56.0),
    },
]
```

```python
data_root = Path(data_root)
output_root = Path(output_root)
bsh_root = Path(bsh_root)

n_age_bins = 220 // age_bin_days   # 22
max_age_bin = n_age_bins - 1        # 21 (bins 0..21)
```

# Regime discovery

```python
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

# Dask cluster

```python
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

# Static inputs (loaded once, shared across all hex_configs)

BSH wet-region polygons and Fucus shapefile are independent of the hex
grid; load them once outside the per-config loop.

```python
wet_gdf = gpd.read_file(data_root / "bsh_coastline/coastline.geojson")
always_wet_gdf = gpd.read_file(
    data_root / "bsh_coastline/coastline_always_wet.geojson"
)
wet_union = unary_union(wet_gdf.geometry)
always_wet_union = unary_union(always_wet_gdf.geometry)

fucus_gdf = (
    gpd.read_file(data_root / "fucus_redlist_shapefile/REDLIST_SIS_Macrophytes.shp")
    .loc[lambda df: df.F_vesiculo != 0, ["geometry"]]
    .to_crs(epsg=4326)
)
fucus_union_3035 = (
    gpd.GeoSeries([unary_union(fucus_gdf.geometry)], crs=4326).to_crs(3035).iloc[0]
)

subbasins = (
    gpd.read_file(
        data_root / "helcom_subbasins/HELCOM_subbasins_2022_level2.shp"
    )
    .to_crs(epsg=4326)
    .rename(columns={"level_2": "subbasin"})
    .reset_index(drop=True)
)
# Stable int8 lookup for the JSON metadata. -1 reserved for "outside".
subbasin_id_to_name = {-1: "_outside"}
for i, name in enumerate(subbasins["subbasin"].tolist()):
    subbasin_id_to_name[int(i)] = str(name)
subbasin_name_to_id = {v: k for k, v in subbasin_id_to_name.items() if k >= 0}
```

# Build all stores

Outer loop over hex_configs. Each iteration:
1. Builds a HexProj and decomposes the BSH H0 grid to get domain hex_ids.
2. If lon_bounds/lat_bounds are set, filters hex_ids to those whose
   centroid falls within the bounding box (keeps DE store small).
3. Builds key.parquet with per-hex geometry and attributes.
4. Builds per-regime counts partitions, filtered to hex_ids in this store.

```python
all_store_summaries = []

for config in hex_configs:
    config_name = config["name"]
    hex_size_meters = config["hex_size_meters"]
    hex_origin_lon = config["hex_origin_lon"]
    hex_origin_lat = config["hex_origin_lat"]
    lon_bounds = config["lon_bounds"]
    lat_bounds = config["lat_bounds"]

    print(f"\n{'='*60}")
    print(f"Building store: {config_name}")
    print(f"  hex_size_meters={hex_size_meters}, origin=({hex_origin_lon}, {hex_origin_lat})")
    if lon_bounds:
        print(f"  bbox filter: lon={lon_bounds}, lat={lat_bounds}")

    dataset_root = output_root / f"HexAggregates/{config_name}"
    dataset_root.mkdir(parents=True, exist_ok=True)

    # Build HexProj for this config.
    hp = HexProj(
        projection_name="laea",
        lon_origin=hex_origin_lon,
        lat_origin=hex_origin_lat,
        hex_size_meters=hex_size_meters,
    )

    # Label every fine+coarse H0 grid point to get the full BSH domain
    # hex set for this projection.
    _domain_ids = set()
    for grid in ("fine", "coarse"):
        h0 = xr.open_dataset(
            bsh_root / f"static_file_{grid}/H0_file_{grid}.nc"
        )
        lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
        labels = hp.label(lon2d.ravel(), lat2d.ravel())
        _domain_ids |= set(int(x) for x in labels if x >= 0)
    hex_ids = np.asarray(sorted(_domain_ids), dtype=np.int32)
    print(f"  BSH-domain hexes before bbox filter: {len(hex_ids):,}")

    # If this config has a bounding-box filter, restrict hex_ids to those
    # whose centroid falls within it. This keeps specialised stores (e.g.
    # German-waters zoom) small without breaking key-completeness for their
    # domain: only (release_hex, target_hex) pairs in-domain get written.
    if lon_bounds is not None and lat_bounds is not None:
        lon_lo, lon_hi = lon_bounds
        lat_lo, lat_hi = lat_bounds
        centroid_gdf = hp.to_geodataframe(hex_ids.tolist())
        centroids = centroid_gdf.geometry.centroid
        in_bbox = (
            (centroids.x >= lon_lo) & (centroids.x <= lon_hi)
            & (centroids.y >= lat_lo) & (centroids.y <= lat_hi)
        )
        hex_ids = np.asarray(
            sorted(int(hid) for hid, keep in zip(hex_ids, in_bbox) if keep),
            dtype=np.int32,
        )
        print(f"  hex_ids after bbox filter: {len(hex_ids):,}")

    hex_ids_set = set(int(x) for x in hex_ids)

    # ------------------------------------------------------------------
    # Key file
    # ------------------------------------------------------------------

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

    hex_gdf["area_m2"] = hex_gdf_3035.geometry.area.astype(np.float32)
    hex_gdf["water_area_m2"] = hex_gdf_3035.geometry.intersection(
        wet_union_3035
    ).area.astype(np.float32)
    hex_gdf["fucus_area_m2"] = hex_gdf_3035.geometry.intersection(
        fucus_union_3035
    ).area.astype(np.float32)

    # Mean depth over always-wet cells (H0 > 0), fine-grid priority.
    def _h0_hex_frame(grid, hp_=hp, bsh_root_=bsh_root):
        h0 = xr.open_dataset(
            bsh_root_ / f"static_file_{grid}/H0_file_{grid}.nc"
        )
        lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
        vals = h0.H0.values
        flat = pd.DataFrame({
            "lon": lon2d.ravel(),
            "lat": lat2d.ravel(),
            "H0": vals.ravel(),
        })
        flat = flat[(flat.H0 > 0) & np.isfinite(flat.H0)]
        flat["hex_id"] = hp_.label(flat.lon.values, flat.lat.values)
        flat = flat[flat.hex_id >= 0]
        flat["grid"] = grid
        return flat[["hex_id", "grid", "H0"]]

    h0_frame = pd.concat(
        [_h0_hex_frame("fine"), _h0_hex_frame("coarse")], ignore_index=True
    )
    hex_has_fine = set(h0_frame.loc[h0_frame.grid == "fine", "hex_id"].unique())
    mask = (h0_frame.grid == "fine") | (~h0_frame.hex_id.isin(hex_has_fine))
    mean_depth = h0_frame[mask].groupby("hex_id")["H0"].mean()
    hex_gdf["mean_depth_m"] = hex_gdf["hex_id"].map(mean_depth).astype(np.float32)
    n_depth_nan = hex_gdf["mean_depth_m"].isna().sum()
    print(f"  hexes without H0 > 0 coverage (mean_depth_m NaN): {n_depth_nan}")

    # Distance from centroid to always-wet coast (m, in EPSG:3035).
    centroids_3035 = hex_gdf.geometry.centroid.to_crs(3035)
    coast_boundary_3035 = always_wet_union_3035.boundary
    hex_gdf["dist_to_coast_m"] = centroids_3035.distance(
        coast_boundary_3035
    ).astype(np.float32).values

    # HELCOM subbasin by centroid.
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

    # Serialise key.parquet as geoparquet with HexProj metadata.
    # Capture config values in locals before writing to avoid closure trap
    # (the loop variable config is rebound on each iteration; writing
    # hex_proj_meta after the loop would use the last config's values).
    _hex_proj_meta = {
        "projection_name": "laea",
        "lon_origin": hex_origin_lon,
        "lat_origin": hex_origin_lat,
        "hex_size_meters": hex_size_meters,
    }
    _key_meta = {
        "hex_proj": _hex_proj_meta,
        "area_crs": "EPSG:3035",
        "subbasin_id_to_name": subbasin_id_to_name,
    }

    key_gdf = hex_gdf[["hex_id", "geometry", "area_m2", "water_area_m2",
                        "fucus_area_m2", "mean_depth_m", "dist_to_coast_m",
                        "helcom_subbasin"]].copy()

    key_path = dataset_root / "key.parquet"
    key_gdf.to_parquet(key_path)

    key_table = pq.read_table(key_path)
    existing_meta = key_table.schema.metadata or {}
    existing_meta[b"hex_aggregate_store"] = json.dumps(_key_meta).encode("utf-8")
    key_table = key_table.replace_schema_metadata(existing_meta)
    pq.write_table(key_table, key_path, compression="zstd")
    print(f"  wrote {key_path} ({key_path.stat().st_size / 1e6:.2f} MB)")

    # ------------------------------------------------------------------
    # Per-regime counts
    # ------------------------------------------------------------------

    def _zarr_for(regime, rl_year=release_year, output_root_=output_root):
        matches = sorted(
            (output_root_ / f"Trajectories/{regime}/{rl_year}").glob("*.zarr")
        )
        assert len(matches) == 1, (regime, matches)
        return matches[0]

    def build_counts(regime, hp_=hp, hids_set=hex_ids_set,
                     dataset_root_=dataset_root,
                     _hex_proj_meta_=_hex_proj_meta,
                     rl_year=release_year):
        zarr_path = _zarr_for(regime, rl_year)
        fn_release_date, fn_regime = parse_zarr_stem(zarr_path)
        assert fn_release_date.year == rl_year, (fn_release_date, rl_year)
        assert fn_regime == regime, (fn_regime, regime)

        stages = {}
        t0 = time.time()

        ds = xr.open_zarr(zarr_path)
        ds, _ = mask_land_seeded(ds)

        release_ts = pd.Timestamp(ds.time.isel(obs=0).compute().values[0])
        release_doy = int(release_ts.dayofyear)
        if release_doy == 366:
            print(f"  [{regime}] leap-year release DOY 366; skipping per plan.")
            return None
        fn_doy = int(fn_release_date.dayofyear)
        assert fn_doy == release_doy, (fn_doy, release_doy, regime)

        # Lazy (trajectory, obs) hex labels.
        target_hex = xr.apply_ufunc(
            hp_.label, ds.lon, ds.lat,
            dask="parallelized", output_dtypes=[np.int64],
        )
        release_hex = xr.apply_ufunc(
            hp_.label, ds.lon.isel(obs=0, drop=True), ds.lat.isel(obs=0, drop=True),
            dask="parallelized", output_dtypes=[np.int64],
        )

        obs_ages_days = ds.obs.values * output_dt_mins / (60 * 24)
        age_bin = (obs_ages_days // age_bin_days).astype(np.int32)
        age_bin_da = xr.DataArray(age_bin, dims=["obs"])

        frame = xr.Dataset({
            "target_hex": target_hex,
            "release_hex": release_hex,
            "age_bin": age_bin_da,
        }).to_dask_dataframe(dim_order=["trajectory", "obs"])

        t1 = time.time()
        valid_frame = frame[
            (frame.target_hex >= 0) & (frame.age_bin >= 0) & (frame.age_bin <= max_age_bin)
        ]
        counts = (
            valid_frame
            .groupby(["release_hex", "age_bin", "target_hex"])
            .size().rename("n_obs").reset_index().compute()
        )
        n_valid = int(counts["n_obs"].sum())
        stages["compute_s"] = time.time() - t1

        counts["release_hex"] = counts["release_hex"].astype(np.int32)
        counts["target_hex"] = counts["target_hex"].astype(np.int32)
        counts["age_bin"] = counts["age_bin"].astype(np.int8)
        counts["n_obs"] = counts["n_obs"].astype(np.int32)
        counts["release_doy"] = np.int16(release_doy)

        # Filter to hex_ids in this store's domain. This preserves
        # key-completeness: every (release_hex, target_hex) stored here
        # is guaranteed to be in key.parquet for this config.
        in_domain = counts["release_hex"].isin(hids_set) & counts["target_hex"].isin(hids_set)
        counts = counts[in_domain].copy()

        t2 = time.time()
        part_dir = dataset_root_ / f"counts/regime={regime}/release_year={rl_year}"
        part_dir.mkdir(parents=True, exist_ok=True)
        part_path = part_dir / "part.parquet"
        part_meta = {
            "hex_proj": _hex_proj_meta_,
            "age_bin_days": age_bin_days,
            "output_dt_mins": output_dt_mins,
            "regime": regime,
            "release_year": int(rl_year),
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

    results = {}
    for regime in regimes:
        print(f"\n  --- {regime} ---")
        results[regime] = build_counts(regime)
        if results[regime] is not None:
            print(f"  {results[regime]['stages']}")

    # ------------------------------------------------------------------
    # Per-config summary
    # ------------------------------------------------------------------

    n_key_hexes = len(key_gdf)
    total_store_bytes = key_path.stat().st_size
    print(f"\n  Summary for {config_name}:")
    print(f"    key hexes: {n_key_hexes:,}")
    for regime in regimes:
        r = results.get(regime)
        if r is None:
            print(f"    {regime}: skipped")
            continue
        rows = len(r["counts"])
        total_obs = int(r["counts"]["n_obs"].sum())
        total_store_bytes += r["part_path"].stat().st_size
        print(f"    {regime}: {rows:,} rows, sum(n_obs)={total_obs:,}")
    print(f"    total store size: {total_store_bytes / 1e6:.2f} MB")

    all_store_summaries.append({
        "name": config_name,
        "n_key_hexes": n_key_hexes,
        "results": results,
    })
```

# Validation

## Key-completeness invariant

Every `release_hex` and `target_hex` in each store's counts must be in
that store's `key.parquet`. Any violation is a data integrity bug.

```python
for summary in all_store_summaries:
    config_name_ = summary["name"]
    dataset_root_ = output_root / f"HexAggregates/{config_name_}"
    key_df_ = gpd.read_parquet(dataset_root_ / "key.parquet")
    key_ids_ = set(int(x) for x in key_df_["hex_id"].values)
    for regime, r in summary["results"].items():
        if r is None:
            continue
        unseen = (
            set(int(x) for x in r["counts"]["release_hex"].values)
            | set(int(x) for x in r["counts"]["target_hex"].values)
        ) - key_ids_
        assert not unseen, (config_name_, regime, sorted(unseen))
    print(f"PASS [{config_name_}]: every release_hex and target_hex is in key.parquet.")
```

## Conservation cross-check

`n_valid` is the pre-domain-filter count (all obs passing age/hex
filters); the stored `sum(n_obs)` may be smaller for stores with a
bbox filter (trajectories that drift outside the store domain are
excluded). Print both for audit; the cross-check within `build_counts`
already confirmed that `n_valid == sum(n_obs)` before the domain filter.

```python
for summary in all_store_summaries:
    print(f"\n[{summary['name']}]")
    for regime, r in summary["results"].items():
        if r is None:
            print(f"  [{regime}] skipped")
            continue
        n_valid = r["n_valid"]
        n_stored = int(r["counts"]["n_obs"].sum())
        print(f"  [{regime}] n_valid (pre-filter)={n_valid:,}  stored sum(n_obs)={n_stored:,}")
```

## Output sizes

```python
for summary in all_store_summaries:
    config_name_ = summary["name"]
    dataset_root_ = output_root / f"HexAggregates/{config_name_}"
    print(f"\n[{config_name_}]")
    total = 0
    key_p = dataset_root_ / "key.parquet"
    size = key_p.stat().st_size
    total += size
    print(f"  {'key.parquet':<60s} {size/1e6:>9.2f} MB")
    for regime in summary["results"]:
        r = summary["results"][regime]
        if r is None:
            continue
        size = r["part_path"].stat().st_size
        total += size
        label = f"counts/{regime}/release_year={release_year}"
        print(f"  {label:<60s} {size/1e6:>9.2f} MB")
    print(f"  {'TOTAL':<60s} {total/1e6:>9.2f} MB")
```
