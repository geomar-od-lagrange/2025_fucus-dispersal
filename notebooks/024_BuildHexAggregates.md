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
store defined in `plans/hex_aggregate_store.md`. One store per run, at
the configured `hex_radius`; multiple radii come from papermill sweeps,
not an in-notebook loop. Each store contains:

- `key.parquet` — one row per hex in the BSH domain (geometry, area,
  water area, depth, coast distance, Fucus area, HELCOM subbasin).
  Built once from static inputs (BSH wet-region geojson, H0 grids,
  Fucus shapefile, HELCOM level-2 polygons). Land hexes are included
  with `water_area_m2 = 0` so the key covers every hex `hp.label` can
  assign to a trajectory position within the domain.
- `counts/regime=…/release_year=…/part.parquet` — one row per
  `(release_hex, age_bin, target_hex)` with `n_obs` aggregate. Built
  per regime × release_year from the zarrs.

Units are SI throughout (m², m). The HexProj configuration travels
as parquet file-level metadata so `hex_id` values can be rematerialised
into geometry downstream.

```python
import json
import os
import re
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import xarray as xr
from dask.distributed import Client
from shapely.ops import unary_union

from hextraj import HexProj
```

```python
# Pattern: Fucus_BSH_YYYYMMDD_{regime}_dt{N}min_seed{S}
# surface_stokes must precede surface so the alternation matches the
# longer form first. Multiple seeds for the same (date, regime) are
# expected — independent reseeded resubmissions for sample-size
# expansion — and aggregated additively on the n_obs counts.
_ZARR_STEM_RE = re.compile(
    r"^Fucus_BSH_(\d{8})_(surface_stokes|surface|bottom)_dt\d+min_seed\d+$"
)


def parse_zarr_stem(path):
    """Parse a trajectory zarr filename into ``(release_date, regime)``.

    Authoritative format (notebook 010):
    ``Fucus_BSH_{YYYYMMDD}_{regime}_dt{N}min_seed{S}.zarr``, where
    ``{regime}`` is ``surface``, ``surface_stokes``, or ``bottom``.
    """
    path = Path(path)
    m = _ZARR_STEM_RE.match(path.stem)
    if m is None:
        raise ValueError(
            f"zarr filename does not match expected pattern "
            f"'Fucus_BSH_YYYYMMDD_<regime>_…_seed<S>': {path.name!r}"
        )
    return pd.Timestamp(m.group(1)), m.group(2)
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (HELCOM polygons, Fucus shapefile, BSH coastline).
data_root = "../data"
# Read root of trajectory zarrs and write root for the hex-aggregate store.
output_root = "../output"
# Read root of BSH H0 / static inputs (defines the hex domain). Demo
# default points at the twin's static subset; NESH runs override to the
# full BSH store.
bsh_root = "../data/bsh_hbmnoku_static"

# Hex radius (corner-to-centre distance) in metres. Sweep this via
# papermill to build stores at multiple radii.
hex_radius = 6000

# Age binning.
age_bin_days = 10        # 22 bins × 10 d = 220 d of drift.
output_dt_mins = 60      # zarr output cadence.

# Dev scope: one release year to build.
release_year = 2019
```

# Derived layout / projection

Layout assumptions encoded by the path-construction below:

- ``bsh_root/static_file_<grid>/H0_file_<grid>.nc`` for ``grid in {fine, coarse}``
- ``output_root/HexAggregates/r{hex_radius}m/`` — this run's store root
- ``output_root/Trajectories/<regime>/<release_year>/*.zarr`` — one zarr per release date
- ``counts/regime=<regime>/release_year=<release_year>/part.parquet`` — counts partition layout

Projection: equal-area (LAEA) centred on the BSH domain centroid
(midpoint of the coarse-grid bounding box, computed below from the H0
files so the centre tracks the actual BSH grid extent).

```python
data_root = Path(data_root)
output_root = Path(output_root)
bsh_root = Path(bsh_root)

n_age_bins = 220 // age_bin_days   # 22
max_age_bin = n_age_bins - 1        # 21 (bins 0..21)

# Compute domain centroid from the coarse-grid H0 (full BSH extent).
_h0_coarse = xr.open_dataset(bsh_root / "static_file_coarse/H0_file_coarse.nc")
domain_lon_origin = float(0.5 * (_h0_coarse.lon.min() + _h0_coarse.lon.max()))
domain_lat_origin = float(0.5 * (_h0_coarse.lat.min() + _h0_coarse.lat.max()))
print(f"BSH domain centroid: lon={domain_lon_origin:.4f}, lat={domain_lat_origin:.4f}")

dataset_root = output_root / "HexAggregates" / f"r{hex_radius}m"
dataset_root.mkdir(parents=True, exist_ok=True)
print(f"Store root: {dataset_root}")
```

# Regime discovery

```python
trajectory_root = output_root / "Trajectories"
regimes = sorted(p.name for p in trajectory_root.iterdir() if p.is_dir())
print(f"Regimes: {regimes}")
```

# Dask cluster

```python
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

# Static inputs

BSH wet-region polygons and Fucus shapefile are independent of the hex
grid; load them once.

```python
wet_gdf = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
always_wet_gdf = gpd.read_file(
    data_root / "bsh_hbmnoku_static/coastline_always_wet.geojson"
)
wet_union = unary_union(wet_gdf.geometry)
always_wet_union = unary_union(always_wet_gdf.geometry)

fucus_gdf = (
    gpd.read_file(data_root / "helcom_fucus_redlist/REDLIST_SIS_Macrophytes.shp")
    .loc[lambda df: df.F_vesiculo != 0, ["geometry"]]
    .to_crs(epsg=4326)
)
fucus_union_3035 = (
    gpd.GeoSeries([unary_union(fucus_gdf.geometry)], crs=4326).to_crs(3035).iloc[0]
)

subbasins = (
    gpd.read_file(
        data_root / "helcom_subbasins_2022/HELCOM_subbasins_2022_level2.shp"
    )
    .to_crs(epsg=4326)
    .rename(columns={"level_2": "subbasin"})
    .reset_index(drop=True)
)
# Stable int8 lookup for the JSON metadata. -1 reserved for "outside".
subbasin_name_to_id = {
    str(name): int(i) for i, name in enumerate(subbasins["subbasin"].tolist())
}
subbasin_id_to_name = {-1: "_outside", **{i: n for n, i in subbasin_name_to_id.items()}}
```

# HexProj and BSH-domain hex set

Build the projection once and label every fine+coarse H0 grid point to
get the full BSH-domain hex set for this radius.

```python
hp = HexProj(
    projection_name="laea",
    lon_origin=domain_lon_origin,
    lat_origin=domain_lat_origin,
    hex_size_meters=hex_radius,
)

_domain_ids = set()
for grid in ("fine", "coarse"):
    h0 = xr.open_dataset(bsh_root / f"static_file_{grid}/H0_file_{grid}.nc")
    lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
    labels = hp.label(lon2d.ravel(), lat2d.ravel())
    _domain_ids |= set(int(x) for x in labels if x >= 0)
hex_ids = np.asarray(sorted(_domain_ids), dtype=np.int32)
print(f"BSH-domain hexes at r={hex_radius} m: {len(hex_ids):,}")
```

# Key file

Per-hex geometry and attributes (area, water area, mean depth, coast
distance, Fucus area, HELCOM subbasin).

```python
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
```

```python
def h0_hex_frame(grid, hp_, bsh_root_):
    """H0 cells (lon, lat, H0 > 0) labelled by hex_id, tagged by grid."""
    h0 = xr.open_dataset(bsh_root_ / f"static_file_{grid}/H0_file_{grid}.nc")
    lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
    flat = pd.DataFrame({
        "lon": lon2d.ravel(),
        "lat": lat2d.ravel(),
        "H0": h0.H0.values.ravel(),
    })
    flat = flat[(flat.H0 > 0) & np.isfinite(flat.H0)]
    flat["hex_id"] = hp_.label(flat.lon.values, flat.lat.values)
    flat = flat[flat.hex_id >= 0]
    flat["grid"] = grid
    return flat[["hex_id", "grid", "H0"]]


# Mean depth over always-wet cells (H0 > 0), fine-grid priority.
h0_frame = pd.concat(
    [h0_hex_frame("fine", hp, bsh_root), h0_hex_frame("coarse", hp, bsh_root)],
    ignore_index=True,
)
hex_has_fine = set(h0_frame.loc[h0_frame.grid == "fine", "hex_id"].unique())
mask = (h0_frame.grid == "fine") | (~h0_frame.hex_id.isin(hex_has_fine))
mean_depth = h0_frame[mask].groupby("hex_id")["H0"].mean()
hex_gdf["mean_depth_m"] = hex_gdf["hex_id"].map(mean_depth).astype(np.float32)
n_depth_nan = hex_gdf["mean_depth_m"].isna().sum()
print(f"hexes without H0 > 0 coverage (mean_depth_m NaN): {n_depth_nan}")

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
hex_proj_meta = {
    "projection_name": "laea",
    "lon_origin": domain_lon_origin,
    "lat_origin": domain_lat_origin,
    "hex_size_meters": hex_radius,
}
key_meta = {
    "hex_proj": hex_proj_meta,
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
existing_meta[b"hex_aggregate_store"] = json.dumps(key_meta).encode("utf-8")
key_table = key_table.replace_schema_metadata(existing_meta)
pq.write_table(key_table, key_path, compression="zstd")
print(f"wrote {key_path} ({key_path.stat().st_size / 1e6:.2f} MB)")
```

# Per-regime counts

```python
def build_counts(regime, hp_, dataset_root_, hex_proj_meta_, rl_year):
    matches = sorted(
        (output_root / f"Trajectories/{regime}/{rl_year}").glob("*.zarr")
    )
    assert len(matches) == 1, (regime, matches)
    zarr_path = matches[0]
    fn_release_date, fn_regime = parse_zarr_stem(zarr_path)
    assert fn_release_date.year == rl_year, (fn_release_date, rl_year)
    assert fn_regime == regime, (fn_regime, regime)

    stages = {}
    t0 = time.time()

    ds = xr.open_zarr(zarr_path)
    # First-step displacement of zero ⇒ trajectory was seeded on land.
    ds = ds.where(~(
        (ds.lon.diff("obs").isel(obs=0, drop=True) == 0)
        & (ds.lat.diff("obs").isel(obs=0, drop=True) == 0)
    ))

    release_ts = pd.Timestamp(ds.time.isel(obs=0).compute().values[0])
    release_doy = int(release_ts.dayofyear)
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

    t2 = time.time()
    part_dir = dataset_root_ / f"counts/regime={regime}/release_year={rl_year}"
    part_dir.mkdir(parents=True, exist_ok=True)
    part_path = part_dir / "part.parquet"
    part_meta = {
        "hex_proj": hex_proj_meta_,
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
    print(f"\n--- {regime} ---")
    results[regime] = build_counts(
        regime, hp, dataset_root, hex_proj_meta, release_year,
    )
    print(f"  {results[regime]['stages']}")
```

# Summary

```python
summary_rows = []
for regime in regimes:
    r = results[regime]
    summary_rows.append({
        "regime": regime,
        "rows": len(r["counts"]),
        "sum_n_obs": int(r["counts"]["n_obs"].sum()),
        "part_size_mb": r["part_path"].stat().st_size / 1e6,
    })
summary_df = pd.DataFrame(summary_rows).set_index("regime")
print(f"key hexes: {len(key_gdf):,}")
print(f"key.parquet size: {key_path.stat().st_size / 1e6:.2f} MB")
summary_df
```

# Validation

## Key-completeness invariant

Every `release_hex` and `target_hex` in the counts must be in
`key.parquet`. Any violation is a data integrity bug.

```python
key_ids = set(int(x) for x in key_gdf["hex_id"].values)
for regime, r in results.items():
    unseen = (
        set(int(x) for x in r["counts"]["release_hex"].values)
        | set(int(x) for x in r["counts"]["target_hex"].values)
    ) - key_ids
    assert not unseen, (regime, sorted(unseen))
print(f"PASS [r{hex_radius}m]: every release_hex and target_hex is in key.parquet.")
```

## Conservation cross-check

`n_valid` is the count of all obs passing age/hex filters; it must
equal `sum(n_obs)` in the stored counts (no domain filter is applied
in the single-domain build).

```python
for regime, r in results.items():
    n_valid = r["n_valid"]
    n_stored = int(r["counts"]["n_obs"].sum())
    print(f"[{regime}] n_valid={n_valid:,}  stored sum(n_obs)={n_stored:,}")
    assert n_valid == n_stored, (regime, n_valid, n_stored)
```

## Output sizes

```python
total_bytes = key_path.stat().st_size
print(f"{'key.parquet':<60s} {key_path.stat().st_size/1e6:>9.2f} MB")
for regime, r in results.items():
    size = r["part_path"].stat().st_size
    total_bytes += size
    label = f"counts/regime={regime}/release_year={release_year}"
    print(f"{label:<60s} {size/1e6:>9.2f} MB")
print(f"{'TOTAL':<60s} {total_bytes/1e6:>9.2f} MB")
```
