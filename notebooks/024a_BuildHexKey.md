---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: md,ipynb
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.19.1
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Build hex-aggregate key file

Per-hex geometry and attributes for the BSH domain at one `hex_radius`.
Single-process, depends only on static inputs (BSH H0, coastline
geojsons, Fucus shapefile, HELCOM polygons). Run once per radius.

Counts partitions are built separately by `024_BuildHexAggregates.md`,
which loads this key as a hard prerequisite (geometry + projection
metadata).

```python
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.ops import unary_union

from hextraj import HexProj
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (HELCOM polygons, Fucus shapefile, BSH
# coastline + static H0/lonlat under bsh_hbmnoku_static/).
data_root = "../data"
# Write root of the hex-aggregate store (key file lives at the top).
output_root = "../output"

# Hex radius (corner-to-centre distance) in metres. Sweep this via
# papermill to build keys at multiple radii.
hex_radius = 6000
```

# Derived layout

```python
data_root = Path(data_root)
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
print(f"key:  {key_path}")
print(f"meta: {meta_path}")
```

# Static inputs

```python
wet_gdf = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
always_wet_gdf = gpd.read_file(
    data_root / "bsh_hbmnoku_static/coastline_always_wet.geojson"
)
wet_union_3035 = (
    gpd.GeoSeries([unary_union(wet_gdf.geometry)], crs=4326).to_crs(3035).iloc[0]
)
always_wet_union_3035 = (
    gpd.GeoSeries([unary_union(always_wet_gdf.geometry)], crs=4326).to_crs(3035).iloc[0]
)

fucus_union_3035 = (
    gpd.read_file(data_root / "helcom_fucus_redlist/REDLIST_SIS_Macrophytes.shp")
    .loc[lambda df: df.F_vesiculo != 0, ["geometry"]]
    .to_crs(epsg=3035)
    .geometry.unary_union
)

subbasins = (
    gpd.read_file(
        data_root / "helcom_subbasins_2022/HELCOM_subbasins_2022_level2.shp"
    )
    .to_crs(epsg=4326)
    .rename(columns={"level_2": "subbasin"})
    .reset_index(drop=True)
)
# Stable int8 lookup for the sidecar metadata. -1 reserved for "outside".
subbasin_name_to_id = {str(name): i for i, name in enumerate(subbasins["subbasin"])}
subbasin_id_to_name = {-1: "_outside", **{i: n for n, i in subbasin_name_to_id.items()}}
```

# HexProj and BSH-domain hex set

Domain centroid from the coarse-grid H0; label every fine+coarse H0
grid point to enumerate every hex `hp.label` can ever assign within the
BSH lat/lon extent.

```python
h0_coarse = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_coarse/H0_file_coarse.nc"
)
h0_fine = xr.open_dataset(
    data_root / "bsh_hbmnoku_static/static_file_fine/H0_file_fine.nc"
)
domain_lon_origin = float(0.5 * (h0_coarse.lon.min() + h0_coarse.lon.max()))
domain_lat_origin = float(0.5 * (h0_coarse.lat.min() + h0_coarse.lat.max()))
print(f"BSH domain centroid: lon={domain_lon_origin:.4f}, lat={domain_lat_origin:.4f}")

hp = HexProj(
    projection_name="laea",
    lon_origin=domain_lon_origin,
    lat_origin=domain_lat_origin,
    hex_size_meters=hex_radius,
)
```

```python
def h0_frame(h0, grid_name):
    """Long-form H0 cells (lon, lat, H0 > 0) labelled by hex_id."""
    lon2d, lat2d = np.meshgrid(h0.lon.values, h0.lat.values)
    df = pd.DataFrame({
        "lon": lon2d.ravel(),
        "lat": lat2d.ravel(),
        "H0": h0.H0.values.ravel(),
    })
    df = df[(df.H0 > 0) & np.isfinite(df.H0)]
    df["hex_id"] = hp.label(df.lon.values, df.lat.values)
    return df.loc[df.hex_id >= 0, ["hex_id", "H0"]].assign(grid=grid_name)


h0_fine_frame = h0_frame(h0_fine, "fine")
h0_coarse_frame = h0_frame(h0_coarse, "coarse")

# Mean depth, fine-grid priority: combine_first fills coarse where fine
# is absent. The hex domain is the union of both.
mean_depth = (
    h0_fine_frame.groupby("hex_id")["H0"].mean()
    .combine_first(h0_coarse_frame.groupby("hex_id")["H0"].mean())
)
hex_ids = mean_depth.index.to_numpy(dtype=np.int32)
print(f"BSH-domain hexes at r={hex_radius} m: {len(hex_ids):,}")
```

# Key file

Per-hex geometry and attributes (area, water area, mean depth, coast
distance, Fucus area, HELCOM subbasin).

```python
key = (
    hp.to_geodataframe(hex_ids.tolist())
    .rename_axis("hex_id").reset_index()
    .assign(hex_id=lambda df: df["hex_id"].astype(np.int32))
    .set_crs(epsg=4326, allow_override=True)
)

geom_3035 = key.geometry.to_crs(3035)
centroids_3035 = key.geometry.centroid.to_crs(3035)

key["area_m2"] = geom_3035.area.astype(np.float32).values
key["water_area_m2"] = (
    geom_3035.intersection(wet_union_3035).area.astype(np.float32).values
)
key["fucus_area_m2"] = (
    geom_3035.intersection(fucus_union_3035).area.astype(np.float32).values
)
key["mean_depth_m"] = (
    key["hex_id"].map(mean_depth).astype(np.float32).values
)
key["dist_to_coast_m"] = (
    centroids_3035.distance(always_wet_union_3035.boundary)
    .astype(np.float32).values
)

# HELCOM subbasin by centroid (in 4326).
joined = gpd.sjoin(
    gpd.GeoDataFrame({"hex_id": key["hex_id"]}, geometry=key.geometry.centroid, crs=4326),
    subbasins[["geometry", "subbasin"]],
    how="left", predicate="within",
).drop_duplicates(subset="hex_id").set_index("hex_id")
key["helcom_subbasin"] = (
    key["hex_id"].map(joined["subbasin"].map(subbasin_name_to_id))
    .fillna(-1).astype(np.int8)
)

key = key[["hex_id", "geometry", "area_m2", "water_area_m2", "fucus_area_m2",
            "mean_depth_m", "dist_to_coast_m", "helcom_subbasin"]]
key.head()
```

```python
key.to_parquet(key_path)

meta = {
    "hex_proj": {
        "projection_name": "laea",
        "lon_origin": domain_lon_origin,
        "lat_origin": domain_lat_origin,
        "hex_size_meters": hex_radius,
    },
    "area_crs": "EPSG:3035",
    "subbasin_id_to_name": {str(k): v for k, v in subbasin_id_to_name.items()},
}
meta_path.write_text(json.dumps(meta, indent=2))

print(f"wrote {key_path} ({key_path.stat().st_size / 1e6:.2f} MB)")
print(f"wrote {meta_path}")
```

# Validation

```python
print(f"hexes:                {len(key):,}")
print(f"  with H0 coverage:   {key['mean_depth_m'].notna().sum():,}")
print(f"  with water area:    {(key['water_area_m2'] > 0).sum():,}")
print(f"  with Fucus area:    {(key['fucus_area_m2'] > 0).sum():,}")
print(f"  inside a subbasin:  {(key['helcom_subbasin'] >= 0).sum():,}")
print(f"geometry valid:       {key.geometry.is_valid.all()}")
```
