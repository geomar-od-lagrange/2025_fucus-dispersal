# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: py:percent,md,ipynb
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Build hex-aggregate key file
#
# Per-hex geometry and attributes for the BSH domain at one `hex_radius`.
# Single-process, depends only on static inputs (BSH H0, coastline
# geojsons, Fucus shapefile, HELCOM polygons). Run once per radius.
#
# Counts partitions are built separately by `024_BuildHexAggregates.md`,
# which loads this key as a hard prerequisite (geometry + projection
# metadata).

# %%
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry import box
from shapely.ops import unary_union

from hextraj import HexProj

# %% [markdown]
# # Parameters

# %% tags=["parameters"]
# Read root of the data twin (HELCOM polygons, Fucus shapefile, BSH
# coastline + static H0/lonlat under bsh_hbmnoku_static/).
data_root = "../data"
# Write root of the hex-aggregate store (key file lives at the top).
output_root = "../output"

# Hex radius (corner-to-centre distance) in metres. Sweep this via
# papermill to build keys at multiple radii.
hex_radius = 6000

# Shoreward snap distance for the HELCOM subbasin assignment, in hex
# radii. A hex that overlaps no subbasin polygon but carries land takes
# the nearest subbasin within this distance; 0 disables the snap and
# leaves the assignment on overlap alone. Expressed in radii rather than
# metres so the cap follows `hex_radius` through a sweep.
subbasin_snap_hex_radii = 1.0

# %% [markdown]
# # Derived layout

# %%
data_root = Path(data_root)
output_root = Path(output_root)
store_root = output_root / "HexAggregates"
store_root.mkdir(parents=True, exist_ok=True)

key_path = store_root / f"HexAgg_key_r{hex_radius}m.parquet"
meta_path = key_path.with_suffix(".json")
subbasin_snap_m = subbasin_snap_hex_radii * hex_radius
print(f"key:  {key_path}")
print(f"meta: {meta_path}")
print(f"subbasin shoreward snap: {subbasin_snap_m:,.0f} m")

# %% [markdown]
# # Static inputs

# %%
wet_gdf = gpd.read_file(data_root / "bsh_hbmnoku_static/coastline.geojson")
always_wet_gdf = gpd.read_file(
    data_root / "bsh_hbmnoku_static/coastline_always_wet.geojson"
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

# %% [markdown]
# # HexProj and BSH-domain hex set
#
# Domain centroid from the coarse-grid H0. The hex domain is every hex
# intersecting the coarse-grid bbox padded by one coarse cell on each
# side — broader than the H0-cell pass because trajectories sample at
# sub-grid positions and can visit hexes that contain no H0 point (yet
# still have water in the wet polygon). H0 is then attached as a per-hex
# mean depth where available, NaN elsewhere.

# %%
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

# %% [markdown]
# # Fine-precedes-coarse coastline merge
#
# The coastline geojsons store fine and coarse staircase polygons in the
# same file with overlap in the German Bight. A plain `unary_union` over
# both wets pixels that fine resolves as land. The right merge is
# `fine_polys ∪ (coarse_polys \ fine_bbox)`, where `fine_bbox`
# is the H0-fine bbox padded by half a fine cell so all fine cells fall
# inside.

# %%
fine_lon_pad = abs(float(h0_fine.lon.diff("lon").mean())) / 2
fine_lat_pad = abs(float(h0_fine.lat.diff("lat").mean())) / 2
fine_bbox = box(
    float(h0_fine.lon.min()) - fine_lon_pad,
    float(h0_fine.lat.min()) - fine_lat_pad,
    float(h0_fine.lon.max()) + fine_lon_pad,
    float(h0_fine.lat.max()) + fine_lat_pad,
)


def fine_first_union(gdf):
    fine = unary_union(gdf.loc[gdf["grid"] == "fine", "geometry"].tolist())
    coarse = unary_union(gdf.loc[gdf["grid"] == "coarse", "geometry"].tolist())
    return unary_union([fine, coarse.difference(fine_bbox)])


wet_union_3035 = (
    gpd.GeoSeries([fine_first_union(wet_gdf)], crs=4326).to_crs(3035).iloc[0]
)
always_wet_union_3035 = (
    gpd.GeoSeries([fine_first_union(always_wet_gdf)], crs=4326).to_crs(3035).iloc[0]
)


# %%
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
# is absent. NaN for hexes without H0 coverage.
mean_depth = (
    h0_fine_frame.groupby("hex_id")["H0"].mean()
    .combine_first(h0_coarse_frame.groupby("hex_id")["H0"].mean())
)

# abs() because `lat` is stored descending — `diff().max()` is negative
# in that case and would shrink the bbox instead of padding it.
lon_step = abs(float(h0_coarse.lon.diff("lon").mean()))
lat_step = abs(float(h0_coarse.lat.diff("lat").mean()))
hex_ids = hp.rectangle_of_hexes(
    float(h0_coarse.lon.min()) - lon_step,
    float(h0_coarse.lon.max()) + lon_step,
    float(h0_coarse.lat.min()) - lat_step,
    float(h0_coarse.lat.max()) + lat_step,
).astype(np.int32)
print(f"BSH-domain hexes at r={hex_radius} m: {len(hex_ids):,}")

# %% [markdown]
# # Key file
#
# Per-hex geometry and attributes (area, water area, mean depth, coast
# distance, Fucus area, HELCOM subbasin).

# %%
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

# HELCOM subbasin per hex, shoreward-extended. The polygons are a marine
# delineation digitised against the coastline, so they stop short of the
# shore: a hex that fronts a basin has most of its water in that basin
# but its centroid on land or in an archipelago gap. Two steps, in
# EPSG:3035 so both the areas and the snap distance are metric.
#
# 1. Largest overlap — the basin holding the most of the hex polygon.
# 2. Shoreward snap — a hex overlapping no basin takes the nearest basin
#    within `subbasin_snap_m`, but only if it carries land. That land
#    test is what keeps the extension landward: the HELCOM outer boundary
#    (Kattegat → Skagerrak, and the whole North Sea side of the BSH
#    domain) runs through open water, so a hex beyond it is all water and
#    never snaps, while a shoreward hex the polygons missed always holds
#    some of the coast that stopped them.
subbasins_3035 = subbasins[["geometry", "subbasin"]].to_crs(3035)
hexes_3035 = gpd.GeoDataFrame({"hex_id": key["hex_id"]}, geometry=geom_3035, crs=3035)

overlap = gpd.overlay(
    hexes_3035, subbasins_3035, how="intersection", keep_geom_type=True
)
largest_overlap = (
    overlap.assign(overlap_m2=overlap.area)
    .sort_values("overlap_m2")
    .groupby("hex_id")["subbasin"].last()
)

carries_land = (key["water_area_m2"] < key["area_m2"]).to_numpy()
snap_pool = hexes_3035[~hexes_3035["hex_id"].isin(largest_overlap.index) & carries_land]
snapped = (
    gpd.sjoin_nearest(
        snap_pool, subbasins_3035, how="inner",
        max_distance=subbasin_snap_m, distance_col="snap_m",
    )
    .sort_values("snap_m")
    .groupby("hex_id")["subbasin"].first()
)

key["helcom_subbasin"] = (
    key["hex_id"]
    .map(largest_overlap.combine_first(snapped).map(subbasin_name_to_id))
    .fillna(-1).astype(np.int8)
)

key = key[["hex_id", "geometry", "area_m2", "water_area_m2", "fucus_area_m2",
            "mean_depth_m", "dist_to_coast_m", "helcom_subbasin"]]
key.head()

# %%
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
    "subbasin_rule": "largest_overlap+shoreward_snap",
    "subbasin_snap_m": subbasin_snap_m,
}
meta_path.write_text(json.dumps(meta, indent=2))

print(f"wrote {key_path} ({key_path.stat().st_size / 1e6:.2f} MB)")
print(f"wrote {meta_path}")

# %% [markdown]
# # Validation

# %%
print(f"hexes:                {len(key):,}")
print(f"  with H0 coverage:   {key['mean_depth_m'].notna().sum():,}")
print(f"  with water area:    {(key['water_area_m2'] > 0).sum():,}")
print(f"  with Fucus area:    {(key['fucus_area_m2'] > 0).sum():,}")
print(f"  inside a subbasin:  {(key['helcom_subbasin'] >= 0).sum():,}")
print(f"    by overlap:       {len(largest_overlap):,}")
print(f"    by shoreward snap:{len(snapped):,}")
print(f"geometry valid:       {key.geometry.is_valid.all()}")

# %% [markdown]
# The residual `-1` class after the shoreward extension is what the HELCOM
# polygons genuinely do not cover: the North Sea and Skagerrak side of the
# BSH domain, plus inland hexes. It is reported split by whether the hex
# holds any model water, because only those can ever carry a particle.

# %%
outside = key["helcom_subbasin"] < 0
print(f"outside (-1):         {int(outside.sum()):,}")
print(f"  with water area:    {int((outside & (key['water_area_m2'] > 0)).sum()):,}")
print(f"  lon range:          {key.loc[outside & (key['water_area_m2'] > 0)].geometry.centroid.x.min():.2f} "
      f"… {key.loc[outside & (key['water_area_m2'] > 0)].geometry.centroid.x.max():.2f}")

# %% [markdown]
# Snap-cap sensitivity: how many further hexes each alternative cap would
# take, over the same land-bearing pool the rule above is restricted to.

# %%
unassigned_pool = hexes_3035[
    ~hexes_3035["hex_id"].isin(largest_overlap.index) & carries_land
]
snap_reach = gpd.sjoin_nearest(
    unassigned_pool, subbasins_3035, how="inner", distance_col="snap_m"
).sort_values("snap_m").drop_duplicates("hex_id")
for radii in (0.5, 1.0, 2.0, 4.0):
    cap = radii * hex_radius
    print(f"  snap ≤ {radii:>3}×r ({cap:>6,.0f} m): "
          f"{int((snap_reach['snap_m'] <= cap).sum()):,} hexes")
