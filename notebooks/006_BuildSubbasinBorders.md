---
jupyter:
  jupytext:
    cell_metadata_filter: tags,-all
    formats: py:percent,md,ipynb
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

# Subbasin borders on the model coastline

The HELCOM level-2 subbasin borders every map of the pipeline draws, built
once from static inputs and written as
`output_root/SubbasinBorders/subbasin_borders.geojson` (one row per pair of
adjacent subbasins). The model coastline is the reference: every border is
continuous and ends exactly on the BSH HBMnoku coastline the maps draw.

The borders are the edges of the **nearest-HELCOM-polygon partition of model
water**. Inside the HELCOM polygons that partition is HELCOM's own; where
the model coast and the HELCOM coast disagree (archipelagos, coarse-grid
coasts) it is the line equidistant from the two polygons, the border's
continuation to the model coast. Construction (docs/visualisations.md
"Subbasin borders"):

1. Rasterise the HELCOM polygons on a `pixel_m` grid in EPSG:3035 (holes,
   the islands, left unlabelled) and unlabel the components smaller than
   `min_component_km2` — archipelago slivers that would otherwise draw
   enclave loops.
2. Give every unlabelled pixel the label of the nearest labelled one, up to
   `reach_m`; beyond it stays unlabelled, so no border is drawn across
   Jutland into the North Sea, where no two subbasins meet.
3. Pixel edges between two labels, merged per pair and Douglas–Peucker
   simplified at `simplify_pixels` (the raster staircase goes).
4. Intersect with the model water polygon in EPSG:4326, so each line stops
   on the model coastline and none runs over model land.

```python
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
import xarray as xr
from PIL import Image, ImageDraw
from scipy import ndimage
from shapely.geometry import box
from shapely.ops import linemerge, unary_union
```

# Parameters

```python tags=["parameters"]
# Read root of the data twin (HELCOM polygons, BSH coastline + static H0).
data_root = "../data"
# Write root; the borders land under SubbasinBorders/.
output_root = "../output"
# Raster cell (m, EPSG:3035) of the partition; the border position is good to
# about this much, far below a line width at page scale.
pixel_m = 200.0
# Label components below this area (km²) are archipelago slivers and are
# refilled from their neighbours.
min_component_km2 = 25.0
# Nearest-polygon fill reach (m): covers every gap between the HELCOM and the
# model coast inside the Baltic, stops short of the North Sea.
reach_m = 25_000.0
# Douglas–Peucker tolerance, in pixels.
simplify_pixels = 1.5
```

```python
data_root = Path(data_root)
output_root = Path(output_root)
out_path = output_root / "SubbasinBorders" / "subbasin_borders.geojson"
out_path.parent.mkdir(parents=True, exist_ok=True)
print(f"borders: {out_path}")
```

# Rasterised HELCOM polygons

```python
subbasins = gpd.read_file(
    data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp"
).to_crs(3035)
names = subbasins["level_2"].tolist()
# Grid over the polygons plus the fill reach, so the fill is never cut by the
# raster edge.
x_min, y_min, x_max, y_max = subbasins.total_bounds + np.array(
    [-reach_m, -reach_m, reach_m, reach_m]
)
nx = int(np.ceil((x_max - x_min) / pixel_m))
ny = int(np.ceil((y_max - y_min) / pixel_m))
print(f"grid {nx} x {ny} at {pixel_m:g} m")


def to_pixels(ring):
    """Ring coordinates (EPSG:3035) as PIL pixel coordinates, row 0 at y_max."""
    xy = np.asarray(ring.coords)[:, :2]
    return list(map(tuple, np.column_stack(
        [(xy[:, 0] - x_min) / pixel_m, (y_max - xy[:, 1]) / pixel_m]
    )))


# Label 0 = unlabelled, i + 1 = subbasin i. A pixel two polygons both cover
# keeps the first.
label = np.zeros((ny, nx), np.uint8)
for i, geom in enumerate(subbasins.geometry, start=1):
    mask = Image.new("1", (nx, ny), 0)
    draw = ImageDraw.Draw(mask)
    for part in shapely.get_parts(geom):
        draw.polygon(to_pixels(part.exterior), fill=1)
    for part in shapely.get_parts(geom):
        for ring in part.interiors:
            draw.polygon(to_pixels(ring), fill=0)
    label[np.array(mask, bool) & (label == 0)] = i

min_component_px = int(min_component_km2 * 1e6 / pixel_m**2)
for i, name in enumerate(names, start=1):
    components, _ = ndimage.label(label == i)
    sizes = np.bincount(components.ravel())
    small = np.flatnonzero(sizes < min_component_px)
    small = small[small > 0]
    label[np.isin(components, small)] = 0
    print(f"{name}: {len(small)} components < {min_component_km2:g} km² refilled")
```

# Nearest-polygon partition

```python
distance, nearest = ndimage.distance_transform_edt(
    label == 0, return_distances=True, return_indices=True
)
partition = label[nearest[0], nearest[1]]
partition[distance * pixel_m > reach_m] = 0
del distance, nearest
```

# Partition edges per subbasin pair

One unit segment per pixel edge between two labels; edges against
unlabelled pixels (outside the reach) are not borders.

```python
segments = {}
rows, cols = np.nonzero(partition[:, 1:] != partition[:, :-1])
for r, c in zip(rows, cols):
    pair = tuple(sorted((int(partition[r, c]), int(partition[r, c + 1]))))
    x = x_min + (c + 1) * pixel_m
    segments.setdefault(pair, []).append(
        ((x, y_max - r * pixel_m), (x, y_max - (r + 1) * pixel_m))
    )
rows, cols = np.nonzero(partition[1:, :] != partition[:-1, :])
for r, c in zip(rows, cols):
    pair = tuple(sorted((int(partition[r, c]), int(partition[r + 1, c]))))
    y = y_max - (r + 1) * pixel_m
    segments.setdefault(pair, []).append(
        ((x_min + c * pixel_m, y), (x_min + (c + 1) * pixel_m, y))
    )
partition_edges = gpd.GeoDataFrame(
    [
        {
            "subbasin_a": names[a - 1],
            "subbasin_b": names[b - 1],
            "geometry": linemerge(shapely.multilinestrings(s)).simplify(
                simplify_pixels * pixel_m
            ),
        }
        for (a, b), s in segments.items()
        if a > 0
    ],
    crs="EPSG:3035",
).to_crs(4326)
print(f"{len(partition_edges)} subbasin pairs meet in the partition")
```

# Clip to model water

The model water polygon is the maps' own: the BSH wet cells
(`bsh_hbmnoku_static/coastline.geojson`), fine nest first inside its
bounding box (024a's merge). Its outline is the coastline every map draws,
so a clipped border ends on it.

```python
h0_fine = xr.open_dataset(
    data_root / "bsh_hbmnoku_static" / "static_file_fine" / "H0_file_fine.nc"
).H0
pad_lon = 0.5 * abs(float(h0_fine.lon.diff("lon").mean()))
pad_lat = 0.5 * abs(float(h0_fine.lat.diff("lat").mean()))
fine_bbox = box(
    float(h0_fine.lon.min()) - pad_lon, float(h0_fine.lat.min()) - pad_lat,
    float(h0_fine.lon.max()) + pad_lon, float(h0_fine.lat.max()) + pad_lat,
)
wet = gpd.read_file(data_root / "bsh_hbmnoku_static" / "coastline.geojson")
water = unary_union([
    unary_union(wet.loc[wet["grid"] == "fine", "geometry"].tolist()),
    unary_union(wet.loc[wet["grid"] == "coarse", "geometry"].tolist()).difference(fine_bbox),
])

borders = partition_edges.assign(
    geometry=partition_edges.geometry.intersection(water).line_merge()
)
borders = borders[~borders.geometry.is_empty].reset_index(drop=True)
print(f"{len(borders)} subbasin pairs share a border in model water")
```

# Checks

Every line end lies on the model coastline or on another border (a
junction), and no border length lies on model land.

```python
parts = shapely.get_parts(borders.geometry.to_numpy())
coastline = water.boundary
dangling = 0
for i, part in enumerate(parts):
    others = unary_union(np.delete(parts, i))
    for end in shapely.get_point(part, [0, -1]):
        if end.distance(coastline) > 1e-9 and end.distance(others) > 1e-9:
            dangling += 1
on_land_m = (
    gpd.GeoSeries([unary_union(parts).difference(water)], crs="EPSG:4326")
    .to_crs(3035).length.sum()
)
print(f"{len(parts)} line parts, {dangling} dangling ends, {on_land_m:.0f} m on model land")
assert dangling == 0, f"{dangling} border ends neither on the coast nor at a junction"
assert on_land_m < 1, f"{on_land_m:.0f} m of border on model land"
print(borders[["subbasin_a", "subbasin_b"]].to_string())
```

```python
borders.to_file(out_path, driver="GeoJSON")
print(f"wrote {out_path}")
```
