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
    display_name: min_data (pixi)
    language: python
    name: min_data
---

# Schleswig-Holstein coast splitting (exploratory)

**Note:** This notebook is exploratory and requires the
`data/Administrative_boundaries/Administrative_boundaries.shp` dataset, which
is not redistributed. It is not on the critical path for the dispersal
pipeline.

Splits the Schleswig-Holstein Baltic coastline into named foerde segments
using hand-placed splitlines, and writes the result to
`output_root/SamplePoints/sh_baltic_coast_split.geojson`.

```python tags=["parameters"]
data_root = "../data"
output_root = "../output"
```

```python
import numpy as np
from pathlib import Path

import geopandas as gpd
import shapely
from shapely.ops import split

import matplotlib.pyplot as plt
import cartopy.crs as ccrs

data_root = Path(data_root)
output_root = Path(output_root)
```

```python
administrative_borders = gpd.read_file(
    data_root / "Administrative_boundaries" / "Administrative_boundaries.shp"
).to_crs(crs=ccrs.Geodetic())
sh_border = administrative_borders.set_index("NAME").loc[["Schleswig-Holstein"]]
```

```python
subbasins_lvl2_file = data_root / "helcom_subbasins_2022" / "HELCOM_subbasins_2022_level2.shp"
gdf_subbasins_area = gpd.read_file(subbasins_lvl2_file).to_crs(crs=ccrs.Geodetic())
```

```python
gdf_baltic = gdf_subbasins_area.clip(mask=(0, 0, 12, 55)).cx[:11.5, :54.6].dissolve()
coast_region_radius_deg = 0.01
gdf_baltic_unbuffed = gdf_baltic.buffer(coast_region_radius_deg).buffer(-coast_region_radius_deg)
gdf_baltic_coast = gdf_baltic.boundary.buffer(coast_region_radius_deg)
gdf_baltic_coast.plot()
```

```python
gdf_baltic_buffed = gdf_baltic.buffer(0.2)
sh_coast = (
    sh_border
    .boundary
    .clip(gdf_baltic_buffed)
    .explode()
)
sh_coast.plot()
```

```python
sh_coast[1]
```

```python
sh_main_coast = sh_coast[1].difference(
    administrative_borders.cx[10:11, 53.5:55].dissolve().buffer(-0.001)
)[0]

splitlines = shapely.MultiLineString([
    ((10.8, 54.14), (10.82, 54.04)),
    ((11.11, 54.38), (11.25, 54.52)),
    ((10.3, 54.4), (10.3, 54.5)),
    ((9.9, 54.53), (10.1, 54.55)),
    ((9.9, 54.7), (10, 54.8)),
    ((10, 54.4), (10.2, 54.5)),
])
sh_coast_split = gpd.GeoDataFrame(
    dict(
        names=[
            "Flensburger Foerde", "Schleimuendung", "Eckernfoerde",
            "Kieler Foerde", "Kieler Bucht", "Mecklenburger Bucht",
            "Luebecker Bucht", "Fehmarn",
        ],
        geometry=np.append(
            shapely.get_parts(split(sh_main_coast, splitlines)),
            sh_coast.geometry[0],
        ),
    ),
).explode().reset_index(drop=True).set_index("names")
```

```python
fig, ax = plt.subplots()
sh_border.to_crs(crs=ccrs.Geodetic()).plot(ax=ax, alpha=0.6)
sh_coast_split.plot(cmap="Paired", ax=ax, lw=3)
plt.show()
```

```python
output_path = output_root / "SamplePoints"
output_path.mkdir(parents=True, exist_ok=True)
sh_coast_split.to_file(output_path / "sh_baltic_coast_split.geojson", driver="GeoJSON")
print(f"Written: {output_path / 'sh_baltic_coast_split.geojson'}")
```
