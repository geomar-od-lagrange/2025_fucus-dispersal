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
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Bake Fucus release points

One-shot preprocess: read the Fucus REDLIST shapefile, filter to
*F. vesiculosus* occurrences, reproject to EPSG:4326, and write the
result next to its source under
`helcom_fucus_redlist/fucus_release_points.geojson`. Run once; the
output is committed to the data twin so downstream notebooks read it
directly without re-running this step.

```python tags=["parameters"]
# Read root of the data twin checkout (HELCOM polygons, Fucus shapefile,
# BSH static + coastlines, CMEMS Stokes sample).
data_root = "../data"
```

```python
from pathlib import Path
import geopandas as gpd

data_root = Path(data_root)
```

```python
# Layout assumption: the Fucus REDLIST shapefile lives under
# data_root/helcom_fucus_redlist/REDLIST_SIS_Macrophytes.shp, and the
# derived release-points geojson is co-located with it in the same dir.
gdf = gpd.read_file(data_root / "helcom_fucus_redlist" / "REDLIST_SIS_Macrophytes.shp")
gdf = gdf.loc[gdf.F_vesiculo != 0].to_crs(epsg=4326)
out_path = data_root / "helcom_fucus_redlist" / "fucus_release_points.geojson"
gdf.to_file(out_path, driver="GeoJSON")
print(f"Written: {out_path}")
```

```python
# Validate the written file
gdf_check = gpd.read_file(out_path)
print("rows:", len(gdf_check))
print("CRS:", gdf_check.crs)
print("bbox:", gdf_check.total_bounds)
```
