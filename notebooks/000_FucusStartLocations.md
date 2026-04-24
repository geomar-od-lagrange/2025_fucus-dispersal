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

# Bake Fucus release points

One-shot preprocess: read the Fucus REDLIST shapefile, filter to
*F. vesiculosus* occurrences, reproject to EPSG:4326, and write the
result to `data/derived/fucus_release_points.geojson`. Run once; the
output is committed to the data twin so downstream notebooks read it
directly without re-running this step.

```python tags=["parameters"]
data_root = "../data"
```

```python
from pathlib import Path
import geopandas as gpd

data_root = Path(data_root)
derived_dir = data_root / "derived"
derived_dir.mkdir(parents=True, exist_ok=True)
```

```python
gdf = gpd.read_file(data_root / "fucus_redlist_shapefile" / "REDLIST_SIS_Macrophytes.shp")
gdf = gdf.loc[gdf.F_vesiculo != 0].to_crs(epsg=4326)
out_path = derived_dir / "fucus_release_points.geojson"
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
