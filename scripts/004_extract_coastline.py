"""Extract staircase coastline polygons from BSH H0 (bathymetry) files.

Reads the fine and coarse H0 files, builds a rectangle for each wet tracer
cell, unions them in integer grid space (to avoid floating-point slivers),
transforms back to lon/lat, and writes the combined MultiPolygon as GeoJSON.

Usage:
    python 004_extract_coastline.py \\
        --bsh-root /path/to/bsh_operationalmodel_data \\
        --output-geojson /path/to/data/bsh_hbmnoku_static/coastline.geojson

    # Exclude tidal flats (H0 < 0) for an "always wet" coastline:
    python 004_extract_coastline.py \\
        --bsh-root /path/to/bsh_operationalmodel_data \\
        --output-geojson /path/to/data/bsh_hbmnoku_static/coastline_always_wet.geojson \\
        --min-h0 0
"""

import argparse
from pathlib import Path

import geopandas as gpd
import numpy as np
import shapely
import xarray as xr
from shapely.geometry import box, MultiPolygon


def coastline_from_h0(h0_path, min_h0=None, scale=6):
    """Build staircase coastline polygon from a BSH H0 file.

    Works in integer grid coordinates (scaled by ``scale``) to get exact
    vertex alignment when unioning cell boxes, then affine-transforms
    back to lon/lat.

    If *min_h0* is given, only cells with ``H0 > min_h0`` are included.
    """
    ds_h0 = xr.open_dataset(h0_path)
    lon = ds_h0.lon.values
    lat = ds_h0.lat.values
    h0 = ds_h0.H0
    if min_h0 is not None:
        wet = (h0.notnull() & (h0 > min_h0)).values
    else:
        wet = h0.notnull().values  # (nlat, nlon)

    dlon = float(lon[1] - lon[0])  # positive
    dlat = float(lat[1] - lat[0])  # negative if lat is descending

    boxes = []
    for j in range(len(lat)):
        for i in range(len(lon)):
            if wet[j, i]:
                boxes.append(box(i * scale, j * scale, (i + 1) * scale, (j + 1) * scale))

    wet_region = shapely.unary_union(boxes)

    # Map integer grid back to lon/lat.
    # j=0 corresponds to lat[0], i=0 corresponds to lon[0].
    # Cell edges sit half a grid step outside the cell centre.
    lon_origin = float(lon[0]) - abs(dlon) / 2
    lat_origin = float(lat[0]) - dlat / 2  # signed: shifts toward lower j

    return shapely.affinity.affine_transform(
        wet_region,
        [abs(dlon) / scale, 0, 0, dlat / scale, lon_origin, lat_origin],
    )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bsh-root", type=Path, required=True,
        help="BSH store root (contains static_file_fine/ and static_file_coarse/)",
    )
    parser.add_argument(
        "--output-geojson", type=Path, required=True,
        help="Output GeoJSON path",
    )
    parser.add_argument(
        "--min-h0", type=float, default=None,
        help="Minimum H0 value to include (e.g. 0 to exclude tidal flats)",
    )
    args = parser.parse_args()

    grids = {
        "fine": args.bsh_root / "static_file_fine" / "H0_file_fine.nc",
        "coarse": args.bsh_root / "static_file_coarse" / "H0_file_coarse.nc",
    }

    all_polys = []
    sources = []
    for label, h0_path in grids.items():
        print(f"Processing {label}: {h0_path}")
        region = coastline_from_h0(h0_path, min_h0=args.min_h0)
        if isinstance(region, MultiPolygon):
            for poly in region.geoms:
                all_polys.append(poly)
                sources.append(label)
        else:
            all_polys.append(region)
            sources.append(label)
        print(f"  {region.geom_type}, {len(getattr(region, 'geoms', [region]))} parts")

    gdf = gpd.GeoDataFrame({"grid": sources}, geometry=all_polys, crs="EPSG:4326")
    args.output_geojson.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(args.output_geojson, driver="GeoJSON")
    print(f"\nWrote {len(gdf)} polygons to {args.output_geojson}")


if __name__ == "__main__":
    main()
