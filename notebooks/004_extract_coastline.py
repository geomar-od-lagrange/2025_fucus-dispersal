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
from shapely.geometry import MultiPolygon


def regular_axis(coord):
    """Step and first centre of a regular axis from its stored centres.

    The BSH static coordinates are stored rounded to 1e-4 deg, so the first
    difference is off by up to 1e-4 deg per cell (0.0834 for 1/12) and an
    origin + i * step built from it drifts across the grid. A linear fit over
    every centre recovers the regular grid the stored values sample; its
    residual is the storage rounding, which the assert bounds at 1 % of a cell.
    """
    index = np.arange(coord.size)
    step, first = np.polyfit(index, coord, 1)
    residual = np.abs(coord - (first + step * index)).max()
    assert residual < 0.01 * abs(step), f"axis not regular: residual {residual:.2e} deg"
    return float(step), float(first)


def coastline_from_h0(h0_path, min_h0=None, scale=6):
    """Build staircase coastline polygon from a BSH H0 file.

    Works in integer grid coordinates (scaled by ``scale``) to get exact
    vertex alignment when unioning cell boxes, then affine-transforms
    back to lon/lat. Cell ``(j, i)`` spans its centre -+ half a step: the
    NE F-points the 003 c_files and 005 wall fields carry are its
    north-east corner.

    If *min_h0* is given, only cells with ``H0 > min_h0`` are included.
    """
    ds_h0 = xr.open_dataset(h0_path)
    h0 = ds_h0.H0
    if min_h0 is not None:
        wet = (h0.notnull() & (h0 > min_h0)).values
    else:
        wet = h0.notnull().values  # (nlat, nlon)

    dlon, lon0 = regular_axis(ds_h0.lon.values)  # dlon > 0
    dlat, lat0 = regular_axis(ds_h0.lat.values)  # dlat < 0: lat descends

    jj, ii = np.nonzero(wet)
    wet_region = shapely.union_all(
        shapely.box(ii * scale, jj * scale, (ii + 1) * scale, (jj + 1) * scale)
    )

    # Integer grid (i, j) -> lon/lat; the cell edge sits half a step before
    # centre 0 along each axis (signed, so it works for descending lat).
    return shapely.affinity.affine_transform(
        wet_region,
        [dlon / scale, 0, 0, dlat / scale, lon0 - dlon / 2, lat0 - dlat / 2],
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
