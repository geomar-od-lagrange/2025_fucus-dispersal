"""Shared helpers for the visualisation notebooks.

Everything here stays lazy on the (trajectory, obs) arrays: no ``.compute()``
on the large arrays. Downstream notebooks trigger one ``dask.compute(*...)``
per scope group instead of re-walking the graph per plot.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import shapely
import xarray as xr
from scipy.ndimage import distance_transform_edt


QUARTER_LABELS = {1: "JFM", 2: "AMJ", 3: "JAS", 4: "OND"}


def load_trajectories(trajectory_path):
    """Lazy-concat all ``*.zarr`` files under ``trajectory_path`` along
    ``trajectory``. Returns ``(ds, zarr_files)``."""
    zarr_files = sorted(Path(trajectory_path).glob("**/*.zarr"))
    ds = xr.concat([xr.open_zarr(z) for z in zarr_files], dim="trajectory")
    return ds, zarr_files


def mask_land_seeded(ds):
    """Return ``(masked_ds, land_seeded)``.

    ``land_seeded`` is a lazy bool along ``trajectory``; trajectories whose
    first step is exactly zero displacement are NaN-filled across all obs.
    """
    dlon0 = ds.lon.diff("obs").isel(obs=0, drop=True)
    dlat0 = ds.lat.diff("obs").isel(obs=0, drop=True)
    land_seeded = (dlon0 == 0) & (dlat0 == 0)
    return ds.where(~land_seeded), land_seeded


def _build_subbasin_raster(
    subbasins, lon_min, lon_max, lat_min, lat_max, resolution, fill_max_pixels
):
    """Raster-Voronoi of HELCOM subbasins.

    Returns ``(filled, lon_edges, lat_edges, id_to_name)``.

    ``filled`` is an int raster on a regular grid; 0 = unassigned (outside the
    fill cap). Interior cells carry the ID of the subbasin they sit in;
    exterior cells within ``fill_max_pixels`` of a polygon inherit from the
    nearest polygon cell — a raster Voronoi that respects polygon shapes and
    is overlap-free by construction.
    """
    lon_edges = np.arange(lon_min, lon_max + resolution, resolution)
    lat_edges = np.arange(lat_min, lat_max + resolution, resolution)
    lon_centres = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_centres = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    xx, yy = np.meshgrid(lon_centres, lat_centres)  # (n_lat, n_lon)

    names = list(subbasins.subbasin.values)
    id_to_name = np.array([None] + list(names), dtype=object)
    raster = np.zeros(xx.shape, dtype=np.int32)
    for i, geom in enumerate(subbasins.geometry, start=1):
        mask = shapely.contains_xy(geom, xx, yy)
        # HELCOM level-2 tiles the Baltic; gaps only (no overlaps expected).
        # If two polygons claim the same pixel, first-come wins.
        raster[mask & (raster == 0)] = i

    dist, (iy, ix) = distance_transform_edt(
        raster == 0, return_distances=True, return_indices=True
    )
    filled = raster[iy, ix]
    filled[dist > fill_max_pixels] = 0
    return filled, lon_edges, lat_edges, id_to_name


def attach_release_metadata(
    ds,
    subbasins,
    lon_min=5,
    lon_max=32,
    lat_min=53,
    lat_max=66,
    raster_resolution=0.02,
    fill_max_pixels=10,
):
    """Attach ``release_quarter`` (int 1-4, JFM/AMJ/JAS/OND) and ``subbasin``
    per trajectory.

    Subbasin comes from a raster-Voronoi lookup on release points
    (``obs=0`` lon/lat). Points inside any HELCOM polygon hit it directly;
    points just outside inherit from the nearest polygon cell (up to
    ``fill_max_pixels`` away ≈ ``fill_max_pixels * raster_resolution * 111`` km).
    Points further out (or NaN-seeded) become ``None``.

    The lookup runs inside ``xr.apply_ufunc(..., dask="parallelized")`` on the
    1-D ``(trajectory,)`` release coords, so shapely only ever sees polygons
    (once, at graph-build time), never trajectories.
    """
    release_time = ds.time.isel(obs=0, drop=True)
    ds = ds.assign(release_quarter=release_time.dt.quarter)

    filled, lon_edges, lat_edges, id_to_name = _build_subbasin_raster(
        subbasins, lon_min, lon_max, lat_min, lat_max,
        raster_resolution, fill_max_pixels,
    )

    def _lookup(lon, lat):
        out = np.full(lon.shape, None, dtype=object)
        valid = ~(np.isnan(lon) | np.isnan(lat))
        if not valid.any():
            return out
        lon_v = lon[valid]
        lat_v = lat[valid]
        lat_idx = np.searchsorted(lat_edges, lat_v) - 1
        lon_idx = np.searchsorted(lon_edges, lon_v) - 1
        in_bounds = (
            (lat_idx >= 0) & (lat_idx < filled.shape[0])
            & (lon_idx >= 0) & (lon_idx < filled.shape[1])
        )
        ids = np.zeros(lat_idx.size, dtype=np.int32)
        ib = in_bounds
        ids[ib] = filled[lat_idx[ib], lon_idx[ib]]
        out[valid] = np.where(ids > 0, id_to_name[ids], None)
        return out

    lon0 = ds.lon.isel(obs=0, drop=True)
    lat0 = ds.lat.isel(obs=0, drop=True)
    subbasin = xr.apply_ufunc(
        _lookup, lon0, lat0,
        dask="parallelized",
        output_dtypes=[object],
    )
    return ds.assign(subbasin=subbasin)


def relabel_quarter(da, dim="release_quarter"):
    """Replace int 1..4 quarter coord values with JFM/AMJ/JAS/OND labels."""
    return da.assign_coords(
        {dim: [QUARTER_LABELS[int(q)] for q in da[dim].values]}
    )


# Pattern: Fucus_BSH_YYYYMMDD_{regime}_dt{N}min
# Regime is one of the three known kernel-forcing variants; surface_stokes
# must be listed before surface so the alternation matches the longer form
# first.
_ZARR_STEM_RE = re.compile(
    r"^Fucus_BSH_(\d{8})_(surface_stokes|surface|bottom)_dt\d+min$"
)


def parse_zarr_stem(path):
    """Parse a trajectory zarr filename into its (release_date, regime) tuple.

    The authoritative filename format (from notebooks/010_FucusDispersal.ipynb) is::

        Fucus_BSH_{YYYYMMDD}_{regime}_dt{N}min.zarr

    where ``{regime}`` is one of ``surface``, ``surface_stokes``, or ``bottom``.

    Returns ``(release_date: pandas.Timestamp, regime: str)``.
    Raises ``ValueError`` with the offending filename on mismatch.
    """
    path = Path(path)
    stem = path.stem  # strips the final .zarr suffix
    m = _ZARR_STEM_RE.match(stem)
    if m is None:
        raise ValueError(
            f"zarr filename does not match expected pattern "
            f"'Fucus_BSH_YYYYMMDD_<regime>_…': {path.name!r}"
        )
    release_date = pd.Timestamp(m.group(1))  # YYYYMMDD → Timestamp
    regime = m.group(2)
    return release_date, regime
