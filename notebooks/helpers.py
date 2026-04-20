"""Shared helpers for the visualisation notebooks.

Everything here stays lazy on the (trajectory, obs) arrays: no ``.compute()``
on the large arrays. Downstream notebooks trigger one ``dask.compute(*...)``
per scope group instead of re-walking the graph per plot.
"""

from pathlib import Path

import dask.array as da
import geopandas as gpd
import numpy as np
import xarray as xr


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


def attach_release_metadata(ds, subbasins):
    """Attach ``release_year``, ``release_quarter`` (int 1-4, JFM/AMJ/JAS/OND),
    and ``subbasin`` per trajectory.

    Subbasin comes from a spatial nearest-neighbour join on release points
    (first-obs lon/lat). NaN release points (land-seeded trajectories) get
    NaN subbasin via reindex.
    """
    release_time = ds.time.isel(obs=0, drop=True)
    ds = ds.assign(
        release_year=release_time.dt.year,
        release_quarter=release_time.dt.quarter,
    )

    release_lon = ds.lon.isel(obs=0, drop=True).to_pandas()
    release_lat = ds.lat.isel(obs=0, drop=True).to_pandas()
    valid = release_lon.notna() & release_lat.notna()
    release_pts = gpd.GeoDataFrame(
        index=release_lon.index[valid],
        geometry=gpd.points_from_xy(release_lon[valid], release_lat[valid]),
        crs=subbasins.crs,
    )
    # HELCOM level-2 subbasins tile the Baltic, so point-in-polygon is the
    # natural join. Keep the first match in case a release point sits on a
    # polygon boundary and hits two subbasins.
    sjoined = gpd.sjoin(release_pts, subbasins, how="left", predicate="within")
    sjoined = sjoined[~sjoined.index.duplicated(keep="first")]
    subbasin_per_traj = sjoined.subbasin.reindex(release_lon.index)
    # Keep the per-trajectory subbasin coord dask-backed and chunk-aligned with
    # ds so downstream ``ds.where(ds.subbasin == sb)`` masks fuse block-wise
    # instead of forcing a client-side fancy index along the concat dim.
    subbasin_dask = da.from_array(
        subbasin_per_traj.values, chunks=ds.chunksizes["trajectory"]
    )
    return ds.assign(subbasin=("trajectory", subbasin_dask))


def relabel_quarter(da, dim="release_quarter"):
    """Replace int 1..4 quarter coord values with JFM/AMJ/JAS/OND labels."""
    return da.assign_coords(
        {dim: [QUARTER_LABELS[int(q)] for q in da[dim].values]}
    )
