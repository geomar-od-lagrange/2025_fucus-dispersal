"""Extract 2D surface, bottom, and surface+stokes fields from a single BSH c_file.

Usage:
    python 003_prepare_2d_fields.py \\
        --c-file /path/to/c_file_fine_2020010100_000_006.nc \\
        --output-root /path/to/outputs

Stokes files (both Baltic high-res and WAVERYS) are read from
<output-root>/stokes/<product>/<year>/ (written by 002_download_stokes.py).
The surface_stokes variant layers Baltic high-res over WAVERYS so the
German Bight strip ~6.2–9 °E (where the Baltic product is undefined)
gets WAVERYS values rather than zero.

2D output files are written to <output-root>/2d_fields/.

Designed to be parallelized: one invocation per c_file.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import xarray as xr


def load_surface(ds):
    """Extract surface layer (layer_number=1) as 2D."""
    u = ds["uvel"].sel(layer_number=1).drop_vars("layer_number")
    v = ds["vvel"].sel(layer_number=1).drop_vars("layer_number")
    return u, v


def load_bottom(ds):
    """Extract deepest fluid layer as 2D, independently for U and V.

    BSH stores `0.0` (not NaN) at C-grid U/V faces blocked by the shallower
    neighbour's bathymetry — a no-flux wall, not fluid. Picking the deepest
    non-NaN would land on that wall at every cell whose east (U) or south
    (V) neighbour is shallower. Instead, for each face, take the deepest
    layer that is both non-NaN AND non-zero. U and V are resolved
    independently because their face bathymetries — min(bathy_self,
    bathy_east) for U, min(bathy_self, bathy_south) for V — can differ.

    Indices are computed from t=0; assumes static bathymetry (no wet/dry).
    """
    n_layers = len(ds.layer_number)

    def deepest_fluid_idx(var):
        a = ds[var].isel(time=0)
        ok = a.notnull() & (a != 0.0)
        ok_rev = ok.isel(layer_number=slice(None, None, -1))
        return (n_layers - 1 - ok_rev.argmax(dim="layer_number")).drop_vars("time")

    u_idx = deepest_fluid_idx("uvel")
    v_idx = deepest_fluid_idx("vvel")

    u = ds["uvel"].isel(layer_number=u_idx).drop_vars("layer_number")
    v = ds["vvel"].isel(layer_number=v_idx).drop_vars("layer_number")
    return u, v


def land_mask_from_surface(ds):
    """Land mask from surface layer NaN pattern. True = land."""
    return ds["uvel"].sel(layer_number=1).isel(time=0).isnull()


def roll_v_to_nemo(v):
    """Roll V from BSH south-face to NEMO north-face convention.

    BSH: V[j,i] on south face of cell (j,i).
    NEMO: V[j,i] on north face of cell (j,i).
    South face of cell j = north face of cell j+1 (descending lat: j+1 = south).
    So V_nemo[j] = V_bsh[j-1], i.e. roll by +1 along lat.

    The wrapped boundary row (j=0, northernmost) gets V_bsh[N-1], which
    was zeroed by the land mask (domain boundary). So V_nemo[0] = 0.
    """
    return v.roll(lat=1, roll_coords=False)


def apply_land_mask_at_edges(u, v, land_mask):
    """Zero velocities at edges adjacent to land.

    U (eastern edge): zero if cell [j,i] or cell [j,i+1] is land.
    V (southern edge): zero if cell [j,i] or cell [j+1,i] is land.
    """
    mask = land_mask.values

    u_mask = mask[:, :-1] | mask[:, 1:]
    u_mask = np.concatenate([u_mask, np.ones((mask.shape[0], 1), dtype=bool)], axis=1)

    v_mask = mask[:-1, :] | mask[1:, :]
    v_mask = np.concatenate([v_mask, np.ones((1, mask.shape[1]), dtype=bool)], axis=0)

    u_out = u.where(~xr.DataArray(u_mask, dims=["lat", "lon"]), 0.0)
    v_out = v.where(~xr.DataArray(v_mask, dims=["lat", "lon"]), 0.0)
    return u_out, v_out


# Iterations of 3x3 rolling-mean spread applied to wave-model fields
# before interpolation onto BSH faces. Determined empirically: N=5
# covers 100% of BSH potentially-wet cells (H0 not NaN) inside each
# wave-model file's bounding box, against both Baltic high-res (2 km)
# and WAVERYS (0.2°). The two grids are fixed, so this never changes;
# raise it only if the wave-model files swap to coarser products.
SPREAD_N_ITER = 5


def _spread_into_nan(arr, n_iter):
    """Iteratively fill NaN cells with the mean of valid 3x3 neighbours.

    Each iteration extends valid values outward by one cell, leaving
    original valid cells untouched.
    """
    for _ in range(n_iter):
        rolled = arr.rolling(lat=3, lon=3, center=True, min_periods=1).mean()
        arr = arr.where(arr.notnull(), rolled)
    return arr


def _interp_stokes_to_bsh(stokes_ds, bsh_lon, bsh_lat, bsh_time, fill_value):
    """Interpolate one Stokes A-grid onto BSH C-grid U/V points and time.

    BSH uses NEMO-like C-grid: lon/lat are cell centres on a regular grid.
    U-points sit on the eastern cell edge:  (lon + 0.5*dlon, lat)
    V-points sit on the southern cell edge: (lon, lat - 0.5*dlat)

    Before the interp, wave-model values spread outward into NaN/land
    cells via ``SPREAD_N_ITER`` iterations of a 3x3 rolling mean — the
    wave-model coastline rarely matches BSH's, and without spreading
    the BSH-side interp at the coastline would pull toward zero-filled
    land and underestimate Stokes there.
    """
    dlon = float(bsh_lon[1] - bsh_lon[0])
    dlat = float(bsh_lat[1] - bsh_lat[0])

    lon_u = bsh_lon + 0.5 * dlon
    lat_v = bsh_lat - 0.5 * abs(dlat)  # abs: works whether lat is ascending or descending

    stokes_flipped = (
        stokes_ds
        .sortby("latitude", ascending=False)
        .rename({"latitude": "lat", "longitude": "lon"})
    )

    vsdx = _spread_into_nan(stokes_flipped["VSDX"], n_iter=SPREAD_N_ITER)
    vsdy = _spread_into_nan(stokes_flipped["VSDY"], n_iter=SPREAD_N_ITER)

    u_stokes = vsdx.interp(
        lon=lon_u, lat=bsh_lat, method="linear",
        kwargs={"fill_value": fill_value},
    ).interp(time=bsh_time, method="nearest")
    # Reassign lon coordinate to match BSH implicit staggering convention
    u_stokes = u_stokes.assign_coords(lon=bsh_lon)

    v_stokes = vsdy.interp(
        lon=bsh_lon, lat=lat_v, method="linear",
        kwargs={"fill_value": fill_value},
    ).interp(time=bsh_time, method="nearest")
    # Reassign lat coordinate to match BSH implicit staggering convention
    v_stokes = v_stokes.assign_coords(lat=bsh_lat)

    return u_stokes, v_stokes


def interpolate_stokes(stokes_baltic_ds, stokes_waverys_ds, bsh_lon, bsh_lat, bsh_time):
    """Layered Stokes onto BSH C-grid: Baltic high-res over WAVERYS over zero.

    Baltic high-res covers ~9–30 °E. Outside its footprint (notably the
    German Bight strip ~6.2–9 °E that BSH fine grid spans), interpolation
    fills NaN; we then substitute the WAVERYS global value at those
    cells. WAVERYS is interpolated with fill_value=0.0 so the rare BSH
    cells outside its (global) footprint default to no Stokes drift.
    """
    u_baltic, v_baltic = _interp_stokes_to_bsh(
        stokes_baltic_ds, bsh_lon, bsh_lat, bsh_time, fill_value=np.nan,
    )
    u_waverys, v_waverys = _interp_stokes_to_bsh(
        stokes_waverys_ds, bsh_lon, bsh_lat, bsh_time, fill_value=0.0,
    )
    u_stokes = u_baltic.where(u_baltic.notnull(), u_waverys)
    v_stokes = v_baltic.where(v_baltic.notnull(), v_waverys)
    return u_stokes, v_stokes


def write_2d(u, v, path, lon_f=None, lat_f=None):
    """Write a 2D UV dataset to netCDF.

    If lon_f/lat_f are provided, replace the coordinate arrays with
    NE F-point values (for Parcels NEMO C-grid convention).
    """
    ds = xr.Dataset({"uvel": u, "vvel": v})
    if lon_f is not None:
        ds = ds.assign_coords(lon=lon_f, lat=lat_f)
    ds.to_netcdf(path)


def derive_stokes_path(c_file: Path, stokes_dir: Path, product: str) -> Path:
    """Derive Stokes file path for a given product from c_file name.

    c_file_fine_2020010100_000_006.nc, product='baltic_highres'
    -> stokes_dir/baltic_highres/2020/stokes_20200101.nc
    """
    match = re.search(r"(\d{4})(\d{4})\d{2}_", c_file.name)
    if not match:
        raise ValueError(f"Cannot parse date from c_file name: {c_file.name}")
    year = match.group(1)
    date = match.group(1) + match.group(2)
    return stokes_dir / product / year / f"stokes_{date}.nc"


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--c-file", type=Path, required=True,
                        help="BSH c_file (3D, with layer_number)")
    parser.add_argument(
        "--output-root", type=Path, default=Path("../output"),
        help="Heavy-outputs root; stokes files are read from <output-root>/stokes/ "
             "and 2D outputs are written to <output-root>/2d_fields/ "
             "(default: %(default)s)",
    )
    args = parser.parse_args()

    stokes_dir = args.output_root / "stokes"
    output_dir = args.output_root / "2d_fields"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.c_file.stem

    out_surf = output_dir / f"{stem}_surface.nc"
    out_bot = output_dir / f"{stem}_bottom.nc"
    out_stokes = output_dir / f"{stem}_surface_stokes.nc"

    # Skip if all outputs exist
    if out_surf.exists() and out_bot.exists() and out_stokes.exists():
        print(f"Skipping {args.c_file.name} (all outputs exist)")
        return

    print(f"Processing {args.c_file.name}")
    ds = xr.open_dataset(args.c_file)
    land_mask = land_mask_from_surface(ds)

    # NE F-point coordinates for Parcels NEMO C-grid.
    # Cell (yi, xi) in Parcels = T-cell (yi+1, xi+1).
    dlon = float(ds.lon[1] - ds.lon[0])
    dlat = float(ds.lat[1] - ds.lat[0])  # negative (descending)
    lon_f = ds.lon.values + dlon / 2
    lat_f = ds.lat.values - dlat / 2      # -dlat/2 = +|dlat|/2, shifts north

    # Fine grids: crop dead boundary rows/columns so particles fall
    # through to the coarse grid via NestedField instead of stalling.
    # Dead boundaries come from:
    #   - lat row 0 (N): V roll wraps boundary zero here
    #   - lat row -1 (S): land mask appends all-True row
    #   - lon col -1 (E): land mask appends all-True column
    is_fine = "fine" in args.c_file.name
    lat_crop = slice(1, -1) if is_fine else slice(None)
    lon_crop = slice(None, -1) if is_fine else slice(None)

    # Surface — mask in BSH convention, then roll V to NEMO north-face
    u_surf, v_surf = load_surface(ds)
    u_surf, v_surf = apply_land_mask_at_edges(u_surf, v_surf, land_mask)
    if not out_surf.exists():
        u_out = u_surf.isel(lat=lat_crop, lon=lon_crop)
        v_out = roll_v_to_nemo(v_surf).isel(lat=lat_crop, lon=lon_crop)
        write_2d(u_out, v_out, out_surf, lon_f[lon_crop], lat_f[lat_crop])
        print(f"  surface: {out_surf}  U=[{float(u_out.min()):.3f}, {float(u_out.max()):.3f}]")

    # Bottom
    if not out_bot.exists():
        u_bot, v_bot = load_bottom(ds)
        u_bot, v_bot = apply_land_mask_at_edges(u_bot, v_bot, land_mask)
        u_out = u_bot.isel(lat=lat_crop, lon=lon_crop)
        v_out = roll_v_to_nemo(v_bot).isel(lat=lat_crop, lon=lon_crop)
        write_2d(u_out, v_out, out_bot, lon_f[lon_crop], lat_f[lat_crop])
        print(f"  bottom:  {out_bot}  U=[{float(u_out.min()):.3f}, {float(u_out.max()):.3f}]")

    # Surface + Stokes — layered Stokes with per-timestep face mask.
    # The spread-and-interp produces Stokes everywhere a wave-model
    # value can be reached; BSH's per-timestep face state then decides
    # whether Stokes is allowed at that face *now*. BSH writes 0.0 at
    # no-slip walls and at faces shut off because a tidal flat is dry —
    # we explicitly zero Stokes there so particles don't advect through
    # blocked faces (this was a likely source of "effective beaching"
    # near shore in earlier runs). Outside-model cells where u_surf is
    # NaN propagate NaN through the sum; Parcels treats NaN as 0.
    if not out_stokes.exists():
        baltic_file = derive_stokes_path(args.c_file, stokes_dir, "baltic_highres")
        waverys_file = derive_stokes_path(args.c_file, stokes_dir, "waverys")
        missing = [f.name for f in (baltic_file, waverys_file) if not f.exists()]
        if missing:
            print(f"  stokes:  SKIPPED (not yet downloaded: {missing})")
            return
        stokes_baltic = xr.open_dataset(baltic_file)
        stokes_waverys = xr.open_dataset(waverys_file)
        u_stokes, v_stokes = interpolate_stokes(
            stokes_baltic, stokes_waverys, ds.lon, ds.lat, ds.time,
        )
        u_stokes = u_stokes.where(u_surf.notnull() & (u_surf != 0), 0)
        v_stokes = v_stokes.where(v_surf.notnull() & (v_surf != 0), 0)
        u_total = u_surf + u_stokes
        v_total = v_surf + v_stokes
        u_total, v_total = apply_land_mask_at_edges(u_total, v_total, land_mask)
        u_out = u_total.isel(lat=lat_crop, lon=lon_crop)
        v_out = roll_v_to_nemo(v_total).isel(lat=lat_crop, lon=lon_crop)
        write_2d(u_out, v_out, out_stokes, lon_f[lon_crop], lat_f[lat_crop])
        print(f"  stokes:  {out_stokes}  U=[{float(u_out.min()):.3f}, {float(u_out.max()):.3f}]")


if __name__ == "__main__":
    main()
