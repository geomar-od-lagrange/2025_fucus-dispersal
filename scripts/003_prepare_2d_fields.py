"""Extract 2D surface, bottom, and surface+stokes fields from a single BSH c_file.

Usage:
    python 003_prepare_2d_fields.py \
        --c-file /path/to/c_file_fine_2020010100_000_006.nc \
        --stokes-dir /path/to/stokes/ \
        --output-dir /path/to/output_2d/

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
    """Extract deepest valid layer as 2D.

    The deepest valid layer varies spatially. We find the index from the
    first timestep and apply to all timesteps. This assumes fixed bathymetry
    (no wetting/drying).
    """
    n_layers = len(ds.layer_number)
    valid_reversed = ds["uvel"].isel(time=0).notnull().isel(
        layer_number=slice(None, None, -1)
    )
    deepest_idx = (
        n_layers - 1 - valid_reversed.argmax(dim="layer_number")
    ).drop_vars("time")

    u = ds["uvel"].isel(layer_number=deepest_idx).drop_vars("layer_number")
    v = ds["vvel"].isel(layer_number=deepest_idx).drop_vars("layer_number")
    return u, v


def land_mask_from_surface(ds):
    """Land mask from surface layer NaN pattern. True = land."""
    return ds["uvel"].sel(layer_number=1).isel(time=0).isnull()


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


def interpolate_stokes(stokes_ds, bsh_lon, bsh_lat, bsh_time):
    """Interpolate Stokes A-grid onto BSH c_file coordinates and time."""
    stokes_flipped = (
        stokes_ds
        .sortby("latitude", ascending=False)
        .rename({"latitude": "lat", "longitude": "lon"})
    )

    u_stokes = stokes_flipped["VSDX"].interp(
        lon=bsh_lon, lat=bsh_lat, method="linear",
        kwargs={"fill_value": 0.0},
    ).interp(time=bsh_time, method="nearest")

    v_stokes = stokes_flipped["VSDY"].interp(
        lon=bsh_lon, lat=bsh_lat, method="linear",
        kwargs={"fill_value": 0.0},
    ).interp(time=bsh_time, method="nearest")

    return u_stokes, v_stokes


def write_2d(u, v, path):
    """Write a 2D UV dataset to netCDF."""
    xr.Dataset({"uvel": u, "vvel": v}).to_netcdf(path)


def derive_stokes_path(c_file: Path, stokes_dir: Path) -> Path:
    """Derive stokes file path from c_file name.

    c_file_fine_2020010100_000_006.nc -> stokes_dir/2020/stokes_20200101.nc
    """
    match = re.search(r"(\d{4})(\d{4})\d{2}_", c_file.name)
    if not match:
        raise ValueError(f"Cannot parse date from c_file name: {c_file.name}")
    year = match.group(1)
    date = match.group(1) + match.group(2)
    return stokes_dir / year / f"stokes_{date}.nc"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--c-file", type=Path, required=True,
                        help="BSH c_file (3D, with layer_number)")
    parser.add_argument("--stokes-dir", type=Path, required=True,
                        help="Stokes directory (file auto-derived from c_file date)")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Output directory for 2D files")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.c_file.stem

    out_surf = args.output_dir / f"{stem}_surface.nc"
    out_bot = args.output_dir / f"{stem}_bottom.nc"
    out_stokes = args.output_dir / f"{stem}_surface_stokes.nc"

    # Skip if all outputs exist
    if out_surf.exists() and out_bot.exists() and out_stokes.exists():
        print(f"Skipping {args.c_file.name} (all outputs exist)")
        return

    print(f"Processing {args.c_file.name}")
    ds = xr.open_dataset(args.c_file)
    land_mask = land_mask_from_surface(ds)

    # Surface
    if not out_surf.exists():
        u_surf, v_surf = load_surface(ds)
        u_surf, v_surf = apply_land_mask_at_edges(u_surf, v_surf, land_mask)
        write_2d(u_surf, v_surf, out_surf)
        print(f"  surface: {out_surf}  U=[{float(u_surf.min()):.3f}, {float(u_surf.max()):.3f}]")
    else:
        # Still need surface for stokes sum
        u_surf, v_surf = load_surface(ds)
        u_surf, v_surf = apply_land_mask_at_edges(u_surf, v_surf, land_mask)

    # Bottom
    if not out_bot.exists():
        u_bot, v_bot = load_bottom(ds)
        u_bot, v_bot = apply_land_mask_at_edges(u_bot, v_bot, land_mask)
        write_2d(u_bot, v_bot, out_bot)
        print(f"  bottom:  {out_bot}  U=[{float(u_bot.min()):.3f}, {float(u_bot.max()):.3f}]")

    # Surface + Stokes
    if not out_stokes.exists():
        stokes_file = derive_stokes_path(args.c_file, args.stokes_dir)
        if not stokes_file.exists():
            print(f"  stokes:  SKIPPED (not yet downloaded: {stokes_file.name})")
            return
        stokes_ds = xr.open_dataset(stokes_file).fillna(0.0)
        u_stokes, v_stokes = interpolate_stokes(
            stokes_ds, ds.lon, ds.lat, ds.time,
        )
        u_total = u_surf + u_stokes
        v_total = v_surf + v_stokes
        u_total, v_total = apply_land_mask_at_edges(u_total, v_total, land_mask)
        write_2d(u_total, v_total, out_stokes)
        print(f"  stokes:  {out_stokes}  U=[{float(u_total.min()):.3f}, {float(u_total.max()):.3f}]")


if __name__ == "__main__":
    main()
