"""Scan 2D bottom/surface fields for a year and flag cells that are always immobile.

For each (lat, lon) the script asks, across every timestep of every file:
    "Did u or v ever differ from 0.0?"

Output (netCDF):
    - bottom_ever_moved : bool mask, True where bottom u or v was != 0 at any time
    - bottom_max_speed  : float, max sqrt(u^2+v^2) seen at the bottom
    - surface_ever_moved: bool mask, True where surface u or v was != 0 at any time
    - surface_max_speed : float, max sqrt(u^2+v^2) seen at the surface
    - bottom_dead_but_wet : bool mask, True where surface moved at least once
                            but bottom never moved

Also prints a summary: counts of land, always-zero surface, always-zero bottom,
and the "wet-but-dead-bottom" cells we actually care about.

Usage:
    python debug_bottom_stationary.py \\
        --input-dir  /path/to/output/2d_fields \\
        --output-dir /path/to/debug_out \\
        --year 2019 \\
        --res fine
"""

import argparse
from pathlib import Path

import numpy as np
import xarray as xr


def scan_pair(files, label):
    ever_moved = None
    max_speed = None
    n_timesteps = 0
    for i, f in enumerate(files):
        with xr.open_dataset(f) as ds:
            u = ds["uvel"].values
            v = ds["vvel"].values
        moved_here = ((u != 0.0) | (v != 0.0)).any(axis=0)
        speed_here = np.sqrt(u * u + v * v).max(axis=0)
        n_timesteps += u.shape[0]
        if ever_moved is None:
            ever_moved = moved_here
            max_speed = speed_here
        else:
            ever_moved |= moved_here
            np.maximum(max_speed, speed_here, out=max_speed)
        if (i + 1) % 200 == 0:
            print(f"  {label}: {i + 1}/{len(files)} files scanned", flush=True)
    print(f"  {label}: done ({len(files)} files, {n_timesteps} timesteps)", flush=True)
    return ever_moved, max_speed


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", type=Path, required=True,
                   help="Directory holding the *_bottom.nc and *_surface.nc files.")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Where to write the summary netCDF.")
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--res", choices=["fine", "coarse"], required=True)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    pattern_bot = f"c_file_{args.res}_{args.year}*_bottom.nc"
    pattern_surf = f"c_file_{args.res}_{args.year}*_surface.nc"
    files_bot = sorted(args.input_dir.glob(pattern_bot))
    files_surf = sorted(args.input_dir.glob(pattern_surf))

    print(f"Found {len(files_bot)} bottom files, {len(files_surf)} surface files")
    if not files_bot or not files_surf:
        raise SystemExit("No files matched — check --input-dir / --year / --res.")

    # Grab a coordinate template from the first file.
    with xr.open_dataset(files_bot[0]) as ds0:
        lat = ds0["lat"].values
        lon = ds0["lon"].values

    print("Scanning bottom...")
    bot_moved, bot_max = scan_pair(files_bot, "bottom")
    print("Scanning surface...")
    surf_moved, surf_max = scan_pair(files_surf, "surface")

    bot_dead = ~bot_moved
    surf_dead = ~surf_moved
    wet_but_dead_bot = surf_moved & bot_dead

    total = bot_moved.size
    print()
    print("=" * 50)
    print(f"Grid: {bot_moved.shape}  total cells: {total}")
    print(f"Surface never moved (likely land + masked edges): {int(surf_dead.sum())}")
    print(f"Bottom  never moved                             : {int(bot_dead.sum())}")
    print(f"Surface moved ∧ bottom never moved  (DEAD BOT) : {int(wet_but_dead_bot.sum())}")
    print("=" * 50)

    out = xr.Dataset(
        data_vars=dict(
            bottom_ever_moved=(("lat", "lon"), bot_moved),
            bottom_max_speed=(("lat", "lon"), bot_max.astype(np.float32)),
            surface_ever_moved=(("lat", "lon"), surf_moved),
            surface_max_speed=(("lat", "lon"), surf_max.astype(np.float32)),
            bottom_dead_but_wet=(("lat", "lon"), wet_but_dead_bot),
        ),
        coords=dict(lat=lat, lon=lon),
        attrs=dict(
            year=args.year,
            res=args.res,
            n_bottom_files=len(files_bot),
            n_surface_files=len(files_surf),
            source_dir=str(args.input_dir),
            description=(
                "Per-cell aggregates across the year. A cell is 'ever_moved' "
                "if u!=0 or v!=0 at any timestep across all files."
            ),
        ),
    )

    out_path = args.output_dir / f"bottom_stationary_{args.res}_{args.year}.nc"
    out.to_netcdf(out_path)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
