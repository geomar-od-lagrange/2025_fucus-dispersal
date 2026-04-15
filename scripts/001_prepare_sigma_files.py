"""Prepare sigma files for Parcels fieldset.

Reads the original sigma files from the shared static data, replaces NaN
values in the first layer with 0, and writes cleaned versions to
<base_path>/output/sigma/.

This is a one-time preprocessing step that must run before 010_FucusDispersal.
"""

import argparse
from pathlib import Path

import xarray as xr


def prepare_sigma(static_path: Path, output_path: Path):
    output_path.mkdir(parents=True, exist_ok=True)

    for resolution in ("fine", "coarse"):
        src = static_path / f"static_file_{resolution}" / f"sigma_file_{resolution}.nc"
        dst = output_path / f"sigma_{resolution}.nc"

        print(f"Reading {src}")
        ds = xr.open_dataset(src)
        ds = ds.assign(sigma=xr.where(ds.layer_number == 1, 0, ds.sigma))
        ds.to_netcdf(dst)
        print(f"Wrote {dst}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-path",
        type=Path,
        default=Path("/gxfs_work/geomar/smomw122/2025_fucus-dispersal"),
        help="Project base path (default: %(default)s)",
    )
    parser.add_argument(
        "--static-path",
        type=Path,
        default=Path("/gxfs_work/geomar/smomw122/bsh_operationalmodel_data"),
        help="Path to shared static data (default: %(default)s)",
    )
    args = parser.parse_args()

    prepare_sigma(
        static_path=args.static_path,
        output_path=args.base_path / "output" / "sigma",
    )
