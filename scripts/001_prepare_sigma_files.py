"""Prepare sigma files for Parcels fieldset.

Reads the sigma files from the BSH store, replaces NaN values in the first
layer with 0, and writes cleaned versions to <output_root>/sigma/.

This is a one-time preprocessing step that must run before 010_FucusDispersal.

Usage:
    python 001_prepare_sigma_files.py \\
        --bsh-root /path/to/bsh_operationalmodel_data \\
        --output-root /path/to/outputs
"""

import argparse
from pathlib import Path

import xarray as xr


def prepare_sigma(bsh_root: Path, output_path: Path):
    output_path.mkdir(parents=True, exist_ok=True)

    for resolution in ("fine", "coarse"):
        src = bsh_root / f"static_file_{resolution}" / f"sigma_file_{resolution}.nc"
        dst = output_path / f"sigma_{resolution}.nc"

        print(f"Reading {src}")
        ds = xr.open_dataset(src)
        ds = ds.assign(sigma=xr.where(ds.layer_number == 1, 0, ds.sigma))
        ds.to_netcdf(dst)
        print(f"Wrote {dst}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--bsh-root",
        type=Path,
        default=Path("../data/bsh_hbmnoku_static"),
        help="BSH store root (contains static_file_fine/, static_file_coarse/, …)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("../output"),
        help="Heavy-outputs root; sigma files are written to <output-root>/sigma/",
    )
    args = parser.parse_args()

    prepare_sigma(
        bsh_root=args.bsh_root,
        output_path=args.output_root / "sigma",
    )
