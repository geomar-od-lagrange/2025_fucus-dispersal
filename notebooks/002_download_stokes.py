"""Download daily Stokes drift files from CMEMS for layered surface_stokes forcing.

Two products are pulled side by side:

- ``baltic_highres`` (``cmems_mod_bal_wav_my_PT1H-i``): 2 km hourly,
  Baltic + Danish Straits + Kattegat (~9–30 °E). Primary source.
- ``waverys`` (``cmems_mod_glo_wav_my_0.2deg_PT3H-i``): 0.2° 3-hourly
  global. Fallback for the German Bight strip ~6.2–9 °E where the
  Baltic high-res product is undefined.

Files are written under ``<output-root>/stokes/<product>/<year>/stokes_<YYYYMMDD>.nc``.
Designed to be resumable — skips files that already exist.

Usage:
    python 002_download_stokes.py --output-root /path/to/outputs
    python 002_download_stokes.py --output-root /path/to/outputs --year 2020
    python 002_download_stokes.py --output-root /path/to/outputs --year 2020 --month 1
"""

import argparse
from datetime import date, timedelta
from pathlib import Path

import copernicusmarine


PRODUCTS = {
    "baltic_highres": "cmems_mod_bal_wav_my_PT1H-i",
    "waverys": "cmems_mod_glo_wav_my_0.2deg_PT3H-i",
}
VARIABLES = ["VSDX", "VSDY"]


def download_day(day: date, product: str, dataset_id: str, stokes_dir: Path):
    """Download one day of Stokes drift for one product. Skips if output exists."""
    filename = f"stokes_{day:%Y%m%d}.nc"
    filepath = stokes_dir / product / f"{day.year}" / filename

    if filepath.exists():
        return False

    filepath.parent.mkdir(parents=True, exist_ok=True)
    next_day = day + timedelta(days=1)
    copernicusmarine.subset(
        dataset_id=dataset_id,
        variables=VARIABLES,
        start_datetime=f"{day}T00:00:00",
        end_datetime=f"{next_day}T00:00:00",
        output_filename=filename,
        output_directory=str(filepath.parent),
    )
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output-root", type=Path, default=Path("../output"),
        help="Heavy-outputs root; Stokes files are written under "
             "<output-root>/stokes/<product>/ (default: %(default)s)",
    )
    parser.add_argument("--start-year", type=int, default=2016)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--year", type=int, default=None,
                        help="Download only this year")
    parser.add_argument("--month", type=int, default=None,
                        help="Download only this month (requires --year)")
    args = parser.parse_args()

    if args.month is not None and args.year is None:
        parser.error("--month requires --year")

    stokes_dir = args.output_root / "stokes"

    if args.year is not None:
        start_year = end_year = args.year
    else:
        start_year, end_year = args.start_year, args.end_year

    start = date(start_year, args.month or 1, 1)
    if args.month:
        if args.month == 12:
            end = date(start_year + 1, 1, 1)
        else:
            end = date(start_year, args.month + 1, 1)
    else:
        end = date(end_year + 1, 1, 1)

    n_downloaded = {p: 0 for p in PRODUCTS}
    n_skipped = {p: 0 for p in PRODUCTS}
    n_failed = {p: 0 for p in PRODUCTS}
    day = start
    while day < end:
        for product, dataset_id in PRODUCTS.items():
            try:
                downloaded = download_day(day, product, dataset_id, stokes_dir)
                if downloaded:
                    n_downloaded[product] += 1
                    print(f"Downloaded {product} {day}")
                else:
                    n_skipped[product] += 1
            except Exception as e:
                n_failed[product] += 1
                print(f"FAILED {product} {day}: {e}")
        day += timedelta(days=1)

    print()
    print("Summary:")
    for product in PRODUCTS:
        print(
            f"  {product}: downloaded={n_downloaded[product]}, "
            f"skipped={n_skipped[product]}, failed={n_failed[product]}"
        )


if __name__ == "__main__":
    main()
