"""Download daily Stokes drift files from CMEMS Baltic Wave Hindcast.

Downloads VSDX/VSDY for 2016-2025, one file per day.
Designed to be resumable — skips files that already exist.

Usage:
    python 002_download_stokes.py --output-dir /path/to/stokes/
    python 002_download_stokes.py --output-dir /path/to/stokes/ --year 2020
    python 002_download_stokes.py --output-dir /path/to/stokes/ --year 2020 --month 1
"""

import argparse
from datetime import date, timedelta
from pathlib import Path

import copernicusmarine


DATASET_ID = "cmems_mod_bal_wav_my_PT1H-i"
VARIABLES = ["VSDX", "VSDY"]


def download_day(day: date, output_dir: Path):
    """Download one day of Stokes drift data. Skips if output exists."""
    filename = f"stokes_{day:%Y%m%d}.nc"
    filepath = output_dir / f"{day.year}" / filename

    if filepath.exists():
        return False  # already downloaded

    filepath.parent.mkdir(parents=True, exist_ok=True)
    next_day = day + timedelta(days=1)
    copernicusmarine.subset(
        dataset_id=DATASET_ID,
        variables=VARIABLES,
        start_datetime=f"{day}T00:00:00",
        end_datetime=f"{next_day}T00:00:00",
        output_filename=filename,
        output_directory=str(filepath.parent),
    )
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Base output directory for Stokes files")
    parser.add_argument("--start-year", type=int, default=2016)
    parser.add_argument("--end-year", type=int, default=2025)
    parser.add_argument("--year", type=int, default=None,
                        help="Download only this year")
    parser.add_argument("--month", type=int, default=None,
                        help="Download only this month (requires --year)")
    args = parser.parse_args()

    if args.month is not None and args.year is None:
        parser.error("--month requires --year")

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

    day = start
    n_downloaded = 0
    n_skipped = 0
    n_failed = 0
    while day < end:
        try:
            downloaded = download_day(day, args.output_dir)
            if downloaded:
                n_downloaded += 1
                print(f"Downloaded {day}")
            else:
                n_skipped += 1
        except Exception as e:
            n_failed += 1
            print(f"FAILED {day}: {e}")
        day += timedelta(days=1)

    print(f"\nDone. Downloaded: {n_downloaded}, skipped: {n_skipped}, failed: {n_failed}")


if __name__ == "__main__":
    main()
