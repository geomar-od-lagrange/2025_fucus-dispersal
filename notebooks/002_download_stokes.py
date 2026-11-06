"""Download daily Stokes drift files from CMEMS for layered surface_stokes forcing.

Two products are pulled side by side:

- ``baltic_highres``: 2 km hourly, Baltic + Danish Straits + Kattegat
  (~9–30 °E). Primary source. Served by the multi-year hindcast
  ``cmems_mod_bal_wav_my_PT1H-i`` up to the hindcast's end and by the
  analysis/forecast ``cmems_mod_bal_wav_anfc_PT1H-i`` beyond it. The two
  share grid, variables, units and hourly step, so a day from either is
  interchangeable for 003 and 024d.
- ``waverys``: 0.2° 3-hourly global (``cmems_mod_glo_wav_my_0.2deg_PT3H-i``).
  Fallback for the German Bight strip ~6.2–9 °E where the Baltic high-res
  product is undefined. It has no NRT continuation here: the global
  analysis/forecast wave product is on a 0.083° grid, so it would change
  the source grid rather than extend this one.

Files are written under ``<output-root>/stokes/<product>/<year>/stokes_<YYYYMMDD>.nc``.
Designed to be resumable — skips files that already exist.

Usage:
    python 002_download_stokes.py --output-root /path/to/outputs
    python 002_download_stokes.py --output-root /path/to/outputs \\
        --start-date 2026-01-01 --end-date 2026-09-14
"""

import argparse
from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import copernicusmarine


@dataclass(frozen=True)
class Product:
    """The dataset(s) covering one product directory, newest source last."""

    multiyear: str
    nrt: str | None


PRODUCTS = {
    "baltic_highres": Product(
        multiyear="cmems_mod_bal_wav_my_PT1H-i",
        nrt="cmems_mod_bal_wav_anfc_PT1H-i",
    ),
    "waverys": Product(
        multiyear="cmems_mod_glo_wav_my_0.2deg_PT3H-i",
        nrt=None,
    ),
}
VARIABLES = ["VSDX", "VSDY"]


def last_full_day(dataset_id: str) -> date:
    """Last day D the dataset covers over the whole ``[D 00:00, D+1 00:00]``.

    A day file spans 25 hourly (or 9 three-hourly) steps ending at the next
    midnight, so the last usable day is one short of the dataset's final
    timestamp whether that timestamp is a midnight or mid-day.
    """
    ds = copernicusmarine.open_dataset(dataset_id=dataset_id, variables=VARIABLES)
    end = ds.time.max().values.astype("datetime64[s]").item()
    return end.date() - timedelta(days=1)


def download_day(day: date, dataset_id: str, out_dir: Path):
    """Download one day of Stokes drift. Skips if output exists."""
    filename = f"stokes_{day:%Y%m%d}.nc"
    filepath = out_dir / filename

    if filepath.exists():
        return False

    out_dir.mkdir(parents=True, exist_ok=True)
    next_day = day + timedelta(days=1)
    copernicusmarine.subset(
        dataset_id=dataset_id,
        variables=VARIABLES,
        start_datetime=f"{day}T00:00:00",
        end_datetime=f"{next_day}T00:00:00",
        output_filename=filename,
        output_directory=str(out_dir),
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
    parser.add_argument(
        "--start-date", type=date.fromisoformat, default=date(2016, 1, 1),
        help="First day to download, YYYY-MM-DD (default: %(default)s)",
    )
    parser.add_argument(
        "--end-date", type=date.fromisoformat, default=None,
        help="Last day to download, inclusive. Default: the last day each "
             "product covers, capped at yesterday so the NRT forecast tail "
             "is not ingested as if it were analysis.",
    )
    args = parser.parse_args()

    stokes_dir = args.output_root / "stokes"
    yesterday = date.today() - timedelta(days=1)

    for name, product in PRODUCTS.items():
        multiyear_end = last_full_day(product.multiyear)
        product_end = (
            last_full_day(product.nrt) if product.nrt else multiyear_end
        )
        end = args.end_date or min(product_end, yesterday)
        nrt_note = (
            f"{product.nrt} through {product_end}"
            if product.nrt else "no NRT continuation"
        )
        print(f"{name}: {product.multiyear} through {multiyear_end}, {nrt_note}")

        n_downloaded = n_skipped = n_failed = 0
        day = args.start_date
        while day <= end:
            dataset_id = product.multiyear if day <= multiyear_end else product.nrt
            if dataset_id is None:
                print(f"  {name}: no NRT continuation past {multiyear_end}, stopping")
                break
            out_dir = stokes_dir / name / f"{day.year}"
            try:
                if download_day(day, dataset_id, out_dir):
                    n_downloaded += 1
                    print(f"Downloaded {name} {day} from {dataset_id}")
                else:
                    n_skipped += 1
            except Exception as e:
                n_failed += 1
                print(f"FAILED {name} {day} from {dataset_id}: {e}")
            day += timedelta(days=1)

        print(
            f"  {name}: downloaded={n_downloaded}, "
            f"skipped={n_skipped}, failed={n_failed}"
        )


if __name__ == "__main__":
    main()
