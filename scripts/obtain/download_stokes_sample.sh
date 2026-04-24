#!/usr/bin/env bash
# Download a one-day Stokes drift sample from the CMEMS Baltic Wave Hindcast.
#
# Covers 2020-01-01 to match the BSH minimal data set so smoke runs in
# notebooks 003, 010, and 025 use temporally consistent inputs.
#
# Authentication: copernicusmarine reads credentials from environment variables
#   COPERNICUSMARINE_SERVICE_USERNAME and COPERNICUSMARINE_SERVICE_PASSWORD,
# or from a stored session created by running `pixi run python -m
# copernicusmarine login` once interactively.  Set the env vars (e.g. via
# a .env file or CI secrets) before running this script unattended.
#
# Output: ./data/cmems_stokes_sample/baltic_stokes_20200101.nc
# Idempotent: skips the download if the file is already present.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/data/cmems_stokes_sample"
OUT_FILE="${OUT_DIR}/baltic_stokes_20200101.nc"

if [ -f "${OUT_FILE}" ]; then
    echo "Stokes sample already present — skipping download."
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "Downloading CMEMS Baltic Stokes sample..."
pixi run python -c "
import copernicusmarine
copernicusmarine.subset(
    dataset_id='cmems_mod_bal_wav_my_PT1H-i',
    variables=['VSDX', 'VSDY'],
    start_datetime='2020-01-01T00:00:00',
    end_datetime='2020-01-02T00:00:00',
    output_filename='baltic_stokes_20200101.nc',
    output_directory='${OUT_DIR}',
)
"

echo "Stokes sample written to: ${OUT_FILE}"
