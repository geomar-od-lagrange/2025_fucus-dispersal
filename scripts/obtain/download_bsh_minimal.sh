#!/usr/bin/env bash
# Download a minimal subset of BSH operational model output from BSH OpenData
# into ./data/bsh_minimal/.
#
# Public source, no authentication: https://gdi.bsh.de/data/OpenData/OperationalModel/
# Files are released under dl-de/by-2-0 — attribution per ATTRIBUTION.md.
#
# Subset: a full day (2020-01-01, 00 / 06 / 12 / 18 UTC) for four file types
# × two grid resolutions = 32 .nc files, ≈ 985 MB total.
#
# File types:
#   c_file — 3D currents (U, V) — primary Parcels input after 003 preprocessing
#   h_file — 3D layer thickness — needed to derive sigma levels at free surface
#   t_file — 3D temperature and salinity — shipped for habitat-suitability checks
#   z_file — 2D sea-surface elevation — needed for dry/wet analysis + sigma derivation
#
# Static files (lonlat, H0, sigma) are NOT downloaded here — they are derived
# artifacts from the per-timestep output and already live under
# data/bsh_minimal/static_file_{fine,coarse}/. If a future rebuild needs to
# regenerate them from the downloaded files, add a separate derivation step.
#
# Idempotent: an individual file is skipped if already present at the destination.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_BASE="${REPO_ROOT}/data/bsh_minimal"
BSH_OPENDATA="https://gdi.bsh.de/data/OpenData/OperationalModel"

DATE="20200101"
YEAR="${DATE:0:4}"
HOURS=(00 06 12 18)

for res in fine coarse; do
    for prefix in c h t z; do
        subdir="${OUT_BASE}/${prefix}_file_${res}_${YEAR}"
        mkdir -p "${subdir}"
        for hh in "${HOURS[@]}"; do
            fn="${prefix}_file_${res}_${DATE}${hh}_000_006.nc"
            dest="${subdir}/${fn}"
            url="${BSH_OPENDATA}/${prefix}_file_${res}_${YEAR}/${fn}"
            if [ -f "${dest}" ]; then
                echo "  skip ${fn}"
            else
                echo "  downloading ${fn} ..."
                curl -fsSL -o "${dest}" "${url}"
            fi
        done
    done
done

echo "BSH minimal subset downloaded to: ${OUT_BASE}"
