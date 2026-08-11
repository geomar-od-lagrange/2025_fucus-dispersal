#!/usr/bin/env bash
# Obtain the BSH coastline shore-type classification (flat_fraction table).
#
# Source: the sidecar repository
#   https://github.com/geomar-od-lagrange/2025_fucus-dispersal_shoreclass
# which classifies the BSH HBMnoku tracer-cell outline from HELCOM BRISK and
# Copernicus CLMS Coastal Zones evidence.  See docs/beaching.md for how this
# study consumes it.
#
# Unlike the other obtain scripts there is no upstream URL to curl: the sidecar
# gitignores its data/processed/, so the table exists only as the output of its
# own pipeline.  This script therefore *derives* the blob from a sidecar
# checkout rather than downloading it.
#
# The full table is 6.8 MB and carries 23 columns, most of them the per-source
# evidence the sidecar needs to revisit its own precedence rule.  This study
# reads five, so the blob committed to the data twin is slimmed — the same
# treatment as the derived release-points geojson baked by 000.
#
# Usage:
#   scripts/obtain/obtain_shoreclass.sh [path/to/2025_fucus-dispersal_shoreclass]
#
# Defaults to a sibling checkout next to this repository.
# Output: ./data/shoreclass_bsh_coastline/bsh_coastline_k2_flatfraction.parquet
# Idempotent: skips if the parquet is already present.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/data/shoreclass_bsh_coastline"
OUT_FILE="${OUT_DIR}/bsh_coastline_k2_flatfraction.parquet"

SIDECAR="${1:-${REPO_ROOT}/../2025_fucus-dispersal_shoreclass}"
SRC="${SIDECAR}/data/processed/bsh_coastline_k2_flatfraction.parquet"

if [ -f "${OUT_FILE}" ]; then
    echo "Shoreclass table already present — skipping."
    exit 0
fi

if [ ! -f "${SRC}" ]; then
    cat >&2 <<EOF
ERROR: no classified coastline at
  ${SRC}

The sidecar does not publish this file; build it there first:

  git clone https://github.com/geomar-od-lagrange/2025_fucus-dispersal_shoreclass
  cd 2025_fucus-dispersal_shoreclass
  pixi install
  pixi run pipeline        # stages 0-5; 0 and 1 need network + a CLMS key

then re-run this script, passing the checkout path if it is not a sibling of
this repository.
EOF
    exit 1
fi

mkdir -p "${OUT_DIR}"

echo "Slimming ${SRC}"
pixi run --manifest-path "${REPO_ROOT}/pixi.toml" python - "${SRC}" "${OUT_FILE}" <<'PY'
import sys

import pandas as pd

src, dst = sys.argv[1], sys.argv[2]
# Midpoint in 3035 (the CRS 024d's beaching raster is built in), the value, the
# length weight, and the provenance/domain flags kept for diagnostics.  The
# geometry column is dropped: it duplicates x_3035/y_3035 as a POINT, and 024d
# snaps coordinates to a raster rather than doing geometry work.
keep = ["x_3035", "y_3035", "flat_fraction", "seg_len_m", "source", "in_baltic"]
df = pd.read_parquet(src, columns=keep)
df.to_parquet(dst, index=False)

baltic = df[df["in_baltic"]]
attributed = baltic["flat_fraction"].notna()
km = baltic["seg_len_m"] / 1e3
print(f"{len(df):,} sub-segments, {df['seg_len_m'].sum() / 1e3:,.0f} km total")
print(f"Baltic: {km.sum():,.0f} km, "
      f"{100 * km[attributed].sum() / km.sum():.1f}% attributed, "
      f"mean flat_fraction {(km[attributed] * baltic['flat_fraction'][attributed]).sum() / km[attributed].sum():.3f}")
PY

echo "Shoreclass table written to: ${OUT_FILE}"
