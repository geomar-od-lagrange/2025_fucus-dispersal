#!/usr/bin/env bash
# Download the MADS/SYKE REDLIST_SIS_Macrophytes shapefile from the HELCOM
# MADS catalog — the source of Fucus vesiculosus release-point geometry.
#
# Source: HELCOM metadata record 5848c347-dd45-4135-bbb0-228be9ddeffb
# The MADS portal calls an ArcGIS GP service to produce a zip on demand; this
# script replicates that flow for non-interactive use.
#
# Output: ./data/fucus_redlist_shapefile/<shapefile components>
# Idempotent: skips the download if the .shp is already present.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/data/fucus_redlist_shapefile"
RECORD_UUID="5848c347-dd45-4135-bbb0-228be9ddeffb"
GP_URL="https://maps.helcom.fi/arcgis/rest/services/MADS/tools/GPServer/getAndCompressDataSource/execute"

if [ -f "${OUT_DIR}/REDLIST_SIS_Macrophytes.shp" ]; then
    echo "Fucus REDLIST shapefile already present — skipping download."
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "Requesting REDLIST_SIS_Macrophytes zip from MADS GP service..."
RESPONSE=$(curl -fsSL "${GP_URL}?id=${RECORD_UUID}&f=json")

DOWNLOAD_URL=$(echo "${RESPONSE}" | python3 -c "
import json, sys
d = json.load(sys.stdin)
for r in d['results']:
    if r['paramName'] == 'output_zip':
        print(r['value']['url'])
        break
")

if [ -z "${DOWNLOAD_URL}" ]; then
    echo "ERROR: could not parse download URL from GP service response:" >&2
    echo "${RESPONSE}" >&2
    exit 1
fi

echo "Downloading: ${DOWNLOAD_URL}"
TMPZIP=$(mktemp /tmp/fucus_redlist_XXXXXX.zip)
trap 'rm -f "${TMPZIP}"' EXIT

curl -fsSL -o "${TMPZIP}" "${DOWNLOAD_URL}"
unzip -o "${TMPZIP}" -d "${OUT_DIR}"

echo "REDLIST_SIS_Macrophytes unpacked to: ${OUT_DIR}"
