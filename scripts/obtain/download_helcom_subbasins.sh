#!/usr/bin/env bash
# Download HELCOM subbasins 2022 level-2 shapefile from the HELCOM MADS catalog.
#
# Source: HELCOM metadata record d4b6296c-fd19-462c-94d2-4c81b9313d77
# The MADS portal doesn't expose a static download URL; instead it calls an
# ArcGIS GP service to package the shapefile on demand, then returns a
# short-lived URL to the resulting zip.  This script replicates that flow so
# it works non-interactively without accepting the browser terms dialog.
#
# Output: ./data/helcom_subbasins/<shapefile components>
# Idempotent: skips the download if the .shp is already present.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/data/helcom_subbasins"
RECORD_UUID="d4b6296c-fd19-462c-94d2-4c81b9313d77"
GP_URL="https://maps.helcom.fi/arcgis/rest/services/MADS/tools/GPServer/getAndCompressDataSource/execute"

if [ -f "${OUT_DIR}/HELCOM_subbasins_2022_level2.shp" ]; then
    echo "HELCOM subbasins shapefile already present — skipping download."
    exit 0
fi

mkdir -p "${OUT_DIR}"

echo "Requesting HELCOM subbasins zip from MADS GP service..."
RESPONSE=$(curl -fsSL "${GP_URL}?id=${RECORD_UUID}&f=json")

# Extract the url field from {"value": {"url": "..."}} inside results[0]
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
TMPZIP=$(mktemp /tmp/helcom_subbasins_XXXXXX.zip)
trap 'rm -f "${TMPZIP}"' EXIT

curl -fsSL -o "${TMPZIP}" "${DOWNLOAD_URL}"
unzip -o "${TMPZIP}" -d "${OUT_DIR}"

echo "HELCOM subbasins unpacked to: ${OUT_DIR}"
