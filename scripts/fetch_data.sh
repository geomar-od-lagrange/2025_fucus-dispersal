#!/usr/bin/env bash
# Upstream rebuild: reconstruct ./data/ from public sources.
#
# This is NOT the primary setup path.  The normal way to get the input data is
# via the git submodule:
#
#   git clone --recurse-submodules <repo-url>
#   git -C data lfs pull
#
# Use THIS script when the submodule is unreachable, when the recipe has
# changed and you want to regenerate the bundled blobs locally, or when
# running the data twin's CI to keep its blobs in sync with the recipe.
#
# Each step lands its output under ./data/ (relative to repo root).
# All sources are public; no NESH mount or other auth required beyond
# git-lfs (for the twin's LFS payload, installed separately) and optionally
# copernicusmarine credentials for the Stokes step.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "=== fetch_data: HELCOM subbasins ==="
bash "${SCRIPT_DIR}/obtain/download_helcom_subbasins.sh"

echo "=== fetch_data: Fucus REDLIST shapefile ==="
bash "${SCRIPT_DIR}/obtain/download_fucus_shapefile.sh"

echo "=== fetch_data: CMEMS Stokes sample ==="
bash "${SCRIPT_DIR}/obtain/download_stokes_sample.sh"

echo "=== fetch_data: BSH HBMnoku demo subset (via BSH OpenData) ==="
bash "${SCRIPT_DIR}/obtain/download_bsh_hbmnoku_demo.sh"

echo "=== fetch_data: bake derived release points ==="
pixi run jupytext --sync --execute notebooks/000_FucusStartLocations.md

echo "=== fetch_data: done ==="
