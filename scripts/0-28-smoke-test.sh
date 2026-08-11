#!/bin/bash
# Full pipeline smoke test (stages 000–028) against the BSH demo subset.
#
# Runs every numbered notebook end-to-end on a single host (login node
# is fine — no SLURM). Designed to flush out integration issues from a
# fresh clone, not to produce scientifically meaningful trajectories:
# uses a 4 h Parcels window, 1 particle per release cell, all three
# regimes (surface, bottom, surface_stokes).
#
# Usage from the repo root:
#     ./scripts/0-28-smoke-test.sh
#
# Override outputs location:
#     OUTPUT_ROOT=/work/<user>/fucus_smoke ./scripts/0-28-smoke-test.sh
#
# Prerequisites (one-time, on the host):
#     git clone --recurse-submodules https://github.com/geomar-od-lagrange/2025_fucus-dispersal.git
#     cd 2025_fucus-dispersal
#     git -C data lfs pull
#     pixi install
#     copernicusmarine login    # for stage 002

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "${REPO_ROOT}"

OUTPUT_ROOT="${OUTPUT_ROOT:-${REPO_ROOT}/output_smoke}"
mkdir -p "${OUTPUT_ROOT}"

REGIMES=(surface bottom surface_stokes)

echo "=== Smoke test: stages 000–028 ==="
echo "Repo:        ${REPO_ROOT}"
echo "Output root: ${OUTPUT_ROOT}"
echo "Regimes:     ${REGIMES[*]}"

if [ ! -f data/helcom_fucus_redlist/REDLIST_SIS_Macrophytes.shp ]; then
    echo "FAIL: data/ submodule not populated. Run:"
    echo "    git submodule update --init data && git -C data lfs pull"
    exit 1
fi

echo
echo "==== 000 FucusStartLocations ===="
pixi run jupytext --sync --execute notebooks/000_FucusStartLocations.md

echo
echo "==== 002 download_stokes (2020-01-01, both products) ===="
pixi run python notebooks/002_download_stokes.py \
    --output-root "${OUTPUT_ROOT}" --year 2020 --month 1 --day 1

echo
echo "==== 003 prepare_2d_fields (first c-file of each grid) ===="
# The 4 h smoke window (01:00–05:00) sits inside the 00 UTC c-file
# (covers 00:00–06:00); only one c-file per grid is needed.
for c in data/bsh_hbmnoku_demo/c_file_fine_2020/c_file_fine_2020010100_*.nc \
         data/bsh_hbmnoku_demo/c_file_coarse_2020/c_file_coarse_2020010100_*.nc; do
    pixi run python notebooks/003_prepare_2d_fields.py \
        --c-file "${c}" --output-root "${OUTPUT_ROOT}"
done

# 004 (extract_coastline): skipped — output already in the data twin.

for regime in "${REGIMES[@]}"; do
    echo
    echo "==== 010 FucusDispersal (${regime}, 4 h, 1 particle/cell) ===="
    pixi run papermill notebooks/010_FucusDispersal.ipynb \
        "${OUTPUT_ROOT}/010_${regime}_smoke.ipynb" \
        -p start_time "2020-01-01T01:00:00" \
        -p end_time "2020-01-01T05:00:00" \
        -p regime "${regime}" \
        -p output_root "${OUTPUT_ROOT}" \
        -p particles_per_cell 1 \
        -p RNG_seed 42 \
        -p allow_time_extrapolation True \
        --cwd notebooks/
done

echo
echo "==== 020 RawTrajectories ===="
pixi run papermill notebooks/020_RawTrajectories.ipynb \
    "${OUTPUT_ROOT}/020_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    --cwd notebooks/

echo
echo "==== 021 TimeStats ===="
pixi run papermill notebooks/021_TimeStats.ipynb \
    "${OUTPUT_ROOT}/021_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    --cwd notebooks/

echo
echo "==== 022 DispersalDistance ===="
pixi run papermill notebooks/022_DispersalDistance.ipynb \
    "${OUTPUT_ROOT}/022_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    --cwd notebooks/

# 023 and 025 are per-regime by design; smoke verifies the code path
# with one regime (surface_stokes — the regime that exercises the new
# WAVERYS-layered Stokes preprocessing end-to-end).
echo
echo "==== 023 Heatmaps (surface_stokes) ===="
pixi run papermill notebooks/023_Heatmaps.ipynb \
    "${OUTPUT_ROOT}/023_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime "surface_stokes" \
    --cwd notebooks/

echo
echo "==== 024a BuildHexKey ===="
pixi run papermill notebooks/024a_BuildHexKey.ipynb \
    "${OUTPUT_ROOT}/024a_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    --cwd notebooks/

echo
echo "==== 024 BuildHexAggregates (surface_stokes) ===="
pixi run papermill notebooks/024_BuildHexAggregates.ipynb \
    "${OUTPUT_ROOT}/024_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime surface_stokes \
    -p release_year 2020 \
    --cwd notebooks/

echo
echo "==== 024b BuildHexDistance (surface_stokes) ===="
pixi run papermill notebooks/024b_BuildHexDistance.ipynb \
    "${OUTPUT_ROOT}/024b_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime surface_stokes \
    -p release_year 2020 \
    --cwd notebooks/

echo
echo "==== 024c BuildHexConnectivity (surface_stokes) ===="
pixi run papermill notebooks/024c_BuildHexConnectivity.ipynb \
    "${OUTPUT_ROOT}/024c_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime surface_stokes \
    -p release_year 2020 \
    --cwd notebooks/

echo
echo "==== 025 HexHeatmaps (surface_stokes) ===="
pixi run papermill notebooks/025_HexHeatmaps.ipynb \
    "${OUTPUT_ROOT}/025_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime "surface_stokes" \
    -p release_year 2020 \
    --cwd notebooks/

# 026/027 pool releases across all months/years; the smoke run only has a
# January release, so pass release_months_csv="" (all months). 026 horizons
# must be age_bin multiples that the 4 h run actually populates — only
# age_bin 0 exists, so map the single "0 d" horizon. 027 has 1 particle per
# cell, so drop the min-trajectory gate.
echo
echo "==== 026 TimeHorizonMaps (surface_stokes) ===="
pixi run papermill notebooks/026_TimeHorizonMaps.ipynb \
    "${OUTPUT_ROOT}/026_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime "surface_stokes" \
    -p release_months_csv "" \
    -p time_horizons_days_csv "0" \
    --cwd notebooks/

echo
echo "==== 027 HexDistanceQuantiles (surface_stokes) ===="
pixi run papermill notebooks/027_HexDistanceQuantiles.ipynb \
    "${OUTPUT_ROOT}/027_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime "surface_stokes" \
    -p release_months_csv "" \
    -p min_traj_per_hex 1 \
    --cwd notebooks/

# 028 pools every release_year for the regime; the smoke run has only a
# January release, so its aug_sep scope is empty and the notebook skips it
# (the all_year scope carries the single release). One regime, like 023–027.
echo
echo "==== 028 SubbasinConnectivityMatrix (surface_stokes) ===="
pixi run papermill notebooks/028_SubbasinConnectivityMatrix.ipynb \
    "${OUTPUT_ROOT}/028_smoke.ipynb" \
    -p output_root "${OUTPUT_ROOT}" \
    -p regime "surface_stokes" \
    --cwd notebooks/

echo
echo "=== Smoke test complete ==="
echo "Outputs:             ${OUTPUT_ROOT}"
echo "Executed notebooks:  ${OUTPUT_ROOT}/0??_*_smoke.ipynb"
