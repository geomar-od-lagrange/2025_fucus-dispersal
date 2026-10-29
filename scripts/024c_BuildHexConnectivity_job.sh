#!/bin/bash
#SBATCH --job-name=024c_BuildHexConnectivity
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/aggregates/%x_%j.out

# Lightweight: reads the 024 counts parquet + key, no Dask cluster.
# One (regime, release_year, hex_radius) per submit, e.g.:
#   sbatch scripts/024c_BuildHexConnectivity_job.sh surface_stokes 2019 6000
#   sbatch scripts/024c_BuildHexConnectivity_job.sh surface_stokes 2020 6000
# 024_BuildHexAggregates_job.sh and 024a_BuildHexKey_job.sh must have run
# first for the matching (regime, release_year, hex_radius).

# Without `set -e` a papermill failure here still reaches `jobinfo` and the
# job exits 0, so an `afterok` dependency downstream would fire on a store
# that was never written.
set -euo pipefail

regime="${1:-surface_stokes}"
year="${2:-2019}"
hex_radius="${3:-6000}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/024c_BuildHexConnectivity.ipynb \
    notebooks_executed/Visualisations/024c_BuildHexConnectivity_${regime}_${year}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p release_year ${year} \
    -p hex_radius ${hex_radius} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
