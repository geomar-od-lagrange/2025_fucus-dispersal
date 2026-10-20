#!/bin/bash
#SBATCH --job-name=036_SeedingMap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads the data twin's release cells and HELCOM subbasins, the
# 024a key, and the first two time steps of one regime-year of trajectory
# zarrs (the land-seeded share per cell), no Dask cluster. One figure per
# hex_radius, no season or member: the release set is the same for every run.
#   sbatch scripts/036_SeedingMap_job.sh 6000
# land_regime (surface_stokes), land_release_year (2019) and
# particles_per_cell (100) fall back to the notebook defaults.
# 024a_BuildHexKey_job.sh must have run first for hex_radius.

hex_radius="${1:-6000}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/036_SeedingMap.ipynb \
    notebooks_executed/Visualisations/036_SeedingMap_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
