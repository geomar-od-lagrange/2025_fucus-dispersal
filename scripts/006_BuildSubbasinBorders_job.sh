#!/bin/bash
#SBATCH --job-name=006_BuildSubbasinBorders
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight, static inputs only (HELCOM polygons, BSH coastline geojson and
# fine H0): writes output_root/SubbasinBorders/subbasin_borders.geojson, the
# subbasin borders every map notebook draws (025–037). Run once, before them.
#   sbatch scripts/006_BuildSubbasinBorders_job.sh

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Inputs/

pixi run papermill --cwd notebooks/ \
    notebooks/006_BuildSubbasinBorders.ipynb \
    notebooks_executed/Inputs/006_BuildSubbasinBorders.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
