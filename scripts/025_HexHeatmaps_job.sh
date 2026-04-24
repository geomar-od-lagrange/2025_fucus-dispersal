#!/bin/bash
#SBATCH --job-name=025_HexHeatmaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates parquet, no Dask cluster needed.
# Submit once per regime, e.g.:
#   sbatch scripts/025_HexHeatmaps_job.sh surface
#   sbatch scripts/025_HexHeatmaps_job.sh surface_stokes
#   sbatch scripts/025_HexHeatmaps_job.sh bottom

regime="${1:-surface}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/025_HexHeatmaps.ipynb \
    notebooks_executed/Visualisations/025_HexHeatmaps_${regime}.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -k python

jobinfo
