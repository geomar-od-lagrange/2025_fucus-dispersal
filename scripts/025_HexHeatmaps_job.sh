#!/bin/bash
#SBATCH --job-name=025_HexHeatmaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates parquet, no Dask cluster needed.
# Submit once per (regime, year), e.g.:
#   sbatch scripts/025_HexHeatmaps_job.sh surface 2019
#   sbatch scripts/025_HexHeatmaps_job.sh surface_stokes 2020
#   sbatch scripts/025_HexHeatmaps_job.sh bottom 2021

regime="${1:-surface}"
year="${2:-2019}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/025_HexHeatmaps.ipynb \
    notebooks_executed/Visualisations/025_HexHeatmaps_${regime}_${year}.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p release_year ${year} \
    -k python

jobinfo
