#!/bin/bash
#SBATCH --job-name=025_HexHeatmaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads from HexAggregates parquet, no Dask cluster needed.
# Pools every release_year partition of the regime and keeps one release
# season. Submit once per (regime, season, hex_radius), e.g.:
#   sbatch scripts/025_HexHeatmaps_job.sh surface_stokes ALL 6000
#   sbatch scripts/025_HexHeatmaps_job.sh bottom SON 6000
# season is one of DJF / MAM / JJA / SON / ALL.

regime="${1:-surface_stokes}"
season="${2:-ALL}"
hex_radius="${3:-6000}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/025_HexHeatmaps.ipynb \
    notebooks_executed/Visualisations/025_HexHeatmaps_${regime}_${season}_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -r season ${season} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
