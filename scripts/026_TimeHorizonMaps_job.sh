#!/bin/bash
#SBATCH --job-name=026_TimeHorizonMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates counts parquet, no Dask cluster.
# Pools every release_year for the regime, so it reads all year partitions
# (hence more memory than 025's single-year read). Submit once per
# (regime, hex_radius), e.g.:
#   sbatch scripts/026_TimeHorizonMaps_job.sh surface 6000
#   sbatch scripts/026_TimeHorizonMaps_job.sh bottom 6000
# age_bin_days falls back to the notebook default (10) and must match the
# value 024 built the counts store with; horizons must be its multiples.

regime="${1:-surface}"
hex_radius="${2:-6000}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/026_TimeHorizonMaps.ipynb \
    notebooks_executed/Visualisations/026_TimeHorizonMaps_${regime}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
