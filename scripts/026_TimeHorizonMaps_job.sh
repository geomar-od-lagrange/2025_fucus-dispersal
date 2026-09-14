#!/bin/bash
#SBATCH --job-name=026_TimeHorizonMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates counts parquet, no Dask cluster.
# Pools every release_year partition of the regime and keeps one release season, then maps cumulative elapsed-time horizons.
# Submit once per (regime, season, hex_radius), e.g.:
#   sbatch scripts/026_TimeHorizonMaps_job.sh surface_stokes ALL 6000
#   sbatch scripts/026_TimeHorizonMaps_job.sh bottom SON 6000
# season is one of DJF / MAM / JJA / SON / ALL. age_bin_days falls back to
# the notebook default (10) and must match the value 024 built the counts
# store with; horizons must be its multiples.

regime="${1:-surface_stokes}"
season="${2:-ALL}"
hex_radius="${3:-6000}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/026_TimeHorizonMaps.ipynb \
    notebooks_executed/Visualisations/026_TimeHorizonMaps_${regime}_${season}_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p season ${season} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
