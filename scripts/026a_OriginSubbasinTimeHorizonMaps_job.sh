#!/bin/bash
#SBATCH --job-name=026a_OriginSubbasinTimeHorizonMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates counts parquet, no Dask cluster.
# Same as 026, emitted once per origin (HELCOM) subbasin.
# Submit once per (regime, season, hex_radius), e.g.:
#   sbatch scripts/026a_OriginSubbasinTimeHorizonMaps_job.sh surface_stokes ALL 6000
#   sbatch scripts/026a_OriginSubbasinTimeHorizonMaps_job.sh bottom SON 6000
# season is one of DJF / MAM / JJA / SON / ALL. age_bin_days falls back to
# the notebook default (10) and must match the value 024 built the counts
# store with; horizons must be its multiples.
# origin_subbasins_csv falls back to the notebook default (all named
# subbasins that seed releases).

regime="${1:-surface_stokes}"
season="${2:-ALL}"
hex_radius="${3:-6000}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/026a_OriginSubbasinTimeHorizonMaps.ipynb \
    notebooks_executed/Visualisations/026a_OriginSubbasinTimeHorizonMaps_${regime}_${season}_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p season ${season} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
