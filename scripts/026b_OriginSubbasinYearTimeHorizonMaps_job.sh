#!/bin/bash
#SBATCH --job-name=026b_OriginSubbasinYearTimeHorizonMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads from HexAggregates counts parquet, no Dask cluster.
# Same as 026a, split a second time by release year.
# Submit once per (regime, season, hex_radius), e.g.:
#   sbatch scripts/026b_OriginSubbasinYearTimeHorizonMaps_job.sh surface_stokes ALL 6000
#   sbatch scripts/026b_OriginSubbasinYearTimeHorizonMaps_job.sh bottom SON 6000
# season is one of DJF / MAM / JJA / SON / ALL. age_bin_days falls back to
# the notebook default (10) and must match the value 024 built the counts
# store with; horizons must be its multiples.
# origin_subbasins_csv and release_years_csv fall back to the notebook
# defaults (all named subbasins that seed releases; all years present).

regime="${1:-surface_stokes}"
season="${2:-ALL}"
hex_radius="${3:-6000}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/026b_OriginSubbasinYearTimeHorizonMaps.ipynb \
    notebooks_executed/Visualisations/026b_OriginSubbasinYearTimeHorizonMaps_${regime}_${season}_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -r season ${season} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
