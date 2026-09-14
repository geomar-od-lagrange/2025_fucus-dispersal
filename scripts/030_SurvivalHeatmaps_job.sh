#!/bin/bash
#SBATCH --job-name=030_SurvivalHeatmaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads the 024f survival-occupancy parquet + the 024a key, no
# Dask cluster. One (regime, hex_radius, member, season) per submit; the
# notebook pools the season's monthly partitions across every release year. E.g.:
#   sbatch scripts/030_SurvivalHeatmaps_job.sh surface_stokes 6000 step_wc0p1_ts3_tcinf ALL
# season is one of DJF / MAM / JJA / SON / ALL (release month).
# 024f_BuildSurvivalOccupancy_job.sh must have run first for the matching
# member; age_bin_days falls back to the notebook default (10) and must match
# the value 024f built with.

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf}"
season="${4:-ALL}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/030_SurvivalHeatmaps.ipynb \
    notebooks_executed/Visualisations/030_SurvivalHeatmaps_${regime}_r${hex_radius}m_${season}_${member}.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -r member ${member} \
    -r season ${season} \
    -k python

jobinfo
