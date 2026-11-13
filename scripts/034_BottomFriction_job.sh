#!/bin/bash
#SBATCH --job-name=034_BottomFriction
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# One counts partition read at a time (~10^8 rows, ~16 GiB peak as 026).
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --partition=base

set -euo pipefail

# Parquet only: the 024 counts partitions of `bottom` and `bottom_ub49` and the
# 024a key. Compares the two on the releases both hold, pooled over the
# release years of one season: area holding 90 % of particle-time per horizon
# (CSV + length-ratio figure) and the map pair at 90 d. Run after 024 for both
# regimes and every release year. Submit once per (season, hex_radius), e.g.:
#   sbatch scripts/034_BottomFriction_job.sh ALL 6000
# season is one of DJF / MAM / JJA / SON / ALL. scripts/025-031_ConsumersAll_job.sh
# runs it for every season when bottom_ub49 is among its regimes.

season="${1:-ALL}"
hex_radius="${2:-6000}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/034_BottomFriction.ipynb \
    notebooks_executed/Visualisations/034_BottomFriction_${season}_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -r season ${season} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
