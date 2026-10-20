#!/bin/bash
#SBATCH --job-name=035_AreaTimeSeries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads 026's area CSVs (Exports/026/026_areas_*.csv), no parquet.
# Draws the area holding 90 % of particle-time against the cumulative horizon,
# one line per regime, for one release season. Run after 026 for every regime
# it draws; a regime whose CSV is missing is reported and left out.
# Submit once per (season, hex_radius), e.g.:
#   sbatch scripts/035_AreaTimeSeries_job.sh ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0
# season is one of DJF / MAM / JJA / SON / ALL. member is the 026 weighting
# tag of surface_stokes (its 024f production member); bottom is always read
# unweighted. hex_radius only labels the figure: the 026 CSV names carry none.

season="${1:-ALL}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/035_AreaTimeSeries.ipynb \
    notebooks_executed/Visualisations/035_AreaTimeSeries_${season}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -r regimes_csv "surface_stokes,bottom" \
    -r tags_csv "${member},unweighted" \
    -r season ${season} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
