#!/bin/bash
#SBATCH --job-name=035_AreaTimeSeries
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads 026's area CSVs (Exports/026/026_areas_*.csv) and 029's
# stranding curves (Exports/029/stranded_curve_*.csv), no parquet. Draws the
# stranded share, the area holding 90 % of particle-time and its square root
# against the cumulative age, one line per regime, for one release season.
# Run after 026 and 029 at `member` for every regime it draws (a regime whose
# area CSV or stranding curve is missing is reported and left out of that
# panel; none at all raises).
# Submit once per (season, hex_radius), e.g.:
#   sbatch scripts/035_AreaTimeSeries_job.sh ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0_b750
# season is one of DJF / MAM / JJA / SON / ALL. member is the rate-model
# member shared by every regime: their 026 weighting tag (the 024f production
# member) and the 024e member of their stranding curves. hex_radius only
# labels the figure: the CSV names carry none.

season="${1:-ALL}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0_b750}"
# "regimes" (one season, one line per regime) or "seasons" (the four release
# seasons as lines, one row per regime; season is then ignored).
layout="${4:-regimes}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/035_AreaTimeSeries.ipynb \
    notebooks_executed/Visualisations/035_AreaTimeSeries_${layout}_${season}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -r regimes_csv "surface_stokes,bottom,bottom_ub49" \
    -r tags_csv "${member},${member},${member}" \
    -r stranded_member ${member} \
    -r season ${season} \
    -r layout ${layout} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
