#!/bin/bash
#SBATCH --job-name=037_MedianAgeGrid
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads 026's median-age GeoJSONs
# (Exports/026/026_medianage_<regime>_<season>_<tag>.geojson), no parquet.
# Draws median-age maps, columns with | without stranding, one row per
# regime, two figures of at most two rows each (one page): surface + Stokes
# drift and bottom; bottom with bed friction. Run after 026 at `member` and
# unweighted for every regime.
# Submit once per (season, hex_radius), e.g.:
#   sbatch scripts/037_MedianAgeGrid_job.sh ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0_b750
# season is one of DJF / MAM / JJA / SON / ALL. member is the 024f rate-model
# member of the "with stranding" column. hex_radius only labels the figure:
# the GeoJSON names carry none.

season="${1:-ALL}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0_b750}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal

mkdir -p notebooks_executed/Visualisations/

for rows in "surface_stokes,bottom" "bottom_ub49"; do
    pixi run papermill --cwd notebooks/ \
        notebooks/037_MedianAgeGrid.ipynb \
        notebooks_executed/Visualisations/037_MedianAgeGrid_${rows//,/_}_${season}_r${hex_radius}m.ipynb \
        -p data_root ${repo_root}/data \
        -p output_root ${output_root} \
        -p hex_radius ${hex_radius} \
        -r rows_csv ${rows} \
        -r member ${member} \
        -r season ${season} \
        -k python
done

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
