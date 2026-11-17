#!/bin/bash
#SBATCH --job-name=038_EmissionGrid
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:15:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads 028's emission-fraction CSVs
# (Exports/028/connectivity_<regime>_<season>_<tag>_T<h>d_emission_fraction.csv),
# no parquet. Draws one 2×2 figure per regime: columns with stranding (the
# 024g member) | without stranding (unweighted), rows age < 30 d | < 180 d.
# Run after 028 at that season for both the member and the unweighted run.
# Submit once per regime, e.g.:
#   sbatch scripts/038_EmissionGrid_job.sh surface_stokes ALL step_wc0p1_ts3_tcinf_tf1_tw0_b750
# season is one of DJF / MAM / JJA / SON / ALL.

regime="${1:-surface_stokes}"
season="${2:-ALL}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0_b750}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/038_EmissionGrid.ipynb \
    notebooks_executed/Visualisations/038_EmissionGrid_${regime}_${season}.ipynb \
    -p output_root ${output_root} \
    -r regime ${regime} \
    -r season ${season} \
    -r horizons_csv "30,180" \
    -r member ${member} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
