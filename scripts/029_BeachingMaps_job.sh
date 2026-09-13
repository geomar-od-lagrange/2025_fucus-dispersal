#!/bin/bash
#SBATCH --job-name=029_BeachingMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Lightweight: reads the 024e beaching parquet + the 024a key, no Dask cluster.
# One (regime, hex_radius, rate-model member) per submit; the notebook pools
# every release year and month of that member. E.g.:
#   sbatch scripts/029_BeachingMaps_job.sh surface_stokes 6000 step_wc0p1_ts3_tcinf
# 024e_BuildBeaching_job.sh must have run first for the matching member, and
# disp_bin_km must match the value 024e built with (notebook default: 10).

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/029_BeachingMaps.ipynb \
    notebooks_executed/Visualisations/029_BeachingMaps_${regime}_r${hex_radius}m_${member}.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -p member ${member} \
    -k python

jobinfo
