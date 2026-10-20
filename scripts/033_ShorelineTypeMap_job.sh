#!/bin/bash
#SBATCH --job-name=033_ShorelineTypeMap
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads the data twin's shoreline-class table, the 024a key and
# the 024e beaching parquet, no Dask cluster. One (regime, hex_radius, member,
# season) per submit; the notebook pools the season's monthly partitions across
# the notebook's release_years and repeats the relationship for `contrast_member`. E.g.:
#   sbatch scripts/033_ShorelineTypeMap_job.sh surface_stokes 6000 \
#       step_wc0p1_ts3_tcinf_tf1_tw0 ALL step_wc0p1_ts3_tcinf
# season is one of DJF / MAM / JJA / SON / ALL (release month). An explicitly
# empty 5th argument ("") switches the contrast off; omitting it takes the
# default contrast member.
# 024e_BuildBeaching_job.sh must have run first for both members.

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0}"
season="${4:-ALL}"
contrast_member="${5-step_wc0p1_ts3_tcinf}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/033_ShorelineTypeMap.ipynb \
    notebooks_executed/Visualisations/033_ShorelineTypeMap_${regime}_r${hex_radius}m_${season}_${member}.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -r member ${member} \
    -r season ${season} \
    -r contrast_member "${contrast_member}" \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
