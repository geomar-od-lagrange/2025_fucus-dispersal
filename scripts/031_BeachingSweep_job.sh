#!/bin/bash
#SBATCH --job-name=031_BeachingSweep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads one 024e beaching store per sweep member + the 024a key,
# no Dask cluster. One (regime, hex_radius, member list) per submit. E.g.:
#   sbatch scripts/031_BeachingSweep_job.sh surface_stokes 6000 \
#       step_wc0p05_ts3_tcinf,step_wc0p1_ts3_tcinf,step_wc0p15_ts3_tcinf
#   sbatch scripts/031_BeachingSweep_job.sh surface_stokes 6000 \
#       "step_wc0p05_ts3_tcinf,step_wc0p1_ts3_tcinf" "w_c 0.05,w_c 0.10" step_wc0p1_ts3_tcinf
# labels_csv (empty = use the member tags) and baseline_member (empty = none)
# are optional. 024e_BuildBeaching_job.sh must have run for every member, and
# disp_bin_km must match the value 024e built with (notebook default: 10).

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
members_csv="${3:-step_wc0p05_ts3_tcinf,step_wc0p1_ts3_tcinf,step_wc0p15_ts3_tcinf}"
labels_csv="${4:-}"
baseline_member="${5:-}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/031_BeachingSweep.ipynb \
    notebooks_executed/Visualisations/031_BeachingSweep_${regime}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -r members_csv "${members_csv}" \
    -r labels_csv "${labels_csv}" \
    -r baseline_member "${baseline_member}" \
    -k python

jobinfo
