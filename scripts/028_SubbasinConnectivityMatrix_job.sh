#!/bin/bash
#SBATCH --job-name=028_SubbasinConnectivityMatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads parquet only (the 024c connectivity store, or the 024g
# survconn store when a member is given), no Dask cluster. Pools every
# available release_year for the regime, restricts to one release season, and
# emits the raw / emission-fraction / interannual-range matrices as print-ready
# PNGs under Figures/028 plus CSVs under Exports/028.
#
# Usage: sbatch scripts/028_SubbasinConnectivityMatrix_job.sh \
#            [regime] [hex_radius] [season] [member]
#   season: DJF | MAM | JJA | SON | ALL   (default ALL)
#   member: "" -> unweighted 024c store (default); a 024g rate-model tag such
#           as step_wc0p1_ts3_tcinf -> survival-weighted w_obs. Only
#           surface_stokes has 024d sidecars, hence 024g members.
# e.g.:
#   sbatch scripts/028_SubbasinConnectivityMatrix_job.sh surface_stokes 6000 SON ""
#   sbatch scripts/028_SubbasinConnectivityMatrix_job.sh surface_stokes 6000 ALL step_wc0p1_ts3_tcinf
#   sbatch scripts/028_SubbasinConnectivityMatrix_job.sh bottom 6000 JJA ""
# age_bin_days and time_horizons_days_csv fall back to the notebook defaults;
# age_bin_days must match the value 024c/024g built the store with.

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
season="${3:-ALL}"
member="${4:-}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

# member is passed with -r (raw string): the empty string is a meaningful
# value (unweighted store) and -p would have papermill parse it away.
pixi run papermill --cwd notebooks/ \
    notebooks/028_SubbasinConnectivityMatrix.ipynb \
    notebooks_executed/Visualisations/028_SubbasinConnectivityMatrix_${regime}_r${hex_radius}m_${season}_${member:-unweighted}.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -r season ${season} \
    -r member "${member}" \
    -k python

jobinfo
