#!/bin/bash
#SBATCH --job-name=028_SubbasinConnectivityMatrix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates connectivity parquet (024c store),
# no Dask cluster. Pools every available release_year for the regime, then
# prints the subbasin→subbasin residence matrix and emits linear + log
# heatmap PNGs for two release-month scopes (all-year and Aug/Sep).
# Submit once per (regime, hex_radius), e.g.:
#   sbatch scripts/028_SubbasinConnectivityMatrix_job.sh surface 6000
#   sbatch scripts/028_SubbasinConnectivityMatrix_job.sh bottom 6000
# age_bin_days falls back to the notebook default (10) and must match the
# value 024c built the connectivity store with.

regime="${1:-surface}"
hex_radius="${2:-6000}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/028_SubbasinConnectivityMatrix.ipynb \
    notebooks_executed/Visualisations/028_SubbasinConnectivityMatrix_${regime}_r${hex_radius}m.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
