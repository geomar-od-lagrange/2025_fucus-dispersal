#!/bin/bash
#SBATCH --job-name=027_HexDistanceQuantiles
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Lightweight: reads from HexAggregates distance parquet, no Dask cluster.
# Pools every release_year for the regime (distance histograms are small).
# Submit once per (regime, hex_radius, season), e.g.:
#   sbatch scripts/027_HexDistanceQuantiles_job.sh surface_stokes 6000 ALL
#   sbatch scripts/027_HexDistanceQuantiles_job.sh bottom 6000 SON
# season is one of DJF / MAM / JJA / SON / ALL (release month).
# distance_bin_km falls back to the notebook default (1.0 km) and must
# match the value 024b built the store with.

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
season="${3:-ALL}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/027_HexDistanceQuantiles.ipynb \
    notebooks_executed/Visualisations/027_HexDistanceQuantiles_${regime}_r${hex_radius}m_${season}.ipynb \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -p hex_radius ${hex_radius} \
    -r season ${season} \
    -k python

jobinfo
