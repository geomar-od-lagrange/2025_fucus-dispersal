#!/bin/bash
#SBATCH --job-name=026_TimeHorizonMaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=02:00:00
#SBATCH --partition=base

set -euo pipefail

# Lightweight: reads HexAggregates counts or survocc parquet, no Dask cluster.
# Pools the regime's partitions in the notebook's release_years and keeps one release season, then maps cumulative elapsed-time horizons.
# Submit once per (regime, season, hex_radius, member), e.g.:
#   sbatch scripts/026_TimeHorizonMaps_job.sh surface_stokes ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0
#   sbatch scripts/026_TimeHorizonMaps_job.sh bottom SON 6000
# season is one of DJF / MAM / JJA / SON / ALL. member "" (default) reads the
# unweighted 024 counts store; a 024f rate-model tag reads that member's
# survival-weighted survocc store (surface_stokes only: 024d sidecars).
# age_bin_days (10), max_age_days (180) and the horizons (30,60,90,180) fall
# back to the notebook defaults; age_bin_days must match the store, horizons
# must be its multiples, and a 024f store must cover max_age_days.

regime="${1:-surface_stokes}"
season="${2:-ALL}"
hex_radius="${3:-6000}"
member="${4:-}"
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

mkdir -p notebooks_executed/Visualisations/

# member is passed with -r (raw string): the empty string is a meaningful
# value (unweighted counts store) and -p would have papermill parse it away.
pixi run papermill --cwd notebooks/ \
    notebooks/026_TimeHorizonMaps.ipynb \
    notebooks_executed/Visualisations/026_TimeHorizonMaps_${regime}_${season}_r${hex_radius}m_${member:-unweighted}.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -r season ${season} \
    -p hex_radius ${hex_radius} \
    -r member "${member}" \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
