#!/bin/bash
#SBATCH --job-name=034_BottomFriction
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# One counts partition read at a time (~10^8 rows, ~16 GiB peak as 026).
#SBATCH --mem-per-cpu=32G
#SBATCH --time=01:00:00
#SBATCH --partition=base

set -euo pipefail

# Parquet only: the 024 counts, 024e beaching and (weighted) 024f survocc
# partitions of `bottom` and `bottom_ub49` and the 024a key. Compares the two on the releases both hold,
# pooled over the release years of one season: area holding 90 % of
# particle-time per horizon (CSV + length-ratio figure), the map pair at 90 d,
# and the stranded % per horizon at the 024e `member` (CSV). Run after 024 and
# 024e (at `member`) for both regimes and every release year; without the
# 024e partitions the stranded table is reported missing and skipped. Submit
# once per (season, hex_radius, weighting), e.g.:
#   sbatch scripts/034_BottomFriction_job.sh ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0_b750 unweighted
#   sbatch scripts/034_BottomFriction_job.sh ALL 6000 step_wc0p1_ts3_tcinf_tf1_tw0_b750 step_wc0p1_ts3_tcinf_tf1_tw0_b750
# season is one of DJF / MAM / JJA / SON / ALL; weighting is "unweighted" (024
# counts) or a 024f member tag (survocc). scripts/025-031_ConsumersAll_job.sh
# runs it for every season and both weightings when bottom_ub49 is among its
# regimes.

season="${1:-ALL}"
hex_radius="${2:-6000}"
member="${3:-step_wc0p1_ts3_tcinf_tf1_tw0_b750}"
weighting="${4:-unweighted}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/034_BottomFriction.ipynb \
    notebooks_executed/Visualisations/034_BottomFriction_${season}_r${hex_radius}m_${weighting}.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -r season ${season} \
    -r member ${member} \
    -r weighting ${weighting} \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
