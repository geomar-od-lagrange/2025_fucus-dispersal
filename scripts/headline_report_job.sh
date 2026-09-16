#!/bin/bash
#SBATCH --job-name=032_HeadlineReport
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/032/%x_%j.out

set -euo pipefail

# Regenerate output_root/Headline/ -- HEADLINE.md and the numbered, downscaled
# copies of the selected manuscript figures.
#
# Lightweight and single-task: the notebook reads the 024x parquet partitions,
# the 024a key and the 025-031 GeoJSON/CSV exports, and draws nothing. Minutes,
# no Dask cluster. Release years, 024e member tags and figure filenames are all
# discovered at run time, so this takes no year and no member list:
#
#   sbatch scripts/headline_report_job.sh
#   sbatch scripts/headline_report_job.sh surface_stokes,bottom step_wc0p1_ts3_tcinf
#
# Last stage of scripts/submit_downstream_chain.sh, which submits it
# afterok on the 025-031 consumer job -- every Exports/ product the tables and
# the figure selection read is written there.

cd "${SLURM_SUBMIT_DIR:-$PWD}"

regimes_csv="${1:-surface_stokes,bottom}"
member="${2:-step_wc0p1_ts3_tcinf}"
hex_radius="${3:-6000}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

# logs/032 is this job's own --output dir; SLURM opens that file before the
# script runs, so it must already exist (scripts/submit_downstream_chain.sh
# pre-creates it). Repeat it here so a direct `sbatch` on a fresh output_root
# fails only once, not silently with exit 53 and no log.
mkdir -p "${output_root}/logs/032" notebooks_executed/Visualisations/

# Blacklisted nodes (see scripts/node-blacklist.txt), for the provenance line.
# The exclusion itself is applied to the allocation by the submitter: an
# #SBATCH directive cannot read a file.
exclude=$(sed "s/#.*//" scripts/node-blacklist.txt | tr -d "[:blank:]" \
          | grep -v "^$" | paste -sd,)
echo "headline report: regimes=${regimes_csv} member=${member} r=${hex_radius}m"
echo "node $(hostname); blacklist=${exclude:-none}"

pixi run papermill --cwd notebooks/ \
    notebooks/032_HeadlineReport.ipynb \
    notebooks_executed/Visualisations/032_HeadlineReport_r${hex_radius}m.ipynb \
    -p output_root "${output_root}" \
    -p hex_radius "${hex_radius}" \
    -r regimes_csv "${regimes_csv}" \
    -r member "${member}" \
    -k python

# The job status must be papermill's, not the accounting tool's.
jobinfo || true
