#!/bin/bash
#SBATCH --job-name=download_stokes
#SBATCH --partition=data
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-stokes-download-%j.out

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

cd ${repo_root}
# No --end-date: each product runs to the last day it covers (capped at
# yesterday, so the NRT forecast tail is not ingested as analysis).
pixi run python notebooks/002_download_stokes.py \
    --output-root ${output_root} \
    --start-date 2016-01-01
