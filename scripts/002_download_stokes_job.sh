#!/bin/bash
#SBATCH --job-name=download_stokes
#SBATCH --partition=data
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --output=slurm-stokes-download-%j.out

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
stokes_dir=/gxfs_work/geomar/smomw122/cmems_mod_bal_wav
# Set up temporary pixi env
pixi init --format pixi /tmp/stokes_env_$$
cd /tmp/stokes_env_$$
pixi add copernicusmarine

# Run download
pixi run python ${base_path}/scripts/002_download_stokes.py \
    --output-dir ${stokes_dir} \
    --start-year 2016 \
    --end-year 2025

# Clean up
rm -rf /tmp/stokes_env_$$
