#!/bin/bash
#SBATCH --job-name=debug_bottom_stationary
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00
#SBATCH --partition=base

# Scan one year's 2D bottom/surface fields and flag cells that never move.
#
# Usage:
#     sbatch debug_bottom_stationary_job.sh                   # default: fine, 2019
#     sbatch debug_bottom_stationary_job.sh coarse 2019
#     sbatch debug_bottom_stationary_job.sh fine   2020

module load gcc12-env/12.3.0
module load singularity/3.11.5

res="${1:-fine}"
year="${2:-2019}"

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
input_dir=${output_root}/2d_fields
output_dir=${output_root}/debug/bottom_stationary
container=parcels-container_2024.10.07-7af7fd0.sif

mkdir -p ${output_dir}

singularity run -B /sfs -B /gxfs_work -B ${repo_root}:/work --pwd /work \
    ${container} bash -c \
    ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
    && python scripts/debug_bottom_stationary.py \
        --input-dir ${input_dir} \
        --output-dir ${output_dir} \
        --year ${year} \
        --res ${res}"

jobinfo
