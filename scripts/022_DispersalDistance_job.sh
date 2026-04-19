#!/bin/bash
#SBATCH --job-name=022_DispersalDistance
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

# Submit from the repo root — relative paths and `-B $PWD:/work` assume it.

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
container=parcels-container_2024.10.07-7af7fd0.sif

mkdir -p notebooks_executed/Visualisations/

singularity run -B /sfs -B /gxfs_work -B $PWD:/work -B $TMPDIR:/tmp --pwd /work \
    ${container} bash -c \
    ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
    && papermill --cwd notebooks/ \
        notebooks/022_DispersalDistance.ipynb \
        notebooks_executed/Visualisations/022_DispersalDistance.ipynb \
        -p base_path ${base_path} \
        -k python"

jobinfo
