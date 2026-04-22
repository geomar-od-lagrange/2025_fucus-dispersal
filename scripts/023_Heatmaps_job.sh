#!/bin/bash
#SBATCH --job-name=023_Heatmaps
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=base

export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128

# Submit from the repo root.

experiment_type="${1:-surface}"

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/023_Heatmaps.ipynb \
    notebooks_executed/Visualisations/023_Heatmaps_${experiment_type}.ipynb \
    -p base_path ${base_path} \
    -p experiment_type ${experiment_type} \
    -k python

jobinfo
