#!/bin/bash
#SBATCH --job-name=021_TimeStats
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=base

export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128

# Submit from the repo root.

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/021_TimeStats.ipynb \
    notebooks_executed/Visualisations/021_TimeStats.ipynb \
    -p base_path ${base_path} \
    -k python

jobinfo
