#!/bin/bash
#SBATCH --job-name=024a_BuildHexKey
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Single-process key file build. Run once per hex_radius, before any
# 024 counts jobs that share the radius — they read this file as a
# hard prerequisite.

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

hex_radius="${1:-6000}"

mkdir -p notebooks_executed/Visualisations/

pixi run papermill --cwd notebooks/ \
    notebooks/024a_BuildHexKey.ipynb \
    notebooks_executed/Visualisations/024a_BuildHexKey_r${hex_radius}m.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p hex_radius ${hex_radius} \
    -k python

jobinfo
