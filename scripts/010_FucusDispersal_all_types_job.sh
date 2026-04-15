#!/bin/bash
#SBATCH --job-name=010_FucusAll
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30G
#SBATCH --time=12:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
start_date=2016-01-01
velocity_factor=1.0

mkdir -p notebooks_executed/TrajectoryCalc/

for experiment_type in surface bottom surface_stokes; do
    srun --ntasks=1 --exact singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work parcels-container_2024.10.07-7af7fd0.sif bash -c \
    ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
    && papermill --cwd notebooks/ \
        notebooks/010_FucusDispersal.ipynb \
        notebooks_executed/TrajectoryCalc/Fucus_${start_date}_${experiment_type}_vf${velocity_factor}.ipynb \
        -p start_date ${start_date} \
        -p experiment_type ${experiment_type} \
        -p max_age_days 50 \
        -p calc_dt_mins 5 \
        -p output_dt_mins 60 \
        -p velocity_factor ${velocity_factor} \
        -p particles_per_cell 10 \
        -p path_2d_fields ${base_path}/output/2d_fields \
        -p path_release_locations ${base_path}/data/Fucus_location_shp \
        -p path_trajectories ${base_path}/output/Trajectories \
        -k python" &
done

wait
jobinfo
