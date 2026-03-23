#!/bin/bash
#SBATCH --job-name=Fucus_dispersal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=48:00:00
#SBATCH --partition=base


# make sure we have Singularity
module load gcc12-env/12.3.0
module load singularity/3.11.5

# to get the image (need to be on a partition which has internet access --> data), run
# $ singularity pull --disable-cache --dir "${PWD}" docker://quay.io/willirath/parcels-container:2024.10.07-7af7fd0

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
start_date=2019-01-01
release_depth=0
relative_particle_speed=0.87
# make sure the output dir exists
mkdir -p notebooks_executed/TrajectoryCalc/
srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work parcels-container_2024.10.07-7af7fd0.sif bash -c \
". /opt/conda/etc/profile.d/conda.sh && conda activate base \
&& papermill --cwd notebooks/ \
    notebooks/010_FucusDispersal.ipynb \
    notebooks_executed/TrajectoryCalc/Fucus_${start_date}_d${release_depth}_s${relative_particle_speed}.ipynb \
    -p start_date ${start_date} \
    -p max_age_days 220 \
    -p calc_dt_mins 5 \
    -p output_dt_mins 60 \
    -p relative_release_depth ${release_depth} \
    -p particles_per_cell 10 \
    -p relative_particle_speed ${relative_particle_speed} \
    -p base_path ${base_path} \
    -k python"

# print resource infos
jobinfo
