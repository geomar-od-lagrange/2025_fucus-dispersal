#!/bin/bash
#SBATCH --job-name=19h1d0s087Fucus_dispersal
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

release_year=2019
first_release_month=1
first_release_day=1
last_release_month=6
last_release_day=30
release_depth=0
relative_particle_speed=0.87
# make sure the output exists
mkdir -p notebooks_executed/
# run for single notebook and put into background
mkdir -p notebooks_executed/TrajectoryCalc/${release_year}/
mkdir -p output/Trajectories/${release_year}/
srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work parcels-container_2024.10.07-7af7fd0.sif bash -c \
". /opt/conda/etc/profile.d/conda.sh && conda activate base \
&& papermill --cwd notebooks/ \
    notebooks/FucusDispersal.ipynb \
    notebooks_executed/TrajectoryCalc/${release_year}/Fucus_y${release_year}_m${first_release_month}-${last_release_month}_d${release_depth}_relspeed${relative_particle_speed}.ipynb \
    -p release_year ${release_year} \
    -p first_release_month ${first_release_month} \
    -p first_release_day ${first_release_day} \
    -p last_release_month ${last_release_month} \
    -p last_release_day ${last_release_day} \
    -p max_age_days 220 \
    -p calc_dt_mins 15 \
    -p output_dt_mins 60 \
    -p relative_release_depth ${release_depth} \
    -p particles_per_cell 10 \
    -p repeated_release True \
    -p repeated_release_dt_days 7 \
    -p relative_particle_speed ${relative_particle_speed} \
    -p is_papermill True \
    -k python" &
    

# wait till background task is done
wait

# print resource infos
jobinfo
