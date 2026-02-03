#!/bin/bash
#SBATCH --job-name=midi1_Fucus_dispersal
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=18:00:00
#SBATCH --partition=base


# make sure we have Singularity
module load gcc12-env/12.3.0
module load singularity/3.11.5

# to get the image (need to be on a partition which has internet access --> data), run
# $ singularity pull --disable-cache --dir "${PWD}" docker://quay.io/willirath/parcels-container:2024.10.07-7af7fd0

release_year=2017
release_depth=1
# make sure the output exists
mkdir -p notebooks_executed/
# run for single notebook and put into background
mkdir -p notebooks_executed/TrajectoryCalc/${release_year}/
mkdir -p output/Trajectories/${release_year}/
srun --ntasks=1 --exclusive singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work parcels-container_2024.10.07-7af7fd0.sif bash -c \
". /opt/conda/etc/profile.d/conda.sh && conda activate base \
&& papermill --cwd notebooks/ \
    notebooks/FucusDispersal.ipynb \
    notebooks_executed/TrajectoryCalc/${release_year}/Fucus_${release_year}_d${release_depth}.ipynb \
    -p release_year ${release_year} \
    -p first_release_month 1 \
    -p first_release_day 1 \
    -p last_release_month 12 \
    -p last_release_day 31 \
    -p max_age_d 220 \
    -p dt_in_minutes 60 \
    -p output_dt_in_minutes 60 \
    -p release_depth_sigma ${release_depth} \
    -p n_particles_per_cell 10 \
    -p repeated_release True \
    -p repeatdt_d 7 \
    -p is_papermill True \
    -k python" &
    

# wait till background task is done
wait

# print resource infos
jobinfo
