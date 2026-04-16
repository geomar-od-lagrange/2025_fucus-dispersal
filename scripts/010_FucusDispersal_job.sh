#!/bin/bash
#SBATCH --job-name=010_FucusProd
#SBATCH --ntasks=50
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=24:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
container=parcels-container_2024.10.07-7af7fd0.sif
year=2019
max_age_days=220
calc_dt_mins=5
output_dt_mins=60
particles_per_cell=100

mkdir -p notebooks_executed/TrajectoryCalc/

run_experiment() {
    local start_date="$1"
    local experiment_type="$2"
    local velocity_factor="$3"
    local rng_seed="$4"

    srun --ntasks=1 --exact \
        singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work \
        ${container} bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
        && papermill --cwd notebooks/ \
            notebooks/010_FucusDispersal.ipynb \
            notebooks_executed/TrajectoryCalc/Fucus_${start_date}_${experiment_type}_vf${velocity_factor}_seed${rng_seed}.ipynb \
            -p start_date ${start_date} \
            -p experiment_type ${experiment_type} \
            -p RNG_seed ${rng_seed} \
            -p max_age_days ${max_age_days} \
            -p calc_dt_mins ${calc_dt_mins} \
            -p output_dt_mins ${output_dt_mins} \
            -p velocity_factor ${velocity_factor} \
            -p particles_per_cell ${particles_per_cell} \
            -p path_2d_fields ${base_path}/output/2d_fields \
            -p path_release_locations ${base_path}/data/Fucus_location_shp \
            -p path_trajectories ${base_path}/output/Trajectories \
            -k python"
}
export -f run_experiment
export base_path container max_age_days calc_dt_mins output_dt_mins particles_per_cell

# Generate all (start_date, experiment_type, velocity_factor) combinations
# and run them in parallel via xargs, respecting SLURM task limit.
# Each record is one null-delimited line: "start_date experiment_type vf seed"
for doy in $(seq 1 5 366); do
    start_date=$(date -d "${year}-01-01 +$(( doy - 1 )) days" +%Y-%m-%d)

    printf '%s\0' "${start_date} surface 1.0 ${RANDOM}"
    printf '%s\0' "${start_date} bottom 1.0 ${RANDOM}"
    printf '%s\0' "${start_date} bottom 0.97 ${RANDOM}"
    printf '%s\0' "${start_date} bottom 0.87 ${RANDOM}"
    printf '%s\0' "${start_date} surface_stokes 1.0 ${RANDOM}"
done | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _

jobinfo
