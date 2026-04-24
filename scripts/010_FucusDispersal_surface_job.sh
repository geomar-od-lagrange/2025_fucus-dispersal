#!/bin/bash
#SBATCH --job-name=010_FucusSurface
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=24:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

# Submit from the repo root — relative paths and `-B $PWD:/work` assume it.

# Override release year via positional arg, e.g. `sbatch <script> 2024`.
year="${1:-2019}"

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
container=parcels-container_2024.10.07-7af7fd0.sif
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
            -p data_root ${repo_root}/data \
            -p output_root ${output_root} \
            -k python"
}
export -f run_experiment
export repo_root output_root container max_age_days calc_dt_mins output_dt_mins particles_per_cell

# Generate all (start_date, experiment_type, velocity_factor) combinations
# and run them in parallel via xargs, respecting SLURM task limit.
for doy in $(seq 1 5 366); do
    start_date=$(date -d "${year}-01-01 +$(( doy - 1 )) days" +%Y-%m-%d)

    printf '%s\0' "${start_date} surface 1.0 $((RANDOM * 32768 + RANDOM))"
    printf '%s\0' "${start_date} surface_stokes 1.0 $((RANDOM * 32768 + RANDOM))"
done | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _

jobinfo
