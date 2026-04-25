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
simulation_days=220
calc_dt_mins=5
output_dt_mins=60
particles_per_cell=100

mkdir -p notebooks_executed/TrajectoryCalc/

run_experiment() {
    local start_time="$1"
    local end_time="$2"
    local regime="$3"
    local rng_seed="$4"
    local start_stem="${start_time//[-:]/}"

    srun --ntasks=1 --exact \
        singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work \
        ${container} bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
        && papermill --cwd notebooks/ \
            notebooks/010_FucusDispersal.ipynb \
            notebooks_executed/TrajectoryCalc/Fucus_${start_stem}_${regime}.ipynb \
            -p start_time ${start_time} \
            -p end_time ${end_time} \
            -p regime ${regime} \
            -p RNG_seed ${rng_seed} \
            -p calc_dt_mins ${calc_dt_mins} \
            -p output_dt_mins ${output_dt_mins} \
            -p particles_per_cell ${particles_per_cell} \
            -p data_root ${repo_root}/data \
            -p output_root ${output_root} \
            -k python"
}
export -f run_experiment
export repo_root output_root container calc_dt_mins output_dt_mins particles_per_cell

# 73 releases per year at doy = 1 + 5n for n ∈ [0, 72] (doys 1..361,
# leap-year-agnostic). One papermill per start_time.
for n in $(seq 0 72); do
    doy=$((1 + 5 * n))
    start_time=$(date -d "${year}-01-01 +$((doy - 1)) days" +%Y-%m-%dT00:15:00)
    end_time=$(date -d "${start_time} +${simulation_days} days" +%Y-%m-%dT00:15:00)

    printf '%s\0' "${start_time} ${end_time} surface $((RANDOM * 32768 + RANDOM))"
done | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _

jobinfo
