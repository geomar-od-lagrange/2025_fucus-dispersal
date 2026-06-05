#!/bin/bash
#SBATCH --job-name=010_FucusGapfill
#SBATCH --ntasks=118
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=60G
#SBATCH --time=24:00:00
#SBATCH --partition=base

# One-off backfill of the 118 trajectory releases missing from the store as
# of the 2026-06-05 audit. --ntasks matches the release count so all 118 run
# concurrently in a single wave; job elapsed then tracks the slowest single
# 220-day release (~17.5h observed), comfortably inside the 24h limit. The
# canonical full-year runs live in
# 010_FucusDispersal_{bottom,surface,surface_stokes}_job.sh; this script only
# exists to avoid re-running the ~750 already-complete releases. Each missing
# date has no existing store, so the fresh per-run seed cannot collide.
# Delete this script once the audit reports 73/73 for every regime-year.
#
# Gap set (regime, year, doys):
#   bottom         2018  201
#   surface        2019  186 361
#   surface_stokes 2017  126 356 361
#   surface_stokes 2018  171..361  (the truncated tail, 39 releases)
#   surface_stokes 2019  1..361    (entire year absent, 73 releases)

module load gcc12-env/12.3.0
module load singularity/3.11.5

# Submit from the repo root — relative paths and `-B $PWD:/work` assume it.

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

# Emit one "year doy regime" line per missing release. doys map to dates the
# same way as the per-regime scripts (date -d "<year>-01-01 +<doy-1> days").
emit_gaps() {
    for doy in 201;              do echo "2018 ${doy} bottom";         done
    for doy in 186 361;          do echo "2019 ${doy} surface";        done
    for doy in 126 356 361;      do echo "2017 ${doy} surface_stokes"; done
    for doy in $(seq 171 5 361); do echo "2018 ${doy} surface_stokes"; done
    for doy in $(seq 1   5 361); do echo "2019 ${doy} surface_stokes"; done
}

emit_gaps | while read -r year doy regime; do
    # Never feed an ISO datetime back into `date -d` — its parser silently
    # drops "+N days" after a T-time.
    start_time=$(date -d "${year}-01-01 +$((doy - 1)) days" +%Y-%m-%dT00:15:00)
    end_time=$(date -d "${year}-01-01 +$((doy - 1 + simulation_days)) days" +%Y-%m-%dT00:15:00)

    printf '%s\0' "${start_time} ${end_time} ${regime} $((RANDOM * 32768 + RANDOM))"
done | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _

jobinfo
