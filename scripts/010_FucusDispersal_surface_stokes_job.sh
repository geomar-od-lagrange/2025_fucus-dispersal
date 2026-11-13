#!/bin/bash
#SBATCH --job-name=010_FucusSurfaceStokes
#SBATCH --ntasks=73
#SBATCH --cpus-per-task=1
# Sized from measured peaks, not guessed. The peak is the in-memory output
# store at the final write, so it scales with how well the trajectories
# compress: pinned (motionless) drifters compress well. Without the wall
# velocity the largest MaxRSS was 11.52 GiB for bottom (job 23836392) and
# 7.00 GiB for surface_stokes (job 23836402, 82 % pinned); with the wall
# velocity surface_stokes pins no longer and is sized like bottom. 16G is
# ~1.4x the old bottom peak (bottom pinned 7 % of particle-hours, so its
# output compresses a little worse with the wall velocity too), and fits ~1,900 one-core tasks on the scattered free
# cores of base where 20G fits ~1,100.
#SBATCH --mem-per-cpu=16G
#SBATCH --time=36:00:00
#SBATCH --partition=base
# One dead node must not cancel every other release's step; the failed
# releases come back as missing zarrs and get resubmitted.
#SBATCH --no-kill
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/010/%x_%j.out

module load gcc12-env/12.3.0
module load singularity/3.11.5

# Submit from the repo root — relative paths and `-B $PWD:/work` assume it.

# Positional args: release year, then how many of the year's releases to run.
# `n_releases` exists because the last forced year is partial -- a release only
# runs if its whole simulation_days window sits inside the available forcing.
# `RELEASE_NS` (env, e.g. "1,19") overrides it with an explicit list of release
# indices, because a lost release is a hole in the middle of the schedule, not
# a missing tail, and a prefix count cannot express one.
# Submit via `scripts/submit_010.sh <regime> <year> [n_releases]`, which sizes
# `--ntasks` from whichever of the two selects the work.
# OUTPUT_ROOT and SIMULATION_DAYS (env) override the output root and the run
# length, for test runs; a test root needs 2d_fields/ and wall_fields/ (links
# to the real ones do).
year="${1:-2019}"
n_releases="${2:-73}"

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=${OUTPUT_ROOT:-/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs}
container=parcels-container_2024.10.07-7af7fd0.sif
simulation_days=${SIMULATION_DAYS:-220}
calc_dt_mins=5
output_dt_mins=60
particles_per_cell=100
# Wall velocity (mm/s) on never-flowing faces, plus the land-step reject
# (docs/wall_velocity.md); the walls come from output_root/wall_fields/ (005).
wall_velocity_mm_s=2.0

# Drop any CPU-bind mask inherited by the batch step; otherwise every one of
# the concurrent srun steps below dies at launch with "CPU binding outside of
# job step allocation, allocated CPUs are: 0x...". --exact plus an explicit
# -N1 and per-step cpu/mem lets each step pick its own binding. Same fix as
# scripts/024d_BuildBeachingForcing_job.sh.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

mkdir -p "${output_root}/logs/010" notebooks_executed/TrajectoryCalc/

run_experiment() {
    local start_time="$1"
    local end_time="$2"
    local regime="$3"
    local rng_seed="$4"
    local start_stem="${start_time//[-:]/}"

    srun -n1 -N1 --exact --no-kill \
        --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
        --mem-per-cpu="${SLURM_MEM_PER_CPU}M" \
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
            -p wall_velocity_mm_s ${wall_velocity_mm_s} \
            -p reject_land_steps True \
            -p data_root ${repo_root}/data \
            -p output_root ${output_root} \
            -k python"
}
export -f run_experiment
export repo_root output_root container calc_dt_mins output_dt_mins particles_per_cell wall_velocity_mm_s

# Either the first `n_releases` of the schedule, or exactly the indices in
# RELEASE_NS. submit_010.sh sizes --ntasks from the same selection, so the
# allocation always matches the number of steps the fan-out launches.
if [ -n "${RELEASE_NS:-}" ]; then
    release_ns=${RELEASE_NS//,/ }
else
    release_ns=$(seq 0 $((n_releases - 1)))
fi

echo "year=${year} releases=[${release_ns}] SLURM_NTASKS=${SLURM_NTASKS}"

# Releases at doy = 1 + 5n over the selected n; a full year is 73 of
# them (doys 1..361, leap-year-agnostic). One papermill per start_time. The surface_stokes
# regime layers BSH surface currents with CMEMS Stokes drift (Baltic
# high-res + WAVERYS fallback for the German Bight strip); the
# preprocessed 2D fieldset is produced by 003_prepare_2d_fields.py.
for n in ${release_ns}; do
    doy=$((1 + 5 * n))
    # Derive both times from year+doy; never feed an ISO datetime back into
    # `date -d` — its parser silently drops "+N days" after a T-time.
    start_time=$(date -d "${year}-01-01 +$((doy - 1)) days" +%Y-%m-%dT00:15:00)
    end_time=$(date -d "${year}-01-01 +$((doy - 1 + simulation_days)) days" +%Y-%m-%dT00:15:00)

    printf '%s\0' "${start_time} ${end_time} surface_stokes $((RANDOM * 32768 + RANDOM))"
done | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _

jobinfo
