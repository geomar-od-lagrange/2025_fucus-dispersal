#!/bin/bash
#SBATCH --job-name=010_FucusBottom
#SBATCH --ntasks=73
#SBATCH --cpus-per-task=1
# Sized from measured peaks, not guessed. The peak is the in-memory output
# store at the final write, so it grows with the compressed trajectory size
# (fixed baseline ~1.8 GiB plus ~3.5x the zarr size, fitted to the
# pre-wall-velocity production maxima: bottom 11.5 GiB at 2.8 GB/zarr,
# surface_stokes 7.0 GiB at 1.5 GB/zarr). The wall velocity grows the
# 180-d bottom zarr by 1.05x (tmp_wall_velocity prototype, eps 0 vs 2 mm/s),
# so ~12 GiB expected; 15G is ~1.25x that.
#SBATCH --mem-per-cpu=15G
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
# The optional third arg is 010's `break_off_speed_mm_s`, the bed-friction
# break-off speed U_b (default 0 = no friction). Off 0 the run writes under its
# own regime, `bottom_ub<U_b>` (49 -> Trajectories/bottom_ub49/), while still
# reading the bottom 2D fields and walls. Not passed by hand: the per-regime
# wrapper scripts/010_FucusDispersal_bottom_ub49_job.sh pins it, so the regime
# name in the script name is the one the run writes (docs/bottom_friction.md).
# Off 0 every release also reuses the RNG seed of the `bottom` run of the same
# start time (read off its zarr name under output_root), so both regimes start
# from identical positions: a paired experiment. A release whose bottom zarr
# is missing or ambiguous is not run and fails the job.
# OUTPUT_ROOT and SIMULATION_DAYS (env) override the output root and the run
# length, for test runs; a test root needs 2d_fields/ and wall_fields/ (links
# to the real ones do).
year="${1:-2019}"
n_releases="${2:-73}"
# %g as 010 formats it, so 49.0 and 49 both name `bottom_ub49`.
break_off_speed_mm_s=$(LC_ALL=C printf '%g' "${3:-0}")
# The output regime, as 010 derives it (for the executed-notebook name only).
if [ "${break_off_speed_mm_s}" = 0 ]; then
    output_regime=bottom
else
    output_regime=bottom_ub${break_off_speed_mm_s//./p}
fi

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

    # Kernel start race: jupyter_client reserves the kernel's ZMQ ports by
    # binding port 0 and closing it, so a kernel starting on the same node in
    # that window can take one ("ZMQError: Address already in use", "Kernel died
    # before replying to kernel_info"). 9 of 1,526 releases died this way on
    # 2026-11-13 when all tasks of a node started in the same second. A random
    # start offset makes the race rare; the retry covers the rest. Only failures
    # in the first 5 min are retried: a late failure (OOM, a numerical error)
    # would fail again hours later.
    sleep $(( RANDOM % 30 ))
    local attempt t0 rc=1
    for attempt in 1 2 3; do
        t0=${SECONDS}
        srun -n1 -N1 --exact --no-kill \
            --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
            --mem-per-cpu="${SLURM_MEM_PER_CPU}M" \
            singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work \
            ${container} bash -c \
            ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
            && papermill --cwd notebooks/ \
                notebooks/010_FucusDispersal.ipynb \
                notebooks_executed/TrajectoryCalc/Fucus_${start_stem}_${output_regime}.ipynb \
                -p start_time ${start_time} \
                -p end_time ${end_time} \
                -p regime ${regime} \
                -p break_off_speed_mm_s ${break_off_speed_mm_s} \
                -p RNG_seed ${rng_seed} \
                -p calc_dt_mins ${calc_dt_mins} \
                -p output_dt_mins ${output_dt_mins} \
                -p particles_per_cell ${particles_per_cell} \
                -p wall_velocity_mm_s ${wall_velocity_mm_s} \
                -p reject_land_steps True \
                -p data_root ${repo_root}/data \
                -p output_root ${output_root} \
                -k python" && return 0
        rc=$?
        if (( SECONDS - t0 > 300 )); then
            echo "FAIL ${start_stem} ${regime} rc=${rc} after $(( SECONDS - t0 )) s, not retried" >&2
            return ${rc}
        fi
        echo "attempt ${attempt} failed at start (rc=${rc}) ${start_stem} ${regime}" >&2
        sleep $(( RANDOM % 20 + 5 ))
    done
    return ${rc}
}
export -f run_experiment
export repo_root output_root container calc_dt_mins output_dt_mins particles_per_cell wall_velocity_mm_s break_off_speed_mm_s output_regime

# Either the first `n_releases` of the schedule, or exactly the indices in
# RELEASE_NS. submit_010.sh sizes --ntasks from the same selection, so the
# allocation always matches the number of steps the fan-out launches.
if [ -n "${RELEASE_NS:-}" ]; then
    release_ns=${RELEASE_NS//,/ }
else
    release_ns=$(seq 0 $((n_releases - 1)))
fi

echo "year=${year} releases=[${release_ns}] break_off_speed_mm_s=${break_off_speed_mm_s} output_regime=${output_regime} simulation_days=${simulation_days} wall_velocity_mm_s=${wall_velocity_mm_s} output_root=${output_root} SLURM_NTASKS=${SLURM_NTASKS}"

# The RNG seed of the `bottom` run released at start stem $1
# (YYYYMMDDTHHMMSS), from its zarr name
# Trajectories/bottom/<year>/Fucus_BSH_<stem>_bottom_dt<N>min_seed<S>.zarr.
# Prints the seed; exactly one such zarr must exist, else it says why on
# stderr and returns 1 -- never a fresh seed, which would unpair the release.
paired_seed() {
    local stem=$1 matches seed
    mapfile -t matches < <(compgen -G \
        "${output_root}/Trajectories/bottom/${year}/Fucus_BSH_${stem}_bottom_dt${output_dt_mins}min_seed*.zarr")
    if [ "${#matches[@]}" -ne 1 ]; then
        echo "FAIL ${stem} ${output_regime}: ${#matches[@]} bottom zarrs to pair with" \
             "under ${output_root}/Trajectories/bottom/${year}/, need exactly 1" >&2
        return 1
    fi
    seed=${matches[0]##*_seed}
    seed=${seed%.zarr}
    if ! [[ "${seed}" =~ ^[0-9]+$ ]]; then
        echo "FAIL ${stem} ${output_regime}: no seed in ${matches[0]}" >&2
        return 1
    fi
    echo "${seed}"
}

# Releases at doy = 1 + 5n over the selected n; a full year is 73 of
# them (doys 1..361, leap-year-agnostic). Each start_time runs the bottom 2D
# fields under 010's wall-velocity RK4 with the bottom walls, with bed
# friction at break_off_speed_mm_s (on U and V, not on the wall velocity).
# The work list is built first, so a release whose paired seed is missing is
# counted here and the job can fail on it after the others have run.
work=()
n_unpaired=0
for n in ${release_ns}; do
    doy=$((1 + 5 * n))
    # Derive both times from year+doy; never feed an ISO datetime back into
    # `date -d` — its parser silently drops "+N days" after a T-time.
    start_time=$(date -d "${year}-01-01 +$((doy - 1)) days" +%Y-%m-%dT00:15:00)
    end_time=$(date -d "${year}-01-01 +$((doy - 1 + simulation_days)) days" +%Y-%m-%dT00:15:00)

    if [ "${break_off_speed_mm_s}" = 0 ]; then
        seed=$((RANDOM * 32768 + RANDOM))
    elif ! seed=$(paired_seed "${start_time//[-:]/}"); then
        n_unpaired=$((n_unpaired + 1))
        continue
    fi
    work+=("${start_time} ${end_time} bottom ${seed}")
done
echo "releases to run: ${#work[@]}, not run for want of a paired bottom seed: ${n_unpaired}"
if [ "${#work[@]}" -gt 0 ]; then
    printf '%s\0' "${work[@]}" | xargs -0 -P ${SLURM_NTASKS} -n 1 bash -c 'run_experiment $1' _
fi

jobinfo
# A release left unpaired is a hole that must not pass as a finished job.
[ "${n_unpaired}" -eq 0 ]
