#!/bin/bash
#SBATCH --job-name=024e_BuildSurvivalOccupancy
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12G
#SBATCH --time=04:00:00
#SBATCH --partition=base
#SBATCH --constraint=sapphire

# Survival-weighted occupancy pass: reads the trajectory zarrs + raw
# baltic_highres Stokes + the 024a key, writes one survocc parquet per
# (regime, release_year, release_month). Same store-then-consume split and
# rate model as 024d, but accumulates residence weighted by the surviving
# (un-beached) fraction S = exp(-A) over a longer horizon (occupancy_max_days,
# notebook default 120 d). Heavier than 024d (~2-3x the obs window + Stokes
# hours), so ~15 min per month cell.
#
# Parallelism fans out the full (year × month) grid with xargs dispatching
# one `srun --ntasks=1 --cpus-per-task --exact` step per cell; njobs is fixed
# (|YEARS|·12 = 48) and concurrency is whatever --ntasks the scheduler grants.
# Pinned to sapphire (srp) nodes via --constraint for reliable IB networking.
#
# Usage: sbatch scripts/024e_BuildSurvivalOccupancy_job.sh [regime] [hex_radius]
# 024a_BuildHexKey_job.sh must have run first for the matching hex_radius.

YEARS=(2016 2017 2018 2019)

# Drop any CPU-bind mask inherited from an outer allocation (present when
# submitted from inside an interactive job); --exact sets each step's own.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
w_half="${3:-0.05}"
# Max random start delay per cell (s), 0.1 s granularity; see 024d job
# script — de-synchronises the Jupyter kernel start-up race.
stagger_max_s="${4:-30}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export output_root regime hex_radius w_half stagger_max_s

mkdir -p notebooks_executed/Visualisations/

for year in "${YEARS[@]}"; do
    for month in $(seq 1 12); do
        printf '%s\0' "${year} ${month}"
    done
done | xargs -0 -P "${SLURM_NTASKS}" -n 1 bash -c '
    read -r year month <<< "$1"
    tenths=$(( RANDOM % (stagger_max_s * 10 + 1) ))
    sleep "$(( tenths / 10 )).$(( tenths % 10 ))"
    # Per-cell Jupyter runtime dir; see 024d job script — shared connection
    # files collide at high concurrency and the kernel never starts.
    export JUPYTER_RUNTIME_DIR="${SLURM_TMPDIR:-/tmp}/jupyter-runtime-$$"
    mkdir -p "${JUPYTER_RUNTIME_DIR}"
    ms=$(printf "_m%02d" "${month}")
    whtag="_wh${w_half//./p}"
    # Retry kernel start-up. jupyter_client picks five free TCP ports by
    # binding to port 0 and CLOSING the socket, then the kernel re-binds them
    # later — a TOCTOU window in which a concurrent kernel on the same host
    # can steal a port. The loser dies with ZMQ "Address already in use"
    # before running any cell. The window is short and ports are re-drawn on
    # each attempt, so a couple of retries removes the failure mode.
    for attempt in 1 2 3; do
        srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --exact \
            pixi run papermill --cwd notebooks/ \
            notebooks/024e_BuildSurvivalOccupancy.ipynb \
            notebooks_executed/Visualisations/024e_BuildSurvivalOccupancy_${regime}_${year}${ms}${whtag}_r${hex_radius}m.ipynb \
            -p output_root ${output_root} \
            -p regime ${regime} \
            -p release_year ${year} \
            -p release_month ${month} \
            -p hex_radius ${hex_radius} \
            -p w_half ${w_half} \
            -k python
        rc=$?
        [ ${rc} -eq 0 ] && break
        echo "cell attempt ${attempt} failed (rc=${rc}); retrying" >&2
        sleep $(( RANDOM % 10 + 1 ))
    done
    # Propagate the final status: without this the loop ends on `sleep`,
    # bash -c exits 0, and xargs reports success for an exhausted cell.
    exit ${rc}
' _
rc=$?

jobinfo
exit ${rc}
