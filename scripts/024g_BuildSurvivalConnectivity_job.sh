#!/bin/bash
#SBATCH --job-name=024g_BuildSurvivalConnectivity
# One cell per (year, month) = |YEARS| x 12 = 48 cells; each cell is a
# seconds-per-sidecar numpy pass over the BeachingForcing sidecar, so a small
# allocation drains the grid in minutes. Raise --ntasks on the command line
# when sweeping several members in one job.
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base
# No --spread-job and no --constraint: like 024f these cells read only the
# compact BeachingForcing sidecar and the key parquet, so they are CPU-bound
# rather than GPFS-bound and packing them tight costs nothing.

set -euo pipefail

# Survival-weighted connectivity reducer: reads the BeachingForcing sidecar
# written by 024d_BuildBeachingForcing plus the 024a key, writes one survconn
# parquet per (regime, release_year, release_month, rate-model member) with
# (origin_subbasin, target_subbasin, release_doy, age_bin) -> (n_obs, w_obs).
# Same rate model as 024e/024f; w_obs is residence weighted by the surviving
# (un-beached) fraction S = exp(-A). 028_SubbasinConnectivityMatrix reads it
# when given a non-empty `member`.
#
# Usage: sbatch scripts/024g_BuildSurvivalConnectivity_job.sh [regime] [hex_radius] [stagger_max_s]
#
# Rate-model parameters come from the environment so a sweep is a loop of
# submissions that differ only in an export:
#
#   W_C         onshore-Stokes threshold, m/s
#   TAU_STRONG_HOURS  e-folding time above W_C
#   TAU_CALM_DAYS    background in-band timescale, d; 0 = off
#   DELTA       linear edge width around W_C, m/s; 0 = hard step
#   BAND_M, CONNECTIVITY_MAX_DAYS  the connectivity horizon must cover the
#                                  largest matrix horizon 028 asks for
#
# 024a_BuildHexKey_job.sh and 024d_BuildBeachingForcing_job.sh must have run
# first for the matching hex_radius / (regime, year).

YEARS=(2016 2017 2018 2019)

# Drop any CPU-bind mask inherited from an outer allocation (present when
# submitted from inside an interactive job); --exact sets each step's own.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"

W_C="${W_C:-0.1}"
TAU_STRONG_HOURS="${TAU_STRONG_HOURS:-3.0}"
TAU_CALM_DAYS="${TAU_CALM_DAYS:-0.0}"
DELTA="${DELTA:-0.0}"
BAND_M="${BAND_M:-2000.0}"
CONNECTIVITY_MAX_DAYS="${CONNECTIVITY_MAX_DAYS:-120}"

# Executed-notebook names carry the same member tag the notebook builds for the
# store filename, so a notebook and the parquet it wrote are matched by eye.
# awk's %g is printf's %g is Python's :g, so the two constructions agree.
fmt_g() { awk -v x="$1" 'BEGIN { printf "%g", x }'; }
tc=$(awk -v x="${TAU_CALM_DAYS}" 'BEGIN { if (x > 0) printf "%g", x; else printf "inf" }')
member="step_wc$(fmt_g "${W_C}")_ts$(fmt_g "${TAU_STRONG_HOURS}")_tc${tc}"
if awk -v x="${DELTA}" 'BEGIN { exit !(x > 0) }'; then
    member="${member}_d$(fmt_g "${DELTA}")"
fi
member="${member//./p}"
echo "member: ${member}"

# Max random start delay per cell (s), 0.1 s granularity; see the 024d job
# script — de-synchronises the Jupyter kernel start-up race.
stagger_max_s="${3:-15}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export output_root regime hex_radius stagger_max_s member
export W_C TAU_STRONG_HOURS TAU_CALM_DAYS DELTA
export BAND_M CONNECTIVITY_MAX_DAYS

mkdir -p notebooks_executed/Visualisations/

# `set -e` would abort on a failing pipeline before `rc` is read, so the
# exit status is captured in the `||` branch instead.
rc=0
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
    # Retry kernel start-up. jupyter_client picks five free TCP ports by
    # binding to port 0 and CLOSING the socket, then the kernel re-binds them
    # later — a TOCTOU window in which a concurrent kernel on the same host
    # can steal a port. The loser dies with ZMQ "Address already in use"
    # before running any cell. The window is short and ports are re-drawn on
    # each attempt, so a couple of retries removes the failure mode.
    for attempt in 1 2 3; do
        srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --exact \
            pixi run papermill --cwd notebooks/ \
            notebooks/024g_BuildSurvivalConnectivity.ipynb \
            notebooks_executed/Visualisations/024g_BuildSurvivalConnectivity_${regime}_${year}${ms}_${member}_r${hex_radius}m.ipynb \
            -p output_root ${output_root} \
            -p regime ${regime} \
            -p release_year ${year} \
            -p release_month ${month} \
            -p hex_radius ${hex_radius} \
            -p w_c ${W_C} \
            -p tau_strong_hours ${TAU_STRONG_HOURS} \
            -p tau_calm_days ${TAU_CALM_DAYS} \
            -p delta ${DELTA} \
            -p band_m ${BAND_M} \
            -p connectivity_max_days ${CONNECTIVITY_MAX_DAYS} \
            -k python
        rc=$?
        [ ${rc} -eq 0 ] && break
        echo "cell attempt ${attempt} failed (rc=${rc}); retrying" >&2
        sleep $(( RANDOM % 10 + 1 ))
    done
    # Propagate the final status: without this the loop ends on `sleep`,
    # bash -c exits 0, and xargs reports success for an exhausted cell.
    exit ${rc}
' _ || rc=$?

jobinfo
exit ${rc}
