#!/bin/bash
#SBATCH --job-name=024f_BuildSurvivalOccupancy
# One cell per (release year, month) that has a BeachingForcing sidecar on
# disk -- 12 per full release year, fewer for a partial one. Each cell is a
# seconds-per-sidecar numpy pass, so a small allocation drains the grid in
# minutes. Raise --ntasks on the command line when sweeping several members in
# one job.
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00
#SBATCH --partition=base
# No --spread-job and no --constraint: unlike 024d, these cells read only the
# compact BeachingForcing sidecar and the key parquet, so they are CPU-bound
# rather than GPFS-bound and packing them tight costs nothing.

# Survival-weighted occupancy reducer: reads the BeachingForcing sidecar
# written by 024d_BuildBeachingForcing plus the 024a key, writes one survocc
# parquet per (regime, release_year, release_month, rate-model member). Same
# rate model as 024e_BuildBeaching, but accumulates residence weighted by the
# surviving (un-beached) fraction S = exp(-A) over a longer horizon
# (OCCUPANCY_MAX_DAYS, notebook default 120 d).
#
# Usage: sbatch scripts/024f_BuildSurvivalOccupancy_job.sh [regime] [hex_radius] [stagger_max_s]
#
# Rate-model parameters come from the environment so a sweep is a loop of
# submissions that differ only in an export:
#
#   W_C         onshore-Stokes threshold, m/s
#   TAU_STRONG_HOURS  e-folding time above W_C
#   TAU_CALM_DAYS    background in-band timescale, d; 0 = off
#   DELTA       linear edge width around W_C, m/s; 0 = hard step
#   TRAP_FLAT, TRAP_WALL  strong-wave rate factor on fully flat /
#               fully walled shore; the realised factor is
#               TRAP_WALL + (TRAP_FLAT - TRAP_WALL)*ff/100
#   BAND_M, OCCUPANCY_MAX_DAYS  the occupancy horizon is this reducer's
#                               analogue of 024e's MAX_FLOAT_DAYS
#
# 024a_BuildHexKey_job.sh and 024d_BuildBeachingForcing_job.sh must have run
# first for the matching hex_radius / (regime, year); the (year, month) grid is
# read off the sidecars they wrote, and $YEARS restricts it to a subset.

# Drop any CPU-bind mask inherited from an outer allocation (present when
# submitted from inside an interactive job); --exact sets each step's own.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

# "<year> <month>" per cell: the distinct release months that have a
# BeachingForcing sidecar on disk for this regime. The year set is not a
# literal -- it is whatever 024d has written, so a release year joins the grid
# as soon as its sidecars land, and a partial year contributes only the months
# it actually has. $YEARS (space-separated) restricts it; unset takes all.
list_cells() {
    for p in "${output_root}/BeachingForcing/${regime}"/*/*.zarr; do
        [ -d "${p}" ] || continue
        year=$(basename "$(dirname "${p}")")
        [ -z "${YEARS:-}" ] || [[ " ${YEARS} " == *" ${year} "* ]] || continue
        stamp=$(basename "${p}")
        stamp=${stamp#*_BSH_}
        month=${stamp:4:2}
        echo "${year} ${month#0}"
    done | sort -u -k1,1n -k2,2n
}

# Work-list preview, for eyeballing what a $YEARS selection resolves to
# without submitting: `LIST_ONLY=1 bash <this script> [regime] [hex_radius]`.
if [ -n "${LIST_ONLY:-}" ]; then
    list_cells
    exit 0
fi

W_C="${W_C:-0.1}"
TAU_STRONG_HOURS="${TAU_STRONG_HOURS:-3.0}"
TAU_CALM_DAYS="${TAU_CALM_DAYS:-0.0}"
DELTA="${DELTA:-0.0}"
BAND_M="${BAND_M:-2000.0}"
TRAP_FLAT="${TRAP_FLAT:-1.0}"
TRAP_WALL="${TRAP_WALL:-0.0}"
OCCUPANCY_MAX_DAYS="${OCCUPANCY_MAX_DAYS:-120}"

# Executed-notebook names carry the same member tag the notebook builds for the
# store filename, so a notebook and the parquet it wrote are matched by eye.
# awk's %g is printf's %g is Python's :g, so the two constructions agree.
fmt_g() { awk -v x="$1" 'BEGIN { printf "%g", x }'; }
tc=$(awk -v x="${TAU_CALM_DAYS}" 'BEGIN { if (x > 0) printf "%g", x; else printf "inf" }')
member="step_wc$(fmt_g "${W_C}")_ts$(fmt_g "${TAU_STRONG_HOURS}")_tc${tc}"
if awk -v x="${DELTA}" 'BEGIN { exit !(x > 0) }'; then
    member="${member}_d$(fmt_g "${DELTA}")"
fi
if awk -v a="${TRAP_FLAT}" -v b="${TRAP_WALL}" 'BEGIN { exit !(a != 1 || b != 1) }'; then
    member="${member}_tf$(fmt_g "${TRAP_FLAT}")_tw$(fmt_g "${TRAP_WALL}")"
fi
member="${member//./p}"
echo "member: ${member}"

# Max random start delay per cell (s), 0.1 s granularity; see the 024d job
# script — de-synchronises the Jupyter kernel start-up race.
stagger_max_s="${3:-15}"

export output_root regime hex_radius stagger_max_s member
export W_C TAU_STRONG_HOURS TAU_CALM_DAYS DELTA
export TRAP_FLAT TRAP_WALL
export BAND_M OCCUPANCY_MAX_DAYS

mkdir -p notebooks_executed/Visualisations/

list_cells | tr '\n' '\0' | xargs -0 -P "${SLURM_NTASKS}" -n 1 bash -c '
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
            notebooks/024f_BuildSurvivalOccupancy.ipynb \
            notebooks_executed/Visualisations/024f_BuildSurvivalOccupancy_${regime}_${year}${ms}_${member}_r${hex_radius}m.ipynb \
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
            -p trap_flat ${TRAP_FLAT} \
            -p trap_wall ${TRAP_WALL} \
            -p occupancy_max_days ${OCCUPANCY_MAX_DAYS} \
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
