#!/bin/bash
#SBATCH --job-name=024d_BuildBeachingForcing
# The grid is |YEARS| x 12 = 48 cells. Asking for more tasks than cells just
# inflates the allocation (idle slots) and makes the job harder to schedule.
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=2
# 8G x 2 CPU = 16 GB/task. The builder holds a handful of (n_real, nobs)
# arrays plus the raster; measured peak is ~10.6 GB, so this is ~50% headroom.
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=base
# Spread tasks over as many nodes as possible. These cells are independent
# single-process papermill runs whose bottleneck is streaming trajectory
# zarrs and hourly Stokes files off GPFS, so what matters is horizontal
# reach into the distributed filesystem, not node locality. Packing them
# tight starves the tail: an 8-node allocation of the same 48 cells ran
# steps min/med/max 11:54/17:55/54:51 against 10:51/16:43/20:54 on 22
# nodes -- same median, 2.6x worse tail. --spread-job disables the
# topology/tree plugin, which costs nothing here since there is no MPI.
#SBATCH --spread-job
# No --constraint: these cells are embarrassingly parallel single-process
# papermill runs with no MPI and no Dask cluster, so the IB-reliability
# rationale for pinning to sapphire (srp) does not apply to them. Leaving
# the whole base partition eligible cuts queue time substantially.

# Beaching forcing sidecar: reads the trajectory zarrs + raw baltic_highres
# Stokes + the 024a key, writes one sidecar zarr per trajectory zarr under
# output_root/BeachingForcing/<regime>/<year>/. A single-process numpy
# notebook (no Dask) — each zarr fits in memory and the bottleneck is Stokes
# I/O, so one papermill run processes a month's ~6 zarrs sequentially
# (~95 s each at window_days=120).
#
# There are no rate-model parameters here: the sidecar caches the
# parameter-free ingredients (onshore Stokes, distance to land, hex id,
# shore type, displacement) and the reducers sweep the rate model over them.
#
# Parallelism fans out the (year x month) grid with xargs dispatching one
# `srun --ntasks=1 -c ${SLURM_CPUS_PER_TASK} --exact` step per cell, like
# 010_FucusDispersal_*. njobs is fixed by the grid but concurrency is
# whatever --ntasks the scheduler grants (`xargs -P ${SLURM_NTASKS}`), so the
# two rescale independently: request fewer tasks under load without editing
# the job. Cells share no output — one zarr per source zarr — and the
# notebook skips sidecars that already carry the current sampling
# parameters, so a partially failed job is simply resubmitted.
#
# Usage: sbatch scripts/024d_BuildBeachingForcing_job.sh [regime] [hex_radius] [stagger_max_s]
#   sbatch scripts/024d_BuildBeachingForcing_job.sh                   # surface_stokes, 48 cells
#   sbatch scripts/024d_BuildBeachingForcing_job.sh surface 6000       # sensitivity regime
#   sbatch --ntasks=8 scripts/024d_BuildBeachingForcing_job.sh         # throttle concurrency
# 024a_BuildHexKey_job.sh must have run first for the matching hex_radius.

YEARS=(2016 2017 2018 2019)

# Drop any CPU-bind mask inherited from an outer allocation (present when
# this is submitted from inside an interactive job); otherwise the
# concurrent srun steps below try to bind to the outer job's CPUs and fail
# with "CPU binding outside of job step allocation". --exact sets each
# step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"
# Max random start delay per cell, in seconds, drawn at 0.1 s granularity.
# This is a de-synchroniser for the Jupyter kernel start-up race (ZMQ
# "Address already in use" when many kernels claim connection files/ports
# at once — it cost us a cell at 48-way and again at 100-way). The race
# window is milliseconds, so tenths of a second are the right unit and a
# few seconds of total spread is plenty; this is NOT bandwidth throttling.
stagger_max_s="${3:-30}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export output_root regime hex_radius stagger_max_s

mkdir -p notebooks_executed/Visualisations/

echo "grid: ${#YEARS[@]} years x 12 months = $((${#YEARS[@]} * 12)) cells," \
     "${SLURM_NTASKS} concurrent"
echo "regime: ${regime}   hex_radius: ${hex_radius}"
echo "stagger: random 0..${stagger_max_s}s per cell, 0.1s granularity"

# Emit every "year month" pair; xargs runs up to ${SLURM_NTASKS} at once,
# each dispatching one srun --ntasks=1 -c N --exact papermill step.
for year in "${YEARS[@]}"; do
    for month in $(seq 1 12); do
        printf '%s\0' "${year} ${month}"
    done
done | xargs -0 -P "${SLURM_NTASKS}" -n 1 bash -c '
    read -r year month <<< "$1"
    # Bash seeds RANDOM per process, so each xargs child gets its own draw.
    tenths=$(( RANDOM % (stagger_max_s * 10 + 1) ))
    sleep "$(( tenths / 10 )).$(( tenths % 10 ))"
    # Per-cell Jupyter runtime dir. Kernels write connection files into a
    # shared runtime dir by default; at high concurrency they collide and the
    # kernel never starts (papermill then writes an output notebook with zero
    # executed cells and no exception). Arrival-time stagger alone did not fix
    # this — the failure rate tracked concurrency, not arrival density — so
    # give every cell its own directory on node-local scratch.
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
            notebooks/024d_BuildBeachingForcing.ipynb \
            notebooks_executed/Visualisations/024d_BuildBeachingForcing_${regime}_${year}${ms}_r${hex_radius}m.ipynb \
            -p output_root ${output_root} \
            -p regime ${regime} \
            -p release_year ${year} \
            -p release_month ${month} \
            -p hex_radius ${hex_radius} \
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
