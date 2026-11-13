#!/bin/bash
#SBATCH --job-name=005_BuildWallFields
#SBATCH --partition=base
# Many one-core tasks: the work is ~2.2 TB of c_file reads per regime, GPFS-
# bound and embarrassingly parallel, and one-core tasks backfill onto the
# scattered free cores of busy nodes. Each worker holds one c_file (<= 0.1 GB)
# plus two boolean face masks.
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH --time=04:00:00
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/005/%x_%j.out

# Wall-velocity and land fields for one regime (notebooks/005_BuildWallFields),
# written to output_root/wall_fields/walls_{fine,coarse}_<regime>.nc. Needs
# every 2d_fields c_file of the regime (003) and is needed by 010.
#
#   sbatch scripts/005_BuildWallFields_job.sh surface_stokes
#   sbatch scripts/005_BuildWallFields_job.sh bottom      # also used by bottom_vs0p8
#
# OUTPUT_ROOT (env) overrides the output root, for tests on a root whose
# 2d_fields/ holds a subset of the c_files.
#
# Multi-task dask layout, one task per SLURM task (as 024):
#   task 0: dask scheduler + one worker
#   task 1: papermill (connects to scheduler via $SCHEDULER_FILE)
#   tasks 2..N-1: extra dask workers

export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128
export no_proxy=localhost,127.0.0.1,0.0.0.0,10.0.0.0/8

regime="${1:-surface_stokes}"

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=${OUTPUT_ROOT:-/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs}

export SCHEDULER_FILE=${repo_root}/.scheduler_${SLURM_JOB_ID}.json

mkdir -p notebooks_executed/Inputs/

cleanup() {
    rm -f "${SCHEDULER_FILE}"
    kill $(jobs -p) 2>/dev/null || true
}
trap cleanup EXIT

# Same step hygiene as scripts/024_BuildHexAggregates_job.sh: drop the
# inherited CPU-bind mask, one node per step, own binding and memory.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

SRUN_STEP="srun -n1 -N1 --exact --cpus-per-task=${SLURM_CPUS_PER_TASK}"
SRUN_STEP="${SRUN_STEP} --mem-per-cpu=${SLURM_MEM_PER_CPU}M"
WORKER_ARGS="--scheduler-file ${SCHEDULER_FILE}"

${SRUN_STEP} pixi run bash -c "
    dask scheduler --scheduler-file ${SCHEDULER_FILE} &
    dask worker ${WORKER_ARGS} &
    wait
" &

# Wait for the scheduler file, bounded (see 024's job script for why).
for _ in $(seq 60); do
    [ -s "${SCHEDULER_FILE}" ] && break
    sleep 5
done
if [ ! -s "${SCHEDULER_FILE}" ]; then
    echo "dask scheduler did not come up within 300 s -- giving up" >&2
    exit 1
fi

${SRUN_STEP} pixi run papermill --cwd notebooks/ \
    notebooks/005_BuildWallFields.ipynb \
    notebooks_executed/Inputs/005_BuildWallFields_${regime}_${SLURM_JOB_ID}.ipynb \
    -p data_root ${repo_root}/data \
    -p output_root ${output_root} \
    -p regime ${regime} \
    -k python &
PAPERMILL_PID=$!

for i in $(seq 3 ${SLURM_NTASKS}); do
    ${SRUN_STEP} pixi run dask worker ${WORKER_ARGS} &
done

wait ${PAPERMILL_PID}
PAPERMILL_EXIT=$?

jobinfo
exit ${PAPERMILL_EXIT}
