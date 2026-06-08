#!/bin/bash
#SBATCH --job-name=021_TimeStats
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --partition=base

# Multi-task dask layout, one task per SLURM task, adjust --ntasks at submit:
#   task 0: dask scheduler + one local worker process set
#   task 1: papermill (connects to scheduler via $SCHEDULER_FILE)
#   tasks 2..N-1: extra dask worker tasks attached to the scheduler
# Dask inter-node traffic is raw TCP; HTTP_PROXY/HTTPS_PROXY are ignored
# by Tornado streams, but we set no_proxy defensively for anything HTTP
# inside the env that might otherwise dial the proxy for node addresses.

export http_proxy=http://10.0.7.235:3128
export https_proxy=http://10.0.7.235:3128
export no_proxy=localhost,127.0.0.1,0.0.0.0,10.0.0.0/8

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs

export SCHEDULER_FILE=${repo_root}/.scheduler_${SLURM_JOB_ID}.json

mkdir -p notebooks_executed/Visualisations/

cleanup() {
    rm -f "${SCHEDULER_FILE}"
    kill $(jobs -p) 2>/dev/null || true
}
trap cleanup EXIT

SRUN_STEP="srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --exact"
WORKER_ARGS="--scheduler-file ${SCHEDULER_FILE} --interface ib0"

# Task 0: scheduler + attached workers (shape auto-detected from cgroup).
# dask worker --scheduler-file polls the file until it appears, so no
# explicit wait loop needed here.
${SRUN_STEP} pixi run bash -c "
    dask scheduler --interface ib0 --scheduler-file ${SCHEDULER_FILE} &
    dask worker ${WORKER_ARGS} &
    wait
" &

sleep 30

# Task 1: papermill.
${SRUN_STEP} pixi run papermill --cwd notebooks/ \
    notebooks/021_TimeStats.ipynb \
    notebooks_executed/Visualisations/021_TimeStats.ipynb \
    -p output_root ${output_root} \
    -k python &
PAPERMILL_PID=$!

# Tasks 2..N-1: additional worker tasks.
for i in $(seq 3 ${SLURM_NTASKS}); do
    ${SRUN_STEP} pixi run dask worker ${WORKER_ARGS} &
done

wait ${PAPERMILL_PID}
PAPERMILL_EXIT=$?

jobinfo
exit ${PAPERMILL_EXIT}
