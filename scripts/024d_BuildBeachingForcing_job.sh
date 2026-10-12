#!/bin/bash
#SBATCH --job-name=024d_forcing
#SBATCH --partition=base
# Concurrency ceiling for a direct `sbatch`. scripts/submit_024d.sh sizes the
# allocation from the discovered work list instead, so a short list books a
# short allocation.
#SBATCH --ntasks=292
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=9G
#SBATCH --time=1-00:00:00
#SBATCH --distribution=cyclic
# Do not abort the whole allocation when one node dies: job 23834688 lost
# nesh-clk542 mid-run and Slurm killed all 292 steps (53 done). With --no-kill
# the surviving steps run on; the stems from the dead node come back as FAIL
# and are resubmitted via ONLY_STEMS_FILE.
#SBATCH --no-kill
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/024d/%x_%j.out

# Beaching forcing sidecar: reads one trajectory zarr + the raw baltic_highres
# Stokes + the 024a key, writes one sidecar zarr under
# output_root/BeachingForcing/<regime>/<year>/.
#
# Unit of work = one trajectory zarr, not a (year, month) cell. At
# window_days=220 (the full Parcels simulation, obs=5280) a cell's ~6 zarrs run
# strictly sequentially inside one papermill process, so a month grid caps
# concurrency at 12 per release year and makes wall time 6x the per-zarr cost.
# One srun step per zarr is one independent single-process numpy run per zarr
# (no MPI, no Dask) whose
# bottleneck is streaming trajectory + hourly Stokes off GPFS, so they want
# horizontal filesystem reach rather than node locality -- but --spread-job is
# the wrong lever at this task count: it asks for one node per task, and a
# 292-node request against ~94 idle nodes pended on (Resources) indefinitely
# (job 23834728). --distribution=cyclic alone round-robins the tasks over
# however many nodes the scheduler grants, which is the spread that matters.
# No --constraint either: it only shrinks the eligible pool.
#
# Memory: 9G x 2 cpu = 18 GB/task. Measured peak at window_days=220 is
# MaxRSS 11556164K = 11.02 GiB (max over the 292 steps of job 23834688), and
# 18 GB is that plus 55% headroom. Steps there ran at 12 GB/task with no
# oom-kill, so the margin is thin rather than absent.
#
# njobs is set by the zarr count on disk, concurrency is whatever --ntasks the
# scheduler grants (`xargs -P ${SLURM_NTASKS}`), so the two rescale
# independently: resubmit with --ntasks=146 under load without editing this.
# Tasks share no output, and the notebook skips sidecars that already carry the
# current sampling parameters, so a partially failed job is just resubmitted.
#
# Submit through scripts/submit_024d.sh: it expands scripts/node-blacklist.txt
# into `sbatch --exclude=` (an #SBATCH directive cannot read a file) and passes
# extra args through. The exclusion applies to the allocation only: a step-level
# `srun --exclude` inside the allocation made SLURM serialise the steps
# ("step creation temporarily disabled ... Requested nodes are busy", job
# 23834727, 1 of 292 steps running), so it is deliberately not repeated here.
#
# Usage: scripts/submit_024d.sh
#   scripts/submit_024d.sh --ntasks=146                            # throttle
#   YEARS="2020 2021" scripts/submit_024d.sh                       # subset
#   ONLY_STEMS_FILE=/path/stems.txt scripts/submit_024d.sh --ntasks=8  # retry
#     (that file holds one "<year> <stem>" line per zarr to rebuild)
#   LIST_ONLY=1 bash scripts/024d_BuildBeachingForcing_job.sh      # print work
# 024a_BuildHexKey_job.sh must have run first for the matching hex_radius.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
regime=surface_stokes
hex_radius=6000
logdir=${output_root}/logs/024d
mkdir -p "${logdir}/executed"

# Drop any CPU-bind mask inherited from an outer allocation; otherwise the
# concurrent srun steps below fail with "CPU binding outside of job step
# allocation". --exact sets each step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

# Blacklisted nodes (see scripts/node-blacklist.txt), for the provenance line.
exclude=$(sed "s/#.*//" scripts/node-blacklist.txt | tr -d "[:blank:]" \
          | grep -v "^$" | paste -sd,)

export output_root regime hex_radius logdir exclude

run_one() {
    year=$1
    stem=$2
    t0=${SECONDS}
    # Per-task Jupyter runtime dir. jupyter_client reserves five ZMQ ports by
    # binding to port 0 and closing them; the kernel re-binds moments later and
    # a concurrent kernel on the same host can steal one in that window. The
    # private dir keeps connection files off GPFS; the 3-attempt retry below is
    # what actually covers the race (port names carry UUIDs).
    rt="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/jupyter-runtime-${year}-${stem}"
    mkdir -p "${rt}"
    rc=1
    for attempt in 1 2 3; do
        srun -n1 -N1 --exact --no-kill \
            --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
            --mem-per-cpu="${SLURM_MEM_PER_CPU}M" \
            --job-name="024d_${stem}" \
            env OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                JUPYTER_RUNTIME_DIR="${rt}" \
            pixi run papermill --cwd notebooks/ \
                notebooks/024d_BuildBeachingForcing.ipynb \
                "${logdir}/executed/024d_${year}_${stem}.ipynb" \
                -p output_root "${output_root}" \
                -p regime "${regime}" \
                -p release_year "${year}" \
                -p hex_radius "${hex_radius}" \
                -p overwrite True \
                -r only_stem "${stem}" \
                -k python && { rc=0; break; }
        rc=$?
        echo "attempt ${attempt} failed (rc=${rc}) ${year} ${stem}" >&2
        sleep $(( RANDOM % 10 + 1 ))
    done
    if [ ${rc} -eq 0 ]; then
        echo "OK ${year} ${stem} $(( SECONDS - t0 ))"
    else
        echo "FAIL ${year} ${stem}"
    fi
    return 0
}
export -f run_one

# "<year> <stem>" per line: every trajectory zarr on disk for this regime, or
# the retry subset. The year set is not a literal -- it is whatever 010 has
# written, so a release year joins the work list as soon as its runs land.
# $YEARS (space-separated) restricts it; an unset $YEARS takes all of them.
list_work() {
    if [ -n "${ONLY_STEMS_FILE:-}" ]; then
        grep -v '^[[:space:]]*$' "${ONLY_STEMS_FILE}"
        return
    fi
    for p in "${output_root}/Trajectories/${regime}"/*/*.zarr; do
        [ -d "${p}" ] || continue
        year=$(basename "$(dirname "${p}")")
        [ -z "${YEARS:-}" ] || [[ " ${YEARS} " == *" ${year} "* ]] || continue
        b=$(basename "${p}")
        echo "${year} ${b%.zarr}"
    done
}

# Work-list preview, for sizing the allocation (scripts/submit_024d.sh) or
# eyeballing what a $YEARS / $ONLY_STEMS_FILE selection resolves to.
if [ -n "${LIST_ONLY:-}" ]; then
    list_work
    exit 0
fi

njobs=$(list_work | wc -l)
# Provenance header: the sidecar zarrs carry git_sha + slurm_job_id as attrs,
# and this pins what that job was, in the .out next to the OK/FAIL tally.
echo "024d provenance: git $(git rev-parse HEAD)"
if [ -n "$(git status --porcelain)" ]; then
    echo "024d WARNING: dirty working tree, this run maps to no commit"
fi
echo "024d slurm: job ${SLURM_JOB_ID}, ntasks ${SLURM_NTASKS}, cpus-per-task"\
     "${SLURM_CPUS_PER_TASK}, mem-per-cpu ${SLURM_MEM_PER_CPU}M, exclude ${exclude:-<none>}"

echo "024d: ${njobs} zarrs, ${SLURM_NTASKS} concurrent, regime=${regime}," \
     "hex_radius=${hex_radius}, logs -> ${logdir}"

# `set -e` must not kill the driver when a task exhausts its retries; run_one
# already swallows its own status and reports FAIL, and the tally below is the
# authoritative exit code.
set +e
list_work | xargs -P "${SLURM_NTASKS}" -L1 bash -c 'run_one "$@"' _ \
    | tee "${logdir}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.tally"
set -e

tally=${logdir}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.tally
nok=$(grep -c '^OK ' "${tally}" || true)
nfail=$(grep -c '^FAIL ' "${tally}" || true)
echo "024d done: ${nok} OK, ${nfail} FAIL of ${njobs}"
grep '^FAIL ' "${tally}" | sed 's/^FAIL //' || true

jobinfo || true
[ "${nfail}" -eq 0 ]
