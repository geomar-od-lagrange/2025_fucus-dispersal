#!/bin/bash
#SBATCH --job-name=024efg_reduce
#SBATCH --partition=base
#SBATCH --ntasks=292
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --distribution=cyclic
# Do not abort the whole allocation when one node dies (see 024d job script):
# the surviving steps run on and the lost cells come back as FAIL.
#SBATCH --no-kill
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/024efg/%x_%j.out

# Re-reduce every beaching store from the rebuilt 024d sidecars, in one job.
#
# Unit of work = one (reducer, member, year, month) cell. Each cell is a
# seconds-per-zarr numpy pass over ~6 sidecar zarrs plus the 024a key, so a
# cell is 1-2 min and the 1200-cell grid drains in a handful of rounds at
# --ntasks=292. Concurrency is whatever --ntasks the scheduler grants
# (`xargs -P ${SLURM_NTASKS}`), independent of the work-list length, so the
# same script can be resubmitted throttled.
#
# Memory: 6G x 2 cpu = 12 GB/task. The single-member 024e/024f/024g jobs ran
# at 8G x 2 = 16 GB with room to spare; these reducers touch only the compact
# sidecar (~30 MB/zarr) and the key parquet.
#
# Submit through scripts/submit_024efg.sh: it expands
# scripts/node-blacklist.txt into `sbatch --exclude=` (an #SBATCH directive
# cannot read a file) and passes extra args through, e.g.
#
#   scripts/submit_024efg.sh --ntasks=292 --dependency=afterok:<024d job id>
#
# The allocation-level exclusion is deliberately NOT repeated as a step-level
# `srun --exclude`: that made SLURM serialise the steps ("step creation
# temporarily disabled ... Requested nodes are busy").
#
# 024a_BuildHexKey_job.sh and 024d_BuildBeachingForcing_job.sh must have run
# first for the matching hex_radius / (regime, year) — the reducers refuse a
# sidecar whose sampling parameters (incl. the shoreclass hash) are stale.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
regime=surface_stokes
hex_radius=6000
band_m=2000.0
max_float_days=60
occupancy_max_days=120
connectivity_max_days=120
YEARS=(2016 2017 2018 2019)
MONTHS=($(seq 1 12))
logdir=${output_root}/logs/024efg
mkdir -p "${logdir}/executed"

# Drop any CPU-bind mask inherited from an outer allocation; otherwise the
# concurrent srun steps below fail with "CPU binding outside of job step
# allocation". --exact sets each step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

# --- the rate-model members, spelled out ---------------------------------
#
# One line per (reducer, member): "reducer w_c tau_strong_h tau_calm_d delta
# trap_flat trap_wall". tau_calm 0 means "no background rate" and tags as
# `tcinf` — that is how the notebooks spell an infinite calm timescale, there
# is no inf sentinel. The member tag is rebuilt below exactly as the notebooks
# build it, so an executed notebook and the parquet it wrote match by eye.
#
# 024e: the 11 step members produced so far plus 3 shore-type (trap) members
# at the production point (step, w_c 0.10, tau_strong 3 h, tau_calm inf).
# The pre-sidecar `sat_t480_wt0p05` member is NOT here: the saturating rate
# form was removed from the notebooks (see plans/done/beaching_review_response.md),
# so it can no longer be produced and papermill would silently write a step
# member instead.
# 024f: the 4 producible survocc members so far plus the same 3 trap members.
# 024g: the production member plus the same 3 trap members.
MEMBERS=$(cat <<'TABLE'
024e 0.05  3  365 0    1.0 1.0
024e 0.05  3  0   0    1.0 1.0
024e 0.075 3  365 0    1.0 1.0
024e 0.075 3  0   0    1.0 1.0
024e 0.15  3  365 0    1.0 1.0
024e 0.15  3  0   0    1.0 1.0
024e 0.1   12 0   0    1.0 1.0
024e 0.1   1  0   0    1.0 1.0
024e 0.1   3  365 0    1.0 1.0
024e 0.1   3  0   0    1.0 1.0
024e 0.1   3  0   0.02 1.0 1.0
024e 0.1   3  0   0    1.0 0.0
024e 0.1   3  0   0    1.0 0.25
024e 0.1   3  0   0    1.0 0.5
024f 0.075 3  0   0    1.0 1.0
024f 0.15  3  0   0    1.0 1.0
024f 0.1   3  365 0    1.0 1.0
024f 0.1   3  0   0    1.0 1.0
024f 0.1   3  0   0    1.0 0.0
024f 0.1   3  0   0    1.0 0.25
024f 0.1   3  0   0    1.0 0.5
024g 0.1   3  0   0    1.0 1.0
024g 0.1   3  0   0    1.0 0.0
024g 0.1   3  0   0    1.0 0.25
024g 0.1   3  0   0    1.0 0.5
TABLE
)

# Blacklisted nodes (see scripts/node-blacklist.txt), for the provenance line.
exclude=$(sed "s/#.*//" scripts/node-blacklist.txt | tr -d "[:blank:]" \
          | grep -v "^$" | paste -sd,)

export output_root regime hex_radius band_m logdir exclude
export max_float_days occupancy_max_days connectivity_max_days

# awk's %g is printf's %g is Python's :g, so this reproduces the notebook's
# member tag character for character.
fmt_g() { awk -v x="$1" 'BEGIN { printf "%g", x }'; }
export -f fmt_g

member_tag() {
    local w_c=$1 ts=$2 tc=$3 delta=$4 tf=$5 tw=$6 m
    m="step_wc$(fmt_g "${w_c}")_ts$(fmt_g "${ts}")"
    m="${m}_tc$(awk -v x="${tc}" 'BEGIN { if (x > 0) printf "%g", x; else printf "inf" }')"
    awk -v x="${delta}" 'BEGIN { exit !(x > 0) }' && m="${m}_d$(fmt_g "${delta}")"
    # Only a shore-type-sensitive member carries the trap weights.
    if awk -v a="${tf}" -v b="${tw}" 'BEGIN { exit !(a != 1.0 || b != 1.0) }'; then
        m="${m}_tf$(fmt_g "${tf}")_tw$(fmt_g "${tw}")"
    fi
    echo "${m//./p}"
}
export -f member_tag

run_one() {
    reducer=$1 w_c=$2 ts=$3 tc=$4 delta=$5 tf=$6 tw=$7 year=$8 month=$9
    t0=${SECONDS}
    tag=$(member_tag "${w_c}" "${ts}" "${tc}" "${delta}" "${tf}" "${tw}")
    ms=$(printf "m%02d" "${month}")

    case "${reducer}" in
        024e) nb=notebooks/024e_BuildBeaching.ipynb
              horizon=(-p max_float_days "${max_float_days}") ;;
        024f) nb=notebooks/024f_BuildSurvivalOccupancy.ipynb
              horizon=(-p occupancy_max_days "${occupancy_max_days}") ;;
        024g) nb=notebooks/024g_BuildSurvivalConnectivity.ipynb
              horizon=(-p connectivity_max_days "${connectivity_max_days}") ;;
        *)    echo "FAIL ${reducer} ? ${year} ${month} (unknown reducer)"; return 0 ;;
    esac

    # Per-task Jupyter runtime dir. jupyter_client reserves five ZMQ ports by
    # binding to port 0 and closing them; the kernel re-binds moments later and
    # a concurrent kernel on the same host can steal one in that window. The
    # private dir keeps connection files off GPFS; the 3-attempt retry below is
    # what actually covers the race (port names carry UUIDs).
    rt="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/jupyter-runtime-${reducer}-${tag}-${year}-${ms}"
    mkdir -p "${rt}"
    rc=1
    for attempt in 1 2 3; do
        srun -n1 -N1 --exact --no-kill \
            --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
            --mem-per-cpu="${SLURM_MEM_PER_CPU}M" \
            --job-name="${reducer}_${tag}" \
            env OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                JUPYTER_RUNTIME_DIR="${rt}" \
            pixi run papermill --cwd notebooks/ \
                "${nb}" \
                "${logdir}/executed/${reducer}_${regime}_${year}_${ms}_${tag}_r${hex_radius}m.ipynb" \
                -p output_root "${output_root}" \
                -r regime "${regime}" \
                -p release_year "${year}" \
                -p release_month "${month}" \
                -p hex_radius "${hex_radius}" \
                -p w_c "${w_c}" \
                -p tau_strong_hours "${ts}" \
                -p tau_calm_days "${tc}" \
                -p delta "${delta}" \
                -p trap_flat "${tf}" \
                -p trap_wall "${tw}" \
                -p band_m "${band_m}" \
                "${horizon[@]}" \
                -k python && { rc=0; break; }
        rc=$?
        echo "attempt ${attempt} failed (rc=${rc}) ${reducer} ${tag} ${year} ${month}" >&2
        sleep $(( RANDOM % 10 + 1 ))
    done
    if [ ${rc} -eq 0 ]; then
        echo "OK ${reducer} ${tag} ${year} ${month} $(( SECONDS - t0 ))"
    else
        echo "FAIL ${reducer} ${tag} ${year} ${month}"
    fi
    return 0
}
export -f run_one

# "reducer w_c ts tc delta tf tw year month" per line.
list_work() {
    while read -r reducer w_c ts tc delta tf tw; do
        [ -n "${reducer}" ] || continue
        for year in "${YEARS[@]}"; do
            for month in "${MONTHS[@]}"; do
                echo "${reducer} ${w_c} ${ts} ${tc} ${delta} ${tf} ${tw} ${year} ${month}"
            done
        done
    done <<< "${MEMBERS}"
}

njobs=$(list_work | wc -l)
# Provenance header: the store parquets carry git_sha + slurm_job_id as attrs,
# and this pins what that job was, in the .out next to the OK/FAIL tally.
echo "024efg provenance: git $(git rev-parse HEAD)"
if [ -n "$(git status --porcelain)" ]; then
    echo "024efg WARNING: dirty working tree, this run maps to no commit"
fi
echo "024efg slurm: job ${SLURM_JOB_ID}, ntasks ${SLURM_NTASKS}, cpus-per-task"\
     "${SLURM_CPUS_PER_TASK}, mem-per-cpu ${SLURM_MEM_PER_CPU}M, exclude ${exclude:-<none>}"
echo "024efg members:"
echo "${MEMBERS}" | while read -r reducer w_c ts tc delta tf tw; do
    [ -n "${reducer}" ] && echo "  ${reducer} $(member_tag "${w_c}" "${ts}" "${tc}" "${delta}" "${tf}" "${tw}")"
done
echo "024efg: ${njobs} cells, ${SLURM_NTASKS} concurrent, regime=${regime}," \
     "hex_radius=${hex_radius}, years ${YEARS[*]}, logs -> ${logdir}"

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
echo "024efg done: ${nok} OK, ${nfail} FAIL of ${njobs}"
grep '^FAIL ' "${tally}" | sed 's/^FAIL //' || true

jobinfo || true
[ "${nfail}" -eq 0 ]
