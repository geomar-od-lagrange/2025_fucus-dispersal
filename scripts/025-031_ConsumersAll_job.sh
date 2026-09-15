#!/bin/bash
#SBATCH --job-name=025-031_consumers
#SBATCH --partition=base
# Concurrency ceiling for a direct `sbatch`. scripts/submit_downstream_chain.sh
# sizes the allocation from the work list instead (capped at the same value).
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=16G
#SBATCH --time=08:00:00
#SBATCH --distribution=cyclic
# Do not abort the whole allocation when one node dies (see the 024d job
# script): the surviving steps run on and the lost units come back as FAIL.
#SBATCH --no-kill
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/025-031/%x_%j.out

# Every parquet-only consumer figure of the study, in one job.
#
# Unit of work = one (stage, regime, season, member) run of an existing
# scripts/0NN_*_job.sh, launched as its own `srun` step. The per-stage scripts
# stay the single source of the papermill contract (notebook path, parameter
# names, `-r` vs `-p`, executed-notebook name); this driver only decides which
# combinations run. The alternative -- ~80 separate `sbatch` submissions --
# means 80 allocations queueing independently behind the same dependency, for
# units that are minutes long.
#
# The matrix:
#
#   025, 026, 026a, 026b, 027  x {surface_stokes, bottom} x 5 seasons   50
#   028  (surface_stokes unweighted, surface_stokes production member,
#         bottom unweighted)                              x 5 seasons   15
#   029, 030  surface_stokes, production member           x 5 seasons   10
#   031  surface_stokes, every 024e member on disk        x 5 seasons    5
#                                                                     ----
#                                                                       80
#
# `surface` (no Stokes) is not produced further (plans/analysis_wrapup.md §0),
# and only surface_stokes has 024d sidecars, hence the beaching/survival stages
# and the weighted 028 are Stokes-only.
#
# Nothing here names a release year: every consumer pools whatever release
# years the HexAggregates store holds, so the figures span the years that
# landed, and a release year that failed upstream is simply absent.
#
# Memory: 8 cpu x 16G = 128 GB/task, the shape 026/026a/026b already ask for
# standalone, and the reason the whole bundle is not sized like the cheap
# stages. Measured `sacct` MaxRSS over the four-release-year stores:
# 026 56.3 GiB, 026b 49.6 GiB, 026a 40.2 GiB, 025 32.8 GiB, 027 0.2 GiB,
# 028 <0.1 GiB -- the counts-store consumers hold the pooled store in memory,
# so their peak grows with the release-year count (four years there, up to ten
# here). 128 GB/task is that peak x ~2.3 and still one task per 192 GB node.
# If a unit OOMs, the fix is `--mem-per-cpu=24G --ntasks=8`, not a bigger
# default: the cheap two thirds of the matrix would pay for it.
#
# Usage:
#   sbatch scripts/025-031_ConsumersAll_job.sh
#   LIST_ONLY=1 bash scripts/025-031_ConsumersAll_job.sh    # print the units
#   STAGES="028 031" sbatch scripts/025-031_ConsumersAll_job.sh  # subset
#   SEASONS="ALL SON" sbatch scripts/025-031_ConsumersAll_job.sh
#
# The 024/024b/024c stores, the 024a key and (for 029/030/031 and the weighted
# 028) the 024e/024f/024g stores must exist for the matching hex_radius --
# scripts/submit_downstream_chain.sh is what wires that ordering.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$PWD}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
hex_radius=${HEX_RADIUS:-6000}
# The rate-model member the study reports (plans/analysis_wrapup.md §0); the
# other members only appear in the 031 sweep.
production_member=${PRODUCTION_MEMBER:-step_wc0p1_ts3_tcinf}
stages=${STAGES:-025 026 026a 026b 027 028 029 030 031}
regimes=${REGIMES:-surface_stokes bottom}
seasons=${SEASONS:-ALL DJF MAM JJA SON}
logdir=${output_root}/logs/025-031
mkdir -p "${logdir}"

# Drop any CPU-bind mask inherited from an outer allocation; otherwise the
# concurrent srun steps below fail with "CPU binding outside of job step
# allocation". --exact sets each step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

# Blacklisted nodes (see scripts/node-blacklist.txt), for the provenance line.
# The exclusion itself is applied to the allocation by the submitter; repeating
# it as a step-level `srun --exclude` makes SLURM serialise the steps.
exclude=$(sed "s/#.*//" scripts/node-blacklist.txt | tr -d "[:blank:]" \
          | grep -v "^$" | paste -sd,)

export output_root repo_root hex_radius production_member logdir exclude

# Every 024e member on disk, as the comma-separated list 031 sweeps over. Read
# at run time, not at submit time: the 024e/f/g reducers write these minutes
# before this job starts, so a hardcoded list would be stale by construction.
beaching_members_csv() {
    local f b
    for f in "${output_root}/HexAggregates"/HexAgg_beaching_r"${hex_radius}"m_surface_stokes_[0-9][0-9][0-9][0-9]_m[0-9][0-9]_*.parquet; do
        [ -e "${f}" ] || continue
        b=$(basename "${f}" .parquet)
        sed -E "s/^HexAgg_beaching_r[0-9]+m_surface_stokes_[0-9]{4}_m[0-9]{2}_//" <<< "${b}"
    done | sort -u | paste -sd,
}
export -f beaching_members_csv

run_one() {
    stage=$1 regime=$2 season=$3 member=$4
    t0=${SECONDS}
    # "-" is how the work list spells "this stage takes no member"; 028 reads
    # the unweighted 024c store for an empty member, so the two are distinct.
    [ "${member}" = "-" ] && member=""

    case "${stage}" in
        025)  job=scripts/025_HexHeatmaps_job.sh
              args=("${regime}" "${season}" "${hex_radius}") ;;
        026)  job=scripts/026_TimeHorizonMaps_job.sh
              args=("${regime}" "${season}" "${hex_radius}") ;;
        026a) job=scripts/026a_OriginSubbasinTimeHorizonMaps_job.sh
              args=("${regime}" "${season}" "${hex_radius}") ;;
        026b) job=scripts/026b_OriginSubbasinYearTimeHorizonMaps_job.sh
              args=("${regime}" "${season}" "${hex_radius}") ;;
        027)  job=scripts/027_HexDistanceQuantiles_job.sh
              args=("${regime}" "${hex_radius}" "${season}") ;;
        028)  job=scripts/028_SubbasinConnectivityMatrix_job.sh
              args=("${regime}" "${hex_radius}" "${season}" "${member}") ;;
        029)  job=scripts/029_BeachingMaps_job.sh
              args=("${regime}" "${hex_radius}" "${member}" "${season}") ;;
        030)  job=scripts/030_SurvivalHeatmaps_job.sh
              args=("${regime}" "${hex_radius}" "${member}" "${season}") ;;
        031)  job=scripts/031_BeachingSweep_job.sh
              args=("${regime}" "${hex_radius}" "$(beaching_members_csv)" "" \
                    "${production_member}" "${season}") ;;
        *)    echo "FAIL ${stage} ${regime} ${season} ${member:-unweighted} (unknown stage)"
              return 0 ;;
    esac

    unit="${stage}_${regime}_${season}_${member:-unweighted}"
    # Per-task Jupyter runtime dir. jupyter_client reserves five ZMQ ports by
    # binding to port 0 and closing them; the kernel re-binds moments later and
    # a concurrent kernel on the same host can steal one in that window. The
    # private dir keeps connection files off GPFS; the 3-attempt retry below is
    # what actually covers the race (port names carry UUIDs).
    rt="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/jupyter-runtime-${unit}"
    mkdir -p "${rt}"
    rc=1
    for attempt in 1 2 3; do
        srun -n1 -N1 --exact --no-kill \
            --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
            --mem-per-cpu="${SLURM_MEM_PER_CPU}M" \
            --job-name="${unit}" \
            env OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}" \
                JUPYTER_RUNTIME_DIR="${rt}" \
            bash "${job}" "${args[@]}" && { rc=0; break; }
        rc=$?
        echo "attempt ${attempt} failed (rc=${rc}) ${unit}" >&2
        sleep $(( RANDOM % 10 + 1 ))
    done
    if [ ${rc} -eq 0 ]; then
        echo "OK ${stage} ${regime} ${season} ${member:-unweighted} $(( SECONDS - t0 ))"
    else
        echo "FAIL ${stage} ${regime} ${season} ${member:-unweighted}"
    fi
    return 0
}
export -f run_one

# "<stage> <regime> <season> <member>" per line, "-" where a stage takes no
# member. $STAGES / $REGIMES / $SEASONS restrict the matrix.
list_work() {
    local stage regime season
    for stage in ${stages}; do
        case "${stage}" in
            025|026|026a|026b|027)
                for regime in ${regimes}; do
                    for season in ${seasons}; do
                        echo "${stage} ${regime} ${season} -"
                    done
                done ;;
            028)
                for season in ${seasons}; do
                    # surface_stokes twice: the unweighted 024c matrix and the
                    # survival-weighted 024g one, on one population.
                    echo "028 surface_stokes ${season} -"
                    echo "028 surface_stokes ${season} ${production_member}"
                    echo "028 bottom ${season} -"
                done ;;
            029|030)
                for season in ${seasons}; do
                    echo "${stage} surface_stokes ${season} ${production_member}"
                done ;;
            031)
                for season in ${seasons}; do
                    echo "031 surface_stokes ${season} -"
                done ;;
            *)  echo "${stage} ? ? -" ;;
        esac
    done
}

# Work-list preview, for sizing the allocation
# (scripts/submit_downstream_chain.sh) or eyeballing a $STAGES selection.
if [ -n "${LIST_ONLY:-}" ]; then
    list_work
    exit 0
fi

njobs=$(list_work | wc -l)
# Provenance header: the figures and exports carry no job id of their own, so
# this .out is what pins a PNG to a commit and a dependency.
echo "025-031 provenance: git $(git rev-parse HEAD)"
if [ -n "$(git status --porcelain)" ]; then
    echo "025-031 WARNING: dirty working tree, this run maps to no commit"
fi
echo "025-031 slurm: job ${SLURM_JOB_ID}, ntasks ${SLURM_NTASKS}, cpus-per-task"\
     "${SLURM_CPUS_PER_TASK}, mem-per-cpu ${SLURM_MEM_PER_CPU}M, exclude ${exclude:-<none>}"
echo "025-031 dependency: $(scontrol show job "${SLURM_JOB_ID}" 2>/dev/null \
      | sed -n 's/.*Dependency=\([^ ]*\).*/\1/p' | head -1)"
echo "025-031: ${njobs} units, ${SLURM_NTASKS} concurrent, hex_radius=${hex_radius},"\
     "production member ${production_member}, logs -> ${logdir}"
echo "025-031 031 sweep members: $(beaching_members_csv)"
ls "${output_root}/HexAggregates"/HexAgg_beaching_r"${hex_radius}"m_surface_stokes_*_"${production_member}".parquet \
    >/dev/null 2>&1 \
    || echo "025-031 WARNING: no 024e store for production member ${production_member}"

# `set -e` must not kill the driver when a unit exhausts its retries; run_one
# already swallows its own status and reports FAIL, and the tally below is the
# authoritative exit code.
set +e
list_work | xargs -P "${SLURM_NTASKS}" -L1 bash -c 'run_one "$@"' _ \
    | tee "${logdir}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.tally"
set -e

tally=${logdir}/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.tally
nok=$(grep -c '^OK ' "${tally}" || true)
nfail=$(grep -c '^FAIL ' "${tally}" || true)
echo "025-031 done: ${nok} OK, ${nfail} FAIL of ${njobs}"
grep '^FAIL ' "${tally}" | sed 's/^FAIL //' || true

jobinfo || true
[ "${nfail}" -eq 0 ]
