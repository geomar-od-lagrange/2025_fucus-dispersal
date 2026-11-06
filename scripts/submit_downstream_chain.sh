#!/bin/bash
# Submit everything downstream of a set of 010 Parcels runs, as one SLURM
# dependency graph, and exit. Nothing here waits, polls or watches: the 010
# runs take 24-48 h, so the chain has to be queued in one go and left alone.
#
#   scripts/submit_downstream_chain.sh <010-job-id> [<010-job-id> ...]
#   scripts/submit_downstream_chain.sh 23836364,23836392,23836393
#
# The stages, and why each dependency type:
#
#   A  024 counts + 024b distance + 024c connectivity, one triple per
#      (regime, year), via scripts/submit_aggregates.sh.
#      afterany on that pair's own 010 job -- afterok would strand the whole
#      chain on one failed release year, and the notebooks discover the
#      trajectory zarrs that exist rather than demanding a full 73.
#   B  024d beaching sidecars, one job over every surface_stokes zarr on disk,
#      afterany on the surface_stokes 010 jobs only. Beaching is Stokes-only:
#      024d and the 024efg reducers pin regime=surface_stokes internally.
#   C  024e/024f/024g reducers (scripts/submit_024efg.sh), afterok on B --
#      a reducer over a half-written sidecar set is not worth having.
#   D  the parquet-only consumers 025..031, one multi-task job
#      (scripts/025-031_ConsumersAll_job.sh), afterok on every A job and on C.
#      afterok is deliberate and strict: a figure pooling a store that failed
#      to build is a wrong figure, not a partial one. The cost is that one
#      failed A job leaves D pending on an unsatisfiable dependency --
#      --kill-on-invalid-dep=yes turns that into a cancelled job rather than a
#      zombie, and the fix is to rerun the failed unit and resubmit D alone.
#   E  the headline regeneration, afterok on D.
#
# Every submission goes through robust_sbatch (scripts/slurm_submit.sh), so a
# controller that times out mid-chain is survivable: each stage carries a job
# name that is unique to it and stable across reruns
# ("<tag>-<regime>-<year>-024b", "<tag>-024d", ...), where <tag> is a hash of
# the 010 id set, and the queue -- not sbatch's exit status -- decides whether
# a stage still needs submitting. That makes the whole chain idempotent: rerun
# the identical command after a partial failure and it adopts what is already
# queued and submits only the rest. The table at the end says which was which.
# A stage that has already *finished* is no longer in the queue and would be
# submitted again, so resume while the chain is still pending.
#
# The (regime, year) pairs are read off the 010 jobs themselves -- `scontrol
# show job <id>` still knows the Command line, which carries the regime (in the
# script name) and the release year (first argument). When scontrol no longer
# knows a job (purged after completion), pass the pairs explicitly:
#
#   REGIME_YEARS="bottom:2020 bottom:2021 surface_stokes:2020" \
#       scripts/submit_downstream_chain.sh 23836364 23836395
#
# With REGIME_YEARS set, every pair depends on the whole 010 id set instead of
# on its own job, since the mapping is no longer known.
#
# Environment:
#   REGIME_YEARS          "<regime>:<year> ..." instead of scontrol discovery
#   HEX_RADIUS            hex radius for every stage (default 6000)
#   MAX_CONSUMER_TASKS    concurrency cap for stage D (default 12)
#   HEADLINE_JOB          stage E job script (default scripts/headline_report_job.sh)
#   CHAIN_TAG             job-name prefix (default: hash of the 010 id set)
#   ADOPT_A ADOPT_B ADOPT_C
#                         ids of stages already queued under names this chain
#                         cannot reconstruct -- an earlier chain, or a hand
#                         submission. Such a stage is not submitted at all and
#                         the later stages depend on these ids. Check them
#                         yourself (scontrol show job <id>): they are taken on
#                         trust.
#   SBATCH_RETRY_DELAYS   submission backoff, see scripts/slurm_submit.sh
#
# A dry run prints the whole graph without touching the queue: put a fake
# `sbatch` that echoes its arguments first on $PATH.
set -euo pipefail
cd "$(dirname "$0")/.."
source "$(dirname "$0")/slurm_submit.sh"

[ $# -ge 1 ] || {
    echo "usage: scripts/submit_downstream_chain.sh <010-job-id> [<010-job-id> ...]" >&2
    exit 1
}

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
hex_radius=${HEX_RADIUS:-6000}
max_consumer_tasks=${MAX_CONSUMER_TASKS:-12}
headline_job=${HEADLINE_JOB:-scripts/headline_report_job.sh}
consumers_job=scripts/025-031_ConsumersAll_job.sh

logroot=${output_root}/logs
# SLURM opens a job's --output before the job runs, so a stage whose #SBATCH
# directive names a log dir cannot create it itself: pre-create every one here.
mkdir -p "${logroot}/aggregates" "${logroot}/024d" "${logroot}/024efg" \
         "${logroot}/025-031" "${logroot}/032" "${logroot}/chain"
chainlog=${logroot}/chain/chain_$(date +%Y%m%dT%H%M%S).log
exec > >(tee -a "${chainlog}") 2>&1

# One id per line, from any mix of repeated arguments, commas and whitespace.
ids=$(printf '%s\n' "$@" | tr ', \t' '\n' | grep -v '^$' | sort -u)
for id in ${ids}; do
    [[ "${id}" =~ ^[0-9]+$ ]] || { echo "not a job id: ${id}" >&2; exit 1; }
done
all_dep=$(paste -sd: <<< "${ids}")
# Stable across reruns of the same chain, distinct between different ones.
tag=${CHAIN_TAG:-fx$(md5sum <<< "${ids}" | cut -c1-8)}

echo "downstream chain: git $(git rev-parse HEAD)"
if [ -n "$(git status --porcelain)" ]; then
    echo "WARNING: dirty working tree, this chain maps to no commit"
fi
echo "010 jobs: $(paste -sd, <<< "${ids}")"
echo "hex_radius=${hex_radius} output_root=${output_root}"
echo "job-name tag: ${tag} (rerun this command to adopt whatever is queued)"
echo

# --- which (regime, year) pairs, and what each waits for ------------------
# "<regime> <year> <dep:dep:...>" per line.
if [ -n "${REGIME_YEARS:-}" ]; then
    pairs=$(printf '%s\n' ${REGIME_YEARS} | tr ',' '\n' | grep -v '^$' \
            | awk -F: -v d="${all_dep}" '{print $1, $2, d}' | sort -u)
    echo "pairs from \$REGIME_YEARS, each waiting on the whole 010 set"
else
    raw=""
    for id in ${ids}; do
        cmd=$(scontrol show job "${id}" 2>/dev/null | sed -n 's/^ *Command=//p' | head -1)
        [ -n "${cmd}" ] || { echo "scontrol no longer knows job ${id}; pass REGIME_YEARS" >&2; exit 1; }
        script=$(basename "$(awk '{print $1}' <<< "${cmd}")")
        year=$(awk '{print $2}' <<< "${cmd}")
        regime=${script#010_FucusDispersal_}
        regime=${regime%_job.sh}
        [ "${regime}" != "${script}" ] && [[ "${year}" =~ ^[0-9]{4}$ ]] || {
            echo "job ${id} is not a 010 run: ${cmd}" >&2; exit 1; }
        raw+="${regime} ${year} ${id}"$'\n'
    done
    # One line per pair, with every 010 job that covers it (a resubmitted year
    # contributes a second id).
    pairs=$(grep -v '^$' <<< "${raw}" \
        | awk '{k=$1" "$2; d[k]=(k in d ? d[k]":"$3 : $3)} END {for (k in d) print k, d[k]}' \
        | sort)
    echo "pairs discovered from scontrol:"
fi
sed 's/^/  /' <<< "${pairs}"
echo

stokes_dep=$(awk '$1=="surface_stokes" {print $3}' <<< "${pairs}" | tr ':' '\n' | sort -u | paste -sd:)
[ -n "${stokes_dep}" ] || { echo "no surface_stokes 010 job: stage B/C/D would have no sidecars" >&2; exit 1; }

# The static hex key is not part of this chain -- it is built once per radius.
[ -f "${output_root}/HexAggregates/HexAgg_key_r${hex_radius}m.parquet" ] || {
    echo "missing 024a key HexAgg_key_r${hex_radius}m.parquet; run 024a_BuildHexKey_job.sh first" >&2
    exit 1; }

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)

graph=""   # "<stage>|<job id>|<how>|<what>|<dependency>" per line
chain_complete=""

# Printed however the script ends, so an aborted run still says what it queued
# and how to finish it.
summary() {
    echo
    printf '%-2s %-10s %-9s %-44s %s\n' stage job how what dependency
    printf '%-2s %-10s %-9s %-44s %s\n' -- ---------- --------- \
        -------------------------------------------- ----------
    printf '%s' "${graph}" \
        | awk -F'|' '{printf "%-2s %-10s %-9s %-44s %s\n", $1, $2, $3, $4, $5}'
    echo
    echo "chain log: ${chainlog}"
    if [ -z "${chain_complete}" ]; then
        echo
        echo "INCOMPLETE: only the stages above are queued. Rerun the identical"
        echo "command to finish it -- queued stages are adopted by job name:"
        echo "  scripts/submit_downstream_chain.sh $(paste -sd' ' <<< "${ids}")"
    fi
}
trap summary EXIT

# --- stage A: per-(regime, year) aggregates -------------------------------
a_ids=()
if [ -n "${ADOPT_A:-}" ]; then
    read -ra a_ids <<< "$(tr ',' ' ' <<< "${ADOPT_A}")"
    for j in "${a_ids[@]}"; do
        graph+="A|${j}|adopted|024/024b/024c (\$ADOPT_A)|as queued"$'\n'
    done
fi
while [ -z "${ADOPT_A:-}" ] && read -r regime year dep; do
    [ -n "${regime}" ] || continue
    # `|| true`: a triple that the controller cut off half-way still prints the
    # ids it did queue, and those go into the table below before this run
    # stops -- a rerun then adopts them by name.
    out=$(ALLOW_MISSING=1 NAME_PREFIX="${tag}-${regime}-${year}-" \
        scripts/submit_aggregates.sh "${regime}" "${year}" "${hex_radius}" \
            --dependency=afterany:"${dep}" \
            --output="${logroot}/aggregates/%x_%j.out") || true
    read -r counts counts_how <<< "$(awk '$1=="024"  {print $3, $4}' <<< "${out}")"
    read -r dist   dist_how   <<< "$(awk '$1=="024b" {print $3, $4}' <<< "${out}")"
    read -r conn   conn_how   <<< "$(awk '$1=="024c" {print $3, $4}' <<< "${out}")"
    [ -n "${counts}" ] && graph+="A|${counts}|${counts_how}|024 counts ${regime} ${year}|afterany:${dep}"$'\n'
    [ -n "${dist}" ] && graph+="A|${dist}|${dist_how}|024b distance ${regime} ${year}|afterany:${dep}"$'\n'
    [ -n "${conn}" ] && graph+="A|${conn}|${conn_how}|024c connectivity ${regime} ${year}|afterok:${counts}"$'\n'
    for j in "${counts}" "${dist}" "${conn}"; do
        [ -n "${j}" ] || { echo "submit_aggregates.sh ${regime} ${year} did not queue its full triple" >&2; exit 1; }
    done
    a_ids+=("${counts}" "${dist}" "${conn}")
done <<< "${pairs}"

# --- stage B: beaching sidecars (surface_stokes only) ---------------------
if [ -n "${ADOPT_B:-}" ]; then
    b_id=${ADOPT_B}
    graph+="B|${b_id}|adopted|024d beaching sidecars (\$ADOPT_B)|as queued"$'\n'
else
    read -r b_id b_how <<< "$(scripts/submit_024d.sh --parsable \
        --job-name="${tag}-024d" \
        --dependency=afterany:"${stokes_dep}" | tail -n1)"
    [ -n "${b_id}" ] || { echo "submit_024d.sh returned no job id" >&2; exit 1; }
    graph+="B|${b_id}|${b_how}|024d beaching sidecars surface_stokes|afterany:${stokes_dep}"$'\n'
fi

# --- stage C: the 024e/024f/024g reducers ---------------------------------
if [ -n "${ADOPT_C:-}" ]; then
    c_id=${ADOPT_C}
    graph+="C|${c_id}|adopted|024efg reducers (\$ADOPT_C)|as queued"$'\n'
else
    read -r c_id c_how <<< "$(scripts/submit_024efg.sh --parsable \
        --job-name="${tag}-024efg" --dependency=afterok:"${b_id}" \
        --kill-on-invalid-dep=yes | tail -n1)"
    [ -n "${c_id}" ] || { echo "submit_024efg.sh returned no job id" >&2; exit 1; }
    graph+="C|${c_id}|${c_how}|024efg reducers, every member|afterok:${b_id}"$'\n'
fi

# --- stage D: the parquet-only consumers ----------------------------------
# The work list is static (the 031 member set is resolved inside the job, at
# run time), so it can be counted now to size the allocation.
nunits=$(LIST_ONLY=1 bash "${consumers_job}" | wc -l)
ntasks=$(( nunits < max_consumer_tasks ? nunits : max_consumer_tasks ))
d_dep=$(printf '%s\n' "${a_ids[@]}" "${c_id}" | paste -sd:)
read -r d_id d_how <<< "$(robust_sbatch --ntasks="${ntasks}" \
        ${exclude:+--exclude="${exclude}"} \
        --job-name="${tag}-consumers" \
        --dependency=afterok:"${d_dep}" --kill-on-invalid-dep=yes \
        "${consumers_job}")"
graph+="D|${d_id}|${d_how}|025-031 consumers, ${nunits} units, ${ntasks} tasks|afterok:${d_dep}"$'\n'

# --- stage E: the headline regeneration -----------------------------------
# ==== HOOK ====  Owned by the headline agent. Wired the moment the script
# exists under this name; until then the chain simply has no stage E and says
# so, rather than queueing a job that cannot run.
if [ -f "${headline_job}" ]; then
    # No --output override: the stage's own #SBATCH directive names its log
    # dir, as every other stage's does, and two names for one stage's logs is
    # how a failed job ends up with no log at all.
    read -r e_id e_how <<< "$(robust_sbatch ${exclude:+--exclude="${exclude}"} \
            --job-name="${tag}-headline" \
            --dependency=afterok:"${d_id}" --kill-on-invalid-dep=yes \
            "${headline_job}")"
    graph+="E|${e_id}|${e_how}|headline (${headline_job})|afterok:${d_id}"$'\n'
else
    echo "WARNING: stage E skipped -- no ${headline_job} on disk."
    echo "WARNING: regenerate the headline by hand, or submit it with"
    echo "WARNING:   sbatch --dependency=afterok:${d_id} ${headline_job}"
    graph+="E|-|-|headline SKIPPED (no ${headline_job})|afterok:${d_id}"$'\n'
fi

chain_complete=1
