#!/bin/bash
# Submit one 010 Parcels run: `scripts/submit_010.sh <regime> <year>
# [n_releases] [extra sbatch args...]`.
#
# It exists for one reason: `--ntasks` and the release count have to agree, and
# an `#SBATCH` directive cannot see the job script's own arguments. The wrapper
# sets both from the same selection, so a partial year (the last forced year,
# where only the releases whose full simulation window fits inside the available
# forcing may run) books exactly as many tasks as it has work. It also expands
# scripts/node-blacklist.txt into `--exclude=`, which `#SBATCH` cannot read from
# a file either.
#
# `RELEASE_NS` picks release indices explicitly instead of a prefix, for
# refilling holes: a release lost mid-schedule (a dead node under `--no-kill`, a
# kernel that failed to start) is not a missing tail, so `n_releases` cannot
# select it. `--ntasks` then comes from the length of the list.
#
#   scripts/submit_010.sh bottom 2022                  # full year, 73 tasks
#   scripts/submit_010.sh surface_stokes 2025 29       # partial year, 29 tasks
#   scripts/submit_010.sh bottom 2023 73 --dependency=afterok:12345
#   RELEASE_NS=1,19 scripts/submit_010.sh surface_stokes 2022   # 2 tasks
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/slurm_submit.sh

regime=${1:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
year=${2:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
n_releases=${3:-73}
shift $(( $# < 3 ? $# : 3 ))

job=scripts/010_FucusDispersal_${regime}_job.sh
[ -f "${job}" ] || { echo "no such regime: ${regime} (${job})" >&2; exit 1; }

# The job name is robust_sbatch's idempotency key, so it has to name the work,
# not the attempt: a rerun of the same releases after a controller timeout must
# hash to the same name and get adopted rather than queued twice.
if [ -n "${RELEASE_NS:-}" ]; then
    release_ns=${RELEASE_NS//,/ }
    for n in ${release_ns}; do
        [ "${n}" -ge 0 ] && [ "${n}" -le 72 ] \
            || { echo "release index must be 0..72, got ${n}" >&2; exit 1; }
    done
    ntasks=$(wc -w <<< "${release_ns}")
    name=010-${regime}-${year}-n${release_ns// /_}
    selection="releases=[${release_ns}]"
else
    [ "${n_releases}" -ge 1 ] && [ "${n_releases}" -le 73 ] \
        || { echo "n_releases must be 1..73, got ${n_releases}" >&2; exit 1; }
    ntasks=${n_releases}
    name=010-${regime}-${year}
    selection="n_releases=${n_releases}"
fi

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)

# #SBATCH --output opens its target when the job starts, before the job
# script's own commands run -- the directory has to exist ahead of
# submission, same as scripts/submit_downstream_chain.sh's logroot mkdir.
mkdir -p /gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/010

echo "regime=${regime} year=${year} ${selection} ntasks=${ntasks}" \
     "excluding: ${exclude:-<none>}"

# RELEASE_NS reaches the job script through the environment sbatch propagates;
# naming it in --export alongside ALL keeps that explicit rather than implied.
robust_sbatch --job-name="${name}" --ntasks="${ntasks}" \
    ${RELEASE_NS:+--export=ALL,RELEASE_NS="${release_ns// /,}"} \
    ${exclude:+--exclude="${exclude}"} "$@" \
    "${job}" "${year}" "${n_releases}"
