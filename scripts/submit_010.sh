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
#
# `VELOCITY_SCALE` (bottom only) runs the bottom fields at a uniform velocity
# factor, written under its own regime (0.8 -> `bottom_vs0p8`); it goes into
# the job name so a scaled run never adopts or blocks the 1x one.
#
#   RELEASE_NS=0,7,12 VELOCITY_SCALE=0.8 scripts/submit_010.sh bottom 2019
set -euo pipefail
cd "$(dirname "$0")/.."
source scripts/slurm_submit.sh

regime=${1:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
year=${2:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
n_releases=${3:-73}
shift $(( $# < 3 ? $# : 3 ))

job=scripts/010_FucusDispersal_${regime}_job.sh
[ -f "${job}" ] || { echo "no such regime: ${regime} (${job})" >&2; exit 1; }
# A job script that never reads VELOCITY_SCALE would run 1x under a scaled name.
if [ -n "${VELOCITY_SCALE:-}" ] && ! grep -q VELOCITY_SCALE "${job}"; then
    echo "${job} does not take VELOCITY_SCALE" >&2; exit 1
fi

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
    name=010-${regime}${VELOCITY_SCALE:+-vs${VELOCITY_SCALE}}-${year}-n${release_ns// /_}
    selection="releases=[${release_ns}]"
else
    [ "${n_releases}" -ge 1 ] && [ "${n_releases}" -le 73 ] \
        || { echo "n_releases must be 1..73, got ${n_releases}" >&2; exit 1; }
    ntasks=${n_releases}
    name=010-${regime}${VELOCITY_SCALE:+-vs${VELOCITY_SCALE}}-${year}
    selection="n_releases=${n_releases}"
fi

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)

# #SBATCH --output opens its target when the job starts, before the job
# script's own commands run -- the directory has to exist ahead of
# submission, same as scripts/submit_downstream_chain.sh's logroot mkdir.
mkdir -p /gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/010

echo "regime=${regime} velocity_scale=${VELOCITY_SCALE:-1.0} year=${year} ${selection} ntasks=${ntasks}" \
     "excluding: ${exclude:-<none>}"

# RELEASE_NS and VELOCITY_SCALE reach the job script through the environment
# sbatch propagates; naming them in --export alongside ALL keeps that explicit
# rather than implied. Names only, no `=value`: --export splits on commas, so
# a RELEASE_NS value of "1,19" would be cut at its first comma.
export_vars=ALL${RELEASE_NS:+,RELEASE_NS}${VELOCITY_SCALE:+,VELOCITY_SCALE}
robust_sbatch --job-name="${name}" --ntasks="${ntasks}" \
    --export="${export_vars}" \
    ${exclude:+--exclude="${exclude}"} "$@" \
    "${job}" "${year}" "${n_releases}"
