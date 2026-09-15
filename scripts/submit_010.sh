#!/bin/bash
# Submit one 010 Parcels run: `scripts/submit_010.sh <regime> <year>
# [n_releases] [extra sbatch args...]`.
#
# It exists for one reason: `--ntasks` and the release count have to agree, and
# an `#SBATCH` directive cannot see the job script's own arguments. The wrapper
# sets both from the same number, so a partial year (the last forced year, where
# only the releases whose full simulation window fits inside the available
# forcing may run) books exactly as many tasks as it has work. It also expands
# scripts/node-blacklist.txt into `--exclude=`, which `#SBATCH` cannot read from
# a file either.
#
#   scripts/submit_010.sh bottom 2022                  # full year, 73 tasks
#   scripts/submit_010.sh surface_stokes 2025 29       # partial year, 29 tasks
#   scripts/submit_010.sh bottom 2023 73 --dependency=afterok:12345
set -euo pipefail
cd "$(dirname "$0")/.."

regime=${1:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
year=${2:?usage: submit_010.sh <regime> <year> [n_releases] [sbatch args...]}
n_releases=${3:-73}
shift $(( $# < 3 ? $# : 3 ))

job=scripts/010_FucusDispersal_${regime}_job.sh
[ -f "${job}" ] || { echo "no such regime: ${regime} (${job})" >&2; exit 1; }
[ "${n_releases}" -ge 1 ] && [ "${n_releases}" -le 73 ] \
    || { echo "n_releases must be 1..73, got ${n_releases}" >&2; exit 1; }

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "regime=${regime} year=${year} n_releases=${n_releases} excluding: ${exclude:-<none>}"
exec sbatch --ntasks="${n_releases}" ${exclude:+--exclude="${exclude}"} "$@" \
    "${job}" "${year}" "${n_releases}"
