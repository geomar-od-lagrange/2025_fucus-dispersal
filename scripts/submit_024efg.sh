#!/bin/bash
# Submit 024efg_ReduceAll_job.sh sized to the work it will actually do.
#
# Same two jobs as scripts/submit_024d.sh: expand scripts/node-blacklist.txt
# into `--exclude=` (an `#SBATCH` directive cannot read a file) and set
# `--ntasks` from the job's own work list (`LIST_ONLY=1`, one line per
# (reducer, member, year, month) cell), capped at $MAX_TASKS. The cap matters
# more here than for 024d: the cell count is the member table times the
# (year, month) grid, thousands of one-to-two-minute cells that a
# `xargs -P ${SLURM_NTASKS}` fan-out drains in rounds. Extra sbatch args are
# passed through after the computed ones, so an explicit --ntasks=N wins.
#
#   scripts/submit_024efg.sh --dependency=afterok:23835103
#   YEARS="2020 2021" scripts/submit_024efg.sh          # restrict the year set
#   MAX_TASKS=48 scripts/submit_024efg.sh               # throttle
set -euo pipefail
cd "$(dirname "$0")/.."

job=scripts/024efg_ReduceAll_job.sh
max_tasks=${MAX_TASKS:-292}

njobs=$(LIST_ONLY=1 bash "${job}" | wc -l)
[ "${njobs}" -gt 0 ] || { echo "nothing to do: empty work list" >&2; exit 1; }
ntasks=$(( njobs < max_tasks ? njobs : max_tasks ))

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "${njobs} cells, --ntasks=${ntasks}, excluding: ${exclude:-<none>}"
exec sbatch --ntasks="${ntasks}" ${exclude:+--exclude="${exclude}"} "$@" "${job}"
