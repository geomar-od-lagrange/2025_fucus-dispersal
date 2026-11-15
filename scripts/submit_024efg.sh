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
#   scripts/submit_024efg.sh <regime> [sbatch args...]
#   scripts/submit_024efg.sh bottom --dependency=afterok:23835103
#   YEARS="2020 2021" scripts/submit_024efg.sh bottom   # restrict the year set
#   MAX_TASKS=48 scripts/submit_024efg.sh bottom        # throttle
set -euo pipefail
cd "$(dirname "$0")/.."
source "$(dirname "$0")/slurm_submit.sh"

job=scripts/024efg_ReduceAll_job.sh
max_tasks=${MAX_TASKS:-292}
regime=${1:?usage: scripts/submit_024efg.sh <regime> [sbatch args...]}
shift

njobs=$(LIST_ONLY=1 bash "${job}" "${regime}" | wc -l)
# The job re-lists its cells at run time. An empty list now means the
# sidecars are not written yet -- this job queued behind the 024d that writes
# them (submit_downstream_chain.sh) -- so book the cap rather than refuse.
if [ "${njobs}" -eq 0 ]; then
    echo "${regime}: no sidecars on disk yet, sizing at the cap" >&2
    njobs=${max_tasks}
fi
ntasks=$(( njobs < max_tasks ? njobs : max_tasks ))

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "${regime}: ${njobs} cells, --ntasks=${ntasks}, excluding: ${exclude:-<none>}"
# robust_sbatch, not sbatch: it retries a controller that times out and adopts
# a timed-out-but-created job by --job-name instead of queueing a second copy.
# Its last line is "<job id> <adopted|submitted>" rather than sbatch's own.
robust_sbatch --ntasks="${ntasks}" ${exclude:+--exclude="${exclude}"} "$@" "${job}" "${regime}"
