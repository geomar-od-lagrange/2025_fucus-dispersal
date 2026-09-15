#!/bin/bash
# Submit 024d_BuildBeachingForcing_job.sh sized to the work it will actually do.
#
# Two things an `#SBATCH` directive cannot do: read scripts/node-blacklist.txt
# into `--exclude=`, and count the trajectory zarrs on disk. Both happen here.
# `--ntasks` is the job's own work list (`LIST_ONLY=1`, one line per zarr),
# capped at $MAX_TASKS so a grown year set does not turn into an allocation that
# pends forever -- `xargs -P ${SLURM_NTASKS}` inside the job decouples
# concurrency from the work count, so the cap costs rounds, not correctness.
# Extra sbatch args are passed through after the computed ones, so an explicit
# --ntasks=N wins.
#
#   scripts/submit_024d.sh
#   YEARS="2020 2021" scripts/submit_024d.sh            # restrict the year set
#   MAX_TASKS=146 scripts/submit_024d.sh                # throttle
#   ONLY_STEMS_FILE=/path/stems.txt scripts/submit_024d.sh   # retry a subset
set -euo pipefail
cd "$(dirname "$0")/.."

job=scripts/024d_BuildBeachingForcing_job.sh
max_tasks=${MAX_TASKS:-292}

njobs=$(LIST_ONLY=1 bash "${job}" | wc -l)
[ "${njobs}" -gt 0 ] || { echo "nothing to do: empty work list" >&2; exit 1; }
ntasks=$(( njobs < max_tasks ? njobs : max_tasks ))

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "${njobs} zarrs, --ntasks=${ntasks}, excluding: ${exclude:-<none>}"
exec sbatch --ntasks="${ntasks}" ${exclude:+--exclude="${exclude}"} "$@" "${job}"
