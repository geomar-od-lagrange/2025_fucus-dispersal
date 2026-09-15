#!/bin/bash
# Submit 024efg_ReduceAll_job.sh with scripts/node-blacklist.txt kept out of
# the allocation. `#SBATCH --exclude=` cannot read a file, so the list is
# expanded here; extra sbatch args are passed through.
#
#   scripts/submit_024efg.sh --ntasks=292 --dependency=afterok:23835103
set -euo pipefail
cd "$(dirname "$0")/.."
exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "excluding: ${exclude:-<none>}"
exec sbatch ${exclude:+--exclude="${exclude}"} "$@" \
    scripts/024efg_ReduceAll_job.sh
