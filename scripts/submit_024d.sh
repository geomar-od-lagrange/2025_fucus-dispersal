#!/bin/bash
# Submit 024d_BuildBeachingForcing_job.sh with scripts/node-blacklist.txt kept
# out of the allocation. `#SBATCH --exclude=` cannot read a file, so the list is
# expanded here; extra sbatch args (e.g. --ntasks=239) are passed through.
#
#   scripts/submit_024d.sh
#   ONLY_STEMS_FILE=/path/stems.txt scripts/submit_024d.sh --ntasks=239
set -euo pipefail
cd "$(dirname "$0")/.."
exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
echo "excluding: ${exclude:-<none>}"
exec sbatch ${exclude:+--exclude="${exclude}"} "$@" \
    scripts/024d_BuildBeachingForcing_job.sh
