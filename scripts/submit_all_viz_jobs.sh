#!/bin/bash
# Submit all visualisation jobs. 020/021/022 run once; 024 runs per regime.
set -euo pipefail

cd "$(dirname "$0")/.."

# Nodes excluded from allocation (ib0 absent or flaky).
EXCLUDE="nesh-clk[385-387,415,430,446,483,557,564,573]"

sbatch --exclude="$EXCLUDE" --ntasks=7 scripts/020_RawTrajectories_job.sh
sbatch --exclude="$EXCLUDE" --ntasks=7 scripts/021_TimeStats_job.sh
sbatch --exclude="$EXCLUDE" --ntasks=7 scripts/022_DispersalDistance_job.sh

sbatch --exclude="$EXCLUDE" --ntasks=9 --cpus-per-task=10 scripts/024_HexHeatmaps_job.sh surface
sbatch --exclude="$EXCLUDE" --ntasks=9 --cpus-per-task=10 scripts/024_HexHeatmaps_job.sh surface_stokes
sbatch --exclude="$EXCLUDE" --ntasks=19 --cpus-per-task=10 scripts/024_HexHeatmaps_job.sh bottom
