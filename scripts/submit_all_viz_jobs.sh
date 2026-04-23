#!/bin/bash
# Submit all visualisation jobs. 020/021/022 run once; 024 runs per regime.
set -euo pipefail

cd "$(dirname "$0")/.."

regimes=(bottom surface surface_stokes)

sbatch scripts/020_RawTrajectories_job.sh
sbatch scripts/021_TimeStats_job.sh
sbatch scripts/022_DispersalDistance_job.sh
for r in "${regimes[@]}"; do
    sbatch scripts/024_HexHeatmaps_job.sh "$r"
done
