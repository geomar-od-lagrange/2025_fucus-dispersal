#!/bin/bash
# Submit all visualisation jobs. 020/021/022 run once; 023/024 run per regime.
set -euo pipefail

cd "$(dirname "$0")/.."

regimes=(bottom surface surface_stokes)

sbatch --cpus-per-task=16 scripts/020_RawTrajectories_job.sh
sbatch --cpus-per-task=16 scripts/021_TimeStats_job.sh
sbatch --cpus-per-task=16 scripts/022_DispersalDistance_job.sh
for r in "${regimes[@]}"; do
    sbatch --cpus-per-task=16 scripts/024_HexHeatmaps_job.sh "$r"
done
