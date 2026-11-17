#!/bin/bash
#SBATCH --job-name=010_FucusBottomUb49
#SBATCH --ntasks=73
#SBATCH --cpus-per-task=1
# The bottom job's sizing: friction only shortens the drift, so the zarr and
# the in-memory output store at the final write stay below the bottom ones.
#SBATCH --mem-per-cpu=15G
#SBATCH --time=36:00:00
#SBATCH --partition=base
#SBATCH --no-kill
#SBATCH --output=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs/logs/010/%x_%j.out

# The `bottom_ub49` regime: the bottom run with bed friction at the lab
# break-off speed U_b = 49 mm/s (docs/bottom_friction.md), written under
# Trajectories/bottom_ub49/. Same arguments and environment as
# scripts/010_FucusDispersal_bottom_job.sh, which does the work; this script
# only pins U_b. One script per output regime means submit_010.sh's
# regime -> script lookup, its job names and submit_downstream_chain.sh's
# regime discovery all key on the script name, with no regime-name parsing,
# while the run logic stays in the bottom script. The #SBATCH sizing above is
# a copy of the bottom script's and does not follow it.
#
#   scripts/submit_010.sh bottom_ub49 2019
#   scripts/submit_010.sh bottom_ub49 2025 29

# Run from the submit dir (the repo root, as submit_010.sh submits it); the
# #SBATCH lines above are this job's, the bottom script's are comments here.
exec bash scripts/010_FucusDispersal_bottom_job.sh "${1:-2019}" "${2:-73}" 49
