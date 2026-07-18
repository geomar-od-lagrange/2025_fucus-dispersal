#!/bin/bash
#SBATCH --job-name=024d_BuildBeaching
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=01:00:00
#SBATCH --partition=base

# Post-simulation beaching pass: reads the trajectory zarrs + raw
# baltic_highres Stokes + the 024a key, writes one beaching parquet per
# (regime, release_year, release_month). A single-process numpy notebook
# (no Dask) — each zarr fits in memory and the bottleneck is Stokes I/O, so
# one papermill run processes a month's ~6 zarrs sequentially (~7 min).
#
# Parallelism fans out the full (year × month) grid with xargs dispatching
# one `srun --ntasks=1 --exact` step per cell, like 010_FucusDispersal_*.
# njobs is fixed (|YEARS|·12 = 48) but concurrency is whatever --ntasks the
# scheduler grants (`xargs -P ${SLURM_NTASKS}`), so it rescales: request
# fewer tasks under load without editing the job. The store is partitioned
# per (regime, year, month) — each cell writes its own `_mMM` file, no shared
# output, no merge; 029 pools the monthly partitions.
#
# Usage: sbatch scripts/024d_BuildBeaching_job.sh [regime] [hex_radius]
#   sbatch scripts/024d_BuildBeaching_job.sh                     # surface_stokes, all cells
#   sbatch --ntasks=8 scripts/024d_BuildBeaching_job.sh          # throttle concurrency
#   sbatch scripts/024d_BuildBeaching_job.sh surface 6000        # surface regime
# 024a_BuildHexKey_job.sh must have run first for the matching hex_radius.

YEARS=(2016 2017 2018 2019)

# Drop any CPU-bind mask inherited from an outer allocation (present when
# this is submitted from inside an interactive job); otherwise the
# concurrent srun steps below try to bind to the outer job's CPUs and fail
# with "CPU binding outside of job step allocation". --exact sets each
# step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

regime="${1:-surface_stokes}"
hex_radius="${2:-6000}"

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
export output_root regime hex_radius

mkdir -p notebooks_executed/Visualisations/

# Emit every "year month" pair; xargs runs up to ${SLURM_NTASKS} at once,
# each dispatching one srun --ntasks=1 --exact papermill step.
for year in "${YEARS[@]}"; do
    for month in $(seq 1 12); do
        printf '%s\0' "${year} ${month}"
    done
done | xargs -0 -P "${SLURM_NTASKS}" -n 1 bash -c '
    read -r year month <<< "$1"
    ms=$(printf "_m%02d" "${month}")
    srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --exact \
        pixi run papermill --cwd notebooks/ \
        notebooks/024d_BuildBeaching.ipynb \
        notebooks_executed/Visualisations/024d_BuildBeaching_${regime}_${year}${ms}_r${hex_radius}m.ipynb \
        -p output_root ${output_root} \
        -p regime ${regime} \
        -p release_year ${year} \
        -p release_month ${month} \
        -p hex_radius ${hex_radius} \
        -k python
' _
rc=$?

jobinfo
exit ${rc}
