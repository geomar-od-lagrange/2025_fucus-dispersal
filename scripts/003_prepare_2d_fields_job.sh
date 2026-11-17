#!/bin/bash
#SBATCH --job-name=prep_2d_fields
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --partition=base

# Which c-files get processed is steered by three optional environment
# variables (`sbatch --export=ALL,YEARS=2025,MISSING_ONLY=1 ...`):
#
#   YEARS="2021 2025"        year loop (default: 2016 … 2025)
#   CFILE_LIST=<path>        process exactly the c-file paths listed in this
#                            file, one per line; YEARS is then unused
#   MISSING_ONLY=1           drop every c-file whose three outputs all exist
#
# 003 itself already skips a c-file whose outputs are all there, but it pays
# a container start-up and an `srun` step per skip. MISSING_ONLY does the
# same test up front from a single directory listing, so refilling a handful
# of holes costs a handful of steps instead of 29 214.

module load gcc12-env/12.3.0
module load singularity/3.11.5

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
bsh_root=/gxfs_work/geomar/smomw400/bsh_operationalmodel_data
container=parcels-container_2024.10.07-7af7fd0.sif

mkdir -p ${output_root}/2d_fields

# Drop any CPU-bind mask inherited from an outer allocation; otherwise the
# concurrent srun steps below fail with "CPU binding outside of job step
# allocation". --exact sets each step's own binding.
unset SLURM_CPU_BIND SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND_VERBOSE

# Process one c_file: called by xargs with the c_file path as argument
process_file() {
    local c_file="$1"
    # -N1: without it a step inherits the job's node count and dies with
    # "can't run 1 processes on N nodes" as soon as the allocation spans
    # more than one node.
    srun -n1 -N1 --exact --no-kill \
        singularity run -B /sfs -B /gxfs_work -B ${repo_root}:/work --pwd /work \
        ${container} bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
        && python notebooks/003_prepare_2d_fields.py \
            --c-file ${c_file} \
            --output-root ${output_root}"
}
export -f process_file
export repo_root output_root container

# Sorted output means early years/timestamps come first.
enumerate_cfiles() {
    if [ -n "${CFILE_LIST:-}" ]; then
        grep -v '^[[:space:]]*$' "${CFILE_LIST}"
        return
    fi
    for year in ${YEARS:-$(seq 2016 2025)}; do
        for res in fine coarse; do
            c_dir="${bsh_root}/c_file_${res}_${year}"
            [ -d "${c_dir}" ] && find "${c_dir}" -name "*.nc" | sort
        done
    done
}

# One directory listing of 2d_fields instead of three stat() calls per
# candidate — kinder to the shared parallel filesystem, and the listing is
# the same cost whether one year or ten is enumerated.
keep_incomplete() {
    awk 'NR == FNR { have[$0] = 1; next }
         {
             stem = $0
             sub(/.*\//, "", stem)
             sub(/\.nc$/, "", stem)
             if (!((stem "_surface.nc") in have) ||
                 !((stem "_bottom.nc") in have) ||
                 !((stem "_surface_stokes.nc") in have)) print
         }' <(ls "${output_root}/2d_fields") -
}

filter=cat
if [ -n "${MISSING_ONLY:-}" ]; then
    filter=keep_incomplete
fi

enumerate_cfiles \
    | ${filter} \
    | xargs -P ${SLURM_NTASKS} -I{} bash -c 'process_file "$@"' _ {}

jobinfo
