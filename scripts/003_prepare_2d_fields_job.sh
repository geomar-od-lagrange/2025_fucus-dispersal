#!/bin/bash
#SBATCH --job-name=prep_2d_fields
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

repo_root=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
bsh_root=/gxfs_work/geomar/smomw400/bsh_operationalmodel_data
container=parcels-container_2024.10.07-7af7fd0.sif

mkdir -p ${output_root}/2d_fields

# Process one c_file: called by xargs with the c_file path as argument
process_file() {
    local c_file="$1"
    srun --ntasks=1 --exact \
        singularity run -B /sfs -B /gxfs_work -B ${repo_root}:/work --pwd /work \
        ${container} bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
        && python notebooks/003_prepare_2d_fields.py \
            --c-file ${c_file} \
            --output-root ${output_root}"
}
export -f process_file
export repo_root output_root container

# Process all years and resolutions in a single xargs pass.
# Sorted find output means early years/timestamps come first.
for year in $(seq 2016 2025); do
    for res in fine coarse; do
        c_dir="${bsh_root}/c_file_${res}_${year}"
        [ -d "${c_dir}" ] && find "${c_dir}" -name "*.nc" | sort
    done
done | xargs -P ${SLURM_NTASKS} -I{} bash -c 'process_file "$@"' _ {}

jobinfo
