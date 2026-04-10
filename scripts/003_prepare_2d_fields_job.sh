#!/bin/bash
#SBATCH --job-name=prep_2d_fields
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --partition=base

module load gcc12-env/12.3.0
module load singularity/3.11.5

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal
bsh_data=/gxfs_work/geomar/smomw400/bsh_operationalmodel_data
stokes_dir=/gxfs_work/geomar/smomw122/cmems_mod_bal_wav
output_dir=${base_path}/output/2d_fields
container=parcels-container_2024.10.07-7af7fd0.sif

mkdir -p ${output_dir}

# Process one c_file: called by xargs with the c_file path as argument
process_file() {
    local c_file="$1"
    srun --ntasks=1 --exact \
        singularity run -B /sfs -B /gxfs_work -B ${base_path}:/work --pwd /work \
        ${container} bash -c \
        ". /opt/conda/etc/profile.d/conda.sh && conda activate base \
        && python scripts/003_prepare_2d_fields.py \
            --c-file ${c_file} \
            --stokes-dir ${stokes_dir} \
            --output-dir ${output_dir}"
}
export -f process_file
export base_path container stokes_dir output_dir

# Loop over years (outer) so each year completes with both resolutions
# before moving on. This lets 010 experiments start as soon as a year is done.
for year in $(seq 2016 2025); do
    for res in fine coarse; do
        c_dir="${bsh_data}/c_file_${res}_${year}"
        if [ ! -d "${c_dir}" ]; then
            echo "Skipping ${c_dir} (not found)"
            continue
        fi
        echo "Processing ${c_dir} ..."
        find "${c_dir}" -name "*.nc" | sort | \
            xargs -P ${SLURM_NTASKS} -I{} bash -c 'process_file "$@"' _ {}
    done
done

jobinfo
