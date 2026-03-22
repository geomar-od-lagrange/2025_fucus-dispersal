#!/bin/bash

# make sure we have Singularity
module load gcc12-env/12.3.0
module load singularity/3.11.5

# to get the image (need to be on a partition which has internet access --> data), run
# $ singularity pull --disable-cache --dir "${PWD}" docker://quay.io/willirath/parcels-container:2024.10.07-7af7fd0

base_path=/gxfs_work/geomar/smomw122/2025_fucus-dispersal/
# make sure the output dir exists
mkdir -p notebooks_executed/
singularity run -B /sfs -B /gxfs_work -B $PWD:/work --pwd /work parcels-container_2024.10.07-7af7fd0.sif bash -c \
". /opt/conda/etc/profile.d/conda.sh && conda activate base \
&& papermill --cwd notebooks/ \
    notebooks/000_FucusStartLocations.ipynb \
    notebooks_executed/000_FucusStartLocations.ipynb \
    -p base_path ${base_path} \
    -k python"
