#!/bin/bash
# Copy a minimal set of BSH operational model files for testing.
# Run from the directory where you want the data to land.
# Copies first 4 timesteps (1 day) of 2020 for all file types.

SRC_ORIG=/gxfs_work/geomar/smomw400/bsh_operationalmodel_data
SRC_STATIC=/gxfs_work/geomar/smomw122/bsh_operationalmodel_data

# Current files (U, V) — first 4 files = first day
for res in fine coarse; do
    mkdir -p bsh_operationalmodel_data/c_file_${res}_2020
    for f in $(ls ${SRC_ORIG}/c_file_${res}_2020/c_file_${res}_202001010* 2>/dev/null); do
        cp -v "$f" bsh_operationalmodel_data/c_file_${res}_2020/
    done
done

# Salt/temp files — first 4 files = first day
for res in fine coarse; do
    mkdir -p bsh_operationalmodel_data/t_file_${res}_2020
    for f in $(ls ${SRC_ORIG}/t_file_${res}_2020/t_file_${res}_202001010* 2>/dev/null); do
        cp -v "$f" bsh_operationalmodel_data/t_file_${res}_2020/
    done
done

# Static files (lonlat, H0)
for res in fine coarse; do
    mkdir -p bsh_operationalmodel_data/static_file_${res}
    cp -v ${SRC_STATIC}/static_file_${res}/lonlat_file_${res}.nc bsh_operationalmodel_data/static_file_${res}/
    cp -v ${SRC_STATIC}/static_file_${res}/H0_file_${res}.nc bsh_operationalmodel_data/static_file_${res}/
done

echo ""
echo "Done. Copied to: $(pwd)/bsh_operationalmodel_data/"
echo "To use: set path_orig_files and path_static_files to $(pwd)/bsh_operationalmodel_data/"
