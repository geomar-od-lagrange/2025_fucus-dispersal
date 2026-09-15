#!/bin/bash
# Submit the three per-(regime, year) hex aggregations of one release year:
#
#   024   counts store        (Dask, reads the year's trajectory zarrs)
#   024b  distance histograms (Dask, reads the same zarrs)
#   024c  connectivity store  (reduces 024's counts -> --dependency=afterok)
#
# This is the thing to fire once per regime-year as the 010 runs land: 024 and
# 024b are independent and start together, 024c waits for 024 to finish clean.
# The beaching branch (024d and its reducers) is a separate, whole-regime pass
# -- see scripts/submit_024d.sh.
#
#   scripts/submit_aggregates.sh <regime> <year> [hex_radius] [sbatch args...]
#
#   scripts/submit_aggregates.sh bottom 2020
#   scripts/submit_aggregates.sh surface_stokes 2025 6000 --time=12:00:00
#
# 024a_BuildHexKey_job.sh must have run first for the matching hex_radius.
set -euo pipefail
cd "$(dirname "$0")/.."

regime=${1:?usage: submit_aggregates.sh <regime> <year> [hex_radius] [sbatch args...]}
year=${2:?usage: submit_aggregates.sh <regime> <year> [hex_radius] [sbatch args...]}
hex_radius=${3:-6000}
shift $(( $# < 3 ? $# : 3 ))

output_root=/gxfs_work/geomar/smomw122/2025_fucus_dispersal_outputs
trajdir=${output_root}/Trajectories/${regime}/${year}
if [ ! -d "${trajdir}" ]; then
    echo "no such trajectory directory: ${trajdir}" >&2
    exit 1
fi
# Under `set -e` a bare count of a year that has not been run yet would abort
# with no message, so the directory is checked before it is counted.
nzarr=$(find "${trajdir}" -maxdepth 1 -name '*.zarr' -type d | wc -l)
if [ "${nzarr}" -eq 0 ]; then
    echo "no trajectory zarrs under ${trajdir}" >&2
    exit 1
fi

exclude=$(sed 's/#.*//' scripts/node-blacklist.txt | tr -d '[:blank:]' \
          | grep -v '^$' | paste -sd,)
sb=(sbatch --parsable)
if [ -n "${exclude}" ]; then
    sb+=(--exclude="${exclude}")
fi
sb+=("$@")

echo "regime=${regime} year=${year} hex_radius=${hex_radius} zarrs=${nzarr}" \
     "excluding: ${exclude:-<none>}"

counts=$("${sb[@]}" scripts/024_BuildHexAggregates_job.sh "${regime}" "${year}" "${hex_radius}")
distance=$("${sb[@]}" scripts/024b_BuildHexDistance_job.sh "${regime}" "${year}" "${hex_radius}")
connectivity=$("${sb[@]}" --dependency=afterok:"${counts}" \
    scripts/024c_BuildHexConnectivity_job.sh "${regime}" "${year}" "${hex_radius}")

echo "024  counts       ${counts}"
echo "024b distance     ${distance}"
echo "024c connectivity ${connectivity}  (afterok:${counts})"
