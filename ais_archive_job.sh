#!/bin/bash
#SBATCH --job-name=ais_archive
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#
# Archive the corrected AIS SMBMIP rerun (job 1807220 output, verified by
# 1807527) to NIRD. New directory; the original taj-era archive
# output_08000m/ is untouched.
#   src: /cluster/work/users/heig/pdd/ais_08/output/  (1359 files, 9 vars)
#   dst: /nird/datalake/NS8002K/heig/PDD-AIS/Output/ais_08/output_08000m_rerun2026-11/
# Note: NS8002K is mounted on compute nodes (checked: the archive-compare
# failure in 1807527 was the flat staging, not the mount; nird is reachable
# from b5201). If the mount is missing on the node, the job fails fast at
# the write test.
# Usage: sbatch ais_archive_job.sh

set -o errexit
set -o nounset

SRC=/cluster/work/users/heig/pdd/ais_08/output
DST=/nird/datalake/NS8002K/heig/PDD-AIS/Output/ais_08/output_08000m_rerun2026-11

echo "== 1. Source check =="
NSRC=$(ls "$SRC" | wc -l)
echo "  source files: $NSRC (expect 1359)"
[ "$NSRC" -eq 1359 ] || { echo "ERROR: unexpected source file count"; exit 1; }

echo "== 2. Destination write test =="
mkdir -p "$DST"
touch "$DST/.wtest" && rm "$DST/.wtest" || { echo "ERROR: cannot write to $DST"; exit 1; }

echo "== 3. rsync =="
rsync -a --info=progress2 "$SRC"/ "$DST"/ | tail -2

echo "== 4. Verify =="
NDST=$(ls "$DST" | wc -l)
echo "  destination files: $NDST (expect 1359)"
[ "$NDST" -eq 1359 ] || { echo "ERROR: destination file count mismatch"; exit 1; }
du -sh "$DST"
cmp "$SRC/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc" \
    "$DST/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc" \
    && echo "  spot-check acabf 1950: identical" \
    || { echo "ERROR: spot-check failed"; exit 1; }
cmp "$SRC/mrro_AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1_2100.nc" \
    "$DST/mrro_AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1_2100.nc" \
    && echo "  spot-check mrro 2100: identical" \
    || { echo "ERROR: spot-check failed"; exit 1; }

echo "== DONE: rerun archived to $DST =="