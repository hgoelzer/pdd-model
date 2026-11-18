#!/bin/bash
# ISMIP7 GrIS 1 km staging - run INTERACTIVELY on a bet LOGIN node.
# Transfers the 1 km CESM2 forcing (~10-15 GB total) and the archive
# comparison files from the NIRD datalake mount to bet work storage, plus
# the 1 km MAR reference and promice mask files. Run once before
# sbatch-ing the 1 km rerun/verify jobs (compute nodes cannot see the
# NIRD mount, so this cannot go into a job).
# Usage: bash gris_stage1km.sh
set -o errexit

N=/nird/datalake/NS8002K/heig/PDD
W=/cluster/work/users/heig/pdd/gris_01
F=$W/Forcing
mkdir -p $F/CESM2/historical/pr_ratio $F/CESM2/historical/tas_anom \
         $F/CESM2/ssp370/pr_ratio $F/CESM2/ssp370/tas_anom \
         $F/MARv3.14LC26 $F/constant \
         $W/archive_cmp/Historical $W/archive_cmp/SSP370

echo "== 1. CESM2 historical 1 km per-year forcing (1950-2014) =="
for Y in $(seq 1950 2014); do
    cp -n $N/CMIP/ismip6/CESM2/historical/pr_ratio/pr_ratio_Amon_CESM2_historical_r11i1p1f1_${Y}_i01000m.nc \
          $F/CESM2/historical/pr_ratio/
    cp -n $N/CMIP/ismip6/CESM2/historical/tas_anom/tas_Amon_CESM2_historical_r11i1p1f1_${Y}_i01000m.nc \
          $F/CESM2/historical/tas_anom/
done

echo "== 2. CESM2 ssp370 1 km per-year forcing (2015-2100) =="
cp -n $N/CMIP/ismip6/CESM2/ssp370/pr_ratio/pr_ratio_Amon_CESM2_ssp370_r11i1p1f1_????_i01000m.nc \
      $F/CESM2/ssp370/pr_ratio/
cp -n $N/CMIP/ismip6/CESM2/ssp370/tas_anom/tas_Amon_CESM2_ssp370_r11i1p1f1_????_i01000m.nc \
      $F/CESM2/ssp370/tas_anom/

echo "== 3. MAR references + promice mask (1 km) =="
cp -n $N/MAR/Lindsey-Clark2026/ismip6/pr_ref_MARv3.14LC26-monthly-ERA5-1960-1989_01000m.nc $F/MARv3.14LC26/
cp -n $N/MAR/MARv3.14/ismip6/tas_ref_MARv3.14-monthly-ERA5-1960-1989_01000m.nc $F/MARv3.14LC26/
cp -n $N/constant/sftgif_promice_01000m.nc $F/constant/

echo "== 4. Archive comparison files (1950/2014/2015/2100, 6 vars) =="
S=$N/Output/gris_01/CESM2_SSP370_r11i1p1f1
for V in acabf mrro pr prra prsn tas snicefreez; do
    for Y in 1950 2014; do
        cp -n $S/${V}_GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_${Y}.nc \
              $W/archive_cmp/Historical/
    done
    for Y in 2015 2100; do
        cp -n $S/${V}_GIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1_${Y}.nc \
              $W/archive_cmp/SSP370/
    done
done

echo "== Counts (expected 65/65/86/86 forcing, 12+12 cmp) =="
ls $F/CESM2/historical/pr_ratio | wc -l
ls $F/CESM2/historical/tas_anom | wc -l
ls $F/CESM2/ssp370/pr_ratio | wc -l
ls $F/CESM2/ssp370/tas_anom | wc -l
ls $W/archive_cmp/Historical | wc -l
ls $W/archive_cmp/SSP370 | wc -l
echo "== STAGING DONE =="
