#!/bin/bash
# Compress nc output

set -e
#set -x

#files=`ls output/ | grep ".nc"`
#files=`ls output/ | grep ".nc" | grep -v acabf`

#for pattern in prsn_GIS_NORCEPDD1_CESM2_Historical prra_GIS_NORCEPDD1_CESM2_SSP370; do
#for pattern in 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_206[89]' 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_20[789]' 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2100'; do

for pattern in 'tas_GIS_NORCEPDD1_MRI-ESM2-0_Historical_r1i1p1f1'  'snicefreez_GIS_NORCEPDD1_MRI-ESM2-0_Historical_r1i1p1f1'; do

#for pattern in 'prra_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2015.nc'    'prsn_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2015.nc'    'rfr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2015.nc'    'snicefreez_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2015.nc'    'tas_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2015.nc'; do
    
    files=`ls output/ | grep ".nc" | grep $pattern`

    for file in $files; do
	#echo $file
	./check_netcdf_compression.sh output/$file
    done

done
