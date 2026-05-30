#!/bin/bash
# Compress nc output

set -e
set -x

#files=`ls output/ | grep ".nc"`
#files=`ls output/ | grep ".nc" | grep -v acabf`

#for pattern in prsn_GIS_NORCEPDD1_CESM2_Historical prra_GIS_NORCEPDD1_CESM2_SSP370; do
#for pattern in 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_206[89]' 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_20[789]' 'pr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1_2100'; do

#for pattern in 'tas_GIS_NORCEPDD1_MRI-ESM2-0_Historical_r1i1p1f1'  'snicefreez_GIS_NORCEPDD1_MRI-ESM2-0_Historical_r1i1p1f1'; do

for pattern in 'prra_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1'    'prsn_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1'    'rfr_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1'    'snicefreez_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1'    'tas_GIS_NORCEPDD1_MRI-ESM2-0_SSP585_r1i1p1f1'; do
    
    files=`ls output/ | grep ".nc" | grep $pattern`

    for file in $files; do
	echo $file
	
	# Use ncdump -hs to show storage format (includes compression info)
	ncdump_output=$(ncdump -hs "output/$file" 2>/dev/null)

	# Check for compression indicators in storage info
	if echo "$ncdump_output" | grep -qi "deflate\|szip\|chunked" >/dev/null 2>&1; then
            echo "  ✓ COMPRESSED"
	else
	    tmpfile="tmp_$(date +%s)_$RANDOM.nc"  # timestamp + random
	    nccopy -d9 output/$file $tmpfile
	    mv $tmpfile output/$file
	fi
	
    done

done
