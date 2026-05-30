#!/bin/bash
# Produce some quick diags

#set -x

# switch from 0 = mm s-1 to 1 = mm yr-1
outputmode=0
#ares=08
ares=01
run=CESM2_r11i1p1f1
#run=MRI-ESM2-0_r1i1p1f1

mkdir -p ${run}

# time series
echo "make time series"
ncrcat -O ../output/${run}/acabf_GIS_NORCEPDD1* ${run}/acabf.nc
ncrcat -O ../output/${run}/pr_GIS_NORCEPDD1* ${run}/pr.nc
#ncrcat -O ../output/${run}/tas_GIS_NORCEPDD1* ${run}/tas.nc
#ncrcat -O ../output/${run}/prsn_GIS_NORCEPDD1* ${run}/prsn.nc
#ncrcat -O ../output/${run}/mrro_GIS_NORCEPDD1* ${run}/mrro.nc

echo "monthly to yearly"
cdo yearmean ${run}/acabf.nc ${run}/acabf_yearly.nc
if (( outputmode == 0 )); then
    ncap2 -O -s "acabf=acabf*31556926" ${run}/acabf_yearly.nc ${run}/acabf_yearly.nc
    ncatted -h -a units,acabf,o,c,"mm w.e. yr-1" ${run}/acabf_yearly.nc
fi

cdo yearmean ${run}/pr.nc ${run}/pr_yearly.nc
if (( outputmode == 0 )); then
    ncap2 -O -s "pr=pr*31556926" ${run}/pr_yearly.nc ${run}/pr_yearly.nc
    ncatted -h -a units,pr,o,c,"mm w.e. yr-1" ${run}/pr_yearly.nc
fi

#cdo yearmean ${run}/tas.nc ${run}/tas_yearly.nc
#
#cdo yearmean ${run}/prsn.nc ${run}/prsn_yearly.nc
#if (( outputmode == 0 )); then
#    ncap2 -O -s "prsn=prsn*31556926" ${run}/prsn_yearly.nc ${run}/prsn_yearly.nc
#    ncatted -h -a units,pr,o,c,"mm w.e. yr-1" ${run}/prsn_yearly.nc
#fi
#
#cdo yearmean ${run}/mrro.nc ${run}/mrro_yearly.nc
#if (( outputmode == 0 )); then
#    ncap2 -O -s "mrro=mrro*31556926" ${run}/mrro_yearly.nc ${run}/mrro_yearly.nc
#    ncatted -h -a units,mrro,o,c,"mm w.e. yr-1" ${run}/mrro_yearly.nc
#fi


## differences
#ncdiff -O tas.nc tas0.nc dtas.nc
#ncdiff -O pr.nc pr0.nc dpr.nc
#ncdiff -O acabf.nc acabf0.nc dacabf.nc

outfile=${run}/smb_gpdd_${ares}000m.nc
# output for scalar calculation
echo "produce summary file"
ncks -O -v acabf ${run}/acabf_yearly.nc ${outfile}
ncks -A -v pr ${run}/pr_yearly.nc ${outfile}
ncap2 -O -s "acabf=float(acabf)" ${outfile} ${outfile}
ncap2 -O -s "pr=float(pr)" ${outfile} ${outfile}
ncap2 -O -s "x=float(x); y=float(y)" ${outfile} ${outfile}

# clean up
#/bin/rm ${run}/acabf.nc ${run}/acabf_yearly.nc ${run}/pr.nc ${run}/pr_yearly.nc
