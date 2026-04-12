# Produce some quick diags

# switch from 0 = mm s-1 to 1 = mm yr-1
outputmode=0
ares=08
#ares=01

run=output_${ares}000m_baseline

# time series
ncrcat -O ../${run}/acabf_GIS_NORCE-PDD* acabf.nc
ncrcat -O ../${run}/pr_GIS_NORCE-PDD* pr.nc
#ncrcat -O ../${run}/tas_GIS_NORCE-PDD* tas.nc
#ncrcat -O ../${run}/prsn_GIS_NORCE-PDD* prsn.nc
#ncrcat -O ../${run}/mrro_GIS_NORCE-PDD* mrro.nc

module load CDO
cdo yearmean acabf.nc acabf_yearly.nc
if (( outputmode == 0 )); then
    ncap2 -O -s "acabf=acabf*31556926" acabf_yearly.nc acabf_yearly.nc
    ncatted -h -a units,acabf,o,c,"mm w.e. yr-1" acabf_yearly.nc
fi

cdo yearmean pr.nc pr_yearly.nc
if (( outputmode == 0 )); then
    ncap2 -O -s "pr=pr*31556926" pr_yearly.nc pr_yearly.nc
    ncatted -h -a units,pr,o,c,"mm w.e. yr-1" pr_yearly.nc
fi

#cdo yearmean tas.nc tas_yearly.nc
#
#cdo yearmean prsn.nc prsn_yearly.nc
#if (( outputmode == 0 )); then
#    ncap2 -O -s "prsn=prsn*31556926" prsn_yearly.nc prsn_yearly.nc
#    ncatted -h -a units,pr,o,c,"mm w.e. yr-1" prsn_yearly.nc
#fi
#
#cdo yearmean mrro.nc mrro_yearly.nc
#if (( outputmode == 0 )); then
#    ncap2 -O -s "mrro=mrro*31556926" mrro_yearly.nc mrro_yearly.nc
#    ncatted -h -a units,mrro,o,c,"mm w.e. yr-1" mrro_yearly.nc
#fi


## differences
#ncdiff -O tas.nc tas0.nc dtas.nc
#ncdiff -O pr.nc pr0.nc dpr.nc
#ncdiff -O acabf.nc acabf0.nc dacabf.nc

# output for scalar calculation
ncks -O -v acabf acabf_yearly.nc ./smb_gpdd_${ares}000m.nc
ncks -A -v pr pr_yearly.nc ./smb_gpdd_${ares}000m.nc
ncap2 -O -s "acabf=float(acabf)" ./smb_gpdd_${ares}000m.nc ./smb_gpdd_${ares}000m.nc
ncap2 -O -s "pr=float(pr)" ./smb_gpdd_${ares}000m.nc ./smb_gpdd_${ares}000m.nc
ncap2 -O -s "x=float(x); y=float(y)" ./smb_gpdd_${ares}000m.nc ./smb_gpdd_${ares}000m.nc

# clean up
/bin/rm acabf.nc acabf_yearly.nc pr.nc pr_yearly.nc
