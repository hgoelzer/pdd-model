#!/bin/bash
# Calculate difference 2070-2099 vs 1960-1989

filename=acabf_yearly.nc

# control time ranges
ncks -O -d time,10,39 ${filename} t0_tmp.nc
ncks -O -d time,120,149 ${filename} t1_tmp.nc
ncdump -t t0_tmp.nc | grep "time ="
ncdump -t t1_tmp.nc | grep "time ="

# calculate differences
ncra -O -d time,10,39 ${filename} t0_tmp.nc
ncra -O -d time,120,149 ${filename} t1_tmp.nc
ncdiff -O t1_tmp.nc t0_tmp.nc dMAP.nc 

ncks -A -v icemask_promice ../../data/icemask_promice_08000m.nc dMAP.nc
ncap2 -O -s "masked_acabf=acabf; where(icemask_promice<0.5)masked_acabf=0" dMAP.nc dMAP.nc

# calculate dMass in Gt
ncap2 -O -s 'dMgt=masked_acabf.total($x,$y)*8000.^2*1e-12' dMAP.nc dMAP.nc

# clean up
#/bin/rm t0_tmp.nc t1_tmp.nc dm_tmp.nc

