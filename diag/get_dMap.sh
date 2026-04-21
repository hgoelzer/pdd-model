#!/bin/bash
# Calculate difference 2070-2099 vs 1960-1989

ares=08

#run=output_MAR
run=output

filename=${run}/acabf_yearly.nc

dmap=${run}/dMAP.nc

# control time ranges
ncks -O -d time,10,39 ${filename} t0_tmp.nc
ncks -O -d time,120,149 ${filename} t1_tmp.nc
ncdump -t t0_tmp.nc | grep "time ="
ncdump -t t1_tmp.nc | grep "time ="

# calculate differences
ncra -O -d time,10,39 ${filename} t0_tmp.nc
ncra -O -d time,120,149 ${filename} t1_tmp.nc
ncdiff -O t1_tmp.nc t0_tmp.nc ${dmap} 

#ncks -A -v sftgif ../../data/sftgif_BM5_${ares}000m.nc ${dmap}
#ncap2 -O -s "masked_acabf=acabf; where(sftgif<0.5)masked_acabf=0" ${dmap} ${dmap}
ncks -A -v icemask_promice ../../data/icemask_promice_${ares}000m.nc ${dmap}
ncap2 -O -s "masked_acabf=acabf; where(icemask_promice<0.5)masked_acabf=0" ${dmap} ${dmap}

# calculate dMass in Gt
ncap2 -O -s 'dMgt=masked_acabf.total($x,$y)*8000.^2*1e-12' ${dmap} ${dmap}

# clean up
#/bin/rm t0_tmp.nc t1_tmp.nc dm_tmp.nc

