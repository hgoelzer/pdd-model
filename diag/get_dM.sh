#!/bin/bash
# Calculate mass change 2070-2099 vs 1960-1989
#cdo expr,'diff = var1 - var2' -merge file1.nc file2.nc output.nc

filename=scalars_mm_08.nc

# control time ranges
ncks -O -d time,110,139 ${filename} t0_tmp.nc
ncks -O -d time,220,249 ${filename} t1_tmp.nc
ncdump -t t0_tmp.nc | grep "time ="
ncdump -t t1_tmp.nc | grep "time ="

# calculate differences
ncra -O -d time,110,139 ${filename} t0_tmp.nc
ncra -O -d time,220,249 ${filename} t1_tmp.nc
ncdiff -O t1_tmp.nc t0_tmp.nc dm_tmp.nc 
ncap2 -O -s "dMgt = acabf*1e-12" dm_tmp.nc dm_tmp.nc 
echo "Mass change in Gt"
ncks -H -C -v dMgt -s "%f\n" dm_tmp.nc
# clean up
/bin/rm t0_tmp.nc t1_tmp.nc dm_tmp.nc

