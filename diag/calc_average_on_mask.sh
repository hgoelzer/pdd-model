# Calculate average on ice mask
#scp nird:/nird/datapeak/NS5011K/ISMIP/ISMIP7/GrIS/Data/FGobs/process/icemask_promice_08000m.nc ./

maskfile=icemask_promice_08000m.nc

infile=pr.nc
tmp_file=pr_tmp.nc
avg_file=pr_avg.nc

cdo yearmean $infile $tmp_file
ncap2 -O -s "pr=pr*31556926" $tmp_file $tmp_file
ncatted -h -a units,pr,o,c,"mm w.e. yr-1" $tmp_file

ncks -A $maskfile $tmp_file

# add missing value
ncatted -a _FillValue,pr,o,d,-99999.0 $tmp_file

ncap2 -O -s 'pr_masked=pr; where(icemask_promice == 0) pr_masked=pr.get_miss()' -v $tmp_file $tmp_file
ncatted -a _FillValue,pr_masked,o,d,-99999.0 $tmp_file 
ncap2 -O -s 'pr_avg=pr_masked.avg($x,$y)' -v $tmp_file $avg_file

# Clean up
/bin/bash pr_tmp.nc
