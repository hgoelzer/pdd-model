# Runs on ISMIP6 grid 

# Data from nird annual and july
scp nird:/nird/datalake/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*yearly*_e01000.nc ./data
scp nird:/nird/datalake/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*monthly-07*_e01000.nc ./data

# Setup
make clean
make gpdd

# Test run
./gpdd.x 

# Run
sbatch runPDD

# cut out 1991-2018 for MAR comparison
ncks -F -v smb -d time,1,28 smb_gpdd.nc smb_CARRAPDD_1991-2018.nc




# Full monthly forcing
# Data from nird 
scp nird:/nird/datalake/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*monthly*_e01000.nc ./data

# Setup
make clean
make gpdd_monthly

acabf_GIS_MARv3.14_CESM2_SSP370_r11i1p1f1_2015.nc
