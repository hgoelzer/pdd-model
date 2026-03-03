# Runs on ISMIP6 grid 

# Data from nird
cp /nird/projects/nird/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*yearly*_e01000.nc ./data
cp /nird/projects/nird/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*monthly-07-*_e01000.nc ./data

# parameters +10%
ddfactorsnow = 0.0033
ddfactorice = 0.0088
sigma = 4.95 


# Setup
make clean
make gpdd

# Test run
./gpdd.x 

# Run
sbatch runPDD

# cut out 1991-2018 for MAR comparison
ncks -F -v smb -d time,1,28 smb_gpdd.nc smb_CARRAPDD_1991-2018.nc


