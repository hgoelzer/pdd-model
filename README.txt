# Runs on ISMIP6 grid 

# Data from nird

cp /nird/projects/nird/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*yearly*_e01000.nc ./data
cp /nird/projects/nird/NS5011K/MIPs/IMBIE3/CARRA/ISMIP6/*monthly-07-*_e01000.nc ./data


# Setup
make clean
make gpdd

# Test run
./gpdd.x 

# Run
sbatch runPDD
