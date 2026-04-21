# Runs on ISMIP6 grid 

# data and Forcing, tyipcally linked in from below
data/
Forcing/

# Setup
edit gpdd_monthly_inout.f90 to set forcing, timing, exp
make clean
make gpdd_monthly_inout

# Test run
./gpdd_monthly_inout.x

# Diagostics, scalar output, mass changes
diag/
