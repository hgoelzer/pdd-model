# PDD model 
Implements PDD model of Huybrechts and De Wolde (1999)

Using ncio by Alexander Robinson (git@github.com:alex-robinson/ncio.git)

# data and Forcing, tyipcally linked in from below
data/

Forcing/

# Setup
edit gpdd_monthly_inout.f90 to set forcing, timing, exp

make clean

make gpdd_monthly_inout


# Run
./gpdd_monthly_inout.x

# Diagostics, scalar output, mass changes
diag/
