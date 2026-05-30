#!/bin/bash
# Source this file to load the NetCDF-Fortran module environment.
# Keeps compile and runtime environments identical.
#
# Interactive compile:  source setenv.sh && make gpdd_monthly_inout
# In SLURM scripts:     source "$(dirname "$0")/setenv.sh"
#
# Same module on bet (Betzy) and oli (Olivia) as of 2026-05.
# Update here when upgrading toolchain — compile and run pick it up automatically.

module purge
module load netCDF-Fortran/4.6.1-gompi-2024a
