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
# Olivia requires NRIS/CPU to expose the module tree before netCDF-Fortran is visible;
# on Betzy this module doesn't exist and the load is a no-op.
module load NRIS/CPU 2>/dev/null || true
module load netCDF-Fortran/4.6.1-gompi-2024a
# Explicitly export library path so srun tasks inherit it even when SLURM
# resets the task environment (some cluster configurations strip LD_LIBRARY_PATH).
export LD_LIBRARY_PATH="${EBROOTNETCDFMINFORTRAN}/lib:${LD_LIBRARY_PATH:-}"
