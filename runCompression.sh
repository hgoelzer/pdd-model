#!/bin/bash
#SBATCH --job-name=gpdd
# CLIM2Ant: nn11016k, ISMIP7: nn5011k, 
#SBATCH --account=nn11016k
#SBATCH --time=02:40:00
#SBATCH --mem-per-cpu=16G
#SBATCH --qos=preproc
#SBATCH --ntasks=1

## Setting variables and prepare runtime environment:
##----------------------------------------------------
## Recommended safety settings:
set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

time srun ./compress_output_with_checks.sh

# Finish the script
exit 0
