#!/bin/bash

##############################################################################
# Script: check_netcdf_compression.sh
# Purpose: Check if NetCDF file is compressed using ncdump
# Usage: ./check_netcdf_compression.sh <file.nc> [file2.nc ...]
##############################################################################

if [ $# -eq 0 ]; then
    echo "Usage: $0 <file.nc> [file2.nc ...]"
    exit 1
fi

for filepath in "$@"; do
    if [ ! -f "$filepath" ]; then
        echo "Error: File not found - $filepath"
        continue
    fi
    
    echo "Checking: $filepath"
    
    # Use ncdump -hs to show storage format (includes compression info)
    ncdump_output=$(ncdump -hs "$filepath" 2>/dev/null)
    
    # Check for compression indicators in storage info
    if echo "$ncdump_output" | grep -qi "deflate\|szip\|chunked" >/dev/null 2>&1; then
        echo "  ✓ COMPRESSED"
    else
        echo "  ✗ NOT COMPRESSED"
    fi
    echo ""
done
