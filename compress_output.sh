#!/bin/bash
# Compress nc output

set -e
set -x

files=`ls output/ | grep ".nc"`

for file in $files; do
    echo $file
    nccopy -d9 output/$file tmp.nc
    mv tmp.nc output/$file
done

