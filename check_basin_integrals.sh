#!/bin/bash
# integrate all basins, to be compared to scalars_mm

set -x
set -e

ncmodel="scalars_mm_01.nc"
ncrignot="scalars_rm_01.nc"
nczwally="scalars_zm_01.nc"

ncap2 -O -s "smb = smb_no+smb_ne+smb_se+smb_sw+smb_cw+smb_nw" -v ${ncrignot} integrals_rignot.nc

ncap2 -O -s "smb = smb_z11+ smb_z12+ smb_z13+ smb_z14 + smb_z21+ smb_z22+ smb_z31+ smb_z32+ smb_z33+ smb_z41+ smb_z42+ smb_z43+ smb_z50+ smb_z61+ smb_z62+ smb_z71+ smb_z72+ smb_z81+ smb_z82" -v ${nczwally} integrals_zwally.nc

ncdump -v smb ${ncmodel}
ncdump -v smb integrals_rignot.nc 
ncdump -v smb integrals_zwally.nc 


#ncdiff -v smb integrals_rignot.nc ${ncmodel} df_rignot.nc
#ncdiff -v smb integrals_zwally.nc ${ncmodel} df_zwally.nc

