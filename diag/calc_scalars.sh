#!/bin/bash
# Calculate scalar variables from 3D model output
# Heiko Goelzer 2019 (heig@norceresearch.no)

# Scalars
# multiply with masks, multiply with area factors integrate spatially

#if [[ $# -ne 3 ]]; then
#    echo "Illegal number of parameters. Need 3 to specify masking for GIC and OBS and resolution"
#    exit 3
#fi

set -x
set -e

# Resolution
#ares=01
ares=08
#ares=${3}

#run=output_MAR
run=output

# Path to mask data
datapath=../data

## Prepare new master file
flg_master=true  # Prepare master
#flg_master=false  # Reuse master

## What output to process
flg_mm=false # Integrals on model mask
flg_rm=false # IMBIE2-Rignot basins

flg_mm=true  # Integrals on model mask
#flg_rm=true  # IMBIE2-Rignot basins


## What masking to apply If true, applied to all output
# Remove GIC contribution? 
flg_GICmask=false # [Default true!]
#flg_GICmask=${1} # [Default true!]
# Remove ice outside observed ice mask (can be combined with GIC masking) 
flg_OBSmask=true # [Default false!]
#flg_OBSmask=${2} # [Default false!]

# Model data
infile=${run}/smb_gpdd_${ares}000m.nc
maskfile=icemasks_${ares}000m.nc

# area factors
af2input=$datapath/af2_ISMIP6_GrIS_${ares}000m.nc
af2file=af2.nc
# af2

# Ellesmere mask
emfile=$datapath/maxmask1_${ares}000m.nc
# maxmask1

# BM3 for observed masks
obsinput=$datapath/BM3_GrIS_nn_e${ares}000m.nc
obsfile=masksOBS.nc

# PROMICE mask
#prominput=$datapath/sftgif_BM5_${ares}000m.nc
prominput=$datapath/icemask_promice_${ares}000m.nc
promfile=maskPROM.nc

# IMBIE2 Rignot extended masks
riginput=$datapath/GrIS_Basins_Rignot_extended_e${ares}000m_v1.nc
rigfile=masksRIG.nc

# GIC area factors
gicinput=$datapath/rgi60_connect01_iaf2_${ares}000m_v1.nc
gicfile=masksGIC.nc
# iaf2

# smb 

# Possible output files
scfile_mm=${run}/scalars_mm_${ares}.nc
scfile_rm=${run}/scalars_rm_${ares}.nc


if $flg_master; then

    # Make a netcdf file with parameters
    ncks -3 -O -v x,y ${af2input} tmp.nc
    # BM3 numbers for density, use Cogley (2012) for ocean area
    #ncap2 -A -s 'rhoi=916.7; rhow=1027.0; rhof=1000.0; oarea=3.618e14' -v tmp.nc params_${ares}.nc
    ncap2 -3 -A -s 'rhoi=916.7; rhow=1027.0; rhof=1000.0; oarea=3.625e14' -v tmp.nc params_${ares}.nc
    # Model specific densities 
    #ncks -3 -A -v rhoi,rhow,rhof,oarea ${infile} params_${ares}.nc

    # Resolution
    #ncap2 -3 -A -s 'dx=x(2)-x(1); dy=y(2)-y(1)' -v tmp.nc params_${ares}.nc
    #ncap2 -3 -A -s 'dx=8000; dy=8000' -v tmp.nc params_${ares}.nc
    #ncap2 -3 -A -s 'dx=1000; dy=1000' -v tmp.nc params_${ares}.nc
    ncap2 -3 -A -s "dx=${ares}000; dy=${ares}000" -v tmp.nc params_${ares}.nc
    ncks -3 -C -O -x -v x,y params_${ares}.nc params_${ares}.nc
    # clean up
    ncatted -h -a history,global,d,, params_${ares}.nc
    ncatted -h -a history_of_appended_files,global,d,, params_${ares}.nc
    ncatted -h -a NCO,global,d,, params_${ares}.nc
    ncatted -h -a CDO,global,d,, params_${ares}.nc
    ncatted -h -a CDI,global,d,, params_${ares}.nc

    /bin/rm tmp.nc

    # prepare BM masks
    ncks -3 -O -v sftgif ${obsinput} ${obsfile} 
    ncks -3 -A -v sftgrf ${obsinput} ${obsfile}
    ncks -3 -A -v sftflf ${obsinput} ${obsfile}
    ncrename -v sftgif,sftgif_BM ${obsfile}
    ncrename -v sftgrf,sftgrf_BM ${obsfile}
    ncrename -v sftflf,sftflf_BM ${obsfile}

    # prepare PROMICE mask
    #ncks -3 -O -v sftgif ${prominput} ${promfile}
    #ncrename -v sftgif,sftgif_PROM ${promfile}
    ncks -3 -O -v icemask_promice ${prominput} ${promfile}
    ncrename -v icemask_promice,sftgif_PROM ${promfile}
    
    # Prepare IMBIE2 Rignot masks, ID: From NO clockwise
    if $flg_rm; then
	ncap2 -3 -O -s 'no=ID==1' -v ${riginput} ${rigfile} 
	ncap2 -3 -A -s 'ne=ID==2' -v ${riginput} ${rigfile} 
	ncap2 -3 -A -s 'se=ID==3' -v ${riginput} ${rigfile} 
	ncap2 -3 -A -s 'sw=ID==4' -v ${riginput} ${rigfile} 
	ncap2 -3 -A -s 'cw=ID==5' -v ${riginput} ${rigfile} 
	ncap2 -3 -A -s 'nw=ID==6' -v ${riginput} ${rigfile} 
    fi


    # Prepare area factors
    ncks -3 -O -v af2 ${af2input} ${af2file}
    # inlcude GIC masking if requested
    if $flg_GICmask; then
	# Prepare GIC masking file
	ncks -3 -O -v iaf2 ${gicinput} ${gicfile}
	# Combine with area factors
	ncks -A -v iaf2 ${gicfile} ${af2file}
	ncap2 -O -s "af2 = af2*iaf2" ${af2file} ${af2file}
    fi

    # Prepare ice masks
    #ncks -3 -O -v sftgif ${infile} ${maskfile}
    #ncks -3 -A -v sftgrf ${infile} ${maskfile}
    #ncks -3 -A -v sftflf ${infile} ${maskfile}
    # inlcude OBS masking if requested
    if $flg_OBSmask; then
	ncks -3 -O -v sftgif_BM ${obsfile} ${maskfile}
	ncks -3 -A -v sftgrf_BM ${obsfile} ${maskfile}
	ncks -3 -A -v sftflf_BM ${obsfile} ${maskfile}
	ncks -3 -A -v sftgif_PROM ${promfile} ${maskfile}
	ncap2 -3 -A -s "sftgif = sftgif_BM" ${maskfile} ${maskfile}
	ncap2 -3 -A -s "sftgrf = sftgrf_BM" ${maskfile} ${maskfile}
	ncap2 -3 -A -s "sftflf = sftflf_BM" ${maskfile} ${maskfile}
    fi

    # Make master netcdf file 
    ncks -3 -O -v sftgif ${maskfile} tmp_mod.nc
    ncks -3 -A -v sftgrf ${maskfile} tmp_mod.nc
    ncks -3 -A -v sftflf ${maskfile} tmp_mod.nc
    ncks -3 -A -v sftgif_PROM ${maskfile} tmp_mod.nc
    ncks -3 -A -v maxmask1 ${emfile} tmp_mod.nc
    ncks -3 -A -v af2 ${af2file} tmp_mod.nc

    if $flg_rm; then
	ncks -3 -A  ${rigfile} tmp_mod.nc
    fi
    # Add porameters
    ncks -3 -A params_${ares}.nc tmp_mod.nc

    # Remember master file
    /bin/cp tmp_mod.nc master_${ares}000m.nc 
else
    /bin/cp master_${ares}000m.nc tmp_mod.nc 
fi

# Add model run specific data
ncks -3 -A -v acabf ${infile} tmp_mod.nc
ncks -3 -A -v pr ${infile} tmp_mod.nc

##################################################################################
# Greenland wide integrals
##################################################################################

if $flg_mm; then

# Make dummy containers for scalar output
ncks -3 -O params_${ares}.nc ${scfile_mm}

# Integrated acabf over grounded ice
/bin/cp params_${ares}.nc tmpaf.nc
ncap2 -A -s 'af=acabf*sftgif_PROM*maxmask1*af2' -v tmp_mod.nc tmpaf.nc
#/bin/cp tmpaf.nc test.nc
ncap2 -O -s 'acabf=af.total($x,$y)*dx^2' -v tmpaf.nc tmpsc.nc
ncks -A -v acabf tmpsc.nc ${scfile_mm}
ncatted -a standard_name,acabf,d,, ${scfile_mm}
ncatted -a long_name,acabf,o,c,"surface mass balance" ${scfile_mm}
ncatted -a units,acabf,o,c,"kg yr-1" ${scfile_mm}
/bin/rm tmpaf.nc tmpsc.nc 

# Integrated pr over grounded ice
/bin/cp params_${ares}.nc tmpaf.nc
ncap2 -A -s 'af=pr*sftgif_PROM*maxmask1*af2' -v tmp_mod.nc tmpaf.nc
ncap2 -O -s 'tp=af.total($x,$y)*dx^2' -v tmpaf.nc tmpsc.nc
ncks -A -v tp tmpsc.nc ${scfile_mm}
ncatted -a standard_name,tp,d,, ${scfile_mm}
ncatted -a long_name,tp,o,c,"total precipitation" ${scfile_mm}
ncatted -a units,tp,o,c,"kg yr-1" ${scfile_mm}
/bin/rm tmpaf.nc tmpsc.nc 

# clean up
ncatted -h -a history,global,d,, ${scfile_mm}
ncatted -h -a history_of_appended_files,global,d,, ${scfile_mm}
ncatted -h -a NCO,global,d,, ${scfile_mm}
ncatted -h -a CDO,global,d,, ${scfile_mm}
ncatted -h -a CDI,global,d,, ${scfile_mm}
ncatted -h -a Description,global,d,, ${scfile_mm}
ncatted -h -a proj4,global,d,, ${scfile_mm}
# Add info
ncatted -h -a Description,global,o,c,"ISMIP6-Greenland recalculated scalar output. Heiko Goelzer 2026, heig@norceresearch.no" ${scfile_mm}

fi

##################################################################################
# IMBIE2 Rignot basins 
##################################################################################
if $flg_rm; then

# Make dummy containers for scalar output
ncks -O params_${ares}.nc ${scfile_rm}

ncks -A params_${ares}.nc tmp.nc

for basin in no ne se sw cw nw; do

# Integrated acabf
ncap2 -O -s "af=acabf*sftgif_PROM*maxmask1*af2*${basin}" tmp_mod.nc tmpaf.nc
ncap2 -O -s 'acabf=af.total($x,$y)*dx^2' -v tmpaf.nc tmpsc.nc
ncrename -v acabf,acabf_${basin} tmpsc.nc
ncks -A -v acabf_${basin} tmpsc.nc ${scfile_rm}
ncatted -a standard_name,acabf_${basin},d,, ${scfile_rm}
ncatted -a long_name,acabf_${basin},o,c,"surface mass balance" ${scfile_rm}
ncatted -a units,acabf_${basin},o,c,"kg yr-1" ${scfile_rm}
/bin/rm tmpaf.nc tmpsc.nc 

# Integrated tp
ncap2 -O -s "af=pr*sftgif_PROM*maxmask1*af2*${basin}" tmp_mod.nc tmpaf.nc
ncap2 -O -s 'tp=af.total($x,$y)*dx^2' -v tmpaf.nc tmpsc.nc
ncrename -v tp,tp_${basin} tmpsc.nc
ncks -A -v tp_${basin} tmpsc.nc ${scfile_rm}
ncatted -a standard_name,tp_${basin},d,, ${scfile_rm}
ncatted -a long_name,tp_${basin},o,c,"total precipitation" ${scfile_rm}
ncatted -a units,tp_${basin},o,c,"kg yr-1" ${scfile_rm}
/bin/rm tmpaf.nc tmpsc.nc 

done

# clean up
ncatted -h -a history,global,d,, ${scfile_rm}
ncatted -h -a history_of_appended_files,global,d,, ${scfile_rm}
ncatted -h -a NCO,global,d,, ${scfile_rm}
ncatted -h -a CDO,global,d,, ${scfile_rm}
ncatted -h -a CDI,global,d,, ${scfile_rm}
ncatted -h -a Description,global,d,, ${scfile_rm}
ncatted -h -a proj4,global,d,, ${scfile_rm}
# Add info
ncatted -h -a Description,global,o,c,"ISMIP6-Greenland recalculated scalar output. Heiko Goelzer 2026, heig@norceresearch.no" ${scfile_rm}

fi


# Clean up
/bin/rm masksOBS.nc maskPROM.nc af2.nc icemasks_08000m.nc 
/bin/rm tmp_mod.nc

