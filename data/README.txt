# input data for diag processing

# Typically create data directory at some level below and link in
mkdir /cluster/work/users/heig/pdd/data
cd /cluster/work/users/heig/pdd/data

rsync -av nird:/nird/datalake/NS8002K/heig/ISMIP6/GrIS/Data/*01000* ./

rsync -av nird:/nird/datapeak/NS5011K/ISMIP/ISMIP7/GrIS/Data/FGobs/process/icemask_promice_01000m.nc ./

#cp ../Forcing/constant/sftgif_BM5_01000m.nc ./

# link back in
#cd ... pdd-model/data
ln -s /cluster/work/users/heig/pdd/data/*.nc ./
