# Get forcing from nird

# Typically create Forcing directory at some level below and link in
mkdir /cluster/work/users/heig/pdd/Forcing
cd /cluster/work/users/heig/pdd/Forcing

mkdir MARv3.14
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/MAR/MARv3.14/ismip6/*_ref_MARv3.14-monthly-ERA5-1960-1989_01000m.nc ./MARv3.14/

mkdir MARv3.14LC26
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/MAR/Lindsey-Clark2026/ismip6/pr_ref_MARv3.14LC26-monthly-ERA5-1960-1989_01000m.nc ./MARv3.14LC26/
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/MAR/MARv3.14/ismip6/tas_ref_MARv3.14-monthly-ERA5-1960-1989_01000m.nc ./MARv3.14LC26/

mkdir -p CESM2/ssp370/tas_anom
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/CMIP/ismip6/CESM2/ssp370/tas_anom/tas_Amon_CESM2_ssp370_r11i1p1f1_*_i01000m.nc CESM2/ssp370/tas_anom/
mkdir -p CESM2/ssp370/pr_ratio
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/CMIP/ismip6/CESM2/ssp370/pr_ratio/pr_ratio_Amon_CESM2_ssp370_r11i1p1f1_*_i01000m.nc CESM2/ssp370/pr_ratio/

mkdir -p CESM2/historical/tas_anom
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/CMIP/ismip6/CESM2/historical/tas_anom/tas_Amon_CESM2_historical_r11i1p1f1_*_i01000m.nc CESM2/historical/tas_anom/
mkdir -p CESM2/historical/pr_ratio
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/CMIP/ismip6/CESM2/historical/pr_ratio/pr_ratio_Amon_CESM2_historical_r11i1p1f1_*_i01000m.nc CESM2/historical/pr_ratio/

mkdir -p constant
rsync -av nird:/nird/datalake/NS8002K/heig/PDD/constant/sftgif_BM5_01000m.nc ./constant/
cp ../data/icemask_promice_01000m.nc ./constant/sftgif_promice_01000m.nc
ncrename -v icemask_promice,sftgif ./constant/sftgif_promice_01000m.nc

# link back in
#cd ... pdd-model/Forcing
ln -s /cluster/work/users/heig/pdd/Forcing/*/ ./
