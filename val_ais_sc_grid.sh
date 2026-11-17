#!/bin/bash
#SBATCH --job-name=ais_val
# Betzy: ISMIP7 nn5011k (quota exhausted) -> nn8002k per commit 178e415
#SBATCH --account=nn8002k
#SBATCH --time=01:00:00
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8G
#SBATCH --qos=preproc

set -o errexit
set -o nounset

RUNDIR=/cluster/projects/nn5011k/heig/PDD/ais_08/pdd-model
cd "$RUNDIR"
source setenv.sh
NCDUMP=$(dirname "${EBROOTNETCDFMINFORTRAN}")/../bin/ncdump || NCDUMP=ncdump
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python

echo "== 1. Verify code version (repo updated outside the job; compute nodes cannot reach GitHub) =="
git merge-base --is-ancestor 90c77b8 HEAD || { echo "ERROR: run dir missing commit 90c77b8 — update the clone first"; exit 1; }
git log --oneline -1

echo "== 2. Compile serial + MPI drivers =="
make clean
make gpdd_monthly_inout gpdd_monthly_inout_mpi

echo "== 3. Write validation namelists (1 yr, historical 1950) =="
COMMON_HEADER='&config
  outputmode = 0
  outputvars = 0
  fmode = 1
  inpathname_pr  = "./Forcing/CESM2/historical/pr_ratio"
  inpathname_tas = "./Forcing/CESM2/historical/tas_anom"
  fileroot_pr    = "pr_ratio_Amon_CESM2_historical_r11i1p1f1"
  fileroot_tas   = "tas_Amon_CESM2_historical_r11i1p1f1"
  year0        = 1950
  nt           = 1
  fileroot_out = "AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1"
  nx         = 761
  ny         = 761
  res        = 8000.0
  res_suffix = "i08000m"
  filename_prref   = "Forcing/MARv3.14/pr_ref_MARv3.14-monthly-ERA5-1980-2009_08000m.nc"
  filename_tasref  = "Forcing/MARv3.14/tas_ref_MARv3.14-monthly-ERA5-1980-2009_08000m.nc"
  filename_defmask = "Forcing/constant/sftgif_BM3_08000m.nc"
  institution   = "NORCE Research"
  contact_name  = "Heiko Goelzer"
  contact_email = "heig@norceresearch.no"
  ddfactorsnow = 0.00297
  ddfactorice  = 0.00791
  sigma        = 4.5
  rainlimit    = 1.0
  pmax = 0.3'

printf '%s\n  seasonalcycle = 1\n  x0 = -3040000.0\n  y0 = -3040000.0\n  grid_information = "ISMIP6 grid, epsg:3031"\n/\n' "$COMMON_HEADER" > params_val_A.nml
printf '%s\n  seasonalcycle = 0\n  x0 = -3040000.0\n  y0 = -3040000.0\n  grid_information = "ISMIP6 grid, epsg:3031"\n/\n' "$COMMON_HEADER" > params_val_B.nml

echo "== 4. Run A (seasonalcycle=1, archive-repro) and B (seasonalcycle=0, SH phase) =="
mkdir -p output_valA output_valB
sed 's|outpathname.*|outpathname = "./output_valA"|' params_val_A.nml > /tmp/na.nml && mv /tmp/na.nml params_val_A.nml
sed 's|outpathname.*|outpathname = "./output_valB"|' params_val_B.nml > /tmp/nb.nml && mv /tmp/nb.nml params_val_B.nml
srun ./gpdd_monthly_inout_mpi.x params_val_A.nml > log_val_A.txt 2>&1
srun ./gpdd_monthly_inout_mpi.x params_val_B.nml > log_val_B.txt 2>&1

echo "== 5. Grid metadata check (run A acabf header) =="
$NCDUMP -h output_valA/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc | grep -E 'x0|x =|actual_range|grid_information' || true

echo "== 6. Annual scalar sums (acabf, mrro, prsn, prra) A vs B =="
$PY - <<'EOF'
import netCDF4 as nc
A = "output_valA/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc"
B = "output_valB/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc"
for v in ["acabf", "mrro", "prsn", "prra"]:
    for tag, f in [("A(sc=1)", A), ("B(sc=0)", B)]:
        d = nc.Dataset(f)
        a = d.variables[v][:]
        print(f"{v:6s} {tag}  annual_mean_Gt/yr = {a.mean().item()*761*761*64e6*31556926/1e12:.3f}")
        d.close()
# month-by-month difference A vs B
dA = nc.Dataset(A); dB = nc.Dataset(B)
diff = (dA.variables["acabf"][:] - dB.variables["acabf"][:]).mean(axis=(0,1))
print("acabf monthly mean diff A-B (kg m-2 s-1):")
print(["%.2e" % x for x in diff.tolist()])
EOF

echo "== DONE =="
exit 0
