#!/bin/bash
#SBATCH --job-name=gris_sc0
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#
# GrIS ISMIP7 sensitivity run: CESM2 historical 1950-2014, 16 km, with
# seasonalcycle=0 (actual monthly temperatures) vs the taj twin
# (seasonalcycle=1 = production configuration). Paired runs quantify the
# effect of the seasonal-cycle treatment on Greenland. Default physics
# (0.00297/0.00791/4.5/1.0, pmax=0.3). Serial driver; 16 km runs are fast.
# Usage: sbatch gris_sc0_job.sh

set -o errexit
set -o nounset

RUNDIR=/cluster/projects/nn5011k/heig/PDD/gris_16/pdd-model
cd "$RUNDIR"
source setenv.sh
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python

NML_SC0=params_CESM2_historical_anom_16km_sc0.nml
NML_SC1=params_CESM2_historical_anom_16km_sc1.nml

echo "== 1. Verify code version and namelist pairing =="
git merge-base --is-ancestor dacfb42 HEAD || { echo "ERROR: run dir missing commit dacfb42 - update the clone first"; exit 1; }
git log --oneline -1
# sc0 and sc1 must differ only in seasonalcycle and outpathname
diff <(grep -v 'seasonalcycle\|outpathname' "$NML_SC0") \
     <(grep -v 'seasonalcycle\|outpathname' "$NML_SC1") >/dev/null \
     || { echo "ERROR: sc0/sc1 namelists differ beyond seasonalcycle/outpathname"; exit 1; }

echo "== 2. Build serial driver =="
make clean
make gpdd_monthly_inout

echo "== 3. Paired runs (sc1 taj twin = production reference, then sc0) =="
rm -rf output_sc1 output_sc0
mkdir -p output_sc1 output_sc0
./gpdd_monthly_inout.x "$NML_SC1" > log_sc1.txt 2>&1
echo "sc1 (taj) done"
./gpdd_monthly_inout.x "$NML_SC0" > log_sc0.txt 2>&1
echo "sc0 (actual temps) done"

echo "== 4. Verification =="
cat > check_gris_sc.py <<'PYEOF'
import glob
import sys
import netCDF4 as nc

CELL_AREA = 16000.0 ** 2
SPM = 30 * 86400.0  # SMBMIP 30-day months, as in test_regression.sh
ROOT = "GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1"
YEARS = range(1950, 2015)

def gtyr(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:].sum().item() * CELL_AREA * SPM / 1e12

fails = []
res = {}
for tag, outdir in [("sc1", "output_sc1"), ("sc0", "output_sc0")]:
    fs = glob.glob(f"{outdir}/acabf_{ROOT}_*.nc")
    if len(fs) != 65:
        fails.append(f"{tag}: {len(fs)} acabf files != 65")
        print(f"FAIL {tag}: {len(fs)} acabf files, expected 65")
        continue
    nv = len(glob.glob(f"{outdir}/*_1950.nc"))
    print(f"{tag}: 65 yr, {nv} vars/yr")
    smb = {}
    for y in YEARS:
        acabf = gtyr(f"{outdir}/acabf_{ROOT}_{y}.nc", "acabf")
        mrro = gtyr(f"{outdir}/mrro_{ROOT}_{y}.nc", "mrro")
        pr = gtyr(f"{outdir}/pr_{ROOT}_{y}.nc", "pr")
        err = abs(acabf + mrro - pr)
        if err > 1e-3:
            fails.append(f"{tag} {y} closure {err:.2e}")
        smb[y] = acabf
    res[tag] = smb
    # grid metadata (GrIS origins, epsg:3413)
    with nc.Dataset(f"{outdir}/acabf_{ROOT}_1950.nc") as d:
        x = d.variables["x"][:]
        yy = d.variables["y"][:]
        gi = getattr(d, "grid_information", "MISSING")
        ok = (abs(x[0] + 720000) < 1 and abs(yy[0] + 3450000) < 1
              and "3413" in gi)
        print(f"{tag} grid: x {x[0]:.0f}..{x[-1]:.0f} y {yy[0]:.0f}..{yy[-1]:.0f} "
              f"{gi} {'OK' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"{tag} grid metadata")

if "sc1" in res:
    b = res["sc1"][1950]
    print(f"sc1 1950 SMB = {b:.2f} Gt/yr (regression baseline 198.06)")
    if abs(b - 198.06) > 1.0:
        fails.append(f"sc1 1950 baseline {b:.2f} != 198.06")

if len(res) == 2:
    print("year   sc1_SMB  sc0_SMB   delta")
    for y in YEARS:
        d = res["sc0"][y] - res["sc1"][y]
        print(f"{y}  {res['sc1'][y]:8.2f} {res['sc0'][y]:8.2f} {d:+8.2f}")
    m1 = sum(res["sc1"].values()) / 65
    m0 = sum(res["sc0"].values()) / 65
    l1 = sum(res["sc1"][y] for y in range(2005, 2015)) / 10
    l0 = sum(res["sc0"][y] for y in range(2005, 2015)) / 10
    print(f"mean 1950-2014: sc1 {m1:.2f}  sc0 {m0:.2f}  delta {m0-m1:+.2f} Gt/yr")
    print(f"mean 2005-2014: sc1 {l1:.2f}  sc0 {l0:.2f}  delta {l0-l1:+.2f} Gt/yr")

if fails:
    print("RESULT: FAILURES:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
PYEOF
$PY check_gris_sc.py

echo "== DONE =="