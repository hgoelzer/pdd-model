#!/bin/bash
#SBATCH --job-name=gris_sc08
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16G
#
# GrIS ISMIP7 seasonalcycle=0 sensitivity at 8 km: CESM2 hist 1950-2014 +
# ssp370 2015-2100, production namelists except seasonalcycle=0 (actual
# monthly temperatures instead of the prescribed taj cycle). The sc1 twin
# is the verified std rerun in ./output_ismip7 - no rerun needed.
# Usage: sbatch gris_sc08_job.sh

set -o errexit
set -o nounset

cd /cluster/projects/nn5011k/heig/PDD/gris_08/pdd-model
source setenv.sh
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python

NML_HIST=params_CESM2_historical_anom_08km_ismip7_sc0.nml
NML_SSP=params_CESM2_ssp370_anom_08km_ismip7_sc0.nml
STD_H=params_CESM2_historical_anom_08km_ismip7.nml
STD_S=params_CESM2_ssp370_anom_08km_ismip7.nml

echo "== 1. Namelist pairing: sc0 differs from std only in seasonalcycle/outpathname =="
for pair in "$NML_HIST $STD_H" "$NML_SSP $STD_S"; do
    set -- $pair
    diff <(grep -vE 'seasonalcycle|outpathname|^ *!' "$1") \
         <(grep -vE 'seasonalcycle|outpathname|^ *!' "$2") >/dev/null \
        || { echo "ERROR: $1 vs $2 differ beyond seasonalcycle/outpathname"; exit 1; }
    grep -q 'seasonalcycle = 0' "$1" || { echo "ERROR: $1 not seasonalcycle=0"; exit 1; }
done
echo "pairing OK"

echo "== 2. Build serial driver =="
make clean >/dev/null
make gpdd_monthly_inout >/dev/null
echo "build done"

echo "== 3. sc0 runs =="
rm -rf output_sc08
mkdir -p output_sc08/Historical output_sc08/SSP370
./gpdd_monthly_inout.x "$NML_HIST" > log_sc08_hist.txt 2>&1
echo "historical 1950-2014 (sc0) done"
./gpdd_monthly_inout.x "$NML_SSP" > log_sc08_ssp.txt 2>&1
echo "ssp370 2015-2100 (sc0) done"

echo "== 4. Verification + sc0 vs std (taj) comparison =="
cat > check_gris_sc08.py <<'PYEOF'
import glob
import sys
import netCDF4 as nc

AREA = 8000.0 ** 2
SPM = 30 * 86400.0  # SMBMIP 30-day months
PHASES = [
    ("Historical", "GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1",
     range(1950, 2015), 585),
    ("SSP370", "GIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1",
     range(2015, 2101), 774),
]

def gtyr(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:].sum().item() * AREA * SPM / 1e12

fails = []
res = {}
for tag, root, years, nexp in PHASES:
    for base in ("output_sc08", "output_ismip7"):
        n = len(glob.glob(f"{base}/{tag}/{root}_*.nc"))
        if n != nexp:
            fails.append(f"{base}/{tag}: {n} files != {nexp}")
            print(f"FAIL {base}/{tag}: {n} files, expected {nexp}")
            continue
        smb = {}
        for y in years:
            acabf = gtyr(f"{base}/{tag}/acabf_{root}_{y}.nc", "acabf")
            mrro = gtyr(f"{base}/{tag}/mrro_{root}_{y}.nc", "mrro")
            pr = gtyr(f"{base}/{tag}/pr_{root}_{y}.nc", "pr")
            if abs(acabf + mrro - pr) > 1e-3:
                fails.append(f"{base}/{tag} {y} closure "
                             f"{abs(acabf+mrro-pr):.2e}")
            smb[y] = acabf
        res[base] = smb
    # grid metadata on the new output
    with nc.Dataset(f"output_sc08/{tag}/acabf_{root}_{years[0]}.nc") as d:
        x = d.variables["x"][:]
        yy = d.variables["y"][:]
        gi = getattr(d, "grid_information", "MISSING")
        ok = (abs(x[0] + 720000) < 1 and abs(yy[0] + 3450000) < 1
              and "3413" in gi)
        print(f"{tag} grid: x {x[0]:.0f} y {yy[0]:.0f} {gi} "
              f"{'OK' if ok else 'FAIL'}")
        if not ok:
            fails.append(f"{tag} grid metadata")

for tag, root, years, _ in PHASES:
    if "output_sc08" not in res or "output_ismip7" not in res:
        continue
    print(f"year  sc0_SMB  taj_SMB   delta ({tag})")
    for y in years:
        d = res["output_sc08"][y] - res["output_ismip7"][y]
        print(f"{y}  {res['output_sc08'][y]:8.2f} "
              f"{res['output_ismip7'][y]:8.2f} {d:+8.2f}")
    yrs = list(years)
    m0 = sum(res["output_sc08"][y] for y in yrs) / len(yrs)
    m1 = sum(res["output_ismip7"][y] for y in yrs) / len(yrs)
    dmax = max(abs(res["output_sc08"][y] - res["output_ismip7"][y])
               for y in yrs)
    print(f"{tag} mean: sc0 {m0:.2f}  taj {m1:.2f}  "
          f"delta {m0-m1:+.2f} Gt/yr  max|yr delta| {dmax:.2f}")

if fails:
    print("RESULT: FAILURES:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
PYEOF
$PY check_gris_sc08.py

echo "== DONE =="
