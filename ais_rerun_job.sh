#!/bin/bash
#SBATCH --job-name=ais_rerun
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#
# AIS SMBMIP rerun 2026-11, corrected configuration:
#   - CESM2 historical 1950-2014 + ssp370 2015-2100, 8 km ISMIP6 AIS grid
#   - original NORCEPDD1-AIS physics (uncalibrated GrIS defaults
#     0.00297/0.00791/4.5/1.0, pmax=0.3)
#   - fixes: seasonalcycle=0 (SH phase), epsg:3031 grid origin,
#     mrro=total runoff, sir=retained meltwater (8 output vars, rfr dropped)
# Serial driver (1-yr 8 km run ~ 15 s, measured job 1807140).
# Usage: sbatch ais_rerun_job.sh

set -o errexit
set -o nounset

RUNDIR=/cluster/projects/nn5011k/heig/PDD/ais_08/pdd-model
cd "$RUNDIR"
source setenv.sh
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python

NML_HIST=params_CESM2_historical_anom_ais_08km.nml
NML_SSP=params_CESM2_ssp370_anom_ais_08km.nml
NML_HIST_VAL=params_CESM2_historical_anom_ais_08km_val1yr.nml
NML_SSP_VAL=params_CESM2_ssp370_anom_ais_08km_val1yr.nml

echo "== 1. Verify code version and namelists =="
git merge-base --is-ancestor 90c77b8 HEAD || { echo "ERROR: run dir missing commit 90c77b8 - update the clone first"; exit 1; }
git log --oneline -1
for nml in "$NML_HIST" "$NML_SSP"; do
    grep -q 'seasonalcycle = 0' "$nml" || { echo "ERROR: $nml lacks seasonalcycle=0 (stale namelist)"; exit 1; }
    grep -q 'x0 = -3040000.0' "$nml" || { echo "ERROR: $nml lacks AIS grid origin x0"; exit 1; }
done

echo "== 2. Build serial driver (guarantees binary matches source) =="
make clean
make gpdd_monthly_inout

echo "== 3. Pre-flight 1-yr runs (hist 1950, ssp370 2015) =="
rm -rf output_val1yr
mkdir -p output_val1yr
./gpdd_monthly_inout.x "$NML_HIST_VAL" > log_val_hist.txt 2>&1
./gpdd_monthly_inout.x "$NML_SSP_VAL" > log_val_ssp.txt 2>&1

echo "== 4. Pre-flight gate =="
cat > check_preflight.py <<'PYEOF'
import sys
import netCDF4 as nc

AREA = 761 * 761 * 64e6
SEC = 31556926.0

def gt(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:].mean().item() * AREA * SEC / 1e12

fails = 0

# Grid metadata (epsg:3031, origin -3040000) in both periods
for f in ["output_val1yr/acabf_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc",
          "output_val1yr/acabf_AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1_2015.nc"]:
    with nc.Dataset(f) as d:
        x = d.variables["x"][:]
        y = d.variables["y"][:]
        gi = getattr(d, "grid_information", "MISSING")
        ok = abs(x[0] + 3040000) < 1 and abs(y[0] + 3040000) < 1 and "3031" in gi
        print(f"{f.split('/')[-1]}: x {x[0]:.0f}..{x[-1]:.0f} "
              f"y {y[0]:.0f}..{y[-1]:.0f} {gi} {'OK' if ok else 'FAIL'}")
        if not ok:
            fails += 1

# 1950 must match validation run B (seasonalcycle=0, default physics,
# job 1807198): acabf 2329.1, mrro 344.4, prsn 2642.4, prra 31.0 Gt/yr
exp = {"acabf": 2329.1, "mrro": 344.4, "prsn": 2642.4, "prra": 31.0}
for v, e in exp.items():
    a = gt(f"output_val1yr/{v}_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc", v)
    ok = abs(a - e) < 1.0
    print(f"1950 {v:6s} {a:9.2f} Gt/yr (expected {e})  {'OK' if ok else 'FAIL'}")
    fails += 0 if ok else 1

# SMB closure acabf = pr - mrro and rain split in both pre-flights
for root, yr in [("AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1", 1950),
                 ("AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1", 2015)]:
    acabf = gt(f"output_val1yr/acabf_{root}_{yr}.nc", "acabf")
    mrro = gt(f"output_val1yr/mrro_{root}_{yr}.nc", "mrro")
    pr = gt(f"output_val1yr/pr_{root}_{yr}.nc", "pr")
    prsn = gt(f"output_val1yr/prsn_{root}_{yr}.nc", "prsn")
    prra = gt(f"output_val1yr/prra_{root}_{yr}.nc", "prra")
    err = abs(acabf + mrro - pr)
    print(f"{yr} closure |acabf+mrro-pr| = {err:.2e} Gt/yr; "
          f"pr {pr:.1f} = prsn {prsn:.1f} + prra {prra:.1f}")
    if err > 1e-3:
        fails += 1

sys.exit(1 if fails else 0)
PYEOF
$PY check_preflight.py

echo "== 5. Full runs (output to work storage via ./output symlink) =="
OUTTGT=$(readlink -f output)
case "$OUTTGT" in
  /cluster/work/users/heig/pdd/ais_08/output) rm -rf "$OUTTGT"/* ;;
  *) echo "ERROR: unexpected output target: $OUTTGT"; exit 1 ;;
esac
./gpdd_monthly_inout.x "$NML_HIST" > log_full_hist.txt 2>&1
echo "historical 1950-2014 done"
./gpdd_monthly_inout.x "$NML_SSP" > log_full_ssp.txt 2>&1
echo "ssp370 2015-2100 done"

echo "== 6. Verification of full output =="
cat > check_rerun.py <<'PYEOF'
import glob
import sys
import netCDF4 as nc

AREA = 761 * 761 * 64e6
SEC = 31556926.0
ROOTS = [("Historical", "AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1", 1950, 2014, 520),
         ("SSP370", "AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1", 2015, 2100, 688)]

def gt(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:].mean().item() * AREA * SEC / 1e12

fails = []
for tag, root, y0, y1, nexp in ROOTS:
    fs = glob.glob(f"output/*_{root}_*.nc")
    if len(fs) != nexp:
        fails.append(f"{tag} file count {len(fs)} != {nexp}")
        print(f"FAIL {tag}: {len(fs)} files, expected {nexp}")
        continue
    print(f"{tag}: {len(fs)} files OK")
    for y in range(y0, y1 + 1):
        acabf = gt(f"output/acabf_{root}_{y}.nc", "acabf")
        mrro = gt(f"output/mrro_{root}_{y}.nc", "mrro")
        pr = gt(f"output/pr_{root}_{y}.nc", "pr")
        err = abs(acabf + mrro - pr)
        if err > 1e-3:
            fails.append(f"{tag} {y} closure {err:.2e}")
        if y == y0 or y == y1 or (y - y0) % 10 == 0:
            print(f"  {y}  pr {pr:8.1f}  acabf {acabf:8.1f}  mrro {mrro:8.1f} Gt/yr")
    for y in (y0, y1):
        with nc.Dataset(f"output/acabf_{root}_{y}.nc") as d:
            x = d.variables["x"][:]
            yax = d.variables["y"][:]
            gi = getattr(d, "grid_information", "MISSING")
            if not (abs(x[0] + 3040000) < 1 and abs(yax[0] + 3040000) < 1
                    and "3031" in gi):
                fails.append(f"{tag} {y} grid metadata")
                print(f"FAIL {tag} {y} grid: {x[0]:.0f} {yax[0]:.0f} {gi}")

# Parameter-independent gate: pr(2014) must match the archive value
pr14 = gt("output/pr_AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_2014.nc", "pr")
print(f"pr(2014) = {pr14:.1f} Gt/yr (expected 2696.6)")
if abs(pr14 - 2696.6) > 1.0:
    fails.append(f"pr(2014) {pr14:.1f} != 2696.6")

if fails:
    print("RESULT: FAILURES:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
PYEOF
$PY check_rerun.py

echo "== DONE =="