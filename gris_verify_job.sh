#!/bin/bash
#SBATCH --job-name=gris_verify
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#
# Verification-only job for the ISMIP7 GrIS corrected rerun (run after
# gris_rerun_job.sh). Checks, per period (Historical/SSP370):
#   1. File counts: 9 vars x 65 yr (hist) / 9 vars x 86 yr (ssp370)
#   2. Cellwise SMB closure acabf = pr - mrro (every year at 8 km, every
#      5th year at 1 km)
#   3. Grid metadata: epsg:3413, origin (-720000, -3450000)
#   4. Archive ground truth at 1950/2014/2015/2100: tas/pr/prsn/prra
#      bit-exact, acabf identical (std only), mrro = arch mrro + arch prra
#      (the fix identity); snicefreez/sir reported informatively
#   5. netCDF4 compression active in the new files
# Comparison files are STAGED under /cluster/work/users/heig/pdd/gris_$RES/
# archive_cmp(_prMAR) because the NS8002K NIRD mount is NOT available on
# compute nodes.
# Env: RES=08|01  VARIANT=std|prMAR   (defaults: 08/std)
# Usage: sbatch --export=RES=08,VARIANT=std gris_verify_job.sh

set -o errexit
set -o nounset

RES=${RES:-08}
VARIANT=${VARIANT:-std}
case "$RES/$VARIANT" in
  08/std)     OUTLINK=output_ismip7;;
  08/prMAR)   OUTLINK=output_ismip7_prMAR;;
  01/std)     OUTLINK=output_ismip7;;
  *) echo "ERROR: bad RES/VARIANT: $RES/$VARIANT"; exit 1;;
esac

RUNDIR=/cluster/projects/nn5011k/heig/PDD/gris_$RES/pdd-model
cd "$RUNDIR"
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python
export RES VARIANT OUTLINK

cat > check_verify_gris.py <<'PYEOF'
import glob
import os
import sys
import netCDF4 as nc

RES, VARIANT, OUT = os.environ["RES"], os.environ["VARIANT"], os.environ["OUTLINK"]
STD = VARIANT == "std"
res_m = 8000.0 if RES == "08" else 1000.0
X0, Y0 = -720000.0, -3450000.0
Gt = res_m * res_m * 31556926.0 / 1e12
STEP = 1 if RES == "08" else 5
ARCH = "/cluster/work/users/heig/pdd/gris_%s/archive_cmp%s" % (
    RES, "" if STD else "_prMAR")
ROOTS = [("Historical", "GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1", 1950, 2014, 585),
         ("SSP370", "GIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1", 2015, 2100, 774)]
CMP_YEARS = {1950, 2014, 2015, 2100}


def arr(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:]


fails = []
for tag, root, y0, y1, nexp in ROOTS:
    fs = glob.glob(f"{OUT}/{tag}/*_{root}_*.nc")
    if len(fs) != nexp:
        fails.append(f"{tag} file count {len(fs)} != {nexp}")
        print(f"FAIL {tag}: {len(fs)} files, expected {nexp}")
        continue
    print(f"{tag}: {len(fs)} files OK (9 vars x {y1-y0+1} yr)")
    for y in range(y0, y1 + 1, STEP):
        a = arr(f"{OUT}/{tag}/acabf_{root}_{y}.nc", "acabf")
        m = arr(f"{OUT}/{tag}/mrro_{root}_{y}.nc", "mrro")
        p = arr(f"{OUT}/{tag}/pr_{root}_{y}.nc", "pr")
        err = abs(a + m - p).max()
        if err > 1e-9:
            fails.append(f"{tag} {y} closure {err:.2e}")
        if y == y0 or y == y1 or (y - y0) % (10 * STEP) == 0:
            print(f"  {y}  pr {p.mean()*Gt:8.1f}  acabf {a.mean()*Gt:8.1f}  "
                  f"mrro {m.mean()*Gt:8.1f}  sir "
                  f"{arr(f'{OUT}/{tag}/sir_{root}_{y}.nc', 'sir').mean()*Gt:6.2f} Gt/yr")
        if arr(f"{OUT}/{tag}/sir_{root}_{y}.nc", "sir").min() < 0:
            fails.append(f"{tag} {y} negative sir")
    for y in (y0, y1):
        with nc.Dataset(f"{OUT}/{tag}/acabf_{root}_{y}.nc") as d:
            x, yax = d.variables["x"][:], d.variables["y"][:]
            gi = getattr(d, "grid_information", "MISSING")
        if not (abs(x[0] - X0) < 1 and abs(yax[0] - Y0) < 1 and "3413" in gi):
            fails.append(f"{tag} {y} grid metadata")
            print(f"FAIL {tag} {y} grid: {x[0]:.0f} {yax[0]:.0f} {gi}")

print("Archive comparison (field max |rerun - archive|):")
for tag, root, y0, y1, _ in ROOTS:
    for y in sorted(CMP_YEARS & set(range(y0, y1 + 1))):
        new = lambda v: arr(f"{OUT}/{tag}/{v}_{root}_{y}.nc", v)
        arc = lambda v: arr(f"{ARCH}/{tag}/{v}_{root}_{y}.nc", v)
        for v in ["tas", "pr", "prsn", "prra"]:
            dmax = abs(new(v) - arc(v)).max()
            if dmax != 0.0:
                fails.append(f"{v} {y} archive diff {dmax:.2e}")
            print(f"  {v:5s} {y}  max|diff| = {dmax:.3e} "
                  f"{'OK' if dmax == 0.0 else 'FAIL'}")
        if STD:
            dmax = abs(new("acabf") - arc("acabf")).max()
            if dmax >= 1e-9:
                fails.append(f"acabf {y} archive diff {dmax:.2e}")
            print(f"  acabf {y}  max|diff| = {dmax:.3e} "
                  f"{'OK' if dmax < 1e-9 else 'FAIL'}")
        dmax = abs(new("mrro") - (arc("mrro") + arc("prra"))).max()
        if dmax >= 1e-9:
            fails.append(f"mrro identity {y} {dmax:.2e}")
        print(f"  mrro identity {y}  max|diff| = {dmax:.3e} "
              f"{'OK' if dmax < 1e-9 else 'FAIL'}")

# compression check: new files must be netCDF4/zlib
f0 = f"{OUT}/Historical/{'acabf'}_{ROOTS[0][1]}_{1950}.nc"
with nc.Dataset(f0) as d:
    flt = d.variables["acabf"].filters()
    print(f"compression: {d.data_model}, zlib={flt.get('zlib')}, "
          f"complevel={flt.get('complevel')}")
if d.data_model != "NETCDF4" or not flt.get("zlib"):
    fails.append("new output not compressed netCDF4")

if fails:
    print("RESULT: FAILURES:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
PYEOF
$PY check_verify_gris.py

echo "== DONE =="
