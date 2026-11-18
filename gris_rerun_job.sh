#!/bin/bash
#SBATCH --job-name=gris_rerun
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8G
#
# ISMIP7 GrIS corrected rerun 2026-11: identical physics to the SMBMIP
# submission (0.00297/0.00791/4.5/1.0, pmax=0.3, taj cycle), driver fixes
# only (mrro = total runoff; snicefreez = actually retained meltwater).
# 1950-2014 + ssp370 2015-2100, all 9 vars, monthly. Serial driver.
# Env: RES=08|01  VARIANT=std|prMAR   (defaults: 08/std)
# Usage:
#   sbatch --export=RES=08,VARIANT=std   gris_rerun_job.sh
#   sbatch --export=RES=08,VARIANT=prMAR gris_rerun_job.sh
#   sbatch --qos=normal --time=08:00:00 --export=RES=01,VARIANT=std gris_rerun_job.sh

set -o errexit
set -o nounset

RES=${RES:-08}
VARIANT=${VARIANT:-std}
case "$RES/$VARIANT" in
  08/std)     NMLSUF=_08km_ismip7;        OUTLINK=output_ismip7;;
  08/prMAR)   NMLSUF=_08km_ismip7_prMAR;  OUTLINK=output_ismip7_prMAR;;
  01/std)     NMLSUF=_01km_ismip7;        OUTLINK=output_ismip7;;
  *) echo "ERROR: bad RES/VARIANT: $RES/$VARIANT"; exit 1;;
esac
NML_HIST=params_CESM2_historical_anom${NMLSUF}.nml
NML_SSP=params_CESM2_ssp370_anom${NMLSUF}.nml

RUNDIR=/cluster/projects/nn5011k/heig/PDD/gris_$RES/pdd-model
cd "$RUNDIR"
source setenv.sh
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python
export RES VARIANT

echo "== 1. Verify code version and namelists =="
git merge-base --is-ancestor 810b265 HEAD || { echo "ERROR: run dir missing commit 810b265 (mrro/sir fix) - update the clone first"; exit 1; }
git log --oneline -1
for nml in "$NML_HIST" "$NML_SSP"; do
    grep -q 'seasonalcycle = 1' "$nml" || { echo "ERROR: $nml lacks seasonalcycle=1 (stale namelist)"; exit 1; }
    grep -q 'x0 = -720000.0' "$nml" || { echo "ERROR: $nml lacks GrIS grid origin x0"; exit 1; }
    grep -q "res_suffix = \"i${RES}000m\"" "$nml" || { echo "ERROR: $nml wrong resolution"; exit 1; }
    grep -q 'pmax = 0.3' "$nml" || { echo "ERROR: $nml lacks pmax=0.3"; exit 1; }
done
case "$VARIANT" in
  std)   grep -q 'MARv3.14LC26/pr_ref' "$NML_HIST" || { echo "ERROR: std nml not using LC26 pr_ref"; exit 1; };;
  prMAR) grep -q 'pr_ref_MARv3.14-oli' "$NML_HIST" && grep -q 'sftgif_BM5' "$NML_HIST" || { echo "ERROR: prMAR nml not using oli raw pr_ref + BM5 mask"; exit 1; };;
esac

echo "== 2. Build serial driver (guarantees binary matches source) =="
mkdir -p obj
make clean
make gpdd_monthly_inout

echo "== 3. Output symlinks (work storage) =="
WORK=/cluster/work/users/heig/pdd/gris_$RES
mkdir -p "$WORK/$OUTLINK/Historical" "$WORK/$OUTLINK/SSP370" "$WORK/output_val1yr"
ln -sfn "$WORK/$OUTLINK" "$OUTLINK"
ln -sfn "$WORK/output_val1yr" output_val1yr

echo "== 4. Pre-flight 1-yr runs (hist 1950, ssp370 2015) =="
sed -E -e 's/^  nt[[:space:]]*= [0-9]+/  nt           = 1/' \
       -e 's#^  outpathname = .*#  outpathname = "./output_val1yr"#' \
       "$NML_HIST" > params_val_hist.nml
sed -E -e 's/^  nt[[:space:]]*= [0-9]+/  nt           = 1/' \
       -e 's#^  outpathname = .*#  outpathname = "./output_val1yr"#' \
       "$NML_SSP" > params_val_ssp.nml
rm -f "$WORK"/output_val1yr/*.nc
./gpdd_monthly_inout.x params_val_hist.nml > log_val_hist.txt 2>&1
./gpdd_monthly_inout.x params_val_ssp.nml > log_val_ssp.txt 2>&1

echo "== 5. Pre-flight gate =="
cat > check_preflight_gris.py <<'PYEOF'
import os
import sys
import netCDF4 as nc

RES, VARIANT = os.environ["RES"], os.environ["VARIANT"]
STD = VARIANT == "std"
res_m = 8000.0 if RES == "08" else 1000.0
NX, NY = (211, 361) if RES == "08" else (1681, 2881)
X0, Y0 = -720000.0, -3450000.0
Gt = res_m * res_m * NX * NY * 31556926.0 / 1e12
ARCH = "/cluster/work/users/heig/pdd/gris_%s/archive_cmp%s" % (
    RES, "" if STD else "_prMAR")
ROOTS = [("Historical", "GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1", 1950),
         ("SSP370", "GIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1", 2015)]


def arr(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:]


fails = 0
for tag, root, yr in ROOTS:
    new = lambda v: arr(f"output_val1yr/{v}_{root}_{yr}.nc", v)
    arc = lambda v: arr(f"{ARCH}/{tag}/{v}_{root}_{yr}.nc", v)
    with nc.Dataset(f"output_val1yr/acabf_{root}_{yr}.nc") as d:
        x, y = d.variables["x"][:], d.variables["y"][:]
        gi = getattr(d, "grid_information", "MISSING")
    ok = abs(x[0] - X0) < 1 and abs(y[0] - Y0) < 1 and "3413" in gi
    print(f"{tag} {yr} grid: x {x[0]:.0f} y {y[0]:.0f} {gi} "
          f"{'OK' if ok else 'FAIL'}")
    fails += 0 if ok else 1
    # tas/pr must match the archive bit-exactly; prsn/prra/acabf carry only
    # float-reassociation noise from the config refactor (<= ~3e-7 observed)
    for v in ["tas", "pr"]:
        dmax = abs(new(v) - arc(v)).max()
        ok = dmax == 0.0
        print(f"{tag} {yr} {v:5s}: max|new-arch| = {dmax:.3e} "
              f"{'OK' if ok else 'FAIL'}")
        fails += 0 if ok else 1
    for v in ["prsn", "prra", "acabf"]:
        dlt = new(v) - arc(v)
        dmax, dmean = abs(dlt).max(), abs(dlt.mean())
        ok = dmax < 1e-5 and dmean < 1e-11
        print(f"{tag} {yr} {v:5s}: max|new-arch| = {dmax:.3e} "
              f"mean = {dmean:.1e} {'OK' if ok else 'FAIL'}")
        fails += 0 if ok else 1
    # mrro fix identity: new mrro (melt + rain) = archive mrro (melt only)
    # + archive prra; float-reassociation noise grows on the 1 km grid
    dlt = new("mrro") - (arc("mrro") + arc("prra"))
    dmax, dmean = abs(dlt).max(), abs(dlt.mean())
    ok = dmax < 1e-5 and dmean < 1e-11
    print(f"{tag} {yr} mrro identity (new = arch + arch prra): "
          f"max|diff| = {dmax:.3e} mean = {dmean:.1e} "
          f"{'OK' if ok else 'FAIL'}")
    fails += 0 if ok else 1
    dmax = abs(new("acabf") + new("mrro") - new("pr")).max()
    ok = dmax < 1e-9
    print(f"{tag} {yr} closure |acabf+mrro-pr|: {dmax:.3e} "
          f"{'OK' if ok else 'FAIL'}")
    fails += 0 if ok else 1
    sf = new("snicefreez")
    print(f"{tag} {yr} means Gt/yr: acabf {new('acabf').mean()*Gt:8.1f} "
          f"mrro {new('mrro').mean()*Gt:8.1f} "
          f"snicefreez {sf.mean()*Gt:6.1f} "
          f"(archive {arc('snicefreez').mean()*Gt:6.1f}); "
          f"sf min {sf.min():.2e}, {(sf < 0).sum()} neg cells")
    if sf.min() < -1e-6:
        print(f"FAIL {tag} {yr}: snicefreez min below -1e-6")
        fails += 1

sys.exit(1 if fails else 0)
PYEOF
$PY check_preflight_gris.py

echo "== 6. Full runs =="
./gpdd_monthly_inout.x "$NML_HIST" > log_full_hist.txt 2>&1
echo "historical 1950-2014 done"
./gpdd_monthly_inout.x "$NML_SSP" > log_full_ssp.txt 2>&1
echo "ssp370 2015-2100 done"

echo "== 7. File counts =="
for tag in Historical SSP370; do
    if [ "$tag" = Historical ]; then nexp=585; else nexp=774; fi
    n=$(ls "$OUTLINK/$tag" | wc -l)
    echo "$OUTLINK/$tag: $n files (expected $nexp)"
done

echo "== DONE =="

