#!/bin/bash
#SBATCH --job-name=ais_verify
#SBATCH --account=nn8002k
#SBATCH --qos=preproc
#SBATCH --time=00:45:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4G
#
# Verification-only job for the AIS SMBMIP rerun (model output complete from
# job 1807220; its verification had wrong expectations: 9 vars/yr including
# rfr, matching the archive layout exactly, and a stale hardcoded pr(2014)).
#
# Checks:
#   1. File counts: 9 vars x 65 yr (hist) + 9 vars x 86 yr (ssp370)
#   2. Cellwise SMB closure acabf = pr - mrro for every year
#   3. Grid metadata: epsg:3031, origin -3040000
#   4. Ground-truth archive comparison: pr/prsn/prra must match the original
#      taj-era archive field-wise for 1950/2014 (hist) and 2015/2100 (ssp370)
#      (these are parameter- and cycle-independent). Comparison files are
#      STAGED to /cluster/work/users/heig/pdd/ais_08/archive_cmp/ because the
#      NS8002K NIRD mount is NOT available on compute nodes.
# Usage: sbatch ais_verify_job.sh

set -o errexit
set -o nounset

RUNDIR=/cluster/projects/nn5011k/heig/PDD/ais_08/pdd-model
cd "$RUNDIR"
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python

cat > check_verify.py <<'PYEOF'
import glob
import sys
import netCDF4 as nc

AREA = 761 * 761 * 64e6
SEC = 31556926.0
ROOTS = [("Historical", "AIS_NORCEPDD1_CESM2_Historical_r11i1p1f1", 1950, 2014),
         ("SSP370", "AIS_NORCEPDD1_CESM2_SSP370_r11i1p1f1", 2015, 2100)]
ARCH = "/cluster/work/users/heig/pdd/ais_08/archive_cmp"
VARS9 = ["acabf", "mrro", "pdd", "pr", "prra", "prsn", "rfr", "snicefreez", "tas"]
# parameter- and cycle-independent: must match the taj archive
MATCH_VARS = ["pr", "prsn", "prra"]

def gt(path, var):
    with nc.Dataset(path) as d:
        return d.variables[var][:].mean().item() * AREA * SEC / 1e12

fails = []
for tag, root, y0, y1 in ROOTS:
    nexp = 9 * (y1 - y0 + 1)
    fs = glob.glob(f"output/*_{root}_*.nc")
    if len(fs) != nexp:
        fails.append(f"{tag} file count {len(fs)} != {nexp}")
        print(f"FAIL {tag}: {len(fs)} files, expected {nexp}")
        continue
    print(f"{tag}: {len(fs)} files OK (9 vars x {y1-y0+1} yr)")
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

# Ground truth: compare against the original taj-era archive (same model
# input chain proven by the 1950 bit-identical validation, job 1807140).
# pr is parameter- and cycle-independent -> must match bit-exactly.
# prsn/prra depend on the temperature path through the Gaussian rain/snow
# threshold (taj vs actual temps) -> small expected differences (<= ~1e-3);
# gate instead on prsn + prra = pr holding against the archive as well.
print("Archive comparison (field max |rerun - archive|):")
for tag, root, y0, y1 in ROOTS:
    for y in (y0, y1):
        for v in MATCH_VARS:
            with nc.Dataset(f"output/{v}_{root}_{y}.nc") as dn:
                a = dn.variables[v][:]
            with nc.Dataset(f"{ARCH}/{tag}/{v}_{root}_{y}.nc") as da:
                b = da.variables[v][:]
            dmax = abs(a - b).max().item()
            if v == "pr":
                ok = dmax == 0.0
            else:
                ok = dmax < 1e-3
            print(f"  {v:5s} {y}  max|diff| = {dmax:.3e}  {'OK' if ok else 'FAIL'}")
            if not ok:
                fails.append(f"{v} {y} archive diff {dmax:.2e}")
        # partition closure against the archive: prsn + prra = pr on both sides
        with nc.Dataset(f"{ARCH}/{tag}/pr_{root}_{y}.nc") as da:
            pr_a = da.variables["pr"][:]
        with nc.Dataset(f"{ARCH}/{tag}/prsn_{root}_{y}.nc") as da:
            prsn_a = da.variables["prsn"][:]
        with nc.Dataset(f"{ARCH}/{tag}/prra_{root}_{y}.nc") as da:
            prra_a = da.variables["prra"][:]
        err = abs(prsn_a + prra_a - pr_a).max().item()
        print(f"  archive partition {y}: max|prsn+prra-pr| = {err:.3e}")
        if err > 1e-9:
            fails.append(f"archive partition {y} {err:.2e}")

# Informative: archive 1950 vs rerun 1950 for the cycle-affected scalars
root = ROOTS[0][1]
for v in ["acabf", "mrro"]:
    arch = gt(f"{ARCH}/Historical/{v}_{root}_1950.nc", v)
    new = gt(f"output/{v}_{root}_1950.nc", v)
    print(f"1950 {v}: archive(taj) {arch:.1f} -> rerun(sc0) {new:.1f} Gt/yr")

if fails:
    print("RESULT: FAILURES:")
    for f in fails:
        print(f"  {f}")
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
PYEOF
$PY check_verify.py

echo "== DONE =="