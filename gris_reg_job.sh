#!/bin/bash
#SBATCH --job-name=gris_reg
#SBATCH --account=nn8002k
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --qos=preproc

set -o errexit
set -o nounset

D=/cluster/projects/nn5011k/heig/PDD/gris_16/pdd-model
PY=/cluster/projects/nn11016k/miniforge3/envs/nc/bin/python
cd "$D"

echo "== 1. Verify version (repo updated outside the job) =="
git merge-base --is-ancestor 90c77b8 HEAD || { echo "ERROR: missing 90c77b8"; exit 1; }
git log --oneline -1

echo "== 2. Compile serial binary =="
source setenv.sh
make gpdd_monthly_inout

echo "== 3. Regression gate (SMB + closure, 16km + 4km) =="
bash test_regression.sh

echo "== 4. Field-level bit-compat vs pre-change baseline (regtest/out16, out04) =="
# Re-run both tests keeping the output, then compare acabf/mrro fields.
for res in 16km 4km; do
    nml=params_CESM2_historical_anom_${res}.nml
    resm=$( [ "$res" = "16km" ] && echo 16000 || echo 4000 )
    tmpnml=/tmp/reg_field_${res}_$$.nml
    outdir=/tmp/reg_field_${res}_$$
    mkdir -p "$outdir"
    sed -e 's/nt\s*=\s*[0-9]*/nt = 1/' -e 's/outputvars\s*=\s*[0-9]*/outputvars = 0/' \
        -e "s|outpathname.*|outpathname = \"${outdir}\"|" "$nml" > "$tmpnml"
    ./gpdd_monthly_inout.x "$tmpnml" > /tmp/reg_field_${res}.log 2>&1
    rm -f "$tmpnml"
    $PY - "$res" "$outdir" "$resm" <<'EOF'
import sys, glob
import numpy as np
from netCDF4 import Dataset
res, outdir, resm = sys.argv[1], sys.argv[2], int(sys.argv[3])
setup = "out16" if res == "16km" else "out04"
base = f"/cluster/work/users/heig/pdd/regtest/{setup}"
for v in ("acabf", "mrro"):
    old = Dataset(f"{base}/{v}_GIS_NORCEPDD1_CESM2_Historical_r11i1p1f1_1950.nc").variables[v][:]
    new = np.concatenate([Dataset(f).variables[v][:] for f in sorted(glob.glob(f"{outdir}/{v}_*.nc"))])
    if old.shape != new.shape:
        print(f"  {res} {v}: SHAPE MISMATCH {old.shape} vs {new.shape}")
        continue
    dmax = float(np.max(np.abs(old - new)))
    rel = dmax / max(1e-30, float(np.max(np.abs(old))))
    print(f"  {res} {v}: max|new-old| = {dmax:.3e}  (rel {rel:.3e})  -> "
          + ("BIT-COMPAT" if dmax == 0.0 else ("near-identical" if rel < 1e-12 else "DIFFERS")))
EOF
    rm -rf "$outdir"
done

echo "== DONE =="
