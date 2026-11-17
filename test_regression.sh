#!/bin/bash
# Regression test: run the serial model for 1 year at 16km and 4km,
# then check integrated SMB against known-good values from Report_serial.txt.
#
# Usage (from pdd-model directory):
#   bash test_regression.sh                 # run both resolutions
#   bash test_regression.sh 16km            # run only 16km
#   RUN_DIR=/some/other/path bash test_regression.sh
#
# Requires:
#   - gpdd_monthly_inout.x (serial binary, built with make gpdd_monthly_inout)
#   - CDO and Python with numpy in PATH / conda env nc
#   - Forcing files for CESM2 historical 1950 at each resolution
#
# On a full queue, run interactively via the preproc QOS (from the repo dir,
# forcing staged under ./Forcing):
#   salloc -A nn5011k --qos=preproc --ntasks=1 --mem=4G -t 0:30:00 \
#          bash test_regression.sh
# Note: namelists are temporarily patched to nt=1 and outputvars=0 (full
# output, needed for the acabf = pr - mrro closure check).

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
NFAIL=0

# ---------- helpers -----------------------------------------------------------

log()  { echo "  $*"; }
pass() { echo "  PASS  $1"; }
fail() { echo "  FAIL  $1  (actual=$2, expected≈$3, tol=$4)"; NFAIL=$((NFAIL+1)); }

# Compute integrated SMB from acabf NetCDF files via CDO + Python.
# conda run does not forward stdin on the clusters, so the Python code is
# written to a temp file first. Usage: smb=$(compute_smb <outdir> <res_m>)
compute_smb() {
    local outdir=$1 res_m=$2
    local pyscript="/tmp/test_regression_smb_$$.py"
    cat > "${pyscript}" <<EOF
import glob, subprocess, numpy as np
files = sorted(glob.glob("${outdir}/acabf_*.nc"))
if not files:
    print("NaN"); exit(0)
cell_area = ${res_m} ** 2
SPM = 30 * 86400
annual = []
for f in files:
    r = subprocess.run(["cdo","-s","output","-fldsum",f], capture_output=True, text=True, check=True)
    vals = list(map(float, r.stdout.split()))
    annual.append(sum(vals) * cell_area * SPM / 1e12)
print("%.2f" % (sum(annual)/len(annual)))
EOF
    conda run -n nc python3 "${pyscript}" || echo "NaN"
    rm -f "${pyscript}"
}

# Run a 1-year test at a given resolution.
# Usage: run_test <setup_name> <nml_file> <res_m> <expected_smb> <tol>
run_test() {
    local setup=$1 nml=$2 res_m=$3 expected=$4 tol=$5
    local work_base="${WORK_BASE:-/cluster/work/users/$USER/pdd}"
    local outdir="${work_base}/${setup}/output"
    local binary="${REPO_DIR}/gpdd_monthly_inout.x"

    echo ""
    echo "--- Test: ${setup} (${res_m} m, 1 year 1950) ---"

    if [ ! -f "${binary}" ]; then
        echo "  SKIP  binary not found: ${binary}"
        return
    fi
    if [ ! -f "${REPO_DIR}/${nml}" ]; then
        echo "  SKIP  namelist not found: ${nml}"
        return
    fi

    # Temporarily set nt=1 and outputvars=0 (full output) by patching a temp namelist
    local tmpnml="/tmp/test_regression_$$.nml"
    sed -e 's/nt\s*=\s*[0-9]*/nt = 1/' -e 's/outputvars\s*=\s*[0-9]*/outputvars = 0/' "${REPO_DIR}/${nml}" > "${tmpnml}"
    # Also redirect output to a temporary directory
    local tmpout="/tmp/test_regression_out_$$"
    mkdir -p "${tmpout}"
    sed -i "s|outpathname.*=.*|outpathname = \"${tmpout}\"|" "${tmpnml}"

    log "running model for 1 year..."
    # Cluster runtime env (netCDF-Fortran lib path); no-op where `module` is absent
    if command -v module > /dev/null 2>&1; then
        # shellcheck disable=SC1091
        source "${REPO_DIR}/setenv.sh"
    fi
    (cd "${REPO_DIR}" && \
        DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}" \
        ./gpdd_monthly_inout.x "${tmpnml}" > /tmp/test_regression_run_$$.log 2>&1) || {
        echo "  FAIL  model run exited non-zero — see /tmp/test_regression_run_$$.log"
        NFAIL=$((NFAIL+1))
        rm -f "${tmpnml}"; rm -rf "${tmpout}"
        return
    }

    log "computing integrated SMB..."
    local smb
    smb=$(compute_smb "${tmpout}" "${res_m}")

    log "checking SMB closure acabf = pr - mrro..."
    local closure
    closure=$(check_closure "${tmpout}")
    if [ "${closure}" = "ok" ]; then
        pass "${setup}: SMB closure acabf = pr - mrro (cellwise)"
    elif [ "${closure}" = "NaN" ]; then
        echo "  FAIL  ${setup}: could not read output files for closure check"
        NFAIL=$((NFAIL+1))
    else
        fail "${setup}: SMB closure acabf = pr - mrro" "${closure}" "max|err| <= 1e-6" "1e-6"
    fi

    rm -f "${tmpnml}"; rm -rf "${tmpout}"

    if [ "${smb}" = "NaN" ]; then
        echo "  FAIL  could not compute SMB (no output files?)"
        NFAIL=$((NFAIL+1))
        return
    fi

    # Check within tolerance
    local ok
    ok=$(python3 -c "print('ok' if abs(${smb} - ${expected}) <= ${tol} else 'fail')" 2>/dev/null)
    if [ "${ok}" = "ok" ]; then
        pass "${setup}: SMB=${smb} Gt/yr  (expected≈${expected} ± ${tol})"
    else
        fail "${setup}: SMB check" "${smb}" "${expected}" "${tol}"
    fi
}

# Check the SMB closure acabf = pr - mrro (mrro = total runoff = melt + rain).
# Cellwise, tolerance is float-rounding level. Usage: check_closure <outdir>
check_closure() {
    local outdir=$1
    local pyscript="/tmp/test_regression_closure_$$.py"
    cat > "${pyscript}" <<'PYEOF'
import glob, sys
import numpy as np
from netCDF4 import Dataset
d = sys.argv[1]
def rd(pat):
    fs = sorted(glob.glob(f"{d}/{pat}"))
    return np.concatenate([Dataset(f).variables[pat.split("_")[0]][:] for f in fs]) if fs else None
acabf, mrro, pr = rd("acabf_*.nc"), rd("mrro_*.nc"), rd("pr_*.nc")
if acabf is None or mrro is None or pr is None or acabf.shape != mrro.shape:
    print("NaN"); sys.exit(0)
err = float(np.max(np.abs(acabf + mrro - pr)))
print("ok" if err <= 1e-6 * max(1.0, float(np.max(np.abs(pr)))) else f"fail {err:.2e}")
PYEOF
    conda run -n nc python3 "${pyscript}" "${outdir}" 2>/dev/null || echo "NaN"
    rm -f "${pyscript}"
}

# ---------- test definitions --------------------------------------------------

FILTER="${1:-all}"

echo "=== test_regression.sh: serial model correctness ==="
echo "    Baselines from Report_serial.txt (CESM2 historical 1950-1954 mean):"
echo "    gris_16: +198 Gt/yr (1950 alone)   gris_04: +200 Gt/yr (1950 alone)"
echo ""

# 16km: expected 1950 value = +198 Gt/yr (from Report_serial.txt annual[0])
# Use ±20 Gt/yr tolerance (inter-annual variability is ~70 Gt/yr std, 1950 is the min year)
if [[ "${FILTER}" == "all" || "${FILTER}" == "16km" ]]; then
    run_test "gris_16" "params_CESM2_historical_anom_16km.nml" "16000" "198" "20"
fi

# 4km: expected 1950 value = +200 Gt/yr
if [[ "${FILTER}" == "all" || "${FILTER}" == "4km" ]]; then
    run_test "gris_04" "params_CESM2_historical_anom_04km.nml" "4000" "200" "20"
fi

# ---------- result ------------------------------------------------------------
echo ""
if [ ${NFAIL} -eq 0 ]; then
    echo "ALL REGRESSION TESTS PASSED"
    exit 0
else
    echo "FAILED: ${NFAIL} test(s)"
    exit 1
fi
