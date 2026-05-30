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

set -euo pipefail

REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
NFAIL=0

# ---------- helpers -----------------------------------------------------------

log()  { echo "  $*"; }
pass() { echo "  PASS  $1"; }
fail() { echo "  FAIL  $1  (actual=$2, expected≈$3, tol=$4)"; NFAIL=$((NFAIL+1)); }

# Compute integrated SMB from acabf NetCDF files via CDO + Python.
# Usage: smb=$(compute_smb <outdir> <res_m>)
compute_smb() {
    local outdir=$1 res_m=$2
    conda run -n nc python3 - <<EOF
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

    # Temporarily set nt=1 by patching a temp namelist
    local tmpnml="/tmp/test_regression_$$.nml"
    sed 's/nt\s*=\s*[0-9]*/nt = 1/' "${REPO_DIR}/${nml}" > "${tmpnml}"
    # Also redirect output to a temporary directory
    local tmpout="/tmp/test_regression_out_$$"
    mkdir -p "${tmpout}"
    sed -i "s|outpathname.*=.*|outpathname = \"${tmpout}\"|" "${tmpnml}"

    log "running model for 1 year..."
    (cd "${REPO_DIR}" && \
        DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:-}" \
        ./gpdd_monthly_inout.x "${tmpnml}" > /tmp/test_regression_run_$$.log 2>&1) || {
        echo "  FAIL  model run exited non-zero — see /tmp/test_regression_run_$$.log"
        NFAIL=$((NFAIL+1))
        rm -f "${tmpnml}"; rm -rf "${tmpout}"
        return
    }

    log "computing integrated SMB..."
    # Inline Python to compute SMB from tmpout
    local smb
    smb=$(conda run -n nc python3 - <<EOF 2>/dev/null || echo "NaN"
import glob, subprocess, numpy as np
files = sorted(glob.glob("${tmpout}/acabf_*.nc"))
if not files:
    print("NaN"); exit(0)
cell_area = ${res_m} ** 2
SPM = 30 * 86400
annual = [sum(map(float, subprocess.run(["cdo","-s","output","-fldsum",f],
    capture_output=True, text=True, check=True).stdout.split()))*cell_area*SPM/1e12
    for f in files]
print("%.2f" % (sum(annual)/len(annual)))
EOF
)

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
