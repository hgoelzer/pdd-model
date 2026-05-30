#!/bin/bash
# Create work-storage directories and replace Forcing/ and output/ with symlinks.
# Run once after cloning/setting up a new configuration directory.
#
# The setup name is derived from the parent directory of the repo:
#   .../gris_16/pdd-model  →  setup = gris_16
#
# Usage:
#   ./setup_work.sh                  # auto-derive work path from setup name
#   ./setup_work.sh /path/on/work    # explicit work base for this setup

set -euo pipefail

REPO_DIR=$(cd "$(dirname "$0")" && pwd)
SETUP=$(basename "$(dirname "${REPO_DIR}")")
WORK_BASE="/cluster/work/users/${USER}/pdd"
WORK_DIR="${1:-${WORK_BASE}/${SETUP}}"

echo "=== setup_work.sh ==="
echo "  Setup:    ${SETUP}"
echo "  Work dir: ${WORK_DIR}"
echo ""

mkdir -p "${WORK_DIR}/Forcing"
mkdir -p "${WORK_DIR}/output"

for dir in Forcing output; do
    src="${REPO_DIR}/${dir}"
    dst="${WORK_DIR}/${dir}"

    if [ -L "${src}" ]; then
        echo "  ${dir}: already a symlink → $(readlink "${src}") — skipping"
        continue
    fi

    if [ -d "${src}" ] && [ -n "$(ls -A "${src}" 2>/dev/null)" ]; then
        echo "  ${dir}: moving contents to ${dst} ..."
        mv "${src}"/* "${dst}"/
    fi

    [ -d "${src}" ] && rmdir "${src}"
    ln -s "${dst}" "${src}"
    echo "  ${dir}: ${src} → ${dst}"
done

echo ""
echo "Done. Verify:"
echo "  ls -la ${REPO_DIR}/"
