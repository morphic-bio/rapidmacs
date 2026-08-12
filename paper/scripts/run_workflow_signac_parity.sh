#!/usr/bin/env bash
# Full-file Signac CallPeaks default parity check. This script records the
# Phase 1 commands; it is not part of `make test` and must be invoked manually.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MACS3="${MACS3:-${REPO_ROOT}/.macs3-reference/macs3-3.0.3-venv/bin/macs3}"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
INPUT="${INPUT:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k_100k/fragments.tsv}"
OUT_ROOT="${OUT_ROOT:-${REPO_ROOT}/build/workflow_compatibility/$(date -u +%Y%m%dT%H%M%SZ)_signac}"

mkdir -p "${OUT_ROOT}/macs3" "${OUT_ROOT}/rapidmacs"
COMMANDS="${OUT_ROOT}/commands.sh"
printf '#!/usr/bin/env bash\nset -euo pipefail\n\n' >"${COMMANDS}"

record_command() {
  printf '%q ' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

macs_cmd=("${MACS3}" callpeak -t "${INPUT}" -f BED -g 2700000000
          -n signac_atac --outdir "${OUT_ROOT}/macs3" -q 0.05 --nomodel
          --extsize 200 --shift -100 --llocal 10000 --keep-dup 1)
rapid_cmd=("${RAPIDMACS}" --macs-profile signac-atac
           --macs-input-format BED --macs-bed-input "${INPUT}"
           --macs-bed-peaks-out "${OUT_ROOT}/rapidmacs/signac_atac_peaks.narrowPeak"
           --macs-bed-summits-out "${OUT_ROOT}/rapidmacs/signac_atac_summits.bed"
           --macs-bed-summary "${OUT_ROOT}/rapidmacs/summary.tsv"
           --outdir "${OUT_ROOT}/rapidmacs" -n signac_atac)

record_command "${macs_cmd[@]}"
"${macs_cmd[@]}" >"${OUT_ROOT}/macs3.stdout" 2>"${OUT_ROOT}/macs3.stderr"
record_command "${rapid_cmd[@]}"
"${rapid_cmd[@]}" >"${OUT_ROOT}/rapidmacs.stdout" 2>"${OUT_ROOT}/rapidmacs.stderr"

macs_np="${OUT_ROOT}/macs3/signac_atac_peaks.narrowPeak"
macs_sm="${OUT_ROOT}/macs3/signac_atac_summits.bed"
rapid_np="${OUT_ROOT}/rapidmacs/signac_atac_peaks.narrowPeak"
rapid_sm="${OUT_ROOT}/rapidmacs/signac_atac_summits.bed"
cmp "${macs_np}" "${rapid_np}"
cmp "${macs_sm}" "${rapid_sm}"
md5sum "${macs_np}" "${rapid_np}" "${macs_sm}" "${rapid_sm}"
echo "PASS: Signac full-file narrowPeak and summits are byte-identical"
echo "Commands: ${COMMANDS}"
