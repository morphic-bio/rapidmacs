#!/usr/bin/env bash
# Copy the controlled ChIP-seq and CUT&RUN BAMPE inputs to the SSD.

set -euo pipefail

STRICT_HOST="${STRICT_HOST:-0}"
guard() {
  [[ "$1" == 0 ]] || return 0
  if [[ "${STRICT_HOST}" == 1 ]]; then echo "ERROR: $2" >&2; exit 2; fi
  echo "WARNING: $2 (timings will not be comparable to the paper)" >&2
}


# ---------------------------------------------------------------------------
# Paths are environment-driven so this runs anywhere. Defaults are relative to
# the repository, so nothing outside it is touched unless you point these
# elsewhere. See paper/README.md.
#   WORK_ROOT      scratch space for run outputs   (default ./work)
#   DATA_DIR       input datasets                  (default ./data)
#   ARCHIVE_ROOT   where completed runs are copied (default ./archive)
#   MACS3_REF_ROOT MACS3 reference installations   (default ./.macs3-reference)
# ---------------------------------------------------------------------------
WORK_ROOT="${WORK_ROOT:-$(pwd)/work}"
DATA_DIR="${DATA_DIR:-$(pwd)/data}"
ARCHIVE_ROOT="${ARCHIVE_ROOT:-$(pwd)/archive}"
MACS3_REF_ROOT="${MACS3_REF_ROOT:-$(pwd)/.macs3-reference}"
CHROMAP_SUITE="${CHROMAP_SUITE:-$(pwd)/chromap-suite}"
CELLRANGER_REF="${CELLRANGER_REF:-${DATA_DIR}/refdata-cellranger-arc-GRCh38-2020-A-2.0.0}"


RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/paper_benchmarks}"
INPUT_ROOT="${RUN_ROOT}/inputs"
PROV_ROOT="${RUN_ROOT}/provenance/staging"
COMMANDS="${PROV_ROOT}/commands.log"

mkdir -p "${INPUT_ROOT}/chipseq" "${INPUT_ROOT}/cutrun" "${PROV_ROOT}"

storage_source="$(findmnt -n -o SOURCE -T "${RUN_ROOT}")"
storage_rotational="$(lsblk -ndo ROTA "${storage_source}" 2>/dev/null | awk 'NR == 1 {print $1}' || true)"
[[ "${storage_rotational}" == 0 ]] || guard 1 "${RUN_ROOT} is not on a confirmed SSD"

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

stage_verified() {
  local source="$1" destination="$2" expected="$3"
  if [[ -s "${destination}" && "$(sha256sum "${destination}" | awk '{print $1}')" == "${expected}" ]]; then
    echo "[reuse] ${destination}"
    return
  fi
  [[ ! -e "${destination}" ]] || { echo "ERROR: wrong existing file: ${destination}" >&2; exit 2; }
  [[ -s "${source}" ]] || { echo "ERROR: missing source: ${source}" >&2; exit 2; }
  [[ "$(sha256sum "${source}" | awk '{print $1}')" == "${expected}" ]] || {
    echo "ERROR: source SHA-256 mismatch: ${source}" >&2; exit 2;
  }
  local partial="${destination}.partial"
  log_command cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  [[ "$(sha256sum "${partial}" | awk '{print $1}')" == "${expected}" ]] || {
    echo "ERROR: staged SHA-256 mismatch: ${partial}" >&2; exit 2;
  }
  mv "${partial}" "${destination}"
}

stage_verified \
  ${MACS3_REF_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF686KKV.bam \
  "${INPUT_ROOT}/chipseq/ENCFF686KKV.bam" \
  37ee55258385f908e139e934dcb542360e91ce0dc4dca4867cc35f5507fbe4ea
stage_verified \
  ${MACS3_REF_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF181CXT.bam \
  "${INPUT_ROOT}/chipseq/ENCFF181CXT.bam" \
  400347e0ebf66079f47bb06593129d0237bf1ceac7579dd9442ace0d2a57ad1f
stage_verified \
  ${MACS3_REF_ROOT}/cutrun/GSE172236/SRR14255054.bam \
  "${INPUT_ROOT}/cutrun/SRR14255054.bam" \
  28f4f45a31b4a90da5a68c44b093b166f0c8bde009b04e00d549636d5fcae95a
stage_verified \
  ${MACS3_REF_ROOT}/cutrun/GSE172236/SRR14255053.bam \
  "${INPUT_ROOT}/cutrun/SRR14255053.bam" \
  33f74032ce65e141bb137582201162f4dddee3950e1983e8963097dc9f57149d

find "${INPUT_ROOT}" -type f -print0 | sort -z | xargs -0 sha256sum >"${PROV_ROOT}/sha256sums.txt"
{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_root\t%s\n' "${RUN_ROOT}"
  printf 'storage_source\t%s\n' "${storage_source}"
  printf 'storage_rotational\t%s\n' "${storage_rotational}"
} >"${PROV_ROOT}/metadata.tsv"

echo "[complete] SSD BAMPE inputs: ${INPUT_ROOT}"
