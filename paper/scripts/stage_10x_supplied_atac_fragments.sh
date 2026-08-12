#!/usr/bin/env bash
# Stage and validate the official 10x Cell Ranger ARC PBMC 3k ATAC fragments.

set -euo pipefail

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


RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/fastq_realign}"
DEST_DIR="${RUN_ROOT}/inputs/atac/10x_supplied"
DEST="${DEST_DIR}/pbmc_unsorted_3k_atac_fragments.tsv.gz"
SOURCE_URL="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_atac_fragments.tsv.gz"
SOURCE_FILE="${SOURCE_FILE:-}"
EXPECTED_SIZE=537308082
PROVENANCE_DIR="${RUN_ROOT}/provenance/10x_supplied_atac"
PARTIAL="${DEST}.partial"

mkdir -p "${DEST_DIR}" "${PROVENANCE_DIR}"

if [[ ! -s "${DEST}" ]]; then
  if [[ -n "${SOURCE_FILE}" ]]; then
    [[ -s "${SOURCE_FILE}" ]] || { echo "ERROR: missing ${SOURCE_FILE}" >&2; exit 2; }
    printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${PROVENANCE_DIR}/commands.log"
    printf '\t%q' cp --reflink=auto --preserve=timestamps \
      "${SOURCE_FILE}" "${PARTIAL}" >>"${PROVENANCE_DIR}/commands.log"
    printf '\n' >>"${PROVENANCE_DIR}/commands.log"
    cp --reflink=auto --preserve=timestamps "${SOURCE_FILE}" "${PARTIAL}"
  else
    printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${PROVENANCE_DIR}/commands.log"
    printf '\t%q' curl -fL --retry 5 --retry-delay 2 --continue-at - \
      -o "${PARTIAL}" "${SOURCE_URL}" >>"${PROVENANCE_DIR}/commands.log"
    printf '\n' >>"${PROVENANCE_DIR}/commands.log"
    curl -fL --retry 5 --retry-delay 2 --continue-at - \
      -o "${PARTIAL}" "${SOURCE_URL}"
  fi
  mv "${PARTIAL}" "${DEST}"
fi

actual_size="$(stat -c %s "${DEST}")"
if [[ "${actual_size}" != "${EXPECTED_SIZE}" ]]; then
  echo "ERROR: expected ${EXPECTED_SIZE} bytes, found ${actual_size}: ${DEST}" >&2
  exit 1
fi
gzip -t "${DEST}"

gzip -cd "${DEST}" | awk -F '\t' '
  /^#/ {next}
  NF != 5 {printf "ERROR: line %d has %d columns\n", NR, NF > "/dev/stderr"; exit 1}
  $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/ {
    printf "ERROR: invalid coordinate/count at line %d\n", NR > "/dev/stderr"; exit 1
  }
  {rows++; if ($5 > max_count) max_count=$5}
  END {
    if (!rows) exit 1
    printf "rows\t%d\ncolumns\t5\nmax_count\t%d\n", rows, max_count
  }
' >"${PROVENANCE_DIR}/schema.tsv"

sha256sum "${DEST}" >"${PROVENANCE_DIR}/sha256sums.txt"
{
  printf 'field\tvalue\n'
  printf 'staged_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'source_url\t%s\n' "${SOURCE_URL}"
  printf 'source_file\t%s\n' "${SOURCE_FILE:-not used}"
  if [[ -n "${SOURCE_FILE}" ]]; then
    printf 'source_file_sha256\t%s\n' "$(sha256sum "${SOURCE_FILE}" | awk '{print $1}')"
  fi
  printf 'destination\t%s\n' "${DEST}"
  printf 'size_bytes\t%s\n' "${actual_size}"
  printf 'storage_source\t%s\n' "$(findmnt -n -o SOURCE -T "${DEST}")"
} >"${PROVENANCE_DIR}/metadata.tsv"

echo "[complete] 10x-supplied ATAC fragments: ${DEST}"
echo "[complete] provenance: ${PROVENANCE_DIR}"
