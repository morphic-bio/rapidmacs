#!/usr/bin/env bash
# Stage ENCSR095QNB inputs or realign its replicate-2 FASTQs, one mode at a time.

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


MODE="${1:-}"
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/fastq_realign}"
ENCODE_ROOT="${RUN_ROOT}/encode_ENCSR095QNB"
INPUT_ROOT="${ENCODE_ROOT}/inputs"
ALIGN_ROOT="${ENCODE_ROOT}/alignments"
PROV_ROOT="${ENCODE_ROOT}/provenance"
TMP_ROOT="${ENCODE_ROOT}/tmp"
REF_ROOT="${RUN_ROOT}/reference"
CHROMAP="${CHROMAP:-${CHROMAP_SUITE}/chromap}"
THREADS="${THREADS:-24}"
COMMANDS="${PROV_ROOT}/commands.log"

BAM="${INPUT_ROOT}/ENCFF646NWY.bam"
R1="${INPUT_ROOT}/ENCFF196VAZ.fastq.gz"
R2="${INPUT_ROOT}/ENCFF195FRB.fastq.gz"
FRAG_PLAIN="${ALIGN_ROOT}/fragments.bulk.bed"
FRAG="${FRAG_PLAIN}.gz"

mkdir -p "${INPUT_ROOT}" "${ALIGN_ROOT}" "${PROV_ROOT}" "${TMP_ROOT}"

case "${MODE}" in
  bam|fastq|align) ;;
  *) echo "Usage: $0 {bam|fastq|align}" >&2; exit 2 ;;
esac

storage_source="$(findmnt -n -o SOURCE -T "${RUN_ROOT}")"
storage_rotational="$(lsblk -ndo ROTA "${storage_source}" 2>/dev/null | awk 'NR == 1 {print $1}' || true)"
if [[ "${storage_rotational}" != 0 ]]; then
  echo "ERROR: ${RUN_ROOT} is not on a confirmed SSD" >&2
  exit 2
fi

for tool in curl md5sum sha256sum gzip bgzip awk findmnt lsblk /usr/bin/time; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: missing ${tool}" >&2; exit 2; }
done

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

md5_matches() {
  local path="$1" expected="$2"
  [[ -s "${path}" ]] && [[ "$(md5sum "${path}" | awk '{print $1}')" == "${expected}" ]]
}

download_verified() {
  local accession="$1" suffix="$2" destination="$3" expected="$4"
  local url="https://www.encodeproject.org/files/${accession}/@@download/${accession}.${suffix}"
  if md5_matches "${destination}" "${expected}"; then
    echo "[reuse] ${destination}"
    return
  fi
  if [[ -e "${destination}" ]]; then
    echo "ERROR: existing file has wrong MD5: ${destination}" >&2
    exit 2
  fi
  local partial="${destination}.partial"
  log_command curl -fL --retry 5 --retry-delay 5 --continue-at - --output "${partial}" "${url}"
  curl -fL --retry 5 --retry-delay 5 --continue-at - --output "${partial}" "${url}"
  if ! md5_matches "${partial}" "${expected}"; then
    echo "ERROR: MD5 validation failed: ${url}" >&2
    exit 2
  fi
  mv "${partial}" "${destination}"
}

case "${MODE}" in
  bam)
    download_verified ENCFF646NWY bam "${BAM}" 63f2d76a79216d623872f66653ec0d41
    ;;
  fastq)
    download_verified ENCFF196VAZ fastq.gz "${R1}" 3909b5296a4f1aa5090274961463c3c0
    download_verified ENCFF195FRB fastq.gz "${R2}" e58702dfb804a94cfb29acfccf825f2f
    ;;
  align)
    md5_matches "${R1}" 3909b5296a4f1aa5090274961463c3c0 || { echo "ERROR: missing/invalid ${R1}" >&2; exit 2; }
    md5_matches "${R2}" e58702dfb804a94cfb29acfccf825f2f || { echo "ERROR: missing/invalid ${R2}" >&2; exit 2; }
    [[ -s "${REF_ROOT}/genome.fa" && -s "${REF_ROOT}/genome.index" ]]
    if [[ ! -s "${FRAG}" ]]; then
      [[ ! -e "${FRAG_PLAIN}" ]] || { echo "ERROR: refusing to overwrite ${FRAG_PLAIN}" >&2; exit 2; }
      log_command "${CHROMAP}" --preset atac -t "${THREADS}" \
        -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa" \
        -1 "${R1}" -2 "${R2}" --summary "${ALIGN_ROOT}/chromap.summary.tsv" \
        --temp-dir "${TMP_ROOT}" -o "${FRAG_PLAIN}"
      /usr/bin/time -v -o "${PROV_ROOT}/chromap.time.txt" \
        "${CHROMAP}" --preset atac -t "${THREADS}" \
        -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa" \
        -1 "${R1}" -2 "${R2}" --summary "${ALIGN_ROOT}/chromap.summary.tsv" \
        --temp-dir "${TMP_ROOT}" -o "${FRAG_PLAIN}" \
        >"${PROV_ROOT}/chromap.stdout.log" 2>"${PROV_ROOT}/chromap.stderr.log"
      partial="${FRAG}.partial"
      printf -v pipeline 'bgzip -@ %q -c %q > %q' "${THREADS}" "${FRAG_PLAIN}" "${partial}"
      log_command bash -o pipefail -c "${pipeline}"
      /usr/bin/time -v -o "${PROV_ROOT}/bgzip.time.txt" \
        bash -o pipefail -c "${pipeline}"
      gzip -t "${partial}"
      mv "${partial}" "${FRAG}"
      log_command rm -f "${FRAG_PLAIN}"
      rm -f "${FRAG_PLAIN}"
    fi
    gzip -t "${FRAG}"
    gzip -cd "${FRAG}" | awk -F '\t' '
      NF != 7 {printf "ERROR: line %d has %d columns\n", NR, NF > "/dev/stderr"; exit 1}
      $4 != "N" || $5 !~ /^[0-9]+$/ || ($6 != "+" && $6 != "-") || $7 !~ /^[0-9]+$/ {
        printf "ERROR: invalid Chromap bulk row at line %d\n", NR > "/dev/stderr"; exit 1
      }
      {rows++; bp += $3-$2; if ($7 > max_count) max_count=$7}
      END {if (!rows) exit 1; printf "rows\t%d\ncolumns\t7\nmean_bp\t%d\nmax_count\t%d\n", rows, bp/rows, max_count}
    ' >"${PROV_ROOT}/fragment_schema.tsv"
    ;;
esac

{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'mode\t%s\n' "${MODE}"
  printf 'experiment\tENCSR095QNB\n'
  printf 'requested_alignment\tENCFF646NWY\n'
  printf 'source_read1\tENCFF196VAZ\n'
  printf 'source_read2\tENCFF195FRB\n'
  printf 'storage_source\t%s\n' "${storage_source}"
  printf 'storage_rotational\t%s\n' "${storage_rotational}"
} >"${PROV_ROOT}/${MODE}.metadata.tsv"

find "${INPUT_ROOT}" "${ALIGN_ROOT}" -maxdepth 1 -type f -print0 | sort -z | \
  xargs -0 -r sha256sum >"${PROV_ROOT}/sha256sums.txt"
echo "[complete] ENCSR095QNB stage mode: ${MODE}"
