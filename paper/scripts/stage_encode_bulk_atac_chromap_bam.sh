#!/usr/bin/env bash
# Create a coordinate-sorted Chromap BAM for ENCSR095QNB from its source FASTQs.

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


RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/fastq_realign}"
ENCODE_ROOT="${RUN_ROOT}/encode_ENCSR095QNB"
INPUT_ROOT="${ENCODE_ROOT}/inputs"
ALIGN_ROOT="${ENCODE_ROOT}/alignments"
PROV_ROOT="${ENCODE_ROOT}/provenance/chromap_bam"
TMP_ROOT="${ENCODE_ROOT}/tmp/chromap_bam"
REF_ROOT="${RUN_ROOT}/reference"
CHROMAP="${CHROMAP:-${CHROMAP_SUITE}/chromap}"
THREADS="${THREADS:-24}"

R1="${INPUT_ROOT}/ENCFF196VAZ.fastq.gz"
R2="${INPUT_ROOT}/ENCFF195FRB.fastq.gz"
BAM="${ALIGN_ROOT}/ENCSR095QNB.chromap.sorted.bam"
FRAG="${ALIGN_ROOT}/ENCSR095QNB.chromap.fragments.tsv.gz"
SUMMARY="${ALIGN_ROOT}/ENCSR095QNB.chromap.summary.tsv"
COMMANDS="${PROV_ROOT}/commands.log"

mkdir -p "${ALIGN_ROOT}" "${PROV_ROOT}" "${TMP_ROOT}"

for tool in "${CHROMAP}" samtools md5sum sha256sum gzip findmnt lsblk /usr/bin/time; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: missing ${tool}" >&2; exit 2; }
done
[[ "$(md5sum "${R1}" | awk '{print $1}')" == 3909b5296a4f1aa5090274961463c3c0 ]]
[[ "$(md5sum "${R2}" | awk '{print $1}')" == e58702dfb804a94cfb29acfccf825f2f ]]
[[ -s "${REF_ROOT}/genome.fa" && -s "${REF_ROOT}/genome.index" ]]
source_device="$(findmnt -n -o SOURCE -T "${ENCODE_ROOT}")"
[[ "$(lsblk -ndo ROTA "${source_device}" | awk 'NR == 1 {print $1}')" == 0 ]] \
  || guard 1 "${ENCODE_ROOT} is not on a confirmed SSD"
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs |[/]chromap / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another benchmark or Chromap process is active" >&2
  exit 2
fi

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

if [[ ! -s "${BAM}" ]]; then
  [[ ! -e "${BAM}" && ! -e "${FRAG}" ]] || {
    echo "ERROR: refusing to overwrite partial Chromap outputs" >&2; exit 2;
  }
  args=("${CHROMAP}" --preset atac -t "${THREADS}"
        -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa"
        -1 "${R1}" -2 "${R2}" --BAM --sort-bam --write-index
        --hts-threads 4 --sort-bam-ram 20G --atac-fragments "${FRAG}"
        --summary "${SUMMARY}" --temp-dir "${TMP_ROOT}" -o "${BAM}")
  log_command /usr/bin/time -v -o "${PROV_ROOT}/time.txt" "${args[@]}"
  /usr/bin/time -v -o "${PROV_ROOT}/time.txt" "${args[@]}" \
    >"${PROV_ROOT}/stdout.log" 2>"${PROV_ROOT}/stderr.log"
fi

samtools quickcheck -v "${BAM}"
[[ "$(samtools view -H "${BAM}" | awk -F '\t' '$1=="@HD" {for(i=2;i<=NF;i++) if($i=="SO:coordinate") ok=1} END{print ok+0}')" == 1 ]]
index_path=""
for candidate in "${BAM}.bai" "${BAM%.bam}.bai"; do
  if [[ -s "${candidate}" ]]; then index_path="${candidate}"; break; fi
done
[[ -n "${index_path}" ]] || { echo "ERROR: Chromap BAM index is missing" >&2; exit 2; }
gzip -t "${FRAG}"

fragment_rows="$(gzip -cd "${FRAG}" | awk 'END{print NR}')"
bam_records="$(samtools view -c "${BAM}")"
[[ "${bam_records}" == "$((fragment_rows * 2))" ]] || {
  echo "ERROR: expected BAM records = 2 * fragment rows" >&2; exit 2;
}

{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'experiment\tENCSR095QNB\n'
  printf 'bam\t%s\n' "${BAM}"
  printf 'bam_index\t%s\n' "${index_path}"
  printf 'fragment_rows\t%s\n' "${fragment_rows}"
  printf 'bam_records\t%s\n' "${bam_records}"
  printf 'storage_source\t%s\n' "${source_device}"
  printf 'chromap_version\t%s\n' "$("${CHROMAP}" --version 2>&1 | head -n 1)"
  printf 'threads\t%s\n' "${THREADS}"
} >"${PROV_ROOT}/metadata.tsv"
sha256sum "${R1}" "${R2}" "${BAM}" "${index_path}" "${FRAG}" \
  "${SUMMARY}" >"${PROV_ROOT}/sha256sums.txt"

echo "[complete] Chromap coordinate-sorted BAM: ${BAM}"
