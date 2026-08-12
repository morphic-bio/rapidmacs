#!/usr/bin/env bash
# Stage public FASTQs/reference/index on the SSD and realign all three assays
# with Chromap. Bulk outputs remain in Chromap's native seven-column schema;
# five-column compatibility projections are also generated for stock MACS3.

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


REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/fastq_realign}"
CHROMAP="${CHROMAP:-${CHROMAP_SUITE}/chromap}"
THREADS="${THREADS:-24}"
REQUIRE_SSD="${REQUIRE_SSD:-1}"
MIN_FREE_BYTES="${MIN_FREE_BYTES:-8000000000}"

REF_SOURCE="${REF_SOURCE:-${CELLRANGER_REF}/fasta/genome.fa}"
INDEX_SOURCE="${INDEX_SOURCE:-${DATA_DIR}/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
ATAC_SOURCE_DIR="${ATAC_SOURCE_DIR:-${DATA_DIR}/atac-seq/benchmarks/pbmc_unsorted_3k_100k/extracted/pbmc_unsorted_3k/atac}"
ATAC_WHITELIST_SOURCE="${ATAC_WHITELIST_SOURCE:-${DATA_DIR}/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
CUTRUN_SOURCE_DIR="${CUTRUN_SOURCE_DIR:-${MACS3_REF_ROOT}/cutrun/GSE172236}"

INPUT_ROOT="${RUN_ROOT}/inputs"
REF_ROOT="${RUN_ROOT}/reference"
ALIGN_ROOT="${RUN_ROOT}/alignments"
PROV_ROOT="${RUN_ROOT}/provenance/alignment"
TMP_ROOT="${RUN_ROOT}/tmp/chromap"
COMMANDS="${PROV_ROOT}/commands.log"

mkdir -p "${INPUT_ROOT}"/{atac,cutrun,chipseq} "${REF_ROOT}" \
  "${ALIGN_ROOT}"/{atac,cutrun,chipseq} "${PROV_ROOT}" \
  "${TMP_ROOT}"/{atac,cutrun,chipseq}

storage_source="$(findmnt -n -o SOURCE -T "${RUN_ROOT}")"
storage_fstype="$(findmnt -n -o FSTYPE -T "${RUN_ROOT}")"
storage_rotational="$(lsblk -ndo ROTA "${storage_source}" 2>/dev/null | awk 'NR == 1 {print $1}' || true)"
if [[ "${REQUIRE_SSD}" == 1 && "${storage_rotational}" != 0 ]]; then
  echo "ERROR: ${RUN_ROOT} resolves to ${storage_source}, not a confirmed SSD" >&2
  exit 2
fi
available_bytes="$(df --output=avail -B1 "${RUN_ROOT}" | awk 'NR==2 {print $1}')"
if (( available_bytes < MIN_FREE_BYTES )); then
  echo "ERROR: at least ${MIN_FREE_BYTES} bytes free are required; found ${available_bytes}" >&2
  exit 2
fi

for tool in "${CHROMAP}" curl md5sum sha256sum gzip bgzip awk findmnt lsblk /usr/bin/time; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: required tool not found: ${tool}" >&2
    exit 2
  fi
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

stage_local() {
  local source="$1" destination="$2" expected="$3"
  if md5_matches "${destination}" "${expected}"; then
    echo "[reuse] ${destination}"
    return
  fi
  if [[ -e "${destination}" ]]; then
    echo "ERROR: existing staged file has the wrong MD5: ${destination}" >&2
    exit 2
  fi
  if ! md5_matches "${source}" "${expected}"; then
    echo "ERROR: source file is missing or has the wrong MD5: ${source}" >&2
    exit 2
  fi
  local partial="${destination}.partial"
  log_command cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  if ! md5_matches "${partial}" "${expected}"; then
    echo "ERROR: staged copy failed MD5 validation: ${partial}" >&2
    exit 2
  fi
  mv "${partial}" "${destination}"
}

download_verified() {
  local url="$1" destination="$2" expected="$3"
  if md5_matches "${destination}" "${expected}"; then
    echo "[reuse] ${destination}"
    return
  fi
  if [[ -e "${destination}" ]]; then
    echo "ERROR: existing download has the wrong MD5: ${destination}" >&2
    exit 2
  fi
  local partial="${destination}.partial"
  log_command curl -fL --retry 5 --retry-delay 5 --continue-at - --output "${partial}" "${url}"
  curl -fL --retry 5 --retry-delay 5 --continue-at - --output "${partial}" "${url}"
  if ! md5_matches "${partial}" "${expected}"; then
    echo "ERROR: download failed MD5 validation: ${url}" >&2
    exit 2
  fi
  mv "${partial}" "${destination}"
}

stage_reference() {
  local source="$1" destination="$2"
  if [[ -s "${destination}" && "$(stat -c%s "${source}")" == "$(stat -c%s "${destination}")" ]]; then
    echo "[reuse] ${destination}"
    return
  fi
  if [[ -e "${destination}" ]]; then
    echo "ERROR: existing reference artifact has an unexpected size: ${destination}" >&2
    exit 2
  fi
  local partial="${destination}.partial"
  log_command cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  cp --reflink=auto --preserve=timestamps "${source}" "${partial}"
  cmp "${source}" "${partial}"
  mv "${partial}" "${destination}"
}

echo "[stage] reference and Chromap index"
stage_reference "${REF_SOURCE}" "${REF_ROOT}/genome.fa"
stage_reference "${INDEX_SOURCE}" "${REF_ROOT}/genome.index"
stage_local "${ATAC_WHITELIST_SOURCE}" "${REF_ROOT}/737K-arc-v1_atac.txt" \
  3f9fc0f6ef9d72540ab010d0b5348aa1

declare -a atac_r1 atac_r2 atac_bc
atac_md5=(
  85cd1dc430a286f8135764ea121aeb68 6352607f8dc573f133e58464ae63abe9
  d338033d82874135bd73abff3a88758d 574800a3b5b353a89a3d102734cb7870
)
atac_bc_md5=(
  17fad147333eb859eadb719ac62203d8 09840390e35429c2223922c45a1c95e8
  544d090ddc3f3270f1f8f55df227dc6d 1240e491374787a14c00692026d91d44
)
atac_r2_md5=(
  c0d12f0e47707d621aef914724a6ff5d 5be9de74085734ed9eda70e744ca869b
  4bb45a4db19eddd66fd89e1d380389ea 967c7061a28c454bc44e92c0937f3b06
)
echo "[stage] ATAC FASTQs"
for lane in 1 2 3 4; do
  lane3="$(printf '%03d' "${lane}")"
  r1="pbmc_unsorted_3k_S3_L${lane3}_R1_001.fastq.gz"
  bc="pbmc_unsorted_3k_S3_L${lane3}_R2_001.fastq.gz"
  r2="pbmc_unsorted_3k_S3_L${lane3}_R3_001.fastq.gz"
  stage_local "${ATAC_SOURCE_DIR}/${r1}" "${INPUT_ROOT}/atac/${r1}" "${atac_md5[$((lane - 1))]}"
  stage_local "${ATAC_SOURCE_DIR}/${bc}" "${INPUT_ROOT}/atac/${bc}" "${atac_bc_md5[$((lane - 1))]}"
  stage_local "${ATAC_SOURCE_DIR}/${r2}" "${INPUT_ROOT}/atac/${r2}" "${atac_r2_md5[$((lane - 1))]}"
  atac_r1+=("${INPUT_ROOT}/atac/${r1}")
  atac_bc+=("${INPUT_ROOT}/atac/${bc}")
  atac_r2+=("${INPUT_ROOT}/atac/${r2}")
done

echo "[stage] CUT&RUN treatment FASTQs"
stage_local "${CUTRUN_SOURCE_DIR}/SRR14255054_1.fastq.gz" \
  "${INPUT_ROOT}/cutrun/SRR14255054_1.fastq.gz" c3a10f1dbadba9e21bd2d6006d968bb2
stage_local "${CUTRUN_SOURCE_DIR}/SRR14255054_2.fastq.gz" \
  "${INPUT_ROOT}/cutrun/SRR14255054_2.fastq.gz" 110e63a1298591f5cad84b942bcebb0a

echo "[stage] ChIP-seq treatment FASTQs from ENCODE"
download_verified \
  'https://www.encodeproject.org/files/ENCFF058OSL/@@download/ENCFF058OSL.fastq.gz' \
  "${INPUT_ROOT}/chipseq/ENCFF058OSL.fastq.gz" d215190816a71c629f2e9bffd482292f
download_verified \
  'https://www.encodeproject.org/files/ENCFF712FTQ/@@download/ENCFF712FTQ.fastq.gz' \
  "${INPUT_ROOT}/chipseq/ENCFF712FTQ.fastq.gz" fba51e5637bac52ca93a7d2bc1a8cfd8

join_comma() {
  local IFS=,
  printf '%s' "$*"
}

run_alignment() {
  local assay="$1" output="$2"
  shift 2
  if [[ -s "${output}" ]]; then
    echo "[reuse] ${assay} alignment output: ${output}"
    return
  fi
  if [[ -e "${output}" ]]; then
    echo "ERROR: refusing to overwrite ${output}; use a new RUN_ROOT" >&2
    exit 2
  fi
  local log_dir="${PROV_ROOT}/${assay}"
  mkdir -p "${log_dir}"
  log_command "$@"
  /usr/bin/time -v -o "${log_dir}/time.txt" "$@" \
    >"${log_dir}/stdout.log" 2>"${log_dir}/stderr.log"
}

compress_fragment_output() {
  local assay="$1" source="$2" destination="$3"
  if [[ -s "${destination}" ]] && gzip -t "${destination}" 2>/dev/null; then
    echo "[reuse] compressed ${assay} fragments: ${destination}"
    return
  fi
  [[ -s "${source}" ]] || { echo "ERROR: missing plain fragment output: ${source}" >&2; exit 2; }
  local partial="${destination}.partial"
  log_command bgzip -@ "${THREADS}" -c "${source}" ">" "${partial}"
  /usr/bin/time -v -o "${PROV_ROOT}/${assay}/bgzip.time.txt" \
    bgzip -@ "${THREADS}" -c "${source}" >"${partial}"
  gzip -t "${partial}"
  mv "${partial}" "${destination}"
  log_command rm -f "${source}"
  rm -f "${source}"
}

atac_plain="${ALIGN_ROOT}/atac/fragments.tsv"
atac_out="${atac_plain}.gz"
if [[ ! -s "${atac_out}" ]] || ! gzip -t "${atac_out}" 2>/dev/null; then
  if [[ ! -s "${atac_plain}" ]]; then
    run_alignment atac "${atac_plain}" \
      "${CHROMAP}" --preset atac -t "${THREADS}" \
      -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa" \
      -1 "$(join_comma "${atac_r1[@]}")" \
      -2 "$(join_comma "${atac_r2[@]}")" \
      -b "$(join_comma "${atac_bc[@]}")" \
      --barcode-whitelist "${REF_ROOT}/737K-arc-v1_atac.txt" \
      --summary "${ALIGN_ROOT}/atac/chromap.summary.tsv" \
      --temp-dir "${TMP_ROOT}/atac" -o "${atac_plain}"
  fi
fi
compress_fragment_output atac "${atac_plain}" "${atac_out}"

cutrun_plain="${ALIGN_ROOT}/cutrun/fragments.bulk.bed"
cutrun_out="${cutrun_plain}.gz"
if [[ ! -s "${cutrun_out}" ]] || ! gzip -t "${cutrun_out}" 2>/dev/null; then
  if [[ ! -s "${cutrun_plain}" ]]; then
    run_alignment cutrun "${cutrun_plain}" \
      "${CHROMAP}" --preset chip --trim-adapters -t "${THREADS}" \
      -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa" \
      -1 "${INPUT_ROOT}/cutrun/SRR14255054_1.fastq.gz" \
      -2 "${INPUT_ROOT}/cutrun/SRR14255054_2.fastq.gz" \
      --summary "${ALIGN_ROOT}/cutrun/chromap.summary.tsv" \
      --temp-dir "${TMP_ROOT}/cutrun" -o "${cutrun_plain}"
  fi
fi
compress_fragment_output cutrun "${cutrun_plain}" "${cutrun_out}"

chipseq_plain="${ALIGN_ROOT}/chipseq/fragments.bulk.bed"
chipseq_out="${chipseq_plain}.gz"
if [[ ! -s "${chipseq_out}" ]] || ! gzip -t "${chipseq_out}" 2>/dev/null; then
  if [[ ! -s "${chipseq_plain}" ]]; then
    run_alignment chipseq "${chipseq_plain}" \
      "${CHROMAP}" --preset chip -t "${THREADS}" \
      -x "${REF_ROOT}/genome.index" -r "${REF_ROOT}/genome.fa" \
      -1 "${INPUT_ROOT}/chipseq/ENCFF058OSL.fastq.gz" \
      -2 "${INPUT_ROOT}/chipseq/ENCFF712FTQ.fastq.gz" \
      --summary "${ALIGN_ROOT}/chipseq/chromap.summary.tsv" \
      --temp-dir "${TMP_ROOT}/chipseq" -o "${chipseq_plain}"
  fi
fi
compress_fragment_output chipseq "${chipseq_plain}" "${chipseq_out}"

project_bulk_for_macs3() {
  local source="$1" destination="$2"
  if [[ -s "${destination}" ]]; then
    echo "[reuse] MACS3 projection: ${destination}"
    return
  fi
  local partial="${destination}.partial"
  local pipeline
  printf -v pipeline \
    "gzip -cd %q | awk -F '\\t' -v OFS='\\t' 'NF == 7 {print \$1,\$2,\$3,\"bulk\",\$7}' | bgzip -@ %q -c > %q" \
    "${source}" "${THREADS}" "${partial}"
  log_command bash -o pipefail -c "${pipeline}"
  gzip -cd "${source}" |
    awk -F '\t' -v OFS='\t' 'NF == 7 {print $1,$2,$3,"bulk",$7}' |
    bgzip -@ "${THREADS}" -c >"${partial}"
  mv "${partial}" "${destination}"
}

project_bulk_for_macs3 "${cutrun_out}" \
  "${ALIGN_ROOT}/cutrun/fragments.macs3.5col.tsv.gz"
project_bulk_for_macs3 "${chipseq_out}" \
  "${ALIGN_ROOT}/chipseq/fragments.macs3.5col.tsv.gz"

validate_columns() {
  local path="$1" expected="$2"
  local observed
  observed="$(gzip -cd "${path}" | awk -F '\t' '!/^#/ && NF {seen[NF]++} END {if (length(seen) != 1) {print "mixed"; exit} for (n in seen) print n}')"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "ERROR: ${path}: expected ${expected} columns, observed ${observed}" >&2
    exit 2
  fi
}
validate_columns "${atac_out}" 5
validate_columns "${cutrun_out}" 7
validate_columns "${chipseq_out}" 7
validate_columns "${ALIGN_ROOT}/cutrun/fragments.macs3.5col.tsv.gz" 5
validate_columns "${ALIGN_ROOT}/chipseq/fragments.macs3.5col.tsv.gz" 5

{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_root\t%s\n' "${RUN_ROOT}"
  printf 'storage_source\t%s\n' "${storage_source}"
  printf 'storage_fstype\t%s\n' "${storage_fstype}"
  printf 'storage_rotational\t%s\n' "${storage_rotational}"
  printf 'chromap_version\t%s\n' "$("${CHROMAP}" --version 2>&1 | head -n 1)"
  printf 'chromap_sha256\t%s\n' "$(sha256sum "${CHROMAP}" | awk '{print $1}')"
  printf 'threads\t%s\n' "${THREADS}"
} >"${PROV_ROOT}/metadata.tsv"

find "${INPUT_ROOT}" "${REF_ROOT}" "${ALIGN_ROOT}" -type f -print0 |
  sort -z | xargs -0 sha256sum >"${PROV_ROOT}/sha256sums.txt"
cp "${REPO_ROOT}/paper/scripts/fastq_sources.tsv" "${PROV_ROOT}/fastq_sources.tsv"

echo "[complete] FASTQ realignment and fragments: ${ALIGN_ROOT}"
echo "[complete] provenance: ${PROV_ROOT}"
