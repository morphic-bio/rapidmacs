#!/usr/bin/env bash
# MACS3 treatment-only baseline for direct or Chromap-realigned ENCSR095QNB.

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


REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/fastq_realign}"
ENCODE_ROOT="${RUN_ROOT}/encode_ENCSR095QNB"
MODE="${MODE:-}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
MACS3="${MACS3:-${MACS3_REF_ROOT}/macs3-3.0.3-no-control-dump-venv/bin/macs3}"
OUT_ROOT="${OUT_ROOT:-${ENCODE_ROOT}/results/$(date -u +%Y%m%dT%H%M%SZ)_${MODE}_macs3}"
RAW_TSV="${OUT_ROOT}/raw.tsv"
COMMANDS="${OUT_ROOT}/provenance/commands.log"

case "${MODE}" in
  bam)
    INPUT="${ENCODE_ROOT}/inputs/ENCFF646NWY.bam"
    FORMAT=BAMPE
    NAME=ENCSR095QNB_ENCFF646NWY_direct
    ;;
  chromap_bam)
    INPUT="${ENCODE_ROOT}/alignments/ENCSR095QNB.chromap.sorted.bam"
    FORMAT=BAMPE
    NAME=ENCSR095QNB_chromap_bam
    ;;
  frag)
    INPUT="${ENCODE_ROOT}/alignments/fragments.macs3.5col.tsv.gz"
    NATIVE_FRAG="${ENCODE_ROOT}/alignments/fragments.bulk.bed.gz"
    FORMAT=FRAG
    NAME=ENCSR095QNB_chromap_realign
    ;;
  *) echo "ERROR: set MODE=bam, MODE=chromap_bam, or MODE=frag" >&2; exit 2 ;;
esac

mkdir -p "${OUT_ROOT}/runs" "${OUT_ROOT}/provenance"
[[ -s "${INPUT}" ]] || { echo "ERROR: missing ${INPUT}" >&2; exit 2; }
source_device="$(findmnt -n -o SOURCE -T "${INPUT}")"
[[ "$(lsblk -ndo ROTA "${source_device}" | awk 'NR == 1 {print $1}')" == 0 ]] \
  || guard 1 "input is not on a confirmed SSD"
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another peak-calling run is active" >&2
  exit 2
fi
[[ "$("${MACS3}" --version 2>&1 | awk 'NR == 1 {print $NF}')" == 3.0.3 ]]

printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

bed3_md5() {
  awk -F '\t' -v OFS='\t' '!/^track/ && NF >= 3 {print $1,$2,$3}' "$1" | md5sum | awk '{print $1}'
}

fragment_mean() {
  gzip -cd "$1" | awk -F '\t' '{sum += $3-$2; n++} END {if (!n) exit 2; printf "%d\n", sum/n}'
}

run_one() {
  local repeat="$1"
  local dir="${OUT_ROOT}/runs/${MODE}.macs3.r${repeat}"
  mkdir -p "${dir}/out"
  local args=("${MACS3}" callpeak -t "${INPUT}" -f "${FORMAT}" -g hs -n "${NAME}" --outdir "${dir}/out" -q 0.05)
  if [[ "${MODE}" == frag ]]; then
    geometry="$(fragment_mean "${NATIVE_FRAG}")"
    args+=(--max-count 1 --min-length "${geometry}" --max-gap "${geometry}")
  fi
  log_command /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}"
  if ! /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}" \
      >"${dir}/stdout.log" 2>"${dir}/stderr.log"; then
    echo "ERROR: MACS3 failed for ${MODE} ${repeat}" >&2
    return 1
  fi
  read -r wall rss <"${dir}/time.tsv"
  np="${dir}/out/${NAME}_peaks.narrowPeak"
  sm="${dir}/out/${NAME}_summits.bed"
  printf '%s\tMACS3\t1\t%s\twarm\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${MODE}" "${repeat}" "${wall}" "${rss}" "$(wc -l <"${np}")" \
    "$(md5sum "${np}" | awk '{print $1}')" "$(bed3_md5 "${np}")" \
    "$(md5sum "${sm}" | awk '{print $1}')" "$(bed3_md5 "${sm}")" >>"${RAW_TSV}"
}

for ((warm=1; warm<=WARMUPS; ++warm)); do
  echo "[warmup] ${MODE} MACS3 ${warm}/${WARMUPS}"
  run_one "warmup${warm}"
  sed -i '$d' "${RAW_TSV}"
done
for ((repeat=1; repeat<=REPEATS; ++repeat)); do
  echo "[benchmark] ${MODE} MACS3 ${repeat}/${REPEATS}"
  run_one "${repeat}"
done

if ! awk -F '\t' 'NR == 2 {np=$9; sm=$11} NR > 2 && ($9 != np || $11 != sm) {bad=1} END {exit bad}' "${RAW_TSV}"; then
  echo "ERROR: MACS3 output changed across repetitions" >&2
  exit 1
fi
python3 "${REPO_ROOT}/paper/scripts/summarize_paper_benchmarks.py" "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'experiment\tENCSR095QNB\n'
  printf 'mode\t%s\n' "${MODE}"
  printf 'input\t%s\n' "${INPUT}"
  printf 'input_sha256\t%s\n' "$(sha256sum "${INPUT}" | awk '{print $1}')"
  printf 'storage_source\t%s\n' "${source_device}"
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'macs3_version\t%s\n' "$("${MACS3}" --version 2>&1 | head -n 1)"
} >"${OUT_ROOT}/provenance/metadata.tsv"
cp "${REPO_ROOT}/paper/scripts/"*.sh "${OUT_ROOT}/provenance/"
echo "[complete] MACS3 baseline: ${OUT_ROOT}/summary.tsv"
