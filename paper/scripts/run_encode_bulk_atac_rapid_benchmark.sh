#!/usr/bin/env bash
# RapidMACS-only parallel benchmark for direct and FASTQ-realigned ENCSR095QNB.

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
THREADS="${THREADS:-4 24}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
OUT_ROOT="${OUT_ROOT:-${ENCODE_ROOT}/results/$(date -u +%Y%m%dT%H%M%SZ)_${MODE}_rapid_only}"
RAW_TSV="${OUT_ROOT}/raw.tsv"
STAGE_TSV="${OUT_ROOT}/stages_raw.tsv"
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
    INPUT="${ENCODE_ROOT}/alignments/fragments.bulk.bed.gz"
    FORMAT=FRAG
    NAME=ENCSR095QNB_chromap_realign
    ;;
  *) echo "ERROR: set MODE=bam, MODE=chromap_bam, or MODE=frag" >&2; exit 2 ;;
esac

mkdir -p "${OUT_ROOT}/runs" "${OUT_ROOT}/provenance"
[[ -s "${INPUT}" ]] || { echo "ERROR: missing ${INPUT}" >&2; exit 2; }
[[ "${INPUT}" == "${RUN_ROOT}"/* ]] || { echo "ERROR: input is not below ${RUN_ROOT}" >&2; exit 2; }
[[ "$(lsblk -ndo ROTA "$(findmnt -n -o SOURCE -T "${INPUT}")" | awk 'NR == 1 {print $1}')" == 0 ]] \
  || guard 1 "input is not on a confirmed SSD"
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another peak-calling experiment is active" >&2
  exit 2
fi
unset RAPIDMACS_LAZY_CONTROL_WINDOWS

printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"
printf 'assay\tthreads\trepeat\tcache_state\tload_s\tlibrary_pipeline_s\tbuild_workspace_s\ttreat_s\tlambda_s\tppois_regions_s\toutput_s\tprofiled_total_s\n' >"${STAGE_TSV}"

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

record_stage() {
  local threads="$1" repeat="$2" file="$3"
  awk -F '\t' -v a="${MODE}" -v t="${threads}" -v r="${repeat}" '
    $1 == "load_fragments" {load=$2}
    $1 == "build_workspace" {build=$2}
    $1 == "finalize_treat_tracks" {treat=$2}
    $1 == "finalize_lambda_tracks" {lambda=$2}
    $1 == "finalize_ppois_tracks" {ppois=$2}
    $1 == "score_build_ppois" {ppois_build=$2}
    $1 == "score_build_pqtable" {ppois_qtable=$2}
    $1 == "score_collect_regions" {regions=$2}
    $1 == "narrowpeak_summits" {output=$2}
    $1 == "total" {total=$2}
    END {
      if (load == "" || total == "") exit 2
      if (ppois == "") ppois=ppois_build+ppois_qtable+regions
      printf "%s\t%s\t%s\twarm\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
             a,t,r,load,total-load,build,treat,lambda,ppois,output,total
    }
  ' "${file}" >>"${STAGE_TSV}"
}

run_one() {
  local threads="$1" repeat="$2"
  local dir="${OUT_ROOT}/runs/${MODE}.rapidmacs.t${threads}.r${repeat}"
  local prefix="${dir}/${NAME}" stage="${dir}/stages.tsv"
  mkdir -p "${dir}"
  local args
  if [[ "${MODE}" == bam || "${MODE}" == chromap_bam ]]; then
    args=("${RAPIDMACS}" --preset cutrun -f BAMPE -t "${INPUT}" -g hs -n "${NAME}"
          -q 0.05 --out-prefix "${prefix}" --peak-caller-threads "${threads}"
          --profile-stages "${stage}")
  else
    geometry="$(fragment_mean "${INPUT}")"
    args=("${RAPIDMACS}" -f FRAG -t "${INPUT}" -g hs -n "${NAME}"
          --macs3-frag-qvalue 0.05 --frag-pileup-max-count 1
          --bdgpeakcall-min-len "${geometry}" --bdgpeakcall-max-gap "${geometry}"
          --macs3-frag-narrowpeak "${prefix}.narrowPeak"
          --macs3-frag-summits "${prefix}.summits.bed"
          --peak-caller-threads "${threads}" --profile-stages "${stage}")
  fi
  log_command /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" env -u RAPIDMACS_LAZY_CONTROL_WINDOWS "${args[@]}"
  if ! /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" \
      env -u RAPIDMACS_LAZY_CONTROL_WINDOWS "${args[@]}" \
      >"${dir}/stdout.log" 2>"${dir}/stderr.log"; then
    echo "ERROR: RapidMACS failed for ${MODE} T=${threads} ${repeat}" >&2
    return 1
  fi
  read -r wall rss <"${dir}/time.tsv"
  np="${prefix}.narrowPeak"; sm="${prefix}.summits.bed"
  printf '%s\tRapidMACS\t%s\t%s\twarm\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${MODE}" "${threads}" "${repeat}" "${wall}" "${rss}" "$(wc -l <"${np}")" \
    "$(md5sum "${np}" | awk '{print $1}')" "$(bed3_md5 "${np}")" \
    "$(md5sum "${sm}" | awk '{print $1}')" "$(bed3_md5 "${sm}")" >>"${RAW_TSV}"
  record_stage "${threads}" "${repeat}" "${stage}"
}

for threads in ${THREADS}; do
  for ((warm=1; warm<=WARMUPS; ++warm)); do
    echo "[warmup] ${MODE} RapidMACS T=${threads} ${warm}/${WARMUPS}"
    run_one "${threads}" "warmup${warm}"
    sed -i '$d' "${RAW_TSV}"
    sed -i '$d' "${STAGE_TSV}"
  done
  for ((repeat=1; repeat<=REPEATS; ++repeat)); do
    echo "[benchmark] ${MODE} RapidMACS T=${threads} ${repeat}/${REPEATS}"
    run_one "${threads}" "${repeat}"
  done
done

if ! awk -F '\t' 'NR == 2 {np=$9; sm=$11} NR > 2 && ($9 != np || $11 != sm) {bad=1} END {exit bad}' "${RAW_TSV}"; then
  echo "ERROR: output changed across runs or thread counts" >&2
  exit 1
fi

python3 "${REPO_ROOT}/paper/scripts/summarize_paper_benchmarks.py" "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
python3 "${REPO_ROOT}/paper/scripts/summarize_stage_benchmarks.py" "${STAGE_TSV}" "${OUT_ROOT}/stages_summary.tsv"
git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'experiment\tENCSR095QNB\n'
  printf 'mode\t%s\n' "${MODE}"
  printf 'input\t%s\n' "${INPUT}"
  printf 'input_sha256\t%s\n' "$(sha256sum "${INPUT}" | awk '{print $1}')"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'storage_source\t%s\n' "$(findmnt -n -o SOURCE -T "${INPUT}")"
  printf 'rapidmacs_sha256\t%s\n' "$(sha256sum "${RAPIDMACS}" | awk '{print $1}')"
  printf 'RAPIDMACS_LAZY_CONTROL_WINDOWS\tunset\n'
} >"${OUT_ROOT}/provenance/metadata.tsv"
cp "${REPO_ROOT}/paper/scripts/"*.sh \
  "${REPO_ROOT}/paper/scripts/fastq_sources.tsv" \
  "${REPO_ROOT}/paper/scripts/fragment_sources.tsv" "${OUT_ROOT}/provenance/"

echo "[complete] summary: ${OUT_ROOT}/summary.tsv"
echo "[complete] provenance: ${OUT_ROOT}/provenance"
