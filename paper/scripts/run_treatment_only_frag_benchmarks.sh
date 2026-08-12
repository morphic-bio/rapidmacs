#!/usr/bin/env bash
# Treatment-only FRAG peak-caller benchmark using SSD-resident Chromap outputs.
# One unrecorded warm-up precedes each set of recorded complete-CLI runs.

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
ALIGN_ROOT="${RUN_ROOT}/alignments"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
MACS3="${MACS3:-${MACS3_REF_ROOT}/macs3-3.0.3-no-control-dump-venv/bin/macs3}"
THREADS="${THREADS:-1 4 24}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
ASSAYS="${ASSAYS:-atac cutrun chipseq}"
TOOLS="${TOOLS:-macs3 rapidmacs}"
OUT_ROOT="${OUT_ROOT:-${RUN_ROOT}/results/$(date -u +%Y%m%dT%H%M%SZ)_treatment_only_frag}"
RAW_TSV="${OUT_ROOT}/raw.tsv"
STAGE_RAW_TSV="${OUT_ROOT}/stages_raw.tsv"
COMMANDS="${OUT_ROOT}/provenance/commands.log"

mkdir -p "${OUT_ROOT}/runs" "${OUT_ROOT}/provenance"

tool_enabled() {
  [[ " ${TOOLS} " == *" $1 "* ]]
}

for tool in awk md5sum sha256sum gzip /usr/bin/time; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: missing ${tool}" >&2; exit 2; }
done
if tool_enabled macs3; then command -v "${MACS3}" >/dev/null 2>&1 || exit 2; fi
if tool_enabled rapidmacs; then command -v "${RAPIDMACS}" >/dev/null 2>&1 || exit 2; fi
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs |[/]chromap / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another alignment or peak-calling experiment is active" >&2
  exit 2
fi

atac_rapid="${ATAC_FRAGMENT_INPUT:-${ALIGN_ROOT}/atac/fragments.tsv.gz}"
atac_macs3="${atac_rapid}"
cutrun_rapid="${ALIGN_ROOT}/cutrun/fragments.bulk.bed.gz"
cutrun_macs3="${ALIGN_ROOT}/cutrun/fragments.macs3.5col.tsv.gz"
chipseq_rapid="${ALIGN_ROOT}/chipseq/fragments.bulk.bed.gz"
chipseq_macs3="${ALIGN_ROOT}/chipseq/fragments.macs3.5col.tsv.gz"

for assay in ${ASSAYS}; do
  case "${assay}" in
    atac) required=("${atac_rapid}") ;;
    cutrun) required=("${cutrun_rapid}" "${cutrun_macs3}") ;;
    chipseq) required=("${chipseq_rapid}" "${chipseq_macs3}") ;;
    *) echo "ERROR: unknown assay: ${assay}" >&2; exit 2 ;;
  esac
  for path in "${required[@]}"; do
    [[ -s "${path}" ]] || { echo "ERROR: missing alignment output: ${path}" >&2; exit 2; }
    [[ "${path}" == "${RUN_ROOT}"/* ]] || { echo "ERROR: input is not under SSD run root: ${path}" >&2; exit 2; }
  done
done

if tool_enabled macs3 && [[ "$("${MACS3}" --version 2>&1 | awk 'NR==1 {print $NF}')" != 3.0.3 ]]; then
  echo "ERROR: MACS3 3.0.3 is required" >&2
  exit 2
fi

printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"
printf 'assay\tthreads\trepeat\tcache_state\tload_s\tlibrary_pipeline_s\tbuild_workspace_s\ttreat_s\tlambda_s\tppois_regions_s\toutput_s\tprofiled_total_s\n' >"${STAGE_RAW_TSV}"

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

bed3_md5() {
  awk -F '\t' -v OFS='\t' '!/^track/ && NF >= 3 {print $1,$2,$3}' "$1" |
    md5sum | awk '{print $1}'
}

fragment_mean() {
  gzip -cd "$1" |
    awk -F '\t' '!/^#/ && NF >= 3 {sum += $3-$2; n++} END {if (!n) exit 2; printf "%d\n", sum/n}'
}

cutrun_geometry=""
chipseq_geometry=""
if [[ " ${ASSAYS} " == *" cutrun "* ]]; then
  cutrun_geometry="$(fragment_mean "${cutrun_rapid}")"
fi
if [[ " ${ASSAYS} " == *" chipseq "* ]]; then
  chipseq_geometry="$(fragment_mean "${chipseq_rapid}")"
fi
{
  printf 'assay\tmin_length\tmax_gap\tbasis\n'
  if [[ " ${ASSAYS} " == *" atac "* ]]; then
    printf 'atac\t200\t30\tpaper FRAG protocol\n'
  fi
  if [[ " ${ASSAYS} " == *" cutrun "* ]]; then
    printf 'cutrun\t%s\t%s\tinteger treatment mean fragment length\n' "${cutrun_geometry}" "${cutrun_geometry}"
  fi
  if [[ " ${ASSAYS} " == *" chipseq "* ]]; then
    printf 'chipseq\t%s\t%s\tinteger treatment mean fragment length\n' "${chipseq_geometry}" "${chipseq_geometry}"
  fi
} >"${OUT_ROOT}/fragment_geometry.tsv"

record_run() {
  local assay="$1" tool="$2" threads="$3" repeat="$4" np="$5" sm="$6" time_file="$7"
  local wall rss peaks np_md5 b3 sm_md5 sm_b3
  read -r wall rss <"${time_file}"
  peaks="$(wc -l <"${np}")"
  np_md5="$(md5sum "${np}" | awk '{print $1}')"
  b3="$(bed3_md5 "${np}")"
  sm_md5="$(md5sum "${sm}" | awk '{print $1}')"
  sm_b3="$(bed3_md5 "${sm}")"
  printf '%s\t%s\t%s\t%s\twarm\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${assay}" "${tool}" "${threads}" "${repeat}" "${wall}" "${rss}" \
    "${peaks}" "${np_md5}" "${b3}" "${sm_md5}" "${sm_b3}" >>"${RAW_TSV}"
}

record_stage_run() {
  local assay="$1" threads="$2" repeat="$3" stage_file="$4"
  awk -F '\t' -v a="${assay}" -v t="${threads}" -v r="${repeat}" '
    $1 == "load_fragments" {load=$2}
    $1 == "build_workspace" {build=$2}
    $1 == "finalize_treat_tracks" {treat=$2}
    $1 == "finalize_lambda_tracks" {lambda=$2}
    $1 == "finalize_ppois_tracks" {ppois=$2}
    $1 == "score_build_ppois" {ppois_build=$2}
    $1 == "score_build_pqtable" {ppois_qtable=$2}
    $1 == "score_collect_regions" {ppois_regions=$2}
    $1 == "narrowpeak_summits" {output=$2}
    $1 == "total" {total=$2}
    END {
      if (load == "" || total == "") exit 2
      if (ppois == "") ppois=ppois_build+ppois_qtable+ppois_regions
      printf "%s\t%s\t%s\twarm\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
             a,t,r,load,total-load,build,treat,lambda,ppois,output,total
    }
  ' "${stage_file}" >>"${STAGE_RAW_TSV}"
}

run_macs3() {
  local assay="$1" repeat="$2" input="$3" geometry="$4" name="$5"
  local dir="${OUT_ROOT}/runs/${assay}.macs3.r${repeat}"
  mkdir -p "${dir}/out"
  local args=("${MACS3}" callpeak -t "${input}" -f FRAG -g hs -n "${name}" --outdir "${dir}/out")
  if [[ "${assay}" == atac ]]; then
    args+=(-p 1e-5 --min-length 200 --max-gap 30)
  else
    args+=(-q 0.05 --max-count 1 --min-length "${geometry}" --max-gap "${geometry}")
  fi
  log_command /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}"
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}" \
    >"${dir}/stdout.log" 2>"${dir}/stderr.log"
  record_run "${assay}" MACS3 1 "${repeat}" \
    "${dir}/out/${name}_peaks.narrowPeak" "${dir}/out/${name}_summits.bed" "${dir}/time.tsv"
}

run_rapidmacs() {
  local assay="$1" threads="$2" repeat="$3" input="$4" geometry="$5" name="$6"
  local dir="${OUT_ROOT}/runs/${assay}.rapidmacs.t${threads}.r${repeat}"
  mkdir -p "${dir}"
  local np="${dir}/${name}.narrowPeak" sm="${dir}/${name}.summits.bed" stage="${dir}/stages.tsv"
  local args=("${RAPIDMACS}" -f FRAG -t "${input}" -g hs -n "${name}"
              --macs3-frag-narrowpeak "${np}" --macs3-frag-summits "${sm}"
              --peak-caller-threads "${threads}" --profile-stages "${stage}")
  if [[ "${assay}" == atac ]]; then
    args+=(-p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30)
  else
    args+=(--macs3-frag-qvalue 0.05 --frag-pileup-max-count 1
           --bdgpeakcall-min-len "${geometry}" --bdgpeakcall-max-gap "${geometry}")
  fi
  log_command /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}"
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}" \
    >"${dir}/stdout.log" 2>"${dir}/stderr.log"
  record_run "${assay}" RapidMACS "${threads}" "${repeat}" "${np}" "${sm}" "${dir}/time.tsv"
  record_stage_run "${assay}" "${threads}" "${repeat}" "${stage}"
}

benchmark_assay() {
  local assay="$1" macs_input="$2" rapid_input="$3" geometry="$4" name="$5"
  local repeat warm threads
  if tool_enabled macs3; then
    for ((warm=1; warm<=WARMUPS; ++warm)); do
      echo "[warmup] ${assay} MACS3 ${warm}/${WARMUPS}"
      run_macs3 "${assay}" "warmup${warm}" "${macs_input}" "${geometry}" "${name}"
      sed -i '$d' "${RAW_TSV}"
    done
    for ((repeat=1; repeat<=REPEATS; ++repeat)); do
      echo "[benchmark] ${assay} MACS3 ${repeat}/${REPEATS}"
      run_macs3 "${assay}" "${repeat}" "${macs_input}" "${geometry}" "${name}"
    done
  fi
  if tool_enabled rapidmacs; then
    for threads in ${THREADS}; do
      for ((warm=1; warm<=WARMUPS; ++warm)); do
        echo "[warmup] ${assay} RapidMACS T=${threads} ${warm}/${WARMUPS}"
        run_rapidmacs "${assay}" "${threads}" "warmup${warm}" "${rapid_input}" "${geometry}" "${name}"
        sed -i '$d' "${RAW_TSV}"
        sed -i '$d' "${STAGE_RAW_TSV}"
      done
      for ((repeat=1; repeat<=REPEATS; ++repeat)); do
        echo "[benchmark] ${assay} RapidMACS T=${threads} ${repeat}/${REPEATS}"
        run_rapidmacs "${assay}" "${threads}" "${repeat}" "${rapid_input}" "${geometry}" "${name}"
      done
    done
  fi

  if tool_enabled macs3 && tool_enabled rapidmacs; then
    local ref_np ref_sm
    ref_np="$(md5sum "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_peaks.narrowPeak" | awk '{print $1}')"
    ref_sm="$(md5sum "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_summits.bed" | awk '{print $1}')"
    if ! awk -F '\t' -v a="${assay}" -v n="${ref_np}" -v s="${ref_sm}" '
      NR > 1 && $1 == a && $2 == "RapidMACS" && ($9 != n || $11 != s) {bad=1}
      END {exit bad}
    ' "${RAW_TSV}"; then
      echo "ERROR: ${assay} RapidMACS output differs from MACS3" >&2
      exit 1
    fi
    echo "[parity] ${assay}: all RapidMACS narrowPeak and summit files are byte-identical to MACS3"
  fi
}

if [[ " ${ASSAYS} " == *" atac "* ]]; then
  benchmark_assay atac "${atac_macs3}" "${atac_rapid}" 200 atac_treatment_only
fi
if [[ " ${ASSAYS} " == *" cutrun "* ]]; then
  benchmark_assay cutrun "${cutrun_macs3}" "${cutrun_rapid}" "${cutrun_geometry}" cutrun_treatment_only
fi
if [[ " ${ASSAYS} " == *" chipseq "* ]]; then
  benchmark_assay chipseq "${chipseq_macs3}" "${chipseq_rapid}" "${chipseq_geometry}" chipseq_treatment_only
fi

python3 "${REPO_ROOT}/paper/scripts/summarize_paper_benchmarks.py" "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
python3 "${REPO_ROOT}/paper/scripts/summarize_stage_benchmarks.py" "${STAGE_RAW_TSV}" "${OUT_ROOT}/stages_summary.tsv"

git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
{
  printf 'field\tvalue\n'
  printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_root\t%s\n' "${RUN_ROOT}"
  printf 'out_root\t%s\n' "${OUT_ROOT}"
  printf 'storage_source\t%s\n' "$(findmnt -n -o SOURCE -T "${RUN_ROOT}")"
  printf 'git_head\t%s\n' "$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  printf 'git_dirty\t%s\n' "$([[ -s "${OUT_ROOT}/provenance/git_status.txt" ]] && echo yes || echo no)"
  printf 'rapidmacs_sha256\t%s\n' "$(sha256sum "${RAPIDMACS}" | awk '{print $1}')"
  printf 'macs3_version\t%s\n' "$("${MACS3}" --version 2>&1 | head -n 1)"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'assays\t%s\n' "${ASSAYS}"
  printf 'atac_fragment_input\t%s\n' "${atac_rapid}"
  printf 'atac_fragment_source\t%s\n' "${ATAC_FRAGMENT_SOURCE:-Chromap realignment from public 10x FASTQs}"
} >"${OUT_ROOT}/provenance/metadata.tsv"
selected_inputs=()
if [[ " ${ASSAYS} " == *" atac "* ]]; then
  selected_inputs+=("${atac_rapid}")
fi
if [[ " ${ASSAYS} " == *" cutrun "* ]]; then
  selected_inputs+=("${cutrun_rapid}" "${cutrun_macs3}")
fi
if [[ " ${ASSAYS} " == *" chipseq "* ]]; then
  selected_inputs+=("${chipseq_rapid}" "${chipseq_macs3}")
fi
sha256sum "${selected_inputs[@]}" >"${OUT_ROOT}/provenance/input_sha256s.txt"
cp "${REPO_ROOT}/paper/scripts/"*.sh "${REPO_ROOT}/paper/scripts/fastq_sources.tsv" \
  "${OUT_ROOT}/provenance/"

echo "[complete] raw results: ${RAW_TSV}"
echo "[complete] summary: ${OUT_ROOT}/summary.tsv"
echo "[complete] provenance: ${OUT_ROOT}/provenance"
