#!/usr/bin/env bash
# Repeated, warm-cache peak-caller benchmark for the Applications Note.
# Alignment is deliberately excluded: each tool receives the same fragments/BAMs.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/.macs3-reference}"
ATAC_TREATMENT="${ATAC_TREATMENT:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k/fragments.tsv.gz}"
CHIP_TREATMENT="${CHIP_TREATMENT:-${DATA_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF686KKV.bam}"
CHIP_CONTROL="${CHIP_CONTROL:-${DATA_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF181CXT.bam}"
CUTRUN_TREATMENT="${CUTRUN_TREATMENT:-${DATA_ROOT}/cutrun/GSE172236/SRR14255054.bam}"
CUTRUN_CONTROL="${CUTRUN_CONTROL:-${DATA_ROOT}/cutrun/GSE172236/SRR14255053.bam}"
THREADS="${THREADS:-1 4 24}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
ASSAYS="${ASSAYS:-atac cutrun chipseq}"
# Select either implementation for incremental reruns. A RapidMACS-only run
# uses REFERENCE_RAW_TSV for its MACS3 parity hashes, so a code-only change
# does not require repeating the multi-minute reference caller.
TOOLS="${TOOLS:-macs3 rapidmacs}"
REFERENCE_RAW_TSV="${REFERENCE_RAW_TSV:-}"
# Defaults to the released build, which is what Table 1 and Supplementary S3
# report. Point MACS3 at macs3-3.0.3-no-control-dump-venv to reproduce the
# other half of Supplementary S4.
MACS3="${MACS3:-${DATA_ROOT}/macs3-3.0.3-venv/bin/macs3}"
OUT_ROOT="${OUT_ROOT:-${REPO_ROOT}/build/paper_benchmarks/$(date -u +%Y%m%dT%H%M%SZ)}"
RAW_TSV="${OUT_ROOT}/raw.tsv"
STAGE_RAW_TSV="${OUT_ROOT}/stages_raw.tsv"

mkdir -p "${OUT_ROOT}/runs"
printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"
printf 'assay\tthreads\trepeat\tcache_state\tload_s\tlibrary_pipeline_s\tbuild_workspace_s\ttreat_s\tlambda_s\tppois_regions_s\toutput_s\tprofiled_total_s\n' >"${STAGE_RAW_TSV}"

tool_enabled() {
  [[ " ${TOOLS} " == *" $1 "* ]]
}

required_tools=(awk md5sum /usr/bin/time)
if tool_enabled rapidmacs; then required_tools+=("${REPO_ROOT}/bin/rapidmacs"); fi
if tool_enabled macs3; then required_tools+=("${MACS3}"); fi
for tool in "${required_tools[@]}"; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: required tool not found: ${tool}" >&2
    exit 2
  fi
done
required_inputs=()
if [[ " ${ASSAYS} " == *" atac "* ]]; then required_inputs+=("${ATAC_TREATMENT}"); fi
if [[ " ${ASSAYS} " == *" chipseq "* ]]; then required_inputs+=("${CHIP_TREATMENT}" "${CHIP_CONTROL}"); fi
if [[ " ${ASSAYS} " == *" cutrun "* ]]; then required_inputs+=("${CUTRUN_TREATMENT}" "${CUTRUN_CONTROL}"); fi
for path in "${required_inputs[@]}"; do
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: benchmark input is missing: ${path}" >&2
    exit 2
  fi
done
if tool_enabled macs3 &&
   [[ "$("${MACS3}" --version 2>&1 | awk 'NR==1{print $NF}')" != 3.0.3 ]]; then
  echo "ERROR: paper parity benchmark requires MACS3 3.0.3" >&2
  exit 2
fi

bed3_md5() {
  awk -F '\t' -v OFS='\t' '!/^track/ && NF >= 3 {print $1,$2,$3}' "$1" | md5sum | awk '{print $1}'
}

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
             a, t, r, load, total-load, build, treat, lambda, ppois, output, total
    }
  ' "${stage_file}" >>"${STAGE_RAW_TSV}"
}

run_macs3() {
  local assay="$1" repeat="$2" treatment="$3" control="$4" format="$5" name="$6"
  local dir="${OUT_ROOT}/runs/${assay}.macs3.r${repeat}"
  mkdir -p "${dir}/out"
  local args=("${MACS3}" callpeak -t "${treatment}" -f "${format}" -g hs -n "${name}" --outdir "${dir}/out")
  if [[ -n "${control}" ]]; then args+=(-c "${control}" -q 0.05); fi
  if [[ "${assay}" == atac ]]; then
    args+=(-p 1e-5 --min-length 200 --max-gap 30)
  elif [[ -z "${control}" ]]; then
    args+=(-q 0.05)
  fi
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}" \
    >"${dir}/stdout" 2>"${dir}/stderr"
  record_run "${assay}" MACS3 1 "${repeat}" \
    "${dir}/out/${name}_peaks.narrowPeak" "${dir}/out/${name}_summits.bed" "${dir}/time.tsv"
}

run_rapidmacs() {
  local assay="$1" threads="$2" repeat="$3" treatment="$4" control="$5" format="$6" preset="$7" name="$8"
  local dir="${OUT_ROOT}/runs/${assay}.rapidmacs.t${threads}.r${repeat}"
  mkdir -p "${dir}"
  local prefix="${dir}/${name}"
  local stage_file="${dir}/stages.tsv"
  local args=("${REPO_ROOT}/bin/rapidmacs" --preset "${preset}" -f "${format}"
              -t "${treatment}" -g hs -n "${name}" --peak-caller-threads "${threads}"
              --profile-stages "${stage_file}")
  if [[ -n "${control}" ]]; then args+=(-c "${control}"); fi
  if [[ "${assay}" == atac ]]; then
    args+=(--macs3-frag-narrowpeak "${prefix}.narrowPeak"
           --macs3-frag-summits "${prefix}.summits.bed" -p 1e-5
           --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30)
  else
    args+=(-q 0.05 --out-prefix "${prefix}")
  fi
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" "${args[@]}" \
    >"${dir}/stdout" 2>"${dir}/stderr"
  record_run "${assay}" RapidMACS "${threads}" "${repeat}" \
    "${prefix}.narrowPeak" "${prefix}.summits.bed" "${dir}/time.tsv"
  record_stage_run "${assay}" "${threads}" "${repeat}" "${stage_file}"
}

benchmark_assay() {
  local assay="$1" treatment="$2" control="$3" format="$4" preset="$5" name="$6"
  local warm repeat threads
  if tool_enabled macs3; then
    echo "[benchmark] ${assay}: warming MACS3 input/cache"
    for ((warm = 1; warm <= WARMUPS; ++warm)); do
      run_macs3 "${assay}" "warmup${warm}" "${treatment}" "${control}" "${format}" "${name}"
      # Warmups are executed through the same checked path, then omitted from raw.tsv.
      sed -i '$d' "${RAW_TSV}"
    done
    for ((repeat = 1; repeat <= REPEATS; ++repeat)); do
      echo "[benchmark] ${assay}: MACS3 repeat ${repeat}/${REPEATS}"
      run_macs3 "${assay}" "${repeat}" "${treatment}" "${control}" "${format}" "${name}"
    done
  fi
  if tool_enabled rapidmacs; then
    for threads in ${THREADS}; do
      echo "[benchmark] ${assay}: warming RapidMACS T=${threads}"
      for ((warm = 1; warm <= WARMUPS; ++warm)); do
        run_rapidmacs "${assay}" "${threads}" "warmup${warm}" "${treatment}" "${control}" "${format}" "${preset}" "${name}"
        sed -i '$d' "${RAW_TSV}"
        sed -i '$d' "${STAGE_RAW_TSV}"
      done
      for ((repeat = 1; repeat <= REPEATS; ++repeat)); do
        echo "[benchmark] ${assay}: RapidMACS T=${threads} repeat ${repeat}/${REPEATS}"
        run_rapidmacs "${assay}" "${threads}" "${repeat}" "${treatment}" "${control}" "${format}" "${preset}" "${name}"
      done
    done
  fi

  if ! tool_enabled rapidmacs; then return; fi
  local ref_narrowpeak ref_bed3 ref_summits ref_summits_bed3
  if tool_enabled macs3; then
    ref_narrowpeak="$(md5sum "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_peaks.narrowPeak" | awk '{print $1}')"
    ref_bed3="$(awk -F '\t' -v OFS='\t' '{print $1,$2,$3}' "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_peaks.narrowPeak" | md5sum | awk '{print $1}')"
    ref_summits="$(md5sum "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_summits.bed" | awk '{print $1}')"
    ref_summits_bed3="$(bed3_md5 "${OUT_ROOT}/runs/${assay}.macs3.r1/out/${name}_summits.bed")"
  elif [[ -s "${REFERENCE_RAW_TSV}" ]]; then
    read -r ref_narrowpeak ref_bed3 ref_summits ref_summits_bed3 < <(
      awk -F '\t' -v a="${assay}" '$1 == a && $2 == "MACS3" {print $9, $10, $11, $12; exit}' \
        "${REFERENCE_RAW_TSV}")
    if [[ -z "${ref_narrowpeak:-}" ]]; then
      echo "ERROR: no MACS3 ${assay} reference in ${REFERENCE_RAW_TSV}" >&2
      exit 2
    fi
  else
    echo "ERROR: RapidMACS-only runs require REFERENCE_RAW_TSV" >&2
    exit 2
  fi
  if awk -F '\t' -v a="${assay}" -v n="${ref_narrowpeak}" -v b="${ref_bed3}" \
      -v s="${ref_summits}" -v sb="${ref_summits_bed3}" '
      NR > 1 && $1 == a && $2 == "RapidMACS" &&
        ($9 != n || $10 != b || $11 != s || $12 != sb) {bad=1}
      END {exit bad}
    ' "${RAW_TSV}"; then
    echo "[parity] ${assay}: every measured RapidMACS run matches MACS3 full narrowPeak + summits"
  else
    echo "ERROR: ${assay} parity changed during benchmarking" >&2
    exit 1
  fi
}

if [[ " ${ASSAYS} " == *" atac "* ]]; then
  benchmark_assay atac "${ATAC_TREATMENT}" "" FRAG shiftTAG atac_full
fi
if [[ " ${ASSAYS} " == *" cutrun "* ]]; then
  benchmark_assay cutrun "${CUTRUN_TREATMENT}" "${CUTRUN_CONTROL}" BAMPE cutrun cutrun_full
fi
if [[ " ${ASSAYS} " == *" chipseq "* ]]; then
  benchmark_assay chipseq "${CHIP_TREATMENT}" "${CHIP_CONTROL}" BAMPE chip chipseq_full
fi

python3 "${REPO_ROOT}/benchmarks/summarize_paper_benchmarks.py" \
  "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
python3 "${REPO_ROOT}/benchmarks/summarize_stage_benchmarks.py" \
  "${STAGE_RAW_TSV}" "${OUT_ROOT}/stages_summary.tsv"
echo "[benchmark] raw: ${RAW_TSV}"
echo "[benchmark] summary: ${OUT_ROOT}/summary.tsv"
echo "[benchmark] stage raw: ${STAGE_RAW_TSV}"
echo "[benchmark] stage summary: ${OUT_ROOT}/stages_summary.tsv"
