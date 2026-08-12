#!/usr/bin/env bash
# Reproduce upstream scATAC peak-calling defaults against a common real input.
# Source revisions audited 2026-08-11:
#   ArchR master:       6feec354ad6c8052ddbc4626a2ca2d858ed465bf
#   Signac master:      12d1ff82774890b08e55b3101f36fab83e674b88
#   pycisTopic main:    787ce422a37f5975b0ebb9e7b19eeaed44847501

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MACS3="${MACS3:-${REPO_ROOT}/.macs3-reference/macs3-3.0.3-venv/bin/macs3}"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
INPUT="${INPUT:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k_100k/fragments.tsv}"
OUT_ROOT="${OUT_ROOT:-${REPO_ROOT}/build/workflow_compatibility/$(date -u +%Y%m%dT%H%M%SZ)_tool_audit}"
TOOLS="${TOOLS:-archr pycistopic}"
# Optional completed Phase 3/4 audit root whose MACS3 outputs can be reused.
# This avoids rerunning an unchanged reference while preserving the exact
# reference command in commands.sh.
MACS_REFERENCE_ROOT="${MACS_REFERENCE_ROOT:-}"

mkdir -p "${OUT_ROOT}"
RESULTS="${OUT_ROOT}/results.tsv"
COMMANDS="${OUT_ROOT}/commands.sh"
printf 'tool\tmacs_peaks\trapid_peaks\tnarrowpeak_md5_equal\tbed3_equal\tsummits_md5_equal\n' >"${RESULTS}"
printf '#!/usr/bin/env bash\nset -euo pipefail\n\n' >"${COMMANDS}"
if [[ -n "${MACS_REFERENCE_ROOT}" ]]; then
  printf 'macs_reference_root\t%s\n' "${MACS_REFERENCE_ROOT}" \
    >"${OUT_ROOT}/reference.tsv"
fi

record_command() {
  printf '%q ' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

bed3_md5() {
  awk -F '\t' -v OFS='\t' '{print $1,$2,$3}' "$1" | md5sum | awk '{print $1}'
}

record_result() {
  local tool="$1" macs_np="$2" rapid_np="$3" macs_sm="$4" rapid_sm="$5"
  local np_equal=no bed3_equal=no sm_equal=no
  cmp -s "${macs_np}" "${rapid_np}" && np_equal=yes
  [[ "$(bed3_md5 "${macs_np}")" == "$(bed3_md5 "${rapid_np}")" ]] && bed3_equal=yes
  cmp -s "${macs_sm}" "${rapid_sm}" && sm_equal=yes
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${tool}" "$(wc -l <"${macs_np}")" \
    "$(wc -l <"${rapid_np}")" "${np_equal}" "${bed3_equal}" "${sm_equal}" >>"${RESULTS}"
}

run_archr() {
  local work="${OUT_ROOT}/archr"
  mkdir -p "${work}/macs3" "${work}/rapidmacs"
  local macs_cmd=("${MACS3}" callpeak -g 2700000000 --name archr
                  --treatment "${INPUT}" --outdir "${work}/macs3" --format BED
                  --call-summits --keep-dup all --nomodel --nolambda
                  --shift -75 --extsize 150 -q 0.1)
  local rapid_cmd=("${RAPIDMACS}" --macs-input-format BED
                   --macs-bed-input "${INPUT}" --outdir "${work}/rapidmacs"
                   --macs-bed-peaks-out "${work}/rapidmacs/archr_peaks.narrowPeak"
                   --macs-bed-summits-out "${work}/rapidmacs/archr_summits.bed"
                   --effective-genome-size 2700000000 --qvalue 0.1 --nomodel
                   --shift -75 --extsize 150 --keep-dup all --call-summits
                   --nolambda --llocal 10000 -n archr)
  record_command "${macs_cmd[@]}"
  local macs_np="${work}/macs3/archr_peaks.narrowPeak"
  local macs_sm="${work}/macs3/archr_summits.bed"
  if [[ -n "${MACS_REFERENCE_ROOT}" ]]; then
    macs_np="${MACS_REFERENCE_ROOT}/archr/macs3/archr_peaks.narrowPeak"
    macs_sm="${MACS_REFERENCE_ROOT}/archr/macs3/archr_summits.bed"
  else
    "${macs_cmd[@]}" >"${work}/macs3.stdout" 2>"${work}/macs3.stderr"
  fi
  record_command "${rapid_cmd[@]}"
  "${rapid_cmd[@]}" >"${work}/rapidmacs.stdout" 2>"${work}/rapidmacs.stderr"
  record_result archr "${macs_np}" \
    "${work}/rapidmacs/archr_peaks.narrowPeak" \
    "${macs_sm}" "${work}/rapidmacs/archr_summits.bed"
}

run_pycistopic() {
  local work="${OUT_ROOT}/pycistopic"
  mkdir -p "${work}/macs3" "${work}/rapidmacs"
  local macs_cmd=("${MACS3}" callpeak --treatment "${INPUT}" --name pycistopic
                  --outdir "${work}/macs3" --format BEDPE --gsize hs
                  --qvalue 0.05 --nomodel --shift 73 --extsize 146
                  --keep-dup all --call-summits --nolambda)
  local rapid_cmd=("${RAPIDMACS}" --preset cutrun -f BEDPE -t "${INPUT}"
                   -g hs -q 0.05 --nomodel --shift 73 --extsize 146
                   --keep-dup all --call-summits --nolambda -n pycistopic
                   --macs3-frag-narrowpeak
                   "${work}/rapidmacs/pycistopic_peaks.narrowPeak"
                   --macs3-frag-summits
                   "${work}/rapidmacs/pycistopic_summits.bed")
  record_command "${macs_cmd[@]}"
  local macs_np="${work}/macs3/pycistopic_peaks.narrowPeak"
  local macs_sm="${work}/macs3/pycistopic_summits.bed"
  if [[ -n "${MACS_REFERENCE_ROOT}" ]]; then
    macs_np="${MACS_REFERENCE_ROOT}/pycistopic/macs3/pycistopic_peaks.narrowPeak"
    macs_sm="${MACS_REFERENCE_ROOT}/pycistopic/macs3/pycistopic_summits.bed"
  else
    "${macs_cmd[@]}" >"${work}/macs3.stdout" 2>"${work}/macs3.stderr"
  fi
  record_command "${rapid_cmd[@]}"
  "${rapid_cmd[@]}" >"${work}/rapidmacs.stdout" 2>"${work}/rapidmacs.stderr"
  record_result pycistopic "${macs_np}" \
    "${work}/rapidmacs/pycistopic_peaks.narrowPeak" \
    "${macs_sm}" \
    "${work}/rapidmacs/pycistopic_summits.bed"
}

for tool in ${TOOLS}; do
  case "${tool}" in
    archr) run_archr ;;
    pycistopic) run_pycistopic ;;
    *) echo "Unsupported tool: ${tool}" >&2; exit 2 ;;
  esac
done

echo "Results: ${RESULTS}"
echo "Commands: ${COMMANDS}"
