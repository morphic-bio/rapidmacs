#!/usr/bin/env bash
# Full-file MACS3 3.0.3 parity matrix for paired-end --keep-dup support.
# This is a compatibility check, not a timing benchmark.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${DATA_ROOT:-${DATA_DIR:-${REPO_ROOT}/data}}"
MACS3="${MACS3:-${REPO_ROOT}/.macs3-reference/macs3-3.0.3-venv/bin/macs3}"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
OUT_ROOT="${OUT_ROOT:-${REPO_ROOT}/build/workflow_compatibility/$(date -u +%Y%m%dT%H%M%SZ)_keepdup}"
ASSAYS="${ASSAYS:-chipseq cutrun}"
FORMATS="${FORMATS:-BAMPE BEDPE}"
KEEP_VALUES="${KEEP_VALUES:-1 2 all}"

CHIP_TREATMENT="${CHIP_TREATMENT:-${DATA_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF686KKV.bam}"
CHIP_CONTROL="${CHIP_CONTROL:-${DATA_ROOT}/chipseq/ENCSR634OAQ_ENCSR409FGC/ENCFF181CXT.bam}"
CUTRUN_TREATMENT="${CUTRUN_TREATMENT:-${DATA_ROOT}/cutrun/GSE172236/SRR14255054.bam}"
CUTRUN_CONTROL="${CUTRUN_CONTROL:-${DATA_ROOT}/cutrun/GSE172236/SRR14255053.bam}"

mkdir -p "${OUT_ROOT}/derived_bedpe" "${OUT_ROOT}/runs"
RESULTS="${OUT_ROOT}/results.tsv"
COMMANDS="${OUT_ROOT}/commands.sh"
if [[ ! -e "${RESULTS}" ]]; then
  printf 'assay\tformat\tkeep_dup\tpeaks\tnarrowpeak_md5\tsummits_md5\tbyte_identical\n' >"${RESULTS}"
fi
if [[ ! -e "${COMMANDS}" ]]; then
  printf '#!/usr/bin/env bash\nset -euo pipefail\n\n' >"${COMMANDS}"
fi

record_command() {
  printf '%q ' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

derive_bedpe() {
  local source_bam="$1" output_bedpe="$2"
  if [[ -s "${output_bedpe}" ]]; then return; fi
  local partial="${output_bedpe}.partial"
  if [[ -e "${partial}" ]]; then
    echo "ERROR: incomplete BEDPE already exists: ${partial}" >&2
    return 1
  fi
  local cmd=("${MACS3}" randsample -i "${source_bam}" -f BAMPE -p 100 -o "${partial}")
  record_command "${cmd[@]}"
  "${cmd[@]}"
  record_command mv "${partial}" "${output_bedpe}"
  mv "${partial}" "${output_bedpe}"
}

run_case() {
  local assay="$1" preset="$2" format="$3" keep="$4" treatment="$5" control="$6"
  local slug="${assay}.${format}.keep_${keep}"
  local run_dir="${OUT_ROOT}/runs/${slug}"
  local name="${assay}_${format}_keep_${keep}"
  mkdir -p "${run_dir}/macs3" "${run_dir}/rapidmacs"

  local macs_np="${run_dir}/macs3/${name}_peaks.narrowPeak"
  local macs_sm="${run_dir}/macs3/${name}_summits.bed"
  local rapid_np="${run_dir}/rapidmacs/${name}.narrowPeak"
  local rapid_sm="${run_dir}/rapidmacs/${name}.summits.bed"
  if [[ -s "${macs_np}" && -s "${macs_sm}" && -s "${rapid_np}" &&
        -s "${rapid_sm}" ]] &&
     awk -F '\t' -v a="${assay}" -v f="${format}" -v k="${keep}" \
       '$1 == a && $2 == f && $3 == k && $7 == "yes" {found=1} END {exit !found}' \
       "${RESULTS}"; then
    echo "SKIP: ${slug} already passed"
    return 0
  fi

  local macs_cmd=("${MACS3}" callpeak -t "${treatment}" -c "${control}"
                  -f "${format}" -g hs -q 0.05 --keep-dup "${keep}"
                  -n "${name}" --outdir "${run_dir}/macs3")
  local rapid_cmd=("${RAPIDMACS}" --preset "${preset}" -f "${format}"
                   -t "${treatment}" -c "${control}" -g hs -q 0.05
                   --keep-dup "${keep}" -n "${name}"
                   --out-prefix "${run_dir}/rapidmacs/${name}")
  if [[ ! -s "${macs_np}" || ! -s "${macs_sm}" ]]; then
    record_command "${macs_cmd[@]}"
    "${macs_cmd[@]}" >"${run_dir}/macs3.stdout" 2>"${run_dir}/macs3.stderr"
  else
    echo "REUSE: existing MACS3 reference for ${slug}"
  fi
  record_command "${rapid_cmd[@]}"
  "${rapid_cmd[@]}" >"${run_dir}/rapidmacs.stdout" 2>"${run_dir}/rapidmacs.stderr"

  local identical=no
  if cmp -s "${macs_np}" "${rapid_np}" && cmp -s "${macs_sm}" "${rapid_sm}"; then
    identical=yes
  fi
  local result_tmp="${RESULTS}.tmp"
  awk -F '\t' -v a="${assay}" -v f="${format}" -v k="${keep}" \
    'NR == 1 || !($1 == a && $2 == f && $3 == k)' "${RESULTS}" >"${result_tmp}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${assay}" "${format}" "${keep}" "$(wc -l <"${macs_np}")" \
    "$(md5sum "${macs_np}" | awk '{print $1}')" \
    "$(md5sum "${macs_sm}" | awk '{print $1}')" "${identical}" >>"${result_tmp}"
  mv "${result_tmp}" "${RESULTS}"
  if [[ "${identical}" != yes ]]; then
    echo "FAIL: ${slug} differs from MACS3" >&2
    return 1
  fi
  echo "PASS: ${slug}"
}

for assay in ${ASSAYS}; do
  case "${assay}" in
    chipseq)
      preset=chip
      bam_treatment="${CHIP_TREATMENT}"
      bam_control="${CHIP_CONTROL}"
      ;;
    cutrun)
      preset=cutrun
      bam_treatment="${CUTRUN_TREATMENT}"
      bam_control="${CUTRUN_CONTROL}"
      ;;
    *)
      echo "Unsupported assay: ${assay}" >&2
      exit 2
      ;;
  esac
  for format in ${FORMATS}; do
    treatment="${bam_treatment}"
    control="${bam_control}"
    if [[ "${format}" == BEDPE ]]; then
      treatment="${OUT_ROOT}/derived_bedpe/${assay}.treatment.bedpe"
      control="${OUT_ROOT}/derived_bedpe/${assay}.control.bedpe"
      derive_bedpe "${bam_treatment}" "${treatment}"
      derive_bedpe "${bam_control}" "${control}"
    elif [[ "${format}" != BAMPE ]]; then
      echo "Unsupported format: ${format}" >&2
      exit 2
    fi
    for keep in ${KEEP_VALUES}; do
      run_case "${assay}" "${preset}" "${format}" "${keep}" "${treatment}" "${control}"
    done
  done
done

echo "Results: ${RESULTS}"
echo "Commands: ${COMMANDS}"
