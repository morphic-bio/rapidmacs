#!/usr/bin/env bash
# Parity benchmark vs upstream MACS3 (Python) and MACS2 (Python).
#
# The validated comparison is:
#   RapidMACS rapidmacs --macs3-frag-narrowpeak ... --bdgpeakcall-cutoff 5
#                       --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30
#   vs
#   macs3 callpeak -f FRAG -p 1e-5 --min-length 200 --max-gap 30
#
# Both consume the same fragments TSV. Threshold semantics are aligned:
#   --bdgpeakcall-cutoff 5  ⇔  -log10(p) > 5  ⇔  -p 1e-5
# Complete narrowPeak and summit files must be byte-identical to MACS3.
# BED3 statistics below also describe the MACS2 comparison.
#
# This script also runs MACS2 with the same -f BEDPE-equivalent path for a
# version-bump noise-floor reference (MACS3 vs MACS2 on the same input).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FORMAT="${FORMAT:-FRAG}"
TREATMENT="${TREATMENT:-${FIXTURE_FRAGS:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k_100k/fragments.tsv}}"
CONTROL="${CONTROL:-}"
OUTDIR="${OUTDIR:-${REPO_ROOT}/build/parity_vs_macs}"
RAPIDMACS="${REPO_ROOT}/bin/rapidmacs"
MACS3="${MACS3:-macs3}"
MACS2="${MACS2:-macs2}"
GENOME_SIZE="${GENOME_SIZE:-hs}"
QVALUE="${QVALUE:-0.05}"
NAME="${NAME:-parity}"
PRESET="${PRESET:-}"

usage() {
  cat >&2 <<EOF
Usage: $0 [--format FRAG|BAMPE|BEDPE] --treatment PATH [--control PATH]
          [--preset shiftTAG|chip|cutrun] [--outdir PATH] [--name NAME]

FRAG retains the historical -p 1e-5/min-length 200/max-gap 30 comparison.
BAMPE/BEDPE use -q ${QVALUE} and MACS3's paired-end dynamic peak geometry.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -f|--format) FORMAT="$2"; shift 2 ;;
    -t|--treatment|--input) TREATMENT="$2"; shift 2 ;;
    -c|--control) CONTROL="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -n|--name) NAME="$2"; shift 2 ;;
    -g|--gsize) GENOME_SIZE="$2"; shift 2 ;;
    -q|--qvalue) QVALUE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done
FORMAT="${FORMAT^^}"
if [[ "${FORMAT}" != FRAG && "${FORMAT}" != BAMPE && "${FORMAT}" != BEDPE ]]; then
  echo "ERROR: --format must be FRAG, BAMPE, or BEDPE" >&2
  exit 2
fi
if [[ ! -s "${TREATMENT}" ]]; then
  echo "ERROR: treatment not found: ${TREATMENT}" >&2
  exit 2
fi
if [[ -n "${CONTROL}" && ! -s "${CONTROL}" ]]; then
  echo "ERROR: control not found: ${CONTROL}" >&2
  exit 2
fi
if [[ -z "${PRESET}" ]]; then
  if [[ "${FORMAT}" == FRAG ]]; then PRESET=shiftTAG; else PRESET=cutrun; fi
fi
TREATMENT="$(realpath "${TREATMENT}")"
if [[ -n "${CONTROL}" ]]; then CONTROL="$(realpath "${CONTROL}")"; fi

if [[ ! -x "${RAPIDMACS}" ]]; then
  echo "ERROR: ${RAPIDMACS} not built. Run 'make' first." >&2
  exit 2
fi
for tool in "${MACS3}" "${MACS2}" bedtools awk sort cmp; do
  if ! command -v "${tool}" >/dev/null 2>&1; then
    echo "ERROR: ${tool} not on PATH" >&2
    exit 2
  fi
done

case "${OUTDIR}" in
  ""|/|"${REPO_ROOT}") echo "ERROR: refusing unsafe OUTDIR=${OUTDIR}" >&2; exit 2 ;;
esac
rm -rf -- "${OUTDIR}"
mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

run_with_time() {
  local label="$1"; shift
  (cd "${OUTDIR}" &&
    /usr/bin/time -v -o "${OUTDIR}/${label}.time" "$@" \
      > "${OUTDIR}/${label}.stdout" 2> "${OUTDIR}/${label}.stderr")
}

# 1. RapidMACS and MACS3 v3.0.3 with aligned paired-end semantics.
lib_prefix="${OUTDIR}/rapidmacs"
lib_args=("${RAPIDMACS}" --preset "${PRESET}" -f "${FORMAT}" -t "${TREATMENT}"
          -g "${GENOME_SIZE}" -n "${NAME}")
macs3_args=("${MACS3}" callpeak -t "${TREATMENT}" -f "${FORMAT}"
            -g "${GENOME_SIZE}" -n "${NAME}")
if [[ -n "${CONTROL}" ]]; then
  lib_args+=(-c "${CONTROL}")
  macs3_args+=(-c "${CONTROL}")
fi
if [[ "${FORMAT}" == FRAG ]]; then
  lib_args+=(--macs3-frag-narrowpeak "${lib_prefix}.narrowPeak"
             --macs3-frag-summits "${lib_prefix}.summits.bed"
             -p 1e-5 --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30)
  macs3_args+=(-p 1e-5 --min-length 200 --max-gap 30)
else
  lib_args+=(-q "${QVALUE}" --out-prefix "${lib_prefix}")
  macs3_args+=(-q "${QVALUE}")
fi
macs3_args+=(--outdir "${OUTDIR}/macs3_out")

echo "[parity] RapidMACS ${FORMAT} treatment${CONTROL:+/control} ..."
run_with_time rapidmacs "${lib_args[@]}"
echo "[parity] MACS3 3.0.3 ${FORMAT} treatment${CONTROL:+/control} ..."
run_with_time macs3 "${macs3_args[@]}"

# 2. MACS2 noise-floor reference on the same paired-end input.
macs2_treatment="${TREATMENT}"
macs2_format="${FORMAT}"
if [[ "${FORMAT}" == FRAG ]]; then
  macs2_treatment="${OUTDIR}/frags.bedpe"
  macs2_format=BEDPE
  awk -v OFS='\t' '{ for (i = 0; i < $5; i++) print $1,$2,$3 }' \
    "${TREATMENT}" >"${macs2_treatment}"
fi
macs2_args=("${MACS2}" callpeak -t "${macs2_treatment}" -f "${macs2_format}"
            -g "${GENOME_SIZE}" -n "${NAME}")
if [[ -n "${CONTROL}" ]]; then macs2_args+=(-c "${CONTROL}"); fi
if [[ "${FORMAT}" == FRAG ]]; then
  macs2_args+=(-p 1e-5 --min-length 200 --max-gap 30 --nomodel)
else
  macs2_args+=(-q "${QVALUE}")
fi
macs2_args+=(--outdir "${OUTDIR}/macs2_out")
echo "[parity] MACS2 ${macs2_format} noise-floor reference ..."
run_with_time macs2 "${macs2_args[@]}"

P_RAPIDMACS="${lib_prefix}.narrowPeak"
P_MACS3="${OUTDIR}/macs3_out/${NAME}_peaks.narrowPeak"
P_MACS2="${OUTDIR}/macs2_out/${NAME}_peaks.narrowPeak"
S_RAPIDMACS="${lib_prefix}.summits.bed"
S_MACS3="${OUTDIR}/macs3_out/${NAME}_summits.bed"

# BED3 extraction for overlap statistics; full-file identity is checked below.
bed3() {
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^#/ && NF>=3 { print $1, $2+0, $3+0 }' "$1" |
    LC_ALL=C sort -k1,1 -k2,2n
}
bed3 "${P_RAPIDMACS}" > "${OUTDIR}/rapidmacs.bed3"
bed3 "${P_MACS3}"    > "${OUTDIR}/macs3.bed3"
bed3 "${P_MACS2}"    > "${OUTDIR}/macs2.bed3"
bed3 "${S_RAPIDMACS}" > "${OUTDIR}/rapidmacs.summits.bed3"
bed3 "${S_MACS3}"    > "${OUTDIR}/macs3.summits.bed3"

parse_time() {
  awk -F': ' '
    /Maximum resident set size/ { rss=$2 }
    /Elapsed \(wall clock\) time/ { wall=$2 }
    END { printf "%s\t%s", (rss?rss:"NA"), (wall?wall:"NA") }
  ' "$1"
}

echo
echo "=== Per-tool stats ==="
{
  printf "tool\tpeaks\trss_kb\twall\n"
  printf "RapidMACS (%s)\t%d\t%s\n" "${FORMAT}" "$(wc -l < "${P_RAPIDMACS}")" "$(parse_time "${OUTDIR}/rapidmacs.time")"
  printf "MACS3 (%s)\t%d\t%s\n" "${FORMAT}" "$(wc -l < "${P_MACS3}")" "$(parse_time "${OUTDIR}/macs3.time")"
  printf "MACS2 (%s)\t%d\t%s\n" "${macs2_format}" "$(wc -l < "${P_MACS2}")" "$(parse_time "${OUTDIR}/macs2.time")"
} | column -t -s $'\t' | tee "${OUTDIR}/per_tool_stats.tsv"

pair_bed3() {
  local label="$1" A="$2" B="$3"
  local a_n b_n jac a_in_b b_in_a a_pct b_pct byte_id
  a_n=$(wc -l < "${A}")
  b_n=$(wc -l < "${B}")
  jac=$(bedtools jaccard -a "${A}" -b "${B}" 2>/dev/null | tail -1 | awk '{print $3}')
  a_in_b=$(bedtools intersect -u -a "${A}" -b "${B}" 2>/dev/null | wc -l)
  b_in_a=$(bedtools intersect -u -a "${B}" -b "${A}" 2>/dev/null | wc -l)
  a_pct=$(awk -v n="${a_in_b}" -v d="${a_n}" 'BEGIN{printf "%.3f", (d>0)?n/d:0}')
  b_pct=$(awk -v n="${b_in_a}" -v d="${b_n}" 'BEGIN{printf "%.3f", (d>0)?n/d:0}')
  byte_id=$(cmp -s "${A}" "${B}" && echo YES || echo NO)
  printf "%s\t%d\t%d\t%s\t%s\t%s\t%s\n" "${label}" "${a_n}" "${b_n}" "${jac}" "${a_pct}" "${b_pct}" "${byte_id}"
}

echo
echo "=== BED3 pairwise comparison (chrom/start/end only) ==="
{
  printf "pair\tA_peaks\tB_peaks\tbp_jaccard\tA_in_B_frac\tB_in_A_frac\tbyte_identical\n"
  pair_bed3 "RapidMACS vs MACS3 (validated)" "${OUTDIR}/rapidmacs.bed3" "${OUTDIR}/macs3.bed3"
  pair_bed3 "MACS3 vs MACS2 (noise floor)"  "${OUTDIR}/macs3.bed3"    "${OUTDIR}/macs2.bed3"
  pair_bed3 "RapidMACS vs MACS2"             "${OUTDIR}/rapidmacs.bed3" "${OUTDIR}/macs2.bed3"
} | column -t -s $'\t' | tee "${OUTDIR}/pairwise_bed3.tsv"

if [[ "${FORMAT}" == FRAG ]]; then
  cmp -s "${P_RAPIDMACS}" "${P_MACS3}" || {
    echo "ERROR: FRAG narrowPeak differs from MACS3" >&2
    exit 1
  }
  cmp -s "${S_RAPIDMACS}" "${S_MACS3}" || {
    echo "ERROR: FRAG summits differ from MACS3" >&2
    exit 1
  }
  parity_note="FRAG full narrowPeak and summits are byte-identical to MACS3"
else
  cmp -s "${P_RAPIDMACS}" "${P_MACS3}" || {
    echo "ERROR: controlled narrowPeak differs from MACS3" >&2
    exit 1
  }
  cmp -s "${S_RAPIDMACS}" "${S_MACS3}" || {
    echo "ERROR: controlled summits differ from MACS3" >&2
    exit 1
  }
  parity_note="controlled narrowPeak and summit files are byte-identical"
fi

echo
echo "Outputs in ${OUTDIR}/"
echo
echo "Parity: ${parity_note}."
