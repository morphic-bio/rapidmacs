#!/usr/bin/env bash
# Reproduce the scATAC 10x row of Table 1 on a single public dataset.
#
# Downloads the 10x Cell Ranger ARC PBMC 3k ATAC fragments, calls peaks with
# released MACS3 3.0.3 and with RapidMACS at several thread counts, and checks
# that every output file is byte-identical to the published result. Exits
# non-zero if any checksum disagrees.
#
# Usage:
#   paper/scripts/run_10x_pbmc_example.sh
#
# Environment:
#   WORK_DIR    scratch + output directory (default ./10x_pbmc_example)
#   MACS3       macs3 3.0.3 executable (default: macs3 on PATH)
#   THREADS     RapidMACS thread counts (default "1 4 24")
#   REPEATS     timed repetitions per configuration (default 1)
#   WARMUPS     unrecorded warm-up runs per configuration (default 0)
#   SKIP_MACS3  set to 1 to time RapidMACS only and check against the
#               published checksums rather than against a local MACS3 run
#   FRAGMENTS   use an already-downloaded fragments file instead of fetching
#
# The paper's protocol is WARMUPS=1 REPEATS=3 and reports medians. The defaults
# here run each configuration once, which is enough to demonstrate parity and
# the order of magnitude but is not a replacement for the full protocol. MACS3
# takes roughly 12 minutes per run on the reference host.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${WORK_DIR:-$(pwd)/10x_pbmc_example}"
MACS3="${MACS3:-macs3}"
THREADS="${THREADS:-1 4 24}"
REPEATS="${REPEATS:-1}"
WARMUPS="${WARMUPS:-0}"
SKIP_MACS3="${SKIP_MACS3:-0}"
RAPIDMACS="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"

SOURCE_URL="https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_unsorted_3k/pbmc_unsorted_3k_atac_fragments.tsv.gz"
SOURCE_SHA256=0189df1de4894f10f48ffb905da4966cefc2466dc5fdd76ee2905244dd5a8dd7
SOURCE_BYTES=537308082

# Published result for this configuration (Table 1, row 1; Supplementary S1).
EXPECTED_PEAKS=50003
EXPECTED_NARROWPEAK_MD5=2c0de7c506514b0afef1e4339c386c9f
EXPECTED_SUMMITS_MD5=c3805d1d3744ed78e3d589a523ad0b09

# Reference wall times in seconds on the paper host (i9-13900KF, NVMe, warm
# cache, medians of three). Printed for orientation only; nothing asserts them,
# because timings are hardware-dependent.
REF_MACS3=747.55
declare -A REF_RAPID=([1]=58.93 [4]=24.44 [24]=16.35)

die() { echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------- preflight --
for tool in curl sha256sum md5sum awk gzip /usr/bin/time; do
  command -v "${tool}" >/dev/null 2>&1 || die "missing required tool: ${tool}"
done
[[ -x "${RAPIDMACS}" ]] || die "rapidmacs not built at ${RAPIDMACS} (run: make -j)"

if [[ "${SKIP_MACS3}" != 1 ]]; then
  command -v "${MACS3}" >/dev/null 2>&1 || die "macs3 not found; set MACS3=/path/to/macs3, or SKIP_MACS3=1"
  macs3_version="$("${MACS3}" --version 2>&1 | awk 'NR == 1 {print $NF}')"
  [[ "${macs3_version}" == 3.0.3 ]] \
    || die "expected MACS3 3.0.3, found ${macs3_version}. See paper/scripts/setup_macs3_reference.sh"
fi

mkdir -p "${WORK_DIR}/runs"

# ------------------------------------------------------------------- input ---
FRAGMENTS="${FRAGMENTS:-${WORK_DIR}/pbmc_unsorted_3k_atac_fragments.tsv.gz}"
if [[ ! -s "${FRAGMENTS}" ]]; then
  echo "[fetch] ${SOURCE_URL}"
  curl -fL --retry 5 --retry-delay 2 --continue-at - -o "${FRAGMENTS}.partial" "${SOURCE_URL}"
  mv "${FRAGMENTS}.partial" "${FRAGMENTS}"
fi
actual_bytes="$(stat -c %s "${FRAGMENTS}")"
[[ "${actual_bytes}" == "${SOURCE_BYTES}" ]] \
  || die "fragments file is ${actual_bytes} bytes, expected ${SOURCE_BYTES}"
echo "[verify] sha256 of input"
printf '%s  %s\n' "${SOURCE_SHA256}" "${FRAGMENTS}" | sha256sum -c - >/dev/null \
  || die "fragments checksum mismatch: not the file the paper used"
gzip -t "${FRAGMENTS}"

# --------------------------------------------------------------------- run ---
RESULTS="${WORK_DIR}/results.tsv"
printf 'tool\tthreads\trepeat\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tsummits_md5\n' >"${RESULTS}"

record() {
  local tool="$1" threads="$2" repeat="$3" dir="$4" narrowpeak="$5" summits="$6"
  local wall rss peaks np_md5 su_md5
  wall="$(awk -F '\t' 'NR == 1 {print $1}' "${dir}/time.tsv")"
  rss="$(awk -F '\t' 'NR == 1 {print $2}' "${dir}/time.tsv")"
  peaks="$(grep -vc '^#' "${narrowpeak}" || true)"
  np_md5="$(md5sum "${narrowpeak}" | awk '{print $1}')"
  su_md5="$(md5sum "${summits}" | awk '{print $1}')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${tool}" "${threads}" "${repeat}" "${wall}" "${rss}" "${peaks}" "${np_md5}" "${su_md5}" >>"${RESULTS}"
  [[ "${repeat}" == warmup* ]] && return 0
  [[ "${peaks}" == "${EXPECTED_PEAKS}" ]] \
    || die "${tool} T=${threads} called ${peaks} peaks, expected ${EXPECTED_PEAKS}"
  [[ "${np_md5}" == "${EXPECTED_NARROWPEAK_MD5}" ]] \
    || die "${tool} T=${threads} narrowPeak MD5 ${np_md5} != published ${EXPECTED_NARROWPEAK_MD5}"
  [[ "${su_md5}" == "${EXPECTED_SUMMITS_MD5}" ]] \
    || die "${tool} T=${threads} summits MD5 ${su_md5} != published ${EXPECTED_SUMMITS_MD5}"
  return 0
}

run_macs3() {
  local repeat="$1" dir="${WORK_DIR}/runs/macs3.${1}"
  rm -rf "${dir}"; mkdir -p "${dir}/out"
  echo "[macs3 ${repeat}] this takes about 12 minutes on the reference host"
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" \
    "${MACS3}" callpeak -t "${FRAGMENTS}" -f FRAG -g hs -n atac \
      --outdir "${dir}/out" -p 1e-5 --min-length 200 --max-gap 30 \
      >"${dir}/stdout" 2>"${dir}/stderr"
  record MACS3 1 "${repeat}" "${dir}" \
    "${dir}/out/atac_peaks.narrowPeak" "${dir}/out/atac_summits.bed"
}

run_rapidmacs() {
  local threads="$1" repeat="$2" dir="${WORK_DIR}/runs/rapidmacs.t${1}.${2}"
  rm -rf "${dir}"; mkdir -p "${dir}"
  echo "[rapidmacs T=${threads} ${repeat}]"
  /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" \
    "${RAPIDMACS}" --preset shiftTAG -f FRAG -t "${FRAGMENTS}" -g hs -n atac \
      --peak-caller-threads "${threads}" -p 1e-5 \
      --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
      --macs3-frag-narrowpeak "${dir}/atac.narrowPeak" \
      --macs3-frag-summits "${dir}/atac.summits.bed" \
      >"${dir}/stdout" 2>"${dir}/stderr"
  record RapidMACS "${threads}" "${repeat}" "${dir}" \
    "${dir}/atac.narrowPeak" "${dir}/atac.summits.bed"
}

if [[ "${SKIP_MACS3}" != 1 ]]; then
  for ((w = 1; w <= WARMUPS; w++)); do run_macs3 "warmup${w}"; done
  for ((r = 1; r <= REPEATS; r++)); do run_macs3 "${r}"; done
fi
for threads in ${THREADS}; do
  for ((w = 1; w <= WARMUPS; w++)); do run_rapidmacs "${threads}" "warmup${w}"; done
  for ((r = 1; r <= REPEATS; r++)); do run_rapidmacs "${threads}" "${r}"; done
done

# ------------------------------------------------------------------ report ---
echo
echo "10x PBMC 3k ATAC fragments (53,015,732 rows; 69,137,946 fragments)"
echo
ref_pairs="MACS3:1=${REF_MACS3}"
for threads in ${THREADS}; do
  [[ -n "${REF_RAPID[${threads}]:-}" ]] && ref_pairs+=",RapidMACS:${threads}=${REF_RAPID[${threads}]}"
done

awk -F '\t' -v ref_macs3="${REF_MACS3}" -v ref_pairs="${ref_pairs}" '
  BEGIN {
    split(ref_pairs, pairs, ",")
    for (i in pairs) { split(pairs[i], kv, "="); ref[kv[1]] = kv[2] }
  }
  NR == 1 { next }
  $3 ~ /^warmup/ { next }
  { key = $1 ":" $2; n[key]++; sum[key] += $4; tool[key] = $1; thr[key] = $2 }
  END {
    printf "-1|  %-11s %7s  %9s  %8s  %11s\n", "tool", "threads", "wall (s)", "speedup", "paper (s)"
    base = ("MACS3:1" in sum) ? sum["MACS3:1"] / n["MACS3:1"] : ref_macs3
    for (key in n) {
      mean = sum[key] / n[key]
      order = (tool[key] == "MACS3") ? 0 : thr[key]
      paper = (key in ref) ? sprintf("%11.2f", ref[key]) : sprintf("%11s", "-")
      printf "%s|  %-11s %7s  %9.2f  %7.1fx  %s\n", order, tool[key], thr[key], mean, base / mean, paper
    }
  }
' "${RESULTS}" | sort -t'|' -k1,1n | cut -d'|' -f2-
echo
if [[ "${SKIP_MACS3}" == 1 ]]; then
  echo "  Speedups are against the paper's MACS3 time, not a local run (SKIP_MACS3=1)."
fi
echo "  \"paper (s)\" is the reference host (i9-13900KF, NVMe, warm cache, median of three)."
echo "PASS: every run produced ${EXPECTED_PEAKS} peaks with the published"
echo "      narrowPeak and summits MD5 sums."
echo "Per-run detail: ${RESULTS}"
