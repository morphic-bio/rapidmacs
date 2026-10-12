#!/usr/bin/env bash
# Check capped counts against original MACS3, including an uncapped control.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
GOLDEN="${REPO_ROOT}/tests/golden/capped_counts"
INPUT="${REPO_ROOT}/tests/fixtures/capped_counts.fragments.tsv"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rapidmacs_max_count.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT
if [[ "${UPDATE_GOLDENS:-0}" == 1 && -z "${MACS3:-}" ]]; then
  echo "ERROR: UPDATE_GOLDENS requires MACS3 pointing to original MACS3 3.0.3" >&2
  exit 2
fi
if [[ -n "${MACS3:-}" ]]; then
  [[ "$("${MACS3}" --version)" == 'macs3 3.0.3' ]]
fi

for cap in 0 1 2 5; do
  candidate="${WORK}/${cap}"
  mkdir -p "${candidate}"
  "${CLI}" -f FRAG -t "${INPUT}" -g 2913022398 -p 1e-5 -n capped \
    --bdgpeakcall-min-len 200 --bdgpeakcall-max-gap 30 \
    --frag-pileup-max-count "${cap}" \
    --macs3-frag-narrowpeak "${candidate}/capped_peaks.narrowPeak" \
    --macs3-frag-summits "${candidate}/capped_summits.bed" \
    >"${candidate}/run.log" 2>&1
  if [[ -n "${MACS3:-}" ]]; then
    reference="${WORK}/reference_${cap}"
    "${MACS3}" callpeak -f FRAG -t "${INPUT}" -g 2913022398 -p 1e-5 \
      --min-length 200 --max-gap 30 --max-count "${cap}" -n capped \
      --outdir "${reference}" >"${candidate}/reference.log" 2>&1
  fi
  for suffix in peaks.narrowPeak summits.bed; do
    file="capped_${suffix}"
    test -s "${candidate}/${file}"
    if [[ "${UPDATE_GOLDENS:-0}" == 1 ]]; then
      mkdir -p "${GOLDEN}/${cap}"
      cp "${reference}/${file}" "${GOLDEN}/${cap}/${file}"
    fi
    if [[ -n "${MACS3:-}" ]]; then
      cmp "${GOLDEN}/${cap}/${file}" "${reference}/${file}"
    fi
    cmp "${GOLDEN}/${cap}/${file}" "${candidate}/${file}"
  done
  echo "PASS: FRAG max-count ${cap}, $(wc -l < "${candidate}/capped_peaks.narrowPeak") peaks; narrowPeak and summits byte-identical"
done
