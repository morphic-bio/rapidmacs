#!/usr/bin/env bash
# Parity smoke for the MACS3 FRAG path (--macs3-frag-narrowpeak +
# --macs3-frag-summits) — the path libchromap actually calls in production
# via RunMacs3FragPeakPipelineFromFragments. Both narrowPeak and summits
# must be byte-identical to the frozen MACS3 3.0.3 output, including peak
# names and every score column.
#
# The CLI uses the sorted streaming FRAG workspace for this invocation; the
# smoke guards that production path against output regressions.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FIXTURE_FRAGS="${FIXTURE_FRAGS:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k_100k/fragments.tsv}"
GOLDEN_DIR="${REPO_ROOT}/tests/golden/100k_macs3_frag"
GOLDEN_NP="${GOLDEN_DIR}/macs3_frag.narrowPeak"
GOLDEN_SM="${GOLDEN_DIR}/macs3_frag.summits.bed"
CLI="${REPO_ROOT}/bin/rapidmacs"

if [[ ! -x "${CLI}" ]]; then
  echo "ERROR: ${CLI} not built. Run 'make' first." >&2
  exit 2
fi
if [[ ! -s "${FIXTURE_FRAGS}" ]]; then
  echo "SKIP: 100K fragments fixture not present at ${FIXTURE_FRAGS}"
  echo "      This fixture is a Chromap realignment and is not redistributed."
  echo "      Point FIXTURE_FRAGS at a fragments TSV to run this test:"
  echo "        FIXTURE_FRAGS=/path/to/fragments.tsv make test-macs3-frag-100k"
  exit 0
fi
if [[ ! -s "${GOLDEN_NP}" || ! -s "${GOLDEN_SM}" ]]; then
  echo "ERROR: MACS3 FRAG golden missing at ${GOLDEN_DIR}" >&2
  exit 2
fi

WORK="$(mktemp -d /tmp/rapidmacs_macs3_frag_smoke.XXXXXX)"
trap 'rm -rf "${WORK}"' EXIT

echo "[smoke:macs3-frag] running rapidmacs --macs3-frag-* on 100K fixture ..."
/usr/bin/time -v -o "${WORK}/time.log" "${CLI}" \
  -i "${FIXTURE_FRAGS}" \
  -n atac_full \
  --macs3-frag-narrowpeak "${WORK}/macs3_frag.narrowPeak" \
  --macs3-frag-summits "${WORK}/macs3_frag.summits.bed" \
  --bin-width 50 --ext-size 150 --local-window 200 \
  --fdr 0.05 --p-value 0.01 --merge-gap 0 > "${WORK}/run.log" 2>&1

fail=0
if ! diff -q "${WORK}/macs3_frag.narrowPeak" "${GOLDEN_NP}" >/dev/null; then
  echo "FAIL: MACS3 FRAG narrowPeak diverged from golden" >&2
  diff "${WORK}/macs3_frag.narrowPeak" "${GOLDEN_NP}" | head -20 >&2
  fail=1
fi
if ! diff -q "${WORK}/macs3_frag.summits.bed" "${GOLDEN_SM}" >/dev/null; then
  echo "FAIL: MACS3 FRAG summits diverged from golden" >&2
  diff "${WORK}/macs3_frag.summits.bed" "${GOLDEN_SM}" | head -20 >&2
  fail=1
fi
if [[ ${fail} -ne 0 ]]; then exit 1; fi

NP_MD5="$(md5sum "${WORK}/macs3_frag.narrowPeak" | awk '{print $1}')"
SM_MD5="$(md5sum "${WORK}/macs3_frag.summits.bed" | awk '{print $1}')"
PEAKS="$(wc -l < "${WORK}/macs3_frag.narrowPeak")"
RSS_KB="$(awk -F': ' '/Maximum resident set size/ {print $2}' "${WORK}/time.log")"
WALL="$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "${WORK}/time.log")"

echo "PASS: MACS3 FRAG output byte-identical to MACS3 3.0.3 (peaks=${PEAKS}, np_md5=${NP_MD5}, sm_md5=${SM_MD5}, RSS=${RSS_KB} kB, wall=${WALL})"
