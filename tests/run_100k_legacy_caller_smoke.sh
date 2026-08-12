#!/usr/bin/env bash
# Parity smoke for the LEGACY custom caller path (--out-prefix).
# This exercises the streaming legacy-caller code path (uint16 bins,
# sliding-window lambda, sparse BH, single-pass FragmentIterator drain).
#
# This is NOT the MACS3 FRAG path that libchromap calls in production.
# For that, see tests/run_100k_macs3_frag_smoke.sh.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCES[0]:-${BASH_SOURCE[0]}}")/.." && pwd)"
FIXTURE_FRAGS="${FIXTURE_FRAGS:-${DATA_DIR:-${REPO_ROOT}/data}/pbmc_unsorted_3k_100k/fragments.tsv}"
GOLDEN="${REPO_ROOT}/tests/golden/100k/golden.narrowPeak"
CLI="${REPO_ROOT}/bin/rapidmacs"

if [[ ! -x "${CLI}" ]]; then
  echo "ERROR: ${CLI} not built. Run 'make' first." >&2
  exit 2
fi
if [[ ! -s "${FIXTURE_FRAGS}" ]]; then
  echo "SKIP: 100K fragments fixture not present at ${FIXTURE_FRAGS}"
  echo "      This fixture is a Chromap realignment and is not redistributed."
  echo "      Point FIXTURE_FRAGS at a fragments TSV to run this test:"
  echo "        FIXTURE_FRAGS=/path/to/fragments.tsv make test-legacy-streaming-100k"
  exit 0
fi
if [[ ! -s "${GOLDEN}" ]]; then
  echo "ERROR: golden baseline missing at ${GOLDEN}" >&2
  exit 2
fi

WORK="$(mktemp -d /tmp/rapidmacs_legacy_smoke.XXXXXX)"
trap 'rm -rf "${WORK}"' EXIT

RSS_CEILING_KB="${RSS_CEILING_KB:-262144}"

echo "[smoke:legacy-streaming] running rapidmacs --out-prefix on 100K fixture ..."
/usr/bin/time -v -o "${WORK}/time.log" "${CLI}" \
  -i "${FIXTURE_FRAGS}" \
  --out-prefix "${WORK}/candidate" \
  --bin-width 50 --ext-size 150 --local-window 200 \
  --fdr 0.05 --p-value 0.01 --merge-gap 0 > "${WORK}/run.log" 2>&1

if ! diff -q "${WORK}/candidate.narrowPeak" "${GOLDEN}" >/dev/null; then
  echo "FAIL: legacy-caller narrowPeak diverged from golden" >&2
  diff "${WORK}/candidate.narrowPeak" "${GOLDEN}" | head -20 >&2
  exit 1
fi

CAND_MD5="$(md5sum "${WORK}/candidate.narrowPeak" | awk '{print $1}')"
PEAKS="$(wc -l < "${WORK}/candidate.narrowPeak")"
RSS_KB="$(awk -F': ' '/Maximum resident set size/ {print $2}' "${WORK}/time.log")"
WALL="$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "${WORK}/time.log")"

if [[ -z "${RSS_KB}" ]] || [[ "${RSS_KB}" -gt "${RSS_CEILING_KB}" ]]; then
  echo "FAIL: peak RSS=${RSS_KB} kB exceeds ceiling ${RSS_CEILING_KB} kB" >&2
  exit 1
fi

echo "PASS: legacy-caller narrowPeak byte-identical (md5=${CAND_MD5}, peaks=${PEAKS}, RSS=${RSS_KB} kB, wall=${WALL})"
