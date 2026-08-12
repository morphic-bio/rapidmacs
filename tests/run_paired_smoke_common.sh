#!/usr/bin/env bash

# Shared BAMPE treatment/control smoke. The caller must define ASSAY, PRESET,
# TREATMENT, CONTROL, GOLDEN_DIR, and NAME before sourcing this file.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${REPO_ROOT}/bin/rapidmacs"

for var in ASSAY PRESET TREATMENT CONTROL GOLDEN_DIR NAME; do
  if [[ -z "${!var:-}" ]]; then
    echo "ERROR: paired smoke variable ${var} is unset" >&2
    exit 2
  fi
done
if [[ ! -x "${CLI}" ]]; then
  echo "ERROR: ${CLI} not built. Run 'make' first." >&2
  exit 2
fi
for path in "${TREATMENT}" "${CONTROL}" \
            "${GOLDEN_DIR}/${NAME}_peaks.narrowPeak" \
            "${GOLDEN_DIR}/${NAME}_summits.bed"; do
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: required ${ASSAY} fixture/golden is missing: ${path}" >&2
    exit 2
  fi
done

WORK="$(mktemp -d "/tmp/rapidmacs_${ASSAY}_smoke.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

echo "[smoke:${ASSAY}] BAMPE treatment/control parity ..."
/usr/bin/time -v -o "${WORK}/time.log" "${CLI}" \
  --preset "${PRESET}" -f BAMPE -t "${TREATMENT}" -c "${CONTROL}" \
  -g hs -q 0.05 -n "${NAME}" --out-prefix "${WORK}/candidate" \
  >"${WORK}/run.stdout" 2>"${WORK}/run.stderr"

candidate_np="${WORK}/candidate.narrowPeak"
candidate_sm="${WORK}/candidate.summits.bed"
golden_np="${GOLDEN_DIR}/${NAME}_peaks.narrowPeak"
golden_sm="${GOLDEN_DIR}/${NAME}_summits.bed"

# The paper claims full-file byte identity, so every narrowPeak field is part
# of the regression contract (coordinates, names, scores, and summit offset).
if ! cmp -s "${candidate_np}" "${golden_np}"; then
  echo "FAIL: ${ASSAY} narrowPeak differs from MACS3 3.0.3" >&2
  diff -u "${golden_np}" "${candidate_np}" | head -40 >&2
  exit 1
fi

if ! cmp -s "${candidate_sm}" "${golden_sm}"; then
  echo "FAIL: ${ASSAY} summits differ from MACS3 3.0.3" >&2
  diff -u "${golden_sm}" "${candidate_sm}" | head -40 >&2
  exit 1
fi

peaks="$(wc -l <"${candidate_np}")"
rss_kb="$(awk -F': ' '/Maximum resident set size/ {print $2}' "${WORK}/time.log")"
wall="$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' "${WORK}/time.log")"
np_md5="$(md5sum "${candidate_np}" | awk '{print $1}')"
sm_md5="$(md5sum "${candidate_sm}" | awk '{print $1}')"
echo "PASS: ${ASSAY} parity (peaks=${peaks}, np_md5=${np_md5}, sm_md5=${sm_md5}, RSS=${rss_kb} kB, wall=${wall})"
