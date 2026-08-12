#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSAY=cutrun
PRESET=cutrun
TREATMENT="${CUTRUN_TREATMENT:-${REPO_ROOT}/tests/fixtures/cutrun_SRR14255054.chr21_20m_30m.bam}"
CONTROL="${CUTRUN_CONTROL:-${REPO_ROOT}/tests/fixtures/cutrun_SRR14255053.chr21_20m_30m.bam}"
GOLDEN_DIR="${CUTRUN_GOLDEN_DIR:-${REPO_ROOT}/tests/golden/cutrun_SRR14255054}"
NAME=cutrun_chr21
source "${REPO_ROOT}/tests/run_paired_smoke_common.sh"

# Regression guard: requesting a diagnostic score track must not change peak
# or summit selection. This exercises the retained-ppois path, while the shared
# smoke above exercises chromosome-at-a-time peak collection.
"${CLI}" --preset "${PRESET}" -f BAMPE -t "${TREATMENT}" -c "${CONTROL}" \
  -g hs -q 0.05 -n "${NAME}" --out-prefix "${WORK}/candidate_with_ppois" \
  --frag-score-ppois-bdg "${WORK}/candidate.ppois.bdg" \
  >"${WORK}/ppois.stdout" 2>"${WORK}/ppois.stderr"
cmp "${WORK}/candidate_with_ppois.narrowPeak" \
  "${GOLDEN_DIR}/${NAME}_peaks.narrowPeak"
cmp "${WORK}/candidate_with_ppois.summits.bed" \
  "${GOLDEN_DIR}/${NAME}_summits.bed"
awk 'NF != 4 { bad=1 } END { exit bad }' "${WORK}/candidate.ppois.bdg"
echo "PASS: cutrun diagnostic ppois export is observational"
