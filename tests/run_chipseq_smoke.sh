#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSAY=chipseq
PRESET=chip
TREATMENT="${CHIPSEQ_TREATMENT:-${REPO_ROOT}/tests/fixtures/chipseq_ENCFF686KKV.chr21_20m_30m.bam}"
CONTROL="${CHIPSEQ_CONTROL:-${REPO_ROOT}/tests/fixtures/chipseq_ENCFF181CXT.chr21_20m_30m.bam}"
GOLDEN_DIR="${CHIPSEQ_GOLDEN_DIR:-${REPO_ROOT}/tests/golden/chipseq_ENCFF686KKV}"
NAME=chip_chr21
source "${REPO_ROOT}/tests/run_paired_smoke_common.sh"
