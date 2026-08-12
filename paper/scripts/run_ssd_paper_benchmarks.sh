#!/usr/bin/env bash
# Run the paper ATAC/ChIP/CUT benchmark matrix sequentially from SSD inputs.

set -euo pipefail

# Host and storage guards. The paper's protocol requires the reference host and
# NVMe placement; those checks are advisory here so the harness runs anywhere.
# Set STRICT_HOST=1 to restore them as hard failures.
STRICT_HOST="${STRICT_HOST:-0}"
guard() {  # guard <ok:0|1> <message>
  [[ "$1" == 0 ]] || return 0
  if [[ "${STRICT_HOST}" == 1 ]]; then echo "ERROR: $2" >&2; exit 2; fi
  echo "WARNING: $2 (timings will not be comparable to the paper)" >&2
}


# ---------------------------------------------------------------------------
# Paths are environment-driven so this runs anywhere. Defaults are relative to
# the repository, so nothing outside it is touched unless you point these
# elsewhere. See paper/README.md.
#   WORK_ROOT      scratch space for run outputs   (default ./work)
#   DATA_DIR       input datasets                  (default ./data)
#   ARCHIVE_ROOT   where completed runs are copied (default ./archive)
#   MACS3_REF_ROOT MACS3 reference installations   (default ./.macs3-reference)
# ---------------------------------------------------------------------------
WORK_ROOT="${WORK_ROOT:-$(pwd)/work}"
DATA_DIR="${DATA_DIR:-$(pwd)/data}"
ARCHIVE_ROOT="${ARCHIVE_ROOT:-$(pwd)/archive}"
MACS3_REF_ROOT="${MACS3_REF_ROOT:-$(pwd)/.macs3-reference}"
CHROMAP_SUITE="${CHROMAP_SUITE:-$(pwd)/chromap-suite}"
CELLRANGER_REF="${CELLRANGER_REF:-${DATA_DIR}/refdata-cellranger-arc-GRCh38-2020-A-2.0.0}"


REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/paper_benchmarks}"
INPUT_ROOT="${RUN_ROOT}/inputs"
ATAC_TREATMENT="${ATAC_TREATMENT:-${WORK_ROOT}/fastq_realign/alignments/atac/fragments.tsv.gz}"
CHIP_TREATMENT="${INPUT_ROOT}/chipseq/ENCFF686KKV.bam"
CHIP_CONTROL="${INPUT_ROOT}/chipseq/ENCFF181CXT.bam"
CUTRUN_TREATMENT="${INPUT_ROOT}/cutrun/SRR14255054.bam"
CUTRUN_CONTROL="${INPUT_ROOT}/cutrun/SRR14255053.bam"
THREADS="${THREADS:-1 4 24}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
ASSAYS="${ASSAYS:-atac cutrun chipseq}"
TOOLS="${TOOLS:-macs3 rapidmacs}"
OUT_ROOT="${OUT_ROOT:-${RUN_ROOT}/results/$(date -u +%Y%m%dT%H%M%SZ)_ssd_paper_matrix}"

mkdir -p "${OUT_ROOT}/provenance"

inputs=("${ATAC_TREATMENT}" "${CHIP_TREATMENT}" "${CHIP_CONTROL}" "${CUTRUN_TREATMENT}" "${CUTRUN_CONTROL}")
for input in "${inputs[@]}"; do
  [[ -s "${input}" ]] || { echo "ERROR: missing ${input}" >&2; exit 2; }
  source_device="$(findmnt -n -o SOURCE -T "${input}")"
  rotational="$(lsblk -ndo ROTA "${source_device}" 2>/dev/null | awk 'NR == 1 {print $1}' || true)"
  [[ "${rotational}" == 0 ]] || guard 1 "input is not on NVMe/SSD: ${input}"
done

if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another peak-calling run is active" >&2
  exit 2
fi

{
  printf 'field\tvalue\n'
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_root\t%s\n' "${RUN_ROOT}"
  printf 'out_root\t%s\n' "${OUT_ROOT}"
  printf 'threads\t%s\n' "${THREADS}"
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'assays\t%s\n' "${ASSAYS}"
  printf 'tools\t%s\n' "${TOOLS}"
} >"${OUT_ROOT}/provenance/metadata.tsv"
sha256sum "${inputs[@]}" >"${OUT_ROOT}/provenance/input_sha256s.txt"
git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
cp "${REPO_ROOT}/paper/scripts/run_paper_benchmarks.sh" \
  "${REPO_ROOT}/paper/scripts/stage_paper_bampe_to_ssd.sh" \
  "${REPO_ROOT}/paper/scripts/run_ssd_paper_benchmarks.sh" \
  "${OUT_ROOT}/provenance/"

ATAC_TREATMENT="${ATAC_TREATMENT}" \
CHIP_TREATMENT="${CHIP_TREATMENT}" CHIP_CONTROL="${CHIP_CONTROL}" \
CUTRUN_TREATMENT="${CUTRUN_TREATMENT}" CUTRUN_CONTROL="${CUTRUN_CONTROL}" \
THREADS="${THREADS}" REPEATS="${REPEATS}" WARMUPS="${WARMUPS}" \
ASSAYS="${ASSAYS}" TOOLS="${TOOLS}" OUT_ROOT="${OUT_ROOT}" \
  "${REPO_ROOT}/paper/scripts/run_paper_benchmarks.sh" 2>&1 | tee "${OUT_ROOT}/benchmark.log"

printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${OUT_ROOT}/provenance/metadata.tsv"
echo "[complete] SSD benchmark summary: ${OUT_ROOT}/summary.tsv"
