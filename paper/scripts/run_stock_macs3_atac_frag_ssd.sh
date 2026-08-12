#!/usr/bin/env bash
# Measure released MACS3 3.0.3 on the supplied 10x PBMC and Chromap-realigned
# ENCSR095QNB bulk-ATAC fragment inputs. No RapidMACS process is run here.

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
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/stock_macs3_atac_frag}"
INPUT_ROOT="${RUN_ROOT}/inputs"
MACS3="${MACS3:-${MACS3_REF_ROOT}/macs3-3.0.3-venv/bin/macs3}"
MACS3_SOURCE="${MACS3_SOURCE:-${MACS3_REF_ROOT}/macs3-3.0.3-venv/lib/python3.10/site-packages/MACS3/Signal/PairedEndTrack.py}"
PBMC_SOURCE="${PBMC_SOURCE:-${ARCHIVE_ROOT}/pbmc_10x_supplied/pbmc_unsorted_3k_atac_fragments.tsv.gz}"
BULK_SOURCE="${BULK_SOURCE:-${ARCHIVE_ROOT}/encode_ENCSR095QNB/alignments/ENCSR095QNB.chromap.fragments.tsv.gz}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
MIN_FREE_BYTES="${MIN_FREE_BYTES:-8000000000}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)_stock_macs3_v3.0.3_atac_frag_ssd}"
OUT_ROOT="${OUT_ROOT:-${RUN_ROOT}/results/${RUN_ID}}"
ARCHIVE_PARENT="${ARCHIVE_PARENT:-${ARCHIVE_ROOT}/stock_macs3_atac_frag}"
ARCHIVE_DIR="${ARCHIVE_DIR:-${ARCHIVE_PARENT}/${RUN_ID}}"

PBMC_INPUT="${INPUT_ROOT}/10x/pbmc_unsorted_3k_atac_fragments.tsv.gz"
BULK_NATIVE_INPUT="${INPUT_ROOT}/bulk_ENCSR095QNB/ENCSR095QNB.chromap.fragments.7col.tsv.gz"
BULK_MACS3_INPUT="${INPUT_ROOT}/bulk_ENCSR095QNB/ENCSR095QNB.chromap.fragments.macs3.5col.tsv.gz"
PBMC_SOURCE_SHA256=0189df1de4894f10f48ffb905da4966cefc2466dc5fdd76ee2905244dd5a8dd7
BULK_SOURCE_SHA256=f77573b6360142143564d95fe342b869b7bad9bf752348129218d42ba44a427a
PBMC_EXPECTED_PEAKS=50003
PBMC_EXPECTED_NARROWPEAK_MD5=2c0de7c506514b0afef1e4339c386c9f
PBMC_EXPECTED_SUMMITS_MD5=c3805d1d3744ed78e3d589a523ad0b09

[[ "${REPEATS}" == 3 && "${WARMUPS}" == 1 ]] || {
  echo "ERROR: protocol requires WARMUPS=1 and REPEATS=3" >&2
  exit 2
}
[[ "$(readlink -f "${MACS3}")" == "${MACS3_REF_ROOT}/macs3-3.0.3-venv/bin/macs3" ]] || {
  echo "ERROR: refusing non-stock MACS3 executable: ${MACS3}" >&2
  exit 2
}
[[ -x "${MACS3}" && -s "${MACS3_SOURCE}" ]] || {
  echo "ERROR: stock MACS3 installation is incomplete" >&2
  exit 2
}
rg -q 'with open\("tmp_b_c\.txt", "w"\)' "${MACS3_SOURCE}" || {
  echo "ERROR: stock source no longer contains the released tmp_b_c.txt write" >&2
  exit 2
}
[[ "$("${MACS3}" --version 2>&1 | awk 'NR == 1 {print $NF}')" == 3.0.3 ]] || {
  echo "ERROR: MACS3 is not version 3.0.3" >&2
  exit 2
}
for tool in rg gzip bgzip awk sha256sum md5sum findmnt lsblk /usr/bin/time; do
  command -v "${tool}" >/dev/null 2>&1 || { echo "ERROR: missing ${tool}" >&2; exit 2; }
done
grep -q "i9-13900KF" /proc/cpuinfo || guard 1 "host CPU is not the i9-13900KF"
[[ -s "${PBMC_SOURCE}" && -s "${BULK_SOURCE}" ]] || { echo "ERROR: archived input is missing" >&2; exit 2; }
[[ "$(sha256sum "${PBMC_SOURCE}" | awk '{print $1}')" == "${PBMC_SOURCE_SHA256}" ]] || {
  echo "ERROR: supplied 10x source checksum mismatch" >&2; exit 2;
}
[[ "$(sha256sum "${BULK_SOURCE}" | awk '{print $1}')" == "${BULK_SOURCE_SHA256}" ]] || {
  echo "ERROR: ENCSR095QNB fragment source checksum mismatch" >&2; exit 2;
}
[[ ! -e "${OUT_ROOT}" ]] || { echo "ERROR: output already exists: ${OUT_ROOT}" >&2; exit 2; }
[[ ! -e "${ARCHIVE_DIR}" ]] || { echo "ERROR: archive already exists: ${ARCHIVE_DIR}" >&2; exit 2; }

storage_source="$(findmnt -n -o SOURCE -T "$(dirname "${RUN_ROOT}")")"
storage_rotational="$(lsblk -ndo ROTA "${storage_source}" 2>/dev/null | awk 'NR == 1 {print $1}')"
[[ "${storage_rotational}" == 0 ]] || guard 1 "${RUN_ROOT} is not on SSD/NVMe"
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another peak-calling process is active" >&2
  exit 2
fi

mkdir -p "$(dirname "${PBMC_INPUT}")" "$(dirname "${BULK_NATIVE_INPUT}")" \
  "${OUT_ROOT}/runs" "${OUT_ROOT}/provenance" "${ARCHIVE_PARENT}"

stage_verified() {
  local source="$1" destination="$2" expected="$3"
  if [[ -s "${destination}" && "$(sha256sum "${destination}" | awk '{print $1}')" == "${expected}" ]]; then
    echo "[reuse] ${destination}"
    return
  fi
  [[ ! -e "${destination}" ]] || { echo "ERROR: wrong staged file: ${destination}" >&2; exit 2; }
  cp --reflink=auto --preserve=timestamps "${source}" "${destination}.partial"
  [[ "$(sha256sum "${destination}.partial" | awk '{print $1}')" == "${expected}" ]] || {
    echo "ERROR: staged checksum mismatch: ${destination}.partial" >&2; exit 2;
  }
  mv "${destination}.partial" "${destination}"
}

stage_verified "${PBMC_SOURCE}" "${PBMC_INPUT}" "${PBMC_SOURCE_SHA256}"
stage_verified "${BULK_SOURCE}" "${BULK_NATIVE_INPUT}" "${BULK_SOURCE_SHA256}"

if [[ ! -s "${BULK_MACS3_INPUT}" ]]; then
  [[ ! -e "${BULK_MACS3_INPUT}" ]] || { echo "ERROR: empty bulk projection exists" >&2; exit 2; }
  gzip -cd "${BULK_NATIVE_INPUT}" |
    awk -F '\t' -v OFS='\t' '
      NF != 7 {printf "ERROR: source line %d has %d fields\n", NR, NF > "/dev/stderr"; exit 1}
      $4 != "N" || $5 !~ /^[0-9]+$/ || ($6 != "+" && $6 != "-") || $7 !~ /^[0-9]+$/ {
        printf "ERROR: invalid seven-column fragment at line %d\n", NR > "/dev/stderr"; exit 1
      }
      {print $1,$2,$3,"bulk",$7}
    ' | bgzip -@ 8 -c >"${BULK_MACS3_INPUT}.partial"
  gzip -t "${BULK_MACS3_INPUT}.partial"
  mv "${BULK_MACS3_INPUT}.partial" "${BULK_MACS3_INPUT}"
fi

gzip -t "${PBMC_INPUT}" "${BULK_NATIVE_INPUT}" "${BULK_MACS3_INPUT}"
gzip -cd "${PBMC_INPUT}" | awk -F '\t' '
  /^#/ || NF == 0 {next}
  NF != 5 {printf "ERROR: 10x line %d has %d fields\n", NR, NF > "/dev/stderr"; exit 1}
  {rows++}
  END {if (!rows) exit 1; printf "rows\t%d\ncolumns\t5\n", rows}
' >"${OUT_ROOT}/provenance/pbmc_schema.tsv"
gzip -cd "${BULK_NATIVE_INPUT}" | awk -F '\t' '
  NF != 7 {exit 1}
  {rows++; bp += $3-$2; if ($7 > max_count) max_count=$7}
  END {
    if (!rows) exit 1
    printf "rows\t%d\ncolumns\t7\nmean_bp\t%d\nmax_count\t%d\n", rows, bp/rows, max_count
  }
' >"${OUT_ROOT}/provenance/bulk_native_schema.tsv"
gzip -cd "${BULK_MACS3_INPUT}" | awk -F '\t' '
  NF != 5 || $4 != "bulk" || $5 !~ /^[0-9]+$/ {exit 1}
  {rows++; bp += $3-$2}
  END {
    if (!rows) exit 1
    printf "rows\t%d\ncolumns\t5\nmean_bp\t%d\n", rows, bp/rows
  }
' >"${OUT_ROOT}/provenance/bulk_macs3_schema.tsv"

bulk_geometry="$(awk -F '\t' '$1 == "mean_bp" {print $2}' "${OUT_ROOT}/provenance/bulk_native_schema.tsv")"
[[ "${bulk_geometry}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid bulk geometry" >&2; exit 2; }

for input in "${PBMC_INPUT}" "${BULK_NATIVE_INPUT}" "${BULK_MACS3_INPUT}"; do
  source_device="$(findmnt -n -o SOURCE -T "${input}")"
  rotational="$(lsblk -ndo ROTA "${source_device}" 2>/dev/null | awk 'NR == 1 {print $1}')"
  [[ "${rotational}" == 0 ]] || guard 1 "staged input is not on NVMe: ${input}"
done

RAW_TSV="${OUT_ROOT}/raw.tsv"
WARMUP_TSV="${OUT_ROOT}/warmup.tsv"
COMMANDS="${OUT_ROOT}/provenance/commands.log"
printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"
printf 'assay\trepeat\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tsummits_md5\n' >"${WARMUP_TSV}"

{
  printf 'field\tvalue\n'
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'out_root\t%s\n' "${OUT_ROOT}"
  printf 'archive_dir\t%s\n' "${ARCHIVE_DIR}"
  printf 'assays\t10x_pbmc_frag bulk_atac_ENCSR095QNB_frag\n'
  printf 'tools\tmacs3\n'
  printf 'threads\t1\n'
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'macs3_executable\t%s\n' "${MACS3}"
  printf 'macs3_version\t%s\n' "$("${MACS3}" --version 2>&1 | head -n 1)"
  printf 'macs3_source\t%s\n' "${MACS3_SOURCE}"
  printf 'tmp_b_c_source_line\t%s\n' "$(rg -n 'with open\("tmp_b_c\.txt", "w"\)' "${MACS3_SOURCE}" | head -n 1)"
  printf 'bulk_geometry_bp\t%s\n' "${bulk_geometry}"
  printf 'storage_source\t%s\n' "${storage_source}"
  printf 'storage_rotational\t%s\n' "${storage_rotational}"
} >"${OUT_ROOT}/provenance/metadata.tsv"

sha256sum "${PBMC_SOURCE}" "${BULK_SOURCE}" >"${OUT_ROOT}/provenance/source_sha256s.txt"
sha256sum "${PBMC_INPUT}" "${BULK_NATIVE_INPUT}" "${BULK_MACS3_INPUT}" \
  >"${OUT_ROOT}/provenance/input_sha256s.txt"
sha256sum "${MACS3}" "${MACS3_SOURCE}" \
  "$(dirname "${MACS3_SOURCE}")/PairedEndTrack.cpython-310-x86_64-linux-gnu.so" \
  >"${OUT_ROOT}/provenance/macs3_sha256s.txt"
git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
cp "${BASH_SOURCE[0]}" "${REPO_ROOT}/paper/scripts/run_encode_bulk_atac_macs3_baseline.sh" \
  "${REPO_ROOT}/paper/scripts/run_treatment_only_frag_benchmarks.sh" \
  "${OUT_ROOT}/provenance/"
{
  uname -a
  lscpu
  free -h
  findmnt -T "${OUT_ROOT}"
  lsblk -o NAME,TYPE,SIZE,ROTA,FSTYPE,MOUNTPOINTS
  /usr/bin/time --version | head -n 1
} >"${OUT_ROOT}/provenance/environment.txt"

log_command() {
  printf '%s' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >>"${COMMANDS}"
  printf '\tcd\t%q\t&&' "$1" >>"${COMMANDS}"
  shift
  printf '\t%q' "$@" >>"${COMMANDS}"
  printf '\n' >>"${COMMANDS}"
}

bed3_md5() {
  awk -F '\t' -v OFS='\t' '!/^track/ && NF >= 3 {print $1,$2,$3}' "$1" |
    md5sum | awk '{print $1}'
}

free_bytes() {
  df -B1 --output=avail "${OUT_ROOT}" | awk 'NR == 2 {print $1}'
}

run_one() {
  local assay="$1" repeat="$2" input="$3" name="$4" expected_peaks="$5"
  local expected_np="$6" expected_sm="$7"
  shift 7
  local dir="${OUT_ROOT}/runs/${assay}.macs3.r${repeat}"
  local np sm wall rss peaks np_md5 np_bed3 sm_md5 sm_bed3
  mkdir -p "${dir}/out" "${dir}/tmp"
  if (( $(free_bytes) < MIN_FREE_BYTES )); then
    echo "ERROR: insufficient free SSD space before ${assay} ${repeat}" >&2
    exit 2
  fi
  args=("${MACS3}" callpeak -t "${input}" -f FRAG -g hs -n "${name}" --outdir "${dir}/out" "$@")
  log_command "${dir}" env "TMPDIR=${dir}/tmp" /usr/bin/time -f '%e\t%M' \
    -o "${dir}/time.tsv" "${args[@]}"
  (
    cd "${dir}"
    TMPDIR="${dir}/tmp" /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" \
      "${args[@]}" >"${dir}/stdout" 2>"${dir}/stderr"
  )
  if [[ -e "${dir}/tmp_b_c.txt" ]]; then
    echo "ERROR: treatment-only FRAG unexpectedly reached the PE control debug path" >&2
    exit 1
  fi
  read -r wall rss <"${dir}/time.tsv"
  np="${dir}/out/${name}_peaks.narrowPeak"
  sm="${dir}/out/${name}_summits.bed"
  peaks="$(wc -l <"${np}")"
  np_md5="$(md5sum "${np}" | awk '{print $1}')"
  np_bed3="$(bed3_md5 "${np}")"
  sm_md5="$(md5sum "${sm}" | awk '{print $1}')"
  sm_bed3="$(bed3_md5 "${sm}")"
  if [[ -n "${expected_peaks}" ]] &&
     [[ "${peaks}" != "${expected_peaks}" || "${np_md5}" != "${expected_np}" || "${sm_md5}" != "${expected_sm}" ]]; then
    echo "ERROR: ${assay} ${repeat} differs from its previous MACS3 baseline" >&2
    echo "observed peaks=${peaks} narrowPeak=${np_md5} summits=${sm_md5}" >&2
    exit 1
  fi
  if [[ "${repeat}" == warmup* ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${assay}" "${repeat}" "${wall}" "${rss}" "${peaks}" "${np_md5}" "${sm_md5}" \
      >>"${WARMUP_TSV}"
  else
    printf '%s\tMACS3\t1\t%s\twarm\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${assay}" "${repeat}" "${wall}" "${rss}" "${peaks}" "${np_md5}" \
      "${np_bed3}" "${sm_md5}" "${sm_bed3}" >>"${RAW_TSV}"
  fi
  echo "[complete] ${assay} ${repeat}: ${wall}s, ${rss} KB RSS, ${peaks} peaks"
}

benchmark_assay() {
  local assay="$1" input="$2" name="$3" expected_peaks="$4" expected_np="$5" expected_sm="$6"
  shift 6
  local warm repeat
  for ((warm=1; warm<=WARMUPS; ++warm)); do
    echo "[warmup] ${assay} ${warm}/${WARMUPS}"
    run_one "${assay}" "warmup${warm}" "${input}" "${name}" \
      "${expected_peaks}" "${expected_np}" "${expected_sm}" "$@"
  done
  for ((repeat=1; repeat<=REPEATS; ++repeat)); do
    echo "[benchmark] ${assay} ${repeat}/${REPEATS}"
    run_one "${assay}" "${repeat}" "${input}" "${name}" \
      "${expected_peaks}" "${expected_np}" "${expected_sm}" "$@"
  done
  if ! awk -F '\t' -v a="${assay}" '
      $1 == a && !seen {peaks=$8; np=$9; sm=$11; seen=1}
      $1 == a && ($8 != peaks || $9 != np || $11 != sm) {bad=1}
      END {exit bad}
    ' "${RAW_TSV}"; then
    echo "ERROR: ${assay} output changed across recorded repetitions" >&2
    exit 1
  fi
}

# Run the smaller bulk sample first as a complete stock-FRAG protocol check,
# then the supplied 10x PBMC dataset. Assays are never run concurrently.
benchmark_assay bulk_atac_ENCSR095QNB_frag "${BULK_MACS3_INPUT}" ENCSR095QNB_chromap_realign \
  '' '' '' -q 0.05 --max-count 1 --min-length "${bulk_geometry}" --max-gap "${bulk_geometry}"
benchmark_assay 10x_pbmc_frag "${PBMC_INPUT}" atac_treatment_only \
  "${PBMC_EXPECTED_PEAKS}" "${PBMC_EXPECTED_NARROWPEAK_MD5}" "${PBMC_EXPECTED_SUMMITS_MD5}" \
  -p 1e-5 --min-length 200 --max-gap 30

python3 "${REPO_ROOT}/paper/scripts/summarize_paper_benchmarks.py" \
  "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
printf 'completed_utc\t%s\nfinal_free_bytes\t%s\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(free_bytes)" >>"${OUT_ROOT}/provenance/metadata.tsv"

(
  cd "${OUT_ROOT}"
  find . -type f ! -path './provenance/archive_manifest.sha256' -print0 |
    sort -z | xargs -0 sha256sum
) >"${OUT_ROOT}/provenance/archive_manifest.sha256"
cp -a "${OUT_ROOT}" "${ARCHIVE_DIR}"
(
  cd "${ARCHIVE_DIR}"
  sha256sum -c provenance/archive_manifest.sha256 >/dev/null
)

echo "[complete] summary: ${OUT_ROOT}/summary.tsv"
echo "[complete] archive: ${ARCHIVE_DIR}"
