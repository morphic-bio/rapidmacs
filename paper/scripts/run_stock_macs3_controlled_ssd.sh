#!/usr/bin/env bash
# Re-measure only the controlled ChIP-seq and CUT&RUN MACS3 baselines with
# released MACS3 3.0.3, including its shipped tmp_b_c.txt diagnostic write.

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
RUN_ROOT="${RUN_ROOT:-${WORK_ROOT}/stock_macs3_controlled}"
INPUT_ROOT="${INPUT_ROOT:-${DATA_DIR}}"
MACS3="${MACS3:-${MACS3_REF_ROOT}/macs3-3.0.3-venv/bin/macs3}"
MACS3_SOURCE="${MACS3_SOURCE:-${MACS3_REF_ROOT}/macs3-3.0.3-venv/lib/python3.10/site-packages/MACS3/Signal/PairedEndTrack.py}"
REPEATS="${REPEATS:-3}"
WARMUPS="${WARMUPS:-1}"
MIN_FREE_BYTES="${MIN_FREE_BYTES:-8000000000}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)_stock_macs3_v3.0.3_controlled_ssd}"
OUT_ROOT="${OUT_ROOT:-${RUN_ROOT}/results/${RUN_ID}}"
ARCHIVE_PARENT="${ARCHIVE_PARENT:-${ARCHIVE_ROOT}/stock_macs3_controlled}"
ARCHIVE_DIR="${ARCHIVE_DIR:-${ARCHIVE_PARENT}/${RUN_ID}}"

CHIP_TREATMENT="${INPUT_ROOT}/chipseq/ENCFF686KKV.bam"
CHIP_CONTROL="${INPUT_ROOT}/chipseq/ENCFF181CXT.bam"
CUTRUN_TREATMENT="${INPUT_ROOT}/cutrun/SRR14255054.bam"
CUTRUN_CONTROL="${INPUT_ROOT}/cutrun/SRR14255053.bam"

EXPECTED_CHIP_TREATMENT_SHA256=37ee55258385f908e139e934dcb542360e91ce0dc4dca4867cc35f5507fbe4ea
EXPECTED_CHIP_CONTROL_SHA256=400347e0ebf66079f47bb06593129d0237bf1ceac7579dd9442ace0d2a57ad1f
EXPECTED_CUTRUN_TREATMENT_SHA256=28f4f45a31b4a90da5a68c44b093b166f0c8bde009b04e00d549636d5fcae95a
EXPECTED_CUTRUN_CONTROL_SHA256=33f74032ce65e141bb137582201162f4dddee3950e1983e8963097dc9f57149d
EXPECTED_CHIP_PEAKS=33985
EXPECTED_CHIP_NARROWPEAK_MD5=5463ccacf4af902c3ad6f388bd141fc4
EXPECTED_CHIP_SUMMITS_MD5=82ae9573bb99733e30635cef96da1f69
EXPECTED_CUTRUN_PEAKS=9097
EXPECTED_CUTRUN_NARROWPEAK_MD5=0fc36102443eae226578201f18df67c3
EXPECTED_CUTRUN_SUMMITS_MD5=1cd24d0c399a5ea0b131f8c6e5fe9c7e

[[ "${REPEATS}" == 3 && "${WARMUPS}" == 1 ]] || {
  echo "ERROR: the controlled protocol requires WARMUPS=1 and REPEATS=3" >&2
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
  echo "ERROR: stock PairedEndTrack.py does not contain the tmp_b_c.txt write" >&2
  exit 2
}
[[ "$("${MACS3}" --version 2>&1 | awk 'NR == 1 {print $NF}')" == 3.0.3 ]] || {
  echo "ERROR: stock MACS3 is not version 3.0.3" >&2
  exit 2
}
command -v /usr/bin/time >/dev/null 2>&1 || { echo "ERROR: GNU time is missing" >&2; exit 2; }
command -v rg >/dev/null 2>&1 || { echo "ERROR: ripgrep is missing" >&2; exit 2; }
grep -q "i9-13900KF" /proc/cpuinfo || guard 1 "host CPU is not the i9-13900KF"
[[ ! -e "${OUT_ROOT}" ]] || { echo "ERROR: output already exists: ${OUT_ROOT}" >&2; exit 2; }
[[ ! -e "${ARCHIVE_DIR}" ]] || { echo "ERROR: archive already exists: ${ARCHIVE_DIR}" >&2; exit 2; }

inputs=("${CHIP_TREATMENT}" "${CHIP_CONTROL}" "${CUTRUN_TREATMENT}" "${CUTRUN_CONTROL}")
expected_input_hashes=(
  "${EXPECTED_CHIP_TREATMENT_SHA256}"
  "${EXPECTED_CHIP_CONTROL_SHA256}"
  "${EXPECTED_CUTRUN_TREATMENT_SHA256}"
  "${EXPECTED_CUTRUN_CONTROL_SHA256}"
)
for i in "${!inputs[@]}"; do
  input="${inputs[$i]}"
  [[ -s "${input}" ]] || { echo "ERROR: missing input: ${input}" >&2; exit 2; }
  source_device="$(findmnt -n -o SOURCE -T "${input}")"
  rotational="$(lsblk -ndo ROTA "${source_device}" 2>/dev/null | awk 'NR == 1 {print $1}')"
  [[ "${rotational}" == 0 ]] || guard 1 "input is not on NVMe/SSD: ${input}"
  observed_hash="$(sha256sum "${input}" | awk '{print $1}')"
  [[ "${observed_hash}" == "${expected_input_hashes[$i]}" ]] || {
    echo "ERROR: input checksum mismatch: ${input}" >&2
    exit 2
  }
done

out_source="$(findmnt -n -o SOURCE -T "$(dirname "${RUN_ROOT}")")"
out_rotational="$(lsblk -ndo ROTA "${out_source}" 2>/dev/null | awk 'NR == 1 {print $1}')"
[[ "${out_rotational}" == 0 ]] || guard 1 "output root is not on NVMe/SSD"
if ps -eo args | awk '/[m]acs3 callpeak|[/]bin\/rapidmacs / {found=1} END {exit found ? 0 : 1}'; then
  echo "ERROR: another peak-calling process is active" >&2
  exit 2
fi

mkdir -p "${OUT_ROOT}/runs" "${OUT_ROOT}/provenance" "${ARCHIVE_PARENT}"
RAW_TSV="${OUT_ROOT}/raw.tsv"
WARMUP_TSV="${OUT_ROOT}/warmup.tsv"
DUMP_TSV="${OUT_ROOT}/tmp_b_c.tsv"
COMMANDS="${OUT_ROOT}/provenance/commands.log"
printf 'assay\ttool\tthreads\trepeat\tcache_state\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tbed3_md5\tsummits_md5\tsummits_bed3_md5\n' >"${RAW_TSV}"
printf 'assay\trepeat\twall_s\trss_kb\tpeaks\tnarrowpeak_md5\tsummits_md5\n' >"${WARMUP_TSV}"
printf 'assay\trepeat\tlogical_bytes\tallocated_bytes\tfree_bytes_before\tfree_bytes_after\tremoved_after_measurement\n' >"${DUMP_TSV}"

{
  printf 'field\tvalue\n'
  printf 'started_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'run_root\t%s\n' "${RUN_ROOT}"
  printf 'out_root\t%s\n' "${OUT_ROOT}"
  printf 'archive_dir\t%s\n' "${ARCHIVE_DIR}"
  printf 'assays\tcutrun chipseq\n'
  printf 'tools\tmacs3\n'
  printf 'threads\t1\n'
  printf 'repeats\t%s\n' "${REPEATS}"
  printf 'warmups\t%s\n' "${WARMUPS}"
  printf 'macs3_executable\t%s\n' "${MACS3}"
  printf 'macs3_version\t%s\n' "$("${MACS3}" --version 2>&1 | head -n 1)"
  printf 'macs3_source\t%s\n' "${MACS3_SOURCE}"
  printf 'tmp_b_c_source_line\t%s\n' "$(rg -n 'with open\("tmp_b_c\.txt", "w"\)' "${MACS3_SOURCE}" | head -n 1)"
  printf 'storage_source\t%s\n' "${out_source}"
  printf 'storage_rotational\t%s\n' "${out_rotational}"
  printf 'minimum_free_bytes\t%s\n' "${MIN_FREE_BYTES}"
} >"${OUT_ROOT}/provenance/metadata.tsv"

sha256sum "${inputs[@]}" >"${OUT_ROOT}/provenance/input_sha256s.txt"
sha256sum "${MACS3}" "${MACS3_SOURCE}" \
  "$(dirname "${MACS3_SOURCE}")/PairedEndTrack.cpython-310-x86_64-linux-gnu.so" \
  >"${OUT_ROOT}/provenance/macs3_sha256s.txt"
git -C "${REPO_ROOT}" status --short >"${OUT_ROOT}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary HEAD >"${OUT_ROOT}/provenance/git_worktree.patch"
cp "${BASH_SOURCE[0]}" "${REPO_ROOT}/paper/scripts/run_paper_benchmarks.sh" \
  "${REPO_ROOT}/paper/scripts/run_ssd_paper_benchmarks.sh" \
  "${REPO_ROOT}/paper/scripts/macs3-3.0.3-no-control-dump.patch" \
  "${OUT_ROOT}/provenance/"
{
  uname -a
  lscpu
  free -h
  findmnt -T "${OUT_ROOT}"
  lsblk -o NAME,TYPE,SIZE,ROTA,FSTYPE,MOUNTPOINTS
  /usr/bin/time --version | head -n 1
} >"${OUT_ROOT}/provenance/environment.txt"
{
  printf 'path\tdevice\trotational\n'
  for input in "${inputs[@]}" "${OUT_ROOT}"; do
    device="$(findmnt -n -o SOURCE -T "${input}")"
    rotational="$(lsblk -ndo ROTA "${device}" 2>/dev/null | awk 'NR == 1 {print $1}')"
    printf '%s\t%s\t%s\n' "${input}" "${device}" "${rotational}"
  done
} >"${OUT_ROOT}/provenance/input_output_storage.tsv"

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
  local assay="$1" repeat="$2" treatment="$3" control="$4" name="$5"
  local expected_peaks="$6" expected_np="$7" expected_sm="$8"
  local dir="${OUT_ROOT}/runs/${assay}.macs3.r${repeat}"
  local free_before free_after dump dump_size dump_allocated np sm
  local wall rss peaks np_md5 np_bed3 sm_md5 sm_bed3
  mkdir -p "${dir}/out" "${dir}/tmp"
  free_before="$(free_bytes)"
  if (( free_before < MIN_FREE_BYTES )); then
    echo "ERROR: only ${free_before} bytes free before ${assay} ${repeat}" >&2
    exit 2
  fi

  args=("${MACS3}" callpeak -t "${treatment}" -f BAMPE -g hs -n "${name}"
        --outdir "${dir}/out" -c "${control}" -q 0.05)
  log_command "${dir}" env "TMPDIR=${dir}/tmp" /usr/bin/time -f '%e\t%M' \
    -o "${dir}/time.tsv" "${args[@]}"
  (
    cd "${dir}"
    TMPDIR="${dir}/tmp" /usr/bin/time -f '%e\t%M' -o "${dir}/time.tsv" \
      "${args[@]}" >"${dir}/stdout" 2>"${dir}/stderr"
  )

  dump="${dir}/tmp_b_c.txt"
  [[ -s "${dump}" ]] || {
    echo "ERROR: stock MACS3 did not produce ${dump}" >&2
    exit 1
  }
  dump_size="$(stat -c '%s' "${dump}")"
  dump_allocated="$(du -B1 "${dump}" | awk '{print $1}')"
  free_after="$(free_bytes)"
  printf '%s\t%s\t%s\t%s\t%s\t%s\tyes\n' \
    "${assay}" "${repeat}" "${dump_size}" "${dump_allocated}" \
    "${free_before}" "${free_after}" >>"${DUMP_TSV}"

  read -r wall rss <"${dir}/time.tsv"
  np="${dir}/out/${name}_peaks.narrowPeak"
  sm="${dir}/out/${name}_summits.bed"
  [[ -s "${np}" && -s "${sm}" ]] || { echo "ERROR: missing MACS3 output" >&2; exit 1; }
  peaks="$(wc -l <"${np}")"
  np_md5="$(md5sum "${np}" | awk '{print $1}')"
  np_bed3="$(bed3_md5 "${np}")"
  sm_md5="$(md5sum "${sm}" | awk '{print $1}')"
  sm_bed3="$(bed3_md5 "${sm}")"
  if [[ "${peaks}" != "${expected_peaks}" || "${np_md5}" != "${expected_np}" ||
        "${sm_md5}" != "${expected_sm}" ]]; then
    echo "ERROR: ${assay} ${repeat} output differs from the existing matrix" >&2
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

  # The released diagnostic file is never consumed. Record its disk footprint,
  # then discard this per-run temporary before starting the next isolated run.
  rm -f -- "${dump}"
  echo "[complete] ${assay} ${repeat}: ${wall}s, ${rss} KB RSS, tmp_b_c=${dump_size} bytes"
}

benchmark_assay() {
  local assay="$1" treatment="$2" control="$3" name="$4"
  local expected_peaks="$5" expected_np="$6" expected_sm="$7"
  local warm repeat
  for ((warm=1; warm<=WARMUPS; ++warm)); do
    echo "[warmup] ${assay} ${warm}/${WARMUPS}"
    run_one "${assay}" "warmup${warm}" "${treatment}" "${control}" "${name}" \
      "${expected_peaks}" "${expected_np}" "${expected_sm}"
  done
  for ((repeat=1; repeat<=REPEATS; ++repeat)); do
    echo "[benchmark] ${assay} ${repeat}/${REPEATS}"
    run_one "${assay}" "${repeat}" "${treatment}" "${control}" "${name}" \
      "${expected_peaks}" "${expected_np}" "${expected_sm}"
  done
}

benchmark_assay cutrun "${CUTRUN_TREATMENT}" "${CUTRUN_CONTROL}" cutrun_full \
  "${EXPECTED_CUTRUN_PEAKS}" "${EXPECTED_CUTRUN_NARROWPEAK_MD5}" "${EXPECTED_CUTRUN_SUMMITS_MD5}"
benchmark_assay chipseq "${CHIP_TREATMENT}" "${CHIP_CONTROL}" chipseq_full \
  "${EXPECTED_CHIP_PEAKS}" "${EXPECTED_CHIP_NARROWPEAK_MD5}" "${EXPECTED_CHIP_SUMMITS_MD5}"

python3 "${REPO_ROOT}/paper/scripts/summarize_paper_benchmarks.py" \
  "${RAW_TSV}" "${OUT_ROOT}/summary.tsv"
awk -F '\t' 'BEGIN {
    OFS="\t";
    print "assay", "macs3_median_s", "rapidmacs_t1_s", "speedup_t1", \
          "rapidmacs_t4_s", "speedup_t4", "rapidmacs_t24_s", "speedup_t24"
  }
  NR > 1 && $1 == "chipseq" {
    printf "chipseq\t%.3f\t96.39\t%.1f\t27.23\t%.1f\t12.11\t%.1f\n", \
           $5, $5/96.39, $5/27.23, $5/12.11
  }
  NR > 1 && $1 == "cutrun" {
    printf "cutrun\t%.3f\t9.81\t%.1f\t2.98\t%.1f\t1.99\t%.1f\n", \
           $5, $5/9.81, $5/2.98, $5/1.99
  }' "${OUT_ROOT}/summary.tsv" >"${OUT_ROOT}/speedups.tsv"

chip_t24_speedup="$(awk -F '\t' '$1 == "chipseq" {print $8}' "${OUT_ROOT}/speedups.tsv")"
if awk -v value="${chip_t24_speedup}" 'BEGIN {exit !(value > 47.1)}'; then
  headline_status="exceeds_47.1x"
else
  headline_status="does_not_exceed_47.1x"
fi
printf 'chipseq_t24_speedup\t%s\nheadline_47.1x_status\t%s\n' \
  "${chip_t24_speedup}" "${headline_status}" >"${OUT_ROOT}/headline_check.tsv"

printf 'completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
  >>"${OUT_ROOT}/provenance/metadata.tsv"
printf 'final_free_bytes\t%s\n' "$(free_bytes)" >>"${OUT_ROOT}/provenance/metadata.tsv"

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
