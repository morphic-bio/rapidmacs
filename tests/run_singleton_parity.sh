#!/usr/bin/env bash
# Frozen original MACS3 3.0.3 references; optionally revalidate or regenerate.
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${RAPIDMACS:-${REPO_ROOT}/bin/rapidmacs}"
GOLDEN="${REPO_ROOT}/tests/golden/singleton_strands"
INPUT="${REPO_ROOT}/tests/fixtures/singleton_strands.bed"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rapidmacs_singleton.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT
if [[ "${UPDATE_GOLDENS:-0}" == 1 && -z "${MACS3:-}" ]]; then
  echo "ERROR: UPDATE_GOLDENS requires MACS3 pointing to original MACS3 3.0.3" >&2
  exit 2
fi
if [[ -n "${MACS3:-}" ]]; then
  [[ "$("${MACS3}" --version)" == 'macs3 3.0.3' ]]
fi

for keep in 1 2 all; do
  for mode in p q; do
    case_name="${mode}_${keep}"
    if [[ "${mode}" == p ]]; then cutoff=1e-5; else cutoff=0.05; fi
    common=(--nomodel --extsize 200 --shift -100 --keep-dup "${keep}"
            "-${mode}" "${cutoff}" -n singleton)
    candidate="${WORK}/${case_name}"
    mkdir -p "${candidate}"
    "${CLI}" --macs-bed-input "${INPUT}" --macs-input-format BED \
      --effective-genome-size 10000 "${common[@]}" --outdir "${candidate}" \
      >"${candidate}/run.log" 2>&1
    if [[ -n "${MACS3:-}" ]]; then
      reference="${WORK}/reference_${case_name}"
      "${MACS3}" callpeak -t "${INPUT}" -f BED -g 10000 "${common[@]}" \
        --outdir "${reference}" >"${candidate}/reference.log" 2>&1
    fi
    for suffix in peaks.narrowPeak summits.bed; do
      file="singleton_${suffix}"
      test -s "${candidate}/${file}"
      if [[ "${UPDATE_GOLDENS:-0}" == 1 ]]; then
        mkdir -p "${GOLDEN}/${case_name}"
        cp "${reference}/${file}" "${GOLDEN}/${case_name}/${file}"
      fi
      if [[ -n "${MACS3:-}" ]]; then
        cmp "${GOLDEN}/${case_name}/${file}" "${reference}/${file}"
      fi
      cmp "${GOLDEN}/${case_name}/${file}" "${candidate}/${file}"
    done
    case "${keep}" in 1) total=73 ;; 2) total=77 ;; all) total=85 ;; esac
    awk -F '\t' -v total="${total}" '
      $1 == "tags_after_filtering" { found = 1; if ($2 != total) exit 1 }
      END { if (!found) exit 1 }
    ' "${candidate}/singleton_chromap_macs_bed_summary.tsv"
    echo "PASS: singleton BED ${case_name}, $(wc -l < "${candidate}/singleton_peaks.narrowPeak") peaks; narrowPeak and summits byte-identical"
  done
done
