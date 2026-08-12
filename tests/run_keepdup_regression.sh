#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLI="${REPO_ROOT}/bin/rapidmacs"
INPUT="${REPO_ROOT}/tests/fixtures/keepdup_three_coordinate_copies.bedpe"
FRAG_INPUT="${REPO_ROOT}/tests/fixtures/keepdup_warning.fragments.tsv"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/rapidmacs_keepdup.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

"${CLI}" -f BEDPE -t "${INPUT}" \
  --frag-span-pileup-bdg "${WORK}/default.bdg" \
  >"${WORK}/default.stdout" 2>"${WORK}/default.stderr"
grep -Fq "4 fragments parsed, 2 retained, 2 coordinate duplicates removed" \
  "${WORK}/default.stderr"

"${CLI}" -f BEDPE -t "${INPUT}" --keep-dup 2 \
  --frag-span-pileup-bdg "${WORK}/two.bdg" \
  >"${WORK}/two.stdout" 2>"${WORK}/two.stderr"
grep -Fq "4 fragments parsed, 3 retained, 1 coordinate duplicates removed" \
  "${WORK}/two.stderr"

"${CLI}" -f BEDPE -t "${INPUT}" --keep-dup all \
  --frag-span-pileup-bdg "${WORK}/all.bdg" \
  >"${WORK}/all.stdout" 2>"${WORK}/all.stderr"
grep -Fq "4 fragments parsed, 4 retained, 0 coordinate duplicates removed" \
  "${WORK}/all.stderr"

for bad in 0 -1 1.5 2147483648 auto junk; do
  if "${CLI}" -f BEDPE -t "${INPUT}" --keep-dup "${bad}" \
      --frag-span-pileup-bdg "${WORK}/bad.bdg" \
      >"${WORK}/bad.stdout" 2>"${WORK}/bad.stderr"; then
    echo "FAIL: --keep-dup ${bad} was accepted" >&2
    exit 1
  fi
  grep -Fq -- "--keep-dup must be 'all' or a positive integer" \
    "${WORK}/bad.stderr"
done

for value in 2 0 -1 auto junk; do
  "${CLI}" -f FRAG -t "${FRAG_INPUT}" --keep-dup "${value}" \
    --frag-span-pileup-bdg "${WORK}/frag-${value}.bdg" \
    >"${WORK}/frag.stdout" 2>"${WORK}/frag.stderr"
  grep -Fq "Since the format is 'FRAG', \`--keep-dup\` will be set as 'all'." \
    "${WORK}/frag.stderr"
done

echo "PASS: BEDPE --keep-dup values and FRAG warn-and-override regression"
