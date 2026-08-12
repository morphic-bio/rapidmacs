#!/usr/bin/env bash
# Compile and run all RapidMACS unit tests in tests/unit/.
# Each .cc is a self-contained binary linking against lib/librapidmacs.a.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIB="${REPO_ROOT}/lib/librapidmacs.a"
OUTDIR="${REPO_ROOT}/build/tests"

if [[ ! -f "${LIB}" ]]; then
  echo "ERROR: ${LIB} not built. Run 'make' first." >&2
  exit 2
fi

mkdir -p "${OUTDIR}"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++11 -Wall -Wextra -O2 -fopenmp}"
HTSLIB_DIR="${HTSLIB_DIR:-}"

CPPFLAGS="-I${REPO_ROOT}/include"
LDFLAGS="-lhts -lm -lz -lpthread -fopenmp"
if [[ -n "${HTSLIB_DIR}" ]]; then
  CPPFLAGS="${CPPFLAGS} -I${HTSLIB_DIR}"
  LDFLAGS="-L${HTSLIB_DIR} ${LDFLAGS}"
fi

failed=0
total=0
for src in "${REPO_ROOT}"/tests/unit/*.cc; do
  name="$(basename "${src}" .cc)"
  out="${OUTDIR}/${name}"
  total=$((total + 1))
  if ! ${CXX} ${CXXFLAGS} ${CPPFLAGS} "${src}" "${LIB}" ${LDFLAGS} -o "${out}" 2> "${OUTDIR}/${name}.build.log"; then
    echo "FAIL build ${name}" >&2
    cat "${OUTDIR}/${name}.build.log" >&2
    failed=$((failed + 1))
    continue
  fi
  if ! "${out}" > "${OUTDIR}/${name}.run.log" 2>&1; then
    echo "FAIL run ${name}" >&2
    cat "${OUTDIR}/${name}.run.log" >&2
    failed=$((failed + 1))
    continue
  fi
  cat "${OUTDIR}/${name}.run.log"
done

if [[ ${failed} -gt 0 ]]; then
  echo "FAIL: ${failed}/${total} unit tests failed"
  exit 1
fi
echo "PASS: ${total}/${total} unit tests"
