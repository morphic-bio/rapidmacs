#!/usr/bin/env bash
# Build the MACS3 3.0.3 reference installations the paper measures against.
#
# Two builds are produced, because the paper reports both:
#
#   macs3-3.0.3-venv                  released 3.0.3, unmodified. This is what
#                                     Table 1 and Supplementary S3 report.
#   macs3-3.0.3-no-control-dump-venv  the same release with the four lines that
#                                     write tmp_b_c.txt removed. Used for
#                                     Supplementary S4, and for the three
#                                     treatment-only configurations, which never
#                                     reach that code.
#
# Both are built from the same PyPI sdist, pinned by SHA-256.
#
# Usage:
#   benchmarks/setup_macs3_reference.sh                      # both builds
#   VARIANT=released benchmarks/setup_macs3_reference.sh
#   VARIANT=no-control-dump benchmarks/setup_macs3_reference.sh
#
# Environment:
#   DATA_ROOT  where to unpack sources and create venvs
#              (default: <repo>/.macs3-reference)
#   VARIANT    released | no-control-dump | both   (default both)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${DATA_ROOT:-${REPO_ROOT}/.macs3-reference}"
VARIANT="${VARIANT:-both}"
SOURCE_ROOT="${DATA_ROOT}/macs3-reference-src"
ARCHIVE="${SOURCE_ROOT}/macs3-3.0.3.tar.gz"
URL='https://files.pythonhosted.org/packages/c1/fe/c557e1170f831876c8673c067c5258326312e6e34b22c228f3e5d72857be/macs3-3.0.3.tar.gz'
SHA256='ee1c892901c4010ff9e201b433c0623cbd747a3058300322386a7185623b1684'

case "${VARIANT}" in
  released|no-control-dump|both) ;;
  *) echo "ERROR: VARIANT must be released, no-control-dump, or both" >&2; exit 2 ;;
esac

mkdir -p "${SOURCE_ROOT}"
if [[ ! -s "${ARCHIVE}" ]]; then curl -L --fail "${URL}" -o "${ARCHIVE}"; fi
printf '%s  %s\n' "${SHA256}" "${ARCHIVE}" | sha256sum -c -

build() {
  local variant="$1" src venv track reported
  src="${SOURCE_ROOT}/macs3-3.0.3-${variant}"
  if [[ "${variant}" == released ]]; then
    venv="${DATA_ROOT}/macs3-3.0.3-venv"
  else
    venv="${DATA_ROOT}/macs3-3.0.3-${variant}-venv"
  fi

  rm -rf "${src}"
  mkdir -p "${src}"
  tar -xzf "${ARCHIVE}" -C "${src}" --strip-components=1

  track="${src}/MACS3/Signal/PairedEndTrack.py"
  if [[ "${variant}" == no-control-dump ]]; then
    grep -q 'tmp_b_c.txt' "${track}" \
      || { echo "ERROR: sdist no longer contains the tmp_b_c.txt write" >&2; exit 1; }
    patch -d "${src}" -p1 <"${REPO_ROOT}/benchmarks/macs3-3.0.3-no-control-dump.patch"
    if grep -q 'tmp_b_c.txt' "${track}"; then
      echo "ERROR: patch did not remove the write" >&2; exit 1
    fi
  else
    grep -q 'tmp_b_c.txt' "${track}" \
      || { echo "ERROR: released source should still contain the write" >&2; exit 1; }
  fi

  rm -rf "${venv}"
  python3 -m venv "${venv}"
  "${venv}/bin/pip" install --disable-pip-version-check "${src}"
  reported="$("${venv}/bin/macs3" --version 2>&1 | awk 'NR == 1 {print $NF}')"
  [[ "${reported}" == 3.0.3 ]] || { echo "ERROR: built ${reported}, expected 3.0.3" >&2; exit 1; }
  echo "[built] ${variant}: ${venv}/bin/macs3"
}

if [[ "${VARIANT}" == released || "${VARIANT}" == both ]]; then build released; fi
if [[ "${VARIANT}" == no-control-dump || "${VARIANT}" == both ]]; then build no-control-dump; fi

echo
echo "Point the benchmark scripts at the build you want with MACS3=<path>."
echo "Table 1 and Supplementary S3 report the released build."
