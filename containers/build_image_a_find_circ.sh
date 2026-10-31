#!/usr/bin/env bash
# Fetch, verify, and patch find_circ for Image A (the paper-reproduction
# image -- see containers/build_notes.md's "Two images" section).
#
# find_circ is GPLv3 (github.com/marvin-jens/find_circ); this repo is MIT.
# Nothing third-party is vendored here -- this script fetches the upstream
# file fresh at build time and applies patches/find_circ_py3.patch to it.
# See that patch file's own header for the full rationale and both sha256
# values this script checks against.
#
# Usage: build_image_a_find_circ.sh <destination_directory>
# Writes <destination_directory>/find_circ.py (patched, executable).
#
# Both sha256 checks below are HARD FAILURES (non-zero exit), not
# warnings -- a warning in a build log is something nobody reads, and
# either of these not matching means the image would ship a find_circ
# that is not actually the one that produced the 19 published datasets.
set -euo pipefail

UPSTREAM_COMMIT="8655dca54970fcf7e92e22fbf57e1188724dda7d"
UPSTREAM_URL="https://raw.githubusercontent.com/marvin-jens/find_circ/${UPSTREAM_COMMIT}/find_circ.py"
UPSTREAM_SHA256="b67a043971c1b8eaa1eb672a1079e2cfc192a0c45ad56124538e8169e67ee174"
PATCHED_SHA256="28a3cb9e99a11a882c26bae250f45235dff3f3fc7594d427ea24dcf553434f00"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH_FILE="${SCRIPT_DIR}/../patches/find_circ_py3.patch"

if [ $# -ne 1 ]; then
    echo "usage: $0 <destination_directory>" >&2
    exit 2
fi
DEST_DIR="$1"
mkdir -p "$DEST_DIR"
TMP_ORIG="$(mktemp)"
trap 'rm -f "$TMP_ORIG"' EXIT

echo "[build_image_a_find_circ] fetching upstream find_circ.py @ ${UPSTREAM_COMMIT}"
curl -fsSL "$UPSTREAM_URL" -o "$TMP_ORIG"

echo "[build_image_a_find_circ] verifying upstream sha256"
ACTUAL_UPSTREAM_SHA256="$(sha256sum "$TMP_ORIG" | cut -d' ' -f1)"
if [ "$ACTUAL_UPSTREAM_SHA256" != "$UPSTREAM_SHA256" ]; then
    echo "[build_image_a_find_circ] FATAL: upstream sha256 mismatch." >&2
    echo "  expected: $UPSTREAM_SHA256" >&2
    echo "  got:      $ACTUAL_UPSTREAM_SHA256" >&2
    echo "  The upstream file at commit ${UPSTREAM_COMMIT} no longer matches" >&2
    echo "  what this patch was written against (or the fetch is corrupt)." >&2
    echo "  Refusing to proceed -- do not vendor a file around this check." >&2
    exit 1
fi
echo "[build_image_a_find_circ] upstream sha256 OK"

echo "[build_image_a_find_circ] applying patches/find_circ_py3.patch"
cp "$TMP_ORIG" "${DEST_DIR}/find_circ.py"
patch "${DEST_DIR}/find_circ.py" < "$PATCH_FILE"

echo "[build_image_a_find_circ] verifying patched sha256"
ACTUAL_PATCHED_SHA256="$(sha256sum "${DEST_DIR}/find_circ.py" | cut -d' ' -f1)"
if [ "$ACTUAL_PATCHED_SHA256" != "$PATCHED_SHA256" ]; then
    echo "[build_image_a_find_circ] FATAL: patched-file sha256 mismatch." >&2
    echo "  expected: $PATCHED_SHA256" >&2
    echo "  got:      $ACTUAL_PATCHED_SHA256" >&2
    echo "  The result of applying this patch does not match the file that" >&2
    echo "  actually produced the 19 published datasets on the server." >&2
    echo "  Refusing to proceed." >&2
    rm -f "${DEST_DIR}/find_circ.py"
    exit 1
fi
echo "[build_image_a_find_circ] patched sha256 OK -- bit-identical to the server's find_circ.py"

chmod +x "${DEST_DIR}/find_circ.py"
echo "[build_image_a_find_circ] done: ${DEST_DIR}/find_circ.py"
