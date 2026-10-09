#!/usr/bin/env bash
# Minimal reproducible example: dual-tool (CIRIquant + DCC) consensus on hand-made toy data.
# Needs only Python 3 (standard library) - no aligner, no genome, no network. Runs in < 1 second.
#
#   bash examples/minimal/run.sh            # writes to examples/minimal/output/ and diffs against expected/
set -euo pipefail
cd "$(dirname "$0")/../.."
OUT=examples/minimal/output
mkdir -p "$OUT"

python3 scripts/consensus_filter.py \
    --cirique examples/minimal/input/SAMPLE.gtf \
    --dcc     examples/minimal/input/DCC/CircCoordinates \
    --output  "$OUT/high_confidence.bed" \
    --summary "$OUT/consensus_summary.tsv" \
    --min-tools 2 --slop 10 --min-bsj 2 \
    --max-junction-ratio 1.0 --qc-bsj-threshold 5 \
    --adaptive --adaptive-ratio 0.1

diff -u examples/minimal/expected/high_confidence.bed   "$OUT/high_confidence.bed"
diff -u examples/minimal/expected/consensus_summary.tsv "$OUT/consensus_summary.tsv"
echo "OK: output matches examples/minimal/expected/"
