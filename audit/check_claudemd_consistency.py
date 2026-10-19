#!/usr/bin/env python3
"""
check_claudemd_consistency.py -- cross-check CLAUDE.md's two dataset tables.

CLAUDE.md states, in a header note directly above the "已分析/已調查資料集"
table, that the "資料集特性稽核校正" table (the correction table) is the
single authoritative source for a dataset's read length and library prep,
and that the comparison table's own Library/Read length columns must be kept
in sync with it. That declaration failed silently three times before anyone
noticed (GSE130078, GSE136569, GSE248612 all drifted back to a stale "150bp"
in the comparison table or narrative text after the correction table was
fixed) -- a declaration that depends on a human remembering to check it by
hand is not a enforcement mechanism. This script is: it parses both tables
directly out of CLAUDE.md and reports every read-length mismatch for a
dataset present in both, with a non-zero exit code on any mismatch, so it
can be run by hand or wired into a pre-commit hook.

Scope: read length only (the one field that both tables actually carry in a
numerically comparable form -- the correction table has no sample-count or
depth column at all, and library-prep text is free-form prose not reliably
diffable field-by-field without false positives). This is deliberately
narrow and exact rather than broad and fuzzy.

Usage: python audit/check_claudemd_consistency.py [path/to/CLAUDE.md]
Exit 0: every dataset present in both tables agrees on read length.
Exit 1: at least one mismatch found (printed to stdout) -- or a parsing
        assumption about either table's row format no longer holds, which
        is itself reported as a mismatch (a parse failure must not be read
        as "nothing wrong").
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

CORRECTION_TABLE_HEADER = "## 資料集特性稽核校正"
COMPARISON_TABLE_HEADER = "### 已分析/已調查資料集"

# Correction table row: | GSE_ID | NNN-NNN | ... | ... | ... | ... |
CORRECTION_ROW_RE = re.compile(
    r'^\|\s*([A-Z]+\d+)\s*\|\s*\**(\d+)-(\d+)\**'
)

# Comparison table row: | **GSE_ID** | ... | ... | ... | <read length cell> | ...
COMPARISON_ROW_RE = re.compile(
    r'^\|\s*\*\*([A-Z]+\d+)\*\*\s*\|'
)
# A read-length cell is either a single "NNNbp" or a range "NNN[-–]NNNbp"
# (hyphen or en-dash, both seen in practice across the two tables -- e.g.
# the correction table uses "138-140", a just-added comparison-table cell
# uses "138–140bp"). Captures both numbers when a range is present, or the
# same number twice for a single value, so range vs range and single vs
# single both compare as plain (low, high) tuples.
READ_LEN_CELL_RE = re.compile(r'(\d+)(?:\s*[-–]\s*(\d+))?\s*bp')


def _extract_range(cell: str) -> tuple[int, int] | None:
    m = READ_LEN_CELL_RE.search(cell)
    if not m:
        return None
    lo = int(m.group(1))
    hi = int(m.group(2)) if m.group(2) else lo
    return (min(lo, hi), max(lo, hi))


def parse_correction_table(lines: list[str]) -> dict[str, tuple[int, int]]:
    """Return {accession: (low_bp, high_bp)} from the correction table's
    "NNN-NNN" cell (both sides equal for ordinary paired-end; a genuine
    range is possible too and is preserved as-is, not averaged)."""
    out = {}
    in_table = False
    for line in lines:
        if CORRECTION_TABLE_HEADER in line:
            in_table = True
            continue
        if in_table and line.startswith("## ") and CORRECTION_TABLE_HEADER not in line:
            break  # left the section
        if not in_table:
            continue
        m = CORRECTION_ROW_RE.match(line)
        if m:
            acc, len1, len2 = m.group(1), int(m.group(2)), int(m.group(3))
            out[acc] = (min(len1, len2), max(len1, len2))
    return out


def parse_comparison_table(lines: list[str]) -> dict[str, tuple[int, int]]:
    """Return {accession: (low_bp, high_bp)} from the comparison table's
    "Read length" column, handling both a single "NNNbp" value and an
    explicit "NNN-NNNbp"/"NNN–NNNbp" range the same way the correction
    table's own cell is parsed, so the two are comparable as equal tuples
    rather than disagreeing purely because of which end of a range each
    table's own formatting convention happens to put right before "bp"."""
    out = {}
    in_table = False
    for line in lines:
        if COMPARISON_TABLE_HEADER in line:
            in_table = True
            continue
        if in_table and line.startswith("### ") and COMPARISON_TABLE_HEADER not in line:
            break
        if not in_table:
            continue
        m = COMPARISON_ROW_RE.match(line)
        if not m:
            continue
        acc = m.group(1)
        cells = line.split("|")
        # Read length is the 5th cell per the known header order:
        # 編號|樣本類型|對比組|Library|Read length|樣本數|...
        if len(cells) < 6:
            continue
        rng = _extract_range(cells[5])
        if rng:
            out[acc] = rng
    return out


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.parent / "CLAUDE.md"
    if not path.exists():
        print(f"STOP: {path} not found")
        return 1
    lines = path.read_text(encoding="utf-8").splitlines()

    correction = parse_correction_table(lines)
    comparison = parse_comparison_table(lines)

    if not correction:
        print("STOP: parsed 0 rows from the correction table -- its format "
              "has likely changed; this script's regex needs updating, "
              "treat this as a failure, not as 'nothing to check'")
        return 1
    if not comparison:
        print("STOP: parsed 0 rows from the comparison table -- same concern")
        return 1

    shared = sorted(set(correction) & set(comparison))
    mismatches = [acc for acc in shared if correction[acc] != comparison[acc]]

    print(f"Correction table: {len(correction)} rows parsed")
    print(f"Comparison table: {len(comparison)} rows parsed")
    print(f"Datasets in both: {len(shared)}")
    print()

    if mismatches:
        print(f"MISMATCH ({len(mismatches)}):")
        for acc in mismatches:
            print(f"  {acc}: correction table says {correction[acc]}bp, "
                  f"comparison table says {comparison[acc]}bp")
        return 1

    print("OK: every dataset present in both tables agrees on read length.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
