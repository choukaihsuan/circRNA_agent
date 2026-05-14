"""
consensus_filter.py — Multi-tool circRNA consensus voting filter.

Reads output from CIRIquant (.gtf), CIRI2 (.txt), and DCC (CircCoordinates),
votes across tools, and outputs high-confidence circRNAs.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


# ── Parsers ──────────────────────────────────────────────────────────────────

def parse_ciriquant(path: str, min_bsj: int) -> set[tuple[str, int, int]]:
    """Parse CIRIquant GTF output. Returns set of (chr, start, end)."""
    coords: set[tuple[str, int, int]] = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "circRNA":
                continue
            chrom = parts[0]
            start = int(parts[3])
            end   = int(parts[4])
            bsj_match = re.search(r'BSJ\s+([\d.]+)', parts[8])
            if not bsj_match:
                continue
            if float(bsj_match.group(1)) >= min_bsj:
                coords.add((chrom, start, end))
    return coords


def parse_ciri2(path: str, min_bsj: int) -> set[tuple[str, int, int]]:
    """Parse CIRI2 output text file. Returns set of (chr, start, end)."""
    coords: set[tuple[str, int, int]] = set()
    with open(path) as fh:
        raw_header = fh.readline().strip().split("\t")
        header = [h.lower().strip() for h in raw_header]

        chr_idx = next((i for i, h in enumerate(header) if h == "chr" or "chrom" in h), 0)
        start_idx = next(
            (i for i, h in enumerate(header) if "start" in h and "circRNA_start" not in h or h == "circrna_start"), 1
        )
        end_idx = next(
            (i for i, h in enumerate(header) if "end" in h and "circrna_end" not in h or h == "circrna_end"), 2
        )
        bsj_idx = next(
            (i for i, h in enumerate(header) if "junction_reads" in h or h == "bsj" or "#junction_reads" in h), 4
        )

        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) <= max(chr_idx, start_idx, end_idx, bsj_idx):
                continue
            try:
                # CIRI2 circ_id format: chr:start|end
                circ_id = parts[0]
                m = re.match(r'^(.+):(\d+)\|(\d+)', circ_id)
                if m:
                    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                else:
                    chrom = parts[chr_idx]
                    start = int(parts[start_idx])
                    end   = int(parts[end_idx])
                bsj = int(parts[bsj_idx])
                if bsj >= min_bsj:
                    coords.add((chrom, start, end))
            except (ValueError, IndexError):
                continue
    return coords


def parse_dcc(path: str, min_bsj: int) -> set[tuple[str, int, int]]:
    """Parse DCC CircCoordinates file. Returns set of (chr, start, end)."""
    coords: set[tuple[str, int, int]] = set()
    with open(path) as fh:
        raw_header = fh.readline().strip().split("\t")
        header = [h.lower().strip() for h in raw_header]

        chr_idx   = next((i for i, h in enumerate(header) if h == "chr"), 0)
        start_idx = next((i for i, h in enumerate(header) if h == "start"), 1)
        end_idx   = next((i for i, h in enumerate(header) if h == "end"), 2)
        count_idx = next((i for i, h in enumerate(header) if h == "c"), 3)

        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) <= max(chr_idx, start_idx, end_idx, count_idx):
                continue
            try:
                chrom = parts[chr_idx]
                start = int(parts[start_idx])
                end   = int(parts[end_idx])
                count = int(parts[count_idx])
                if count >= min_bsj:
                    coords.add((chrom, start, end))
            except (ValueError, IndexError):
                continue
    return coords


# ── Voting ───────────────────────────────────────────────────────────────────

def coords_match(
    c1: tuple[str, int, int],
    c2: tuple[str, int, int],
    slop: int,
) -> bool:
    return c1[0] == c2[0] and abs(c1[1] - c2[1]) <= slop and abs(c1[2] - c2[2]) <= slop


def vote(
    sets: list[set[tuple[str, int, int]]],
    min_tools: int,
    slop: int,
) -> list[tuple[str, int, int, int]]:
    """Return circRNAs supported by >= min_tools tools. Fourth element is n_tools."""
    all_coords: set[tuple[str, int, int]] = set()
    for s in sets:
        all_coords.update(s)

    results: list[tuple[str, int, int, int]] = []
    for coord in all_coords:
        n_tools = sum(
            any(coords_match(coord, c, slop) for c in s)
            for s in sets
        )
        if n_tools >= min_tools:
            results.append((*coord, n_tools))  # type: ignore[arg-type]

    return sorted(results, key=lambda x: (x[0], x[1]))


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Multi-tool circRNA consensus filter (CIRIquant + CIRI2 + DCC)"
    )
    parser.add_argument("--cirique",   required=True,  help="CIRIquant GTF output")
    parser.add_argument("--ciri2",     default=None,   help="CIRI2 output text file (optional)")
    parser.add_argument("--dcc",       required=True,  help="DCC CircCoordinates file")
    parser.add_argument("--output",    required=True, help="Output BED file")
    parser.add_argument("--summary",   required=True, help="Output summary TSV file")
    parser.add_argument("--min-tools", type=int, default=2, dest="min_tools",
                        help="Minimum number of tools that must agree (default: 2)")
    parser.add_argument("--slop",      type=int, default=10,
                        help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--min-bsj",   type=int, default=2, dest="min_bsj",
                        help="Minimum BSJ reads per tool before voting (default: 2)")
    args = parser.parse_args()

    cirique_set = parse_ciriquant(args.cirique, args.min_bsj)
    ciri2_set   = parse_ciri2(args.ciri2, args.min_bsj) if args.ciri2 else set()
    dcc_set     = parse_dcc(args.dcc, args.min_bsj)

    tool_sets = [cirique_set, dcc_set]
    if ciri2_set:
        tool_sets.append(ciri2_set)

    all_count = len(set().union(*tool_sets))

    print(f"[consensus] CIRIquant: {len(cirique_set)} circRNAs", file=sys.stderr)
    if ciri2_set:
        print(f"[consensus] CIRI2:     {len(ciri2_set)} circRNAs", file=sys.stderr)
    print(f"[consensus] DCC:       {len(dcc_set)} circRNAs",     file=sys.stderr)

    results = vote(tool_sets, args.min_tools, args.slop)

    print(
        f"[consensus] ≥{args.min_tools} 工具共識: "
        f"{len(results)} circRNAs（過濾掉 {all_count - len(results)} 個）",
        file=sys.stderr,
    )

    # BED output: chr start end . n_tools .
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        for chrom, start, end, n_tools in results:
            fh.write(f"{chrom}\t{start}\t{end}\t.\t{n_tools}\t.\n")

    # Summary TSV
    rows = []
    for chrom, start, end, n_tools in results:
        coord = (chrom, start, end)
        rows.append({
            "circ_id":      f"{chrom}:{start}|{end}",
            "chr":          chrom,
            "start":        start,
            "end":          end,
            "n_tools":      n_tools,
            "in_ciriquant": int(any(coords_match(coord, c, args.slop) for c in cirique_set)),
            "in_ciri2":     int(any(coords_match(coord, c, args.slop) for c in ciri2_set)),
            "in_dcc":       int(any(coords_match(coord, c, args.slop) for c in dcc_set)),
        })

    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
    print(f"[consensus] Summary written → {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
