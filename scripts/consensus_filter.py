"""
consensus_filter.py — Multi-tool circRNA consensus voting filter with adaptive scoring.

Reads output from CIRIquant (.gtf) and DCC (CircCoordinates),
votes across tools, computes per-circRNA confidence scores, and outputs
high-confidence circRNAs.

Confidence score formula (per circRNA):
    score = Σ over agreeing tools: log2(BSJ+1) × (1 − coord_dist/slop)
            ÷ n_tools_total
Higher score = higher BSJ expression AND more precise coordinate agreement.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd


# ── Type alias ────────────────────────────────────────────────────────────────
# Maps (chr, start, end) → BSJ count
CoordMap = dict[tuple[str, int, int], float]


# ── Parsers ───────────────────────────────────────────────────────────────────

def parse_ciriquant(path: str, min_bsj: int, max_junction_ratio: float = 1.0) -> CoordMap:
    """Parse CIRIquant GTF. Returns {(chr, start, end): bsj_count}.

    Pseudo-circ QC: removes entries where BSJ/FSJ > max_junction_ratio.
    Real circRNAs almost always have fewer BSJ than FSJ reads; a ratio > 1.0
    is a strong indicator of mis-mapping or genomic repeat artifacts.
    FSJ = 0 entries are skipped from ratio filtering (no linear evidence present).
    """
    coords: CoordMap = {}
    n_ratio_filtered = 0
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
            bsj_m = re.search(r'BSJ\s+([\d.]+)', parts[8])
            fsj_m = re.search(r'FSJ\s+([\d.]+)', parts[8])
            if not bsj_m:
                continue
            bsj = float(bsj_m.group(1))
            fsj = float(fsj_m.group(1)) if fsj_m else 0.0
            if bsj < min_bsj:
                continue
            if fsj > 0 and (bsj / fsj) > max_junction_ratio:
                n_ratio_filtered += 1
                continue
            coords[(chrom, start, end)] = bsj
    if n_ratio_filtered:
        print(
            f"[consensus] CIRIquant junction_ratio filter (BSJ/FSJ > {max_junction_ratio}): "
            f"{n_ratio_filtered} removed",
            file=sys.stderr,
        )
    return coords


def parse_dcc(path: str, min_bsj: int) -> CoordMap:
    """Parse DCC CircCoordinates. Returns {(chr, start, end): count}."""
    coords: CoordMap = {}
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
                count = float(parts[count_idx])
                if count >= min_bsj:
                    coords[(chrom, start, end)] = count
            except (ValueError, IndexError):
                continue
    return coords


def parse_ciri2(path: str, min_bsj: int) -> CoordMap:
    """Parse CIRI2 output text file. Returns {(chr, start, end): bsj}."""
    coords: CoordMap = {}
    with open(path) as fh:
        raw_header = fh.readline().strip().split("\t")
        header = [h.lower().strip() for h in raw_header]
        chr_idx   = next((i for i, h in enumerate(header) if h == "chr" or "chrom" in h), 0)
        start_idx = next((i for i, h in enumerate(header) if "start" in h), 1)
        end_idx   = next((i for i, h in enumerate(header) if "end" in h), 2)
        bsj_idx   = next((i for i, h in enumerate(header)
                          if "junction_reads" in h or h == "bsj" or "#junction_reads" in h), 4)
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) <= max(chr_idx, start_idx, end_idx, bsj_idx):
                continue
            try:
                circ_id = parts[0]
                m = re.match(r'^(.+):(\d+)\|(\d+)', circ_id)
                if m:
                    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                else:
                    chrom = parts[chr_idx]; start = int(parts[start_idx]); end = int(parts[end_idx])
                bsj = float(parts[bsj_idx])
                if bsj >= min_bsj:
                    coords[(chrom, start, end)] = bsj
            except (ValueError, IndexError):
                continue
    return coords


# ── Adaptive voting ───────────────────────────────────────────────────────────

def _best_match(
    ref: tuple[str, int, int],
    pool: CoordMap,
    slop: int,
) -> tuple[float, float] | None:
    """Find best matching coordinate in pool within slop. Returns (bsj, dist)."""
    best_bsj, best_dist = None, float("inf")
    for c2, bsj in pool.items():
        if c2[0] != ref[0]:
            continue
        dist = max(abs(c2[1] - ref[1]), abs(c2[2] - ref[2]))
        if dist <= slop and dist < best_dist:
            best_dist = dist
            best_bsj = bsj
    if best_bsj is None:
        return None
    return best_bsj, best_dist


def vote(
    tool_maps: list[CoordMap],
    min_tools: int,
    slop: int,
) -> list[tuple[str, int, int, int, float]]:
    """
    Return circRNAs supported by >= min_tools tools.
    Output: (chr, start, end, n_tools, confidence_score)

    confidence_score = Σ[log2(bsj+1) × (1 − dist/slop)] / n_supporting_tools
    Denominator is the number of tools that actually support this circRNA,
    so the score reflects mean per-tool evidence rather than dilution by absent tools.
    Higher score = more expression + better coordinate agreement among supporting tools.
    """
    all_coords: set[tuple[str, int, int]] = set()
    for m in tool_maps:
        all_coords.update(m.keys())

    results: list[tuple[str, int, int, int, float]] = []

    for coord in all_coords:
        matches = [_best_match(coord, m, slop) for m in tool_maps]
        supporting = [(bsj, dist) for m in matches if m is not None for bsj, dist in [m]]
        n_tools = len(supporting)
        if n_tools < min_tools:
            continue

        score = sum(
            math.log2(bsj + 1) * (1.0 - dist / max(slop, 1))
            for bsj, dist in supporting
        ) / n_tools

        results.append((*coord, n_tools, round(score, 4)))  # type: ignore[arg-type]

    return sorted(results, key=lambda x: (-x[4], x[0], x[1]))


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Multi-tool circRNA consensus filter with adaptive confidence scoring"
    )
    parser.add_argument("--cirique",   default=None, help="CIRIquant GTF output (optional)")
    parser.add_argument("--ciri2",     default=None, help="CIRI2 output text file (optional)")
    parser.add_argument("--dcc",       default=None, help="DCC CircCoordinates file (optional)")
    parser.add_argument("--output",    required=True, help="Output BED file")
    parser.add_argument("--summary",   required=True, help="Output summary TSV file")
    parser.add_argument("--min-tools", type=int, default=2, dest="min_tools",
                        help="Minimum tools that must agree (default: 2)")
    parser.add_argument("--slop",      type=int, default=10,
                        help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--min-bsj",   type=int, default=2, dest="min_bsj",
                        help="Minimum BSJ reads per tool before voting (default: 2)")
    parser.add_argument("--max-junction-ratio", type=float, default=1.0,
                        dest="max_junction_ratio",
                        help="Max BSJ/FSJ ratio for CIRIquant pseudo-circ QC (default: 1.0)")
    args = parser.parse_args()

    if not args.cirique and not args.dcc:
        parser.error("At少需要提供 --cirique 或 --dcc 其中一個")

    tool_maps: list[CoordMap] = []
    tool_names: list[str] = []
    cirique_map: CoordMap = {}
    ciri2_map:   CoordMap = {}
    dcc_map:     CoordMap = {}

    if args.cirique:
        cirique_map = parse_ciriquant(args.cirique, args.min_bsj, args.max_junction_ratio)
        tool_maps.append(cirique_map)
        tool_names.append("ciriquant")
        print(f"[consensus] CIRIquant: {len(cirique_map)} circRNAs", file=sys.stderr)

    if args.ciri2:
        ciri2_map = parse_ciri2(args.ciri2, args.min_bsj)
        tool_maps.append(ciri2_map)
        tool_names.append("ciri2")
        print(f"[consensus] CIRI2:     {len(ciri2_map)} circRNAs", file=sys.stderr)

    if args.dcc:
        dcc_map = parse_dcc(args.dcc, args.min_bsj)
        tool_maps.append(dcc_map)
        tool_names.append("dcc")
        print(f"[consensus] DCC:       {len(dcc_map)} circRNAs", file=sys.stderr)

    effective_min_tools = min(args.min_tools, len(tool_maps))
    if effective_min_tools < args.min_tools:
        print(
            f"[consensus] min_tools 調整為 {effective_min_tools}（只有 {len(tool_maps)} 個工具）",
            file=sys.stderr,
        )

    all_count = len(set().union(*[set(m.keys()) for m in tool_maps]))
    results = vote(tool_maps, effective_min_tools, args.slop)

    print(
        f"[consensus] ≥{effective_min_tools} 工具共識: "
        f"{len(results)} circRNAs（過濾掉 {all_count - len(results)} 個）",
        file=sys.stderr,
    )
    if results:
        scores = [r[4] for r in results]
        print(
            f"[consensus] confidence score: min={min(scores):.3f} "
            f"median={sorted(scores)[len(scores)//2]:.3f} max={max(scores):.3f}",
            file=sys.stderr,
        )

    # ── BED output: chr start end . n_tools . confidence_score ───────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as fh:
        for chrom, start, end, n_tools, score in results:
            fh.write(f"{chrom}\t{start}\t{end}\t.\t{n_tools}\t.\t{score}\n")

    # ── Summary TSV ───────────────────────────────────────────────────────────
    def _in(coord, pool, slop):
        return int(_best_match(coord, pool, slop) is not None)

    rows = []
    for chrom, start, end, n_tools, score in results:
        coord = (chrom, start, end)
        row = {
            "circ_id":          f"{chrom}:{start}|{end}",
            "chr":              chrom,
            "start":            start,
            "end":              end,
            "n_tools":          n_tools,
            "confidence_score": score,
        }
        if cirique_map:
            m = _best_match(coord, cirique_map, args.slop)
            row["in_ciriquant"] = 1 if m else 0
            row["bsj_ciriquant"] = m[0] if m else 0
        if ciri2_map:
            m = _best_match(coord, ciri2_map, args.slop)
            row["in_ciri2"] = 1 if m else 0
            row["bsj_ciri2"] = m[0] if m else 0
        if dcc_map:
            m = _best_match(coord, dcc_map, args.slop)
            row["in_dcc"] = 1 if m else 0
            row["bsj_dcc"] = m[0] if m else 0
        rows.append(row)

    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
    print(f"[consensus] Summary written → {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
