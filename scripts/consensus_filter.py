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
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


# ── Type alias ────────────────────────────────────────────────────────────────
# Maps (chr, start, end) → BSJ count  (typing.Dict for Python 3.7 compat)
CoordMap = Dict[Tuple[str, int, int], float]


# ── Parsers ───────────────────────────────────────────────────────────────────

def parse_ciriquant(
    path: str,
    min_bsj: int,
    max_junction_ratio: float = 1.0,
    qc_bsj_threshold: int = 5,
) -> CoordMap:
    """Parse CIRIquant GTF. Returns {(chr, start, end): bsj_count}.

    Pseudo-circ QC: removes entries where BSJ/FSJ > max_junction_ratio AND BSJ < qc_bsj_threshold.
    High-BSJ circRNAs (BSJ >= qc_bsj_threshold) are exempt from ratio filtering —
    they are likely real even if the host gene has low linear expression (FSJ << BSJ).
    FSJ = 0 entries are always skipped from ratio filtering.
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
            bsj_m = re.search(r'bsj\s+([\d.]+)', parts[8], re.IGNORECASE)
            fsj_m = re.search(r'fsj\s+([\d.]+)', parts[8], re.IGNORECASE)
            if not bsj_m:
                continue
            bsj = float(bsj_m.group(1))
            fsj = float(fsj_m.group(1)) if fsj_m else 0.0
            if bsj < min_bsj:
                continue
            # Apply QC only to low-BSJ circRNAs; high-BSJ ones are kept regardless of ratio
            if fsj > 0 and bsj < qc_bsj_threshold and (bsj / fsj) > max_junction_ratio:
                n_ratio_filtered += 1
                continue
            coords[(chrom, start, end)] = bsj
    if n_ratio_filtered:
        print(
            f"[consensus] CIRIquant junction_ratio filter "
            f"(BSJ<{qc_bsj_threshold} AND BSJ/FSJ>{max_junction_ratio}): "
            f"{n_ratio_filtered} removed",
            file=sys.stderr,
        )
    return coords


def parse_dcc(path: str, min_bsj: int) -> CoordMap:
    """Parse DCC output. Prefers CircRNACount (has read counts) over CircCoordinates.
    DCC writes coordinates and counts in separate files; CircRNACount has columns:
    Chr, Start, End, <sample_name> (the junction file name, index 3).
    """
    coords: CoordMap = {}
    # CircRNACount sits in the same directory as CircCoordinates
    count_path = str(Path(path).parent / "CircRNACount")
    parse_path = count_path if Path(count_path).exists() else path

    with open(parse_path) as fh:
        raw_header = fh.readline().strip().split("\t")
        header = [h.lower().strip() for h in raw_header]
        chr_idx   = next((i for i, h in enumerate(header) if h in ("chr", "chrom")), 0)
        start_idx = next((i for i, h in enumerate(header) if h == "start"), 1)
        end_idx   = next((i for i, h in enumerate(header) if h == "end"), 2)
        # CircRNACount: count is 4th column (index 3), header name is the junction file
        # CircCoordinates: no count column; fall back to nominal 5 (passed -Nr 5 1 filter)
        use_count_col = parse_path == count_path
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) <= 2:
                continue
            try:
                chrom = parts[chr_idx]
                start = int(parts[start_idx])
                end   = int(parts[end_idx])
                if use_count_col and len(parts) > 3:
                    count = float(parts[3])
                else:
                    count = 5.0  # DCC already filtered by -Nr 5 1
                if count >= min_bsj:
                    coords[(chrom, start, end)] = count
            except (ValueError, IndexError):
                continue
    return coords


def parse_circexplorer2(path: str, min_bsj: int) -> CoordMap:
    """Parse CIRCexplorer2 output: auto-detects format.

    Handles two formats:
    - back_spliced_junction.bed (≤6 cols): col 4 = BSJ count (score)
    - known_circ.txt / circularRNA_known.txt (≥13 cols): col 12 = readNumber

    CIRCexplorer2 uses BED 0-based start; we add 1 to match CIRIquant GTF (1-based).
    """
    coords: CoordMap = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            n = len(parts)
            if n < 5:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1]) + 1  # BED 0-based → 1-based (match CIRIquant)
                end   = int(parts[2])
                bsj   = float(parts[12]) if n >= 13 else float(parts[4])
                if bsj >= min_bsj:
                    coords[(chrom, start, end)] = bsj
            except (ValueError, IndexError):
                continue
    return coords


def parse_find_circ(path: str, min_bsj: int) -> CoordMap:
    """Parse find_circ splice_sites.bed output.

    Format (BED 0-based): chrom start end name n_reads strand n_uniq uniq_bridges
      best_qual_left best_qual_right tissues tiss_counts edits anchor_overlap
      breakpoints signal strandmatch category
    Column 4 (n_reads) = BSJ read count.
    Column 17 (category) must contain 'CIRCULAR' — filters out LINEAR/AMBIGUOUS junctions.
    Start is converted to 1-based to match CIRIquant GTF coordinates.
    """
    coords: CoordMap = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            # Filter to CIRCULAR category only (col 17); skip if column absent
            if len(parts) >= 18 and "CIRCULAR" not in parts[17]:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1]) + 1  # BED 0-based → 1-based
                end   = int(parts[2])
                bsj   = float(parts[4])    # n_reads column
                if bsj >= min_bsj:
                    coords[(chrom, start, end)] = bsj
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


def _cluster_coords(
    coords: list[tuple[str, int, int]],
    slop: int,
    primary: CoordMap,
) -> list[tuple[str, int, int]]:
    """
    Collapse coordinates lying within `slop` of one another (same chromosome)
    into a single representative each.

    Different detection tools use different aligners (CIRIquant: HISAT2/BWA,
    DCC: STAR), so the *same* back-splice junction is frequently reported at
    coordinates a few bp apart. Without this step, vote() would emit one row per
    tool-specific coordinate for the same biological circRNA, inflating the
    consensus count.

    The representative for a cluster prefers a coordinate that exists in
    `primary` (the CIRIquant map, always tool_maps[0] when present), because the
    downstream count matrix keys on CIRIquant GTF coordinates via an exact-string
    match — picking the CIRIquant coordinate keeps that join intact.
    """
    by_chr: dict[str, list[tuple[str, int, int]]] = {}
    for c in coords:
        by_chr.setdefault(c[0], []).append(c)

    reps: list[tuple[str, int, int]] = []
    for chrom, members in by_chr.items():
        members.sort(key=lambda c: (c[1], c[2]))
        clusters: list[list[tuple[str, int, int]]] = []
        for c in members:
            for cl in clusters:
                anchor = cl[0]  # leftmost coord; members are start-sorted
                if max(abs(c[1] - anchor[1]), abs(c[2] - anchor[2])) <= slop:
                    cl.append(c)
                    break
            else:
                clusters.append([c])
        for cl in clusters:
            reps.append(next((c for c in cl if c in primary), cl[0]))
    return reps


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

    Near-duplicate coordinates (same junction reported a few bp apart by different
    tools) are first clustered to a single representative, so each biological
    circRNA is voted on — and emitted — exactly once.
    """
    all_coords: set[tuple[str, int, int]] = set()
    for m in tool_maps:
        all_coords.update(m.keys())

    primary = tool_maps[0] if tool_maps else {}
    canonical = _cluster_coords(sorted(all_coords), slop, primary)

    results: list[tuple[str, int, int, int, float]] = []

    for coord in canonical:
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
    parser.add_argument("--cirique",       default=None, help="CIRIquant GTF output (optional)")
    parser.add_argument("--ciri2",         default=None, help="CIRI2 output text file (optional)")
    parser.add_argument("--dcc",           default=None, help="DCC CircCoordinates file (optional)")
    parser.add_argument("--circexplorer2", default=None,
                        help="CIRCexplorer2 back_spliced_junction.bed or known_circ.txt (optional)")
    parser.add_argument("--find-circ",     default=None, dest="find_circ",
                        help="find_circ splice_sites.bed output (optional)")
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
    parser.add_argument("--qc-bsj-threshold", type=int, default=5,
                        dest="qc_bsj_threshold",
                        help="BSJ count threshold: only apply ratio QC when BSJ < this value "
                             "(default: 5). High-BSJ circRNAs are exempt from ratio filtering.")
    parser.add_argument("--adaptive", action="store_true", default=False,
                        help="Enable adaptive fallback: if one tool detects far fewer "
                             "circRNAs than another (< --adaptive-ratio), automatically "
                             "lower min_tools to 1 to avoid near-zero recall")
    parser.add_argument("--adaptive-ratio", type=float, default=0.1,
                        dest="adaptive_ratio",
                        help="Threshold for adaptive fallback: if min(tool_counts) / "
                             "max(tool_counts) < ratio, trigger single-tool mode (default: 0.1)")
    parser.add_argument("--min-confidence", type=float, default=0.0,
                        dest="min_confidence",
                        help="Minimum confidence score threshold after voting (default: 0.0). "
                             "Removes low-evidence circRNAs to improve precision.")
    args = parser.parse_args()

    if not any([args.cirique, args.ciri2, args.dcc, args.circexplorer2, args.find_circ]):
        parser.error("至少需要提供 --cirique、--ciri2、--dcc、--circexplorer2 或 --find-circ 其中一個")

    tool_maps: list[CoordMap] = []
    tool_names: list[str] = []
    cirique_map:       CoordMap = {}
    ciri2_map:         CoordMap = {}
    dcc_map:           CoordMap = {}
    circexplorer2_map: CoordMap = {}
    find_circ_map:     CoordMap = {}

    if args.cirique:
        cirique_map = parse_ciriquant(
            args.cirique, args.min_bsj, args.max_junction_ratio, args.qc_bsj_threshold
        )
        tool_maps.append(cirique_map)
        tool_names.append("ciriquant")
        print(f"[consensus] CIRIquant:     {len(cirique_map)} circRNAs", file=sys.stderr)

    if args.ciri2:
        ciri2_map = parse_ciri2(args.ciri2, args.min_bsj)
        tool_maps.append(ciri2_map)
        tool_names.append("ciri2")
        print(f"[consensus] CIRI2:         {len(ciri2_map)} circRNAs", file=sys.stderr)

    if args.dcc:
        dcc_map = parse_dcc(args.dcc, args.min_bsj)
        tool_maps.append(dcc_map)
        tool_names.append("dcc")
        print(f"[consensus] DCC:           {len(dcc_map)} circRNAs", file=sys.stderr)

    if args.circexplorer2:
        circexplorer2_map = parse_circexplorer2(args.circexplorer2, args.min_bsj)
        tool_maps.append(circexplorer2_map)
        tool_names.append("circexplorer2")
        print(f"[consensus] CIRCexplorer2: {len(circexplorer2_map)} circRNAs", file=sys.stderr)

    if args.find_circ:
        find_circ_map = parse_find_circ(args.find_circ, args.min_bsj)
        tool_maps.append(find_circ_map)
        tool_names.append("find_circ")
        print(f"[consensus] find_circ:     {len(find_circ_map)} circRNAs", file=sys.stderr)

    effective_min_tools = min(args.min_tools, len(tool_maps))
    if effective_min_tools < args.min_tools:
        print(
            f"[consensus] min_tools 調整為 {effective_min_tools}（只有 {len(tool_maps)} 個工具）",
            file=sys.stderr,
        )

    # ── Adaptive fallback: if tool detection counts are severely imbalanced ───
    if args.adaptive and effective_min_tools > 1 and len(tool_maps) >= 2:
        counts = [len(m) for m in tool_maps]
        imbalance = min(counts) / max(counts) if max(counts) > 0 else 0.0
        if imbalance < args.adaptive_ratio:
            print(
                f"[consensus] ⚠ Adaptive fallback triggered: "
                f"tool counts = {dict(zip(tool_names, counts))} "
                f"(imbalance ratio = {imbalance:.3f} < {args.adaptive_ratio}). "
                f"Lowering min_tools 2 → 1 to prevent near-zero recall.",
                file=sys.stderr,
            )
            effective_min_tools = 1

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

    # ── Confidence score filter ────────────────────────────────────────────────
    if args.min_confidence > 0.0:
        before = len(results)
        results = [r for r in results if r[4] >= args.min_confidence]
        print(
            f"[consensus] confidence_score ≥ {args.min_confidence}: "
            f"{len(results)} circRNAs 保留（移除 {before - len(results)} 個）",
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
        if circexplorer2_map:
            m = _best_match(coord, circexplorer2_map, args.slop)
            row["in_circexplorer2"] = 1 if m else 0
            row["bsj_circexplorer2"] = m[0] if m else 0
        if find_circ_map:
            m = _best_match(coord, find_circ_map, args.slop)
            row["in_find_circ"] = 1 if m else 0
            row["bsj_find_circ"] = m[0] if m else 0
        rows.append(row)

    pd.DataFrame(rows).to_csv(args.summary, sep="\t", index=False)
    print(f"[consensus] Summary written → {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
