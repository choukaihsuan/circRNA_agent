"""
assign_isoforms.py
──────────────────
Map circRNA BSJ coordinates to host genes from the GTF, producing an
isoform grouping table for downstream IUI (Isoform Usage Index) analysis.

circ_id format expected in count_matrix row names: chrom:start|end
  (1-based coordinates, matching CIRIquant GTF output)

Usage:
    python scripts/assign_isoforms.py \
        --count-matrix results/circRNA/count_matrix.tsv \
        --gtf          /path/to/hg19.gtf \
        --out          results/circRNA/isoform_groups.tsv
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


# ── GTF parsing ───────────────────────────────────────────────────────────────

def _parse_attrs(attr_str: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for token in attr_str.split(";"):
        token = token.strip()
        if not token:
            continue
        kv = token.split(" ", 1)
        if len(kv) == 2:
            attrs[kv[0]] = kv[1].strip('"')
    return attrs


def parse_gtf_genes(gtf_path: str) -> pd.DataFrame:
    """
    Parse gene-level intervals from a GTF.
    Tries 'gene' feature first; falls back to deriving boundaries from
    'transcript' entries grouped by gene_id (handles UCSC-style GTFs).
    Returns DataFrame with columns:
        chrom, g_start (1-based), g_end (1-based), strand, gene_id, gene_name
    """
    gene_rows: List[dict] = []
    # Accumulate transcript boundaries per gene as fallback
    tx_bounds: Dict[str, dict] = defaultdict(lambda: {
        "chrom": "", "g_start": 10**12, "g_end": 0,
        "strand": ".", "gene_name": "",
    })

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("gene", "transcript"):
                continue

            chrom  = parts[0]
            start  = int(parts[3])   # 1-based inclusive
            end    = int(parts[4])   # 1-based inclusive
            strand = parts[6]
            attrs  = _parse_attrs(parts[8])
            gene_id   = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", gene_id)

            if feat == "gene":
                gene_rows.append({
                    "chrom": chrom, "g_start": start, "g_end": end,
                    "strand": strand, "gene_id": gene_id, "gene_name": gene_name,
                })
            else:
                rec = tx_bounds[gene_id]
                rec["chrom"]     = chrom
                rec["g_start"]   = min(rec["g_start"], start)
                rec["g_end"]     = max(rec["g_end"],   end)
                rec["strand"]    = strand
                rec["gene_name"] = gene_name

    if gene_rows:
        df = pd.DataFrame(gene_rows)
        print(f"[assign_isoforms] {len(df):,} gene features from GTF", file=sys.stderr)
        return df

    # Fallback: build gene intervals from transcript entries
    rows = []
    for gid, rec in tx_bounds.items():
        if rec["chrom"] and rec["g_end"] > 0:
            rows.append({"gene_id": gid, **rec})
    df = pd.DataFrame(rows)
    print(
        f"[assign_isoforms] {len(df):,} genes derived from transcript features (fallback)",
        file=sys.stderr,
    )
    return df


# ── Exon parsing & exon-span annotation ──────────────────────────────────────

def parse_gtf_exons(gtf_path: str) -> Dict[str, List[Tuple[int, int, str, str]]]:
    """
    Parse exon features from GTF.
    Returns dict: gene_id → list of (start, end, transcript_id, exon_number).
    exon_number is '' when the attribute is absent.
    """
    gene_exons: Dict[str, List[Tuple[int, int, str, str]]] = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            start  = int(parts[3])
            end    = int(parts[4])
            attrs  = _parse_attrs(parts[8])
            gid    = attrs.get("gene_id", "")
            tx_id  = attrs.get("transcript_id", "")
            exon_n = attrs.get("exon_number", "")
            if gid:
                gene_exons[gid].append((start, end, tx_id, exon_n))
    n_genes = len(gene_exons)
    print(f"[assign_isoforms] Exon index built: {n_genes:,} genes", file=sys.stderr)
    return gene_exons


def _exon_span(circ_start: int, circ_end: int, gene_id: str,
               gene_exons: Dict[str, List[Tuple[int, int, str, str]]],
               slop: int = 10) -> str:
    """
    Find the exon numbers that form the circRNA backsplice boundaries.

    For a circRNA (start, end):
      - circ_start should match an exon's *start* coordinate (upstream donor)
      - circ_end   should match an exon's *end*   coordinate (downstream acceptor)

    Tries every transcript of the host gene; picks the transcript where both
    endpoints are explained.  Reports 'eN-eM' using exon_number attributes, or
    positional rank within the transcript if exon_number is absent.
    Returns '' when no matching transcript is found.
    """
    exons = gene_exons.get(gene_id, [])
    if not exons:
        return ""

    # Group by transcript
    by_tx: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for s, e, tx, en in exons:
        by_tx[tx].append((s, e, en))

    for tx, tx_exons in by_tx.items():
        tx_sorted = sorted(tx_exons, key=lambda x: x[0])  # genomic order

        start_hit = None
        end_hit   = None
        for rank, (s, e, en) in enumerate(tx_sorted, start=1):
            if abs(s - circ_start) <= slop:
                start_hit = en if en else str(rank)
            if abs(e - circ_end) <= slop:
                end_hit = en if en else str(rank)

        if start_hit is not None and end_hit is not None:
            if start_hit == end_hit:
                return f"e{start_hit}"
            return f"e{start_hit}-e{end_hit}"

    return ""


# ── circRNA coordinate parsing ────────────────────────────────────────────────

def parse_circ_ids(count_matrix_path: str) -> pd.DataFrame:
    """
    Extract circRNA coordinates from count_matrix row names.
    Supported formats: chrom:start|end  or  chrom:start-end  (1-based)
    """
    mat = pd.read_csv(count_matrix_path, sep="\t", index_col=0)
    circ_ids = list(mat.index)

    rows = []
    skipped = 0
    for circ_id in circ_ids:
        try:
            chrom, coord = circ_id.rsplit(":", 1)
            if "|" in coord:
                start_s, end_s = coord.split("|", 1)
            elif "-" in coord:
                start_s, end_s = coord.split("-", 1)
            else:
                skipped += 1
                continue
            rows.append({
                "circ_id": circ_id,
                "chrom":   chrom,
                "start":   int(start_s),
                "end":     int(end_s),
            })
        except (ValueError, IndexError):
            skipped += 1

    if skipped:
        print(f"[assign_isoforms] {skipped} circ_ids could not be parsed", file=sys.stderr)
    return pd.DataFrame(rows)


# ── Host gene assignment ──────────────────────────────────────────────────────

def assign_host_gene(circ_df: pd.DataFrame, gene_df: pd.DataFrame) -> pd.DataFrame:
    """
    For each circRNA, find the smallest gene body that fully contains both
    BSJ endpoints (containment = g_start <= c_start AND g_end >= c_end).
    When multiple genes qualify, the narrowest span wins (most specific).
    """
    # Index genes by chromosome for O(n_genes/n_chr) per lookup
    by_chrom: Dict[str, pd.DataFrame] = {
        chrom: grp.reset_index(drop=True)
        for chrom, grp in gene_df.groupby("chrom")
    }

    results = []
    for _, circ in circ_df.iterrows():
        c_chrom, c_start, c_end = circ["chrom"], circ["start"], circ["end"]
        gene_id = gene_name = "intergenic"

        if c_chrom in by_chrom:
            cands = by_chrom[c_chrom]
            mask  = (cands["g_start"] <= c_start) & (cands["g_end"] >= c_end)
            hits  = cands[mask]
            if not hits.empty:
                hits = hits.copy()
                hits["_span"] = hits["g_end"] - hits["g_start"]
                best      = hits.loc[hits["_span"].idxmin()]
                gene_id   = best["gene_id"]
                gene_name = best["gene_name"]

        strand = "."
        if c_chrom in by_chrom:
            cands = by_chrom[c_chrom]
            mask  = (cands["g_start"] <= c_start) & (cands["g_end"] >= c_end)
            hits  = cands[mask]
            if not hits.empty:
                hits = hits.copy()
                hits["_span"] = hits["g_end"] - hits["g_start"]
                best   = hits.loc[hits["_span"].idxmin()]
                strand = best.get("strand", ".")

        results.append({
            "circ_id":    circ["circ_id"],
            "chrom":      c_chrom,
            "start":      c_start,
            "end":        c_end,
            "strand":     strand,
            "gene_id":    gene_id,
            "gene_name":  gene_name,
            "isoform_id": f"{gene_id}|{c_chrom}:{c_start}-{c_end}",
        })


    df = pd.DataFrame(results)

    # Count isoforms per gene and merge back
    iso_counts = (
        df.groupby("gene_id")["circ_id"]
        .count()
        .rename("n_isoforms")
    )
    df = df.merge(iso_counts, on="gene_id")
    return df


# ── Entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description="Assign circRNA isoforms to host genes using a GTF"
    )
    p.add_argument(
        "--count-matrix", required=True,
        help="BSJ count matrix TSV (row names: chrom:start|end)",
    )
    p.add_argument("--gtf", required=True, help="Genome annotation GTF")
    p.add_argument("--out", required=True, help="Output isoform groups TSV")
    args = p.parse_args()

    circ_df = parse_circ_ids(args.count_matrix)
    print(f"[assign_isoforms] {len(circ_df):,} circRNAs to assign", file=sys.stderr)

    gene_df = parse_gtf_genes(args.gtf)

    result = assign_host_gene(circ_df, gene_df)

    # Add exon-span annotation
    print("[assign_isoforms] Parsing exons for exon-span annotation…", file=sys.stderr)
    gene_exons = parse_gtf_exons(args.gtf)
    result["exon_span"] = result.apply(
        lambda r: _exon_span(r["start"], r["end"], r["gene_id"], gene_exons)
        if r["gene_id"] != "intergenic" else "",
        axis=1,
    )
    n_annotated = (result["exon_span"] != "").sum()
    print(
        f"[assign_isoforms] Exon span annotated: {n_annotated}/{len(result)} circRNAs",
        file=sys.stderr,
    )

    # Region classification
    # exonic: both BSJ boundaries match exon edges (ecircRNA)
    # intronic: within gene body but boundaries don't match exon edges (ciRNA)
    # intergenic: outside any annotated gene
    def _classify(row):
        if row["gene_id"] == "intergenic":
            return "intergenic"
        return "exonic" if row["exon_span"] else "intronic"

    result["region"] = result.apply(_classify, axis=1)
    counts = result["region"].value_counts().to_dict()
    print(
        f"[assign_isoforms] Region: exonic={counts.get('exonic',0)} "
        f"intronic={counts.get('intronic',0)} intergenic={counts.get('intergenic',0)}",
        file=sys.stderr,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)

    n_total = len(result)
    n_multi = int((result["n_isoforms"] > 1).sum())
    n_genes = result[result["n_isoforms"] > 1]["gene_id"].nunique()
    pct     = n_multi / n_total * 100 if n_total > 0 else 0.0
    print(
        f"[assign_isoforms] {n_multi}/{n_total} circRNAs ({pct:.1f}%) "
        f"in {n_genes} multi-isoform genes → {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
