"""
annotate_circbase.py – Annotate pipeline circRNAs against circBase coordinates.

Downloads (or reads a local copy of) the circBase hg19 database and matches
each circRNA from the consensus summary TSV by coordinate overlap within slop.

Output TSV adds columns:
    circbase_id   – circBase identifier (e.g. hsa_circ_0000001), or "novel"
    circbase_gene – host gene name from circBase
    in_circbase   – 1 / 0

Usage:
    python scripts/annotate_circbase.py \
        --summary  results/.../consensus_summary.tsv \
        --output   results/.../circbase_annotated.tsv \
        --slop     10 \
        [--circbase-file /path/to/hsa_hg19_circRNA.txt]
        # omit --circbase-file to auto-download
"""

from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import pandas as pd

CIRCBASE_URL = (
    "http://www.circbase.org/download/hsa_hg19_circRNA.txt"
)
CIRCBASE_FALLBACK_URL = (
    "https://raw.githubusercontent.com/circBase/circBase/master/"
    "data/hsa_hg19_circRNA.txt"
)


# ── circBase loader ───────────────────────────────────────────────────────────

def load_circbase(path: str | None) -> pd.DataFrame:
    """Load circBase TSV; download if path is None or file absent."""
    if path and Path(path).exists():
        print(f"[circBase] loading local file: {path}", file=sys.stderr)
        return _read_circbase(path)

    cache = Path(path) if path else Path("/tmp/circbase_hg19.txt")
    if not cache.exists():
        print(f"[circBase] downloading to {cache} …", file=sys.stderr)
        try:
            urllib.request.urlretrieve(CIRCBASE_URL, str(cache))
        except Exception as e:
            print(f"[circBase] primary URL failed ({e}), trying fallback …", file=sys.stderr)
            try:
                urllib.request.urlretrieve(CIRCBASE_FALLBACK_URL, str(cache))
            except Exception as e2:
                sys.exit(
                    f"[circBase] ERROR: could not download circBase file.\n"
                    f"  Primary:  {CIRCBASE_URL}\n"
                    f"  Fallback: {CIRCBASE_FALLBACK_URL}\n"
                    f"  Last error: {e2}\n"
                    f"  Provide a local copy with --circbase-file."
                )
    return _read_circbase(str(cache))


def _read_circbase(path: str) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        first_line = fh.readline().strip()

    # Strip leading '#' that circBase uses for its header line
    if first_line.startswith("#"):
        first_line = first_line.lstrip("#").strip()
    raw_header = first_line.split("\t")
    header = [h.lower().strip() for h in raw_header]

    # Use header=None + skiprows=1 so the actual first data row is NOT lost.
    # (comment="#" + header=0 would silently drop the first data entry.)
    df = pd.read_csv(
        path, sep="\t",
        header=None, names=header, skiprows=1,
        dtype=str,
        compression="gzip" if path.endswith(".gz") else None,
    )

    # Normalise column names across circBase versions.
    # Handles "# chrom", "circRNA ID", "gene symbol", etc.
    rename_map = {}
    for col in header:
        if "chrom" in col or col == "chr":
            rename_map[col] = "chr"
        elif col in ("start", "chromstart"):
            rename_map[col] = "start"
        elif col in ("end", "chromend"):
            rename_map[col] = "end"
        elif ("name" in col or "circ_id" in col or col == "id"
              or ("circrna" in col and "id" in col)):
            rename_map[col] = "circbase_id"
        elif "gene" in col or "symbol" in col:
            rename_map[col] = "gene"
    df = df.rename(columns=rename_map)

    for c in ("start", "end"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["chr", "start", "end"])
    print(f"[circBase] loaded {len(df):,} entries", file=sys.stderr)
    return df


# ── Annotation ────────────────────────────────────────────────────────────────

def _annotate_row(chrom: str, start: int, end: int,
                  cb: pd.DataFrame, slop: int) -> tuple[str, str]:
    """Return (circbase_id, gene) for the best matching circBase entry, or ('novel', '')."""
    sub = cb[cb["chr"] == chrom]
    if sub.empty:
        return "novel", ""
    dist = sub.apply(
        lambda r: max(abs(r["start"] - start), abs(r["end"] - end)), axis=1
    )
    best_idx = dist.idxmin()
    if dist[best_idx] <= slop:
        row = sub.loc[best_idx]
        cb_id = row.get("circbase_id", "")
        gene  = row.get("gene", "")
        return str(cb_id), str(gene)
    return "novel", ""


def annotate(summary_files: list[str], circbase_df: pd.DataFrame,
             slop: int) -> pd.DataFrame:
    """
    Load one or more consensus_summary.tsv files, deduplicate by circ_id,
    then annotate all unique circRNAs against circBase.
    Accepting multiple files ensures every circRNA in the merged count matrix
    is annotated, not just those from a single sample.
    """
    frames = [pd.read_csv(f, sep="\t") for f in summary_files]
    df = pd.concat(frames, ignore_index=True)

    # Keep the row with the highest confidence_score when duplicated
    if "confidence_score" in df.columns:
        df = (df.sort_values("confidence_score", ascending=False)
                .drop_duplicates(subset="circ_id")
                .reset_index(drop=True))
    else:
        df = df.drop_duplicates(subset="circ_id").reset_index(drop=True)

    print(f"[circBase] {len(df)} unique circRNAs to annotate "
          f"(from {len(summary_files)} sample(s))", file=sys.stderr)

    # Build per-chromosome index for O(n_circ × n_chr) instead of O(n_circ × N)
    cb_by_chr: dict[str, pd.DataFrame] = {
        chrom: grp for chrom, grp in circbase_df.groupby("chr")
    }

    cb_ids, cb_genes, in_cb = [], [], []
    for _, row in df.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end   = int(row["end"])
        sub   = cb_by_chr.get(chrom, pd.DataFrame())

        if sub.empty:
            cid, gene = "novel", ""
        else:
            dist      = sub.apply(
                lambda r: max(abs(r["start"] - start), abs(r["end"] - end)), axis=1
            )
            best_idx  = dist.idxmin()
            if dist[best_idx] <= slop:
                best  = sub.loc[best_idx]
                cid   = str(best.get("circbase_id", ""))
                gene  = str(best.get("gene", ""))
            else:
                cid, gene = "novel", ""

        cb_ids.append(cid)
        cb_genes.append(gene)
        in_cb.append(0 if cid == "novel" else 1)

    df["circbase_id"]   = cb_ids
    df["circbase_gene"] = cb_genes
    df["in_circbase"]   = in_cb

    n_known = sum(in_cb)
    print(
        f"[circBase] annotated {n_known}/{len(df)} circRNAs "
        f"({100*n_known/max(len(df),1):.1f}% known)",
        file=sys.stderr,
    )
    return df


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Annotate circRNAs against circBase (hg19)"
    )
    parser.add_argument("--summary", nargs="+", required=True,
                        help="One or more consensus_summary.tsv files from consensus_filter.py")
    parser.add_argument("--output",        required=True, help="Output annotated TSV")
    parser.add_argument("--slop",          type=int, default=10, help="Coordinate tolerance in bp (default: 10)")
    parser.add_argument("--circbase-file", default=None, dest="circbase_file",
                        help="Local circBase hg19 file (auto-downloaded if absent)")
    args = parser.parse_args()

    cb = load_circbase(args.circbase_file)
    df = annotate(args.summary, cb, args.slop)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[circBase] written → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
