"""
cross_dataset_biomarkers.py – Cross-dataset reproducibility of biomarker candidates.

Motivation: every dimension of the biomarker composite score (rank_biomarkers.py)
comes from one pipeline run, so a top rank is not *independently* corroborated.
The cheapest orthogonal evidence available is reproducibility across independent
cohorts: a circRNA that ranks highly as a biomarker in several independent
datasets of the same cancer type is validated by data no single run can provide.

This script pools multiple `biomarker_candidates.tsv` files, matches circRNAs
across them, and reports how many datasets each circRNA recurs in (with its
rank/score in each). It makes no statistical claim — it is a recurrence table,
the raw material for a "recurrent across N cohorts" figure/statement.

Matching key:
  * Known circRNAs  → circBase id (hsa_circ_XXXXXXX): stable across cohorts,
    assigned by annotate_circbase with slop tolerance — the correct cross-cohort
    identifier.
  * Novel circRNAs  → genomic coordinates clustered within --slop bp across the
    pooled set, so the same junction detected at slightly different coordinates
    in different cohorts still matches.

Usage:
  # explicit files (name inferred from the *_results/ path component):
  python scripts/cross_dataset_biomarkers.py \
      --inputs ~/GSE113230_results/de/biomarker_candidates.tsv \
               ~/GSE133998_results/de/biomarker_candidates.tsv \
               ~/SRP156355_results/de/biomarker_candidates.tsv \
      --out breast_biomarker_recurrence.tsv

  # or a results root + dataset list:
  python scripts/cross_dataset_biomarkers.py \
      --results-root ~ --datasets GSE113230 GSE133998 SRP156355 GSE58135 \
      --out breast_biomarker_recurrence.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

_COORD_RE = re.compile(r"^(chr[^:]+):(\d+)[|\-](\d+)$")


def _infer_name(path: Path) -> str:
    """Pull a dataset name from a .../{NAME}_results/de/biomarker_candidates.tsv path."""
    for part in path.parts:
        m = re.match(r"([A-Za-z]+\d+)_results", part)
        if m:
            return m.group(1)
    # fallback: grandparent dir name (…/de/file → grandparent)
    return path.parent.parent.name or path.stem


def _parse_coord(circ_id: str) -> tuple[str, int, int] | None:
    m = _COORD_RE.match(str(circ_id).strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def _valid_circbase(v) -> bool:
    s = str(v).strip().lower()
    return bool(s) and s not in ("novel", "nan", "none", "na", "")


def _cluster_novel(coords: list[tuple[str, int, int]], slop: int) -> dict[tuple[str, int, int], str]:
    """Greedy same-chromosome clustering of novel coordinates → {coord: cluster_key}."""
    by_chr: dict[str, list[tuple[str, int, int]]] = {}
    for c in coords:
        by_chr.setdefault(c[0], []).append(c)
    mapping: dict[tuple[str, int, int], str] = {}
    for chrom, members in by_chr.items():
        members = sorted(set(members), key=lambda c: (c[1], c[2]))
        clusters: list[list[tuple[str, int, int]]] = []
        for c in members:
            for cl in clusters:
                a = cl[0]
                if max(abs(c[1] - a[1]), abs(c[2] - a[2])) <= slop:
                    cl.append(c)
                    break
            else:
                clusters.append([c])
        for cl in clusters:
            key = f"novel:{chrom}:{cl[0][1]}-{cl[0][2]}"  # anchor coord names the cluster
            for c in cl:
                mapping[c] = key
    return mapping


def load_all(inputs: list[tuple[str, Path]], slop: int, top_n: int) -> pd.DataFrame:
    frames = []
    for name, path in inputs:
        if not path.exists():
            print(f"[recur] skip (missing): {path}", file=sys.stderr)
            continue
        df = pd.read_csv(path, sep="\t")
        if "circ_id" not in df.columns:
            print(f"[recur] skip (no circ_id): {path}", file=sys.stderr)
            continue
        if "rank" in df.columns:
            df = df.sort_values("rank")
        if top_n and top_n > 0:
            df = df.head(top_n)
        df = df.copy()
        df["dataset"] = name
        frames.append(df)
        print(f"[recur] {name}: {len(df)} candidates from {path}", file=sys.stderr)
    if not frames:
        sys.exit("[recur] no usable biomarker_candidates.tsv inputs")

    pool = pd.concat(frames, ignore_index=True)

    # Build the cross-dataset match key.
    novel_coords = []
    for _, r in pool.iterrows():
        if not _valid_circbase(r.get("circbase_id")):
            c = _parse_coord(r["circ_id"])
            if c:
                novel_coords.append(c)
    novel_map = _cluster_novel(novel_coords, slop)

    def _key(r) -> str:
        if _valid_circbase(r.get("circbase_id")):
            return f"circbase:{str(r['circbase_id']).strip()}"
        c = _parse_coord(r["circ_id"])
        if c and c in novel_map:
            return novel_map[c]
        return f"raw:{r['circ_id']}"  # unparseable coord: match only on exact id

    pool["match_key"] = pool.apply(_key, axis=1)
    return pool


def summarise(pool: pd.DataFrame) -> pd.DataFrame:
    def _gene(g: pd.DataFrame) -> str:
        _skip = ("nan", "intergenic", "", "=", "-", "novel", "none", "na")
        for col in ("circbase_gene", "gene_name"):
            if col in g.columns:
                vals = [str(v) for v in g[col] if str(v).strip().lower() not in _skip]
                if vals:
                    return max(set(vals), key=vals.count)
        return ""

    def _circbase(g: pd.DataFrame) -> str:
        if "circbase_id" in g.columns:
            vals = [str(v).strip() for v in g["circbase_id"] if _valid_circbase(v)]
            if vals:
                return max(set(vals), key=vals.count)
        return ""

    has_score = "biomarker_score" in pool.columns
    has_rank  = "rank" in pool.columns
    has_lfc   = "log2FC" in pool.columns

    rows = []
    for key, g in pool.groupby("match_key"):
        datasets = sorted(g["dataset"].unique())
        row = {
            "match_key":   key,
            "circbase_id": _circbase(g),
            "gene":        _gene(g),
            "n_datasets":  len(datasets),
            "datasets":    ",".join(datasets),
        }
        if has_score:
            row["mean_score"] = round(g["biomarker_score"].mean(), 4)
            row["best_score"] = round(g["biomarker_score"].max(), 4)
        if has_rank:
            row["best_rank"] = int(g["rank"].min())
            # per-dataset rank, e.g. "GSE113230:3;GSE133998:11"
            row["ranks"] = ";".join(
                f"{d}:{int(g.loc[g['dataset'] == d, 'rank'].min())}" for d in datasets
            )
        if has_lfc:
            # direction consistency: do all cohorts agree on up/down?
            signs = {"+" if v > 0 else "-" for v in g["log2FC"] if pd.notna(v)}
            row["lfc_dir"] = "".join(sorted(signs)) if signs else ""
            row["dir_consistent"] = int(len(signs) == 1)
        rows.append(row)

    out = pd.DataFrame(rows)
    sort_cols = ["n_datasets"] + (["mean_score"] if has_score else [])
    out = out.sort_values(sort_cols, ascending=[False] * len(sort_cols)).reset_index(drop=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Cross-dataset biomarker recurrence")
    ap.add_argument("--inputs", nargs="*", default=[],
                    help="biomarker_candidates.tsv paths (name inferred from *_results/ path); "
                         "or NAME=path pairs to set the name explicitly")
    ap.add_argument("--results-root", default=None,
                    help="Root holding {DATASET}_results/de/biomarker_candidates.tsv")
    ap.add_argument("--datasets", nargs="*", default=[],
                    help="Dataset names to combine with --results-root")
    ap.add_argument("--slop", type=int, default=10,
                    help="bp tolerance when matching novel circRNAs across datasets (default 10)")
    ap.add_argument("--top-n", type=int, default=0,
                    help="Only consider each dataset's top-N ranked candidates (0 = all, default)")
    ap.add_argument("--min-datasets", type=int, default=2,
                    help="Recurrence threshold to flag/report (default 2)")
    ap.add_argument("--out", default=None, help="Output recurrence TSV")
    args = ap.parse_args()

    inputs: list[tuple[str, Path]] = []
    for item in args.inputs:
        if "=" in item and not item.startswith("/") and "/" not in item.split("=")[0]:
            name, p = item.split("=", 1)
            inputs.append((name, Path(p).expanduser()))
        else:
            p = Path(item).expanduser()
            inputs.append((_infer_name(p), p))
    if args.results_root:
        root = Path(args.results_root).expanduser()
        for ds in args.datasets:
            inputs.append((ds, root / f"{ds}_results" / "de" / "biomarker_candidates.tsv"))

    if not inputs:
        sys.exit("[recur] provide --inputs or (--results-root + --datasets)")

    pool = load_all(inputs, args.slop, args.top_n)
    summary = summarise(pool)

    n_ds = pool["dataset"].nunique()
    recurrent = summary[summary["n_datasets"] >= args.min_datasets]
    print(f"\n[recur] pooled {len(summary)} distinct circRNAs across {n_ds} dataset(s)",
          file=sys.stderr)
    print(f"[recur] recurrent in ≥{args.min_datasets} datasets: {len(recurrent)}",
          file=sys.stderr)

    if not recurrent.empty:
        show = recurrent.head(30)
        cols = [c for c in ["circbase_id", "gene", "n_datasets", "datasets",
                            "mean_score", "best_rank", "dir_consistent"] if c in show.columns]
        print("\n[recur] top recurrent biomarker candidates:", file=sys.stderr)
        with pd.option_context("display.max_colwidth", 40, "display.width", 200):
            print(show[cols].to_string(index=False), file=sys.stderr)

    if args.out:
        summary.to_csv(args.out, sep="\t", index=False)
        print(f"\n[recur] full recurrence table → {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
