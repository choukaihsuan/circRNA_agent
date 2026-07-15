"""
validate_condition_labels.py – QC for automatic tumor/normal label assignment.

Re-runs the condition auto-detection (`prepare_metadata._detect_condition`) on
every dataset's sample metadata and compares the result against the condition
label already stored in `sample_groups.csv` (the label the pipeline actually
used). The point is to surface samples where auto-detection *disagrees* with the
stored label — those are the rows a human must eyeball, because a mis-labelled
tumor/normal sample silently flips the direction of every downstream DE result.

This is meant to be run on the server, where all datasets live under
`metadata/{DATASET}/{library_info.csv, sample_groups.csv}`.

Verdicts per sample:
  MATCH       detected == stored label                → fine
  MISMATCH    detected != stored, detected not None    → ⚠ inspect (possible flip)
  UNDETECTED  detection returned None                  → auto couldn't tell;
                                                          pipeline relied on the
                                                          stored/manual label (ok,
                                                          but not independently
                                                          confirmed)

Usage:
    python scripts/validate_condition_labels.py \
        [--metadata-root metadata] [--out condition_qc.tsv]

Exit status is 0 regardless of findings (this is a report, not a gate); the
count of MISMATCH rows is printed to stderr so it can be grepped in CI/logs.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Import the *live* detection functions from prepare_metadata so this QC always
# tests exactly what the pipeline runs (no copy-paste drift).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_metadata import (  # noqa: E402
    _detect_condition,
    _match_case_kw,
    _match_control_kw,
)

# Descriptive columns (in priority order) whose text carries the tumor/normal
# signal. sample_groups.csv usually only has sample_name; library_info.csv adds
# tissue / source_name / disease_state, which for many GEO datasets is where the
# "adjacent normal" / "tumor tissue" wording actually lives.
_TEXT_COLS = ["sample_name", "description", "disease_state", "source_name", "tissue", "title"]


def _labels_from_conditions(conds: list[str]) -> tuple[str, str] | None:
    """
    Given the distinct condition strings stored for a dataset, decide which is
    the case label and which is the control label, so detection can be run with
    the dataset's own vocabulary (tumor/normal, treated/vehicle, …).

    Returns (case_label, control_label), or None if it can't be determined
    (e.g. two treatment codes like EPZ6438/DMSO with no keyword cue — those
    samples are simply not auto-detectable and will show up as UNDETECTED).
    """
    labels = sorted({c.strip() for c in conds if str(c).strip()})
    if len(labels) != 2:
        return None
    a, b = labels
    a_case, b_case = _match_case_kw(a) is not None, _match_case_kw(b) is not None
    a_ctrl, b_ctrl = _match_control_kw(a) is not None, _match_control_kw(b) is not None
    if a_case and not b_case:
        return a, b
    if b_case and not a_case:
        return b, a
    if a_ctrl and not b_ctrl:
        return b, a
    if b_ctrl and not a_ctrl:
        return a, b
    return None  # ambiguous — treatment codes etc.


def _build_text(row: pd.Series) -> str:
    parts = []
    for c in _TEXT_COLS:
        v = row.get(c)
        if v is not None and str(v).strip() and str(v).lower() != "nan":
            parts.append(str(v))
    return " ".join(parts)


def _load_dataset(ds_dir: Path) -> pd.DataFrame | None:
    """Merge sample_groups (stored condition) with library_info (descriptive text)."""
    sg_path = ds_dir / "sample_groups.csv"
    if not sg_path.exists():
        return None
    sg = pd.read_csv(sg_path)
    if "condition" not in sg.columns or "srr_id" not in sg.columns:
        return None

    li_path = ds_dir / "library_info.csv"
    if li_path.exists():
        li = pd.read_csv(li_path)
        extra = [c for c in _TEXT_COLS if c in li.columns and c not in sg.columns]
        if "srr_id" in li.columns and extra:
            sg = sg.merge(li[["srr_id", *extra]], on="srr_id", how="left")
    return sg


def validate_dataset(name: str, df: pd.DataFrame) -> pd.DataFrame:
    conds = df["condition"].dropna().astype(str).tolist()
    labels = _labels_from_conditions(conds)
    if labels is None:
        case_label, control_label = "tumor", "normal"  # best-effort defaults
        label_note = "defaulted(tumor/normal)"
    else:
        case_label, control_label = labels
        label_note = f"{case_label}/{control_label}"

    rows = []
    for _, r in df.iterrows():
        text = _build_text(r)
        stored = str(r.get("condition", "")).strip()
        detected = _detect_condition(text, case_label, control_label)
        if detected is None:
            verdict = "UNDETECTED"
        elif stored and detected == stored:
            verdict = "MATCH"
        else:
            verdict = "MISMATCH"
        rows.append({
            "dataset":    name,
            "srr_id":     r.get("srr_id", ""),
            "text":       text[:80],
            "stored":     stored,
            "detected":   detected if detected is not None else "",
            "verdict":    verdict,
            "labels":     label_note,
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="QC tumor/normal auto-detection vs stored labels")
    ap.add_argument("--metadata-root", default="metadata",
                    help="Root dir holding per-dataset subdirs (default: metadata)")
    ap.add_argument("--out", default=None, help="Optional TSV of every sample's verdict")
    args = ap.parse_args()

    root = Path(args.metadata_root)
    if not root.is_dir():
        sys.exit(f"[qc] metadata root not found: {root}")

    # Each subdir with a sample_groups.csv is a dataset. Also accept a
    # sample_groups.csv directly at the root (single-project layout).
    ds_dirs = sorted(d for d in root.iterdir() if d.is_dir() and (d / "sample_groups.csv").exists())
    if (root / "sample_groups.csv").exists():
        ds_dirs.append(root)

    if not ds_dirs:
        sys.exit(f"[qc] no datasets with sample_groups.csv under {root}")

    all_rows = []
    print(f"{'dataset':<16} {'n':>4} {'match':>6} {'mismatch':>9} {'undet':>6}  labels", file=sys.stderr)
    print("-" * 70, file=sys.stderr)
    total_mismatch = 0
    for d in ds_dirs:
        df = _load_dataset(d)
        if df is None or df.empty:
            continue
        res = validate_dataset(d.name if d != root else root.name, df)
        all_rows.append(res)
        vc = res["verdict"].value_counts()
        n_match = int(vc.get("MATCH", 0))
        n_mis   = int(vc.get("MISMATCH", 0))
        n_und   = int(vc.get("UNDETECTED", 0))
        total_mismatch += n_mis
        flag = " ⚠" if n_mis else ""
        print(f"{res['dataset'].iloc[0]:<16} {len(res):>4} {n_match:>6} {n_mis:>9} {n_und:>6}  "
              f"{res['labels'].iloc[0]}{flag}", file=sys.stderr)

    if not all_rows:
        sys.exit("[qc] no usable datasets")

    full = pd.concat(all_rows, ignore_index=True)

    mism = full[full["verdict"] == "MISMATCH"]
    if not mism.empty:
        print("\n[qc] ⚠ MISMATCH rows (auto-detection disagrees with stored label — inspect):",
              file=sys.stderr)
        with pd.option_context("display.max_colwidth", 60, "display.width", 200):
            print(mism[["dataset", "srr_id", "text", "stored", "detected"]].to_string(index=False),
                  file=sys.stderr)
    else:
        print("\n[qc] ✅ No MISMATCH rows: auto-detection agrees with every stored label "
              "(where detectable).", file=sys.stderr)

    if args.out:
        full.to_csv(args.out, sep="\t", index=False)
        print(f"\n[qc] full per-sample report → {args.out}", file=sys.stderr)

    print(f"\n[qc] total MISMATCH samples across all datasets: {total_mismatch}", file=sys.stderr)


if __name__ == "__main__":
    main()
