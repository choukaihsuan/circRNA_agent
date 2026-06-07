"""
prepare_metadata.py – Assign case/control group labels to samples.

Modes:
  1. Auto-detect: match common keywords in sample_name
  2. Interactive: prompt user via CLI
  3. CSV import: load pre-made sample_groups.csv
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

# Keywords for auto-detection (case-insensitive).
# Matches are mapped to the configured case/control labels (default: tumor/normal).
_CASE_KEYWORDS    = ["tumor", "cancer", "carcinoma", "malignant", "case", "patient",
                     "hypoxia", "treated", "treatment", "drug", "knockdown", "ko",
                     "overexpression", "oe", "stimulated"]
_CONTROL_KEYWORDS = ["normal", "control", "ctrl", "healthy", "adjacent", "non-tumor",
                     "benign", "normoxia", "untreated", "vehicle", "wildtype", "wt",
                     "scramble", "mock"]


def _detect_condition(name: str, case_label: str = "tumor",
                      control_label: str = "normal") -> str | None:
    """Return case_label or control_label if a keyword matches, else None."""
    name_lower = name.lower()
    for kw in _CASE_KEYWORDS:
        if re.search(kw, name_lower):
            return case_label
    for kw in _CONTROL_KEYWORDS:
        if re.search(kw, name_lower):
            return control_label
    return None


def auto_assign(metadata: pd.DataFrame,
                case_label: str = "tumor",
                control_label: str = "normal") -> pd.DataFrame:
    """
    Try to assign conditions from the 'sample_name' column.
    Returns a DataFrame with columns [srr_id, sample_name, condition].
    Rows where condition cannot be inferred are left as empty string.
    """
    if "sample_name" not in metadata.columns:
        raise ValueError("metadata must contain a 'sample_name' column for auto-detection")

    rows = []
    for _, row in metadata.iterrows():
        condition = _detect_condition(
            str(row.get("sample_name", "")), case_label, control_label
        )
        rows.append({
            "srr_id":      row["srr_id"],
            "sample_name": row.get("sample_name", ""),
            "condition":   condition or "",
        })

    df = pd.DataFrame(rows)
    detected   = (df["condition"] != "").sum()
    undetected = (df["condition"] == "").sum()
    print(f"[OK] Auto-assigned: {detected} samples | unresolved: {undetected}")
    return df


def interactive_assign(metadata: pd.DataFrame,
                       case_label: str = "tumor",
                       control_label: str = "normal") -> pd.DataFrame:
    """
    Walk through each sample and ask user to assign a condition.
    Pre-fills auto-detected values so user can press Enter to accept.
    Accepts any non-empty string (not limited to tumor/normal).
    """
    groups = auto_assign(metadata, case_label, control_label)

    print("\n── Manual group assignment ──────────────────────────")
    print("Press Enter to accept the auto-detected value, or type a new one.")
    print(f"Typical values: {case_label} | {control_label}\n")

    updated = []
    for _, row in groups.iterrows():
        suggestion = row["condition"] or "?"
        user_input = input(
            f"  {row['srr_id']}  [{row['sample_name']}]  condition [{suggestion}]: "
        ).strip()
        condition = user_input if user_input else suggestion
        if not condition or condition == "?":
            condition = "unknown"
        updated.append({
            "srr_id":      row["srr_id"],
            "sample_name": row["sample_name"],
            "condition":   condition,
        })

    return pd.DataFrame(updated)


def load_or_create(
    metadata_file:  str  = "metadata/library_info.csv",
    groups_file:    str  = "metadata/sample_groups.csv",
    interactive:    bool = False,
    case_label:     str  = "tumor",
    control_label:  str  = "normal",
) -> pd.DataFrame:
    """
    Return sample groups DataFrame.
    If groups_file already exists, load it.
    Otherwise attempt auto-detection, with optional interactive override.
    """
    if Path(groups_file).exists():
        df = pd.read_csv(groups_file)
        print(f"[OK] Loaded existing groups from {groups_file}")
        return df

    metadata = pd.read_csv(metadata_file)
    groups   = (interactive_assign(metadata, case_label, control_label)
                if interactive
                else auto_assign(metadata, case_label, control_label))

    Path(groups_file).parent.mkdir(parents=True, exist_ok=True)
    groups.to_csv(groups_file, index=False)
    print(f"[OK] Groups saved → {groups_file}")
    return groups


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Assign case/control labels to samples")
    parser.add_argument("--metadata",      default="metadata/library_info.csv")
    parser.add_argument("--groups",        default="metadata/sample_groups.csv")
    parser.add_argument("--interactive",   action="store_true",
                        help="Prompt user for each sample's condition")
    parser.add_argument("--case-label",    default="tumor",
                        help="Label for the case/treatment group (default: tumor)")
    parser.add_argument("--control-label", default="normal",
                        help="Label for the control group (default: normal)")
    args = parser.parse_args()

    df = load_or_create(args.metadata, args.groups, args.interactive,
                        args.case_label, args.control_label)
    print(df.to_string(index=False))
