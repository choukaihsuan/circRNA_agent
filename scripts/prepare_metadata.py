"""
prepare_metadata.py – Assign case/control group labels to samples.

Modes:
  1. Auto-detect: suffix pattern (T/N/APN) then keyword matching in sample_name
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


def _normalize_label(text: str) -> str:
    """Clean a raw metadata value into a compact label (lowercase, underscores)."""
    text = re.sub(r"[^a-zA-Z0-9\s_-]", "", text.strip().lower())
    text = re.sub(r"\s+", "_", text)
    return text[:30]  # cap at 30 chars


def _match_case_kw(text: str) -> str | None:
    """Return the matched case keyword (the actual keyword string) or None."""
    text_lower = text.lower()
    for kw in _CASE_KEYWORDS:
        if re.search(r"\b" + kw + r"\b", text_lower):
            return kw
    return None


def _match_control_kw(text: str) -> str | None:
    """Return the matched control keyword (the actual keyword string) or None."""
    text_lower = text.lower()
    for kw in _CONTROL_KEYWORDS:
        if re.search(r"\b" + kw + r"\b", text_lower):
            return kw
    return None


def detect_group_labels(df: pd.DataFrame) -> tuple[str, str]:
    """
    Infer case/control labels from GEO metadata columns.

    Strategy (in priority order):
    1. Structured fields (disease_state, source_name, tissue): if exactly 2
       distinct non-empty values exist, classify each as case/control and
       use the actual value text as the label.
    2. Keyword extraction from all text fields: find the dominant case/control
       keyword across all samples and return those keyword strings.
    3. Fallback: ("tumor", "normal").

    Returns (case_label, control_label).
    """
    # ── Strategy 1: structured fields ────────────────────────────────────────
    for col in ("disease_state", "source_name", "tissue"):
        if col not in df.columns:
            continue
        vals = df[col].dropna().astype(str)
        vals = vals[vals.str.strip() != ""]
        unique_vals = vals.str.strip().str.lower().unique()
        if len(unique_vals) != 2:
            continue

        # Classify each unique value as case or control
        label_a, label_b = unique_vals[0], unique_vals[1]
        a_is_case    = _match_case_kw(label_a) is not None
        b_is_case    = _match_case_kw(label_b) is not None
        a_is_control = _match_control_kw(label_a) is not None
        b_is_control = _match_control_kw(label_b) is not None

        if a_is_case and b_is_control:
            return _normalize_label(label_a), _normalize_label(label_b)
        if b_is_case and a_is_control:
            return _normalize_label(label_b), _normalize_label(label_a)
        # Ambiguous but exactly 2 values: still useful — pick shorter as label
        if not (a_is_case or a_is_control or b_is_case or b_is_control):
            # Can't classify by keyword; skip this field
            continue
        # One side matches, other doesn't — treat matched side as case
        if a_is_case or (not b_is_case and a_is_control is False):
            return _normalize_label(label_a), _normalize_label(label_b)
        return _normalize_label(label_b), _normalize_label(label_a)

    # ── Strategy 2: keyword extraction across all text fields ────────────────
    scan_cols = [c for c in ("disease_state", "source_name", "tissue", "sample_name")
                 if c in df.columns]
    case_kw_counts:    dict[str, int] = {}
    control_kw_counts: dict[str, int] = {}

    for _, row in df.iterrows():
        combined = " ".join(str(row.get(c, "")) for c in scan_cols)
        kw = _match_case_kw(combined)
        if kw:
            case_kw_counts[kw] = case_kw_counts.get(kw, 0) + 1
        kw = _match_control_kw(combined)
        if kw:
            control_kw_counts[kw] = control_kw_counts.get(kw, 0) + 1

    case_label    = max(case_kw_counts,    key=case_kw_counts.__getitem__)    if case_kw_counts    else "tumor"
    control_label = max(control_kw_counts, key=control_kw_counts.__getitem__) if control_kw_counts else "normal"

    if case_label == control_label:
        return "tumor", "normal"

    return case_label, control_label


def _detect_condition_by_suffix(name: str, case_label: str = "tumor",
                                 control_label: str = "normal") -> str | None:
    """
    Detect tumor/normal from SRA SampleName trailing suffix conventions.

    Patterns handled (case-insensitive, must appear at end of name):
      Tumor  : …T,  …_T,  …-T
      Normal : …N,  …_N,  …-N,  …APN,  …CN,  …WT (wildtype control)

    The suffix token must be anchored to a DIGIT (the SRA sample-code
    convention, e.g. M269T, 138T, 87APN, 6WT). This prevents free-text
    metadata such as "Bioengineered construct" or "adjacent normal tissue"
    from matching a bare trailing "t"/"n" and being mis-classified.

    Control tokens (WT/APN/CN/N) are checked BEFORE the tumor `T` rule so that
    "6WT" (wildtype control) is not mis-read as tumor via the trailing "T".
    """
    s = name.strip()
    # Explicit control suffixes first, each anchored to a preceding digit.
    if re.search(r'(?<=[0-9])[-_]?(WT|APN|CN|N)$', s, re.IGNORECASE):
        return control_label
    if re.search(r'(?<=[0-9])[-_]?T$', s, re.IGNORECASE):
        return case_label
    return None


def _detect_condition(name: str, case_label: str = "tumor",
                      control_label: str = "normal") -> str | None:
    """Return case_label or control_label — tries suffix pattern first, then keywords."""
    # Priority 1: SRA-style T/N/APN suffix (e.g. M269T, M269N, 138T, 87APN)
    result = _detect_condition_by_suffix(name, case_label, control_label)
    if result is not None:
        return result
    # Priority 2: keyword search (GEO-style names).
    # Use word-boundary matching (via the shared helpers) so short keywords like
    # "oe"/"ko" don't match inside unrelated words ("Bioengineered", "Tokyo").
    case_kw    = _match_case_kw(name)
    control_kw = _match_control_kw(name)
    if case_kw and control_kw:
        # Both cue types present, e.g. "tumor-adjacent normal tissue".
        # "adjacent"/"non-tumor"/"para-tumor" explicitly mark a normal control
        # sample even though the word "tumor" appears — control wins here.
        if re.search(r'\b(adjacent|non[-_ ]?tumou?r|para[-_ ]?tumou?r)\b',
                     name.lower()):
            return control_label
        # Genuinely ambiguous: leave blank so a human confirms, rather than
        # silently picking the wrong side (old code always returned case).
        return None
    if case_kw:
        return case_label
    if control_kw:
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
