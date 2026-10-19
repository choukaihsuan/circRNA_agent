"""
download_geo.py – Fetch SRR metadata from GEO/SRA.

Priority:
  1. RunInfo.csv (offline, most reliable on HPC)
  2. GSE* → pysradb SRAweb: gse_to_srp → srp_to_srr
  3. PRJNA*/SRP*/SRX* → NCBI eUtils API (no pysradb needed)
  4. Raises with instructions for manual fallback
"""

from __future__ import annotations

import io
import json
import urllib.request
import urllib.error
from pathlib import Path

import pandas as pd


# ── Column name normalisation map ────────────────────────────────────────────
_RUNINFO_COLS = {
    "Run": "srr_id",
    "SampleName": "sample_name",
    "SRX_accession": "gsm_id",
    "Experiment": "gsm_id",     # alternative column name
    "BioSample": "biosample",
    "LibraryLayout": "layout",
    "ScientificName": "organism",
    "LibrarySelection": "library_selection",
    "LibrarySource": "library_source",
    "Disease": "disease_state",
    "Tumor": "tumor_flag",
    "CenterName": "center",
    "SRAStudy": "srp_id",
    "BioProject": "bioproject",
}

_SRAWEB_COLS = {
    "run_accession":      "srr_id",
    "experiment_alias":   "gsm_id",       # GSM accession lives here
    "experiment_title":   "sample_name",
    "library_layout":     "layout",
    "organism_name":      "organism",
    "tissue":             "tissue",
    "disease state":      "disease_state",
    "source_name":        "source_name",
}

# Accession prefixes routed to eUtils (not GEO)
_SRA_PREFIXES = ("PRJNA", "SRP", "SRX", "SRS", "SRR", "ERP", "DRP")


def from_runinfo(csv_path: str) -> pd.DataFrame:
    """Parse a RunInfo.csv downloaded from NCBI SRA Run Selector."""
    df = pd.read_csv(csv_path)

    # Normalise known columns; keep extras
    rename = {k: v for k, v in _RUNINFO_COLS.items() if k in df.columns}
    df = df.rename(columns=rename)

    if "srr_id" not in df.columns:
        raise ValueError(
            f"Could not find SRR IDs in {csv_path}.\n"
            "Expected a 'Run' column (NCBI SRA RunInfo format)."
        )

    df = df[["srr_id"] + [c for c in ["sample_name", "gsm_id", "layout", "organism"] if c in df.columns]]
    print(f"[OK] Loaded {len(df)} samples from {csv_path}")
    return df


def from_sra_accession(accession: str) -> pd.DataFrame:
    """
    Fetch SRR metadata for PRJNA/SRP/SRX accessions via NCBI eUtils API.
    No external dependencies beyond urllib (stdlib).

    Supports: PRJNA* (BioProject), SRP* (SRA Study), SRX* (SRA Experiment)
    """
    acc = accession.strip()
    print(f"[..] Querying NCBI eUtils for {acc} ...")

    # Step 1: esearch → internal SRA ID list
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={acc}&retmax=500&retmode=json"
    )
    try:
        with urllib.request.urlopen(search_url, timeout=30) as r:
            data = json.loads(r.read())
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"NCBI eUtils connection failed: {exc}\n"
            "On HPC without internet access, download RunInfo.csv manually:\n"
            "  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<PRJNA_or_SRP>\n"
            "Then run: python scripts/agent.py --runinfo RunInfo.csv"
        ) from exc

    ids = data["esearchresult"]["idlist"]
    count = int(data["esearchresult"]["count"])
    if not ids:
        raise ValueError(f"No SRA records found for {acc}")
    print(f"[OK] {count} runs found for {acc}")

    # Step 2: efetch → RunInfo CSV (same format as manually-downloaded RunInfo)
    fetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=sra&id={','.join(ids)}&rettype=runinfo&retmode=text"
    )
    with urllib.request.urlopen(fetch_url, timeout=60) as r:
        csv_text = r.read().decode("utf-8")

    df_raw = pd.read_csv(io.StringIO(csv_text))

    # Normalise column names (reuse _RUNINFO_COLS)
    rename = {k: v for k, v in _RUNINFO_COLS.items() if k in df_raw.columns}
    df = df_raw.rename(columns=rename)

    if "srr_id" not in df.columns:
        raise ValueError(
            f"Unexpected RunInfo schema – 'Run' column missing.\n"
            f"Columns available: {list(df_raw.columns)}"
        )

    keep = ["srr_id"] + [
        c for c in ["sample_name", "gsm_id", "layout", "organism",
                    "disease_state", "source_name", "tissue",
                    "library_selection", "library_source", "bioproject", "srp_id"]
        if c in df.columns
    ]
    df = df[keep].drop_duplicates(subset="srr_id")
    print(f"[OK] Retrieved {len(df)} samples for {acc}")

    # Warn if library looks like poly-A selection (bad for circRNA)
    if "library_selection" in df.columns:
        sel = df["library_selection"].dropna().str.upper().unique()
        if any("POLY" in s or s == "CDNA" for s in sel):
            print(
                f"[WARN] LibrarySelection={list(sel)} — may be poly-A selection.\n"
                "       circRNAs lack poly-A tails; detection yield could be very low."
            )

    return df


def from_geo(gse_id: str) -> pd.DataFrame:
    """
    Fetch SRR metadata from SRA via pysradb SRAweb.
    Flow: GSE → SRP (gse_to_srp) → SRR list (srp_to_srr, detailed=True).
    """
    try:
        from pysradb.sraweb import SRAweb  # type: ignore[import-untyped]
    except ImportError:
        raise ImportError("pysradb not installed.  Run: pip install pysradb")

    print(f"[..] Querying SRA for {gse_id} ...")
    try:
        db     = SRAweb()
        srp_df = db.gse_to_srp(gse_id)
        if srp_df is None or srp_df.empty:
            raise ValueError(f"No SRP accession found for {gse_id}")
        srp_id = srp_df["study_accession"].iloc[0]
        print(f"[OK] Resolved {gse_id} → {srp_id}")

        df = db.srp_to_srr(srp_id, detailed=True)
    except Exception as exc:
        raise RuntimeError(
            f"pysradb query failed for {gse_id}: {exc}\n"
            "On HPC without internet access, download RunInfo.csv manually:\n"
            "  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=<GSE>&o=acc_s%3Aa\n"
            "Then run: python scripts/agent.py --runinfo RunInfo.csv"
        ) from exc

    if df is None or df.empty:
        raise ValueError(f"No SRR records found under {srp_id}")

    rename = {k: v for k, v in _SRAWEB_COLS.items() if k in df.columns}
    df = df.rename(columns=rename)

    if "srr_id" not in df.columns:
        raise ValueError(
            f"Unexpected pysradb schema – 'run_accession' column missing.\n"
            f"Available columns: {list(df.columns)}"
        )

    keep = ["srr_id"] + [
        c for c in ["sample_name", "gsm_id", "layout", "organism",
                    "tissue", "disease_state", "source_name"]
        if c in df.columns
    ]
    df = df[keep].drop_duplicates(subset="srr_id")
    print(f"[OK] Found {len(df)} samples for {gse_id} ({srp_id})")
    return df


def preflight_accession(accession: str) -> dict:
    """Check a dataset accession's existence and scale WITHOUT any side
    effects (no file writes, no directory creation, no config/queue state) --
    safe to call before any job is created. Reuses exactly the same network
    calls prepare_metadata() would make (pysradb for GSE*, NCBI eUtils for
    PRJNA*/SRP*/SRX*/SRS*/SRR*/ERP*/DRP*), just reading fields
    (library_layout, total/public size, total_spots/run_total_bases for GSE;
    LibraryLayout/size_MB/avgLength for the eUtils RunInfo path) that
    prepare_metadata()'s own column-narrowing step discards -- not issuing an
    extra query beyond what the real download would already make.

    Returns a dict, always (never raises): {"exists", "n_samples",
    "total_bytes", "layout", "avg_read_length_bp", "error"}. "exists": False
    with a populated "error" covers both "genuinely doesn't exist" and "query
    failed for some other reason" (network/API trouble, unexpected schema) --
    callers must treat both as "cannot proceed," never silently default to
    proceeding on an inconclusive answer."""
    acc = (accession or "").strip().upper()
    empty = {"exists": False, "n_samples": None, "total_bytes": None,
              "layout": None, "avg_read_length_bp": None, "error": None}
    try:
        if acc.startswith("GSE"):
            from pysradb.sraweb import SRAweb  # type: ignore[import-untyped]
            db = SRAweb()
            srp_df = db.gse_to_srp(acc)
            if srp_df is None or srp_df.empty:
                return {**empty, "error": f"No SRP study found for {acc}"}
            srp_id = srp_df["study_accession"].iloc[0]
            df = db.srp_to_srr(srp_id, detailed=True)
            if df is None or df.empty:
                return {**empty, "error": f"No runs found under {srp_id}"}

            layout = None
            if "library_layout" in df.columns:
                vals = df["library_layout"].dropna().unique()
                layout = vals[0] if len(vals) == 1 else "/".join(sorted(set(str(v) for v in vals)))

            total_bytes = None
            if "total_size" in df.columns:
                sizes = pd.to_numeric(df["total_size"], errors="coerce").dropna()
                if len(sizes):
                    total_bytes = int(sizes.sum())

            avg_len = None
            if "run_total_bases" in df.columns and "total_spots" in df.columns:
                bases = pd.to_numeric(df["run_total_bases"], errors="coerce")
                spots = pd.to_numeric(df["total_spots"], errors="coerce")
                ratio = (bases / spots).dropna()
                if len(ratio):
                    divisor = 2 if layout and "PAIRED" in str(layout).upper() else 1
                    avg_len = int(round(ratio.mean() / divisor))

            return {"exists": True, "n_samples": len(df), "total_bytes": total_bytes,
                    "layout": layout, "avg_read_length_bp": avg_len, "error": None}

        elif any(acc.startswith(p) for p in _SRA_PREFIXES):
            search_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                f"?db=sra&term={acc}&retmax=500&retmode=json"
            )
            with urllib.request.urlopen(search_url, timeout=15) as r:
                data = json.loads(r.read())
            ids = data["esearchresult"]["idlist"]
            if not ids:
                return {**empty, "error": f"No SRA records found for {acc}"}

            fetch_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=sra&id={','.join(ids)}&rettype=runinfo&retmode=text"
            )
            with urllib.request.urlopen(fetch_url, timeout=30) as r:
                csv_text = r.read().decode("utf-8")
            df = pd.read_csv(io.StringIO(csv_text))
            if df.empty or "Run" not in df.columns:
                return {**empty, "error": f"Unexpected or empty RunInfo for {acc}"}

            layout = None
            if "LibraryLayout" in df.columns:
                vals = df["LibraryLayout"].dropna().unique()
                layout = vals[0] if len(vals) == 1 else "/".join(sorted(set(str(v) for v in vals)))

            total_bytes = None
            if "size_MB" in df.columns:
                mb = pd.to_numeric(df["size_MB"], errors="coerce").dropna()
                if len(mb):
                    total_bytes = int(mb.sum() * 1024 * 1024)

            avg_len = None
            if "avgLength" in df.columns:
                lens = pd.to_numeric(df["avgLength"], errors="coerce").dropna()
                if len(lens):
                    divisor = 2 if layout and "PAIRED" in str(layout).upper() else 1
                    avg_len = int(round(lens.mean() / divisor))

            return {"exists": True, "n_samples": len(df), "total_bytes": total_bytes,
                    "layout": layout, "avg_read_length_bp": avg_len, "error": None}

        else:
            return {**empty, "error": f"Unrecognised accession format: {acc!r}"}
    except Exception as exc:
        return {**empty, "error": f"{type(exc).__name__}: {exc}"}


def prepare_metadata(
    gse_id: str | None = None,
    runinfo_csv: str | None = None,
    output_file: str = "metadata/library_info.csv",
) -> pd.DataFrame:
    """
    Resolve SRR metadata from any source and save to output_file.

    Accession routing:
      GSE*              → GEO via pysradb
      PRJNA*/SRP*/SRX*  → NCBI eUtils API (no pysradb)
      runinfo_csv       → local CSV file
    """
    if runinfo_csv:
        df = from_runinfo(runinfo_csv)
    elif gse_id:
        acc = gse_id.strip().upper()
        if acc.startswith("GSE"):
            df = from_geo(gse_id)
        elif any(acc.startswith(p) for p in _SRA_PREFIXES):
            df = from_sra_accession(gse_id)
        else:
            raise ValueError(
                f"Unrecognised accession format: {gse_id!r}\n"
                "Expected: GSE* (GEO), PRJNA* (BioProject), SRP* (SRA Study), "
                "SRX* (SRA Experiment), or --runinfo RunInfo.csv"
            )
    else:
        raise ValueError("Provide either gse_id/accession or runinfo_csv")

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)
    print(f"[OK] Metadata saved → {output_file}")
    return df


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Fetch SRR metadata from GEO or RunInfo.csv")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gse",      help="GEO accession (e.g. GSE113230)")
    group.add_argument("--runinfo",  help="Path to NCBI SRA RunInfo.csv")
    parser.add_argument("--out", default="metadata/library_info.csv", help="Output CSV path")
    args = parser.parse_args()

    prepare_metadata(gse_id=args.gse, runinfo_csv=args.runinfo, output_file=args.out)

    # Quick test:
    #   python scripts/download_geo.py --gse PRJNA808398
    #   python scripts/download_geo.py --gse SRP156355
    #   python scripts/download_geo.py --gse GSE113230
