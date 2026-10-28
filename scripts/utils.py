"""
utils.py – Shared helpers for the circRNA agent.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import yaml


# ── Config I/O ───────────────────────────────────────────────────────────────

def load_config(path: str = "config.yaml") -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def save_config(cfg: dict, path: str = "config.yaml") -> None:
    with open(path, "w", encoding="utf-8") as fh:
        yaml.dump(cfg, fh, default_flow_style=False, allow_unicode=True)


# ── CIRIquant config generation ──────────────────────────────────────────────

def generate_ciriquant_config(cfg: dict, output: str = "config/ciriquant.yaml") -> None:
    """
    Write a CIRIquant-compatible YAML from the main config's genome section.
    CIRIquant expects keys: name, tools.{bwa,hisat2,stringtie,samtools}, database.{gtf,genome,...}
    """
    genome = cfg["genome"]
    cirique = {
        "name": genome["species"],
        "tools": {
            "bwa":       "bwa",
            "hisat2":    "hisat2",
            "stringtie": "stringtie",
            "samtools":  "samtools",
        },
        "database": {
            "gtf":           genome["gtf"],
            "genome":        genome["fasta"],
            "hisat2_genome": genome["hisat2_index"],
        },
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as fh:
        yaml.dump(cirique, fh, default_flow_style=False)
    print(f"[OK] CIRIquant config written → {output}")


# ── Tool availability check ──────────────────────────────────────────────────

_REQUIRED_TOOLS = [
    "fastqc", "fastp", "multiqc",
    "prefetch", "fasterq-dump", "pigz",
    "CIRIquant",
    "bwa", "hisat2", "samtools", "stringtie",
]

def check_tools(tools: list[str] | None = None) -> bool:
    """Return True if all tools are found on PATH; print missing ones."""
    tools = tools or _REQUIRED_TOOLS
    missing = []
    for tool in tools:
        result = subprocess.run(
            ["which", tool], capture_output=True, text=True
        )
        if result.returncode != 0:
            missing.append(tool)

    if missing:
        print(f"[WARN] Missing tools: {', '.join(missing)}")
        print("       Install them via conda (see envs/circrna.yaml) or module load.")
        return False
    print("[OK] All required tools found on PATH")
    return True


# ── Disk resource gate ───────────────────────────────────────────────────────

def disk_gb_resource(path: str = ".", safety_factor: float = 0.7) -> int:
    """Return a conservative --resources disk_gb ceiling: whatever is
    ACTUALLY free under `path` right now, times safety_factor, floored to an
    integer GB. This replaces a hardcoded GB number with one that reflects
    live disk state at launch time -- a fixed ceiling goes stale the moment
    real usage changes, in either direction (too low wastes capacity once
    something is cleaned up elsewhere; too high lets Snakemake schedule more
    concurrent disk-heavy jobs than the disk can actually hold).

    safety_factor (default 0.7) leaves headroom for: disk use by anything
    else on the server outside this pipeline's control, and error in any
    single rule's own disk_gb resource claim (itself currently an estimate,
    not a measurement -- see workflow/rules/circrna.smk's
    _DISK_GB_PER_SAMPLE_ESTIMATE comment). Not yet tuned against a measured
    per-sample peak; revisit once logs/disk_usage_{srr}.tsv (written by
    scripts/sample_disk_usage.py during a real run) gives one."""
    free_bytes = shutil.disk_usage(path).free
    return int(free_bytes / (1024 ** 3) * safety_factor)


# ── Snakemake runner ─────────────────────────────────────────────────────────

def run_snakemake(
    snakefile: str = "workflow/Snakefile",
    cores: int = 8,
    extra_args: list[str] | None = None,
    dry_run: bool = False,
    disk_check_path: str = ".",
) -> int:
    """Run Snakemake as a subprocess; return exit code."""
    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", snakefile,
        "--cores", str(cores),
        "--directory", ".",
        "--resources", f"disk_gb={disk_gb_resource(disk_check_path)}",
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if dry_run:
        cmd.extend(["--dry-run", "--quiet", "rules"])
    if extra_args:
        cmd.extend(extra_args)

    print(f"[CMD] {' '.join(cmd)}")
    result = subprocess.run(cmd)
    return result.returncode


# ── Report generation-era detection ──────────────────────────────────────────

# The columns rank_biomarkers.py's v2 output format added; generate_report.py
# (as of commit dd1da2b, 2026-10-04) requires every one of these to be present
# in biomarker_candidates_v2_<method>.tsv before it will build a report. A
# dataset's own report.html predates that requirement whenever its de/
# directory never got these files -- there is no other way to tell from the
# platform side, since report.html itself carries no version marker.
_V2_BIOMARKER_REQUIRED_COLS = {
    "rank", "rank_default", "pass_default",
    "three_method_sig", "effect_stable", "biomarker_score",
}
_V2_BIOMARKER_METHODS = ("edgeR", "DESeq2", "limma")


def _v2_biomarker_header_ok(tsv_path: Path) -> bool:
    """True if tsv_path exists and its header line has every v2 column.
    Reads only the first line -- these files can run into the thousands of
    rows, and the header is all this check needs."""
    if not tsv_path.exists():
        return False
    try:
        with open(tsv_path, encoding="utf-8") as fh:
            header = set(fh.readline().rstrip("\n").split("\t"))
    except OSError:
        return False
    return _V2_BIOMARKER_REQUIRED_COLS.issubset(header)


def is_legacy_biomarker_generation(results_dir: str) -> bool:
    """True if this dataset's de/ directory predates the rank_biomarkers.py
    v2 / generate_report.py dd1da2b migration (2026-10-04) -- i.e. the
    current generate_report.py would raise SystemExit if asked to rebuild
    this dataset's report.html right now, because at least one of the three
    per-method biomarker_candidates_v2_<method>.tsv files (or its
    .status.json sidecar) is missing, or lacks the v2 columns.

    Deliberately file-feature-based, not a hardcoded accession list: the
    2026-10-04 migration itself (~/run_rank_v2.sh, server-only, not in this
    repo) WAS a hardcoded list of exactly 15 accessions, transcribed by hand
    from CLAUDE.md's already-documented 15-dataset cross-reproducibility
    inclusion criteria (資料集特性稽核校正) rather than read from any single
    canonical list file -- see docs/platform_security_notes.md. A function
    keyed on accession names would need to be remembered and updated by hand
    on every future migration; checking the files this function actually
    needs does not."""
    R = Path(results_dir)
    for m in _V2_BIOMARKER_METHODS:
        tsv = R / f"de/biomarker_candidates_v2_{m}.tsv"
        status_json = R / f"de/biomarker_candidates_v2_{m}.status.json"
        if not _v2_biomarker_header_ok(tsv) or not status_json.exists():
            return True
    return False


# ── Deterministic top-N selection ────────────────────────────────────────────

def select_top_n(df, score_col: str, n: int, id_col: str = "circ_id"):
    """Return the first n rows of df sorted ascending by score_col, with
    id_col as a deterministic tie-break secondary key -- so two rows that
    share an exact score_col value (common with p-values/padj on sparse
    datasets, and not rare: GSE148036's DESeq2 top-50 cutoff ties 131 rows
    at once) always resolve to the same n rows on every run, regardless of
    the order df's rows arrived in or which sort algorithm pandas's default
    (unstable) quicksort happens to pick.

    Before this helper existed, the same `.sort_values(score_col).head(n)`
    pattern, with no secondary key, was independently written at four
    separate call sites in generate_report.py (the main report's heatmap
    pool selection) plus predict_interactions.py's miRNA/RBP query-set
    selection -- the latter feeds n_mirna/n_rbp, two of rank_biomarkers.py's
    six biomarker_score dimensions, so a tie there changes which dimensions
    a candidate's score is even computed from (see
    rank_biomarkers.py's interaction_queried / Layer-2 renormalization),
    not just which row displays where. One shared function, not five
    independent copies of the same unguarded pattern.

    id_col must be unique across df's rows (circ_id and gene_name both are,
    in this codebase's usage) for the tie-break itself to be unambiguous;
    this is the caller's responsibility, not checked here.

    kind="stable" is redundant once id_col fully disambiguates the sort key
    (no tie survives a (score_col, id_col) tuple that is itself unique per
    row) but is passed anyway as defense in depth, in case a future caller
    passes a non-unique id_col by mistake -- stable sort at least falls
    back to preserving df's own row order instead of an arbitrary one."""
    return df.sort_values([score_col, id_col], kind="stable").head(n)


# ── Directory scaffolding ────────────────────────────────────────────────────

def ensure_dirs() -> None:
    for d in [
        "metadata", "config", "logs",
        "results/qc/raw", "results/qc/fastp",
        "results/circRNA",
        "results/de", "results/plots",
    ]:
        Path(d).mkdir(parents=True, exist_ok=True)
