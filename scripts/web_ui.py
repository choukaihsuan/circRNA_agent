"""
Web UI for circRNA pipeline — tool and DE method selection.

Usage:
    cd ~/circRNA_agent
    conda activate ciriquant
    python scripts/web_ui.py          # opens on http://localhost:5000
    python scripts/web_ui.py --port 8080 --host 0.0.0.0
"""

from __future__ import annotations

import argparse
import json
import random
import re
import string
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml
import csv
import io
import sys as _sys
_sys.path.insert(0, str(Path(__file__).parent))

from flask import Flask, jsonify, redirect, render_template, request, url_for

app = Flask(__name__, template_folder="templates")
BASE_DIR = Path(__file__).parent.parent
CONFIG_PATH = BASE_DIR / "config.yaml"
LOG_PATH = BASE_DIR / "logs" / "pipeline_run.log"
REGISTRY_PATH = BASE_DIR / "jobs" / "registry.json"


# ── Job registry ──────────────────────────────────────────────────────────────

def generate_job_id(gse_id: str) -> str:
    suffix = "".join(random.choices(string.ascii_uppercase + string.digits, k=4))
    return f"{gse_id}-{suffix}"


def load_registry() -> dict:
    if not REGISTRY_PATH.exists():
        return {}
    try:
        return json.loads(REGISTRY_PATH.read_text())
    except (json.JSONDecodeError, OSError):
        return {}


def save_registry(registry: dict) -> None:
    REGISTRY_PATH.parent.mkdir(parents=True, exist_ok=True)
    REGISTRY_PATH.write_text(json.dumps(registry, ensure_ascii=False, indent=2))


def register_job(job_id: str, gse_id: str, log_path: Path, cores: int) -> None:
    registry = load_registry()
    registry[job_id] = {
        "job_id":     job_id,
        "gse_id":     gse_id,
        "start_time": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "log":        str(log_path.relative_to(BASE_DIR)),
        "cores":      cores,
    }
    save_registry(registry)


def job_log_path(job_id: str) -> Path:
    return BASE_DIR / "logs" / f"pipeline_{job_id}.log"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def save_config(cfg: dict) -> None:
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True, sort_keys=False)


def _project_config_path(gse_id: str) -> Path:
    return BASE_DIR / "config" / "projects" / f"{gse_id}.yaml"


def _project_meta_dir(gse_id: str) -> Path:
    return BASE_DIR / "metadata" / gse_id


def _configfile_for(gse_id: str) -> str:
    """Return configfile path (relative to BASE_DIR) for snakemake --configfile."""
    p = _project_config_path(gse_id)
    if p.exists():
        return str(p.relative_to(BASE_DIR))
    return "config.yaml"


def save_project_snapshot(cfg: dict) -> None:
    """Write config.yaml AND a per-project snapshot in config/projects/{id}.yaml."""
    save_config(cfg)
    gse_id = cfg.get("project_id", "")
    if gse_id:
        proj_path = _project_config_path(gse_id)
        proj_path.parent.mkdir(parents=True, exist_ok=True)
        with open(proj_path, "w") as f:
            yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True, sort_keys=False)


def load_project_config(gse_id: str) -> dict:
    """Load project-specific config if snapshot exists, else fall back to main config."""
    p = _project_config_path(gse_id)
    if p.exists():
        with open(p) as f:
            return yaml.safe_load(f)
    return load_config()


def _snake_bin() -> str:
    """Return absolute path to snakemake in the current conda env."""
    import sys as _sys
    return str(Path(_sys.executable).parent / "snakemake")


def _snake_env() -> dict:
    """Return os.environ copy with working SRA-tool paths prepended to PATH.

    Priority order (front → back):
      1. sra_env  — has fasterq-dump that works on CentOS 7 glibc 2.17
      2. circrna  — has aria2c for fast multi-connection S3 downloads
      3. base conda bin — has prefetch / srapath
    """
    import os, sys as _sys
    env = dict(os.environ)
    # envs/ciriquant/bin/python → miniconda3
    conda_root = Path(_sys.executable).parents[3]
    prepend = [
        str(conda_root / "envs" / "sra_env"  / "bin"),
        str(conda_root / "envs" / "circrna"  / "bin"),
        str(conda_root / "bin"),
    ]
    # Also check miniforge3 if present
    miniforge = Path(_sys.executable).parents[3].parent / "miniforge3"
    if miniforge.is_dir():
        prepend += [
            str(miniforge / "envs" / "circrna" / "bin"),
            str(miniforge / "bin"),
        ]
    path = env.get("PATH", "")
    extra = ":".join(p for p in prepend if p not in path)
    if extra:
        env["PATH"] = extra + ":" + path
    return env


def _update_paths_for_project(cfg: dict, new_pid: str) -> dict:
    """Replace the old project_id embedded in raw/trimmed/results paths,
    and ensure the download section exists with server-appropriate defaults."""
    old_pid = cfg.get("project_id", "")
    if old_pid and old_pid != new_pid:
        for key in ("raw_dir", "trimmed_dir", "results_dir"):
            if key in cfg and old_pid in cfg[key]:
                cfg[key] = cfg[key].replace(old_pid, new_pid)
        # Also update sra_cache_dir / tmp_dir inside the download section
        dl = cfg.get("download", {})
        for key in ("sra_cache_dir", "tmp_dir"):
            if key in dl and old_pid in dl[key]:
                dl[key] = dl[key].replace(old_pid, new_pid)
        if dl:
            cfg["download"] = dl

    # Ensure download section exists so download.smk doesn't KeyError
    if "download" not in cfg:
        # Derive a sensible cache dir from raw_dir (sibling directory)
        raw = cfg.get("raw_dir", "")
        parent = str(Path(raw).parent) if raw else str(Path.home())
        cfg["download"] = {
            "sra_cache_dir": parent + "/sra_cache",
            "tmp_dir":       parent + "/sra_tmp",
            "retry":         3,
            "ascp_speed":    "500m",
        }
    return cfg


PIPELINE_STAGES = [
    ("download_fastq",         "下載 SRA → FASTQ"),
    ("fastqc_raw",             "FastQC (raw)"),
    ("fastp_trim",             "fastp 品質過濾"),
    ("multiqc",                "MultiQC 報告"),
    ("check_ciriquant_config", "驗證 CIRIquant 設定"),
    ("ciriquant",              "CIRIquant"),
    ("star_align",             "STAR paired"),
    ("star_align_mate1",       "STAR mate1"),
    ("star_align_mate2",       "STAR mate2"),
    ("dcc",                    "DCC"),
    ("consensus_filter",       "Consensus filter"),
    ("merge_counts",           "合併計數矩陣"),
    ("annotate_circbase",      "circBase 注釋"),
    ("de_analysis",            "DE 分析"),
    ("rank_biomarkers",        "Biomarker 排序"),
    ("assign_isoforms",        "Isoform 分組"),
    ("isoform_switching",      "Isoform switching"),
    ("generate_report",        "HTML 報告"),
]


def _parse_one_run(run_text: str) -> tuple:
    """Parse one Snakemake run section.
    Returns (job_to_rule, rule_started, rule_done, rule_failed, finished_count, total_count).
    """
    job_to_rule: dict[str, str] = {}
    rule_started: dict[str, int] = {}
    rule_done:    dict[str, int] = {}
    rule_failed:  set[str]       = set()
    current_rule: Optional[str] = None
    finished_count = total_count = 0

    for line in run_text.splitlines():
        m = re.match(r"\s*(?:local)?rule\s+(\w+)\s*:", line)
        if m:
            current_rule = m.group(1)
            continue

        if current_rule:
            m = re.match(r"\s+jobid:\s*(\d+)", line)
            if m:
                jid = m.group(1)
                job_to_rule[jid] = current_rule
                rule_started[current_rule] = rule_started.get(current_rule, 0) + 1
                current_rule = None
                continue
            if line and not line[0].isspace():
                current_rule = None

        m = re.search(r"Finished job (\d+)\.", line)
        if m:
            rule = job_to_rule.get(m.group(1))
            if rule:
                rule_done[rule] = rule_done.get(rule, 0) + 1

        m = re.match(r"\s*Error in rule\s+(\w+)\s*:", line)
        if m:
            rule_failed.add(m.group(1))

        m = re.search(r"(\d+) of (\d+) steps", line)
        if m:
            finished_count = int(m.group(1))
            total_count    = int(m.group(2))

    return job_to_rule, rule_started, rule_done, rule_failed, finished_count, total_count


def parse_log_progress(log_text: str) -> dict:
    """Parse a Snakemake log and return per-rule status + overall progress.

    Strategy: parse ALL runs to build the set of rules that ever completed
    (handles steps done in early runs), then use the LAST run to determine
    current active status (error/running).  This prevents re-run retries from
    inflating started counts while still showing historical completions.
    """
    dag_marker = "Building DAG of jobs"
    parts = log_text.split(dag_marker)
    run_texts = [dag_marker + p for p in parts[1:]] if len(parts) > 1 else [log_text]

    # Rules that ever had at least one Finished job across all runs
    rule_ever_done: set[str] = set()
    rule_ever_done_count: dict[str, int] = {}
    for rt in run_texts:
        _, _, rd, _, _, _ = _parse_one_run(rt)
        for rule, cnt in rd.items():
            rule_ever_done.add(rule)
            rule_ever_done_count[rule] = rule_ever_done_count.get(rule, 0) + cnt

    # Find the last "main" Snakemake run start.
    # Subprocess jobs also print "Building DAG..." but have no Snakemake timestamp
    # entries ([Day Mon DD HH:MM:SS YYYY]) in their segment. We want the last
    # Building DAG segment that IS followed by a timestamp — that is the real
    # latest main-process run. We then use everything from that point onward as
    # the authoritative "current run" text so that rules logged before subprocess
    # Building DAG noise are still visible.
    ts_re = re.compile(r'^\[(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) ', re.MULTILINE)
    dag_positions = [m.start() for m in re.finditer(re.escape(dag_marker), log_text)]
    last_main_dag_pos = dag_positions[0] if dag_positions else 0
    for i, pos in enumerate(dag_positions):
        end = dag_positions[i + 1] if i + 1 < len(dag_positions) else len(log_text)
        if ts_re.search(log_text[pos:end]):
            last_main_dag_pos = pos
    last_run_text = log_text[last_main_dag_pos:]

    # Last run: authoritative for current started/done/failed/progress counts
    _, rule_started, rule_done, rule_failed, finished_count, total_count = \
        _parse_one_run(last_run_text)

    # If generate_report ever completed, all prerequisite steps must be done
    # (they were skipped by Snakemake because outputs already existed)
    pipeline_finished = "generate_report" in rule_ever_done

    stages = []
    for rule_id, label in PIPELINE_STAGES:
        started = rule_started.get(rule_id, 0)
        done    = rule_done.get(rule_id, 0)
        failed  = rule_id in rule_failed

        if failed and done < started:
            # Last run had an error and not all jobs finished → genuine failure
            status = "error"
        elif started > 0 and done >= started:
            # Last run completed all its jobs for this rule
            status = "done"
        elif started > 0:
            # Last run started but not done yet
            status = "running"
        elif rule_id in rule_ever_done:
            # Completed in an earlier run; not touched in last run
            status = "done"
            done    = rule_ever_done_count.get(rule_id, 1)
            started = done
        elif pipeline_finished:
            # Never appeared in any log (Snakemake skipped because outputs exist)
            # but the final step completed → infer this step was done
            status = "done"
            done    = 1
            started = 1
        else:
            status = "pending"

        stages.append({
            "id":      rule_id,
            "label":   label,
            "status":  status,
            "started": started,
            "done":    done,
        })

    # Second pass: if a later stage is running/done, all earlier pending stages
    # must have completed (Snakemake skipped them because outputs already existed).
    last_active_idx = -1
    for i, s in enumerate(stages):
        if s["status"] in ("running", "done"):
            last_active_idx = i
    for i in range(last_active_idx):
        if stages[i]["status"] == "pending":
            stages[i]["status"]  = "done"
            stages[i]["started"] = 1
            stages[i]["done"]    = 1

    return {
        "stages":         stages,
        "finished_count": finished_count,
        "total_count":    total_count,
    }


def _infer_stages_from_files(results_dir_str: str, stages: list,
                              cfg: Optional[dict] = None) -> list:
    """Supplement log-parsed stage statuses by checking actual output files.

    Called after parse_log_progress so file evidence can upgrade 'pending'
    stages to 'done' when web_ui was restarted after a completed pipeline.
    Also performs per-sample file counting to show N/M progress.
    """
    results_dir = Path(results_dir_str) if results_dir_str else None
    if not results_dir or not results_dir.exists():
        return stages

    stage_map = {s["id"]: s for s in stages}

    def mark_done(*ids):
        for sid in ids:
            if sid in stage_map and stage_map[sid]["status"] == "pending":
                stage_map[sid].update({"status": "done", "started": 1, "done": 1})

    def mark_partial(sid: str, n_done: int, n_total: int):
        """Update a stage's N/M counts from filesystem; upgrade pending→running/done."""
        if sid not in stage_map or n_done == 0:
            return
        s = stage_map[sid]
        s["started"] = n_total
        s["done"]    = n_done
        if n_done >= n_total:
            s["status"] = "done"
        elif s["status"] == "pending":
            s["status"] = "running"

    # ── Per-sample counting (requires sample list from config metadata) ──────
    sample_ids: list = []
    if cfg:
        meta_rel = cfg.get("metadata", "")
        if meta_rel:
            meta_path = Path(meta_rel)
            if not meta_path.is_absolute():
                meta_path = BASE_DIR / meta_path
            if meta_path.exists():
                try:
                    import csv as _csv
                    with open(meta_path, newline="") as f:
                        sample_ids = [r["srr_id"] for r in _csv.DictReader(f)
                                      if "srr_id" in r]
                except Exception:
                    pass

    n = len(sample_ids)
    if n > 0:
        raw_dir = Path(cfg.get("raw_dir", "")) if cfg else Path("")

        # download_fastq
        mark_partial("download_fastq",
                     sum(1 for s in sample_ids
                         if (raw_dir / f"{s}_1.fastq.gz").exists()), n)

        # FastQC (raw) — check for any fastqc zip in per-sample subdir or flat dir
        fqc_dir = results_dir / "qc" / "fastqc"
        mark_partial("fastqc_raw",
                     sum(1 for s in sample_ids
                         if list(fqc_dir.glob(f"{s}*_fastqc.html")) if fqc_dir.exists()), n)

        # fastp
        mark_partial("fastp_trim",
                     sum(1 for s in sample_ids
                         if (results_dir / "qc" / "fastp" / f"{s}.json").exists()), n)

        # CIRIquant
        mark_partial("ciriquant",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / f"{s}.gtf").exists()), n)

        # STAR paired
        mark_partial("star_align",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / "Chimeric.out.junction").exists()), n)

        # STAR mate1
        mark_partial("star_align_mate1",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / "mate1" / "Chimeric.out.junction").exists()), n)

        # STAR mate2
        mark_partial("star_align_mate2",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / "mate2" / "Chimeric.out.junction").exists()), n)

        # DCC
        mark_partial("dcc",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / "DCC" / "CircCoordinates").exists()), n)

        # consensus_filter
        mark_partial("consensus_filter",
                     sum(1 for s in sample_ids
                         if (results_dir / "circRNA" / s / "high_confidence.bed").exists()), n)

    # ── Milestone-based cascade (single-output rules) ─────────────────────────
    # check_ciriquant_config: independent of results_dir; uses project-level touch file
    if (BASE_DIR / "config" / ".ciriquant_ready").exists():
        mark_done("check_ciriquant_config")

    if (results_dir / "qc" / "multiqc_report.html").exists():
        mark_done("download_fastq", "fastqc_raw", "fastp_trim", "multiqc")

    if (results_dir / "circRNA" / "count_matrix.tsv").exists():
        mark_done("download_fastq", "fastqc_raw", "fastp_trim", "multiqc",
                  "check_ciriquant_config", "ciriquant",
                  "star_align", "star_align_mate1", "star_align_mate2",
                  "dcc", "consensus_filter", "merge_counts")

    if (results_dir / "circRNA" / "circbase_annotated.tsv").exists():
        mark_done("annotate_circbase")

    if (results_dir / "circRNA" / "isoform_groups.tsv").exists():
        mark_done("assign_isoforms")

    if (results_dir / "de" / "de_results.tsv").exists():
        mark_done("de_analysis")

    if (results_dir / "de" / "biomarker_candidates.tsv").exists():
        mark_done("rank_biomarkers")

    if (results_dir / "de" / "isoform_switching.tsv").exists():
        mark_done("isoform_switching")

    if (results_dir / "report.html").exists():
        for s in stages:
            if s["status"] == "pending":
                s.update({"status": "done", "started": 1, "done": 1})

    return stages


def tail_log(n: int = 80, path: Optional[Path] = None) -> str:
    p = path or LOG_PATH
    if not p.exists():
        return "（尚無 log 檔）"
    lines = p.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


def pipeline_is_running() -> bool:
    """Check if snakemake is actually running the pipeline (not monitoring scripts)."""
    try:
        result = subprocess.run(
            ["pgrep", "-f", "snakemake.*--snakefile.*workflow/Snakefile"],
            capture_output=True, text=True,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def kill_pipeline() -> None:
    """Kill all snakemake-related processes robustly using pgrep + kill -9."""
    import signal
    try:
        result = subprocess.run(
            ["pgrep", "-f", "snakemake.*workflow/Snakefile"],
            capture_output=True, text=True,
        )
        pids = [int(p) for p in result.stdout.split() if p.strip().isdigit()]
        for pid in pids:
            try:
                import os
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    except Exception:
        pass


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    cfg = load_config()
    tools = cfg.get("consensus", {}).get("tools", ["ciriquant", "dcc"])
    de_method = cfg.get("de", {}).get("method", "deseq2")
    registry = load_registry()
    recent_jobs = list(reversed(list(registry.items())))[:5]
    return render_template(
        "index.html",
        config=cfg,
        selected_tools=tools,
        de_method=de_method,
        running=pipeline_is_running(),
        recent_jobs=recent_jobs,
    )


@app.route("/update", methods=["POST"])
def update():
    cfg = load_config()

    # circRNA tools
    tools = request.form.getlist("tools")
    if not tools:
        tools = ["ciriquant"]
    cfg.setdefault("consensus", {})
    cfg["consensus"]["tools"] = tools
    cfg["consensus"]["min_tools"] = max(1, int(request.form.get("min_tools", len(tools))))

    # DE method
    cfg.setdefault("de", {})
    cfg["de"]["method"] = request.form.get("de_method", "edgeR_ciriquant")

    # Numeric params
    cfg["consensus"]["min_bsj_reads"]    = int(request.form.get("min_bsj", 2))
    cfg["consensus"]["slop"]             = int(request.form.get("slop", 10))
    cfg["consensus"]["max_junction_ratio"] = float(request.form.get("max_junction_ratio", 1.0))
    cfg["de"]["fdr_cutoff"]              = float(request.form.get("fdr", 0.05))
    cfg["de"]["log2fc_cutoff"]           = float(request.form.get("log2fc", 1.0))
    cfg["threads"]                       = int(request.form.get("threads", 8))
    cfg["study_title"]                   = request.form.get("study_title", "").strip()

    save_project_snapshot(cfg)

    if request.form.get("action") == "run":
        kill_pipeline()
        cores = cfg["threads"] * 4
        gse_id = cfg.get("project_id", "pipeline")
        job_id = generate_job_id(gse_id)
        log_path = job_log_path(job_id)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as log_f:
            subprocess.Popen(
                [
                    _snake_bin(),
                    "--snakefile", "workflow/Snakefile",
                    "--configfile", _configfile_for(gse_id),
                    "--cores", str(cores),
                    "--resources", f"mem_gb={cfg.get('resources', {}).get('mem_gb', 300)}",
                    "--keep-going",
                    "--rerun-incomplete",
                ],
                cwd=str(BASE_DIR),
                stdout=log_f,
                stderr=subprocess.STDOUT,
                env=_snake_env(),
            )
        register_job(job_id, gse_id, log_path, cores)
        return redirect(url_for("status_job", job_id=job_id))

    return redirect(url_for("index"))


@app.route("/run_gse", methods=["POST"])
def run_gse():
    gse_id = request.form.get("gse_id", "").strip()
    cores  = int(request.form.get("cores", 8))
    if not gse_id:
        return redirect(url_for("index"))

    job_id = generate_job_id(gse_id)
    log_path = job_log_path(job_id)

    cfg = load_project_config(gse_id)
    cfg["project_id"] = gse_id
    cfg["threads"]    = cores
    cfg = _update_paths_for_project(cfg, gse_id)

    # Pre-detect labels if metadata is already available (re-run scenario)
    _meta_path = BASE_DIR / cfg.get("metadata", f"metadata/{gse_id}/library_info.csv")
    if not _meta_path.exists():
        _meta_path = BASE_DIR / "metadata" / gse_id / "library_info.csv"
    if _meta_path.exists():
        try:
            import pandas as _pd
            from prepare_metadata import detect_group_labels as _detect
            _meta_df = _pd.read_csv(str(_meta_path))
            _case, _ctrl = _detect(_meta_df)
            cfg.setdefault("de", {})["tumor_label"]  = _case
            cfg.setdefault("de", {})["normal_label"]  = _ctrl
        except Exception:
            pass

    save_project_snapshot(cfg)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_f:
        subprocess.Popen(
            ["python", "scripts/agent.py", "--gse", gse_id, "--cores", str(cores)],
            cwd=str(BASE_DIR),
            stdout=log_f,
            stderr=subprocess.STDOUT,
        )
    register_job(job_id, gse_id, log_path, cores)
    return redirect(url_for("status_job", job_id=job_id))


@app.route("/run_manual", methods=["POST"])
def run_manual():
    """Start pipeline from manual SRR list or uploaded CSV."""
    cores      = int(request.form.get("cores", 8))
    project_id = request.form.get("project_id", "CUSTOM").strip() or "CUSTOM"
    tumor_label  = request.form.get("tumor_label",  "tumor").strip()  or "tumor"
    normal_label = request.form.get("normal_label", "normal").strip() or "normal"

    # ── Build library_info rows ────────────────────────────────────────────────
    rows: list[dict] = []

    # Option A: uploaded CSV file
    csv_file = request.files.get("csv_file")
    if csv_file and csv_file.filename:
        content = csv_file.read().decode("utf-8-sig")
        reader  = csv.DictReader(io.StringIO(content))
        for r in reader:
            srr = (r.get("srr_id") or r.get("SRR_ID") or r.get("run") or "").strip()
            if srr:
                rows.append({
                    "srr_id":    srr,
                    "condition": (r.get("condition") or r.get("group") or "tumor").strip(),
                    "paired":    (r.get("paired") or "true").strip(),
                })
    else:
        # Option B: manual SRR entries from the form
        srr_ids    = request.form.getlist("srr_id[]")
        conditions = request.form.getlist("condition[]")
        for srr, cond in zip(srr_ids, conditions):
            srr  = srr.strip()
            cond = cond.strip()
            if srr:
                rows.append({"srr_id": srr, "condition": cond, "paired": "true"})

    if not rows:
        return redirect(url_for("index"))

    # ── Write metadata files (shared + per-project copy) ──────────────────────
    meta_dir     = BASE_DIR / "metadata"
    proj_meta    = _project_meta_dir(project_id)
    meta_dir.mkdir(exist_ok=True)
    proj_meta.mkdir(parents=True, exist_ok=True)

    lib_path      = meta_dir     / "library_info.csv"
    grp_path      = meta_dir     / "sample_groups.csv"
    proj_lib_path = proj_meta    / "library_info.csv"
    proj_grp_path = proj_meta    / "sample_groups.csv"

    for lp in (lib_path, proj_lib_path):
        with open(lp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["srr_id", "paired"])
            w.writeheader()
            for r in rows:
                w.writerow({"srr_id": r["srr_id"], "paired": r.get("paired", "true")})

    for gp in (grp_path, proj_grp_path):
        with open(gp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["srr_id", "condition"])
            w.writeheader()
            for r in rows:
                w.writerow({"srr_id": r["srr_id"], "condition": r["condition"]})

    # ── Update config ───────────────────────────────────────────────────────────
    cfg = load_project_config(project_id)
    cfg = _update_paths_for_project(cfg, project_id)
    cfg["project_id"]          = project_id
    cfg["threads"]             = cores
    cfg["metadata"]            = str(proj_lib_path.relative_to(BASE_DIR))
    cfg["groups"]              = str(proj_grp_path.relative_to(BASE_DIR))
    cfg["de"]["tumor_label"]   = tumor_label
    cfg["de"]["normal_label"]  = normal_label
    notify_email = request.form.get("notify_email", "").strip()
    if notify_email:
        cfg.setdefault("notify", {})["email_to"] = notify_email
    # Auto-fetch GEO title if not already set
    if not cfg.get("study_title"):
        cfg["study_title"] = _fetch_geo_title(project_id)
    save_project_snapshot(cfg)

    # ── Kill any leftover pipeline before launching ────────────────────────────
    kill_pipeline()

    # ── Launch Snakemake via same Python env as web_ui ──────────────────────────
    _snake = _snake_bin()

    job_id   = generate_job_id(project_id)
    log_path = job_log_path(job_id)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_f:
        subprocess.Popen(
            [_snake,
             "--snakefile", "workflow/Snakefile",
             "--configfile", _configfile_for(project_id),
             "--cores", str(cores),
             "--keep-going", "--rerun-incomplete"],
            cwd=str(BASE_DIR),
            stdout=log_f,
            stderr=subprocess.STDOUT,
            env=_snake_env(),
        )
    register_job(job_id, project_id, log_path, cores)
    return redirect(url_for("status_job", job_id=job_id))


@app.route("/api/scan_fastq")
def api_scan_fastq():
    """Scan a server-side directory for paired FASTQ files."""
    path = request.args.get("path", "").strip()
    if not path:
        return jsonify({"error": "未指定路徑"}), 400
    p = Path(path)
    if not p.exists():
        return jsonify({"error": f"路徑不存在：{path}"}), 404
    if not p.is_dir():
        return jsonify({"error": f"非目錄：{path}"}), 400

    PAIR_PATTERNS = [
        ("_1.fastq.gz",      "_2.fastq.gz"),
        ("_R1.fastq.gz",     "_R2.fastq.gz"),
        ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
        ("_1.fq.gz",         "_2.fq.gz"),
        ("_R1.fq.gz",        "_R2.fq.gz"),
        ("_1.fastq",         "_2.fastq"),
        ("_R1.fastq",        "_R2.fastq"),
    ]

    seen: set = set()
    samples: list = []
    for r1_suf, r2_suf in PAIR_PATTERNS:
        for r1 in sorted(p.glob(f"*{r1_suf}")):
            r2 = p / (r1.name[: -len(r1_suf)] + r2_suf)
            if not r2.exists():
                continue
            base = r1.name[: -len(r1_suf)]
            if base in seen:
                continue
            seen.add(base)
            samples.append({
                "name":      base,
                "r1":        str(r1),
                "r2":        str(r2),
                "r1_name":   r1.name,
                "r2_name":   r2.name,
                "condition": "tumor",
            })

    return jsonify({"samples": samples, "count": len(samples), "path": str(p)})


@app.route("/run_local", methods=["POST"])
def run_local():
    """Start pipeline with local FASTQ files (server-side paths via symlinks)."""
    project_id   = request.form.get("project_id",   "LOCAL").strip() or "LOCAL"
    tumor_label  = request.form.get("tumor_label",  "tumor").strip() or "tumor"
    normal_label = request.form.get("normal_label", "normal").strip() or "normal"
    cores        = int(request.form.get("cores", 8))
    samples_json = request.form.get("samples_json", "[]")

    try:
        samples = json.loads(samples_json)
    except json.JSONDecodeError:
        return "Invalid samples JSON", 400

    if not samples:
        return redirect(url_for("index"))

    cfg     = load_project_config(project_id)
    cfg     = _update_paths_for_project(cfg, project_id)
    raw_dir = Path(cfg["raw_dir"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    rows: list = []
    for s in samples:
        name = (s.get("name") or "").strip()
        r1   = (s.get("r1")   or "").strip()
        r2   = (s.get("r2")   or "").strip()
        cond = (s.get("condition") or "tumor").strip()
        if not name or not r1:
            continue
        for link, target in [(raw_dir / f"{name}_1.fastq.gz", r1),
                             (raw_dir / f"{name}_2.fastq.gz", r2)]:
            if link.exists() or link.is_symlink():
                link.unlink()
            if target:
                link.symlink_to(target)
        rows.append({"srr_id": name, "condition": cond})

    if not rows:
        return redirect(url_for("index"))

    meta_dir  = BASE_DIR / "metadata"
    proj_meta = _project_meta_dir(project_id)
    meta_dir.mkdir(exist_ok=True)
    proj_meta.mkdir(parents=True, exist_ok=True)

    for lp in (meta_dir / "library_info.csv", proj_meta / "library_info.csv"):
        with open(lp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["srr_id", "paired"])
            w.writeheader()
            for r in rows:
                w.writerow({"srr_id": r["srr_id"], "paired": "true"})

    for gp in (meta_dir / "sample_groups.csv", proj_meta / "sample_groups.csv"):
        with open(gp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["srr_id", "condition"])
            w.writeheader()
            for r in rows:
                w.writerow({"srr_id": r["srr_id"], "condition": r["condition"]})

    cfg["project_id"]                        = project_id
    cfg["threads"]                           = cores
    cfg["metadata"]                          = str((proj_meta / "library_info.csv").relative_to(BASE_DIR))
    cfg["groups"]                            = str((proj_meta / "sample_groups.csv").relative_to(BASE_DIR))
    cfg.setdefault("de", {})["tumor_label"]  = tumor_label
    cfg["de"]["normal_label"]                = normal_label
    notify_email = request.form.get("notify_email", "").strip()
    if notify_email:
        cfg.setdefault("notify", {})["email_to"] = notify_email
    save_project_snapshot(cfg)

    kill_pipeline()
    _snake   = _snake_bin()
    job_id   = generate_job_id(project_id)
    log_path = job_log_path(job_id)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w") as log_f:
        subprocess.Popen(
            [_snake,
             "--snakefile", "workflow/Snakefile",
             "--configfile", _configfile_for(project_id),
             "--cores", str(cores),
             "--keep-going", "--rerun-incomplete"],
            cwd=str(BASE_DIR),
            stdout=log_f,
            stderr=subprocess.STDOUT,
            env=_snake_env(),
        )
    register_job(job_id, project_id, log_path, cores)
    return redirect(url_for("status_job", job_id=job_id))


@app.route("/status")
def status():
    registry = load_registry()
    if registry:
        latest = list(registry.keys())[-1]
        return redirect(url_for("status_job", job_id=latest))
    return render_template(
        "status.html",
        log=tail_log(),
        running=pipeline_is_running(),
        job=None,
        job_id=None,
    )


@app.route("/status/<job_id>")
def status_job(job_id: str):
    registry = load_registry()
    job = registry.get(job_id)
    if job is None:
        return render_template(
            "status.html",
            log=f"找不到任務編號：{job_id}",
            running=False,
            job=None,
            job_id=job_id,
        ), 404
    log_path = BASE_DIR / job["log"]
    return render_template(
        "status.html",
        log=tail_log(path=log_path),
        running=pipeline_is_running(),
        job=job,
        job_id=job_id,
    )


@app.route("/status-query")
def status_query():
    job_id = request.args.get("job_id", "").strip()
    if not job_id:
        return redirect(url_for("index"))
    return redirect(url_for("status_job", job_id=job_id))


@app.route("/api/log")
def api_log():
    job_id = request.args.get("job")
    path: Optional[Path] = None
    if job_id:
        registry = load_registry()
        job = registry.get(job_id)
        if job:
            path = BASE_DIR / job["log"]
    return jsonify({"log": tail_log(50, path=path), "running": pipeline_is_running()})


@app.route("/api/progress")
def api_progress():
    job_id = request.args.get("job")
    path: Optional[Path] = None
    gse_id = None
    if job_id:
        registry = load_registry()
        job = registry.get(job_id)
        if job:
            path = BASE_DIR / job["log"]
            gse_id = job.get("gse_id")
    log_text = (path or LOG_PATH).read_text(errors="replace") \
               if (path or LOG_PATH).exists() else ""
    data = parse_log_progress(log_text)
    data["running"] = pipeline_is_running()
    # Supplement with filesystem evidence so status survives web_ui restarts
    report_exists = False
    qc_exists = False
    if gse_id:
        cfg = load_project_config(gse_id)
        results_dir = cfg.get("results_dir", "")
        if results_dir:
            data["stages"] = _infer_stages_from_files(results_dir, data["stages"], cfg)
            report_exists = (Path(results_dir) / "report.html").exists()
            qc_exists = (Path(results_dir) / "qc" / "multiqc_report.html").exists()
    data["report_exists"] = report_exists
    data["qc_exists"] = qc_exists
    return jsonify(data)


def _fetch_geo_title(gse_id: str) -> str:
    """Fetch GEO study title via NCBI eUtils. Returns empty string on failure."""
    import re as _re, urllib.request as _ur, json as _js
    m = _re.match(r"GSE(\d+)", gse_id.upper())
    if not m:
        return ""
    uid = "200" + m.group(1)
    url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
           f"?db=gds&id={uid}&retmode=json")
    try:
        with _ur.urlopen(url, timeout=8) as resp:
            data = _js.loads(resp.read())
        return data.get("result", {}).get(uid, {}).get("title", "")
    except Exception:
        return ""


@app.route("/api/detect_labels")
def api_detect_labels():
    """Detect case/control labels and GEO title for a given GSE."""
    gse_id = request.args.get("gse", "").strip()
    if not gse_id:
        return jsonify({"error": "missing gse"}), 400

    # Fetch GEO title from NCBI eUtils (non-blocking; empty string on failure)
    geo_title = _fetch_geo_title(gse_id)

    # Look for metadata in standard locations
    candidates = [
        BASE_DIR / "metadata" / gse_id / "library_info.csv",
        BASE_DIR / "metadata" / "library_info.csv",
    ]
    meta_path = next((p for p in candidates if p.exists()), None)
    if meta_path is None:
        return jsonify({"detected": False, "case": "tumor", "control": "normal",
                        "geo_title": geo_title, "note": "metadata not found yet"})

    try:
        import pandas as _pd
        from prepare_metadata import detect_group_labels as _detect
        df = _pd.read_csv(str(meta_path))
        case_label, ctrl_label = _detect(df)
        return jsonify({"detected": True, "case": case_label, "control": ctrl_label,
                        "rows": len(df), "geo_title": geo_title})
    except Exception as exc:
        return jsonify({"detected": False, "case": "tumor", "control": "normal",
                        "geo_title": geo_title, "error": str(exc)})


@app.route("/report/<job_id>")
def serve_report(job_id: str):
    """Serve the HTML report for a finished job."""
    from flask import send_file, abort
    registry = load_registry()
    job = registry.get(job_id)
    if not job:
        abort(404)
    gse_id = job["gse_id"]
    cfg = load_project_config(gse_id)
    results_dir = cfg.get("results_dir", "")
    report_path = Path(results_dir) / "report.html"
    if not report_path.exists():
        abort(404)
    return send_file(str(report_path), mimetype="text/html",
                     as_attachment=False, download_name=f"{gse_id}_report.html")


@app.route("/download/<job_id>")
def download_report(job_id: str):
    """Force-download the HTML report."""
    from flask import send_file, abort
    registry = load_registry()
    job = registry.get(job_id)
    if not job:
        abort(404)
    gse_id = job["gse_id"]
    cfg = load_project_config(gse_id)
    results_dir = cfg.get("results_dir", "")
    report_path = Path(results_dir) / "report.html"
    if not report_path.exists():
        abort(404)
    return send_file(str(report_path), mimetype="text/html",
                     as_attachment=True, download_name=f"{gse_id}_report.html")


@app.route("/qc/<job_id>")
def serve_qc(job_id: str):
    """Serve the MultiQC report for a finished job."""
    from flask import send_file
    registry = load_registry()
    job = registry.get(job_id)
    if not job:
        return "<h2>找不到此 job。</h2>", 404
    gse_id = job["gse_id"]
    cfg = load_project_config(gse_id)
    results_dir = cfg.get("results_dir", "")
    qc_path = Path(results_dir) / "qc" / "multiqc_report.html"
    if not qc_path.exists():
        return (
            f"<html><head><meta charset='utf-8'></head><body style='font-family:sans-serif;padding:40px'>"
            f"<h2>⏳ MultiQC 報告尚未產生</h2>"
            f"<p>MultiQC 需要所有 sample 的 FastQC 完成後才會執行。<br>"
            f"目前仍有 sample 正在 QC / 下載中，請等 pipeline 全部完成後再試。</p>"
            f"<p style='color:#888'>預期路徑：{qc_path}</p>"
            f"<a href='javascript:history.back()'>← 返回</a>"
            f"</body></html>"
        ), 202
    return send_file(str(qc_path), mimetype="text/html",
                     as_attachment=False, download_name=f"{gse_id}_multiqc_report.html")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=5000)
    args = parser.parse_args()
    print(f"  circRNA Pipeline UI  →  http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False)
