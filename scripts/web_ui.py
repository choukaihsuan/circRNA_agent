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
import re
import subprocess
from pathlib import Path

import yaml
from flask import Flask, jsonify, redirect, render_template, request, url_for

app = Flask(__name__, template_folder="templates")
BASE_DIR = Path(__file__).parent.parent
CONFIG_PATH = BASE_DIR / "config.yaml"
LOG_PATH = BASE_DIR / "logs" / "pipeline_run.log"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def save_config(cfg: dict) -> None:
    with open(CONFIG_PATH, "w") as f:
        yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True, sort_keys=False)


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


def parse_log_progress(log_text: str) -> dict:
    """Parse a Snakemake log and return per-rule status + overall progress."""
    job_to_rule: dict[str, str] = {}
    rule_started: dict[str, int] = {}
    rule_done:    dict[str, int] = {}
    rule_failed:  set[str]       = set()
    current_rule: str | None     = None
    finished_count = total_count = 0

    for line in log_text.splitlines():
        # "rule ciriquant:" or "localrule all:"
        m = re.match(r"\s*(?:local)?rule\s+(\w+)\s*:", line)
        if m:
            current_rule = m.group(1)
            continue

        # "    jobid: 5"  (indented, follows immediately after rule block)
        if current_rule:
            m = re.match(r"\s+jobid:\s*(\d+)", line)
            if m:
                jid = m.group(1)
                job_to_rule[jid]      = current_rule
                rule_started[current_rule] = rule_started.get(current_rule, 0) + 1
                current_rule = None
                continue
            # Non-indented non-empty line → lost rule context
            if line and not line[0].isspace():
                current_rule = None

        # "Finished job 5."
        m = re.search(r"Finished job (\d+)\.", line)
        if m:
            rule = job_to_rule.get(m.group(1))
            if rule:
                rule_done[rule] = rule_done.get(rule, 0) + 1

        # "Error in rule ciriquant:"
        m = re.match(r"\s*Error in rule\s+(\w+)\s*:", line)
        if m:
            rule_failed.add(m.group(1))

        # "3 of 19 steps (16%) done"
        m = re.search(r"(\d+) of (\d+) steps", line)
        if m:
            finished_count = int(m.group(1))
            total_count    = int(m.group(2))

    stages = []
    for rule_id, label in PIPELINE_STAGES:
        started = rule_started.get(rule_id, 0)
        done    = rule_done.get(rule_id, 0)
        failed  = rule_id in rule_failed
        if failed:
            status = "error"
        elif started > 0 and done >= started:
            status = "done"
        elif started > 0:
            status = "running"
        else:
            status = "pending"
        stages.append({
            "id":      rule_id,
            "label":   label,
            "status":  status,
            "started": started,
            "done":    done,
        })

    return {
        "stages":         stages,
        "finished_count": finished_count,
        "total_count":    total_count,
    }


def tail_log(n: int = 80) -> str:
    if not LOG_PATH.exists():
        return "（尚無 log 檔）"
    lines = LOG_PATH.read_text(errors="replace").splitlines()
    return "\n".join(lines[-n:])


def pipeline_is_running() -> bool:
    try:
        result = subprocess.run(
            ["pgrep", "-f", "snakemake.*circRNA_agent"],
            capture_output=True, text=True,
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    cfg = load_config()
    tools = cfg.get("consensus", {}).get("tools", ["ciriquant", "dcc"])
    de_method = cfg.get("de", {}).get("method", "deseq2")
    return render_template(
        "index.html",
        config=cfg,
        selected_tools=tools,
        de_method=de_method,
        running=pipeline_is_running(),
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

    save_config(cfg)

    if request.form.get("action") == "run":
        cores = cfg["threads"] * 4
        LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
        with open(LOG_PATH, "w") as log_f:
            subprocess.Popen(
                [
                    "snakemake",
                    "--snakefile", "workflow/Snakefile",
                    "--configfile", "config.yaml",
                    "--cores", str(cores),
                    "--resources", f"mem_gb={cfg.get('resources', {}).get('mem_gb', 300)}",
                    "--keep-going",
                    "--rerun-incomplete",
                ],
                cwd=str(BASE_DIR),
                stdout=log_f,
                stderr=subprocess.STDOUT,
            )
        return redirect(url_for("status"))

    return redirect(url_for("index"))


@app.route("/run_gse", methods=["POST"])
def run_gse():
    gse_id = request.form.get("gse_id", "").strip()
    cores  = int(request.form.get("cores", 8))
    if not gse_id:
        return redirect(url_for("index"))

    cfg = load_config()
    cfg["project_id"] = gse_id
    cfg["threads"]    = cores
    save_config(cfg)

    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "w") as log_f:
        subprocess.Popen(
            ["python", "scripts/agent.py", "--gse", gse_id, "--cores", str(cores)],
            cwd=str(BASE_DIR),
            stdout=log_f,
            stderr=subprocess.STDOUT,
        )
    return redirect(url_for("status"))


@app.route("/status")
def status():
    return render_template(
        "status.html",
        log=tail_log(),
        running=pipeline_is_running(),
    )


@app.route("/api/log")
def api_log():
    return jsonify({"log": tail_log(50), "running": pipeline_is_running()})


@app.route("/api/progress")
def api_progress():
    log_text = LOG_PATH.read_text(errors="replace") if LOG_PATH.exists() else ""
    data = parse_log_progress(log_text)
    data["running"] = pipeline_is_running()
    return jsonify(data)


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=5000)
    args = parser.parse_args()
    print(f"  circRNA Pipeline UI  →  http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False)
