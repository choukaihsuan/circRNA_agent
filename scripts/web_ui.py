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
    cfg["consensus"]["min_bsj_reads"] = int(request.form.get("min_bsj", 2))
    cfg["consensus"]["slop"]          = int(request.form.get("slop", 10))
    cfg["de"]["fdr_cutoff"]           = float(request.form.get("fdr", 0.05))
    cfg["de"]["log2fc_cutoff"]        = float(request.form.get("log2fc", 1.0))
    cfg["threads"]                    = int(request.form.get("threads", 8))

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


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=5000)
    args = parser.parse_args()
    print(f"  circRNA Pipeline UI  →  http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False)
