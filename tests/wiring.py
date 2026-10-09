"""Helpers for the config -> command-line wiring tests.

We render the *real* shell commands Snakemake would run by doing a dry-run
(``snakemake -n -p --forceall``) inside a throw-away fixture project. Nothing is executed, no
published results are touched, no real data is needed (genome files are empty placeholders).

Why a dry-run and not a lighter parser: the commands are assembled from module-level Python in
the .smk files (``adaptive_flag``, ``USE_CIRIQUANT``-conditional inputs, lambdas, ``expand``,
``shell.prefix``) and Snakemake's own wildcard/params resolution. Re-implementing that in a
regex would test our reimplementation, not the pipeline - and silent wiring gaps are exactly
the bugs that live in that glue. A dry-run takes ~1 s and is read-only.
"""
from __future__ import annotations

import copy
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

import yaml

ROOT = Path(__file__).resolve().parent.parent
SNAKEFILE = ROOT / "workflow" / "Snakefile"
BASE_CONFIG = ROOT / "config.yaml"


def snakemake_bin() -> Optional[str]:
    cand = Path(sys.executable).parent / "snakemake"
    return str(cand) if cand.exists() else shutil.which("snakemake")


def build_fixture(fx: Path, cfg: dict) -> None:
    """Write config.yaml + placeholder inputs under ``fx`` so the DAG can be built."""
    cfg = copy.deepcopy(cfg)
    (fx / "metadata").mkdir(parents=True, exist_ok=True)
    (fx / "config").mkdir(exist_ok=True)
    (fx / "ref").mkdir(exist_ok=True)
    (fx / "metadata" / "library_info.csv").write_text(
        "srr_id,paired\nSRR1,true\nSRR2,true\nSRR3,true\nSRR4,true\n")
    (fx / "metadata" / "sample_groups.csv").write_text(
        "srr_id,condition\nSRR1,tumor\nSRR2,tumor\nSRR3,normal\nSRR4,normal\n")
    shutil.copy(ROOT / "config" / "ciriquant.yaml", fx / "config" / "ciriquant.yaml")
    shutil.copy(ROOT / "config" / "ciriquant_container.yaml", fx / "config" / "ciriquant_container.yaml")
    for k, v in list(cfg.get("genome", {}).items()):
        if k == "species":
            continue
        p = fx / "ref" / k
        p.write_text("")
        cfg["genome"][k] = str(p)
    cfg["metadata"] = "metadata/library_info.csv"
    cfg["groups"] = "metadata/sample_groups.csv"
    for k in ("raw_dir", "trimmed_dir", "results_dir"):
        cfg[k] = str(fx / "out" / k)
    cfg["use_container"] = False
    cfg["circbase_file"] = ""
    (fx / "config.yaml").write_text(yaml.safe_dump(cfg))


_HDR = re.compile(r"^(?:local)?rule (\w+):\s*$")


def parse_dry_run(text: str, fx: Path) -> Dict[str, List[str]]:
    """{rule_name: sorted unique normalised shell commands}. Rules without shell -> []."""
    out: Dict[str, set] = {}
    rule = None
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        m = _HDR.match(lines[i])
        if m:
            rule = m.group(1)
            out.setdefault(rule, set())
        elif lines[i].startswith("Shell command:") and rule:
            buf = []
            j = i + 1
            while j < len(lines) and lines[j].strip() != "" and not lines[j].startswith("["):
                buf.append(lines[j].strip())
                j += 1
            cmd = " ".join(" ".join(buf).split())
            if cmd and cmd != "None":
                cmd = cmd.replace(str(fx), "<FX>")
                cmd = re.sub(r"\bSRR\d\b", "SRRn", cmd)         # per-sample wildcards collapse
                out[rule].add(cmd)
            i = j
            continue
        i += 1
    return {k: sorted(v) for k, v in out.items()}


def render(fx: Path, cfg: dict) -> Dict[str, List[str]]:
    exe = snakemake_bin()
    if not exe:
        raise RuntimeError("snakemake not installed")
    build_fixture(fx, cfg)
    proc = subprocess.run(
        [exe, "-n", "-p", "--forceall", "--cores", "4", "--quiet", "rules", "--snakefile",
         str(SNAKEFILE), "--directory", str(fx)],
        capture_output=True, text=True, timeout=120)
    # --quiet rules suppresses the rule table only; re-run without it to keep the command text
    proc = subprocess.run(
        [exe, "-n", "-p", "--forceall", "--cores", "4", "--snakefile", str(SNAKEFILE),
         "--directory", str(fx)],
        capture_output=True, text=True, timeout=120)
    text = proc.stdout + "\n" + proc.stderr
    if proc.returncode != 0:
        raise RuntimeError("dry-run failed:\n" + text[-2000:])
    return parse_dry_run(text, fx)


def load_base_config() -> dict:
    return yaml.safe_load(BASE_CONFIG.read_text())


def leaf_paths(d: dict, prefix=()) -> List[tuple]:
    out = []
    for k, v in d.items():
        if isinstance(v, dict):
            out += leaf_paths(v, prefix + (k,))
        else:
            out.append(prefix + (k,))
    return out


def get_path(d: dict, path: tuple):
    for p in path:
        d = d[p]
    return d


def set_path(d: dict, path: tuple, value) -> None:
    for p in path[:-1]:
        d = d.setdefault(p, {})
    d[path[-1]] = value


def deviate(value):
    """A value clearly different from the default and unlikely to collide with other numbers."""
    if isinstance(value, bool):
        return not value
    if isinstance(value, int):
        return value + 7 if value != 7 else value + 5     # 10 -> 17 (distinct, easy to grep)
    if isinstance(value, float):
        return round(value * 1.7 + 0.0123, 4)
    if isinstance(value, str):
        return value + "_X"
    return value
