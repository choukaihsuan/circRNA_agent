#!/usr/bin/env python3
"""Lightweight import-chain preflight check.

Catches a silently-missing hard dependency between scripts/*.py modules
BEFORE any pipeline stage or web_ui.py route that would call it is reached.

Named case this check exists because of: scripts/consensus_filter.py does a
module-level `import run_params` (line 26, not inside any function, so it
fails at import time, not at call time) -- and so do scripts/rank_biomarkers.py
and scripts/predict_interactions.py, independently. In 2026-10, run_params.py
was found missing from the server's deployed copy of scripts/ while every one
of these three callers was present and otherwise intact. Nothing failed until
someone actually ran consensus_filter / rank_biomarkers / predict_interactions
mid-pipeline, at which point the traceback names run_params but gives no hint
that the real cause is a repo/server sync gap rather than a code bug. This
script's only job is to turn that class of failure into a startup-time
failure instead of a mid-pipeline one.

Scope: every scripts/*.py module that web_ui.py calls (as a subprocess, or
via a lazy `from X import Y` inside a function) or that a Snakemake
`script:`/shell: rule invokes (workflow/rules/*.smk) -- EXCLUDING web_ui.py
itself. If web_ui.py's own top-level imports were broken, the Python
interpreter would already refuse to start web_ui.py before this check could
ever run, so re-checking it here would be redundant.

For R scripts (analysis.R / isoform_switching.R / loo_stability.R), "import"
has no direct equivalent; this checks for `Rscript` on PATH and that each
script's own top-level `library(...)` calls name R packages that are
actually installed (`requireNamespace(pkg, quietly=TRUE)`) -- not that the
script itself is syntactically valid or produces correct output.

Severity is NOT uniform across the two halves, and that is deliberate:

- A missing PYTHON module is a hard failure (exit code 1, web_ui.py refuses
  to start). Those modules are web_ui.py's own runtime dependencies -- if
  one is missing, the web service cannot do its job, so refusing to come
  up at all is strictly better than coming up broken.
- A missing R PACKAGE is a WARNING ONLY -- printed, never affects the exit
  code. The R scripts run inside Snakemake rules, not inside the web
  process; a missing R package fails that one Snakemake rule the next time
  it runs, it does not fail the web interface. The R environment may also
  legitimately differ from the process running web_ui.py (the project's
  ciriquant conda env is a known, deliberately separate environment -- see
  CLAUDE.md's existing env-isolation rules), so "the web process's R
  toolchain doesn't see env X's packages" is an expected, not an alarming,
  state. Gating the whole platform on it would take down a working service
  over a dependency that service doesn't use.

This script only imports modules / checks package presence. It never calls
any module's main()/CLI entry point, never touches any results_dir,
config.yaml, or network resource, and never executes Snakemake.
"""
from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent

# Every scripts/*.py module invoked by web_ui.py (subprocess or lazy
# import) or by a Snakemake script:/shell: rule. Import-only below is
# safe because each of these already guards its main logic behind
# `if __name__ == "__main__":` (generate_report.py additionally behind
# `if "snakemake" in dir():`) -- importing the module cannot trigger any
# pipeline work.
PYTHON_MODULES = [
    "utils",
    "run_params",
    "consensus_filter",
    "merge_counts",
    "annotate_circbase",
    "assign_isoforms",
    "rank_biomarkers",
    "predict_interactions",
    "generate_report",
    "notify",
    "prepare_metadata",
    "download_geo",
    "agent",
    "sample_disk_usage",
]

# R scripts in scope, and the R packages their own library(...) calls name.
R_SCRIPTS = {
    "analysis.R": ["ggplot2", "pheatmap", "RColorBrewer", "dplyr",
                   "qvalue", "edgeR", "DESeq2", "limma"],
    "isoform_switching.R": ["dplyr", "tidyr", "tibble"],
    "loo_stability.R": ["edgeR", "DESeq2", "limma"],
}


def check_python_modules():
    sys.path.insert(0, str(SCRIPTS_DIR))
    failures = []
    for name in PYTHON_MODULES:
        path = SCRIPTS_DIR / f"{name}.py"
        if not path.is_file():
            failures.append((f"{name}.py", "file not found in scripts/"))
            continue
        # Drop any cached entry first so a stale sys.modules import from
        # earlier in the same process can't mask a real failure here.
        sys.modules.pop(name, None)
        try:
            importlib.import_module(name)
        except Exception as exc:  # ImportError/ModuleNotFoundError, and anything
            # else a module-level statement could raise during import.
            failures.append((f"{name}.py", f"{type(exc).__name__}: {exc}"))
    return failures


def check_r_scripts():
    failures = []
    rscript = shutil.which("Rscript")
    if not rscript:
        for script_name in R_SCRIPTS:
            failures.append((script_name, "Rscript not found on PATH -- cannot check R package dependencies"))
        return failures

    for script_name, pkgs in R_SCRIPTS.items():
        script_path = SCRIPTS_DIR / script_name
        if not script_path.is_file():
            failures.append((script_name, "file not found in scripts/"))
            continue
        pkg_vec = ", ".join(f'"{p}"' for p in pkgs)
        r_expr = (
            f"pkgs <- c({pkg_vec}); "
            "missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; "
            "cat(paste(missing, collapse=','))"
        )
        result = subprocess.run([rscript, "--vanilla", "-e", r_expr],
                                 capture_output=True, text=True, timeout=120)
        missing = [p for p in result.stdout.strip().split(",") if p]
        if missing:
            failures.append((script_name, f"missing R package(s): {', '.join(missing)}"))
        if result.returncode != 0 and not missing:
            failures.append((script_name, f"Rscript check itself failed (rc={result.returncode}): {result.stderr.strip()[:300]}"))
    return failures


def main() -> int:
    py_failures = check_python_modules()
    r_failures = check_r_scripts()

    if r_failures:
        print("preflight_imports: WARNING -- the R script(s) below have a "
              "missing dependency. This does NOT affect the exit code and "
              "does NOT stop web_ui.py from starting: these scripts run "
              "inside Snakemake rules, not inside the web process, and the "
              "R toolchain visible here may legitimately differ from the "
              "one used to run the actual pipeline (e.g. the ciriquant "
              "conda env). Install the listed package(s) in whichever R "
              "environment actually runs these Snakemake rules before "
              "relying on them.", file=sys.stderr)
        for name, detail in r_failures:
            print(f"  - {name}: {detail}", file=sys.stderr)

    if not py_failures:
        print("preflight_imports: Python OK -- all "
              f"{len(PYTHON_MODULES)} modules in scope imported cleanly.")
        return 0

    print("preflight_imports: FAILED -- the Python module(s) below have a "
          "missing dependency. This is the same error web_ui.py's own code "
          "would hit the first time it reaches that module; this check "
          "only moves that failure to startup time instead of into a "
          "running, partially-broken service.", file=sys.stderr)
    for name, detail in py_failures:
        print(f"  - {name}: {detail}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
