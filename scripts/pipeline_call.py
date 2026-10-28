#!/usr/bin/env python3
"""Thin CLI dispatcher so web_ui.py can call a handful of pipeline functions
(prepare_metadata.detect_group_labels, download_geo.preflight_accession) as a
subprocess instead of importing them in-process.

Why this exists: both of those modules `import pandas as pd` at module
level. Importing them directly inside a Flask route handler means pandas
runs inside the web process itself -- any ImportError (e.g. once web_ui.py
runs from an env that doesn't carry pandas) crashes that route, and any
crash/OOM inside pandas's CSV parsing takes the web worker down with it,
not just the one request. Going through a subprocess means the web process
never imports pandas at all, and a crash on the pipeline side is just a
failed subprocess call.

Contract with the caller: on success, prints exactly one line that is the
sentinel `===CIRCDEX_RESULT===`, followed by exactly one line of JSON, to
stdout, then exits 0. The caller must only look at the line immediately
after the sentinel line -- anything before it (library warnings, stray
prints) is not part of the contract and must be ignored. On failure, the
sentinel is never printed; a traceback goes to stderr and the process exits
non-zero. Callers that need to show the user why something failed should
show (a truncated) stderr, not just "failed".
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

SENTINEL = "===CIRCDEX_RESULT==="

# Same derivation as web_ui.py's BASE_DIR (Path(__file__).parent.parent),
# computed independently here rather than imported from web_ui.py -- this
# script must stay importable/runnable without Flask on the path (see the
# module docstring and scripts/preflight_imports.py's own note on keeping
# this standalone), so it can't just `from web_ui import BASE_DIR`.
ALLOWED_ROOT = (Path(__file__).resolve().parent.parent / "metadata").resolve()


def _emit(result) -> None:
    print(SENTINEL)
    print(json.dumps(result))


def cmd_preflight_accession(accession: str) -> None:
    from download_geo import preflight_accession
    _emit(preflight_accession(accession))


def cmd_detect_group_labels(csv_path: str) -> None:
    # Containment check BEFORE anything touches the filesystem with this
    # path. This belongs here, in the subprocess, not only in web_ui.py:
    # the parent-side check is a UX nicety (fail fast, nice error), this one
    # is the actual security boundary -- this script can in principle be
    # invoked by something other than web_ui.py later, and a boundary that
    # only exists in one specific caller isn't a boundary.
    #
    # Path.relative_to(), not str(p).startswith(str(root)) -- the prefix
    # form is the exact confinement bug fixed in 6ac3278
    # (_is_within_root()): it also accepts a sibling directory whose name
    # merely starts with the same characters as the allowed root (e.g.
    # root=.../metadata wrongly matching .../metadata_evil/...).
    #
    # resolve() first so a traversal sequence (../../etc/passwd) or a
    # symlink pointing outside ALLOWED_ROOT is caught by relative_to()
    # rather than by string shape. The rejection message says nothing about
    # whether the path exists, is a file vs directory, or why resolution
    # failed -- only that it's outside the allowed directory -- so this
    # can't be used to probe the filesystem via the stderr tail the caller
    # surfaces to the user.
    resolved = Path(csv_path).resolve()
    try:
        resolved.relative_to(ALLOWED_ROOT)
    except ValueError:
        raise PermissionError(
            "csv_path is outside the allowed metadata directory; rejected")

    import pandas as pd
    from prepare_metadata import detect_group_labels
    df = pd.read_csv(resolved)
    case_label, ctrl_label = detect_group_labels(df)
    _emit({"case": case_label, "control": ctrl_label, "rows": len(df)})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p1 = sub.add_parser("preflight-accession")
    p1.add_argument("accession")

    p2 = sub.add_parser("detect-group-labels")
    p2.add_argument("csv_path")

    args = parser.parse_args()

    try:
        if args.command == "preflight-accession":
            cmd_preflight_accession(args.accession)
        elif args.command == "detect-group-labels":
            cmd_detect_group_labels(args.csv_path)
    except Exception:
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
