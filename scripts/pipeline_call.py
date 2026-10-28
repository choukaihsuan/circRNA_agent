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

SENTINEL = "===CIRCDEX_RESULT==="


def _emit(result) -> None:
    print(SENTINEL)
    print(json.dumps(result))


def cmd_preflight_accession(accession: str) -> None:
    from download_geo import preflight_accession
    _emit(preflight_accession(accession))


def cmd_detect_group_labels(csv_path: str) -> None:
    import pandas as pd
    from prepare_metadata import detect_group_labels
    df = pd.read_csv(csv_path)
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
