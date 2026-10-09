#!/usr/bin/env python3
"""Read-only audit: does each project snapshot's study_title match GEO?

Run on the server (it needs network access to NCBI and the real config/projects/*.yaml):

    python scripts/audit_study_titles.py            # table of OK / MISMATCH / EMPTY / N-A
    python scripts/audit_study_titles.py --json     # machine readable

It never writes anything. Background: before the fix in `web_ui.py`, a new dataset inherited the
previous project's `study_title`, and 5-digit GSE ids were looked up under the wrong GEO UID, so
some already-published reports may carry another dataset's title.
"""
import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from geo_title import fetch_geo_title, geo_uid  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


def _norm(t: str) -> str:
    return " ".join((t or "").lower().split())


def audit(projects_dir: Path, fetch=fetch_geo_title):
    rows = []
    for f in sorted(projects_dir.glob("*.yaml")):
        pid = f.stem
        cfg = yaml.safe_load(f.read_text()) or {}
        stored = str(cfg.get("study_title", "") or "")
        if geo_uid(pid) is None:                      # SRP*/PRJNA*/custom: no GEO series record
            status, geo = "N-A", ""
        else:
            geo = fetch(pid)
            if not geo:
                status = "GEO-UNREACHABLE"
            elif not stored:
                status = "EMPTY"
            elif _norm(stored) == _norm(geo):
                status = "OK"
            else:
                status = "MISMATCH"
        rows.append({"project": pid, "status": status, "stored": stored, "geo": geo})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--projects-dir", default=str(ROOT / "config" / "projects"))
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args(argv)
    rows = audit(Path(a.projects_dir))
    if a.json:
        print(json.dumps(rows, ensure_ascii=False, indent=2))
    else:
        for r in rows:
            print("%-12s %-16s stored=%r" % (r["project"], r["status"], r["stored"][:70]))
            if r["status"] == "MISMATCH":
                print("%-12s %-16s geo   =%r" % ("", "", r["geo"][:70]))
    return 1 if any(r["status"] == "MISMATCH" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
