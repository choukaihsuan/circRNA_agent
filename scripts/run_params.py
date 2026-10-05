"""Per-step execution records (run_params), written by the step itself at the end of its work.

Layout:
  {results_dir}/run_params/{rule}.yaml                 rules without a wildcard
  {results_dir}/run_params/{rule}.{wildcard}.yaml      one file per wildcard value (e.g. method, srr)
Each step writes only its own file, so parallel jobs never read-modify-write a shared file, and no lock
is needed (NFS locks are unreliable here). Files are written to a temporary name and renamed, so a
reader never sees a partial file.

Consolidation: `python scripts/run_params.py consolidate <results_dir> <project_id>` merges the per-step
files into {results_dir}/run_params.yaml. It is a separate, single-process step, run after the rules.
It is NOT a Snakemake output: declaring it as one would make every published dataset "missing output".

Provenance: records written here have provenance "executed". Values reconstructed afterwards from config,
code defaults or dry-run expansion must go to {results_dir}/run_params.reconstructed.yaml with provenance
"reconstructed", never into run_params.yaml.

canonical_key(record) gives the cache key for item L3 (see its docstring).
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

SCHEMA_VERSION = 1
HASH_LIMIT_BYTES = 1 << 30  # files above 1 GB are recorded by size and mtime only (no sha256)
REPO = Path(__file__).resolve().parent.parent


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _git(args: list[str]) -> str | None:
    try:
        # cwd, not "git -C": the CentOS 7 server has git too old for -C and would silently record null
        out = subprocess.run(["git", *args], cwd=str(REPO), capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired):
        return None
    return out.stdout.strip() if out.returncode == 0 else None


def code_info(script: str) -> dict:
    """Commit of the repository and whether the code that defines this step has uncommitted edits.

    git_dirty looks only at scripts/ and workflow/ (the code that runs the step), not at audit notes.
    Values that cannot be read are recorded as null, never guessed.
    """
    commit = _git(["rev-parse", "HEAD"])
    dirty = _git(["status", "--porcelain", "--", "scripts", "workflow"])
    return {
        "git_commit": commit,
        "git_dirty": None if dirty is None else bool(dirty),
        "script": script,
    }


def environment_info() -> dict:
    return {
        "hostname": socket.gethostname(),
        "python": platform.python_version(),
        "pandas": _pkg_version("pandas"),
        "yaml": _pkg_version("yaml"),
    }


def _pkg_version(name: str) -> str | None:
    try:
        mod = __import__(name)
        return getattr(mod, "__version__", None)
    except Exception:
        return None


def fingerprint(path: str | os.PathLike) -> dict:
    """path / size / mtime_ns / sha256. Files over 1 GB keep sha256 null with an explicit reason."""
    p = Path(path)
    st = p.stat()
    rec = {"path": str(p), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    if st.st_size > HASH_LIMIT_BYTES:
        rec["sha256"] = None
        rec["reason"] = "size>1GB"
        return rec
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    rec["sha256"] = h.hexdigest()
    return rec


def _safe(v):
    """Make values YAML- and JSON-serialisable without changing their meaning."""
    if v is None or isinstance(v, (bool, int, str)):
        return v
    if isinstance(v, float):
        return v
    if isinstance(v, (list, tuple)):
        return [_safe(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _safe(x) for k, x in v.items()}
    if isinstance(v, os.PathLike):
        return str(v)
    return str(v)


def record_path(results_dir: str | os.PathLike, rule: str, wildcard: str | None) -> Path:
    name = f"{rule}.{wildcard}.yaml" if wildcard else f"{rule}.yaml"
    return Path(results_dir) / "run_params" / name


def write_record(results_dir, rule, *, wildcard=None, wildcards=None, started_at, effective_params,
                 inputs=None, script, external_queries=None, extra=None) -> Path:
    """Write this step's record. Called at the END of the step, after its outputs exist.

    Raises on failure: a step must not report success without its record.
    """
    doc = {
        "schema_version": SCHEMA_VERSION,
        "provenance": "executed",
        "rule": rule,
        "wildcards": _safe(wildcards or ({"value": wildcard} if wildcard else {})),
        "started_at": started_at,
        "finished_at": now_iso(),
        "effective_params": _safe(effective_params),
        "code": code_info(script),
        "inputs": [fingerprint(p) for p in (inputs or []) if p and Path(p).exists()],
        "environment": environment_info(),
    }
    if external_queries is not None:
        doc["external_queries"] = _safe(external_queries)
    if extra:
        doc.update(_safe(extra))
    path = record_path(results_dir, rule, wildcard)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, allow_unicode=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return path


def _normalise(v):
    if isinstance(v, float):
        return float(f"{v:.12g}")
    if isinstance(v, dict):
        return {k: _normalise(x) for k, x in sorted(v.items())}
    if isinstance(v, (list, tuple)):
        return [_normalise(x) for x in v]
    return v


def canonical_key(record: dict, include_code_commit: bool = False) -> str:
    """Cache key for item L3. Only fields that change the step's output enter the key.

    In the key: rule, wildcards, effective_params, and each input's sha256 (or size and mtime_ns when the
    sha256 is null). Not in the key: started_at, finished_at, hostname, environment, external_queries
    dates, and the code commit unless include_code_commit is set (a repository commit changes on every
    commit, which would invalidate every cache entry; decide once real records exist).
    """
    key = {
        "rule": record.get("rule"),
        "wildcards": record.get("wildcards"),
        "effective_params": record.get("effective_params"),
        "inputs": [
            {"path": i.get("path"),
             "sha256": i.get("sha256"),
             "size": i.get("size") if i.get("sha256") is None else None,
             "mtime_ns": i.get("mtime_ns") if i.get("sha256") is None else None}
            for i in record.get("inputs", [])
        ],
    }
    if include_code_commit:
        key["code_commit"] = (record.get("code") or {}).get("git_commit")
    blob = json.dumps(_normalise(key), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def consolidate(results_dir: str | os.PathLike, project_id: str) -> Path:
    base = Path(results_dir) / "run_params"
    rules: dict[str, list] = {}
    for f in sorted(base.glob("*.yaml")):
        if f.name.startswith("."):
            continue
        doc = yaml.safe_load(open(f, encoding="utf-8"))
        if doc.get("provenance") != "executed":
            raise SystemExit(f"{f}: provenance is not 'executed'; refusing to consolidate")
        rules.setdefault(doc["rule"], []).append({"file": f.name, **doc})
    out = {
        "schema_version": SCHEMA_VERSION,
        "provenance": "executed",
        "project_id": project_id,
        "generated_at": now_iso(),
        "rules": rules,
    }
    path = Path(results_dir) / "run_params.yaml"
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        yaml.safe_dump(out, fh, sort_keys=False, allow_unicode=True)
    os.replace(tmp, path)
    return path


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "consolidate":
        print(consolidate(sys.argv[2], sys.argv[3]))
    else:
        raise SystemExit("usage: run_params.py consolidate <results_dir> <project_id>")
