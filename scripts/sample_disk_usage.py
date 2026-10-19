"""
sample_disk_usage.py -- lightweight background disk-usage sampler.

Appends one row every 60s to a TSV: timestamp, free space (GB) under a
reference path (e.g. the data_root filesystem), and total size (bytes) of
a specific sample's own working directory. Meant to run in the background
(`&`) alongside a real pipeline rule, started just before and killed just
after the rule's main command -- near-zero cost (one `shutil.disk_usage`
call and one directory walk per minute), and NOT something that should ever
be run as its own standalone analysis. Its only purpose is to turn the
current disk_gb resource estimate (workflow/rules/circrna.smk's
_DISK_GB_PER_SAMPLE_ESTIMATE = 200, itself reverse-engineered from an old
disk-full incident, not measured -- see docs/platform_security_notes.md)
into a real measurement, for free, the next time a real sample happens to
run through this rule.

Usage (from within a rule's shell: block, backgrounded and killed after the
real command finishes):
    python scripts/sample_disk_usage.py <work_dir> <log_tsv> [--ref-path PATH] [--interval 60] &
    SAMPLER_PID=$!
    <real command>
    kill $SAMPLER_PID 2>/dev/null || true
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


def _du_bytes(path: Path) -> int:
    """Total size of path's contents, via `du -sb` (fast, uses directory
    metadata rather than Python's own slower recursive stat() walk for
    large trees)."""
    try:
        out = subprocess.run(
            ["du", "-sb", str(path)], capture_output=True, text=True, timeout=30
        )
        if out.returncode == 0:
            return int(out.stdout.split()[0])
    except Exception:
        pass
    return -1  # sentinel: du failed or timed out, not "zero bytes"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work_dir", help="Directory to measure with du -sb (e.g. a sample's own outdir)")
    ap.add_argument("log_tsv", help="Path to append rows to")
    ap.add_argument("--ref-path", default=None,
                     help="Path to check free space under (default: work_dir's own filesystem)")
    ap.add_argument("--interval", type=int, default=60, help="Seconds between samples")
    args = ap.parse_args()

    work_dir = Path(args.work_dir)
    ref_path = Path(args.ref_path) if args.ref_path else work_dir
    log_path = Path(args.log_tsv)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    new_file = not log_path.exists()
    with open(log_path, "a") as f:
        if new_file:
            f.write("timestamp\tfree_gb\twork_dir_bytes\n")
        while True:
            free_gb = shutil.disk_usage(str(ref_path)).free / (1024 ** 3)
            size_b = _du_bytes(work_dir)
            f.write(f"{datetime.now().isoformat(timespec='seconds')}\t{free_gb:.1f}\t{size_b}\n")
            f.flush()
            time.sleep(args.interval)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(0)
