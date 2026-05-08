"""
Download rules: SRA → paired-end FASTQ (gzipped)
Uses Python run: block for cross-platform compatibility (Windows dev + Linux HPC).
"""

import shutil
import subprocess
import time
from pathlib import Path


rule download_fastq:
    """Fetch SRA accession and convert to gzipped FASTQ pair."""
    output:
        r1 = protected("data/raw_fastq/{srr}_1.fastq.gz"),
        r2 = protected("data/raw_fastq/{srr}_2.fastq.gz"),
    params:
        sra_cache = "data/sra",
        tmp_dir   = "data/tmp/{srr}",
        out_dir   = "data/raw_fastq",
        retry     = config["download"]["retry"],
    threads: min(config["threads"], 8)
    log: "logs/download/{srr}.log"
    run:
        srr       = wildcards.srr
        out_dir   = Path(params.out_dir)
        tmp_dir   = Path(params.tmp_dir)
        sra_cache = Path(params.sra_cache)
        log_path  = Path(log[0])

        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        sra_cache.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        def run_cmd(cmd, logf):
            logf.write(f"[CMD] {' '.join(str(c) for c in cmd)}\n")
            logf.flush()
            result = subprocess.run(cmd, stdout=logf, stderr=logf, text=True)
            return result.returncode

        with open(log_path, "w") as logf:
            # ── prefetch with retry ──────────────────────────────
            for attempt in range(1, int(params.retry) + 1):
                logf.write(f"[Attempt {attempt}] prefetch {srr}\n")
                logf.flush()
                rc = run_cmd(["prefetch", srr, "-O", str(sra_cache)], logf)
                if rc == 0:
                    break
                if attempt == int(params.retry):
                    raise RuntimeError(
                        f"prefetch failed after {params.retry} attempts. "
                        f"See {log_path}"
                    )
                time.sleep(30)

            # ── fasterq-dump ─────────────────────────────────────
            rc = run_cmd(
                ["fasterq-dump", srr,
                 "-O", str(out_dir),
                 "-t", str(tmp_dir),
                 "--split-files",
                 "--skip-technical",
                 "--threads", str(threads)],
                logf,
            )
            if rc != 0:
                raise RuntimeError(f"fasterq-dump failed for {srr}. See {log_path}")

            # ── Compress with pigz (fallback: gzip) ──────────────
            for suffix in ("_1.fastq", "_2.fastq"):
                fastq = out_dir / f"{srr}{suffix}"
                if not fastq.exists():
                    raise RuntimeError(f"Expected file not found: {fastq}")
                # try pigz first (parallel), fall back to gzip
                rc = run_cmd(["pigz", "-p", str(threads), str(fastq)], logf)
                if rc != 0:
                    logf.write("[WARN] pigz not found, falling back to gzip\n")
                    rc = run_cmd(["gzip", str(fastq)], logf)
                if rc != 0:
                    raise RuntimeError(f"Compression failed for {fastq}")

            # ── Cleanup ──────────────────────────────────────────
            shutil.rmtree(tmp_dir, ignore_errors=True)
            shutil.rmtree(sra_cache / srr, ignore_errors=True)
            logf.write(f"[OK] {srr} downloaded and compressed\n")
