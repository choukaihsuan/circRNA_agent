"""
Download rules: SRA → paired-end FASTQ (gzipped)
Download priority:
  1. srapath + aria2c  (AWS S3 direct, fastest)
  2. ascp              (Aspera, fast but requires valid key)
  3. prefetch          (NCBI HTTPS, fallback)
"""

import shutil
import subprocess
import time
from pathlib import Path


def _find_tool(name):
    """Find a binary, preferring known-working envs over system PATH.

    Priority: sra_env > circrna > base conda > which() > /usr/local
    sra_env has fasterq-dump that works on CentOS 7 glibc 2.17.
    circrna has aria2c for fast multi-connection S3 downloads.
    """
    import shutil as sh, os, glob
    home = os.path.expanduser("~")
    bases = [
        f"{home}/miniconda3", f"{home}/miniforge3",
        f"{home}/miniconda",  f"{home}/anaconda3",
    ]
    # Check high-priority named envs first
    for prio_env in ("sra_env", "circrna", "ciriquant"):
        for base in bases:
            c = f"{base}/envs/{prio_env}/bin/{name}"
            if os.path.isfile(c) and os.access(c, os.X_OK):
                return c
    # Fall back to which() (uses current PATH)
    found = sh.which(name)
    if found:
        return found
    # Last resort: base conda bins + all envs
    candidates = []
    for base in bases:
        candidates.append(f"{base}/bin/{name}")
        for d in glob.glob(f"{base}/envs/*/bin/{name}"):
            candidates.append(d)
    candidates.append(f"/usr/local/bin/{name}")
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return name  # best-effort fallback


def _find_sra_tool(name):
    """Alias for SRA-toolkit tools."""
    return _find_tool(name)


def _find_ascp():
    """Return (ascp_path, key_path) or (None, None) if Aspera is not installed."""
    import shutil as sh
    ascp = sh.which("ascp")
    if not ascp:
        return None, None
    key = Path(ascp).parent.parent / "etc" / "asperaweb_id_dsa.openssh"
    if key.exists():
        return ascp, str(key)
    return None, None


def _sra_ascp_path(srr):
    """Build NCBI FTP path for ascp."""
    prefix = srr[:6]
    return f"/sra/sra-instant/reads/ByRun/sra/SRR/{prefix}/{srr}/{srr}.sra"


rule download_fastq:
    """Fetch SRA accession and convert to gzipped FASTQ pair."""
    output:
        r1 = protected(RAW_DIR + "/{srr}_1.fastq.gz"),
        r2 = protected(RAW_DIR + "/{srr}_2.fastq.gz"),
    params:
        sra_cache       = config.get("download", {}).get("sra_cache_dir",
                              str(Path(RAW_DIR).parent / "sra_cache")),
        tmp_dir         = config.get("download", {}).get("tmp_dir",
                              str(Path(RAW_DIR).parent / "sra_tmp")) + "/{srr}",
        out_dir         = RAW_DIR,
        retry           = config.get("download", {}).get("retry", 3),
        ascp_speed      = config.get("download", {}).get("ascp_speed", "500m"),
        s3_url_retry    = config.get("download", {}).get("s3_url_retry", 5),
        s3_retry_wait   = config.get("download", {}).get("s3_retry_wait_sec", 300),
    threads: 4   # 24 cores ÷ 4 = 6 parallel jobs × 8 conn = 48 S3 connections max
    log: "logs/download/{srr}.log"
    run:
        srr       = wildcards.srr
        out_dir   = Path(params.out_dir)
        tmp_dir   = Path(params.tmp_dir)
        sra_cache = Path(params.sra_cache)
        sra_dir   = sra_cache / srr
        sra_file  = sra_dir / f"{srr}.sra"
        log_path  = Path(log[0])

        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        sra_dir.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        def run_cmd(cmd, logf):
            logf.write(f"[CMD] {' '.join(str(c) for c in cmd)}\n")
            logf.flush()
            result = subprocess.run(cmd, stdout=logf, stderr=logf, text=True)
            return result.returncode

        with open(log_path, "w") as logf:
            downloaded = False

            # ── 1. srapath + aria2c/curl (AWS S3 direct) ─────────
            import shutil as sh
            aria2c  = _find_tool("aria2c")   # searches miniconda3 + miniforge3 envs
            srapath = _find_sra_tool("srapath")
            curl    = _find_tool("curl")
            aria2_meta = sra_dir / f"{srr}.sra.aria2"
            sra_partial = aria2_meta.exists()  # incomplete aria2c download
            # Short-circuit: complete .sra already present — skip all network downloads
            if sra_file.exists() and not sra_partial:
                downloaded = True
                logf.write(f"[skip] SRA already downloaded: {sra_file}\n")
            if not downloaded and srapath and (not sra_file.exists() or sra_partial):
                # Retry S3 URL lookup: S3 is often temporarily unavailable but recovers quickly.
                # Retry up to s3_url_retry times with s3_retry_wait seconds between attempts
                # before giving up and falling through to ascp/prefetch.
                s3_max  = int(params.s3_url_retry)
                s3_wait = int(params.s3_retry_wait)
                s3_url  = ""
                for s3_attempt in range(s3_max + 1):
                    logf.write(f"[s3] Getting S3 URL for {srr}"
                               f" (attempt {s3_attempt + 1}/{s3_max + 1})\n")
                    logf.flush()
                    result = subprocess.run(
                        [srapath, "--location", "s3", srr],
                        capture_output=True, text=True
                    )
                    s3_url = result.stdout.strip()
                    if s3_url.startswith("http"):
                        logf.write(f"[s3] URL: {s3_url}\n")
                        break
                    logf.write(f"[WARN] srapath returned no S3 URL: '{s3_url}'\n")
                    if s3_attempt < s3_max:
                        logf.write(f"[s3] Waiting {s3_wait}s before retry…\n")
                        logf.flush()
                        time.sleep(s3_wait)
                    else:
                        logf.write(f"[WARN] S3 URL unavailable after {s3_max + 1} attempts,"
                                   f" falling back to ascp/prefetch\n")

                if s3_url.startswith("http"):
                    logf.flush()
                    if aria2c:
                        rc = run_cmd([
                            aria2c,
                            "-x", "8", "-s", "8", "-k", "10M",
                            "--file-allocation=none",
                            "--retry-wait=10", "--max-tries=5",
                            "-c",  # resume partial download via .aria2 metadata
                            "-d", str(sra_dir),
                            "-o", f"{srr}.sra",
                            s3_url,
                        ], logf)
                    elif curl:
                        rc = run_cmd([
                            curl, "-L", "--retry", "3", "--retry-delay", "5",
                            "-o", str(sra_file), s3_url,
                        ], logf)
                    else:
                        rc = 1
                    # Trust rc=0 from aria2c/curl; NFS stat() cache can lag
                    if rc == 0:
                        downloaded = True
                        logf.write(f"[s3] Download complete\n")
                    else:
                        logf.write(f"[WARN] S3 download failed (rc={rc}), trying ascp\n")

            # ── 2. ascp (Aspera) ─────────────────────────────────
            if not downloaded:
                ascp, key = _find_ascp()
                if ascp and not sra_file.exists():
                    remote = _sra_ascp_path(srr)
                    logf.write(f"[ascp] Downloading {srr} via Aspera\n")
                    logf.flush()
                    rc = run_cmd([
                        ascp, "-i", key, "-k", "1", "-T",
                        "-l", str(params.ascp_speed),
                        f"anonftp@ftp.ncbi.nlm.nih.gov:{remote}",
                        str(sra_dir) + "/",
                    ], logf)
                    if rc == 0 and sra_file.exists():
                        downloaded = True
                        logf.write(f"[ascp] Download complete\n")
                    else:
                        logf.write(f"[WARN] ascp failed (rc={rc}), falling back to prefetch\n")

            # ── 3. prefetch (NCBI HTTPS fallback) ────────────────
            if not downloaded:
                # Remove stale lock files from previous interrupted downloads
                # Use try/except for Python 3.7 compatibility (missing_ok added in 3.8)
                for lock in sra_dir.glob("*.lock"):
                    try:
                        lock.unlink()
                    except OSError:
                        pass
                for tmp in sra_dir.glob("*.tmp"):
                    try:
                        tmp.unlink()
                    except OSError:
                        pass

                prefetch_bin = _find_sra_tool("prefetch")
                for attempt in range(1, int(params.retry) + 1):
                    logf.write(f"[prefetch] Attempt {attempt} for {srr} ({prefetch_bin})\n")
                    logf.flush()
                    rc = run_cmd([prefetch_bin, srr, "-O", str(sra_cache)], logf)
                    if rc == 0:
                        downloaded = True
                        break
                    if attempt < int(params.retry):
                        time.sleep(30)
                if not downloaded:
                    logf.write(f"[WARN] prefetch failed all attempts; "
                               f"will try fasterq-dump direct download\n")

            # ── fasterq-dump: use SRA file if available, else direct SRR download ──
            fastq_stage = tmp_dir / "fastq_stage"
            fastq_stage.mkdir(parents=True, exist_ok=True)
            fasterq_input = str(sra_file) if sra_file.exists() else srr
            fasterq_bin = _find_tool("fasterq-dump")
            rc = run_cmd(
                [fasterq_bin, fasterq_input,
                 "-O", str(fastq_stage),
                 "-t", str(tmp_dir),
                 "--split-files",
                 "--skip-technical",
                 "--threads", str(threads)],
                logf,
            )
            if rc != 0:
                raise RuntimeError(f"fasterq-dump failed for {srr}. See {log_path}")

            # ── Compress with pigz (fallback: gzip), then move to out_dir ──
            pigz_bin = _find_tool("pigz")
            for suffix in ("_1.fastq", "_2.fastq"):
                fastq = fastq_stage / f"{srr}{suffix}"
                if not fastq.exists():
                    raise RuntimeError(f"Expected file not found: {fastq}")
                rc = 1
                try:
                    rc = run_cmd([pigz_bin, "-p", str(threads), str(fastq)], logf)
                except (FileNotFoundError, OSError):
                    logf.write("[WARN] pigz not available, falling back to gzip\n")
                if rc != 0:
                    rc = run_cmd(["gzip", str(fastq)], logf)
                if rc != 0:
                    raise RuntimeError(f"Compression failed for {fastq}")
                gz = fastq_stage / f"{srr}{suffix}.gz"
                shutil.move(str(gz), str(out_dir / gz.name))
                logf.write(f"[move] {gz.name} → {out_dir}\n")

            # ── Cleanup ──────────────────────────────────────────
            shutil.rmtree(tmp_dir, ignore_errors=True)
            shutil.rmtree(sra_dir, ignore_errors=True)
            logf.write(f"[OK] {srr} downloaded and compressed\n")
