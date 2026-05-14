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
        sra_cache  = config["download"]["sra_cache_dir"],
        tmp_dir    = config["download"]["tmp_dir"] + "/{srr}",
        out_dir    = RAW_DIR,
        retry      = config["download"]["retry"],
        ascp_speed = config["download"].get("ascp_speed", "500m"),
    threads: min(config["threads"], 8)
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

            # ── 1. srapath + aria2c (AWS S3 direct) ─────────────
            import shutil as sh
            aria2c = sh.which("aria2c")
            srapath = sh.which("srapath")
            if aria2c and srapath and not sra_file.exists():
                logf.write(f"[aria2c] Getting S3 URL for {srr}\n")
                logf.flush()
                result = subprocess.run(
                    [srapath, srr], capture_output=True, text=True
                )
                s3_url = result.stdout.strip()
                if s3_url.startswith("http"):
                    logf.write(f"[aria2c] Downloading from {s3_url}\n")
                    logf.flush()
                    rc = run_cmd([
                        aria2c,
                        "-x", "16", "-s", "16", "-k", "10M",
                        "--file-allocation=none",
                        "--retry-wait=5", "--max-tries=3",
                        "-d", str(sra_dir),
                        "-o", f"{srr}.sra",
                        s3_url,
                    ], logf)
                    if rc == 0 and sra_file.exists():
                        downloaded = True
                        logf.write(f"[aria2c] Download complete\n")
                    else:
                        logf.write(f"[WARN] aria2c failed (rc={rc}), trying ascp\n")
                else:
                    logf.write(f"[WARN] srapath returned no S3 URL: {s3_url}\n")

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
                for lock in sra_dir.glob("*.lock"):
                    lock.unlink(missing_ok=True)
                for tmp in sra_dir.glob("*.tmp"):
                    tmp.unlink(missing_ok=True)

                for attempt in range(1, int(params.retry) + 1):
                    logf.write(f"[prefetch] Attempt {attempt} for {srr}\n")
                    logf.flush()
                    rc = run_cmd(["prefetch", srr, "-O", str(sra_cache)], logf)
                    if rc == 0:
                        downloaded = True
                        break
                    if attempt == int(params.retry):
                        raise RuntimeError(
                            f"All download methods failed after {params.retry} attempts. "
                            f"See {log_path}"
                        )
                    time.sleep(30)

            # ── fasterq-dump (output to Linux fs to avoid NTFS space check) ──
            fastq_stage = tmp_dir / "fastq_stage"
            fastq_stage.mkdir(parents=True, exist_ok=True)
            fasterq_input = str(sra_file) if sra_file.exists() else srr
            rc = run_cmd(
                ["fasterq-dump", fasterq_input,
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
            for suffix in ("_1.fastq", "_2.fastq"):
                fastq = fastq_stage / f"{srr}{suffix}"
                if not fastq.exists():
                    raise RuntimeError(f"Expected file not found: {fastq}")
                rc = run_cmd(["pigz", "-p", str(threads), str(fastq)], logf)
                if rc != 0:
                    logf.write("[WARN] pigz not found, falling back to gzip\n")
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
