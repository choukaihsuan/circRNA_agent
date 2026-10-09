"""
Input-validation and path-safety helpers for the Web UI.

Pure functions, no Flask dependency, so they can be unit-tested without a server.
Must stay Python 3.7 compatible (the HPC server's conda env runs 3.7):
no walrus operator, no ``Path.is_relative_to`` (3.9).
"""

from __future__ import annotations

import os
import re
import secrets
from pathlib import Path
from typing import Iterable, Optional

# Accession whitelists (OWASP ASVS V5.1: positive validation)
ACCESSION_RE = re.compile(r"^(GSE|SRP|PRJNA)\d+\Z")
SRR_RE = re.compile(r"^SRR\d+\Z")
# Free-form project ids used by manual/local runs ("CUSTOM", "LOCAL", user labels).
# Must start with a letter so it can never look like a CLI option ("-x").
CUSTOM_ID_RE = re.compile(r"^[A-Z][A-Z0-9_-]{0,31}\Z")
# Local FASTQ sample names (become file names under raw_dir)
SAMPLE_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}\Z")
# Job ids: "<project>-<suffix>"; old ids were 4 chars, new ones are url-safe tokens
JOB_ID_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}\Z")
# Group / condition labels
LABEL_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _.+-]{0,63}\Z")

DE_METHODS = ("edgeR_ciriquant", "deseq2", "limma")
TOOL_NAMES = ("ciriquant", "dcc")

# Directories a user-supplied symlink target must never resolve into.
FORBIDDEN_ROOTS = tuple(
    Path(p) for p in ("/etc", "/proc", "/sys", "/dev", "/root", "/boot",
                      "/var/lib", "/var/log", "/run", "/bin", "/sbin", "/usr")
)


class ValidationError(ValueError):
    """Raised on rejected input. Message is deliberately generic: never echo input."""


def validate_project_id(value: str) -> str:
    """Return the normalised (upper-case) project id or raise ValidationError."""
    v = (value or "").strip().upper()
    if ACCESSION_RE.match(v) or CUSTOM_ID_RE.match(v):
        return v
    raise ValidationError("invalid project id")


def validate_accession(value: str) -> str:
    """Strict GSE/SRP/PRJNA accession only (no free-form ids)."""
    v = (value or "").strip().upper()
    if ACCESSION_RE.match(v):
        return v
    raise ValidationError("invalid accession")


def validate_srr(value: str) -> str:
    v = (value or "").strip().upper()
    if SRR_RE.match(v):
        return v
    raise ValidationError("invalid run accession")


def validate_sample_name(value: str) -> str:
    v = (value or "").strip()
    if SAMPLE_NAME_RE.match(v) and ".." not in v:
        return v
    raise ValidationError("invalid sample name")


def validate_label(value: str, default: str) -> str:
    v = (value or "").strip() or default
    if LABEL_RE.match(v):
        return v
    raise ValidationError("invalid label")


def validate_job_id(value: str) -> str:
    v = (value or "").strip()
    if JOB_ID_RE.match(v):
        return v
    raise ValidationError("invalid job id")


def bounded_int(value, default: int, lo: int, hi: int) -> int:
    """Parse an int from user input, clamp to [lo, hi]; fall back to default on garbage."""
    try:
        n = int(str(value).strip())
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, n))


def bounded_float(value, default: float, lo: float, hi: float) -> float:
    try:
        x = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    if x != x or x in (float("inf"), float("-inf")):  # NaN / inf
        return default
    return max(lo, min(hi, x))


def _within(child: Path, root: Path) -> bool:
    try:
        child.relative_to(root)
        return True
    except ValueError:
        return False


def safe_join(root: os.PathLike, *parts: str) -> Path:
    """Join ``parts`` under ``root`` and prove the result stays inside ``root``.

    Rejects absolute components, NUL bytes, '..' and symlinks that escape root.
    """
    root_p = Path(root).resolve()
    for part in parts:
        s = str(part)
        if "\x00" in s or os.path.isabs(s) or ".." in Path(s).parts:
            raise ValidationError("path rejected")
    candidate = root_p.joinpath(*[str(p) for p in parts]).resolve()
    if not _within(candidate, root_p):
        raise ValidationError("path rejected")
    return candidate


def safe_target(target: str, allowed_roots: Iterable[os.PathLike]) -> Path:
    """Resolve a user-supplied absolute path and require it to live under one of
    ``allowed_roots`` (compared by path components, not string prefix) and not under
    a system directory. Returns the resolved path."""
    if not target or "\x00" in target:
        raise ValidationError("path rejected")
    resolved = Path(target).resolve()
    for bad in FORBIDDEN_ROOTS:
        if _within(resolved, bad):
            raise ValidationError("path rejected")
    for root in allowed_roots:
        if root and _within(resolved, Path(root).resolve()):
            return resolved
    raise ValidationError("path rejected")


def new_job_id(project_id: str) -> str:
    """Unguessable job id: project id + 72 bits of CSPRNG (url-safe, 12 chars)."""
    return "{}-{}".format(project_id, secrets.token_urlsafe(9))


def mask_token(token: str) -> str:
    """Show only first/last 4 chars of a secret for logs."""
    if not token or len(token) <= 8:
        return "****"
    return "{}…{}".format(token[:4], token[-4:])


def mask_link(link: str) -> str:
    """Mask the token segment of a magic link (/auth/<token>) for log output."""
    return re.sub(r"(/auth/)([^/?#]+)", lambda m: m.group(1) + mask_token(m.group(2)), link)


def allowed_emails(env: Optional[str]) -> Optional[set]:
    """Parse PIPELINE_ALLOWED_EMAILS. None => not configured (no enforcement)."""
    if not env or not env.strip():
        return None
    return {e.strip().lower() for e in env.split(",") if e.strip()}


EMAIL_RE = re.compile(r"^[A-Za-z0-9._%+-]{1,64}@[A-Za-z0-9.-]{1,190}\.[A-Za-z]{2,24}\Z")


def validate_email(value: str) -> str:
    """Single, header-injection-safe address (no CR/LF/space/comma)."""
    v = (value or "").strip()
    if len(v) <= 254 and EMAIL_RE.match(v):
        return v.lower()
    raise ValidationError("invalid email")
