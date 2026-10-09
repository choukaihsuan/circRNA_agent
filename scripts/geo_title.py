"""GEO series title lookup (NCBI E-utilities esummary, db=gds).

Kept free of Flask/pandas so both the Web UI and a stand-alone audit script can use it.
"""
from __future__ import annotations

import json
import re
import urllib.request
from typing import Optional

_GSE_RE = re.compile(r"^GSE(\d+)\Z")


def geo_uid(accession: str) -> Optional[str]:
    """GEO DataSets UID for a GSE accession: ``200000000 + <number>``.

    The old code used ``"200" + <number>``, which is only right for 6-digit series
    (GSE113230 -> 200113230). For GSE58135 it produced 20058135 - a *different* record -
    so the report got somebody else's study title. Correct: GSE58135 -> 200058135.
    """
    m = _GSE_RE.match((accession or "").strip().upper())
    if not m:
        return None
    return str(200000000 + int(m.group(1)))


def fetch_geo_title(accession: str, timeout: int = 8) -> str:
    """Return the study title, or "" on any failure or when the record is not that accession."""
    uid = geo_uid(accession)
    if not uid:
        return ""
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
           "?db=gds&id=%s&retmode=json" % uid)
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read())
        rec = data.get("result", {}).get(uid, {})
        # Guard against attaching another series' title: the record must be this accession.
        acc = str(rec.get("accession", "")).upper()
        if acc and acc != accession.strip().upper():
            return ""
        return str(rec.get("title", ""))
    except Exception:
        return ""
