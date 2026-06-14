"""
predict_interactions.py
-----------------------
Query two sources for miRNA/RBP interactions of top DE circRNAs:
  1. CircInteractome (NIH)  — TargetScan predictions, uses circBase ID, hg19
  2. ENCORI (CLIP-seq)      — experimental validation, uses gene symbol, hg38

Generates interactions.json embedded into the HTML report.

Usage:
    python scripts/predict_interactions.py \
        --de    results/de/de_results.tsv \
        --iso   results/circRNA/isoform_groups.tsv \
        --out   results/de/interactions.json \
        --top-n 50 --clip-exp-num 1 --program-num 2
"""

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

try:
    import requests
    _REQUESTS = True
except ImportError:
    _REQUESTS = False
    print("[predict] 'requests' not installed — API queries skipped", file=sys.stderr)

CIRCINTERACTOME_CIRC  = "https://circinteractome.nia.nih.gov/api/v2/circsearch"
CIRCINTERACTOME_MIRNA = "https://circinteractome.nia.nih.gov/api/v2/mirnasearch"
CIRCINTERACTOME_RBP   = "https://circinteractome.nia.nih.gov/api/v2/rbpsearch"

ENCORI_MIRNA = "https://rnasysu.com/encori/api/miRNATarget/"
ENCORI_RBP   = "https://rnasysu.com/encori/api/RBPTarget/"

_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Referer":    "https://circinteractome.nia.nih.gov/circular_rna.html",
}
_ENCORI_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "Accept":     "application/json",
}

# CircInteractome rbpsearch accepts uppercase RBP names
_RBP_NAME_MAP = {
    "HuR": "HUR", "huR": "HUR", "hur": "HUR",
    "Lin28a": "LIN28A", "Lin28b": "LIN28B",
    "Fmrp": "FMRP", "fmrp": "FMRP",
}
_VALID_RBPS = {
    "AGO1","AGO2","AGO3","ALKBH5","AUF1","C17ORF85","C22ORF28","CAPRIN1",
    "DGCR8","EIF4A3","EWSR1","FMRP","FOX2","FUS","FXR1","FXR2","HNRNPC",
    "HUR","IGF2BP1","IGF2BP2","IGF2BP3","LIN28A","LIN28B","METTL3",
    "MOV10","PTB","PUM2","QKI","SFRS1","TAF15","TDP43","TIA1","TIAL1",
    "TNRC6","U2AF65","UPF1","WTAP","YTHDF1","YTHDF2","YTHDF3","ZC3H7B",
}


# ── HTML parsing helpers ──────────────────────────────────────────────────────

def _strip_tags(html: str) -> str:
    return re.sub(r"<[^>]+>", "", html).replace("&nbsp;", " ").strip()


def _parse_spliced_length(html: str) -> int:
    """Parse spliced sequence length (nt) from CircInteractome circsearch HTML."""
    m = re.search(
        r"Spliced Seq Length.*?<TD[^>]*><FONT[^>]*>\s*(\d+)\s*bp", html, re.DOTALL
    )
    return int(m.group(1)) if m else 0


def _parse_rbp_counts(html: str) -> List[dict]:
    """Extract RBP binding site counts from CircInteractome circsearch HTML."""
    rbp_list = []

    # RBP inside circRNA — from Google Visualization JS setCell calls
    internal = {}
    for m in re.finditer(
        r"data\.setCell\s*\(\s*(\d+)\s*,\s*0\s*,\s*'([^']+)'\s*\)",
        html,
    ):
        internal[m.group(1)] = m.group(2)
    for m in re.finditer(
        r"data\.setCell\s*\(\s*(\d+)\s*,\s*1\s*,\s*(\d+)\s*\)",
        html,
    ):
        idx = m.group(1)
        if idx in internal:
            rbp_list.append({
                "RBPName":      internal[idx],
                "bindingSites": int(m.group(2)),
                "location":     "internal",
                "sites":        [],  # will be filled by _fetch_rbp_positions
                "source":       "CircInteractome",
            })

    # RBP in flanking regions — from HTML table
    flanking_match = re.search(
        r"RNA-binding protein sites matching flanking regions.*?</TABLE>",
        html, re.DOTALL,
    )
    if flanking_match:
        rows = re.findall(
            r"<TR[^>]*>(.*?)</TR>", flanking_match.group(0), re.DOTALL
        )
        for row in rows[2:]:  # skip header rows
            cells = re.findall(r"<TD[^>]*>(.*?)</TD>", row, re.DOTALL)
            if len(cells) >= 2:
                rbp  = _strip_tags(cells[0])
                tags = _strip_tags(cells[1])
                if rbp:
                    rbp_list.append({
                        "RBPName":      rbp,
                        "bindingSites": int(tags) if tags.isdigit() else 0,
                        "location":     "flanking",
                        "sites":        [],
                        "source":       "CircInteractome",
                    })

    return rbp_list


def _fetch_rbp_positions(circbase_id: str, rbp_name: str) -> List[dict]:
    """Query rbpsearch API for individual binding site positions within the circRNA."""
    if not _REQUESTS:
        return []
    # Normalise name to API enum format
    api_name = _RBP_NAME_MAP.get(rbp_name, rbp_name.upper())
    if api_name not in _VALID_RBPS:
        api_name = rbp_name.upper()
    if api_name not in _VALID_RBPS:
        return []

    h = dict(_HEADERS)
    h["Referer"] = "https://circinteractome.nia.nih.gov/rna_binding_protein.html"
    try:
        r = requests.get(
            CIRCINTERACTOME_RBP,
            params={
                "circular_rna_query": circbase_id,
                "rbp_query":          api_name,
                "output_type":        "Web",
            },
            headers=h, timeout=30,
        )
        if r.status_code != 200:
            return []
    except Exception:
        return []

    sites = []
    html = r.text
    rows = re.findall(r"<TR[^>]*>(.*?)</TR>", html, re.DOTALL)
    for row in rows[1:]:  # skip header
        cells = re.findall(r"<TD[^>]*>(.*?)</TD>", row, re.DOTALL)
        cleaned = [_strip_tags(c) for c in cells]
        if len(cleaned) >= 10:
            try:
                circ_s = int(cleaned[8])
                circ_e = int(cleaned[9])
                sites.append({
                    "circ_start": circ_s,
                    "circ_end":   circ_e,
                    "circ_pos":   f"{circ_s}–{circ_e}",
                })
            except (ValueError, IndexError):
                pass
    time.sleep(0.3)
    return sites


def _parse_mirna(html: str) -> List[dict]:
    """Extract TargetScan miRNA predictions from CircInteractome mirnasearch HTML."""
    mirna_list = []

    idx = html.find("TargetScan miRNA predictions")
    if idx < 0:
        return mirna_list
    end = html.find("</TABLE>", idx)
    segment = html[idx : end + 8]

    rows = re.findall(
        r"<TR[^>]*>\s*((?:<TD[^>]*>.*?</TD>\s*)+)</TR>", segment, re.DOTALL
    )
    for row in rows:
        cells = re.findall(r"<TD[^>]*>(.*?)</TD>", row, re.DOTALL)
        if len(cells) < 5:
            continue
        mirna_match = re.search(r"(hsa-miR-[\w-]+|hsa-let-[\w-]+)", cells[0])
        if not mirna_match:
            continue
        mirna_name = mirna_match.group(1)
        site_type  = _strip_tags(cells[2])
        pos_start  = _strip_tags(cells[3])
        pos_end    = _strip_tags(cells[4])
        circ_pos   = f"{pos_start}–{pos_end}" if pos_start and pos_end else "—"
        mirna_list.append({
            "miRNAName": mirna_name,
            "siteType":  site_type,
            "circ_pos":  circ_pos,
            "source":    "CircInteractome",
        })

    return mirna_list


# ── CircInteractome API ───────────────────────────────────────────────────────

def _query_circinteractome(
    circbase_id: str,
) -> Tuple[List[dict], List[dict], int]:
    """Return (mirna_list, rbp_list_with_positions, spliced_length)."""
    if not _REQUESTS or not circbase_id:
        return [], [], 0

    mirna: List[dict] = []
    rbp_list: List[dict] = []
    spliced_length = 0

    # ── 1. circsearch: RBP counts + spliced_length ──
    try:
        r = requests.get(
            CIRCINTERACTOME_CIRC,
            params={"circular_rna_query": circbase_id},
            headers=_HEADERS, timeout=30,
        )
        if r.status_code == 200:
            rbp_list      = _parse_rbp_counts(r.text)
            spliced_length = _parse_spliced_length(r.text)
    except Exception as exc:
        print(f"[predict] circsearch error: {exc}", file=sys.stderr)
    time.sleep(0.4)

    # ── 2. rbpsearch: fetch per-site positions for internal RBPs ──
    for item in rbp_list:
        if item["location"] == "internal":
            sites = _fetch_rbp_positions(circbase_id, item["RBPName"])
            item["sites"] = sites
            if sites:
                # build readable circ_pos summary
                item["circ_pos"] = "; ".join(s["circ_pos"] for s in sites[:5])
            else:
                item["circ_pos"] = "—"
        else:
            item["circ_pos"] = "flanking"

    # ── 3. mirnasearch ──
    try:
        h2 = dict(_HEADERS)
        h2["Referer"] = "https://circinteractome.nia.nih.gov/mirna_target_sites.html"
        r = requests.get(
            CIRCINTERACTOME_MIRNA,
            params={"circular_rna_query": circbase_id, "mirna_query": ""},
            headers=h2, timeout=30,
        )
        if r.status_code == 200:
            mirna = _parse_mirna(r.text)
    except Exception as exc:
        print(f"[predict] mirnasearch error: {exc}", file=sys.stderr)
    time.sleep(0.4)

    return mirna, rbp_list, spliced_length


# ── ENCORI API ───────────────────────────────────────────────────────────────

def _parse_encori_tsv(text: str) -> List[dict]:
    """Parse ENCORI TSV response (strips #-comment lines, returns list of dicts)."""
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    result = []
    for line in lines[1:]:
        cells = line.split("\t")
        result.append(dict(zip(header, cells)))
    return result


def _fetch_encori_mirna(gene_name: str, clip_exp_num: int = 1,
                        program_num: int = 2,
                        lo=None, coord_idx: Optional[dict] = None,
                        exon_nums: Optional[List[int]] = None) -> List[dict]:
    """Query ENCORI miRNATarget (TSV) for circRNA host gene. Returns one entry per unique miRNA."""
    if not _REQUESTS or not gene_name or gene_name in ("", "nan"):
        return []
    try:
        r = requests.get(
            ENCORI_MIRNA,
            params={
                "assembly":       "hg38",
                "geneType":       "circRNA",
                "miRNA":          "all",
                "target":         gene_name,
                "clipExpNum":     clip_exp_num,
                "programNum":     program_num,
                "pancancerNum":   0,
                "cellType":       0,
                "interactionNum": 0,
            },
            headers=_ENCORI_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        rows = _parse_encori_tsv(r.text)
    except Exception as exc:
        print(f"[predict] ENCORI miRNA error ({gene_name}): {exc}", file=sys.stderr)
        return []

    # aggregate by miRNA name: keep max clipExpNum, collect distinct cell types + positions
    agg: dict = {}
    for row in rows:
        name = row.get("miRNAname", "").strip()
        if not name:
            continue
        clip  = int(row.get("clipExpNum") or 0)
        cell  = row.get("cellline/tissue", "").strip()
        site  = row.get("merClass", "CLIP-seq").strip()
        chrom = row.get("chromosome", "").strip()
        start = row.get("start", "").strip()
        end   = row.get("end", "").strip()
        pos   = f"{chrom}:{start}-{end}" if chrom and start and end else "—"
        # hg38 coords as tuple for liftover (0-based half-open)
        hg38  = (chrom, int(start) - 1, int(end)) if (chrom and start.isdigit() and end.isdigit()) else None
        if name not in agg:
            agg[name] = {"clip": clip, "cells": set(), "site": site, "pos": pos, "hg38": hg38}
        else:
            agg[name]["clip"] = max(agg[name]["clip"], clip)
        if cell:
            agg[name]["cells"].add(cell)

    _circ_pos_re = re.compile(r'^\d+[–\-]\d+$')
    mirna_list = []
    for name, info in agg.items():
        circ_pos = info["pos"]
        if lo is not None and coord_idx is not None and exon_nums and info["hg38"]:
            hg38 = info["hg38"]
            lifted = _liftover_interval(lo, hg38[0], hg38[1], hg38[2])
            if lifted:
                mapped = _map_to_circ_pos(lifted[0], lifted[1], lifted[2],
                                          gene_name, exon_nums, coord_idx)
                if mapped:
                    circ_pos = mapped
        in_circ = bool(_circ_pos_re.match(circ_pos))
        mirna_list.append({
            "miRNAName":  name,
            "siteType":   info["site"],
            "circ_pos":   circ_pos,
            "clipExpNum": info["clip"],
            "cellType":   ", ".join(sorted(info["cells"])[:3]),
            "source":     "ENCORI",
            "in_circ":    in_circ,
        })
    time.sleep(0.3)
    return mirna_list


def _fetch_encori_rbp(gene_name: str, clip_exp_num: int = 1,
                      lo=None, coord_idx: Optional[dict] = None,
                      exon_nums: Optional[List[int]] = None) -> List[dict]:
    """Query ENCORI RBPTarget (TSV) for circRNA host gene. Returns one entry per unique RBP."""
    if not _REQUESTS or not gene_name or gene_name in ("", "nan"):
        return []
    try:
        r = requests.get(
            ENCORI_RBP,
            params={
                "assembly":       "hg38",
                "geneType":       "circRNA",
                "RBP":            "all",
                "target":         gene_name,
                "clipExpNum":     clip_exp_num,
                "programNum":     0,
                "pancancerNum":   0,
                "cellType":       0,
                "interactionNum": 0,
            },
            headers=_ENCORI_HEADERS,
            timeout=30,
        )
        if r.status_code != 200:
            return []
        rows = _parse_encori_tsv(r.text)
    except Exception as exc:
        print(f"[predict] ENCORI RBP error ({gene_name}): {exc}", file=sys.stderr)
        return []

    # aggregate by RBP name: max clipExpNum, collect cell types + positions
    agg: dict = {}
    for row in rows:
        name  = row.get("RBP", "").strip()
        if not name:
            continue
        clip  = int(row.get("clipExpNum") or 0)
        cell  = row.get("cellline/tissue", "").strip()
        chrom = row.get("chromosome", "").strip()
        start = row.get("narrowStart", "").strip()
        end   = row.get("narrowEnd",   "").strip()
        pos   = f"{chrom}:{start}-{end}" if chrom and start and end else "—"
        hg38  = (chrom, int(start) - 1, int(end)) if (chrom and start.isdigit() and end.isdigit()) else None
        if name not in agg:
            agg[name] = {"clip": clip, "cells": set(), "pos": pos, "hg38": hg38}
        else:
            agg[name]["clip"] = max(agg[name]["clip"], clip)
        if cell:
            agg[name]["cells"].add(cell)

    _circ_pos_re2 = re.compile(r'^(\d+)[–\-](\d+)$')
    rbp_list = []
    for name, info in agg.items():
        circ_pos = info["pos"]
        sites: List[dict] = []
        if lo is not None and coord_idx is not None and exon_nums and info["hg38"]:
            hg38   = info["hg38"]
            lifted = _liftover_interval(lo, hg38[0], hg38[1], hg38[2])
            if lifted:
                mapped = _map_to_circ_pos(lifted[0], lifted[1], lifted[2],
                                          gene_name, exon_nums, coord_idx)
                if mapped:
                    circ_pos = mapped
                    m2 = _circ_pos_re2.match(mapped)
                    if m2:
                        cs, ce = int(m2.group(1)), int(m2.group(2))
                        sites = [{"circ_start": cs, "circ_end": ce, "circ_pos": mapped}]
        in_circ = bool(_circ_pos_re2.match(circ_pos))
        rbp_list.append({
            "RBPName":      name,
            "bindingSites": info["clip"],
            "location":     "internal",
            "sites":        sites,
            "circ_pos":     circ_pos,
            "clipExpNum":   info["clip"],
            "cellType":     ", ".join(sorted(info["cells"])[:3]),
            "source":       "ENCORI",
            "in_circ":      in_circ,
        })
    time.sleep(0.3)
    return rbp_list


# ── Exon structure SVG (linear, kept as fallback) ─────────────────────────────

def _exon_svg(exon_span: str, gene_name: str,
              strand: str, region: str) -> str:
    """Return a linear SVG string of the circRNA exon diagram (fallback)."""
    start_e = end_e = None
    if exon_span:
        m = re.match(r"e(\d+)(?:-e(\d+))?", exon_span.strip())
        if m:
            start_e = int(m.group(1))
            end_e   = int(m.group(2)) if m.group(2) else start_e

    gene_label = f"{gene_name} ({strand})" if gene_name else ""

    if start_e is None:
        label = region or "intergenic"
        return (
            f'<svg width="320" height="80" style="font-family:sans-serif">'
            f'<line x1="10" y1="40" x2="310" y2="40" stroke="#ccc" stroke-width="2"/>'
            f'<rect x="60" y="26" width="200" height="28" fill="#f0ad4e" rx="4" stroke="#c87f00"/>'
            f'<text x="160" y="44" text-anchor="middle" font-size="11" fill="white">{label}</text>'
            f'<text x="160" y="16" text-anchor="middle" font-size="10" fill="#555">{gene_label}</text>'
            f'</svg>'
        )

    n_before = min(2, start_e - 1)
    first_e  = start_e - n_before
    last_e   = end_e + 2
    exons    = list(range(first_e, last_e + 1))

    EW, GAP, EH, PAD = 36, 12, 24, 24
    total_w = PAD * 2 + len(exons) * (EW + GAP) - GAP
    SVG_H   = 110
    y_ex    = 70

    defs = (
        '<defs><marker id="bsj-arr" markerWidth="8" markerHeight="8" '
        'refX="0" refY="3" orient="auto">'
        '<path d="M0,0 L0,6 L6,3 z" fill="#d62728"/></marker></defs>'
    )
    parts = [
        f'<svg width="{total_w}" height="{SVG_H}" viewBox="0 0 {total_w} {SVG_H}" '
        f'style="font-family:sans-serif;overflow:visible">',
        defs,
        f'<line x1="{PAD}" y1="{y_ex + EH//2}" x2="{total_w - PAD}" '
        f'y2="{y_ex + EH//2}" stroke="#ccc" stroke-width="2"/>',
    ]
    ay = y_ex + EH // 2
    if strand == "-":
        ax = PAD - 6
        parts.append(
            f'<polygon points="{ax},{ay} {ax+9},{ay-5} {ax+9},{ay+5}" fill="#aaa"/>'
        )
    else:
        ax = total_w - PAD + 6
        parts.append(
            f'<polygon points="{ax},{ay} {ax-9},{ay-5} {ax-9},{ay+5}" fill="#aaa"/>'
        )
    circ_xs = []
    for i, en in enumerate(exons):
        x    = PAD + i * (EW + GAP)
        circ = start_e <= en <= end_e
        fill, tf, stroke = (
            ("#d62728", "white", "#a01010") if circ else ("#c8dff0", "#445", "#8ab0c8")
        )
        parts += [
            f'<rect x="{x}" y="{y_ex}" width="{EW}" height="{EH}" '
            f'fill="{fill}" rx="4" stroke="{stroke}" stroke-width="1"/>',
            f'<text x="{x + EW//2}" y="{y_ex + EH//2 + 4}" text-anchor="middle" '
            f'font-size="9" fill="{tf}">e{en}</text>',
        ]
        if circ:
            circ_xs.append(x + EW // 2)
    if len(circ_xs) >= 2:
        x1, x2 = circ_xs[0], circ_xs[-1]
        mid    = (x1 + x2) / 2
        arc_h  = max(28, int((x2 - x1) * 0.45))
        cy     = y_ex - arc_h
        parts += [
            f'<path d="M {x2},{y_ex} C {x2},{cy} {x1},{cy} {x1},{y_ex}" '
            f'fill="none" stroke="#d62728" stroke-width="1.8" '
            f'stroke-dasharray="5,3" marker-end="url(#bsj-arr)"/>',
            f'<text x="{mid}" y="{cy - 5}" text-anchor="middle" '
            f'font-size="9" fill="#d62728">back-splice</text>',
        ]
    elif len(circ_xs) == 1:
        parts.append(
            f'<text x="{circ_xs[0]}" y="{y_ex - 8}" text-anchor="middle" '
            f'font-size="9" fill="#d62728">single-exon circRNA</text>'
        )
    span_label = f"exon {exon_span}" if exon_span else region
    parts.append(
        f'<text x="{total_w//2}" y="14" text-anchor="middle" '
        f'font-size="11" font-weight="bold" fill="#333">'
        f'{gene_label}  |  {span_label}</text>'
    )
    parts.append("</svg>")
    return "".join(parts)


# ── GTF exon helpers ─────────────────────────────────────────────────────────

def _build_exon_index(gtf_file):
    # type: (str) -> Tuple[Dict, Dict]
    """
    Stream GTF once → (length_index, coord_index).
    length_index : {gene_name: {exon_num: max_length}}
    coord_index  : {gene_name: {exon_num: (chrom, start0, end0, strand)}}
                   0-based half-open coords; longest exon per number.
    """
    len_idx   = {}  # type: Dict[str, Dict[int, int]]
    coord_idx = {}  # type: Dict[str, Dict[int, tuple]]
    print("[predict] Indexing exon lengths/coords from GTF …", file=sys.stderr, flush=True)
    try:
        with open(gtf_file) as fh:
            for line in fh:
                if line.startswith('#') or '\texon\t' not in line:
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9:
                    continue
                chrom  = parts[0]
                start0 = int(parts[3]) - 1   # GTF 1-based → 0-based
                end0   = int(parts[4])        # GTF end is inclusive → exclusive
                strand = parts[6]
                attr   = parts[8]
                gm = re.search(r'gene_name "([^"]+)"', attr)
                em = re.search(r'exon_number "?(\d+)"?', attr)
                if not gm or not em:
                    continue
                gene   = gm.group(1)
                en     = int(em.group(1))
                length = end0 - start0
                if gene not in len_idx:
                    len_idx[gene]   = {}
                    coord_idx[gene] = {}
                if en not in len_idx[gene] or length > len_idx[gene][en]:
                    len_idx[gene][en]   = length
                    coord_idx[gene][en] = (chrom, start0, end0, strand)
    except Exception as exc:
        print(f"[predict] GTF index error: {exc}", file=sys.stderr)
    print(f"[predict] Indexed {len(len_idx)} genes from GTF", file=sys.stderr, flush=True)
    return len_idx, coord_idx


def _init_liftover():
    """Initialize pyliftover LiftOver(hg38→hg19). Returns None on failure."""
    try:
        from pyliftover import LiftOver  # type: ignore
        lo = LiftOver("hg38", "hg19")
        print("[predict] pyliftover hg38→hg19 ready", file=sys.stderr, flush=True)
        return lo
    except Exception as exc:
        print(f"[predict] LiftOver unavailable: {exc}", file=sys.stderr)
        return None


def _liftover_interval(lo, chrom: str, start0: int, end0: int) -> Optional[Tuple[str, int, int]]:
    """Liftover a 0-based half-open interval from hg38→hg19. Returns None on failure."""
    if lo is None:
        return None
    try:
        r1 = lo.convert_coordinate(chrom, start0)
        r2 = lo.convert_coordinate(chrom, end0 - 1)
        if not r1 or not r2:
            return None
        c1, p1 = r1[0][0], int(r1[0][1])
        c2, p2 = r2[0][0], int(r2[0][1])
        if c1 != c2:
            return None
        return (c1, min(p1, p2), max(p1, p2) + 1)
    except Exception:
        return None


def _get_exon_nums_ordered(exon_span: str, strand: str) -> List[int]:
    """Return exon numbers in mRNA (5′→3′) order for the circRNA."""
    m = re.match(r'e(\d+)-e(\d+)', exon_span.strip())
    if not m:
        return []
    e1, e2 = int(m.group(1)), int(m.group(2))
    nums = list(range(min(e1, e2), max(e1, e2) + 1))
    return nums[::-1] if strand == '-' else nums


def _map_to_circ_pos(hg19_chrom: str, hg19_start: int, hg19_end: int,
                     gene_name: str, exon_nums: List[int],
                     coord_idx: dict) -> Optional[str]:
    """
    Map an hg19 0-based half-open interval to a 1-based circRNA-relative position.
    Returns '1234–1256' string, or None if no overlap found.
    """
    gene_exons = coord_idx.get(gene_name, {})
    cum = 0
    for en in exon_nums:
        coords = gene_exons.get(en)
        if not coords:
            continue
        ec_chrom, ec_start, ec_end, strand = coords
        ex_len = ec_end - ec_start
        if ec_chrom == hg19_chrom:
            ov_s = max(hg19_start, ec_start)
            ov_e = min(hg19_end,   ec_end)
            if ov_s < ov_e:
                if strand == '-':
                    cs = cum + (ec_end - ov_e)   + 1
                    ce = cum + (ec_end - ov_s)
                else:
                    cs = cum + (ov_s - ec_start) + 1
                    ce = cum + (ov_e - ec_start)
                return f"{cs}–{ce}"
        cum += ex_len
    return None


def _get_exon_boundaries(gene_name, exon_span, strand, exon_index, spliced_length):
    # type: (str, str, str, Dict, int) -> List[dict]
    """Return [{label, cum_start, cum_end}] proportional to actual exon lengths."""
    if not gene_name or not exon_span or not exon_index:
        return []
    m = re.match(r'e(\d+)-e(\d+)', exon_span.strip())
    if not m:
        return []
    e1, e2 = int(m.group(1)), int(m.group(2))
    exon_nums = list(range(min(e1, e2), max(e1, e2) + 1))
    if strand == '-':
        exon_nums = exon_nums[::-1]   # minus strand: higher exon# is 5' in mRNA

    gene_exons = exon_index.get(gene_name, {})
    if not gene_exons:
        return []

    lengths = [gene_exons.get(en, 100) for en in exon_nums]
    raw_total = sum(lengths)
    if raw_total == 0:
        return []

    scale = spliced_length / raw_total if spliced_length > 0 else 1.0
    result = []
    cum = 0
    for i, en in enumerate(exon_nums):
        scaled = max(1, int(round(lengths[i] * scale)))
        result.append({'label': f'e{en}', 'cum_start': cum, 'cum_end': cum + scaled})
        cum += scaled
    if result and spliced_length > 0:
        result[-1]['cum_end'] = spliced_length   # exact match
    return result


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Pre-fetch CircInteractome miRNA/RBP interactions for top DE circRNAs"
    )
    p.add_argument("--de",            required=True)
    p.add_argument("--iso",           required=True)
    p.add_argument("--out",           required=True)
    p.add_argument("--circbase",      default=None,
                   help="circbase_annotated.tsv for circBase ID lookup")
    p.add_argument("--gtf",           default=None,
                   help="Reference GTF for exon boundary lengths (optional)")
    p.add_argument("--de-edger",  default=None, help="edgeR DE results TSV")
    p.add_argument("--de-deseq2", default=None, help="DESeq2 DE results TSV")
    p.add_argument("--de-limma",  default=None, help="limma DE results TSV")
    p.add_argument("--top-n",         type=int, default=50)
    p.add_argument("--clip-exp-num",  type=int, default=1,
                   help="ENCORI: min CLIP experiments (1=all, 2=moderate, 3=high)")
    p.add_argument("--program-num",   type=int, default=2,
                   help="ENCORI: min prediction programs agreeing")
    args = p.parse_args()

    de  = pd.read_csv(args.de,  sep="\t")
    iso = pd.read_csv(args.iso, sep="\t")

    # Build method-independent metadata from iso + circbase (covers all circRNAs)
    iso_cols = [c for c in ("circ_id", "gene_name", "strand", "region", "exon_span")
                if c in iso.columns]
    meta = iso[iso_cols].copy() if iso_cols else pd.DataFrame(columns=["circ_id"])

    if args.circbase:
        try:
            cb = pd.read_csv(args.circbase, sep="\t",
                             usecols=lambda c: c in ("circ_id", "circbase_id", "in_circbase"))
            meta = meta.merge(cb, on="circ_id", how="left")
        except Exception as exc:
            print(f"[predict] circbase merge warning: {exc}", file=sys.stderr)

    # Also merge iso/circbase into primary de for backward-compat fallback
    if "circ_id" in de.columns and iso_cols:
        de = de.merge(iso[iso_cols], on="circ_id", how="left")

    exon_index = {}
    coord_index = {}
    if args.gtf and os.path.exists(args.gtf):
        exon_index, coord_index = _build_exon_index(args.gtf)

    lo = _init_liftover()

    def _top_ids_from(path: Optional[str], top_n: int) -> set:
        """Return top_n up + top_n down circ_ids from a DE TSV."""
        if not path or not os.path.exists(path):
            return set()
        df = pd.read_csv(path, sep="\t")
        p_col = "pvalue" if "pvalue" in df.columns else "padj"
        if p_col not in df.columns or "circ_id" not in df.columns:
            return set()
        df = df.dropna(subset=[p_col])
        lfc_col = "log2FC" if "log2FC" in df.columns else "log2FoldChange"
        if lfc_col not in df.columns:
            return set(df.sort_values(p_col).head(top_n * 2)["circ_id"])
        up = set(df[df[lfc_col] > 0].sort_values(p_col).head(top_n)["circ_id"])
        dn = set(df[df[lfc_col] < 0].sort_values(p_col).head(top_n)["circ_id"])
        return up | dn

    # Union mode: collect top-N from each provided method file
    extra_paths = [args.de_edger, args.de_deseq2, args.de_limma]
    if any(extra_paths):
        combined_ids: set = set()
        n_methods = 0
        for path in extra_paths:
            ids = _top_ids_from(path, args.top_n)
            if ids:
                combined_ids |= ids
                n_methods += 1
        print(f"[predict] union from {n_methods} DE methods: {len(combined_ids)} unique circRNAs "
              f"(top-{args.top_n} up + top-{args.top_n} down each)", file=sys.stderr)
    else:
        # Fallback: use primary --de only
        p_col = "pvalue" if "pvalue" in de.columns else "padj"
        de_sig = de.dropna(subset=[p_col])
        top_up = set(de_sig[de_sig["log2FC"] > 0].sort_values(p_col)
                     .head(args.top_n)["circ_id"]) if "log2FC" in de_sig.columns else set()
        top_dn = set(de_sig[de_sig["log2FC"] < 0].sort_values(p_col)
                     .head(args.top_n)["circ_id"]) if "log2FC" in de_sig.columns else set()
        combined_ids = top_up | top_dn
        print(f"[predict] querying {len(combined_ids)} circRNAs "
              f"(top-{args.top_n} up + top-{args.top_n} down by p-value, primary DE only)",
              file=sys.stderr)

    # Build query dataframe from method-independent metadata
    if not meta.empty and "circ_id" in meta.columns:
        top = meta[meta["circ_id"].isin(combined_ids)].drop_duplicates("circ_id").copy()
    else:
        top = de[de["circ_id"].isin(combined_ids)].copy()
    print(f"[predict] metadata resolved for {len(top)}/{len(combined_ids)} circRNAs",
          file=sys.stderr)

    results = {}
    total = len(top)

    for i, (_, row) in enumerate(top.iterrows(), 1):
        circ_id     = str(row["circ_id"])
        circbase_id = str(row.get("circbase_id", "") or "")
        gene_name   = str(row.get("gene_name",  "") or "")
        strand      = str(row.get("strand",     "+") or "+")
        region      = str(row.get("region",     "") or "")
        exon_span   = str(row.get("exon_span",  "") or "")

        can_query = (
            circbase_id
            and circbase_id not in ("", "nan", "novel")
            and circbase_id.startswith("hsa_circ_")
        )

        print(
            f"[predict] [{i}/{total}] {circ_id}  cb={circbase_id or 'novel'}",
            end="  ", file=sys.stderr, flush=True,
        )

        mirna: List[dict] = []
        rbp:   List[dict] = []
        spliced_length = 0

        # ── Source 1: CircInteractome (requires circBase ID) ──
        if can_query:
            mirna, rbp, spliced_length = _query_circinteractome(circbase_id)

        # ── Source 2: ENCORI (requires gene_name, hg38) ──
        can_query_encori = bool(gene_name and gene_name not in ("", "nan"))
        exon_nums = _get_exon_nums_ordered(exon_span, strand) if exon_span else []
        if can_query_encori:
            enc_mirna = _fetch_encori_mirna(gene_name, args.clip_exp_num, args.program_num,
                                            lo=lo, coord_idx=coord_index, exon_nums=exon_nums)
            enc_rbp   = _fetch_encori_rbp(gene_name, args.clip_exp_num,
                                          lo=lo, coord_idx=coord_index, exon_nums=exon_nums)
            mirna = mirna + enc_mirna
            rbp   = rbp   + enc_rbp

        exon_bds = _get_exon_boundaries(gene_name, exon_span, strand, exon_index, spliced_length)  # uses len_index

        n_ci  = sum(1 for x in mirna if x.get("source") == "CircInteractome")
        n_enc = sum(1 for x in mirna if x.get("source") == "ENCORI")
        n_rbp_pos = sum(len(r.get("sites", [])) > 0 for r in rbp)
        print(f"miRNA={len(mirna)}(CI:{n_ci}/ENCORI:{n_enc}) RBP={len(rbp)}(pos:{n_rbp_pos})",
              file=sys.stderr)

        results[circ_id] = {
            "info": {
                "gene_name":       gene_name,
                "strand":          strand,
                "region":          region,
                "exon_span":       exon_span,
                "circbase_id":     circbase_id,
                "spliced_length":  spliced_length,
                "exon_boundaries": exon_bds,
            },
            "exon_svg": _exon_svg(exon_span, gene_name, strand, region),
            "mirna":    mirna,
            "rbp":      rbp,
        }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=2)

    n_mi_ci  = sum(any(x.get("source")=="CircInteractome" for x in v["mirna"]) for v in results.values())
    n_mi_enc = sum(any(x.get("source")=="ENCORI"          for x in v["mirna"]) for v in results.values())
    n_rb_ci  = sum(any(x.get("source")=="CircInteractome" for x in v["rbp"])   for v in results.values())
    n_rb_enc = sum(any(x.get("source")=="ENCORI"          for x in v["rbp"])   for v in results.values())
    print(
        f"[predict] {len(results)} circRNAs → {args.out}\n"
        f"  miRNA: CircInteractome={n_mi_ci}  ENCORI={n_mi_enc}\n"
        f"  RBP:   CircInteractome={n_rb_ci}  ENCORI={n_rb_enc}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
