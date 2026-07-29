"""
generate_report.py – Build a self-contained HTML summary report.

Called as a Snakemake script.  Embeds PDF plots as base64-encoded images
and includes top DE results as a table.
"""

import base64
import html as _html_mod
import math
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

try:
    import numpy as np
    import plotly.graph_objects as go  # type: ignore[import]
    _PLOTLY = True
except ImportError:
    _PLOTLY = False


# ── Helpers ──────────────────────────────────────────────────────────────────

def _pdf_to_b64(pdf_path: str) -> str:
    with open(pdf_path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def _embed_pdf(pdf_path: str, height: int = 480) -> str:
    b64 = _pdf_to_b64(pdf_path)
    return (
        f'<embed src="data:application/pdf;base64,{b64}" '
        f'type="application/pdf" width="100%" height="{height}px" />'
    )


def _fmt_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Format float columns to 3 significant figures to improve readability."""
    out = df.copy()
    for col in out.select_dtypes(include="float").columns:
        out[col] = out[col].apply(
            lambda v: f"{v:.3g}" if pd.notna(v) else "—"
        )
    return out


def _dl_wrap(table_html: str, table_id: str, csv_filename: str) -> str:
    """Inject an id into the first <table> tag and prepend a CSV download button."""
    html = table_html.replace('<table', f'<table id="{table_id}"', 1)
    btn  = (f'<div class="tbl-dl-bar no-print">'
            f'<button class="dl-btn" onclick="dlCSV(\'{table_id}\',\'{csv_filename}\')">'
            f'⬇ CSV</button></div>')
    return btn + html


def _df_to_html(df: pd.DataFrame, max_rows: int = 50, escape: bool = True) -> str:
    return (
        _fmt_floats(df.head(max_rows))
          .to_html(index=False, classes="table", border=0, na_rep="—", escape=escape)
    )


def _make_label(row: pd.Series) -> str:
    """Build a human-readable primary ID: circbase_id or gene+exon, fallback to circ_id."""
    gene  = str(row.get("gene_name", "") or "").strip()
    exon  = str(row.get("exon_span",  "") or "").strip()
    cb_id = str(row.get("circbase_id", "") or "").strip()
    if cb_id and cb_id.lower() not in ("", "nan", "novel"):
        label = cb_id
        if gene:
            label += f" ({gene}"
            if exon:
                label += f" {exon}"
            label += ")"
    elif gene:
        label = gene
        if exon:
            label += f" {exon}"
    else:
        label = str(row.get("circ_id", ""))
    return label


def _eff_sig(de: pd.DataFrame, de_sig_by: str, fdr: float):
    """Return (column, threshold, label) for the effective significance criterion."""
    if de_sig_by == "auto":
        if "qvalue" in de.columns and (de["qvalue"] < 0.2).any():
            return "qvalue", 0.2, "Storey q-value"
        return "pvalue", 0.05, "p-value (nominal)"
    if de_sig_by == "qvalue":
        return "qvalue", 0.2, "Storey q-value"
    if de_sig_by == "pvalue":
        return "pvalue", fdr, "p-value (nominal)"
    return "padj", fdr, "adjusted p-value"


def _sig_ids_from_de(de: pd.DataFrame, p_col: str, sig_thr: float, lfc: float) -> set:
    """Return set of significant circRNA IDs."""
    if p_col not in de.columns or "log2FC" not in de.columns or "circ_id" not in de.columns:
        return set()
    mask = (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)
    return set(de.loc[mask, "circ_id"].dropna().astype(str))


def _venn_3_svg(sig_sets: dict, lfc: float, de_lookup: dict = None) -> str:
    """Generate a 3-circle Venn SVG with clickable regions showing circRNA lists."""
    import json as _json2
    methods = [k for k in ("edgeR_ciriquant", "deseq2", "limma") if k in sig_sets]
    if len(methods) < 2:
        return ""
    sets = [sig_sets[m] for m in methods]
    A, B = sets[0], sets[1]
    C = sets[2] if len(sets) > 2 else set()
    ab = A & B; ac = A & C; bc = B & C; abc = A & B & C
    regions_sets = {
        "A":   A - B - C,  "B":  B - A - C,  "C":  C - A - B,
        "AB":  ab - C,     "AC": ac - B,      "BC": bc - A,
        "ABC": abc,
    }
    counts = {k: len(v) for k, v in regions_sets.items()}
    total_union = len(A | B | C) if len(sets) > 2 else len(A | B)

    # Build per-region circRNA detail rows from de_lookup
    method_labels = {"edgeR_ciriquant": "edgeR", "deseq2": "DESeq2", "limma": "limma"}
    region_labels_zh = {
        "A":   "edgeR only",        "B":   "DESeq2 only",
        "C":   "limma only",        "AB":  "edgeR ∩ DESeq2 only",
        "AC":  "edgeR ∩ limma only","BC":  "DESeq2 ∩ limma only",
        "ABC": "三方法均顯著",
    }
    region_labels_en = {
        "A":   "edgeR only",        "B":   "DESeq2 only",
        "C":   "limma only",        "AB":  "edgeR ∩ DESeq2 only",
        "AC":  "edgeR ∩ limma only","BC":  "DESeq2 ∩ limma only",
        "ABC": "All 3 methods",
    }
    region_labels = region_labels_zh
    venn_region_data: dict = {}
    if de_lookup:
        # Build circ→{gene, lfc, pval, cb, methods[]} — collect ALL significant methods
        circ_info: dict = {}
        for m, df in de_lookup.items():
            if df is None or df.empty:
                continue
            pcol = "pvalue" if "pvalue" in df.columns else ("PValue" if "PValue" in df.columns else None)
            if pcol is None:
                continue
            mlabel = method_labels.get(m, m)
            m_sig_set = sig_sets.get(m, set())
            for _, row in df.iterrows():
                cid = str(row.get("circ_id", ""))
                if not cid:
                    continue
                pv = float(row.get(pcol, 1.0))
                is_sig = cid in m_sig_set
                if cid not in circ_info:
                    circ_info[cid] = {
                        "gene": str(row.get("gene_name", "")),
                        "lfc":  round(float(row.get("log2FC", 0.0)), 3),
                        "pval": pv,
                        "cb":   str(row.get("circbase_id", "novel") or "novel"),
                        "methods": [mlabel] if is_sig else [],
                    }
                else:
                    # update best pval / lfc
                    if pv < circ_info[cid]["pval"]:
                        circ_info[cid]["pval"] = pv
                        circ_info[cid]["lfc"]  = round(float(row.get("log2FC", 0.0)), 3)
                    if is_sig and mlabel not in circ_info[cid]["methods"]:
                        circ_info[cid]["methods"].append(mlabel)
        for region, cid_set in regions_sets.items():
            rows = []
            for cid in cid_set:
                info = circ_info.get(cid, {})
                rows.append({
                    "id":   cid,
                    "gene": info.get("gene", ""),
                    "lfc":  info.get("lfc", 0.0),
                    "pval": info.get("pval", 1.0),
                    "cb":   info.get("cb", "novel"),
                    "m":    ", ".join(info.get("methods", [])),
                })
            rows.sort(key=lambda r: r["pval"])
            venn_region_data[region] = {
                "label_zh": f"{region_labels_zh[region]} ({counts[region]} 個)",
                "label_en": f"{region_labels_en[region]} ({counts[region]})",
                "circs": rows,
            }

    W, H = 460, 345
    cx = [170, 290, 230]; cy = [128, 128, 210]; r = 90
    colors  = ["rgba(214,39,40,0.20)", "rgba(44,160,44,0.20)", "rgba(44,119,214,0.20)"]
    strokes = ["#d62728", "#2CA02C", "#2c6fad"]
    labels  = {"edgeR_ciriquant": "edgeR (FSJ offset)", "deseq2": "DESeq2", "limma": "limma-voom"}
    lpos    = [(170, 22), (290, 22), (230, 316)]
    rpos    = {
        "A":  (118, 118), "B":  (342, 118), "C":  (230, 263),
        "AB": (230, 88),  "AC": (170, 190), "BC": (280, 182),
        "ABC":(230, 163),
    }
    svg_parts = [f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
                 f'id="venn-svg" xmlns="http://www.w3.org/2000/svg" '
                 f'style="font-family:sans-serif;max-width:480px">']
    for i in range(min(len(methods), 3)):
        svg_parts.append(f'<circle cx="{cx[i]}" cy="{cy[i]}" r="{r}" '
                         f'fill="{colors[i]}" stroke="{strokes[i]}" stroke-width="2.5"/>')
    for i, m in enumerate(methods[:3]):
        lx, ly = lpos[i]
        svg_parts.append(f'<text x="{lx}" y="{ly}" font-size="11" font-weight="bold" '
                         f'fill="{strokes[i]}" text-anchor="middle">{labels.get(m,m)}</text>')
        svg_parts.append(f'<text x="{lx}" y="{ly+14}" font-size="10" '
                         f'fill="{strokes[i]}" text-anchor="middle">(n={len(sig_sets[m])})</text>')
    for region, (rx, ry) in rpos.items():
        if len(methods) < 3 and region in ("C", "AC", "BC", "ABC"):
            continue
        v = counts[region]
        if v == 0:
            continue
        has_data = bool(de_lookup) and region in venn_region_data
        cursor    = "cursor:pointer" if has_data else ""
        hover_tip = f'title="點擊查看 {v} 個 circRNA / Click to view {v} circRNAs"' if has_data else ""
        onclick   = f'onclick="showVennRegion(\'{region}\')"' if has_data else ""
        # larger transparent click target
        svg_parts.append(f'<rect x="{rx-16}" y="{ry-14}" width="32" height="28" '
                         f'fill="transparent" rx="5" style="{cursor}" {onclick} {hover_tip}/>')
        svg_parts.append(f'<text x="{rx}" y="{ry}" font-size="15" font-weight="bold" '
                         f'text-anchor="middle" dominant-baseline="middle" '
                         f'stroke="white" stroke-width="4" paint-order="stroke" fill="#333" '
                         f'style="{cursor};text-decoration:{"underline" if has_data else "none"}" '
                         f'{onclick} {hover_tip}>{v}</text>')
    svg_parts.append('</svg>')
    note = (f'<p style="font-size:12px;color:#666;margin:4px 0 0">'
            f'<span data-en="Numbers = significant circRNAs per region (|log₂FC| &gt; {lfc}); union {total_union}.'
            f'{"&lt;b&gt;Click numbers&lt;/b&gt; to view circRNA lists." if de_lookup else ""}">'
            f'數字 = 各區域顯著 circRNA 數量（閾值：|log₂FC| &gt; {lfc}）；'
            f'三方法聯集共 {total_union} 個。'
            f'{"<b>點擊數字</b>可查看對應 circRNA 清單。" if de_lookup else ""}'
            f'</span></p>')

    venn_js = _json2.dumps(venn_region_data, ensure_ascii=False)
    detail_panel = f"""
<div id="venn-detail" style="display:none;margin:12px auto;max-width:860px;
     border:1px solid #dde;border-radius:8px;padding:14px;background:#fafbff">
  <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px">
    <span id="venn-detail-title" style="font-weight:bold;font-size:14px;color:#2c6fad"></span>
    <span>
      <button id="venn-dl-btn" onclick="vennDownloadCSV()"
        style="font-size:12px;padding:3px 10px;border:1px solid #2c6fad;border-radius:4px;
               background:#fff;color:#2c6fad;cursor:pointer;margin-right:6px">&#11015; CSV</button>
      <button onclick="document.getElementById('venn-detail').style.display='none';
                       document.getElementById('venn-detail').dataset.active=''"
        style="font-size:12px;padding:3px 8px;border:1px solid #ccc;border-radius:4px;
               background:#fff;cursor:pointer" data-en="&#x2715; Close">&#x2715; 關閉</button>
    </span>
  </div>
  <div id="venn-detail-body" style="max-height:420px;overflow-y:auto">
    <table id="venn-detail-tbl" class="table" style="font-size:12px;width:100%">
      <thead><tr>
        <th>#</th><th>circ_id</th><th>gene_name</th>
        <th>log2FC</th><th>p-value</th><th>circbase_id</th><th data-en="Methods">方法</th>
      </tr></thead>
      <tbody id="venn-detail-tbody"></tbody>
    </table>
  </div>
</div>
<script>
const VENN_REGION_DATA = {venn_js};
let _vennCSVRows = [];
function showVennRegion(region) {{
  const panel = document.getElementById('venn-detail');
  if (panel.dataset.active === region) {{
    panel.style.display = 'none'; panel.dataset.active = ''; return;
  }}
  panel.dataset.active = region;
  const d = VENN_REGION_DATA[region];
  if (!d) return;
  document.getElementById('venn-detail-title').textContent = d['label_'+(_LANG||'zh')] || d.label_zh;
  const tbody = document.getElementById('venn-detail-tbody');
  tbody.innerHTML = '';
  _vennCSVRows = [['#','circ_id','gene_name','log2FC','p-value','circbase_id','method']];
  d.circs.forEach((c, i) => {{
    const tr = document.createElement('tr');
    const cb = (c.cb && c.cb !== 'novel') ?
      `<span style="color:#e07b39;font-weight:bold">${{_cbLink(c.cb)}}</span>` : 'novel';
    tr.innerHTML = `<td>${{i+1}}</td>
      <td><a class="circ-link" onclick="showCircDetail('${{c.id}}')" title="View detail">${{c.id}}</a></td>
      <td>${{c.gene||'—'}}</td>
      <td style="color:${{c.lfc>0?'#d62728':'#2CA02C'}};font-weight:bold">${{c.lfc.toFixed(2)}}</td>
      <td>${{c.pval < 0.001 ? c.pval.toExponential(2) : c.pval.toFixed(4)}}</td>
      <td>${{cb}}</td><td style="color:#888">${{c.m}}</td>`;
    tbody.appendChild(tr);
    _vennCSVRows.push([i+1, c.id, c.gene||'', c.lfc, c.pval, c.cb||'novel', c.m]);
  }});
  panel.style.display = 'block';
  panel.scrollIntoView({{behavior:'smooth', block:'nearest'}});
}}
function vennDownloadCSV() {{
  const panel = document.getElementById('venn-detail');
  const region = panel.dataset.active || 'venn';
  const d = VENN_REGION_DATA[region];
  const filename = 'GSE133998_venn_' + region + '.csv';
  const csv = _vennCSVRows.map(r => r.map(v => '"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\\n');
  const a = document.createElement('a');
  a.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
  a.download = filename; a.click();
}}
</script>"""

    return (f'<div style="text-align:center">{"".join(svg_parts)}{note}</div>'
            + detail_panel)


def _build_method_js_data(
    de: pd.DataFrame,
    matrix: pd.DataFrame,
    p_col: str,
    sig_thr: float,
    lfc: float,
    groups_file: Optional[str],
    normal_label: str,
    hm_pool: int = 50,
) -> dict:
    """Build {stats, volcano, heatmap} for one DE method (for ALL_DE_METHODS JS)."""
    import json as _json
    # Stats
    n_sig = n_up = n_dn = 0
    if p_col in de.columns and "log2FC" in de.columns:
        mask = (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)
        sig  = de[mask]
        n_sig = len(sig)
        n_up  = int((sig["log2FC"] > 0).sum()) if n_sig else 0
        n_dn  = int((sig["log2FC"] < 0).sum()) if n_sig else 0
    # Volcano rows
    vol_rows = []
    if "circ_id" in de.columns and p_col in de.columns and "log2FC" in de.columns:
        for _, r in de.iterrows():
            if pd.isna(r["log2FC"]) or pd.isna(r.get(p_col)):
                continue
            x  = round(float(r["log2FC"]), 3)
            pv = float(r[p_col])
            y  = round(-math.log10(max(pv, 1e-300)), 3)
            s  = ("U" if pv < sig_thr and r["log2FC"] > lfc else
                  "D" if pv < sig_thr and r["log2FC"] < -lfc else "N")
            vol_rows.append([x, y, s, str(r["circ_id"])])
    # Heatmap data (mirrors FULL_HEATMAP_DATA logic)
    hm_data = None
    if not de.empty and not matrix.empty:
        try:
            if {"circ_id", "log2FC", p_col}.issubset(de.columns):
                _sig_hm = de[(de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)] if {p_col, "log2FC"}.issubset(de.columns) else de
                up_p = _sig_hm[_sig_hm["log2FC"] > 0].dropna(subset=[p_col]).sort_values(p_col).head(hm_pool)
                dn_p = _sig_hm[_sig_hm["log2FC"] < 0].dropna(subset=[p_col]).sort_values(p_col).head(hm_pool)
                hm_ids  = list(up_p["circ_id"]) + list(dn_p["circ_id"])
                hm_avail = [i for i in hm_ids if i in matrix.index]
                if hm_avail:
                    samps  = matrix.columns.tolist()
                    sub    = matrix.loc[hm_avail].astype(float)
                    log_s  = (sub + 1).apply(lambda col: col.apply(lambda v: math.log2(v) if v > 0 else 0.0))
                    cmap   = {}
                    if groups_file and Path(groups_file).exists():
                        try:
                            grp_df = pd.read_csv(groups_file)
                            cmap   = dict(zip(grp_df["srr_id"].astype(str), grp_df["condition"].astype(str)))
                        except Exception:
                            pass
                    ncols  = [c for c in samps if cmap.get(c, "") == normal_label]
                    rcols  = ncols if ncols else samps
                    rmean  = log_s[rcols].mean(axis=1)
                    rstd   = log_s.std(axis=1).clip(lower=0.1)
                    z_hm   = log_s.sub(rmean, axis=0).div(rstd, axis=0)
                    # Reorder columns: tumor first, then normal
                    _t_c = [c for c in samps if cmap.get(c, "") not in (normal_label, "")]
                    _n_c = [c for c in samps if cmap.get(c, "") == normal_label]
                    _o_c = [c for c in samps if c not in set(_t_c) and c not in set(_n_c)]
                    if _t_c or _n_c:
                        samps = _t_c + _n_c + _o_c
                        z_hm  = z_hm[samps]
                    cb_map = {}
                    if "circbase_id" in de.columns and "in_circbase" in de.columns:
                        kn = de[de["in_circbase"] == 1].dropna(subset=["circbase_id"])
                        cb_map = dict(zip(kn["circ_id"].astype(str), kn["circbase_id"].astype(str)))
                    pv_map = dict(zip(de["circ_id"].astype(str), de[p_col]))
                    lf_map = dict(zip(de["circ_id"].astype(str), de["log2FC"]))
                    rows_hm = {}
                    for cid in hm_avail:
                        cb  = cb_map.get(cid, "")
                        lbl = cb if cb and cb not in ("", "novel") else cid
                        rows_hm[cid] = {
                            "z":      [round(v, 3) for v in z_hm.loc[cid].tolist()],
                            "pval":   float(pv_map.get(cid, 1.0)),
                            "log2fc": float(lf_map.get(cid, 0.0)),
                            "label":  lbl,
                        }
                    hm_data = {
                        "samples":    samps,
                        "conditions": cmap,
                        "up_order":   [i for i in up_p["circ_id"] if i in set(hm_avail)],
                        "dn_order":   [i for i in dn_p["circ_id"] if i in set(hm_avail)],
                        "rows":       rows_hm,
                    }
        except Exception:
            pass
    # DE table data (top 50 up + 50 down) for dynamic re-render on method switch
    de_table = {"cols": [], "up": [], "dn": []}
    sig_ids: list = []
    if p_col in de.columns and "log2FC" in de.columns and "circ_id" in de.columns:
        _mask = (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)
        sig_ids = [str(x) for x in de.loc[_mask, "circ_id"].tolist()]
        _sig_de = de[_mask].sort_values(p_col)
        _want = ["circ_id", "gene_name", "strand", "region", "exon_span",
                 "circbase_id", "log2FC", p_col, "Type"]
        _avail = [c for c in _want if c in de.columns]
        de_table["cols"] = [("p-value" if c == p_col else c) for c in _avail]
        def _fmt_row(r):
            row = []
            for c in _avail:
                v = r[c]
                if isinstance(v, float) and math.isnan(v):
                    row.append(None)
                elif c == "log2FC":
                    row.append(round(float(v), 3))
                elif c == p_col:
                    row.append(float(v))
                else:
                    row.append(str(v))
            return row
        de_table["up"] = [_fmt_row(r) for _, r in _sig_de[_sig_de["log2FC"] > 0].head(50).iterrows()]
        de_table["dn"] = [_fmt_row(r) for _, r in _sig_de[_sig_de["log2FC"] < 0].head(50).iterrows()]

    return {"stats": {"n_sig": n_sig, "n_up": n_up, "n_dn": n_dn},
            "volcano": vol_rows, "heatmap": hm_data,
            "de_table": de_table, "sig_ids": sig_ids}


def _load_interactions(json_file: Optional[str]) -> dict:
    """Load pre-fetched ENCORI interactions JSON, return {} on any error."""
    if not json_file:
        return {}
    try:
        import json as _json
        return _json.loads(Path(json_file).read_text(encoding="utf-8"))
    except Exception:
        return {}


_CIRCBASE_URL = "https://www.circbase.org/cgi-bin/singlerecord.cgi?id="


def _cb_link(cb_id) -> str:
    """Render a circbase_id value as a link to its circBase record page, or
    a plain '—'/'novel' fallback for missing/non-circBase entries."""
    v = str(cb_id or "").strip()
    if not v or v.lower() in ("nan", "novel", "none", "—"):
        return v if v else "—"
    from urllib.parse import quote as _urlquote
    return (f'<a href="{_CIRCBASE_URL}{_urlquote(v)}" target="_blank" '
            f'rel="noopener" class="cb-link">{v}</a>')


def _de_table_clickable(df: pd.DataFrame, interactions: dict) -> str:
    """Render DE table with circ_position column as clickable links."""
    if df.empty:
        return "<p class='no-data'>No data.</p>"
    disp = df.copy()
    id_col = "circ_position" if "circ_position" in disp.columns else (
             "circ_id" if "circ_id" in disp.columns else None)
    if id_col:
        def _link(v: str) -> str:
            tip = "in interaction data" if str(v) in interactions else "no interaction data pre-fetched"
            return (f'<a class="circ-link" onclick="showCircDetail(\'{v}\')" '
                    f'title="{tip}">{v}</a>')
        disp[id_col] = disp[id_col].astype(str).apply(_link)
    if "circbase_id" in disp.columns:
        disp["circbase_id"] = disp["circbase_id"].apply(_cb_link)
    disp = _fmt_floats(disp)
    return disp.head(50).to_html(index=False, classes="table", border=0,
                                  na_rep="—", escape=False)


def _de_split_tables(sig: pd.DataFrame, tumor_label: str = "tumor",
                     normal_label: str = "normal",
                     interactions: Optional[dict] = None) -> str:
    """Return two HTML tables: up-regulated and down-regulated in tumor."""
    if sig.empty or "log2FC" not in sig.columns:
        return _df_to_html(sig)
    ixn  = interactions or {}
    up   = sig[sig["log2FC"] > 0].copy()
    down = sig[sig["log2FC"] < 0].copy()
    html_parts = []
    if not up.empty:
        html_parts.append(
            f'<h3 style="color:#d62728">&#8593; Up-regulated in {tumor_label} '
            f'(log₂FC &gt; 0) — {len(up)} circRNAs</h3>'
        )
        html_parts.append(_dl_wrap(_de_table_clickable(up, ixn),
                                   "tbl_de_up", "de_upregulated.csv"))
    if not down.empty:
        html_parts.append(
            f'<h3 style="color:#2CA02C">&#8595; Down-regulated in {tumor_label} / '
            f'Up-regulated in {normal_label} (log₂FC &lt; 0) — {len(down)} circRNAs</h3>'
        )
        html_parts.append(_dl_wrap(_de_table_clickable(down, ixn),
                                   "tbl_de_dn", "de_downregulated.csv"))
    return "\n".join(html_parts) if html_parts else _df_to_html(sig)


# ── Report template ──────────────────────────────────────────────────────────

_STYLE = """
<style>
  body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; color: #222; }
  h1   { color: #2c6fad; border-bottom: 2px solid #2c6fad; padding-bottom: 8px; }
  h2   { color: #444; margin-top: 40px; }
  .table { border-collapse: collapse; width: 100%; font-size: 13px; }
  .table th { background: #2c6fad; color: #fff; padding: 6px 10px; text-align: left;
              cursor: pointer; user-select: none; white-space: nowrap; }
  .table th:hover { background: #1a5090; }
  .table td { border-bottom: 1px solid #ddd; padding: 5px 10px; }
  .table tr:hover td { background: #f0f7ff; }
  .stat-box { display:inline-block; background:#f4f8ff; border:1px solid #b6d0f0;
              border-radius:6px; padding:12px 24px; margin:6px; min-width:150px; text-align:center; }
  .stat-box .num { font-size:2em; font-weight:bold; color:#2c6fad; }
  .stat-box .lbl { font-size:0.85em; color:#666; }
  .type-bar { display:flex; align-items:center; gap:12px; margin: 8px 0; }
  .type-bar .bar-wrap { flex:1; background:#e8e8e8; border-radius:4px; height:22px; }
  .type-bar .bar-fill { height:22px; border-radius:4px; display:flex; align-items:center;
                        padding-left:8px; color:#fff; font-size:12px; font-weight:bold; }
  .bar-type1 { background: #2c6fad; }
  .bar-type2 { background: #e07b39; }
  .badge { display:inline-block; padding:2px 9px; border-radius:10px; font-size:0.82em;
           font-weight:bold; color:#fff; }
  .badge-type1 { background:#2c6fad; }
  .badge-type2 { background:#e07b39; }
  .method-tag { display:inline-block; background:#eaf3ff; border:1px solid #99c2f0;
                border-radius:4px; padding:2px 10px; font-size:0.9em; color:#2c6fad; }
  /* Clickable circRNA link */
  .circ-link { color:#2c6fad; text-decoration:none; cursor:pointer; font-family:monospace; font-size:12px; }
  .circ-link:hover { text-decoration:underline; }
  /* Modal overlay */
  .circ-modal { display:none; position:fixed; z-index:2000; inset:0;
                background:rgba(0,0,0,0.55); justify-content:center; align-items:center; }
  .circ-modal-box { background:#fff; border-radius:10px; padding:28px 32px;
                    max-width:760px; width:94%; max-height:85vh; overflow-y:auto;
                    position:relative; box-shadow:0 8px 32px rgba(0,0,0,0.25); }
  .circ-modal-close { position:absolute; top:14px; right:18px; font-size:22px;
                      cursor:pointer; color:#888; line-height:1; }
  .circ-modal-close:hover { color:#333; }
  .circ-modal-title { font-size:1.15em; font-weight:bold; color:#333; margin-bottom:2px; }
  .circ-modal-sub { font-size:12px; color:#888; margin-bottom:14px; }
  /* Tabs */
  .ctab-bar { display:flex; gap:0; border-bottom:2px solid #e0e8f0; margin-bottom:14px; }
  .ctab-btn { padding:7px 20px; border:none; background:none; cursor:pointer;
              color:#666; font-size:13px; border-bottom:3px solid transparent;
              margin-bottom:-2px; }
  .ctab-btn.active { color:#2c6fad; border-bottom-color:#2c6fad; font-weight:bold; }
  .ctab-content { display:none; }
  .ctab-content.active { display:block; }
  /* Exon SVG container */
  .exon-wrap { overflow-x:auto; padding:8px 0; }
  /* Interaction tables */
  .itable { border-collapse:collapse; width:100%; font-size:12px; }
  .itable th { background:#f0f5ff; color:#335; padding:5px 10px; border:1px solid #dde; text-align:left;
              cursor:pointer; user-select:none; white-space:nowrap; }
  .itable th:hover { background:#dce8ff; }
  .itable td { padding:4px 10px; border:1px solid #eee; }
  .itable tr:hover td { background:#f8faff; }
  .no-data { color:#aaa; font-size:13px; font-style:italic; padding:8px 0; }
  .conf-high { color:#2c6fad; font-weight:bold; }

  /* Download / print toolbar */
  .dl-btn { display:inline-block; margin:4px 4px 4px 0; padding:4px 10px;
            font-size:11px; font-weight:bold; border:1px solid #2c6fad;
            border-radius:4px; color:#2c6fad; background:#fff; cursor:pointer;
            text-decoration:none; transition:background .15s; }
  .dl-btn:hover { background:#eaf3ff; }
  .tbl-dl-bar { text-align:right; margin-bottom:2px; }
  /* ── CircDEX Report Header ── */
  .print-bar  { position:sticky; top:0; z-index:99;
                background:#0F2137; color:#fff;
                padding:0 24px; display:flex; align-items:center;
                min-height:70px; gap:0; flex-wrap:nowrap; }
  .print-bar::after { content:''; position:absolute; bottom:0; left:0; right:0; height:1px;
                      background:linear-gradient(90deg,transparent,#00B4C6,transparent);
                      opacity:.4; }
  .cd-rpt-brand { display:flex; align-items:center; gap:12px; text-decoration:none;
                  padding:10px 0; flex-shrink:0; }
  .cd-rpt-meta { display:flex; flex-direction:column; }
  .cd-rpt-wm   { font-size:26px; line-height:1; }
  .cd-rpt-circ { font-weight:300; color:rgba(255,255,255,.82); }
  .cd-rpt-dex  { font-weight:800; letter-spacing:.07em; color:#00B4C6; }
  .cd-rpt-sub  { font-size:14px; color:rgba(255,255,255,.52); margin-top:3px; }
  .cd-rpt-chips { display:flex; gap:5px; margin-top:5px; }
  .cd-rpt-chip { font-size:11px; font-weight:700; letter-spacing:.04em; text-transform:uppercase;
                 color:#00B4C6; background:rgba(0,180,198,.14);
                 border:1px solid rgba(0,180,198,.30); border-radius:4px; padding:3px 9px; }
  .cd-rpt-nav  { margin-left:auto; display:flex; align-items:center; gap:6px; flex-shrink:0; }
  .cd-rpt-proj { font-size:12px; font-weight:600; color:rgba(255,255,255,.65);
                 background:rgba(255,255,255,.08); border:1px solid rgba(255,255,255,.15);
                 border-radius:5px; padding:3px 10px; }
  .print-btn  { padding:5px 13px; border:1.5px solid rgba(255,255,255,.35);
                border-radius:5px; background:transparent; color:rgba(255,255,255,.75);
                font-size:12px; font-weight:600; cursor:pointer; transition:all .15s; }
  .print-btn:hover { background:rgba(255,255,255,.10); color:#fff;
                     border-color:rgba(255,255,255,.6); }

  /* DE method switcher */
  .msw-btn { padding:5px 14px; border:1px solid #2c6fad; border-radius:4px;
             background:white; color:#2c6fad; cursor:pointer; font-size:13px;
             margin-right:6px; transition:all 0.15s; }
  .msw-btn.active { background:#2c6fad; color:white; }
  .msw-btn:hover:not(.active) { background:#eaf3ff; }

  @media print {
    .print-bar, .dl-btn, .no-print { display:none !important; }
    body { max-width:100%; margin:0; padding:0 12px; font-size:11px; }

    /* Keep heading + immediately following content together */
    h1, h2, h3 {
      break-after: avoid;
      page-break-after: avoid;
      color:#2c6fad !important;
      -webkit-print-color-adjust: exact;
      print-color-adjust: exact;
    }
    /* Avoid breaking inside tables, stat boxes, type-bar, biomarker blocks */
    .table { break-inside: avoid; page-break-inside: avoid; font-size:10px; }
    .table th { background:#2c6fad !important; -webkit-print-color-adjust:exact;
                print-color-adjust:exact; }
    .stat-box { break-inside: avoid; page-break-inside: avoid; }
    .type-bar { break-inside: avoid; page-break-inside: avoid; }

    /* Section wrappers: keep heading + first element together */
    #de-tables-section, #biomarker-section, #isoform-section {
      break-inside: auto;
    }
    #de-tables-section h2, #biomarker-section h2, #isoform-section h2 {
      break-after: avoid;
      page-break-after: avoid;
    }

    /* Plotly charts: scale to fit page width and avoid splitting */
    .plotly-graph-div, [id$="-plot"], [id="main-heatmap-plot"],
    [id="main-volcano-plot"] {
      break-inside: avoid;
      page-break-inside: avoid;
      max-width: 100% !important;
      width: 100% !important;
    }
    /* Volcano and Heatmap: try to keep on one page */
    h2 + p + .plotly-graph-div,
    h2 + .plotly-graph-div { break-before: avoid; page-break-before: avoid; }

    /* Venn diagram: keep together */
    #biomarker-section + * div[style*="text-align:center"] {
      break-inside: avoid;
    }

    /* Force page break before major sections for clean layout */
    #de-tables-section { break-before: page; page-break-before: always; }

    /* Two-column distribution plots: stack vertically for print */
    div[style*="display:flex"] { display: block !important; }
  }
  .lang-btn-rpt { padding:2px 9px; border:1px solid rgba(255,255,255,.22); border-radius:4px;
                  font-size:11px; font-weight:700; cursor:pointer;
                  background:transparent; color:rgba(255,255,255,.50);
                  transition:all .15s; }
  .lang-btn-rpt:hover { background:rgba(255,255,255,.10); color:#fff; }
  .lang-btn-rpt.active { background:#00B4C6; border-color:#00B4C6; color:#fff; }
  #qc-section > summary span:first-child { transition:transform .2s; display:inline-block; }
  #qc-section[open] > summary span:first-child { transform:rotate(90deg); }
</style>
"""

_SCRIPT = """
<script>
// ── Language switcher ─────────────────────────────────────────────────────────
let _LANG = localStorage.getItem('circrna_report_lang') || 'zh';
const _LS = {
  zh: {
    hmStatus:   (n,u,d)=>`（已顯示 ${n} 筆；資料庫最多 ${u} up + ${d} down）`,
    swConc:     p=>p<0.05?'✗ 非常態 / Non-normal (α=0.05)':'✓ 常態 / Normal (α=0.05)',
    obsDist:    '觀測分布',
    bmAll:      n=>`全部（${n}）`,
    bm2:        n=>`≥ 2 方法顯著（${n}）`,
    bm3:        n=>`3 方法均顯著（${n}）`,
    bmNSig:     '在所選方法下不顯著',
    toggleShow: '全顯示',
    toggleHide: '全隱藏',
    toggleTip:  '點擊顯示/隱藏',
  },
  en: {
    hmStatus:   (n,u,d)=>`(showing ${n}; pool: ${u} up + ${d} down)`,
    swConc:     p=>p<0.05?'✗ Non-normal (α=0.05)':'✓ Normal (α=0.05)',
    obsDist:    'Observed distribution',
    bmAll:      n=>`All (${n})`,
    bm2:        n=>`≥ 2 methods (${n})`,
    bm3:        n=>`All 3 methods (${n})`,
    bmNSig:     'Not significant under selected method',
    toggleShow: 'Show all',
    toggleHide: 'Hide all',
    toggleTip:  'Click to show/hide',
  },
};

function switchReportLang(lang) {
  _LANG = lang;
  localStorage.setItem('circrna_report_lang', lang);
  document.querySelectorAll('[data-en]').forEach(function(el) {
    if (lang === 'en') {
      if (!el.dataset.zh) el.dataset.zh = el.innerHTML;
      el.innerHTML = el.dataset.en;
    } else {
      if (el.dataset.zh !== undefined) el.innerHTML = el.dataset.zh;
    }
  });
  document.querySelectorAll('.lang-btn-rpt').forEach(function(b) {
    b.classList.toggle('active', b.dataset.lang === lang);
  });
  const statusEl = document.getElementById('heatmap-status');
  if (statusEl && statusEl.dataset.n) {
    statusEl.textContent = _LS[_LANG].hmStatus(statusEl.dataset.n, statusEl.dataset.u, statusEl.dataset.d);
  }
  _updateBMFilterBtnLabels();
  const vennPanel = document.getElementById('venn-detail');
  if (vennPanel && vennPanel.dataset.active) {
    const d = VENN_REGION_DATA[vennPanel.dataset.active];
    if (d) document.getElementById('venn-detail-title').textContent = d['label_' + _LANG] || d.label_zh;
  }
  _refreshBMHistLang();
}

function _applyLangToContainer(el) {
  /* Re-apply current report language to dynamically injected HTML (e.g., modal tabs). */
  if (!el || typeof _LANG === 'undefined') return;
  el.querySelectorAll('[data-en]').forEach(function(node) {
    if (_LANG === 'en') {
      if (!node.dataset.zh) node.dataset.zh = node.innerHTML;
      node.innerHTML = node.dataset.en;
    } else {
      if (node.dataset.zh !== undefined) node.innerHTML = node.dataset.zh;
    }
  });
}

function _updateBMFilterBtnLabels() {
  const btns = document.querySelectorAll('.bm-filter-btn');
  const tbody = document.querySelector('#tbl_biomarker tbody');
  if (!btns.length || !tbody) return;
  const allRows = [...tbody.querySelectorAll('tr')];
  const n_all = allRows.length;
  const n2 = allRows.filter(r => (parseInt(r.getAttribute('data-nsig')||'1')) >= 2).length;
  const n3 = allRows.filter(r => (parseInt(r.getAttribute('data-nsig')||'1')) >= 3).length;
  btns.forEach((b, i) => {
    if (i === 0) b.textContent = _LS[_LANG].bmAll(n_all);
    else if (i === 1) b.textContent = _LS[_LANG].bm2(n2);
    else if (i === 2) b.textContent = _LS[_LANG].bm3(n3);
  });
}

function _refreshBMHistLang() {
  if (typeof Plotly === 'undefined' || !document.getElementById('bm-hist-plot')) return;
  try {
    Plotly.restyle('bm-hist-plot', {name: [_LS[_LANG].obsDist]}, [0]);
  } catch(e) {}
}

(function() {
  var saved = localStorage.getItem('circrna_report_lang');
  if (saved === 'en') switchReportLang('en');
})();

function _cbLink(cb) {
  // Render a circbase_id value as a link to its circBase record page.
  var v = (cb == null ? '' : String(cb)).trim();
  if (!v || v.toLowerCase() === 'nan' || v.toLowerCase() === 'novel' || v === '—') return v || '—';
  return '<a href="https://www.circbase.org/cgi-bin/singlerecord.cgi?id=' + encodeURIComponent(v) +
         '" target="_blank" rel="noopener" class="cb-link">' + v + '</a>';
}
function _parseGenCoord(s) {
  // Parse chrN:start|end → {c: chrNum, p: startPos} for genomic sort
  var m = s.match(/^chr([0-9]+|[XYMxym])[:|\s](\d+)/i);
  if (!m) return null;
  var n = m[1].toUpperCase();
  var cn = /^\d+$/.test(n) ? parseInt(n) : (n==='X'?23:n==='Y'?24:n==='M'?26:25);
  return {c: cn, p: parseInt(m[2])};
}
function _makeSortable(tableId) {
  var tbl = document.getElementById(tableId);
  if (!tbl) return;
  var ths = Array.from(tbl.querySelectorAll('thead th'));
  var _sc = null, _asc = true;
  ths.forEach(function(th, i) {
    th.addEventListener('click', function() {
      if (_sc === i) { _asc = !_asc; } else { _sc = i; _asc = true; }
      ths.forEach(function(h, j) {
        h.innerHTML = h.innerHTML.replace(/\s*[▲▼]$/, '');
        if (j === i) h.innerHTML += _asc ? ' ▲' : ' ▼';
      });
      var tbody = tbl.querySelector('tbody');
      if (!tbody) return;
      var rows = Array.from(tbody.querySelectorAll('tr'));
      rows.sort(function(a, b) {
        var av = (a.cells[i] ? a.cells[i].textContent : '').trim();
        var bv = (b.cells[i] ? b.cells[i].textContent : '').trim();
        // Genomic coordinate sort (chrN:start format)
        var ca = _parseGenCoord(av), cb = _parseGenCoord(bv);
        if (ca && cb) {
          if (ca.c !== cb.c) return _asc ? ca.c - cb.c : cb.c - ca.c;
          return _asc ? ca.p - cb.p : cb.p - ca.p;
        }
        var an = parseFloat(av);
        var bn = parseFloat(bv);
        if (!isNaN(an) && !isNaN(bn)) return _asc ? an - bn : bn - an;
        return _asc ? av.localeCompare(bv) : bv.localeCompare(av);
      });
      rows.forEach(function(r) { tbody.appendChild(r); });
    });
  });
}

function dlCSV(tid, fname) {
  var tbl = document.getElementById(tid);
  if (!tbl) return;
  var rows = tbl.querySelectorAll('tr');
  var lines = [];
  rows.forEach(function(r) {
    var cells = r.querySelectorAll('th,td');
    var cols = [];
    cells.forEach(function(c) {
      cols.push('"' + c.innerText.replace(/"/g,'""') + '"');
    });
    lines.push(cols.join(','));
  });
  var blob = new Blob([lines.join('\\n')], {type:'text/csv;charset=utf-8;'});
  var a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = fname;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
}
</script>
"""


def _sample_overview_section(
    groups_file: Optional[str],
    tumor_label: str,
    normal_label: str,
    results_dir: Optional[str] = None,
    study_title: str = "",
) -> str:
    """Return an HTML section describing samples + fastp QC stats."""
    if not groups_file or not Path(groups_file).exists():
        return ""
    try:
        grp = pd.read_csv(groups_file)
    except Exception:
        return ""
    if "srr_id" not in grp.columns or "condition" not in grp.columns:
        return ""

    cond_counts = grp["condition"].value_counts().to_dict()
    n_case  = cond_counts.get(tumor_label, 0)
    n_ctrl  = cond_counts.get(normal_label, 0)
    n_other = sum(v for k, v in cond_counts.items() if k not in (tumor_label, normal_label))

    # Try to read fastp JSON per sample
    import json as _json
    rows = []
    for _, row in grp.iterrows():
        srr  = str(row["srr_id"])
        cond = str(row.get("condition", ""))
        pid  = str(row["patient_id"]) if "patient_id" in grp.columns else ""
        total_reads = q30 = mean_len = ""
        if results_dir:
            jpath = Path(results_dir) / "qc" / "fastp" / f"{srr}.json"
            if jpath.exists():
                try:
                    with open(jpath) as fh:
                        qc = _json.load(fh)
                    bf = qc.get("summary", {}).get("before_filtering", {})
                    total_reads = f'{bf.get("total_reads", 0) / 1e6:.1f}M'
                    mean_len    = str(int(bf.get("read1_mean_length", 0))) + " bp"
                    q30         = f'{bf.get("q30_rate", 0) * 100:.1f}%'
                except Exception:
                    pass
        rows.append({"srr": srr, "cond": cond, "pid": pid,
                     "reads": total_reads, "len": mean_len, "q30": q30})

    # Sort: case first then control, then by srr
    order = {tumor_label: 0, normal_label: 1}
    rows.sort(key=lambda r: (order.get(r["cond"], 2), r["srr"]))

    has_qc  = any(r["reads"] for r in rows)
    has_pid = any(r["pid"] for r in rows)

    def _cond_badge(c: str) -> str:
        color = "#2c6fad" if c == tumor_label else "#16a34a" if c == normal_label else "#888"
        return (f'<span style="background:{color};color:white;padding:1px 8px;'
                f'border-radius:10px;font-size:11px">{c}</span>')

    tbody = ""
    for i, r in enumerate(rows):
        bg = "background:#f9fafb;" if i % 2 == 0 else ""
        pid_cell = f'<td style="padding:5px 12px">{r["pid"]}</td>' if has_pid else ""
        qc_cells = (
            f'<td style="padding:5px 12px;text-align:right">{r["reads"] or "—"}</td>'
            f'<td style="padding:5px 12px;text-align:right">{r["len"] or "—"}</td>'
            f'<td style="padding:5px 12px;text-align:right">{r["q30"] or "—"}</td>'
        ) if has_qc else ""
        tbody += (
            f'<tr style="{bg}">'
            f'<td style="padding:5px 12px;font-family:monospace;font-size:12px">{r["srr"]}</td>'
            f'<td style="padding:5px 12px">{_cond_badge(r["cond"])}</td>'
            f'{pid_cell}{qc_cells}</tr>\n'
        )

    pid_th  = '<th style="padding:6px 12px;text-align:left;border-bottom:1px solid #e0e8f0">Patient</th>' if has_pid else ""
    qc_ths  = (
        '<th style="padding:6px 12px;text-align:right;border-bottom:1px solid #e0e8f0">Total Reads</th>'
        '<th style="padding:6px 12px;text-align:right;border-bottom:1px solid #e0e8f0">Avg Length</th>'
        '<th style="padding:6px 12px;text-align:right;border-bottom:1px solid #e0e8f0">Q30</th>'
    ) if has_qc else ""

    other_box = (
        f'<div class="stat-box"><div class="num" style="color:#888">{n_other}</div>'
        f'<div class="lbl">Other</div></div>'
    ) if n_other else ""

    title_html = (
        f'<p style="font-size:14px;color:#374151;margin:0 0 12px;font-style:italic">'
        f'{study_title}</p>'
    ) if study_title else ""

    return f"""
<h2>Samples</h2>
{title_html}<div style="display:flex;gap:16px;flex-wrap:wrap;margin:8px 0 16px">
  <div class="stat-box"><div class="num" style="color:#2c6fad">{n_case}</div><div class="lbl">{tumor_label.title()}</div></div>
  <div class="stat-box"><div class="num" style="color:#16a34a">{n_ctrl}</div><div class="lbl">{normal_label.title()}</div></div>
  {other_box}
</div>
<details style="margin:0 0 20px">
  <summary style="cursor:pointer;color:#2c6fad;font-size:13px;user-select:none">
    <span data-en='&#9654; Show sample list ({len(rows)} samples){"&nbsp;&nbsp;<span style=&apos;color:#888;font-size:11px&apos;>with fastp QC stats</span>" if has_qc else ""}'>&#9654; 顯示樣本清單（{len(rows)} samples）{"&nbsp;&nbsp;<span style='color:#888;font-size:11px'>含 fastp QC 統計</span>" if has_qc else ""}</span>
  </summary>
  <div style="overflow-x:auto;margin-top:10px">
  <table style="border-collapse:collapse;font-size:13px;width:100%;max-width:800px">
    <thead>
      <tr style="background:#f1f5f9;font-weight:600">
        <th style="padding:6px 12px;text-align:left;border-bottom:1px solid #e0e8f0">SRR ID</th>
        <th style="padding:6px 12px;text-align:left;border-bottom:1px solid #e0e8f0">Condition</th>
        {pid_th}{qc_ths}
      </tr>
    </thead>
    <tbody>{tbody}</tbody>
  </table>
  </div>
</details>
"""


def _type_section(sig: pd.DataFrame) -> str:
    """Return HTML for Type I/II breakdown; empty string if no Type column."""
    if "Type" not in sig.columns or sig.empty:
        return ""
    n_type1 = int((sig["Type"] == "Type_I").sum())
    n_type2 = int((sig["Type"] == "Type_II").sum())
    total   = n_type1 + n_type2
    if total == 0:
        return ""
    pct1 = round(100 * n_type1 / total)
    pct2 = 100 - pct1
    return f"""
  <h2>Type I / II Classification <small style="font-size:0.7em;color:#888">(edgeR + per-locus FSJ offset)</small></h2>
  <p>
    <span class="badge badge-type1">Type I</span>&nbsp;circRNA-specific splicing change (BSJ/FSJ ratio shift, FSJ stable) &nbsp;|&nbsp;
    <span class="badge badge-type2">Type II</span>&nbsp;Gene-level change (FSJ also DE, concordant direction)
  </p>
  <div style="display:flex; gap:32px; align-items:flex-start; margin:12px 0;">
    <div class="stat-box">
      <div class="num" style="color:#2c6fad">{n_type1}</div>
      <div class="lbl">Type I (circRNA-specific)</div>
    </div>
    <div class="stat-box">
      <div class="num" style="color:#e07b39">{n_type2}</div>
      <div class="lbl">Type II (gene-level)</div>
    </div>
  </div>
  <div class="type-bar" style="max-width:500px">
    <span style="width:55px;font-size:12px">Type I</span>
    <div class="bar-wrap">
      <div class="bar-fill bar-type1" style="width:{pct1}%">{pct1}%</div>
    </div>
  </div>
  <div class="type-bar" style="max-width:500px">
    <span style="width:55px;font-size:12px">Type II</span>
    <div class="bar-wrap">
      <div class="bar-fill bar-type2" style="width:{pct2}%">{pct2}%</div>
    </div>
  </div>
"""


def _plot_isoform_usage(sig: pd.DataFrame, top_n: int = 10,
                        case_label: str = "tumor",
                        control_label: str = "normal") -> str:
    """Plotly stacked-bar of IUI for top switching genes; falls back to empty string."""
    try:
        import plotly.graph_objects as go  # type: ignore[import]
    except ImportError:
        return (
            "<p style='font-size:13px;color:#888'>"
            "<em>Install plotly for the interactive IUI chart: "
            "<code>pip install plotly</code></em></p>"
        )
    if sig.empty or "gene_name" not in sig.columns:
        return ""

    # Detect IUI column names (dynamic based on labels, fall back to tumor/normal)
    iui_case_col    = f"iui_{case_label}"    if f"iui_{case_label}"    in sig.columns else "iui_tumor"
    iui_control_col = f"iui_{control_label}" if f"iui_{control_label}" in sig.columns else "iui_normal"

    top_genes = (
        sig.groupby("gene_name")["padj_global"]
        .min()
        .sort_values()
        .head(top_n)
        .index.tolist()
    )

    colors = [
        "#1D9E75", "#534AB7", "#BA7517", "#D85A30", "#378ADD",
        "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#17BECF",
    ]
    fig = go.Figure()
    for gene in top_genes:
        rows = sig[sig["gene_name"] == gene].reset_index(drop=True)
        for i, row in rows.iterrows():
            short = str(row["circ_id"]).split(":")[-1]
            fig.add_trace(go.Bar(
                name        = f"{gene} | {short}",
                x           = [[gene, gene], [control_label.title(), case_label.title()]],
                y           = [row.get(iui_control_col, 0), row.get(iui_case_col, 0)],
                marker_color= colors[int(i) % len(colors)],
                legendgroup = gene,
            ))

    fig.update_layout(
        barmode       = "stack",
        title         = f"Isoform Usage Index — {case_label.title()} vs {control_label.title()} (significant switching only)",
        yaxis_title   = "Isoform Usage Index (IUI)",
        height        = 500,
        legend_title  = "Gene | BSJ coords",
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
        margin        = dict(b=80, t=60),
        xaxis         = dict(tickfont=dict(size=11)),
    )
    return fig.to_html(include_plotlyjs=False, full_html=False)


def _isoform_section(switching_file: Optional[str],
                     isoform_file:   Optional[str] = None,
                     circbase_file:  Optional[str] = None,
                     case_label:     str = "tumor",
                     control_label:  str = "normal") -> str:
    """Return HTML block for isoform switching results; empty string if unavailable."""
    if not switching_file or not Path(switching_file).exists():
        return ""
    try:
        sw = pd.read_csv(switching_file, sep="\t")
    except Exception:
        return ""
    if sw.empty or "is_switching" not in sw.columns:
        return ""

    sig = sw[sw["is_switching"] == True].copy()
    n_genes_tested  = sw["gene_id"].nunique()  if "gene_id"  in sw.columns else 0
    n_genes_sig     = sig["gene_id"].nunique() if not sig.empty else 0
    n_events        = len(sig)

    # Merge exon structure annotations (exon_span, region, strand)
    if isoform_file and "circ_id" in sig.columns:
        try:
            iso = pd.read_csv(isoform_file, sep="\t",
                              usecols=lambda c: c in ("circ_id", "gene_name", "strand", "region", "exon_span"))
            new_cols = ["circ_id"] + [c for c in iso.columns if c not in sig.columns and c != "circ_id"]
            if len(new_cols) > 1:
                sig = sig.merge(iso[new_cols], on="circ_id", how="left")
        except Exception:
            pass

    # Merge circBase annotation
    if circbase_file and "circ_id" in sig.columns:
        try:
            cb = pd.read_csv(circbase_file, sep="\t",
                             usecols=lambda c: c in ("circ_id", "circbase_id", "in_circbase"))
            new_cols = ["circ_id"] + [c for c in cb.columns if c not in sig.columns and c != "circ_id"]
            if len(new_cols) > 1:
                sig = sig.merge(cb[new_cols], on="circ_id", how="left")
        except Exception:
            pass

    # Build compact exon_structure column: "e3-e7 · exonic" / "intronic" / "intergenic"
    if "region" in sig.columns:
        def _struct(row: pd.Series) -> str:
            span   = str(row.get("exon_span", "") or "").strip()
            region = str(row.get("region",    "") or "").strip()
            strand = str(row.get("strand",    "") or "").strip()
            parts  = []
            if span:
                parts.append(span)
            if region:
                parts.append(region)
            if strand:
                parts.append(strand)
            return " · ".join(parts) if parts else "—"
        sig["exon_structure"] = sig.apply(_struct, axis=1)

    plot_html = _plot_isoform_usage(sig, case_label=case_label, control_label=control_label)

    iui_case_col    = f"iui_{case_label}"
    iui_control_col = f"iui_{control_label}"
    show_cols = [c for c in
                 ["gene_name", "circ_id", "exon_structure", "circbase_id",
                  iui_control_col, iui_case_col, "delta_iui", "p_value"]
                 if c in sig.columns]
    sort_keys = [c for c in ["gene_name", "p_value"] if c in sig.columns]
    disp = sig.sort_values(sort_keys)[show_cols].copy()
    if "circbase_id" in disp.columns:
        disp["circbase_id"] = disp["circbase_id"].apply(_cb_link)
    table_html = (
        _dl_wrap(_df_to_html(disp, max_rows=40, escape=False),
                 "tbl_isoform", "isoform_switching.csv")
        if show_cols and sort_keys else ""
    )

    return f"""
  <h2>Isoform Switching Analysis</h2>
  <p style="font-size:13px;color:#555;margin-bottom:8px;">
    <b>IUI</b> (Isoform Usage Index) = BSJ<sub>i</sub> / Σ BSJ for all isoforms of the same gene.
    Wilcoxon rank-sum test per isoform; BH-corrected FDR at two levels
    (within-gene and global). Threshold: FDR &lt; 0.05 &amp; |ΔIUI| &gt; 0.10.
  </p>
  <div>
    <div class="stat-box">
      <div class="num">{n_genes_tested}</div>
      <div class="lbl">Multi-isoform genes tested</div>
    </div>
    <div class="stat-box">
      <div class="num" style="color:#d62728">{n_genes_sig}</div>
      <div class="lbl">Genes with switching</div>
    </div>
    <div class="stat-box">
      <div class="num" style="color:#d62728">{n_events}</div>
      <div class="lbl">Switching isoforms</div>
    </div>
  </div>
  {plot_html}
  {table_html}
"""


def _compute_score_dist_data(
    de: pd.DataFrame,
    p_col: str,
    sig_thr: float,
    lfc: float,
    bm_lookup: dict,
    method_label: str = "",
) -> Optional[dict]:
    """Pre-compute score distribution data for one DE method (for JS-side Plotly.react)."""
    if p_col not in de.columns or "log2FC" not in de.columns or "circ_id" not in de.columns:
        return None
    mask = de[p_col].notna() & (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)
    sig = de[mask].copy()
    if sig.empty:
        return None
    n = len(sig)
    sig = sig.copy()
    sig["_sig"] = sig[p_col].apply(lambda p: min(-math.log10(max(p, 1e-10)), 10))
    sig["_fc"]  = sig["log2FC"].abs().clip(upper=5)
    def _mn(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn + 1e-10)
    sig["_sig_n"] = _mn(sig["_sig"])
    sig["_fc_n"]  = _mn(sig["_fc"])
    use_ixn = any(bm_lookup.get(str(r), {}).get("mirna_n", 0) > 0 for r in sig["circ_id"])
    ndim = 6.0 if use_ixn else 4.0
    scores = []
    for _, r in sig.iterrows():
        lu = bm_lookup.get(str(r["circ_id"]), {})
        base = r["_sig_n"] + r["_fc_n"] + lu.get("conf_n", 0) + lu.get("known", 0)
        if use_ixn:
            base += lu.get("mirna_n", 0) + lu.get("rbp_n", 0)
        scores.append(round(float(base) / ndim, 4))
    scores.sort(reverse=True)
    try:
        import numpy as _np
        arr = _np.array(scores)
        mu = float(arr.mean()); sd = float(arr.std())
        x_n = _np.linspace(max(0, arr.min()-0.05), min(1, arr.max()+0.05), 200)
        try:
            from scipy.stats import norm as _norm, shapiro as _shap
            y_n = _norm.pdf(x_n, mu, sd).tolist()
            sw_w, sw_p = (_shap(arr) if n <= 5000 else (None, None))
            sw_w = float(sw_w) if sw_w is not None else None
            sw_p = float(sw_p) if sw_p is not None else None
        except Exception:
            y_n = [(1/(sd*(2*3.14159)**0.5))*math.exp(-0.5*((x-mu)/sd)**2) for x in x_n.tolist()]
            sw_w = sw_p = None
        return {"n": n, "label": method_label, "scores": scores,
                "mu": round(mu,4), "sd": round(sd,4),
                "x_norm": [round(float(x),4) for x in x_n],
                "y_norm": [round(float(y),4) for y in y_n],
                "sw_w": round(sw_w,4) if sw_w else None, "sw_p": sw_p}
    except Exception:
        return None


def _compute_bm_table_data(
    de: pd.DataFrame,
    p_col: str,
    sig_thr: float,
    lfc: float,
    bm_lookup: dict,
    sig_sets_all: Optional[dict] = None,
) -> Optional[dict]:
    """Compute all significant biomarker rows for one DE method (stored in ALL_DE_METHODS[m].bm_table); JS slices to displayed N."""
    if p_col not in de.columns or "log2FC" not in de.columns or "circ_id" not in de.columns:
        return None
    mask = de[p_col].notna() & (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc)
    sig = de[mask].copy()
    if sig.empty:
        return None
    n_total = len(sig)
    sig["_sig_v"] = sig[p_col].apply(lambda p: min(-math.log10(max(p, 1e-10)), 10))
    sig["_fc_v"]  = sig["log2FC"].abs().clip(upper=5)
    def _mn(s):
        mn, mx = s.min(), s.max()
        return (s - mn) / (mx - mn + 1e-10)
    sig["_sig_n"] = _mn(sig["_sig_v"])
    sig["_fc_n"]  = _mn(sig["_fc_v"])
    use_ixn = any(bm_lookup.get(str(r), {}).get("n_mirna", 0) > 0 for r in sig["circ_id"])
    ndim = 6.0 if use_ixn else 4.0
    # Normalize mirna/rbp within the current sig set to avoid cross-set scaling artifacts
    _sig_mirna_mx = max((float(bm_lookup.get(str(c), {}).get("n_mirna", 0)) for c in sig["circ_id"]), default=1) or 1
    _sig_rbp_mx   = max((float(bm_lookup.get(str(c), {}).get("n_rbp",   0)) for c in sig["circ_id"]), default=1) or 1
    scores_list = []
    for _, r in sig.iterrows():
        lu = bm_lookup.get(str(r["circ_id"]), {})
        mirna_n = float(lu.get("n_mirna", 0)) / _sig_mirna_mx if use_ixn else 0.0
        rbp_n   = float(lu.get("n_rbp",   0)) / _sig_rbp_mx   if use_ixn else 0.0
        base = r["_sig_n"] + r["_fc_n"] + lu.get("conf_n", 0) + lu.get("known", 0)
        if use_ixn:
            base += mirna_n + rbp_n
        scores_list.append(round(float(base) / ndim, 4))
    sig["_score"] = scores_list
    all_sig = sig_sets_all or {}
    all_methods = list(all_sig.keys())
    if all_methods:
        sig["_n_sig"] = sig["circ_id"].apply(
            lambda cid: sum(1 for m in all_methods if str(cid) in all_sig[m])
        )
    else:
        sig["_n_sig"] = 1
    sig_top = sig.sort_values("_score", ascending=False).reset_index(drop=True)
    cols = ["rank", "circ_id", "log2FC", "n_mirna", "n_rbp",
            "biomarker_score", "n_sig_methods", "in_circbase", "circbase_id", "circbase_gene", "Type"]
    rows = []
    for i, r in sig_top.iterrows():
        cid = str(r["circ_id"])
        lu = bm_lookup.get(cid, {})
        rows.append([
            int(i) + 1,
            cid,
            round(float(r.get("log2FC", 0)), 2),
            int(lu.get("n_mirna", 0)),
            int(lu.get("n_rbp", 0)),
            round(float(r.get("_score", 0)), 3),
            int(r.get("_n_sig", 1)),
            int(lu.get("in_circbase", 0)),
            lu.get("circbase_id", ""),
            (lambda _g, _n: _g if _g and _g.lower() not in ("none","nan","novel","intergenic","") else (_n if _n and _n.lower() not in ("none","nan","novel","intergenic","") else ""))(lu.get("circbase_gene",""), lu.get("gene_name","")),
            (lambda _tv, _lu: (
                str(_tv) if (_tv is not None and _tv == _tv and str(_tv).lower() not in ("nan","none","na",""))
                else (_lu.get("type_edger") or "—")
            ))(r.get("Type") if "Type" in r.index else None, lu),
        ])
    return {"cols": cols, "rows": rows, "n_total": n_total}


def _biomarker_score_dist(bm: pd.DataFrame) -> str:
    """Ranked scatter plot of biomarker scores for all DE circRNAs."""
    if not _PLOTLY or "biomarker_score" not in bm.columns or "circ_id" not in bm.columns:
        return ""
    import plotly.graph_objects as go

    bm2 = bm.reset_index(drop=True).copy()
    bm2["_rank"] = range(1, len(bm2) + 1)
    top_n = 30

    def _hover(r):
        parts = [f"<b>#{int(r['_rank'])}</b>  {r['circ_id']}",
                 f"Score: {r['biomarker_score']:.4f}"]
        if "log2FC" in r and pd.notna(r["log2FC"]):
            parts.append(f"log₂FC: {r['log2FC']:.3f}")
        if "Type" in r and pd.notna(r["Type"]):
            parts.append(f"Type: {r['Type']}")
        if "in_circbase" in r and r["in_circbase"]:
            cb = r.get("circbase_id", "")
            parts.append(f"circBase: {cb}")
        return "<br>".join(parts)

    bm2["_hover"] = bm2.apply(_hover, axis=1)
    top = bm2.iloc[:top_n]
    rest = bm2.iloc[top_n:]

    fig = go.Figure()
    if not rest.empty:
        fig.add_trace(go.Scatter(
            x=rest["_rank"], y=rest["biomarker_score"],
            mode="markers", name=f"Rank {top_n+1}–{len(bm2)}",
            marker=dict(color="rgba(120,120,120,0.45)", size=5),
            hovertext=rest["_hover"], hoverinfo="text",
        ))
    fig.add_trace(go.Scatter(
        x=top["_rank"], y=top["biomarker_score"],
        mode="markers", name=f"Top {top_n} (table)",
        marker=dict(color="#d62728", size=7, line=dict(color="white", width=0.8)),
        hovertext=top["_hover"], hoverinfo="text",
    ))
    # Elbow line
    fig.add_shape(type="line", x0=top_n + 0.5, x1=top_n + 0.5,
                  y0=0, y1=1, yref="paper",
                  line=dict(dash="dot", color="#d62728", width=1.2))
    fig.add_annotation(x=top_n + 1.5, y=0.98, yref="paper",
                       text=f"Top {top_n}", showarrow=False,
                       font=dict(size=10, color="#d62728"), xanchor="left")

    score_min = float(bm2["biomarker_score"].min())
    score_max = float(bm2["biomarker_score"].max())
    fig.update_layout(
        xaxis=dict(title="Rank", tickmode="auto"),
        yaxis=dict(title="Biomarker Score",
                   range=[max(0, score_min - 0.05), min(1, score_max + 0.05)]),
        height=320,
        margin=dict(t=20, b=60, l=70, r=40),
        legend=dict(x=0.75, y=0.95, bgcolor="rgba(255,255,255,0.8)"),
        plot_bgcolor="white", paper_bgcolor="white",
        hovermode="closest",
    )
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id="bm-scatter-plot")


def _biomarker_normality_plot(bm: pd.DataFrame) -> str:
    """Histogram + fitted normal curve + Shapiro-Wilk test for biomarker score distribution."""
    if not _PLOTLY or "biomarker_score" not in bm.columns:
        return ""
    import plotly.graph_objects as go
    import numpy as np

    scores = bm["biomarker_score"].dropna().values
    n = len(scores)
    if n == 0:
        return ""
    mu, sd = float(scores.mean()), float(scores.std())

    # Shapiro-Wilk (n ≤ 5000); fallback to KS for larger sets
    sw_text = ""
    conclusion = ""
    try:
        from scipy import stats as _stats
        if n < 3:
            sw_text = f"樣本數不足（n={n}），無法進行常態檢定 / Insufficient samples (n={n}) for normality test"
            conclusion = ""
        elif n <= 5000:
            W, p_sw = _stats.shapiro(scores)
            sw_text = f"Shapiro-Wilk: W = {W:.4f}, p = {p_sw:.4e}"
            conclusion = ("✗ 非常態 / Non-normal" if p_sw < 0.05 else "✓ 常態 / Normal") + " (α = 0.05)"
        else:
            D, p_sw = _stats.kstest(scores, "norm", args=(mu, sd))
            sw_text = f"KS test: D = {D:.4f}, p = {p_sw:.4e}"
            conclusion = ("✗ 非常態 / Non-normal" if p_sw < 0.05 else "✓ 常態 / Normal") + " (α = 0.05)"
    except ImportError:
        sw_text = "scipy 未安裝，無法執行常態檢定 / scipy not installed"

    # Histogram
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=scores, nbinsx=30, name="Observed distribution",
        marker_color="rgba(44,119,214,0.55)",
        marker_line=dict(color="rgba(44,119,214,0.9)", width=0.8),
        histnorm="probability density",
    ))

    # Fitted normal curve
    x_norm = np.linspace(scores.min() - 0.05, scores.max() + 0.05, 300)
    try:
        from scipy.stats import norm as _norm
        y_norm = _norm.pdf(x_norm, mu, sd)
    except ImportError:
        y_norm = (1 / (sd * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_norm - mu) / sd) ** 2)

    fig.add_trace(go.Scatter(
        x=x_norm.tolist(), y=y_norm.tolist(),
        mode="lines", name=f"Normal(μ={mu:.3f}, σ={sd:.3f})",
        line=dict(color="#d62728", width=2),
    ))

    # Mean ± 1SD / ± 2SD reference lines
    ref_lines = [
        (mu,        "μ",    "solid", "#555", 1.2),
        (mu - sd,   "μ−σ",  "dot",   "#888", 1.0),
        (mu + sd,   "μ+σ",  "dot",   "#888", 1.0),
        (mu - 2*sd, "μ−2σ", "dash",  "#bbb", 0.8),
        (mu + 2*sd, "μ+2σ", "dash",  "#bbb", 0.8),
    ]
    for xv, lbl, dash, col, lw in ref_lines:
        fig.add_shape(type="line", x0=xv, x1=xv, y0=0, y1=1, yref="paper",
                      line=dict(color=col, width=lw, dash=dash))
        fig.add_annotation(x=xv, y=1.02, yref="paper", text=lbl, showarrow=False,
                           font=dict(size=9, color=col), xanchor="center")

    sw_color = "#d62728" if "Non-normal" in conclusion else "#2CA02C"
    fig.update_layout(
        xaxis=dict(title="Biomarker Score"),
        yaxis=dict(title="Probability Density"),
        height=320,
        margin=dict(t=40, b=90, l=70, r=40),
        legend=dict(x=0.65, y=0.95, bgcolor="rgba(255,255,255,0.85)"),
        plot_bgcolor="white", paper_bgcolor="white",
        annotations=[dict(
            x=0.5, y=-0.28, xref="paper", yref="paper",
            text=f"{sw_text}    <b style='color:{sw_color}'>{conclusion}</b>",
            showarrow=False, align="center", xanchor="center", yanchor="top",
            bgcolor="rgba(255,255,255,0.88)", bordercolor="#ccc", borderwidth=1,
            font=dict(size=10),
        )],
    )
    chart_html = fig.to_html(include_plotlyjs=False, full_html=False, div_id="bm-hist-plot")
    caption = (
        "<p style='font-size:11px;color:#888;margin:2px 0 0;line-height:1.6'>"
        "Vertical lines: <b>solid</b> = μ ({mu:.3f}); "
        "<b>dotted</b> = μ ± σ ({lo1:.3f} – {hi1:.3f}, ~68% of data); "
        "<b>dashed</b> = μ ± 2σ ({lo2:.3f} – {hi2:.3f}, ~95% of data)"
        "</p>"
    ).format(mu=mu, lo1=mu - sd, hi1=mu + sd, lo2=mu - 2*sd, hi2=mu + 2*sd)
    return chart_html + caption


def _biomarker_section(biomarker_file: Optional[str],
                       interactions: Optional[dict] = None,
                       iso_lookup: Optional[dict] = None) -> str:
    if not biomarker_file or not Path(biomarker_file).exists():
        return ""
    try:
        bm = pd.read_csv(biomarker_file, sep="\t")
    except Exception:
        return ""
    if bm.empty:
        return ""
    show_cols = [c for c in ["rank", "circ_id", "log2FC", "n_mirna", "n_rbp",
                              "biomarker_score", "n_sig_methods",
                              "in_circbase", "circbase_id", "circbase_gene", "Type"]
                 if c in bm.columns]
    has_interactions = "n_mirna" in bm.columns and (bm["n_mirna"].max() + bm["n_rbp"].max()) > 0
    if has_interactions:
        score_desc = "Score = mean of: −log₁₀(p-value), |log₂FC|, confidence score, circBase bonus, #miRNA binders, #RBP binders (each normalised 0–1; 6D)."
    else:
        score_desc = "Score = mean of: −log₁₀(p-value), |log₂FC|, confidence score, circBase known bonus (each normalised 0–1; 4D)."

    dist_html = _biomarker_score_dist(bm)
    norm_html = _biomarker_normality_plot(bm)

    n_total = len(bm)
    bm_top = bm.head(30).copy()
    # Fix circbase_gene: replace "None"/"nan"/"novel" with gene_name from iso_lookup
    if "circbase_gene" in bm_top.columns and iso_lookup:
        _BAD_GENE = {"none", "nan", "novel", "intergenic", ""}
        def _fix_cbg(row):
            cg = str(row.get("circbase_gene", "") or "")
            if cg.lower() in _BAD_GENE:
                gn = iso_lookup.get(str(row.get("circ_id", "")), {}).get("gene_name", "") if iso_lookup else ""
                return gn if gn and gn.lower() not in _BAD_GENE else "—"
            return cg
        bm_top["circbase_gene"] = bm_top.apply(_fix_cbg, axis=1)
    disp = _fmt_floats(bm_top[show_cols].copy())
    if "circ_id" in disp.columns and interactions is not None:
        def _link(v: str) -> str:
            tip = "in interaction data" if str(v) in interactions else "no interaction data pre-fetched"
            return (f'<a class="circ-link" onclick="showCircDetail(\'{v}\')" '
                    f'title="{tip}">{v}</a>')
        disp["circ_id"] = disp["circ_id"].astype(str).apply(_link)
    if "circbase_id" in disp.columns:
        disp["circbase_id"] = disp["circbase_id"].apply(_cb_link)

    raw_html = disp.to_html(index=False, classes="table", border=0, na_rep="—", escape=False)
    # Inject data-nsig into body <tr> tags.
    # pandas renders the header as <tr style="text-align: right;"> which does NOT
    # match plain <tr>, so every matched <tr> is a body row — no header skip needed.
    if "n_sig_methods" in bm_top.columns:
        import re as _re
        nsig_vals = bm_top["n_sig_methods"].tolist()
        _idx = [0]
        def _inject(m):
            v = nsig_vals[_idx[0]] if _idx[0] < len(nsig_vals) else 1
            _idx[0] += 1
            fw = 'font-weight:bold;' if int(v) >= 3 else ''
            return f'<tr data-nsig="{v}" style="{fw}">'
        raw_html = _re.sub(r'<tr>', _inject, raw_html)

    table_html = _dl_wrap(raw_html, "tbl_biomarker", "biomarker_candidates.csv")

    # Build filter UI (only when n_sig_methods column exists)
    filter_ui = ""
    if "n_sig_methods" in bm.columns:
        n3 = int((bm.head(30)["n_sig_methods"] == 3).sum())
        n2 = int((bm.head(30)["n_sig_methods"] >= 2).sum())
        filter_ui = f"""
  <div id="bm-filter-bar" style="display:flex;align-items:center;gap:8px;margin:10px 0 6px;flex-wrap:wrap">
    <span style="font-size:12px;color:#666" data-en="Show:">顯示：</span>
    <button class="bm-filter-btn active" onclick="filterBiomarker(0)"
      style="font-size:12px;padding:3px 12px;border:1px solid #bbb;border-radius:12px;
             background:#2c6fad;color:white;cursor:pointer" data-en="All ({min(n_total,30)})" data-count="{min(n_total,30)}" data-filter="0">全部（{min(n_total,30)}）</button>
    <button class="bm-filter-btn" onclick="filterBiomarker(2)"
      style="font-size:12px;padding:3px 12px;border:1px solid #bbb;border-radius:12px;
             background:white;color:#555;cursor:pointer" data-en="≥ 2 methods ({n2})" data-count="{n2}" data-filter="2">≥ 2 方法顯著（{n2}）</button>
    <button class="bm-filter-btn" onclick="filterBiomarker(3)"
      style="font-size:12px;padding:3px 12px;border:1px solid #bbb;border-radius:12px;
             background:white;color:#555;cursor:pointer" data-en="All 3 methods ({n3})" data-count="{n3}" data-filter="3">3 方法均顯著（{n3}）</button>
    <span style="font-size:11px;color:#999;margin-left:4px" data-en="n_sig_methods = significant in 1/2/3 DE methods">n_sig_methods = 1/2/3 種 DE 方法中顯著</span>
  </div>
  <script>
  function filterBiomarker(minN) {{
    document.querySelectorAll('.bm-filter-btn').forEach(function(b) {{
      b.style.background = 'white'; b.style.color = '#555';
    }});
    event.currentTarget.style.background = '#2c6fad';
    event.currentTarget.style.color = 'white';
    var rows = document.querySelectorAll('#tbl_biomarker tbody tr');
    rows.forEach(function(r) {{
      var n = parseInt(r.getAttribute('data-nsig') || '1');
      r.style.display = (minN === 0 || n >= minN) ? '' : 'none';
    }});
  }}
  </script>"""

    dist_section = ""
    if dist_html or norm_html:
        dist_section = f"<h3 id='bm-dist-title' style='font-size:14px;color:#444;margin:16px 0 4px'>Score Distribution — all {n_total} significant DE circRNAs</h3>"
        if dist_html and norm_html:
            dist_section += (
                "<div style='display:flex;gap:16px;flex-wrap:wrap'>"
                f"<div style='flex:1;min-width:320px'><p style='font-size:12px;color:#666;margin:0 0 4px'>Ranked Score</p>{dist_html}</div>"
                f"<div style='flex:1;min-width:320px'><p style='font-size:12px;color:#666;margin:0 0 4px'>Histogram + Normal Fit</p>{norm_html}</div>"
                "</div>"
            )
        else:
            dist_section += dist_html or norm_html

    return f"""
  <h2 id="biomarker-section-title">Biomarker Candidates (top 30 by composite score)</h2>
  <p style="font-size:13px;color:#555;margin-bottom:8px;">{score_desc}</p>
  {dist_section}
  <div style="display:flex;align-items:center;gap:10px;margin:20px 0 4px;flex-wrap:wrap">
    <h3 style="font-size:14px;color:#444;margin:0">Biomarker Candidates</h3>
    <div style="display:flex;align-items:center;gap:6px;background:#f4f8ff;padding:6px 12px;border-radius:6px;border:1px solid #d0e4f7">
      <span style="font-size:13px" data-en="Show top">顯示前</span>
      <input type="number" id="bm-n-input" value="30" min="1" max="{n_total}"
             style="width:64px;padding:3px 6px;border:1px solid #bbb;border-radius:4px;font-size:13px"
             onkeydown="if(event.key==='Enter')updateBiomarkerN()">
      <button onclick="updateBiomarkerN()"
              style="background:#2c6fad;color:white;border:none;border-radius:4px;padding:4px 14px;cursor:pointer;font-size:13px"
              data-en="Update">更新</button>
      <span id="bm-n-status" style="font-size:12px;color:#888"></span>
    </div>
  </div>
  <p style="font-size:12px;color:#666">&#128204; Click a <strong>circ_id</strong> to view exon diagram, miRNA and RBP binding sites.</p>
  {filter_ui}
  {table_html}
"""


def _plotly_volcano(de: pd.DataFrame, fdr: float, lfc: float, de_method: str,
                    p_col: str = "padj", sig_thr: float = 0.05, p_label: str = "adjusted p-value",
                    heatmap_ids: Optional[set] = None,
                    div_id: str = "main-volcano-plot") -> str:
    """Interactive Plotly volcano; returns '' when plotly is unavailable or DE is empty."""
    if not _PLOTLY:
        return ""
    if p_col not in de.columns:
        p_col = "padj" if "padj" in de.columns else None
    if p_col is None:
        return ""
    df = de[de[p_col].notna() & de["log2FC"].notna()].copy()
    if df.empty:
        return ""

    df["nlp"] = df[p_col].clip(lower=1e-300).apply(lambda p: -math.log10(p))

    def _sig(row):
        if row[p_col] < sig_thr and row["log2FC"] > lfc:
            return "Up"
        if row[p_col] < sig_thr and row["log2FC"] < -lfc:
            return "Down"
        return "NS"

    df["sig"] = df.apply(_sig, axis=1)
    has_type = "Type" in df.columns

    colors = {"Up": "#d62728", "Down": "#2CA02C", "NS": "rgba(150,150,150,0.35)"}
    sizes  = {"Up": 6, "Down": 6, "NS": 4}

    traces = []
    for grp in ["NS", "Down", "Up"]:
        sub = df[df["sig"] == grp]
        if sub.empty:
            continue
        hover = sub.apply(
            lambda r: (
                f"<b>{r.get('circ_id', '')}</b><br>"
                f"log₂FC: {r['log2FC']:.3f}<br>"
                f"{p_label}: {r[p_col]:.2e}"
                + (f"<br>Type: {r['Type']}" if has_type and pd.notna(r.get("Type")) else "")
            ), axis=1,
        )
        traces.append(go.Scatter(
            x=sub["log2FC"].tolist(), y=sub["nlp"].tolist(),
            mode="markers", name=grp,
            marker=dict(color=colors[grp], size=sizes[grp], opacity=0.8,
                        line=dict(width=0)),
            text=hover.tolist(),
            hovertemplate="%{text}<extra></extra>",
        ))

    # Overlay heatmap top circRNAs with open circle markers + labels
    if heatmap_ids:
        heat_df = df[df["circ_id"].isin(heatmap_ids)]
        if not heat_df.empty:
            heat_hover = heat_df.apply(
                lambda r: (
                    f"<b>{r.get('circ_id','')}</b><br>"
                    f"★ Heatmap top {int(len(heatmap_ids)//2)} up+down<br>"
                    f"log₂FC: {r['log2FC']:.3f}<br>"
                    f"{p_label}: {r[p_col]:.2e}"
                ), axis=1,
            )
            traces.append(go.Scatter(
                x=heat_df["log2FC"].tolist(), y=heat_df["nlp"].tolist(),
                mode="markers+text",
                name="Heatmap top",
                visible="legendonly",
                marker=dict(symbol="circle-open", size=13, color="black",
                            line=dict(width=2, color="black")),
                text=heat_df["circ_id"].tolist(),
                textposition="top center",
                textfont=dict(size=8, color="black"),
                hovertext=heat_hover.tolist(),
                hovertemplate="%{hovertext}<extra></extra>",
            ))

    fig = go.Figure(traces)
    fig.add_vline(x=lfc,   line_dash="dot", line_color="#aaa", line_width=1)
    fig.add_vline(x=-lfc,  line_dash="dot", line_color="#aaa", line_width=1)
    fig.add_hline(y=-math.log10(sig_thr), line_dash="dot", line_color="#aaa", line_width=1)

    # toggle button for heatmap top circles (only when trace exists)
    heatmap_trace_idx = len(traces) - 1 if heatmap_ids else None
    updatemenus = []
    if heatmap_trace_idx is not None:
        updatemenus = [dict(
            type="buttons", direction="left", showactive=True,
            x=1.0, y=1.13, xanchor="right", yanchor="top",
            bgcolor="#f4f8ff", bordercolor="#99c2f0", borderwidth=1,
            font=dict(size=11),
            buttons=[
                dict(label="○ Heatmap circles: OFF",
                     method="restyle",
                     args=[{"visible": "legendonly"}, [heatmap_trace_idx]]),
                dict(label="◉ Heatmap circles: ON",
                     method="restyle",
                     args=[{"visible": True}, [heatmap_trace_idx]]),
            ],
        )]

    fig.update_layout(
        title=dict(text=f"Volcano Plot [{de_method}]", font_size=14),
        xaxis_title="log₂ Fold Change (Tumor / Normal)",
        yaxis_title=f"−log₁₀({p_label})",
        height=500, plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(title="", orientation="h", y=1.02, x=0),
        margin=dict(t=80),
        updatemenus=updatemenus,
    )
    fig.update_xaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=False)
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id=div_id)


def _plotly_heatmap(de: pd.DataFrame, matrix: pd.DataFrame, top_n: int = 10,
                    p_col: str = "padj",
                    groups_file: Optional[str] = None,
                    normal_label: str = "normal") -> str:
    """Interactive Plotly heatmap: top N up + N down DE circRNAs (normal-centered z-score)."""
    if not _PLOTLY:
        return ""
    if p_col not in de.columns:
        p_col = "padj" if "padj" in de.columns else "pvalue"
    if "circ_id" not in de.columns or p_col not in de.columns or "log2FC" not in de.columns:
        return ""

    # circbase label map (de already enriched by _enrich_de)
    cb_map = {}
    if "circbase_id" in de.columns and "in_circbase" in de.columns:
        known = de[de["in_circbase"] == 1].dropna(subset=["circbase_id"])
        cb_map = dict(zip(known["circ_id"], known["circbase_id"]))

    def _label(cid: str) -> str:
        cb = cb_map.get(cid, "")
        return cb if cb and cb not in ("", "novel") else cid

    up_ids = (
        de[de["log2FC"] > 0].dropna(subset=[p_col])
        .sort_values(p_col).head(top_n)["circ_id"].tolist()
    )
    dn_ids = (
        de[de["log2FC"] < 0].dropna(subset=[p_col])
        .sort_values(p_col).head(top_n)["circ_id"].tolist()
    )
    top_ids = up_ids + dn_ids

    avail = [i for i in top_ids if i in matrix.index]
    if len(avail) < 2:
        return ""

    sub = matrix.loc[avail].astype(float)
    log_sub = np.log2(sub + 1)

    # Normal-centered z-score: use normal group mean/SD as reference
    normal_cols = []
    cond_map_hm: dict = {}
    if groups_file and Path(groups_file).exists():
        try:
            grp = pd.read_csv(groups_file)
            cond_map_hm = dict(zip(grp["srr_id"], grp["condition"]))
            normal_cols = [c for c in log_sub.columns if cond_map_hm.get(c, "") == normal_label]
        except Exception:
            pass
    ref_cols = normal_cols if normal_cols else log_sub.columns.tolist()
    row_mean = log_sub[ref_cols].mean(axis=1)       # center on normal mean
    row_std  = log_sub.std(axis=1).clip(lower=0.1)  # scale on all-sample SD
    z = log_sub.sub(row_mean, axis=0).div(row_std, axis=0)

    # Reorder columns: tumor (non-normal) first, then normal
    _t_cols = [c for c in z.columns if cond_map_hm.get(c, "") not in (normal_label, "")]
    _n_cols = [c for c in z.columns if cond_map_hm.get(c, "") == normal_label]
    _o_cols = [c for c in z.columns if c not in set(_t_cols) and c not in set(_n_cols)]
    if _t_cols or _n_cols:
        z = z[_t_cols + _n_cols + _o_cols]

    n_up = sum(i in up_ids for i in avail)
    n_dn = sum(i in dn_ids for i in avail)
    y_labels = [_label(i) for i in avail]

    fig = go.Figure(go.Heatmap(
        z=z.values.tolist(), x=z.columns.tolist(), y=y_labels,
        colorscale=[[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']], zmid=0,
        colorbar=dict(title=f"z-score<br>({normal_label}-<br>centered)"),
        hovertemplate="<b>%{y}</b><br>%{x}<br>z-score: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=""),
        yaxis=dict(tickfont=dict(size=8), autorange="reversed"),
        height=max(420, len(avail) * 22 + 130),
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(t=70, l=300),
    )
    return fig.to_html(include_plotlyjs=False, full_html=False, div_id="main-heatmap-plot")


def _plotly_pca(matrix: pd.DataFrame, groups_file: Optional[str] = None) -> str:
    """Interactive Plotly PCA coloured by condition (tumor/normal)."""
    if not _PLOTLY or matrix.shape[1] < 2:
        return ""

    condition_map = {}  # type: ignore[var-annotated]
    if groups_file and Path(groups_file).exists():
        try:
            grp = pd.read_csv(groups_file)
            if {"srr_id", "condition"} <= set(grp.columns):
                condition_map = dict(zip(grp["srr_id"], grp["condition"]))
        except Exception:
            pass

    log_mat = np.log2(matrix.values.astype(float) + 1)
    X = log_mat.T
    X -= X.mean(axis=0)
    try:
        U, s, _ = np.linalg.svd(X, full_matrices=False)
    except np.linalg.LinAlgError:
        return ""

    scores  = U * s
    var_pct = (s ** 2 / max((s ** 2).sum(), 1e-10)) * 100

    samples    = matrix.columns.tolist()
    conditions = [condition_map.get(smp, "unknown") for smp in samples]
    fallback   = ["#d62728", "#2CA02C", "#2563eb", "#e07b39", "#9333ea", "#dc2626"]
    uniq_conds = sorted(set(conditions) - {"unknown"})
    color_map: dict = {"unknown": "#888"}
    for i, c in enumerate(uniq_conds):
        color_map[c] = fallback[i % len(fallback)]

    fig = go.Figure()
    for cond in sorted(set(conditions)):
        idx = [i for i, c in enumerate(conditions) if c == cond]
        fig.add_trace(go.Scatter(
            x=[float(scores[i, 0]) for i in idx],
            y=[float(scores[i, 1]) for i in idx],
            mode="markers+text", name=cond,
            text=[samples[i] for i in idx], textposition="top center",
            textfont=dict(size=10),
            marker=dict(color=color_map[cond], size=12,
                        line=dict(color="white", width=1.5)),
            hovertemplate="%{text}<extra></extra>",
        ))
    fig.update_layout(
        title=dict(text="PCA (log₂ normalized)", font_size=14),
        xaxis_title=f"PC1: {var_pct[0]:.1f}%",
        yaxis_title=f"PC2: {var_pct[1]:.1f}%",
        height=440, plot_bgcolor="white", paper_bgcolor="white",
        legend_title="Condition", margin=dict(t=60),
    )
    fig.update_xaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=True, zerolinecolor="#ccc")
    fig.update_yaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=True, zerolinecolor="#ccc")
    return fig.to_html(include_plotlyjs=False, full_html=False)


def _enrich_de(de: pd.DataFrame,
               isoform_file: Optional[str],
               circbase_file: Optional[str]) -> pd.DataFrame:
    """Merge host gene name, circBase ID, and exon span into the DE table."""
    if isoform_file:
        try:
            iso = pd.read_csv(isoform_file, sep="\t",
                              usecols=lambda c: c in
                              ("circ_id", "gene_name", "strand", "region", "exon_span"))
            de = de.merge(iso, on="circ_id", how="left")
        except Exception:
            pass
    if circbase_file:
        try:
            cb = pd.read_csv(circbase_file, sep="\t",
                             usecols=lambda c: c in
                             ("circ_id", "circbase_id", "circbase_gene", "in_circbase"))
            de = de.merge(cb, on="circ_id", how="left")
        except Exception:
            pass
    return de


def build_report(
    de_file:        str,
    matrix_file:    str,
    volcano_pdf:    str,
    heatmap_pdf:    str,
    pca_pdf:        str,
    output_file:    str,
    project_id:     str,
    fdr:            float = 0.05,
    lfc:            float = 1.0,
    de_method:      str   = "deseq2",
    biomarker_file: Optional[str] = None,
    switching_file: Optional[str] = None,
    groups_file:    Optional[str] = None,
    de_sig_by:      str   = "auto",
    tumor_label:    str   = "tumor",
    normal_label:   str   = "normal",
    isoform_file:      Optional[str] = None,
    circbase_file:     Optional[str] = None,
    heatmap_top_n:     int   = 10,
    interactions_file: Optional[str] = None,
    multiqc_file:      Optional[str] = None,
    de_files:          dict  = {},
    study_title:       str   = "",
) -> None:
    de     = pd.read_csv(de_file, sep="\t")
    matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)

    # Normalise column names — analysis.R always outputs log2FC
    if "log2FC" not in de.columns and "log2FoldChange" in de.columns:
        de = de.rename(columns={"log2FoldChange": "log2FC"})

    # Enrich with host gene, circBase ID, exon span
    de = _enrich_de(de, isoform_file, circbase_file)

    p_col, sig_thr, _sig_label_str = _eff_sig(de, de_sig_by, fdr)
    sig_mask = (de[p_col] < sig_thr) & (de["log2FC"].abs() > lfc) if p_col in de.columns else pd.Series(False, index=de.index)
    sig: pd.DataFrame = de.loc[sig_mask]

    n_total  = len(matrix)
    n_sig    = len(sig)
    n_up     = int((sig["log2FC"] > 0).sum()) if len(sig) else 0
    n_dn     = int((sig["log2FC"] < 0).sum()) if len(sig) else 0
    n_sample = matrix.shape[1]
    sig_label = f"{_sig_label_str}&lt;{sig_thr}"

    # Rename circ_id to circ_position for display
    sig = sig.copy()
    if "circ_id" in sig.columns:
        sig = sig.rename(columns={"circ_id": "circ_position"})

    top_cols = [c for c in [
        "circ_position", "gene_name", "strand", "region", "exon_span", "circbase_id",
        "log2FC", "pvalue", "Type",
    ] if c in sig.columns]
    top_table = sig.sort_values(p_col)[top_cols] if top_cols else sig.head(20)

    # Load pre-fetched interaction data
    interactions = _load_interactions(interactions_file)

    # Sync isoform/circbase metadata from the enriched DE table into interactions so that
    # the Circular Structure modal displays the same gene name as the table.
    # circbase_gene (e.g. "KIAA0182") takes priority over gene_name (e.g. "GSE1 / KIAA0182
    # renamed in 2019) because the biomarker/DE tables display circbase_gene in their columns.
    _iso_lookup: dict = {}
    if "circ_id" in de.columns:
        for _, _r in de.iterrows():
            _cid = str(_r["circ_id"])
            if _cid not in _iso_lookup:
                _cb_gene = str(_r.get("circbase_gene", "") or "")
                _gname   = str(_r.get("gene_name",    "") or "")
                _iso_lookup[_cid] = {
                    "gene_name":    _cb_gene if _cb_gene and _cb_gene not in ("nan", "None", "novel") else _gname,
                    "strand":       str(_r.get("strand",    "") or ""),
                    "region":       str(_r.get("region",    "") or ""),
                    "exon_span":    str(_r.get("exon_span", "") or ""),
                    "circbase_gene": _cb_gene,
                }
    for _cid, _entry in interactions.items():
        if _cid in _iso_lookup:
            _entry.setdefault("info", {}).update(_iso_lookup[_cid])

    import json as _json
    interactions_js = _json.dumps(interactions, ensure_ascii=False)

    # Build compact volcano data for mini-plot in modal
    _p_col_v = p_col
    _vol_rows = []
    if "log2FC" in de.columns and _p_col_v in de.columns and "circ_id" in de.columns:
        for _, _r in de.iterrows():
            if pd.isna(_r["log2FC"]) or pd.isna(_r.get(_p_col_v)):
                continue
            _x = round(float(_r["log2FC"]), 3)
            _pv = float(_r[_p_col_v])
            _y = round(-math.log10(max(_pv, 1e-300)), 3)
            _s = ("U" if (_pv < fdr and _r["log2FC"] > lfc) else
                  "D" if (_pv < fdr and _r["log2FC"] < -lfc) else "N")
            _vol_rows.append([_x, _y, _s, str(_r["circ_id"])])
    volcano_data_js = _json.dumps(_vol_rows, ensure_ascii=False)
    _fdr_js  = sig_thr
    _lfc_js  = lfc
    _sig_label_js = _sig_label_str

    # Build FULL_HEATMAP_DATA for dynamic main heatmap + modal mini-heatmap (pool: top 50 up + 50 down)
    _HM_POOL = 50
    full_heatmap_data_js = "null"
    if not de.empty and not matrix.empty:
        try:
            _pc_hm = p_col
            if {"circ_id", "log2FC", _pc_hm}.issubset(de.columns):
                _sig_hm = de[sig_mask] if "log2FC" in de.columns else de
                _up_pool = _sig_hm[_sig_hm["log2FC"] > 0].dropna(subset=[_pc_hm]).sort_values(_pc_hm).head(_HM_POOL)
                _dn_pool = _sig_hm[_sig_hm["log2FC"] < 0].dropna(subset=[_pc_hm]).sort_values(_pc_hm).head(_HM_POOL)
                _hm_ids  = list(_up_pool["circ_id"]) + list(_dn_pool["circ_id"])
                _hm_avail = [i for i in _hm_ids if i in matrix.index]
                if _hm_avail:
                    _samps  = matrix.columns.tolist()
                    _sub_hm = matrix.loc[_hm_avail].astype(float)
                    _log_hm = (_sub_hm + 1).apply(lambda col: col.apply(lambda v: math.log2(v) if v > 0 else 0.0))
                    _cmap_hm = {}
                    if groups_file and Path(groups_file).exists():
                        try:
                            _grp_df = pd.read_csv(groups_file)
                            _cmap_hm = dict(zip(_grp_df["srr_id"].astype(str), _grp_df["condition"].astype(str)))
                        except Exception:
                            pass
                    _ncols = [c for c in _samps if _cmap_hm.get(c, "") == normal_label]
                    _rcols = _ncols if _ncols else _samps
                    _rmean = _log_hm[_rcols].mean(axis=1)
                    _rstd  = _log_hm.std(axis=1).clip(lower=0.1)
                    _z_hm  = _log_hm.sub(_rmean, axis=0).div(_rstd, axis=0)
                    # Reorder: tumor (non-normal) first, then normal
                    _t_samps = [c for c in _samps if _cmap_hm.get(c, "") not in (normal_label, "")]
                    _n_samps = [c for c in _samps if _cmap_hm.get(c, "") == normal_label]
                    _o_samps = [c for c in _samps if c not in set(_t_samps) and c not in set(_n_samps)]
                    if _t_samps or _n_samps:
                        _samps = _t_samps + _n_samps + _o_samps
                        _z_hm  = _z_hm[_samps]
                    # circbase label map
                    _cb_map_hm = {}
                    if "circbase_id" in de.columns and "in_circbase" in de.columns:
                        _known_cb = de[de["in_circbase"] == 1].dropna(subset=["circbase_id"])
                        _cb_map_hm = dict(zip(_known_cb["circ_id"].astype(str), _known_cb["circbase_id"].astype(str)))
                    _pval_map_hm = dict(zip(de["circ_id"].astype(str), de[_pc_hm]))
                    _lfc_map_hm  = dict(zip(de["circ_id"].astype(str), de["log2FC"]))
                    _rows_hm = {}
                    for _cid in _hm_avail:
                        _cb = _cb_map_hm.get(_cid, "")
                        _lbl = _cb if _cb and _cb not in ("", "novel") else _cid
                        _rows_hm[_cid] = {
                            "z":      [round(v, 3) for v in _z_hm.loc[_cid].tolist()],
                            "pval":   float(_pval_map_hm.get(_cid, 1.0)),
                            "log2fc": float(_lfc_map_hm.get(_cid, 0.0)),
                            "label":  _lbl,
                        }
                    _up_avail = [i for i in _up_pool["circ_id"] if i in set(_hm_avail)]
                    _dn_avail = [i for i in _dn_pool["circ_id"] if i in set(_hm_avail)]
                    full_heatmap_data_js = _json.dumps({
                        "samples":    _samps,
                        "conditions": _cmap_hm,
                        "up_order":   _up_avail,
                        "dn_order":   _dn_avail,
                        "rows":       _rows_hm,
                    }, ensure_ascii=False)
        except Exception as _hm_exc:
            import sys as _sys
            print(f"[report] FULL_HEATMAP_DATA build failed: {_hm_exc}", file=_sys.stderr)

    # Build CLUST_HEATMAP_DATA: all significant circRNAs (primary method), hierarchically clustered rows
    clust_heatmap_data_js = "null"
    if not de.empty and not matrix.empty:
        try:
            _sig_all = de[sig_mask].copy() if "log2FC" in de.columns else de.copy()
            _cids_all = [c for c in _sig_all["circ_id"].tolist() if c in matrix.index]
            if len(_cids_all) >= 3:
                _sub_c = matrix.loc[_cids_all].astype(float)
                _log_c = (_sub_c + 1).apply(lambda col: col.apply(lambda v: math.log2(v) if v > 0 else 0.0))
                # Use same normal-sample mean centering as main heatmap
                _rcols_c = [c for c in _samps if _cmap_hm.get(c, "") == normal_label] if "_samps" in dir() and "_cmap_hm" in dir() else _log_c.columns.tolist()
                if not _rcols_c:
                    _rcols_c = _log_c.columns.tolist()
                _rmean_c = _log_c[_rcols_c].mean(axis=1)
                _rstd_c  = _log_c.std(axis=1).clip(lower=0.1)
                _z_c     = _log_c.sub(_rmean_c, axis=0).div(_rstd_c, axis=0)
                # Apply same sample column order (tumor first, then normal)
                _samps_c = _samps if "_samps" in dir() and set(_samps).issubset(set(_z_c.columns)) else _z_c.columns.tolist()
                _z_c = _z_c[_samps_c]
                # Hierarchical clustering on rows (ward linkage, euclidean distance)
                from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram as _scipy_dendro
                _Z_mat = _z_c.fillna(0).values
                _link  = linkage(_Z_mat, method="ward", metric="euclidean")
                _order = leaves_list(_link).tolist()
                _cids_ord = [_cids_all[i] for i in _order]
                # Compute dendrogram segment coordinates for JS rendering
                _dend = _scipy_dendro(_link, no_plot=True)
                _max_dist = float(max(max(d) for d in _dend['dcoord']))
                # Drop all segments whose merge distance >= 80% of max_dist.
                # Near-root U-shapes span huge y-gaps and create detached "stray" arms.
                _clip_d = _max_dist * 0.80
                _pairs = [(ic, dc) for ic, dc in zip(_dend['icoord'], _dend['dcoord'])
                          if max(dc) < _clip_d]
                _dendro_icoord = [p[0] for p in _pairs]
                _dendro_dcoord = [p[1] for p in _pairs]
                # Build label map
                _cb_map_c = {}
                if "circbase_id" in de.columns and "in_circbase" in de.columns:
                    _known = de[de["in_circbase"] == 1].dropna(subset=["circbase_id"])
                    _cb_map_c = dict(zip(_known["circ_id"].astype(str), _known["circbase_id"].astype(str)))
                _lfc_map_c = dict(zip(de["circ_id"].astype(str), de["log2FC"]))
                _pval_map_c = dict(zip(de["circ_id"].astype(str), de[p_col]))
                _rows_c = {}
                for _cid in _cids_ord:
                    _cb = _cb_map_c.get(_cid, "")
                    _lbl = _cb if _cb and _cb not in ("", "novel") else _cid
                    _rows_c[_cid] = {
                        "z":      [round(float(v), 3) for v in _z_c.loc[_cid].tolist()],
                        "lfc":    round(float(_lfc_map_c.get(_cid, 0.0)), 2),
                        "pval":   float(_pval_map_c.get(_cid, 1.0)),
                        "label":  _lbl,
                    }
                clust_heatmap_data_js = _json.dumps({
                    "samples":    _samps_c,
                    "conditions": _cmap_hm if "_cmap_hm" in dir() else {},
                    "order":      _cids_ord,
                    "rows":       _rows_c,
                    "n_total":    len(_cids_ord),
                    "dendro": {
                        "icoord":    _dendro_icoord,
                        "dcoord":    _dendro_dcoord,
                        "max_dist":  _max_dist,
                    },
                }, ensure_ascii=False)
        except Exception as _clust_exc:
            import sys as _sys
            print(f"[report] CLUST_HEATMAP_DATA build failed: {_clust_exc}", file=_sys.stderr)

    # Compute heatmap top IDs for volcano annotation
    heatmap_ids: set = set()
    if p_col in de.columns and "log2FC" in de.columns and "circ_id" in de.columns:
        heatmap_ids = set(
            de[de["log2FC"] > 0].dropna(subset=[p_col])
            .sort_values(p_col).head(heatmap_top_n)["circ_id"].tolist()
        ) | set(
            de[de["log2FC"] < 0].dropna(subset=[p_col])
            .sort_values(p_col).head(heatmap_top_n)["circ_id"].tolist()
        )

    # ── Build ALL_DE_METHODS data (method switcher + Venn diagram) ──────────────
    all_de_data: dict = {}
    sig_sets_venn: dict = {}
    de_lookup_venn: dict = {}   # method → enriched DataFrame for Venn clickable detail
    for _mkey, _mfile in (de_files or {}).items():
        if not _mfile or not Path(_mfile).exists():
            continue
        try:
            _m_de = pd.read_csv(_mfile, sep="\t")
            if "log2FC" not in _m_de.columns and "log2FoldChange" in _m_de.columns:
                _m_de = _m_de.rename(columns={"log2FoldChange": "log2FC"})
            _m_de = _enrich_de(_m_de, isoform_file, circbase_file)
            _m_pcol, _m_thr, _m_lbl = _eff_sig(_m_de, de_sig_by, fdr)
            all_de_data[_mkey] = _build_method_js_data(
                _m_de, matrix, _m_pcol, _m_thr, lfc, groups_file, normal_label
            )
            sig_sets_venn[_mkey] = _sig_ids_from_de(_m_de, _m_pcol, _m_thr, lfc)
            de_lookup_venn[_mkey] = _m_de
        except Exception as _exc:
            import sys as _sys
            print(f"[report] DE file for {_mkey} failed: {_exc}", file=_sys.stderr)

    # Ensure primary method is in venn sets
    if de_method not in sig_sets_venn:
        sig_sets_venn[de_method] = _sig_ids_from_de(de, p_col, sig_thr, lfc)
    if de_method not in de_lookup_venn:
        de_lookup_venn[de_method] = de

    # ── Build bm_lookup for per-method score distribution ──────────────────────
    _bm_lookup: dict = {}
    if biomarker_file and Path(biomarker_file).exists():
        try:
            _bm_df = pd.read_csv(biomarker_file, sep="\t")
            _conf_vals = pd.to_numeric(_bm_df.get("confidence_score", pd.Series()), errors="coerce").fillna(0)
            _mirna_vals = pd.to_numeric(_bm_df.get("n_mirna", pd.Series(0)), errors="coerce").fillna(0)
            _rbp_vals   = pd.to_numeric(_bm_df.get("n_rbp",   pd.Series(0)), errors="coerce").fillna(0)
            _conf_mn, _conf_mx = _conf_vals.min(), _conf_vals.max()
            _mirna_mx = _mirna_vals.max(); _rbp_mx = _rbp_vals.max()
            for _, _br in _bm_df.iterrows():
                _cid = str(_br.get("circ_id", ""))
                if not _cid: continue
                _c = float(pd.to_numeric(_br.get("confidence_score", 0), errors="coerce") or 0)
                _type_raw = _br.get("Type", "")
                _type_edger = "" if (
                    _type_raw is None or _type_raw != _type_raw or
                    str(_type_raw).lower() in ("nan", "none", "na", "")
                ) else str(_type_raw)
                _bm_lookup[_cid] = {
                    "conf_n":  (_c - _conf_mn) / (_conf_mx - _conf_mn + 1e-10),
                    "known":   float(int(_br.get("in_circbase", 0) or 0)),
                    "mirna_n": float(_br.get("n_mirna", 0) or 0) / (_mirna_mx + 1e-10),
                    "rbp_n":   float(_br.get("n_rbp",   0) or 0) / (_rbp_mx   + 1e-10),
                    # Raw display values for bm_table
                    "n_mirna":       int(float(_br.get("n_mirna",       0) or 0)),
                    "n_rbp":         int(float(_br.get("n_rbp",         0) or 0)),
                    "in_circbase":   int(float(_br.get("in_circbase",   0) or 0)),
                    "circbase_id":   str(_br.get("circbase_id",   "") or ""),
                    "circbase_gene": str(_br.get("circbase_gene", "") or ""),
                    "type_edger":    _type_edger,
                }
        except Exception:
            pass

    # Override n_mirna / n_rbp from current interactions.json (may be newer than biomarker_candidates.tsv)
    # rank_biomarkers.py counts only entries where in_circ=True (or CircInteractome source).
    # After predict_interactions re-run with --gtf, more ENCORI entries have in_circ=True.
    def _is_in_circ_entry(e: dict) -> bool:
        return e.get("source") == "CircInteractome" or bool(e.get("in_circ", False))
    if interactions:
        for _oc, _olu in _bm_lookup.items():
            if _oc in interactions:
                _ixn = interactions[_oc]
                _nm = len(set(m.get("miRNAName", "") for m in _ixn.get("mirna", [])
                              if m.get("miRNAName") and _is_in_circ_entry(m)))
                _nr = len(set(r.get("RBPName", "") for r in _ixn.get("rbp", [])
                              if r.get("RBPName") and _is_in_circ_entry(r)))
                _olu["n_mirna"] = _nm
                _olu["n_rbp"]   = _nr
        # Recompute pre-normalized mirna_n / rbp_n in _bm_lookup after override
        _ixn_mirna_mx = max((_bm_lookup[c].get("n_mirna", 0) for c in _bm_lookup), default=1) or 1
        _ixn_rbp_mx   = max((_bm_lookup[c].get("n_rbp",   0) for c in _bm_lookup), default=1) or 1
        for _olu in _bm_lookup.values():
            _olu["mirna_n"] = _olu.get("n_mirna", 0) / _ixn_mirna_mx
            _olu["rbp_n"]   = _olu.get("n_rbp",   0) / _ixn_rbp_mx

    # Supplement _bm_lookup from circbase_annotated.tsv + interactions
    # for circRNAs significant in DESeq2/limma but not in edgeR biomarker_candidates.tsv
    _mirna_mx_lk = max((_bm_lookup[c].get("n_mirna", 0) for c in _bm_lookup), default=1) or 1
    _rbp_mx_lk   = max((_bm_lookup[c].get("n_rbp",   0) for c in _bm_lookup), default=1) or 1
    try:
        if circbase_file and Path(circbase_file).exists():
            _cb_df = pd.read_csv(circbase_file, sep="\t")
            _cb_conf_vals = pd.to_numeric(_cb_df.get("confidence_score", pd.Series()), errors="coerce").fillna(0)
            _cb_conf_mn = _cb_conf_vals.min(); _cb_conf_mx = _cb_conf_vals.max()
            for _, _cbr in _cb_df.iterrows():
                _cid2 = str(_cbr.get("circ_id", "") or "")
                if not _cid2 or _cid2 in _bm_lookup:
                    continue
                _cc = float(pd.to_numeric(_cbr.get("confidence_score", 0), errors="coerce") or 0)
                _nm = 0; _nr = 0
                if interactions and _cid2 in interactions:
                    _nm = len(interactions[_cid2].get("mirna", []))
                    _nr = len(interactions[_cid2].get("rbp",   []))
                _bm_lookup[_cid2] = {
                    "conf_n":  (_cc - _cb_conf_mn) / (_cb_conf_mx - _cb_conf_mn + 1e-10),
                    "known":   float(int(_cbr.get("in_circbase", 0) or 0)),
                    "mirna_n": _nm / (_mirna_mx_lk + 1e-10),
                    "rbp_n":   _nr / (_rbp_mx_lk   + 1e-10),
                    "n_mirna": _nm, "n_rbp": _nr,
                    "in_circbase":   int(float(_cbr.get("in_circbase",   0) or 0)),
                    "circbase_id":   str(_cbr.get("circbase_id",   "") or ""),
                    "circbase_gene": str(_cbr.get("circbase_gene", "") or ""),
                }
    except Exception:
        pass

    # Backfill gene_name into _bm_lookup from _iso_lookup so circbase_gene column
    # can fall back to GTF host gene when circBase has no gene annotation (e.g. "None").
    for _bid in _bm_lookup:
        if _bid in _iso_lookup:
            _bm_lookup[_bid].setdefault("gene_name", _iso_lookup[_bid].get("gene_name", ""))

    # Add score_dist to each method's data (primary + alternates)
    _msw_labels_sd = {"edgeR_ciriquant": "edgeR (FSJ offset)", "deseq2": "DESeq2", "limma": "limma-voom"}
    # Primary method
    _primary_sd = _compute_score_dist_data(de, p_col, sig_thr, lfc, _bm_lookup, _msw_labels_sd.get(de_method, de_method))
    if de_method not in all_de_data:
        all_de_data[de_method] = {}
    all_de_data[de_method]["score_dist"] = _primary_sd
    # Alternate methods
    for _mkey, _mfile in (de_files or {}).items():
        if _mkey not in all_de_data or not _mfile or not Path(_mfile).exists():
            continue
        try:
            _m_de2 = pd.read_csv(_mfile, sep="\t")
            if "log2FC" not in _m_de2.columns and "log2FoldChange" in _m_de2.columns:
                _m_de2 = _m_de2.rename(columns={"log2FoldChange": "log2FC"})
            _m_pcol2, _m_thr2, _ = _eff_sig(_m_de2, de_sig_by, fdr)
            all_de_data[_mkey]["score_dist"] = _compute_score_dist_data(
                _m_de2, _m_pcol2, _m_thr2, lfc, _bm_lookup, _msw_labels_sd.get(_mkey, _mkey))
        except Exception:
            pass

    # Compute per-method bm_table (top-30 biomarker ranking) for dynamic table update
    # Primary method (de already enriched)
    if de_method not in all_de_data:
        all_de_data[de_method] = {}
    all_de_data[de_method]["bm_table"] = _compute_bm_table_data(
        de, p_col, sig_thr, lfc, _bm_lookup, sig_sets_venn
    )
    # Alternate methods (use de_lookup_venn which stores enriched DFs from the first loop)
    for _mkey2, _m_enr_de in de_lookup_venn.items():
        if _mkey2 == de_method or _mkey2 not in all_de_data:
            continue
        try:
            _m_pcol3, _m_thr3, _ = _eff_sig(_m_enr_de, de_sig_by, fdr)
            all_de_data[_mkey2]["bm_table"] = _compute_bm_table_data(
                _m_enr_de, _m_pcol3, _m_thr3, lfc, _bm_lookup, sig_sets_venn
            )
        except Exception:
            pass

    all_de_methods_js = _json.dumps(all_de_data, ensure_ascii=False)
    venn_html = _venn_3_svg(sig_sets_venn, lfc, de_lookup_venn) if len(sig_sets_venn) >= 2 else ""

    # Build method switcher HTML
    _msw_methods = [de_method] + [m for m in all_de_data if m != de_method]
    _msw_labels  = {"edgeR_ciriquant": "edgeR (FSJ offset)", "deseq2": "DESeq2", "limma": "limma-voom"}
    _msw_btns    = [
        f'<button class="msw-btn{" active" if m == de_method else ""}" '
        f'onclick="switchDEMethod(\'{m}\')">{_msw_labels.get(m, m)}</button>'
        for m in _msw_methods if m in all_de_data or m == de_method
    ]
    _msw_html = (
        '<div style="margin:10px 0 8px;padding:8px 14px;background:#f8f9fa;'
        'border:1px solid #e0e8f0;border-radius:6px;display:flex;align-items:center;gap:6px;flex-wrap:wrap">'
        '<span style="font-size:13px;color:#555;margin-right:6px" data-en="DE Method:">分析方法：</span>'
        + "".join(_msw_btns)
        + '<span style="font-size:11px;color:#999;margin-left:8px" data-en="Switching updates Volcano / Heatmap / stats instantly">切換後 Volcano / Heatmap / 統計即時更新</span>'
        + '</div>'
    ) if len(_msw_btns) > 1 else ""

    results_dir    = str(Path(output_file).parent)
    sample_html    = _sample_overview_section(
        groups_file, tumor_label, normal_label,
        results_dir=results_dir, study_title=study_title)
    type_html      = _type_section(sig)
    biomarker_html = _biomarker_section(biomarker_file, interactions=interactions, iso_lookup=_iso_lookup)
    isoform_html   = _isoform_section(switching_file,
                                       isoform_file=isoform_file,
                                       circbase_file=circbase_file,
                                       case_label=tumor_label,
                                       control_label=normal_label)

    # Interactive Plotly charts; fall back to static PDF embeds when unavailable
    p_volcano = _plotly_volcano(de, fdr, lfc, de_method, p_col=p_col, sig_thr=sig_thr,
                                p_label=_sig_label_str, heatmap_ids=heatmap_ids)
    p_pca     = _plotly_pca(matrix, groups_file)
    volcano_html = p_volcano if p_volcano else _embed_pdf(volcano_pdf)
    pca_html     = p_pca     if p_pca     else _embed_pdf(pca_pdf)

    n_ixn = len(interactions)
    n_ixn_mirna = sum(len(v.get("mirna", [])) > 0 for v in interactions.values())

    _modal_js = f"""
<script>
const CIRC_DATA          = {interactions_js};
const VOLCANO_DATA       = {volcano_data_js};
const FULL_HEATMAP_DATA  = {full_heatmap_data_js};
const CLUST_HEATMAP_DATA = {clust_heatmap_data_js};
const ALL_DE_METHODS     = {all_de_methods_js};
const _FDR = {_fdr_js};
let _HEATMAP_DATA_CACHE  = FULL_HEATMAP_DATA;
let _CURRENT_VOLCANO_DATA = VOLCANO_DATA;
const _LFC = {_lfc_js};
const _SIG_LABEL = "{_sig_label_js}";

function showCircDetail(circId) {{
  try {{
    const d    = CIRC_DATA[circId] || {{}};
    const info = d.info || {{}};
    document.getElementById('cm-title').textContent = circId;
    document.getElementById('cm-sub').innerHTML =
      [
        (info.gene_name && info.gene_name !== 'intergenic') ? info.gene_name : null,
        (info.strand && info.strand !== '.') ? info.strand : null,
        info.region || null,
        (info.exon_span && info.exon_span !== 'nan') ? info.exon_span : null
      ].filter(Boolean).join(' &nbsp;|&nbsp; ');
    const circEl = document.getElementById('cm-circle-wrap');
    circEl.innerHTML = '';
    _drawCircleRNA(circId, circEl);
    const _mirnaEl = document.getElementById('cm-mirna');
    _mirnaEl.innerHTML = _buildInteractionTable(
      d.mirna||[], ['_priority','miRNAName','siteType','circ_pos','_seq_logo','clipExpNum','cellType','source'],
      ['Priority','miRNA','Site Type','Chr Position','Binding Seq','CLIP Exp.','Cell Type','Source'], circId, 'mirna');
    _applyLangToContainer(_mirnaEl);
    const _rbpEl = document.getElementById('cm-rbp');
    _rbpEl.innerHTML = _buildInteractionTable(
      d.rbp||[], ['_priority','RBPName','bindingSites','_mapped','circ_pos','_seq_logo','location','clipExpNum','cellType','source'],
      ['Priority','RBP','Sites','Mapped','Site Positions (hg19)','Binding Seq','Location','CLIP Exp.','Cell Type','Source'], circId, 'rbp');
    _applyLangToContainer(_rbpEl);
    _fetchSeqsInTable('cm-mirna');
    _fetchSeqsInTable('cm-rbp');
    const vEl = document.getElementById('cm-volcano');
    vEl.dataset.circId = circId;
    vEl.innerHTML = '<p style="color:#aaa;font-size:12px;padding:8px">Click tab to load volcano.</p>';
    vEl._plotlyLoaded = false;
    const hEl = document.getElementById('cm-heatmap');
    hEl.dataset.circId = circId;
    hEl.innerHTML = '<p style="color:#aaa;font-size:12px;padding:8px">Click tab to load heatmap.</p>';
    hEl._plotlyLoaded = false;
    document.getElementById('cm-tab-heatmap').style.display = '';
    _switchTab('exon');
    _updateDlBar('exon');
    document.getElementById('circ-modal').style.display = 'flex';
  }} catch(err) {{
    console.error('[showCircDetail]', err);
    document.getElementById('cm-title').textContent = circId;
    document.getElementById('cm-circle-wrap').innerHTML =
      '<p style="color:#c0392b;padding:12px"><b>JS Error:</b> '+err.message+'<br>'
      +'<small style="color:#888">Check browser console (F12) for details.</small></p>';
    document.getElementById('circ-modal').style.display = 'flex';
  }}
}}

function closeCircModal() {{
  document.getElementById('circ-modal').style.display = 'none';
  ['cm-volcano','cm-heatmap'].forEach(id => {{
    const el = document.getElementById(id);
    if (el && el._plotlyLoaded) {{ Plotly.purge(el); el._plotlyLoaded = false; }}
  }});
}}

function _switchTab(name) {{
  document.querySelectorAll('.ctab-content').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.ctab-btn').forEach(b => b.classList.remove('active'));
  document.getElementById('cm-' + name).classList.add('active');
  document.querySelector('.ctab-btn[data-tab="' + name + '"]').classList.add('active');
  if (name === 'volcano') {{
    const vEl = document.getElementById('cm-volcano');
    if (!vEl._plotlyLoaded) _buildMiniVolcano(vEl.dataset.circId);
  }}
  if (name === 'heatmap') {{
    const hEl = document.getElementById('cm-heatmap');
    if (!hEl._plotlyLoaded) _buildMiniHeatmap(hEl.dataset.circId);
  }}
  _updateDlBar(name);
}}

function _updateDlBar(name) {{
  const bar    = document.getElementById('cm-dl-bar');
  const circId = document.getElementById('cm-title').textContent.trim();
  const safe   = circId.replace(/[:|]/g,'_');
  const btn = (label, fn) =>
    `<button class="dl-btn" onclick="${{fn}}">${{label}}</button>`;
  const map = {{
    exon:    btn('⬇ SVG', `_dlSVG('${{safe}}')`),
    mirna:   btn('⬇ CSV', `_dlTabCSV('cm-mirna','${{safe}}_mirna.csv')`),
    rbp:     btn('⬇ CSV', `_dlTabCSV('cm-rbp','${{safe}}_rbp.csv')`),
    volcano: btn('⬇ PNG', `_dlPlotly('cm-volcano','${{safe}}_volcano')`),
    heatmap: btn('⬇ PNG', `_dlPlotly('cm-heatmap','${{safe}}_heatmap')`),
  }};
  bar.innerHTML = map[name] || '';
}}

function _dlSVG(safeName) {{
  const svg = document.querySelector('#cm-circle-wrap svg');
  if (!svg) return;
  const blob = new Blob([svg.outerHTML], {{type:'image/svg+xml;charset=utf-8'}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = safeName + '_circle.svg';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

function _dlTabCSV(containerId, fname) {{
  const tbl = document.querySelector('#' + containerId + ' table');
  if (!tbl) return;
  const rows = tbl.querySelectorAll('tr');
  const lines = [];
  rows.forEach(function(r) {{
    const cells = r.querySelectorAll('th,td');
    const cols = [];
    cells.forEach(function(c) {{ cols.push('"' + c.innerText.replace(/"/g,'""') + '"'); }});
    lines.push(cols.join(','));
  }});
  const blob = new Blob([lines.join('\\n')], {{type:'text/csv;charset=utf-8;'}});
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = fname;
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
}}

function _dlPlotly(elId, fname) {{
  const el = document.getElementById(elId);
  if (!el || !el._plotlyLoaded || typeof Plotly === 'undefined') return;
  Plotly.downloadImage(el, {{format:'png', filename:fname, width:900, height:550}});
}}

// ── Circular RNA diagram ──────────────────────────────────────────────────────
function _drawCircleRNA(circId, container) {{
  const d        = CIRC_DATA[circId] || {{}};
  const info     = d.info || {{}};
  const mirnaList= d.mirna || [];
  const rbpList  = d.rbp   || [];
  const exonBds  = info.exon_boundaries || [];
  const totalLen = parseInt(info.spliced_length) ||
                   (exonBds.length > 0 ? (exonBds[exonBds.length-1].cum_end || 0) : 0);

  const W=480, H=440, cx=240, cy=210;
  const ROUT=148, RIN=115, MI_OUT=176, MI_IN=153, RBP_OUT=113, RBP_IN=90;

  const MI_COLORS  = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628',
                      '#c4007a','#0077b6','#c96a00','#2d6a4f','#6d4c91','#b5470b',
                      '#1d7a3a','#7b2d8b','#0059a6','#8b6900','#c1121f','#3d7a00',
                      '#1d91c0','#6a3d9a'];
  const RBP_COLORS = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e',
                      '#e6ab02','#a6761d','#333333','#1f78b4','#b2df8a','#fb8072','#80b1d3'];
  // Distinct dash patterns for overlapping arc boundaries (indexed by badge num)
  const BOUND_DASH=['5,3','2,2','7,2,2,2','4,2,1,2','1,2','9,3','3,1,1,1,1,1'];
  const SITE_LABELS='abcdefghijklmnopqrstuvwxyz'.split('');

  function posToAngle(pos) {{
    if (!totalLen) return -Math.PI/2;
    return -Math.PI/2 + (pos/totalLen)*2*Math.PI;
  }}
  function polar(r,a) {{ return [cx+r*Math.cos(a), cy+r*Math.sin(a)]; }}
  function arcPath(ro,ri,a1,a2) {{
    const span=((a2-a1)%(2*Math.PI)+2*Math.PI)%(2*Math.PI);
    const large=span>Math.PI?1:0, f=n=>n.toFixed(2);
    const [x1,y1]=polar(ro,a1),[x2,y2]=polar(ro,a2);
    const [x3,y3]=polar(ri,a2),[x4,y4]=polar(ri,a1);
    return `M${{f(x1)}},${{f(y1)}} A${{ro}},${{ro}} 0 ${{large}},1 ${{f(x2)}},${{f(y2)}} L${{f(x3)}},${{f(y3)}} A${{ri}},${{ri}} 0 ${{large}},0 ${{f(x4)}},${{f(y4)}}Z`;
  }}

  let svg=`<svg width="${{W}}" height="${{H}}" viewBox="0 0 ${{W}} ${{H}}" style="font-family:sans-serif">`;
  svg+=`<circle cx="${{cx}}" cy="${{cy}}" r="${{(ROUT+RIN)/2}}" fill="none" stroke="#fce8e8" stroke-width="${{ROUT-RIN}}"/>`;
  if(totalLen>0&&mirnaList.length>0)
    svg+=`<circle cx="${{cx}}" cy="${{cy}}" r="${{(MI_OUT+MI_IN)/2}}" fill="none" stroke="#f0f0f0" stroke-width="${{MI_OUT-MI_IN}}"/>`;
  if(totalLen>0&&rbpList.filter(r=>r.location==='internal').length>0)
    svg+=`<circle cx="${{cx}}" cy="${{cy}}" r="${{(RBP_OUT+RBP_IN)/2}}" fill="none" stroke="#e8f0e8" stroke-width="${{RBP_OUT-RBP_IN}}"/>`;

  // ── Proportional exon boundaries (divider lines only here; labels drawn last) ──
  let exonLabelSvg='';
  if(exonBds.length>0&&totalLen>0){{
    exonBds.forEach(eb=>{{
      // divider line at exon junction
      if(eb.cum_start>0){{
        const a=posToAngle(eb.cum_start);
        const [lx1,ly1]=polar(RIN-2,a),[lx2,ly2]=polar(ROUT+2,a);
        svg+=`<line x1="${{lx1.toFixed(1)}}" y1="${{ly1.toFixed(1)}}" x2="${{lx2.toFixed(1)}}" y2="${{ly2.toFixed(1)}}" stroke="#999" stroke-width="1.5"/>`;
      }}
      // exon label — collected and drawn after arcs so labels are on top
      const arcFrac=(eb.cum_end-eb.cum_start)/totalLen;
      const arcPx=arcFrac*2*Math.PI*(RIN+ROUT)/2;
      if(arcPx>12){{  // only label if enough space (12px ≈ 1.5% of ring for small exons)
        const aMid=posToAngle((eb.cum_start+eb.cum_end)/2);
        const [lx,ly]=polar((RIN+ROUT)/2,aMid);
        const isLow=(aMid>0&&aMid<Math.PI);
        const tdeg=isLow?(aMid*180/Math.PI-90):(aMid*180/Math.PI+90);
        const exonLabel=eb.label.replace(/^e(\\d+)$/,'exon $1');
        const fsize=arcPx>25?11:9;
        exonLabelSvg+=`<text transform="translate(${{lx.toFixed(1)}},${{ly.toFixed(1)}}) rotate(${{tdeg.toFixed(1)}})" text-anchor="middle" dominant-baseline="central" font-size="${{fsize}}" fill="#444" font-weight="600">${{exonLabel}}</text>`;
      }}
    }});
  }}

  // helper: normalize angle to [0, 2π)
  function normA(a){{ return((a%(2*Math.PI))+2*Math.PI)%(2*Math.PI); }}
  // helper: de-overlap badge angles — forward bump, min gap in radians
  function deOverlap(list, gap){{
    list.sort((a,b)=>a.ang-b.ang);
    for(let i=1;i<list.length;i++){{
      if(list[i].ang-list[i-1].ang<gap) list[i].ang=list[i-1].ang+gap;
    }}
    return list;
  }}

  // ── miRNA arcs (outer ring) — group arcs per name, then de-overlapped badges ──
  window._circ_miData={{}};window._circ_rbpData={{}};  // reset on each modal open
  const miMap={{}};let miIdx=0;
  const miLegend=[];
  const miBadgeAngs={{}};
  const miArcGroups={{}};
  const miArcSegsFlat=[];
  if(totalLen>0){{
    mirnaList.forEach(item=>{{
      const name=item.miRNAName||'';
      const rawCp=String(item.circ_pos||'');
      // Skip absolute chromosome coords (ENCORI liftover failed) — cannot draw arc
      if(/^chr/.test(rawCp))return;
      const m=rawCp.match(/(\\d+)[–-](\\d+)/);
      if(!m)return;
      if(!(name in miMap)){{
        miIdx++;
        miMap[name]={{num:miIdx,color:MI_COLORS[(miIdx-1)%MI_COLORS.length]}};
        miLegend.push({{num:miIdx,name,color:miMap[name].color,st:item.siteType||'',src:item.source||''}});
        miArcGroups[name]=[];
      }}
      const e=miMap[name];
      const a1=posToAngle(parseInt(m[1])),a2r=posToAngle(Math.max(parseInt(m[2]),parseInt(m[1])+1));
      const minA=2*Math.PI/totalLen*3, a2=Math.max(a2r,a1+minA);
      const _mfi=miArcSegsFlat.length;
      miArcSegsFlat.push({{name,num:e.num,color:e.color,a1,a2}});
      miArcGroups[name].push({{d:arcPath(MI_OUT,MI_IN,a1,a2),title:`${{name}} · ${{item.siteType||''}} · ${{item.circ_pos}}`,a1,a2,fi:_mfi}});
      if(!(name in miBadgeAngs)) miBadgeAngs[name]=normA((a1+a2)/2);
    }});
    // Compute overlaps; assign site letters sorted by angle
    const miOvlSet=new Set();
    for(let i=0;i<miArcSegsFlat.length;i++){{
      for(let j=i+1;j<miArcSegsFlat.length;j++){{
        const s=miArcSegsFlat[i],t=miArcSegsFlat[j];
        if(s.a1<t.a2&&s.a2>t.a1){{miOvlSet.add(i);miOvlSet.add(j);}}
      }}
    }}
    const miOvlSorted=[...miOvlSet].sort((a,b)=>miArcSegsFlat[a].a1-miArcSegsFlat[b].a1);
    miOvlSorted.forEach((fi,li)=>{{miArcSegsFlat[fi].label=SITE_LABELS[li%26];}});
    // Build per-molecule arc data for dynamic right-panel site toggle
    const _miArcData={{}};
    Object.entries(miArcGroups).forEach(([name,arcs])=>{{
      const e=miMap[name];if(!e)return;
      _miArcData[e.num]={{name,color:e.color,arcs:arcs.map((a,ai)=>
        ({{letter:SITE_LABELS[ai%26],id:`_mi_site_${{e.num}}_${{ai}}`,pos:a.title||''}}))
      }};
    }});
    // Draw arc groups (each name wrapped in <g id="_miarc_N"> for L1 show/hide;
    //   each arc wrapped in <g id="_mi_site_N_i"> for L2 per-site toggle)
    Object.entries(miArcGroups).forEach(([name,arcs])=>{{
      const e=miMap[name];
      svg+=`<g id="_miarc_${{e.num}}" class="_miarc">`;
      arcs.forEach((a,ai)=>{{
        const arcLetter=SITE_LABELS[ai%26];
        svg+=`<g id="_mi_site_${{e.num}}_${{ai}}" class="_mi_site">`;
        svg+=`<path d="${{a.d}}" fill="${{e.color}}" opacity="0.85"><title>${{a.title}}</title></path>`;
        const fi=(a.fi!==undefined)?a.fi:-1;
        if(fi>=0&&miOvlSet.has(fi)){{
          const seg=miArcSegsFlat[fi];
          const dp=BOUND_DASH[(seg.num-1)%BOUND_DASH.length];
          if(seg.label)svg+=`<g id="_ann_mi_${{seg.label}}" class="_miann">`;
          [a.a1,a.a2].forEach(ang=>{{
            const [lx1,ly1]=polar(MI_IN-1,ang),[lx2,ly2]=polar(MI_OUT+1,ang);
            svg+=`<line x1="${{lx1.toFixed(1)}}" y1="${{ly1.toFixed(1)}}" x2="${{lx2.toFixed(1)}}" y2="${{ly2.toFixed(1)}}" stroke="${{e.color}}" stroke-width="2.5" stroke-dasharray="${{dp}}" opacity="0.95"/>`;
          }});
          // overlap letter badge removed — per-arc site labels (_arc_lbl) used instead
          if(seg.label)svg+=`</g>`;
        }}
        // hidden per-arc letter badge — only made visible when this molecule is the sole visible one
        if(Math.abs(a.a2-a.a1)*(MI_IN+MI_OUT)/2>8){{
          const midA=(a.a1+a.a2)/2,midR=(MI_IN+MI_OUT)/2;
          const [lx,ly]=polar(midR,midA);
          svg+=`<g class="_arc_lbl" data-arc-type="mi" data-arc-num="${{e.num}}" style="visibility:hidden">`;
          svg+=`<circle cx="${{lx.toFixed(1)}}" cy="${{ly.toFixed(1)}}" r="6" fill="${{e.color}}" opacity="0.93" stroke="white" stroke-width="1.2"/>`;
          svg+=`<text x="${{lx.toFixed(1)}}" y="${{ly.toFixed(1)}}" text-anchor="middle" dominant-baseline="central" font-size="9" fill="white" font-weight="bold">${{arcLetter}}</text>`;
          svg+=`</g>`;
        }}
        svg+=`</g>`; // end _mi_site
      }});
      svg+=`</g>`;
    }});
    // place badges immediately outside arc; radially stagger only when angularly too close
    const miBadgesSorted=Object.entries(miBadgeAngs)
      .map(([n,a])=>{{return{{name:n,ang:a,e:miMap[n]}}}})
      .filter(b=>b.e);
    miBadgesSorted.sort((a,b)=>a.ang-b.ang);
    const MI_R0=MI_OUT+7, MI_R1=MI_OUT+19, MI_ANG_MIN=0.14;
    let miLastAng=-Infinity, miRIdx=0;
    miBadgesSorted.forEach(b=>{{
      if(b.ang-miLastAng<MI_ANG_MIN){{ miRIdx=1-miRIdx; }}
      else{{ miRIdx=0; }}
      b.r=miRIdx===0?MI_R0:MI_R1;
      miLastAng=b.ang;
    }});
    miBadgesSorted.forEach(b=>{{
      const [nx,ny]=polar(b.r,b.ang);
      svg+=`<g id="_mib_${{b.e.num}}" class="_mibadge">`;
      if(b.r>MI_R0){{
        const [ox,oy]=polar(MI_OUT+2,b.ang);
        svg+=`<line x1="${{ox.toFixed(1)}}" y1="${{oy.toFixed(1)}}" x2="${{nx.toFixed(1)}}" y2="${{ny.toFixed(1)}}" stroke="${{b.e.color}}" stroke-width="0.8" opacity="0.5" stroke-dasharray="2,2"/>`;
      }}
      svg+=`<circle cx="${{nx.toFixed(1)}}" cy="${{ny.toFixed(1)}}" r="6" fill="${{b.e.color}}"/>`;
      svg+=`<text x="${{nx.toFixed(1)}}" y="${{ny.toFixed(1)}}" text-anchor="middle" dominant-baseline="central" font-size="7.5" fill="white" font-weight="bold">${{b.e.num}}</text>`;
      svg+=`</g>`;
    }});
    window._circ_miData=_miArcData;  // store before block closes
  }}

  // ── RBP arcs (inner ring) — draw arcs first, then de-overlapped badges ──
  const rbpMap={{}};let rbpIdx=0;
  const rbpLegend=[];
  const rbpBadgeAngs={{}};
  const rbpArcGroups={{}};
  const rbpArcSegsFlat=[];
  if(totalLen>0){{
    rbpList.forEach(item=>{{
      if(item.location!=='internal')return;
      const name=item.RBPName||'';
      if(!(name in rbpMap)){{
        rbpIdx++;
        rbpMap[name]={{num:rbpIdx,color:RBP_COLORS[(rbpIdx-1)%RBP_COLORS.length]}};
        const rbpTotal=parseInt(item.bindingSites||0)||(item.sites||[]).length;
        const rbpMapped=(item.sites||[]).length;
        rbpLegend.push({{num:rbpIdx,name,color:rbpMap[name].color,ns:rbpTotal,nm:rbpMapped,src:item.source||''}});
        rbpArcGroups[name]=[];
      }}
      const e=rbpMap[name];
      const sites=item.sites||[];
      if(sites.length>0){{
        sites.forEach(s=>{{
          const s0=Math.max(0,s.circ_start);
          const s1=Math.min(Math.max(s.circ_end,s0+1),totalLen); // clamp to [0, totalLen]
          const a1=posToAngle(s0),a2=posToAngle(s1);
          const _rfi=rbpArcSegsFlat.length;
          rbpArcSegsFlat.push({{name,num:e.num,color:e.color,a1,a2}});
          rbpArcGroups[name].push({{d:arcPath(RBP_OUT,RBP_IN,a1,a2),title:`${{name}} · ${{s.circ_pos}}`,a1,a2,fi:_rfi}});
          if(!(name in rbpBadgeAngs)) rbpBadgeAngs[name]=normA((a1+a2)/2);
        }});
      }} else {{
        if(!(name in rbpBadgeAngs)){{
          const fallbackAng=normA(-Math.PI/2+(rbpIdx-1)/(Math.max(rbpList.filter(r=>r.location==='internal').length,1))*2*Math.PI);
          rbpBadgeAngs[name]=fallbackAng;
        }}
      }}
    }});
    // Compute overlaps; assign site letters sorted by angle
    const rbpOvlSet=new Set();
    for(let i=0;i<rbpArcSegsFlat.length;i++){{
      for(let j=i+1;j<rbpArcSegsFlat.length;j++){{
        const s=rbpArcSegsFlat[i],t=rbpArcSegsFlat[j];
        if(s.a1<t.a2&&s.a2>t.a1){{rbpOvlSet.add(i);rbpOvlSet.add(j);}}
      }}
    }}
    const rbpOvlSorted=[...rbpOvlSet].sort((a,b)=>rbpArcSegsFlat[a].a1-rbpArcSegsFlat[b].a1);
    rbpOvlSorted.forEach((fi,li)=>{{rbpArcSegsFlat[fi].label=SITE_LABELS[li%26];}});
    // Build per-molecule arc data for dynamic right-panel site toggle (RBP)
    const _rbpArcData={{}};
    Object.entries(rbpArcGroups).forEach(([name,arcs])=>{{
      const e=rbpMap[name];if(!e)return;
      _rbpArcData[e.num]={{name,color:e.color,arcs:arcs.map((a,ai)=>
        ({{letter:SITE_LABELS[ai%26],id:`_rbp_site_${{e.num}}_${{ai}}`,pos:a.title||''}}))
      }};
    }});
    // Draw RBP arc groups (L1: _rbparc_N; L2: _rbp_site_N_i per arc)
    Object.entries(rbpArcGroups).forEach(([name,arcs])=>{{
      const e=rbpMap[name];
      svg+=`<g id="_rbparc_${{e.num}}" class="_rbparc">`;
      arcs.forEach((a,ai)=>{{
        const arcLetter=SITE_LABELS[ai%26];
        svg+=`<g id="_rbp_site_${{e.num}}_${{ai}}" class="_rbp_site">`;
        svg+=`<path d="${{a.d}}" fill="${{e.color}}" opacity="0.9" stroke="rgba(0,0,0,0.3)" stroke-width="0.8"><title>${{a.title}}</title></path>`;
        const fi=(a.fi!==undefined)?a.fi:-1;
        if(fi>=0&&rbpOvlSet.has(fi)){{
          const seg=rbpArcSegsFlat[fi];
          const dp=BOUND_DASH[(seg.num-1)%BOUND_DASH.length];
          if(seg.label)svg+=`<g id="_ann_rbp_${{seg.label}}" class="_rbpann">`;
          [a.a1,a.a2].forEach(ang=>{{
            const [lx1,ly1]=polar(RBP_IN-1,ang),[lx2,ly2]=polar(RBP_OUT+1,ang);
            svg+=`<line x1="${{lx1.toFixed(1)}}" y1="${{ly1.toFixed(1)}}" x2="${{lx2.toFixed(1)}}" y2="${{ly2.toFixed(1)}}" stroke="${{e.color}}" stroke-width="2.5" stroke-dasharray="${{dp}}" opacity="0.95"/>`;
          }});
          // overlap letter badge removed — per-arc site labels (_arc_lbl) used instead
          if(seg.label)svg+=`</g>`;
        }}
        // hidden per-arc letter badge — placed on exon ring (midR between RIN and ROUT)
        if(Math.abs(a.a2-a.a1)*(RBP_IN+RBP_OUT)/2>8){{
          const midA=(a.a1+a.a2)/2,midR=(RIN+ROUT)/2;
          const [lx,ly]=polar(midR,midA);
          svg+=`<g class="_arc_lbl" data-arc-type="rbp" data-arc-num="${{e.num}}" style="visibility:hidden">`;
          svg+=`<circle cx="${{lx.toFixed(1)}}" cy="${{ly.toFixed(1)}}" r="6" fill="${{e.color}}" opacity="0.93" stroke="white" stroke-width="1.2"/>`;
          svg+=`<text x="${{lx.toFixed(1)}}" y="${{ly.toFixed(1)}}" text-anchor="middle" dominant-baseline="central" font-size="9" fill="white" font-weight="bold">${{arcLetter}}</text>`;
          svg+=`</g>`;
        }}
        svg+=`</g>`; // end _rbp_site
      }});
      svg+=`</g>`;
    }});
    // build + de-overlap badge list, place inside RBP ring (in center-hole)
    const rbpBadgesSorted=Object.entries(rbpBadgeAngs)
      .map(([n,a])=>{{return{{name:n,ang:a,e:rbpMap[n]}}}})
      .filter(b=>b.e);
    rbpBadgesSorted.sort((a,b)=>a.ang-b.ang);
    const RBP_R0=RBP_IN-8, RBP_R1=RBP_IN-20, RBP_ANG_MIN=0.14;
    let rbpLastAng=-Infinity, rbpRIdx=0;
    rbpBadgesSorted.forEach(b=>{{
      if(b.ang-rbpLastAng<RBP_ANG_MIN){{ rbpRIdx=1-rbpRIdx; }}
      else{{ rbpRIdx=0; }}
      b.r=rbpRIdx===0?RBP_R0:RBP_R1;
      rbpLastAng=b.ang;
    }});
    rbpBadgesSorted.forEach(b=>{{
      const [nx,ny]=polar(b.r,b.ang);
      svg+=`<g id="_rbpb_${{b.e.num}}" class="_rbpbadge">`;
      if(b.r<RBP_R0){{
        const [ox,oy]=polar(RBP_IN-2,b.ang);
        svg+=`<line x1="${{ox.toFixed(1)}}" y1="${{oy.toFixed(1)}}" x2="${{nx.toFixed(1)}}" y2="${{ny.toFixed(1)}}" stroke="${{b.e.color}}" stroke-width="0.8" opacity="0.5" stroke-dasharray="2,2"/>`;
      }}
      svg+=`<circle cx="${{nx.toFixed(1)}}" cy="${{ny.toFixed(1)}}" r="6" fill="${{b.e.color}}"/>`;
      svg+=`<text x="${{nx.toFixed(1)}}" y="${{ny.toFixed(1)}}" text-anchor="middle" dominant-baseline="central" font-size="7.5" fill="white" font-weight="bold">${{b.e.num}}</text>`;
      svg+=`</g>`;
    }});
    window._circ_rbpData=_rbpArcData;  // store before block closes
  }}

  // ── Exon labels drawn last (on top of all arc layers) ──
  svg+=exonLabelSvg;

  // BSJ junction line — middle ring only, same weight as exon dividers
  {{
    const bsjA=-Math.PI/2;
    const [bl1x,bl1y]=polar(RIN-2,bsjA),[bl2x,bl2y]=polar(ROUT+2,bsjA);
    svg+=`<line x1="${{bl1x.toFixed(1)}}" y1="${{bl1y.toFixed(1)}}" x2="${{bl2x.toFixed(1)}}" y2="${{bl2y.toFixed(1)}}" stroke="#999" stroke-width="1.5"/>`;
  }}
  // BSJ marker triangle + label
  const [bx,by]=polar(ROUT+12,-Math.PI/2);
  svg+=`<polygon points="${{bx.toFixed(1)}},${{(by-10).toFixed(1)}} ${{(bx-7).toFixed(1)}},${{(by+4).toFixed(1)}} ${{(bx+7).toFixed(1)}},${{(by+4).toFixed(1)}}" fill="#d62728"/>`;
  svg+=`<text x="${{bx.toFixed(1)}}" y="${{(by-13).toFixed(1)}}" text-anchor="middle" font-size="9" fill="#d62728" font-weight="bold">BSJ</text>`;

  // Center text
  const _gn=(info.gene_name && info.gene_name!=='intergenic')?info.gene_name:(info.region||'—');
  const _st=(info.strand && info.strand!=='.')?` (${{info.strand}})`:'';
  const gname=`${{_gn}}${{_st}}`;
  const eline=(info.exon_span && info.exon_span!=='nan')?`exon ${{info.exon_span}}`:(info.region&&info.region!=='intergenic'?info.region:'');
  const lline=totalLen?`${{totalLen}} nt`:'(length unknown)';
  svg+=`<text x="${{cx}}" y="${{cy-18}}" text-anchor="middle" font-size="13" font-weight="bold" fill="#333">${{gname}}</text>`;
  svg+=`<text x="${{cx}}" y="${{cy}}"    text-anchor="middle" font-size="11" fill="#666">${{eline}}</text>`;
  svg+=`<text x="${{cx}}" y="${{cy+16}}" text-anchor="middle" font-size="10" fill="#999">${{lline}}</text>`;
  svg+='</svg>';

  // Legend as HTML below SVG — clickable chips to toggle badge visibility
  const badge=(c,n)=>`<span style="background:${{c}};color:white;border-radius:50%;width:15px;height:15px;display:inline-flex;align-items:center;justify-content:center;font-size:8px;font-weight:bold;margin-right:3px;flex-shrink:0">${{n}}</span>`;
  const srcBadge=(src)=>{{
    const col=src==='ENCORI'?'#0077b6':'#6c757d';
    return `<span style="background:${{col}};color:white;border-radius:3px;padding:0 3px;font-size:7px;margin-left:2px;vertical-align:middle">${{src}}</span>`;
  }};
  const toggleBtn=(label,onclick)=>`<button onclick="${{onclick}}" style="font-size:9px;padding:1px 5px;border:1px solid #ccc;border-radius:3px;background:#f8f8f8;cursor:pointer;margin-left:4px">${{label}}</button>`;

  let legHtml='';
  if(miLegend.length>0){{
    legHtml+=`<div style="margin-top:8px;font-size:10px">
      <span style="font-weight:bold">miRNA (outer ring):</span>
      ${{toggleBtn(_LS[_LANG||'zh'].toggleShow,'_toggleAll(\\'mi\\',true,this.closest(\\'div\\'))')}}
      ${{toggleBtn(_LS[_LANG||'zh'].toggleHide,'_toggleAll(\\'mi\\',false,this.closest(\\'div\\'))')}}
      <div style="display:flex;flex-wrap:wrap;gap:3px 8px;margin-top:4px">`;
    miLegend.forEach(e=>{{
      legHtml+=`<span id="_leg_mi_${{e.num}}" onclick="_toggleBadge('mi',${{e.num}},this)"
        title="${{_LS[_LANG||'zh'].toggleTip}}" style="display:inline-flex;align-items:center;gap:2px;cursor:pointer;padding:1px 3px;border-radius:3px;border:1px solid #eee">
        ${{badge(e.color,e.num)}}${{e.name}}${{e.src?srcBadge(e.src):''}}</span>`;
    }});
    legHtml+='</div></div>';
  }}
  if(rbpLegend.length>0){{
    legHtml+=`<div style="margin-top:6px;font-size:10px">
      <span style="font-weight:bold">RBP (inner ring):</span>
      ${{toggleBtn(_LS[_LANG||'zh'].toggleShow,'_toggleAll(\\'rbp\\',true,this.closest(\\'div\\'))')}}
      ${{toggleBtn(_LS[_LANG||'zh'].toggleHide,'_toggleAll(\\'rbp\\',false,this.closest(\\'div\\'))')}}
      <div style="display:flex;flex-wrap:wrap;gap:3px 8px;margin-top:4px">`;
    rbpLegend.forEach(e=>{{
      // nm=0 + ENCORI → fallback badge (no arc drawn); show as CLIP count not sites
      const encoriFallback=(e.nm===0&&(e.src||'').includes('ENCORI'));
      const siteStr=encoriFallback
        ?`<span style="color:#aaa" title="ENCORI 染色體絕對座標，位置未對應至 spliced 序列">${{e.ns}} CLIP</span>`
        :(e.nm>0&&e.nm<e.ns?`${{e.nm}}/${{e.ns}} sites`:`${{e.ns}} sites`);
      const borderStyle=encoriFallback?'border:1px dashed #ccc;opacity:.7':'border:1px solid #eee';
      legHtml+=`<span id="_leg_rbp_${{e.num}}" onclick="_toggleBadge('rbp',${{e.num}},this)"
        title="${{encoriFallback?'ENCORI (染色體絕對座標，badge 為均分位置)':_LS[_LANG||'zh'].toggleTip}}" style="display:inline-flex;align-items:center;gap:2px;cursor:pointer;padding:1px 3px;border-radius:3px;${{borderStyle}}">
        ${{badge(e.color,e.num)}}${{e.name}} (${{siteStr}})${{e.src?srcBadge(e.src):''}}</span>`;
    }});
    legHtml+='</div></div>';
  }}
  // Dynamic right panel (populated by _showSitePanel when exactly one molecule is visible)
  const sitePanelDiv=`<div id="circ-site-panel" style="display:none;font-size:10px;padding:6px 0 0 10px;border-left:2px solid #d8e8f4;min-width:150px;max-width:210px;align-self:flex-start;margin-top:8px"></div>`;
  container.innerHTML=`<div style="display:flex;justify-content:center;align-items:flex-start;gap:6px"><div style="text-align:center">${{svg}}</div>${{sitePanelDiv}}</div>`+legHtml;
}}

// ── Mini heatmap ──────────────────────────────────────────────────────────────
function _buildMiniHeatmap(circId) {{
  const el=document.getElementById('cm-heatmap');
  if(!el)return;
  if(typeof Plotly==='undefined'){{el.innerHTML='<p class="no-data">Plotly not available.</p>';return;}}
  const cd=CLUST_HEATMAP_DATA;
  if(!cd||!cd.order||!cd.order.length){{
    el.innerHTML='<p class="no-data">Clustering heatmap data not available.</p>';return;
  }}
  const order=cd.order, rows=cd.rows||{{}}, samps=cd.samples||[];
  const idx=order.indexOf(circId);
  if(idx<0){{
    el.innerHTML=`<p class="no-data">${{circId}} is not in the primary method's significant set (not clustered).</p>`;
    return;
  }}
  // Zoom to a window of ~20 neighboring rows in cluster order (Ward-linkage
  // neighbors = similar expression pattern), target centered, clamped/shifted
  // at array bounds so the window always has up to N rows.
  const N=20, HALF=Math.floor(N/2);
  let lo=idx-HALF, hi=idx+(N-HALF);
  if(lo<0){{hi+=(-lo);lo=0;}}
  if(hi>order.length){{lo-=(hi-order.length);hi=order.length;}}
  lo=Math.max(0,lo);
  const winIds=order.slice(lo,hi);
  const tIdx=winIds.indexOf(circId);
  const winRows={{}}; winIds.forEach(id=>{{winRows[id]=rows[id];}});
  const winDendro=_sliceDendro(cd.dendro, lo, hi);
  const titleText=`${{circId}} — clustering neighborhood (${{winIds.length}} of ${{order.length}})`;
  _renderClustPanel(el, winIds, winRows, samps, cd.conditions||{{}}, winDendro,
                     {{highlightIdx: tIdx>=0?tIdx:null, title: titleText}});
  el._plotlyLoaded=true;
}}

// ── Mini volcano ──────────────────────────────────────────────────────────────
function _buildMiniVolcano(circId) {{
  const el=document.getElementById('cm-volcano');
  if(!circId||typeof Plotly==='undefined'){{el.innerHTML='<p class="no-data">Plotly not available.</p>';return;}}
  const _vd=_CURRENT_VOLCANO_DATA||VOLCANO_DATA;
  const pt=_vd.find(d=>d[3]===circId);
  if(!pt){{el.innerHTML='<p class="no-data">circRNA not found in volcano data.</p>';return;}}
  const ns_x=[],ns_y=[],up_x=[],up_y=[],dn_x=[],dn_y=[];
  _vd.forEach(d=>{{
    if(d[2]==='N'){{ns_x.push(d[0]);ns_y.push(d[1]);}}
    else if(d[2]==='U'){{up_x.push(d[0]);up_y.push(d[1]);}}
    else{{dn_x.push(d[0]);dn_y.push(d[1]);}}
  }});
  const thr_y=-Math.log10(_FDR);
  const yLab=`−log₁₀(${{_SIG_LABEL}})`;
  Plotly.newPlot(el,[
    {{x:ns_x,y:ns_y,mode:'markers',name:'NS',marker:{{color:'rgba(150,150,150,0.25)',size:3,line:{{width:0}}}},hoverinfo:'skip'}},
    {{x:up_x,y:up_y,mode:'markers',name:'Up',marker:{{color:'rgba(214,39,40,0.55)',size:4}},hoverinfo:'skip'}},
    {{x:dn_x,y:dn_y,mode:'markers',name:'Down',marker:{{color:'rgba(44,160,44,0.55)',size:4}},hoverinfo:'skip'}},
    {{x:[pt[0]],y:[pt[1]],mode:'markers+text',name:circId,
      text:[circId.split(':').slice(-1)[0]],textposition:'top center',textfont:{{size:9}},
      marker:{{color:'#ff8c00',size:16,symbol:'star',line:{{color:'#000',width:1.5}}}},
      hovertemplate:`<b>${{circId}}</b><br>log₂FC: ${{pt[0].toFixed(3)}}<br>${{yLab}}: ${{pt[1].toFixed(3)}}<extra></extra>`}},
  ],{{
    height:320,margin:{{t:36,b:50,l:60,r:20}},
    title:{{text:'Volcano Plot Position',font:{{size:13}}}},
    xaxis:{{title:'log₂ Fold Change (Tumor / Normal)',showgrid:true,gridcolor:'#f0f0f0',zeroline:false}},
    yaxis:{{title:yLab,showgrid:true,gridcolor:'#f0f0f0',zeroline:false}},
    showlegend:false,plot_bgcolor:'white',paper_bgcolor:'white',
    shapes:[
      {{type:'line',x0:_LFC,x1:_LFC,y0:0,y1:1,yref:'paper',line:{{dash:'dot',color:'#bbb',width:1}}}},
      {{type:'line',x0:-_LFC,x1:-_LFC,y0:0,y1:1,yref:'paper',line:{{dash:'dot',color:'#bbb',width:1}}}},
      {{type:'line',x0:0,x1:1,xref:'paper',y0:thr_y,y1:thr_y,line:{{dash:'dot',color:'#bbb',width:1}}}},
    ],
  }},{{responsive:true,displayModeBar:false}});
  el._plotlyLoaded=true;
}}

function _calcPriority(r, type) {{
  let s=0;
  if(type==='mirna') {{
    const st=(r.siteType||'').toLowerCase();
    if(st.includes('8mer')) s+=4;
    else if(st.includes('7mer-m8')) s+=3;
    else if(st.includes('7mer-1a')) s+=2;
    else if(st.includes('6mer')) s+=1;
  }} else {{
    const n=parseInt(r.bindingSites||0);
    s+=Math.min(3, Math.log2(n+1));
    if((r.location||'').toLowerCase().includes('internal')) s+=2;
  }}
  if(parseInt(r.clipExpNum||0)>0) s+=3;
  if(r.source==='ENCORI') s+=2;
  const ic=r.in_circ;
  if(ic===true||ic==='true') s+=1;
  return Math.round(s*10)/10;
}}

function _sortITable(th) {{
  const tbl=th.closest('table');
  const col=th.cellIndex;
  const asc=th.dataset.asc!=='true';
  th.dataset.asc=asc;
  th.closest('tr').querySelectorAll('th').forEach(h=>{{
    h.textContent=h.textContent.replace(' ▲','').replace(' ▼','');
  }});
  th.textContent+=(asc?' ▲':' ▼');
  const tbody=tbl.querySelector('tbody');
  Array.from(tbody.querySelectorAll('tr'))
    .sort((a,b)=>{{
      const av=a.cells[col].dataset.val!==undefined?parseFloat(a.cells[col].dataset.val):a.cells[col].innerText.trim();
      const bv=b.cells[col].dataset.val!==undefined?parseFloat(b.cells[col].dataset.val):b.cells[col].innerText.trim();
      if(typeof av==='number'&&typeof bv==='number') return asc?av-bv:bv-av;
      return asc?String(av).localeCompare(String(bv)):String(bv).localeCompare(String(av));
    }})
    .forEach(r=>tbody.appendChild(r));
}}

function _buildInteractionTable(rows, keys, headers, circId, tableType) {{
  if(!rows||!rows.length) return '<p class="no-data">No data available (novel circRNA or source returned no results).</p>';

  // Drop entries that failed to liftover onto this circRNA (in_circ === false)
  // rather than showing an unusable/unverifiable row.
  rows=rows.filter(r=>{{const ic=r.in_circ; return ic===undefined||ic===true||ic==='true';}});
  if(!rows.length) return '<p class="no-data">No data available (no binding sites could be mapped onto this circRNA).</p>';

  // compute priority and pre-sort descending
  rows=rows.map(r=>Object.assign({{}},r,{{_priority:_calcPriority(r,tableType||'mirna')}}));
  rows.sort((a,b)=>b._priority-a._priority);

  // parse chr and genomic start from circId (e.g. "chr2:56813056|56820808")
  let chrom='', chromStart=0;
  if(circId) {{
    const cm=circId.match(/^(.+):(\\d+)[|](\\d+)$/);
    if(cm) {{ chrom=cm[1]; chromStart=parseInt(cm[2]); }}
  }}

  const _srcBadge=s=>{{
    if(!s)return'—';
    const col=s==='ENCORI'?'#0077b6':'#6c757d';
    return`<span style="background:${{col}};color:white;border-radius:3px;padding:1px 5px;font-size:10px">${{s}}</span>`;
  }};

  const _absPos=v=>{{
    if(!v||v==='—') return v;
    const sv=String(v);
    // ENCORI: already absolute "chrN:start-end" → pass through as-is
    if(/^chr[^\s:]+:\d/.test(sv)) return sv.replace(/[–-]/,'-');
    // CircInteractome: relative "start–end" → add chromStart
    if(!chrom) return v;
    const pm=sv.match(/(\\d+)[–-](\\d+)/);
    if(!pm) return v;
    const a=chromStart+parseInt(pm[1])-1;
    const b=chromStart+parseInt(pm[2])-1;
    return chrom+':'+a.toLocaleString()+'-'+b.toLocaleString();
  }};

  const priTitle='Priority score (click to sort):\\n'
    +(tableType==='mirna'
      ?'8mer=+4, 7mer-m8=+3, 7mer-1A=+2, 6mer=+1\\nCLIP exp>0=+3, ENCORI=+2, in_circ=+1'
      :'bindingSites log2×2 (max+3), internal=+2\\nCLIP exp>0=+3, ENCORI=+2, in_circ=+1');

  // info banner: miRNA seed types OR RBP scoring guide
  const _infoDiv='background:#f0f6ff;border:1px solid #c5d9f0;border-radius:6px;padding:8px 12px;margin-bottom:10px;font-size:12px;line-height:1.6;color:#444';
  const _th='padding:1px 8px 1px 0;font-weight:bold;white-space:nowrap;vertical-align:top';
  const infoNote = tableType==='mirna' ? `
  <div style="${{_infoDiv}}">
    <span style="font-weight:bold;color:#2c6fad" data-en="ℹ miRNA Seed Type Guide">ℹ miRNA Seed Type 說明</span>
    &nbsp;&mdash;&nbsp;<span data-en="Binding strength is determined by the seed region (miRNA nt 2–8) complementarity:">miRNA 與靶點結合的強度由 seed 區（miRNA 2–8 nt）的配對完整度決定：</span>
    <table style="margin-top:5px;border-collapse:collapse;font-size:11.5px">
      <tr>
        <td style="${{_th}};color:#1a6e3c">8mer</td>
        <td><span data-en="seed (pos 2–8) + pos 8 match + pos 1 = A → &lt;b style='color:#1a6e3c'&gt;Strongest&lt;/b&gt;, highest TargetScan confidence">seed（位置 2–8）+ 位置 8 配對 + 位置 1 為 A &nbsp;→&nbsp;<b style="color:#1a6e3c">最強</b>，TargetScan 最高可信度</span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#3a7ebf">7mer-m8</td>
        <td><span data-en="seed (pos 2–8) + pos 8 match, no pos 1 requirement → &lt;b style='color:#3a7ebf'&gt;Strong&lt;/b&gt;">seed（位置 2–8）+ 位置 8 配對，無位置 1 限制 &nbsp;→&nbsp;<b style="color:#3a7ebf">強</b></span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#e07b39">7mer-1A</td>
        <td><span data-en="seed (pos 2–7) + pos 1 = A, no pos 8 requirement → &lt;b style='color:#e07b39'&gt;Moderate&lt;/b&gt;">seed（位置 2–7）+ 位置 1 為 A，無位置 8 限制 &nbsp;→&nbsp;<b style="color:#e07b39">中等</b></span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#888">6mer</td>
        <td><span data-en="seed (pos 2–7) match only → &lt;b style='color:#888'&gt;Weak&lt;/b&gt;, higher false-positive rate">僅 seed（位置 2–7）配對，最短 seed &nbsp;→&nbsp;<b style="color:#888">弱</b>，誤報率較高</span></td>
      </tr>
    </table>
    <div style="margin-top:5px;font-size:11px;color:#888">
      <span data-en="Priority score: seed strength (8mer=+4, 7mer-m8=+3, 7mer-1A=+2, 6mer=+1) + CLIP exp&gt;0 (+3) + ENCORI (+2) + in_circ (+1)">Priority score：seed 強度（8mer=+4, 7mer-m8=+3, 7mer-1A=+2, 6mer=+1）+ CLIP 實驗支持（+3）+ ENCORI 收錄（+2）+ 位於 circRNA 內（+1）</span>
    </div>
  </div>` : tableType==='rbp' ? `
  <div style="${{_infoDiv}}">
    <span style="font-weight:bold;color:#2c6fad" data-en="ℹ RBP Priority Score Guide">ℹ RBP Priority Score 說明</span>
    &nbsp;&mdash;&nbsp;<span data-en="Score reflects the strength of RBP–circRNA binding evidence:">分數反映 RBP 與 circRNA 結合的實驗證據強度：</span>
    <table style="margin-top:5px;border-collapse:collapse;font-size:11.5px">
      <tr>
        <td style="${{_th}};color:#2c6fad"><span data-en="Binding Sites">結合位點數</span></td>
        <td><span data-en="log₂(sites+1)×2, capped at +3 (≥ 3 sites = max)">log₂(sites+1)×2，上限 +3（≥ 3 個位點即達最高分）</span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#1a6e3c">internal</td>
        <td><span data-en="Binding within exon body, not at junction boundary → +2">結合位點位於 exon 內部（非 junction 邊界）→ +2</span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#e07b39"><span data-en="CLIP Exp.">CLIP 實驗</span></td>
        <td><span data-en="CLIP experiment count &gt; 0 → +3 (binding validated in cell lines)">CLIP 實驗數量 &gt; 0 → +3（細胞株中實驗驗證的結合）</span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#2c6fad">ENCORI</td>
        <td><span data-en="Listed in ENCORI database → +2">收錄於 ENCORI 資料庫 → +2</span></td>
      </tr>
      <tr>
        <td style="${{_th}};color:#1a6e3c"><span data-en="In circRNA">In circRNA</span></td>
        <td><span data-en="Binding site falls within the circRNA span → +1">結合位點位於 circRNA 範圍內 → +1</span></td>
      </tr>
    </table>
    <div style="margin-top:5px;font-size:11px;color:#888">
      <span data-en="Priority score: log₂(sites+1)×2 (max +3) + internal (+2) + CLIP exp&gt;0 (+3) + ENCORI (+2) + in_circ (+1)">Priority score：log₂(sites+1)×2（上限 +3）+ internal（+2）+ CLIP 實驗（+3）+ ENCORI（+2）+ in circRNA（+1）</span>
    </div>
  </div>` : '';

  let html= infoNote + '<table class="itable"><thead><tr>';
  headers.forEach((h,i)=>{{
    const tip=i===0?` title="${{priTitle}}"`:' title="Click to sort"';
    html+=`<th style="cursor:pointer;user-select:none" onclick="_sortITable(this)"${{tip}}>${{h}}${{i===0?' ▼':''}}</th>`;
  }});
  html+='</tr></thead><tbody>';
  rows.forEach(r=>{{
    html+='<tr>';
    keys.forEach(k=>{{
      let v=r[k]!==undefined&&r[k]!==''?r[k]:'—';
      if(k==='_priority') {{
        const sc=r._priority;
        const col=sc>=7?'#2ca02c':sc>=4?'#e07b39':'#aaa';
        html+=`<td data-val="${{sc}}" style="white-space:nowrap">`
             +`<span style="background:${{col}};color:#fff;font-weight:bold;border-radius:4px;`
             +`padding:2px 7px;font-size:11px">${{sc.toFixed(1)}}</span></td>`;
        return;
      }}
      if(k==='bindingSites') {{
        // ENCORI: bindingSites = clipExpNum (not position count); show max(bindingSites, sites.length)
        const posCount=(r.sites||[]).length;
        const v2=Math.max(parseInt(v)||0, posCount);
        html+=`<td data-val="${{v2}}">${{v2}}</td>`;
        return;
      }}
      if(k==='_mapped') {{
        const mapped=(r.sites||[]).length;
        const total=Math.max(parseInt(r.bindingSites||0)||0, mapped);
        // ENCORI entries use absolute genomic coords — cannot map to spliced position
        const encoriAbs=(r.source||'').includes('ENCORI') && mapped===0;
        if(encoriAbs) {{
          html+=`<td data-val="0" title="ENCORI 提供染色體絕對座標，無法直接對應至 spliced 序列位置"><span style="color:#bbb;font-size:11px">—</span></td>`;
        }} else if(mapped===0) {{
          html+=`<td data-val="0"><span style="color:#bbb;font-size:11px">0</span></td>`;
        }} else if(mapped<total) {{
          html+=`<td data-val="${{mapped}}" title="${{mapped}} / ${{total}} sites have hg19 position">`
               +`<span style="color:#e07b39;font-weight:bold">${{mapped}}</span>`
               +`<span style="color:#aaa;font-size:10px"> /${{total}}</span></td>`;
        }} else {{
          html+=`<td data-val="${{mapped}}">${{mapped}}</td>`;
        }}
        return;
      }}
      if(k==='circ_pos') {{
        if(tableType==='rbp' && r.sites && r.sites.length>0) {{
          // Per-site absolute hg19 positions with letter labels matching Structure tab
          const _lbl='abcdefghijklmnopqrstuvwxyz';
          let _posHtml=r.sites.map((s,i)=>{{
            const _letter=_lbl[i%26];
            const _rawPos=s.circ_pos||(s.circ_start+'–'+s.circ_end);
            const _absP=_absPos(_rawPos);
            return `<span style="white-space:nowrap;display:block">`
              +`<b style="background:#555;color:#fff;border-radius:3px;padding:0 4px;`
              +`font-size:10px;margin-right:3px;font-family:sans-serif">${{_letter}}</b>`
              +`<span style="font-family:monospace;font-size:10px">${{_absP}}</span></span>`;
          }}).join('');
          html+=`<td style="line-height:1.7;vertical-align:top">${{_posHtml}}</td>`;
          return;
        }}
        v=_absPos(v);
      }}
      if(k==='_seq_logo') {{
        const rawPos=String(r.circ_pos||'');
        // ENCORI: absolute coord "chrN:start-end"
        const absM=rawPos.match(/^(chr[^\s:]+):(\d+)[–-](\d+)/);
        // CircInteractome: relative "start–end"
        const relM=!absM&&rawPos.match(/(\\d+)[–-](\\d+)/);
        let seqChrom='', seqS=0, seqE=0;
        if(absM) {{
          seqChrom=absM[1]; seqS=parseInt(absM[2]); seqE=parseInt(absM[3]);
        }} else if(relM&&chrom) {{
          seqChrom=chrom;
          seqS=chromStart+parseInt(relM[1])-1;
          seqE=chromStart+parseInt(relM[2]);
        }}
        if(seqChrom) {{
          const len=seqE-seqS;
          if(len>200)
            v=`<span style="color:#aaa;font-size:10px">${{len}}bp</span>`;
          else
            v=`<span class="_seqCell" data-chrom="${{seqChrom}}" data-s="${{seqS}}" data-e="${{seqE}}"
                 style="color:#aaa;font-size:10px;font-family:monospace">⟳</span>`;
        }} else {{ v='—'; }}
      }}
      if(k==='source') v=_srcBadge(String(r[k]||''));
      html+='<td>'+v+'</td>';
    }});
    html+='</tr>';
  }});
  html+='</tbody></table>';
  const nCI=rows.filter(r=>r.source==='CircInteractome').length;
  const nEN=rows.filter(r=>r.source==='ENCORI').length;
  html+=`<p style="font-size:11px;color:#aaa;margin-top:6px">
    <span style="background:#6c757d;color:white;border-radius:3px;padding:1px 5px;font-size:10px">CircInteractome</span>
    ${{nCI}} records · TargetScan predictions · hg19
    &nbsp;&nbsp;
    <span style="background:#0077b6;color:white;border-radius:3px;padding:1px 5px;font-size:10px">ENCORI</span>
    ${{nEN}} records · CLIP-seq validated · hg38
    &nbsp;&nbsp; Binding Seq: hg19 UCSC (⟳ = loading)
    &nbsp;&nbsp;
    <strong>Priority</strong>: `
    +(tableType==='mirna'
      ?'seed strength (8mer=4…6mer=1) + CLIP>0(+3) + ENCORI(+2) + in_circ(+1)'
      :'sites log₂×2 max3 + internal(+2) + CLIP>0(+3) + ENCORI(+2) + in_circ(+1)')
    +(` · <span style="background:#2ca02c;color:#fff;border-radius:3px;padding:1px 5px;font-size:10px">≥7 high</span>`
    +` <span style="background:#e07b39;color:#fff;border-radius:3px;padding:1px 5px;font-size:10px">4–7 medium</span>`
    +` <span style="background:#aaa;color:#fff;border-radius:3px;padding:1px 5px;font-size:10px">&lt;4 low</span>`)
    +`</p>`;
  return html;
}}

function _seqLogoSpan(seq) {{
  const C={{A:'#2ca02c',T:'#d62728',G:'#e07b39',C:'#1a5c96',
            U:'#d62728',R:'#8c564b',Y:'#9467bd',K:'#bcbd22',
            M:'#17becf',S:'#7f7f7f',W:'#e377c2',N:'#aaa'}};
  const MAX=42;
  const s=seq.toUpperCase();
  const disp=s.length>MAX?s.slice(0,MAX)+'…':s;
  return '<span title="'+s+'" style="font-family:monospace;font-size:11px;letter-spacing:0.3px">'
    +disp.split('').map(c=>{{
      const col=C[c]||'#888';
      return `<span style="color:${{col}};font-weight:bold;`
            +`border-bottom:2.5px solid ${{col}};padding-bottom:1px;margin:0 0.3px">${{c}}</span>`;
    }}).join('')+'</span>';
}}

function _fetchSeqsInTable(containerId) {{
  const cells=Array.from(document.querySelectorAll('#'+containerId+' ._seqCell'));
  if(!cells.length) return;
  Promise.all(cells.map(cell=>{{
    const {{chrom,s,e}}=cell.dataset;
    const url='https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom='+chrom+';start='+s+';end='+e;
    return fetch(url)
      .then(r=>r.ok?r.json():null)
      .then(json=>{{
        if(json&&json.dna) cell.outerHTML=_seqLogoSpan(json.dna);
        else cell.textContent='N/A';
      }})
      .catch(()=>{{cell.textContent='✗';}});
  }}));
}}

function _spotlight(type, num, chip) {{
  const wrap = document.getElementById('cm-circle-wrap');
  if (!wrap) return;
  const isOn = chip.dataset.spotlight === '1';
  // Reset all arcs and chips for this type
  wrap.querySelectorAll(`._${{type}}arc`).forEach(g => {{
    g.style.opacity = '';
    g.querySelectorAll('path').forEach(p => {{ p.style.stroke=''; p.style.strokeWidth=''; }});
  }});
  wrap.querySelectorAll(`._${{type}}badge`).forEach(g => g.style.opacity = '');
  document.querySelectorAll(`[id^="_leg_${{type}}_"]`).forEach(c => {{
    c.dataset.spotlight = '0';
    c.style.boxShadow = '';
    c.style.fontWeight = '';
    c.style.background = '';
  }});
  if (!isOn) {{
    // Dim all arcs and badges
    wrap.querySelectorAll(`._${{type}}arc`).forEach(g => g.style.opacity = '0.04');
    wrap.querySelectorAll(`._${{type}}badge`).forEach(g => g.style.opacity = '0.15');
    // Spotlight selected arc
    const selArc = wrap.querySelector(`#_${{type}}arc_${{num}}`);
    if (selArc) {{
      selArc.style.opacity = '1';
      selArc.querySelectorAll('path').forEach(p => {{
        p.style.stroke = 'rgba(0,0,0,0.45)';
        p.style.strokeWidth = '1.5';
      }});
    }}
    // Spotlight selected badge
    const selBadge = wrap.querySelector(`#_${{type}}b_${{num}}`);
    if (selBadge) selBadge.style.opacity = '1';
    chip.dataset.spotlight = '1';
    chip.style.boxShadow = '0 0 0 2px #333';
    chip.style.fontWeight = 'bold';
    chip.style.background = '#f5f5f5';
  }}
}}

// ── Level 2: site-panel (appears when exactly 1 molecule visible) ────────────
function _updateSitePanel(wrap) {{
  const visMi=[...wrap.querySelectorAll('._miarc')].filter(g=>g.style.display!=='none');
  const visRbp=[...wrap.querySelectorAll('._rbparc')].filter(g=>g.style.display!=='none');
  const total=visMi.length+visRbp.length;
  if(total===1){{
    const g=visMi.length===1?visMi[0]:visRbp[0];
    const type=visMi.length===1?'mi':'rbp';
    const num=parseInt(g.id.replace(/^.*_(\d+)$/,'$1'));
    _showSitePanel(wrap,type,num);
  }}else{{
    _hideSitePanel(wrap);
  }}
}}
function _showSitePanel(wrap,type,num){{
  const data=type==='mi'?(window._circ_miData||{{}})[num]:(window._circ_rbpData||{{}})[num];
  const panel=document.getElementById('circ-site-panel');
  if(!panel||!data){{_hideSitePanel(wrap);return;}}
  wrap.querySelectorAll('._arc_lbl').forEach(g=>g.style.visibility='hidden');
  wrap.querySelectorAll(`._arc_lbl[data-arc-type="${{type}}"][data-arc-num="${{num}}"]`).forEach(g=>g.style.visibility='');
  let h=`<div style="font-weight:600;margin-bottom:4px;color:#446;font-size:10px;word-break:break-all">${{data.name}}</div>`;
  h+='<div style="font-size:9px;color:#888;margin-bottom:5px">點擊字母切換顯示 / Click to toggle</div>';
  h+='<div style="max-height:320px;overflow-y:auto">';
  data.arcs.forEach(arc=>{{
    h+=`<div onclick="_toggleSite('${{arc.id}}',this)" style="display:flex;align-items:center;gap:5px;cursor:pointer;padding:2px 4px;border-radius:4px;margin-bottom:3px;transition:opacity .15s">
      <span style="background:${{data.color}};color:white;border-radius:50%;width:16px;height:16px;display:inline-flex;align-items:center;justify-content:center;font-size:9px;font-weight:bold;flex-shrink:0">${{arc.letter}}</span>
      <span style="flex:1;font-size:9px;color:#333">site ${{arc.letter}}</span>
      <span class="_eye" style="font-size:11px">👁</span>
    </div>`;
  }});
  h+='</div>';
  panel.innerHTML=h;panel.style.display='';
  // Highlight arcs so the letter badge is clearly on the arc
  const aGrp=wrap.querySelector(`#_${{type}}arc_${{num}}`);
  if(aGrp)aGrp.querySelectorAll('path').forEach(p=>{{
    p.style.stroke='white';p.style.strokeWidth='1.5';p.style.strokeLinejoin='round';p.style.opacity='1';
  }});
}}
function _hideSitePanel(wrap){{
  const panel=document.getElementById('circ-site-panel');
  if(panel){{panel.style.display='none';panel.innerHTML='';}}
  if(wrap){{
    wrap.querySelectorAll('._arc_lbl').forEach(g=>g.style.visibility='hidden');
    wrap.querySelectorAll('._miarc path,._rbparc path').forEach(p=>{{
      p.style.stroke='';p.style.strokeWidth='';p.style.opacity='';
    }});
  }}
}}
function _toggleSite(arcId,row){{
  const g=document.getElementById(arcId);if(!g)return;
  const nowHidden=g.style.display==='none';
  g.style.display=nowHidden?'':'none';
  const eye=row?row.querySelector('._eye'):null;
  if(eye)eye.textContent=nowHidden?'👁':'🙈';
  if(row)row.style.opacity=nowHidden?'1':'0.4';
}}

function _toggleAnn(type, letter, row) {{
  // Toggle overlap-annotation arc group (boundary lines + letter badge)
  const wrap = document.getElementById('cm-circle-wrap');
  if (!wrap) return;
  const g = wrap.querySelector(`#_ann_${{type}}_${{letter}}`);
  if (!g) return;
  const nowHidden = g.style.display === 'none';
  g.style.display = nowHidden ? '' : 'none';
  const eye = row ? row.querySelector('._eye') : null;
  if (eye) eye.textContent = nowHidden ? '👁' : '🙈';
  if (row) row.style.opacity = nowHidden ? '1' : '0.4';
}}

function _toggleBadge(type, num, chip) {{
  const wrap = document.getElementById('cm-circle-wrap');
  if (!wrap) return;
  const badge = wrap.querySelector(`#_${{type}}b_${{num}}`);
  const arc   = wrap.querySelector(`#_${{type}}arc_${{num}}`);
  const nowHidden = badge ? (badge.style.display === 'none') : false;
  if (badge) badge.style.display = nowHidden ? '' : 'none';
  if (arc)   arc.style.display   = nowHidden ? '' : 'none';
  if (chip) {{
    chip.style.opacity        = nowHidden ? '1' : '0.3';
    chip.style.textDecoration = nowHidden ? '' : 'line-through';
  }}
  _updateSitePanel(wrap);
}}

function _toggleAll(type, show, container) {{
  const wrap = document.getElementById('cm-circle-wrap');
  if (!wrap) return;
  // Clear any spotlight state first
  wrap.querySelectorAll(`._${{type}}arc`).forEach(g => {{
    g.style.opacity = '';
    g.querySelectorAll('path').forEach(p => {{ p.style.stroke=''; p.style.strokeWidth=''; }});
  }});
  wrap.querySelectorAll(`._${{type}}badge`).forEach(g => g.style.opacity = '');
  if (container) {{
    container.querySelectorAll(`[id^="_leg_${{type}}_"]`).forEach(c => {{
      c.dataset.spotlight='0'; c.style.boxShadow=''; c.style.fontWeight=''; c.style.background='';
    }});
  }}
  wrap.querySelectorAll(`._${{type}}badge`).forEach(el => el.style.display = show ? '' : 'none');
  wrap.querySelectorAll(`._${{type}}arc`).forEach(el   => el.style.display = show ? '' : 'none');
  if (container) {{
    container.querySelectorAll(`[id^="_leg_${{type}}_"]`).forEach(chip => {{
      chip.style.opacity        = show ? '1' : '0.3';
      chip.style.textDecoration = show ? '' : 'line-through';
    }});
  }}
  _updateSitePanel(wrap);
}}

// ── Main heatmap dynamic update ───────────────────────────────────────────────
function updateMainHeatmap() {{
  // Top-N heatmap section was removed from the report; this is now only
  // called (harmlessly) from switchDEMethod(). No-op if its DOM is gone.
  const nInput=document.getElementById('heatmap-n-input');
  if(!nInput)return;
  const _hmData=_HEATMAP_DATA_CACHE||FULL_HEATMAP_DATA;
  if(!_hmData||typeof Plotly==='undefined')return;
  const n=Math.max(1,Math.min(50,parseInt(nInput.value)||10));
  nInput.value=n;
  const upIds=(_hmData.up_order||[]).slice(0,n);
  const dnIds=(_hmData.dn_order||[]).slice(0,n);
  const allIds=[...upIds,...dnIds];
  if(allIds.length<2)return;
  const samps=_hmData.samples||[];
  const conds=_hmData.conditions||{{}};
  const rows=_hmData.rows||{{}};
  const validIds=allIds.filter(id=>rows[id]);
  const zMatrix=validIds.map(id=>rows[id].z);
  const yLabels=validIds.map(id=>rows[id].label||id);
  const nUp=upIds.filter(id=>rows[id]).length;
  const nDn=dnIds.filter(id=>rows[id]).length;
  const TUMOR_COL='#d62728',NORMAL_COL='#2CA02C';
  const grps=[];let cur=null;
  samps.forEach((s,i)=>{{
    const c=conds[s]||'';
    const col=c==='{tumor_label}'?TUMOR_COL:c==='{normal_label}'?NORMAL_COL:'#888';
    if(!cur||cur.c!==c){{cur={{c,col,s:i,e:i}};grps.push(cur);}}else cur.e=i;
  }});
  const groupShapes=grps.map(g=>{{return{{type:'rect',xref:'x',yref:'paper',
    x0:g.s-0.45,x1:g.e+0.45,y0:1.04,y1:1.09,fillcolor:g.col,line:{{width:0}}}}}});
  const groupAnno=grps.map(g=>{{return{{xref:'x',yref:'paper',x:(g.s+g.e)/2,y:1.065,
    yanchor:'middle',xanchor:'center',text:g.c,showarrow:false,
    font:{{size:13,color:'white',family:'sans-serif'}}}}}});
  Plotly.react('main-heatmap-plot',[{{
    type:'heatmap',z:zMatrix,x:samps,y:yLabels,
    colorscale:[[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']],
    zmid:0,colorbar:{{title:'z-score<br>({normal_label}-<br>centered)'}},
    hovertemplate:'<b>%{{y}}</b><br>%{{x}}<br>z-score: %{{z:.2f}}<extra></extra>',
  }}],{{
    title:{{text:''}},
    yaxis:{{tickfont:{{size:8}},autorange:'reversed'}},
    height:Math.max(420,validIds.length*22+130),
    plot_bgcolor:'white',paper_bgcolor:'white',
    margin:{{t:70,l:300}},
    shapes:groupShapes,
    annotations:groupAnno,
  }});
  const titleEl=document.getElementById('heatmap-section-title');
  if(titleEl)titleEl.textContent=`Heatmap (top ${{nUp}} significant up + ${{nDn}} significant down DE circRNAs)`;
  const maxUp=(_hmData.up_order||[]).length;
  const maxDn=(_hmData.dn_order||[]).length;
  const statusEl=document.getElementById('heatmap-status');
  if(statusEl){{statusEl.dataset.n=validIds.length;statusEl.dataset.u=maxUp;statusEl.dataset.d=maxDn;statusEl.textContent=_LS[_LANG||'zh'].hmStatus(validIds.length,maxUp,maxDn);}}
}}

// ── Biomarker top-N update ────────────────────────────────────────────────────
function updateBiomarkerN() {{
  const inp = document.getElementById('bm-n-input');
  if (!inp) return;
  const v = Math.max(1, parseInt(inp.value) || 30);
  inp.value = v;
  // Find the currently active DE method
  let curMethod = '{de_method}';
  document.querySelectorAll('.msw-btn').forEach(b => {{
    if (b.classList.contains('active')) {{
      const m = (b.getAttribute('onclick') || '').match(/'([^']+)'/);
      if (m) curMethod = m[1];
    }}
  }});
  const md = ALL_DE_METHODS && ALL_DE_METHODS[curMethod];
  _renderBiomarkerTable(curMethod, md);
}}

// ── DE method switcher ────────────────────────────────────────────────────────
function switchDEMethod(method) {{
  const md=ALL_DE_METHODS&&ALL_DE_METHODS[method];
  if(!md)return;
  // Update stat boxes
  const s=md.stats||{{}};
  const sigEl=document.getElementById('stat-n-sig');
  const upEl=document.getElementById('stat-n-up');
  const dnEl=document.getElementById('stat-n-dn');
  if(sigEl)sigEl.textContent=s.n_sig!=null?s.n_sig:'—';
  if(upEl) upEl.textContent =s.n_up !=null?s.n_up :'—';
  if(dnEl) dnEl.textContent =s.n_dn !=null?s.n_dn :'—';
  // Update switcher button state
  document.querySelectorAll('.msw-btn').forEach(b=>{{
    const m=b.getAttribute('onclick')||'';
    b.classList.toggle('active',m.includes("'"+method+"'"));
  }});
  // Update volcano via Plotly.react
  const volData=md.volcano||[];
  _CURRENT_VOLCANO_DATA=volData;
  const ns_x=[],ns_y=[],up_x=[],up_y=[],dn_x=[],dn_y=[];
  volData.forEach(d=>{{
    if(d[2]==='N'){{ns_x.push(d[0]);ns_y.push(d[1]);}}
    else if(d[2]==='U'){{up_x.push(d[0]);up_y.push(d[1]);}}
    else{{dn_x.push(d[0]);dn_y.push(d[1]);}}
  }});
  if(typeof Plotly!=='undefined'&&document.getElementById('main-volcano-plot')){{
    const mLabels={{'edgeR_ciriquant':'edgeR (FSJ offset)','deseq2':'DESeq2','limma':'limma-voom'}};
    Plotly.react('main-volcano-plot',[
      {{x:ns_x,y:ns_y,mode:'markers',name:'NS',
        marker:{{color:'rgba(150,150,150,0.35)',size:4,line:{{width:0}}}},hoverinfo:'skip'}},
      {{x:up_x,y:up_y,mode:'markers',name:'Up',
        marker:{{color:'#d62728',size:6,opacity:0.8}}}},
      {{x:dn_x,y:dn_y,mode:'markers',name:'Down',
        marker:{{color:'#2CA02C',size:6,opacity:0.8}}}},
    ],{{
      title:{{text:'Volcano Plot ['+(mLabels[method]||method)+']',font:{{size:14}}}},
      xaxis:{{title:'log₂ Fold Change (Tumor / Normal)',showgrid:true,gridcolor:'#f0f0f0',zeroline:false}},
      yaxis:{{title:'−log₁₀('+_SIG_LABEL+')',showgrid:true,gridcolor:'#f0f0f0',zeroline:false}},
      height:500,plot_bgcolor:'white',paper_bgcolor:'white',
      legend:{{title:'',orientation:'h',y:1.02,x:0}},
      margin:{{t:60}},
      shapes:[
        {{type:'line',x0:_LFC,x1:_LFC,y0:0,y1:1,yref:'paper',line:{{dash:'dot',color:'#aaa',width:1}}}},
        {{type:'line',x0:-_LFC,x1:-_LFC,y0:0,y1:1,yref:'paper',line:{{dash:'dot',color:'#aaa',width:1}}}},
        {{type:'line',x0:0,x1:1,xref:'paper',y0:-Math.log10(_FDR),y1:-Math.log10(_FDR),line:{{dash:'dot',color:'#aaa',width:1}}}},
      ],
    }});
  }}
  // Update heatmap data cache and re-render
  // Conditions (sample→group mapping) are identical across all DE methods;
  // fall back to FULL_HEATMAP_DATA.conditions if the alternate method's is empty.
  if(md.heatmap){{
    const _hm=Object.assign({{}},md.heatmap);
    if((!_hm.conditions||!Object.keys(_hm.conditions).length)&&FULL_HEATMAP_DATA){{
      _hm.conditions=FULL_HEATMAP_DATA.conditions;
    }}
    _HEATMAP_DATA_CACHE=_hm;
    updateMainHeatmap();
  }}
  // Re-render Top DE tables
  _renderDETables(method, md);
  // Re-render Biomarker table (per-method ranking; fallback to graying if no bm_table)
  _renderBiomarkerTable(method, md);
  // Update biomarker score distribution plots
  _updateScoreDist(method, md);
}}

function _renderDETables(method, md) {{
  const sec=document.getElementById('de-tables-section');
  if(!sec)return;
  const dt=md.de_table||{{}};
  const cols=dt.cols||[];
  const mLabels={{'edgeR_ciriquant':'edgeR (FSJ offset)','deseq2':'DESeq2','limma':'limma-voom'}};
  const mName=mLabels[method]||method;
  const s=md.stats||{{}};
  const heading=document.getElementById('de-tables-heading');
  if(heading)heading.innerHTML=`Top Differentially Expressed circRNAs <span class="method-tag">${{mName}}</span>`;

  function _mkTable(rows,direction) {{
    if(!rows||!rows.length)return '<p style="color:#aaa;font-size:13px">No significant circRNAs.</p>';
    const color=direction==='up'?'#d62728':'#2CA02C';
    const arrow=direction==='up'?'&#8593;':'&#8595;';
    const label=direction==='up'?`Up-regulated — ${{rows.length}} circRNAs`:`Down-regulated — ${{rows.length}} circRNAs`;
    const tid=`de-${{direction}}-table`;
    let html=`<h3 style="color:${{color}}">${{arrow}} ${{label}}</h3><div style="overflow-x:auto"><table id="${{tid}}" class="table" border="0"><thead><tr>`;
    cols.forEach(c=>{{ html+=`<th>${{c}}</th>`; }});
    html+='</tr></thead><tbody>';
    rows.forEach(row=>{{
      html+='<tr>';
      row.forEach((v,i)=>{{
        const col=cols[i]||'';
        if(col==='circ_id'&&v){{
          html+=`<td><a class="circ-link" onclick="showCircDetail('${{v}}')">${{v}}</a></td>`;
        }} else if(col==='circbase_id'){{
          html+=`<td>${{_cbLink(v)}}</td>`;
        }} else if(col==='log2FC'&&v!=null){{
          html+=`<td>${{(+v).toFixed(3)}}</td>`;
        }} else if(col==='p-value'&&v!=null){{
          html+=`<td>${{(+v).toExponential(3)}}</td>`;
        }} else if(col==='gene_name'){{
          html+=`<td>${{(v==null||v==='intergenic')?'—':v}}</td>`;
        }} else {{
          html+=`<td>${{v==null?'—':v}}</td>`;
        }}
      }});
      html+='</tr>';
    }});
    html+='</tbody></table></div>';
    return html;
  }}
  const upHtml=_mkTable(dt.up,'up');
  const dnHtml=_mkTable(dt.dn,'dn');
  // Replace only the table content, keep the heading div
  const tableWrap=sec.querySelector('#de-dynamic-tables');
  if(tableWrap){{ tableWrap.innerHTML=upHtml+dnHtml; }}
  else {{
    const wrap=document.createElement('div');
    wrap.id='de-dynamic-tables';
    wrap.innerHTML=upHtml+dnHtml;
    // Remove static tables, headings, and CSV download bar divs
    const staticTbls=sec.querySelectorAll('table,h3,.tbl-dl-bar');
    staticTbls.forEach(el=>el.remove());
    sec.appendChild(wrap);
  }}
  _makeSortable('de-up-table');
  _makeSortable('de-dn-table');
}}

function _updateScoreDist(method, md) {{
  const sd=md.score_dist;
  if(!sd||typeof Plotly==='undefined')return;
  const n=sd.scores.length;
  const mLabels={{'edgeR_ciriquant':'edgeR (FSJ offset)','deseq2':'DESeq2','limma':'limma-voom'}};
  // Update title
  const t=document.getElementById('bm-dist-title');
  if(t)t.textContent=`Score Distribution — all ${{n}} significant DE circRNAs (${{mLabels[method]||method}})`;
  // Scatter plot
  const top30=sd.scores.slice(0,30), rest=sd.scores.slice(30);
  const xTop=top30.map((_,i)=>i+1), xRest=rest.map((_,i)=>i+31);
  if(document.getElementById('bm-scatter-plot')){{
    Plotly.react('bm-scatter-plot',[
      {{x:xRest,y:rest,mode:'markers',name:`Rank 31–${{n}}`,
        marker:{{color:'rgba(120,120,120,0.45)',size:5}}}},
      {{x:xTop,y:top30,mode:'markers',name:'Top 30 (table)',
        marker:{{color:'#d62728',size:7,line:{{color:'white',width:0.8}}}}}},
    ],{{
      xaxis:{{title:'Rank'}},
      yaxis:{{title:'Biomarker Score',range:[Math.max(0,sd.scores[n-1]-0.05),Math.min(1,sd.scores[0]+0.05)]}},
      height:320,margin:{{t:20,b:60,l:70,r:40}},
      plot_bgcolor:'white',paper_bgcolor:'white',
      shapes:[{{type:'line',x0:30.5,x1:30.5,y0:0,y1:1,yref:'paper',
                line:{{dash:'dot',color:'#d62728',width:1.2}}}}],
      annotations:[{{x:31.5,y:0.98,yref:'paper',text:'Top 30',showarrow:false,
                     font:{{size:10,color:'#d62728'}},xanchor:'left'}}],
      legend:{{x:0.75,y:0.95,bgcolor:'rgba(255,255,255,0.8)'}},
    }});
  }}
  // Histogram
  if(document.getElementById('bm-hist-plot')){{
    const swColor=(sd.sw_p&&sd.sw_p<0.05)?'#d62728':'#2CA02C';
    const swText=sd.sw_w?`Shapiro-Wilk: W = ${{sd.sw_w.toFixed(4)}}, p = ${{sd.sw_p.toExponential(2)}}`:'';
    const swConc=sd.sw_p?_LS[_LANG||'zh'].swConc(sd.sw_p):'';
    Plotly.react('bm-hist-plot',[
      {{x:sd.scores,type:'histogram',nbinsx:30,name:_LS[_LANG||'zh'].obsDist,histnorm:'probability density',
        marker:{{color:'rgba(44,119,214,0.55)',line:{{color:'rgba(44,119,214,0.9)',width:0.8}}}}}},
      {{x:sd.x_norm,y:sd.y_norm,mode:'lines',
        name:`Normal(μ=${{sd.mu.toFixed(3)}}, σ=${{sd.sd.toFixed(3)}})`,
        line:{{color:'#d62728',width:2}}}},
    ],{{
      xaxis:{{title:'Biomarker Score'}},
      yaxis:{{title:'Probability Density'}},
      height:320,margin:{{t:40,b:90,l:70,r:40}},
      legend:{{x:0.65,y:0.95,bgcolor:'rgba(255,255,255,0.85)'}},
      plot_bgcolor:'white',paper_bgcolor:'white',
      annotations:[{{x:0.5,y:-0.28,xref:'paper',yref:'paper',
        text:`${{swText}}    <b style="color:${{swColor}}">${{swConc}}</b>`,
        showarrow:false,align:'center',xanchor:'center',yanchor:'top',
        bgcolor:'rgba(255,255,255,0.88)',bordercolor:'#ccc',borderwidth:1,
        font:{{size:10}}}}],
    }});
  }}
}}

function _renderBiomarkerTable(method, md) {{
  const bt = md && md.bm_table;
  const tbody = document.querySelector('#tbl_biomarker tbody');
  if (!tbody || !bt || !bt.rows || !bt.rows.length) {{
    // Fallback: gray out rows not significant under this method
    _updateBiomarkerHighlight(md && md.sig_ids || []);
    return;
  }}
  const cols = bt.cols || [];
  const mLabels = {{'edgeR_ciriquant':'edgeR (FSJ offset)','deseq2':'DESeq2','limma':'limma-voom'}};
  // Respect top-N input
  const _bmInp = document.getElementById('bm-n-input');
  const _bmN = _bmInp ? Math.max(1, parseInt(_bmInp.value) || 30) : 30;
  const displayRows = bt.rows.slice(0, _bmN);
  // Renumber rank column (index 0 = 'rank')
  const rankIdx = cols.indexOf('rank');
  const rows = displayRows.map((row, i) => {{
    const r = [...row];
    if (rankIdx >= 0) r[rankIdx] = i + 1;
    return r;
  }});
  tbody.innerHTML = rows.map(row => {{
    const nsigIdx = cols.indexOf('n_sig_methods');
    const nsig = nsigIdx >= 0 ? (parseInt(row[nsigIdx]) || 1) : 1;
    const fw = nsig >= 3 ? 'font-weight:bold;' : '';
    const cells = row.map((v, i) => {{
      const col = cols[i] || '';
      if (col === 'circ_id' && v) {{
        const hasDat = CIRC_DATA && CIRC_DATA[v];
        const tip = hasDat ? 'in interaction data' : 'no interaction data pre-fetched';
        return `<td><a class="circ-link" onclick="showCircDetail('${{v}}')" title="${{tip}}">${{v}}</a></td>`;
      }}
      if (col === 'log2FC' && v != null) return `<td>${{(+v).toFixed(2)}}</td>`;
      if (col === 'biomarker_score' && v != null) return `<td>${{(+v).toFixed(3)}}</td>`;
      if (col === 'circbase_id') return `<td>${{_cbLink(v)}}</td>`;
      if (v === '' || v == null) return '<td>—</td>';
      return `<td>${{v}}</td>`;
    }}).join('');
    return `<tr data-nsig="${{nsig}}" style="${{fw}}">${{cells}}</tr>`;
  }}).join('');
  // Reset filter buttons to 全部 + update counts
  const allRows = [...tbody.querySelectorAll('tr')];
  const n_all = allRows.length;
  const n2 = allRows.filter(r => (parseInt(r.getAttribute('data-nsig')||'1')) >= 2).length;
  const n3 = allRows.filter(r => (parseInt(r.getAttribute('data-nsig')||'1')) >= 3).length;
  document.querySelectorAll('.bm-filter-btn').forEach((b, i) => {{
    b.style.background = 'white'; b.style.color = '#555';
    if (i === 0) {{ b.textContent = _LS[_LANG||'zh'].bmAll(n_all); b.style.background = '#2c6fad'; b.style.color = 'white'; }}
    else if (i === 1) b.textContent = _LS[_LANG||'zh'].bm2(n2);
    else if (i === 2) b.textContent = _LS[_LANG||'zh'].bm3(n3);
  }});
  // Update section headings
  const bm_h2 = document.getElementById('biomarker-section-title');
  if (bm_h2) bm_h2.textContent = `Biomarker Candidates (top ${{n_all}} by composite score) [${{mLabels[method]||method}}]`;
  const bmStatus = document.getElementById('bm-n-status');
  if (bmStatus) {{ const tot=bt.rows.length; bmStatus.textContent=n_all<tot?`${{n_all}} / ${{tot}} shown`:''; }}
  // Re-attach sort listeners (tbody was rebuilt)
  _makeSortable('tbl_biomarker');
}}

function _updateBiomarkerHighlight(sigIds) {{
  const sigSet=new Set(sigIds);
  const bmSec=document.getElementById('biomarker-section');
  if(!bmSec)return;
  const rows=bmSec.querySelectorAll('table tbody tr');
  rows.forEach(row=>{{
    const link=row.querySelector('a.circ-link');
    const m=link&&link.getAttribute('onclick')&&link.getAttribute('onclick').match(/'([^']+)'/);
    const cid=m?m[1]:null;
    row.style.opacity=(cid&&!sigSet.has(cid))?'0.25':'1';
    row.title=(cid&&!sigSet.has(cid))?_LS[_LANG||'zh'].bmNSig:'';
  }});
}}

document.addEventListener('keydown',e=>{{if(e.key==='Escape')closeCircModal();}});
// Make static tables sortable on page load
_makeSortable('tbl_isoform');
_makeSortable('tbl_biomarker');
// Replace static DE tables with dynamic (sortable) ones on first load
(function() {{
  var initMethod = '{de_method}';
  if (ALL_DE_METHODS && ALL_DE_METHODS[initMethod]) {{
    switchDEMethod(initMethod);
  }}
}})();

// ── Clustering Heatmap (always expanded, renders immediately) ────────────────
(function() {{
  const sec = document.getElementById('clust-heatmap-section');
  if (!sec || !CLUST_HEATMAP_DATA || typeof Plotly === 'undefined') return;
  _drawClustHeatmap();
}})();

// Shared renderer: dendrogram (left) + heatmap (right) two-panel layout.
// Used by both the full-report clustering heatmap and the per-circRNA
// modal's windowed "clustering neighborhood" mini-heatmap.
//   order/rows/samps/conds: same shape as CLUST_HEATMAP_DATA fields, but
//     order/rows may be a windowed subset (mini-heatmap use case).
//   dendro: {{icoord, dcoord, max_dist}} already local to `order` (x=0 is
//     the first row, x=10*(order.length) the last) -- see _sliceDendro().
//   opts.highlightIdx: index into `order` to draw an orange highlight box
//     around (mini-heatmap's clicked circRNA); null/undefined = no box.
//   opts.title: optional Plotly chart title.
function _renderClustPanel(div, order, rows, samps, conds, dendro, opts) {{
  opts = opts || {{}};
  const n = order.length;
  if (n === 0 || samps.length === 0) return;

  // Scipy leaf convention: row i sits at y = 10*i+5 (5,15,25,...)
  const y_positions = order.map((_, i) => 10 * i + 5);
  const ylbls       = order.map(id => rows[id] ? rows[id].label || id : id);

  const fontSz = n > 200 ? 4 : n > 100 ? 6 : n > 50 ? 8 : n > 25 ? 9 : 10;
  const rowH   = n > 200 ? 4 : n > 100 ? 6 : n > 50 ? 8 : n > 25 ? 14 : 20;
  const plotH  = opts.height || Math.min(700, Math.max(280, n * rowH + 140));

  const TUMOR_COL='#d62728', NORMAL_COL='#2CA02C';
  const grps=[]; let cur=null;
  samps.forEach((s,i) => {{
    const c = conds[s] || '';
    const col = c==='{tumor_label}' ? TUMOR_COL : c==='{normal_label}' ? NORMAL_COL : '#888';
    if (!cur || cur.c !== c) {{ cur={{c,col,s:i,e:i}}; grps.push(cur); }} else cur.e = i;
  }});
  const groupShapes = grps.map(g => ({{
    type:'rect', xref:'x2', yref:'paper',
    x0:g.s-0.45, x1:g.e+0.45, y0:1.02, y1:1.07,
    fillcolor:g.col, line:{{width:0}}
  }}));
  const groupAnno = grps.map(g => ({{
    xref:'x2', yref:'paper', x:(g.s+g.e)/2, y:1.045,
    yanchor:'middle', xanchor:'center', text:g.c, showarrow:false,
    font:{{size:11, color:'white', family:'sans-serif'}}
  }}));

  const z = order.map(id => rows[id] ? rows[id].z : samps.map(() => 0));
  // customdata MUST be 2D [n_rows × n_cols] for heatmap hovertemplate to work
  const customdata = order.map((id, i) =>
    samps.map(() => [
      rows[id] ? rows[id].lfc  : 0,
      rows[id] ? rows[id].pval : 1,
      ylbls[i],
    ])
  );

  const traces = [];

  // ── Dendrogram traces (left panel, xaxis) ──
  const hasDendro = !!(dendro && dendro.icoord && dendro.icoord.length);
  const maxD = hasDendro ? (dendro.max_dist || 1) : 1;
  if (hasDendro) {{
    dendro.icoord.forEach((ic, k) => {{
      traces.push({{
        type: 'scatter',
        x: dendro.dcoord[k].map(v => -v),   // negate: root far-left, leaves near 0
        y: ic,
        mode: 'lines',
        line: {{color:'#555', width: n > 150 ? 0.5 : 0.8}},
        showlegend: false,
        hoverinfo: 'none',
        xaxis: 'x',
        yaxis: 'y',
      }});
    }});
  }}

  // ── Heatmap trace (main panel, xaxis2) ──
  // Colorbar sits on the LEFT (in margin.l) to avoid overlap with right-side ID labels
  traces.push({{
    type: 'heatmap',
    z: z,
    x: samps,
    y: y_positions,
    colorscale: [[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']],
    zmid: 0,
    colorbar: {{
      title: {{text:'z-score', side:'right', font:{{size:10}}}},
      thickness: 10, len: 0.4,
      x: 0.0, xanchor: 'right',   // left of the dendrogram panel
      y: 0.5, yanchor: 'middle',
    }},
    customdata: customdata,
    hovertemplate: '<b>%{{customdata[2]}}</b><br>Sample: %{{x}}<br>z-score: %{{z:.2f}}<br>log2FC: %{{customdata[0]:.2f}}<br>p-value: %{{customdata[1]:.3g}}<extra></extra>',
    xaxis: 'x2',
    yaxis: 'y',
  }});

  const dendroFrac = hasDendro ? 0.13 : 0;
  const gapFrac    = hasDendro ? 0.01  : 0;
  // Right margin: accommodate circRNA ID labels on the right
  const maxLblLen = Math.max(...ylbls.map(s => s.length));
  const rMargin   = Math.min(300, Math.max(130, maxLblLen * fontSz * 0.65));

  const hlShape = (opts.highlightIdx != null) ? [{{
    type:'rect', xref:'x2', yref:'y',
    x0:-0.5, x1:samps.length-0.5,
    y0:10*opts.highlightIdx, y1:10*opts.highlightIdx+10,
    line:{{color:'#ff8c00', width:2.5}}, fillcolor:'rgba(255,140,0,0.07)'
  }}] : [];

  Plotly.newPlot(div, traces, {{
    height: plotH,
    title: opts.title ? {{text: opts.title, font:{{size:11}}}} : undefined,
    plot_bgcolor: 'white',
    paper_bgcolor: 'white',
    margin: {{t: opts.title ? 85 : 60, l:72, r:rMargin, b:60}},   // l:72 reserves space for the left-side colorbar
    shapes: [...groupShapes, ...hlShape],
    annotations: groupAnno,
    dragmode: 'zoom',
    // ── Dendrogram axis (left narrow panel) ──
    xaxis: {{
      domain: [0, dendroFrac],
      anchor: 'y',
      range: [-(maxD * 0.82), maxD * 0.05],   // clip near-root arms (matches Python 80% filter)
      showticklabels: false,
      showgrid: false,
      zeroline: false,
      fixedrange: false,
    }},
    // ── Heatmap axis (main right panel) ──
    xaxis2: {{
      domain: [dendroFrac + gapFrac, 1.0],
      anchor: 'y',
      showgrid: false,
      zeroline: false,
      showline: false,
      fixedrange: false,
    }},
    // ── Shared y-axis: labels appear on the RIGHT of the heatmap ──
    yaxis: {{
      anchor: 'x2',     // anchor to heatmap's right edge
      side:   'right',  // labels to the right of the heatmap panel
      tickmode: 'array',
      tickvals: y_positions,
      ticktext: ylbls,
      tickfont: {{size: fontSz}},
      autorange: 'reversed',
      fixedrange: false,
      showgrid: false,
      zeroline: false,   // hide y=0 line (appears above row-1 in reversed axis)
      showline: false,
    }},
  }}, {{
    responsive: true,
    displayModeBar: true,
    scrollZoom: true,                              // mouse-wheel / pinch to zoom
    modeBarButtonsToRemove: ['select2d','lasso2d'],
  }});
  // ── Zoom: dynamically scale y-axis tick labels as row height changes ──
  let _chBusy = false;
  div.on('plotly_relayout', function(ev) {{
    if (_chBusy) return;
    const yl = div._fullLayout && div._fullLayout.yaxis;
    if (!yl || !yl.range) return;
    const span = Math.abs(yl.range[1] - yl.range[0]);
    const vis  = Math.max(1, Math.round(span / 10));   // number of visible rows
    const newSz = vis > 120 ? 4 : vis > 60 ? 6 : vis > 30 ? 8 : vis > 15 ? 10 : vis > 6 ? 12 : 14;
    const curSz = (yl.tickfont && yl.tickfont.size) || 0;
    if (newSz === curSz) return;
    _chBusy = true;
    Plotly.relayout(div, {{'yaxis.tickfont.size': newSz}})
      .then( () => {{ _chBusy = false; }} )
      .catch( () => {{ _chBusy = false; }} );
  }});
}}

// Extract the local sub-dendrogram for a contiguous leaf-index window
// [lo, hi) of the full cluster order. Standard dendrograms never cross
// branches, so any merge segment whose x-endpoints both fall inside the
// window's x-range belongs entirely to leaves within that window; segments
// reaching outside are simply dropped (their partner leaf isn't visible).
function _sliceDendro(dendro, lo, hi) {{
  if (!dendro || !dendro.icoord || !dendro.icoord.length) return null;
  const xlo = 10*lo, xhi = 10*hi;
  const icoord=[], dcoord=[];
  dendro.icoord.forEach((ic, k) => {{
    if (ic.every(x => x >= xlo - 1e-6 && x <= xhi + 1e-6)) {{
      icoord.push(ic.map(x => x - xlo));
      dcoord.push(dendro.dcoord[k]);
    }}
  }});
  if (!icoord.length) return null;
  const maxD = Math.max(...dcoord.map(d => Math.max(...d)));
  return {{icoord, dcoord, max_dist: maxD || dendro.max_dist || 1}};
}}

function _drawClustHeatmap() {{
  const cd = CLUST_HEATMAP_DATA;
  if (!cd || typeof Plotly === 'undefined') return;
  const div = document.getElementById('clust-heatmap-plot');
  if (!div) return;
  _renderClustPanel(div, cd.order || [], cd.rows || {{}}, cd.samples || [],
                     cd.conditions || {{}}, cd.dendro, {{}});
}}
</script>"""

    _modal_html = """
<div id="circ-modal" class="circ-modal" onclick="if(event.target===this)closeCircModal()">
  <div class="circ-modal-box">
    <span class="circ-modal-close" onclick="closeCircModal()">&#10005;</span>
    <div class="circ-modal-title" id="cm-title"></div>
    <div class="circ-modal-sub" id="cm-sub"></div>
    <div class="ctab-bar">
      <button class="ctab-btn active" data-tab="exon"    onclick="_switchTab('exon')"   >&#11835; Circular Structure</button>
      <button class="ctab-btn"        data-tab="mirna"   onclick="_switchTab('mirna')"  >&#128250; miRNA Sponge</button>
      <button class="ctab-btn"        data-tab="rbp"     onclick="_switchTab('rbp')"    >&#129520; RBP Binding</button>
      <button class="ctab-btn"        data-tab="volcano" onclick="_switchTab('volcano')">&#128200; Volcano</button>
      <button class="ctab-btn" id="cm-tab-heatmap" data-tab="heatmap" onclick="_switchTab('heatmap')">&#128293; Heatmap</button>
    </div>
    <div id="cm-exon"    class="ctab-content active"><div id="cm-circle-wrap"></div></div>
    <div id="cm-mirna"   class="ctab-content"></div>
    <div id="cm-rbp"     class="ctab-content"></div>
    <div id="cm-volcano" class="ctab-content" style="min-height:340px"></div>
    <div id="cm-heatmap" class="ctab-content" style="min-height:340px"></div>
    <div id="cm-dl-bar" style="text-align:right;margin-top:10px;border-top:1px solid #eee;padding-top:8px"></div>
  </div>
</div>"""

    # Build self-contained MultiQC section (srcdoc).
    # Fix 1: inject nav-guard so empty/relative hrefs (MultiQC logo, etc.) cannot navigate
    #         the iframe to the parent report URL.
    # Fix 2: on toggle-open, reflow Highcharts at 200/700/1500ms to catch all charts that
    #         initialised with 0-width containers while <details> was collapsed.
    if multiqc_file and os.path.exists(multiqc_file):
        with open(multiqc_file, encoding="utf-8", errors="replace") as _mf:
            _mqc_raw = _mf.read()
        _nav_guard = (
            '<script>'
            '(function(){'
            'document.addEventListener("click",function(e){'
            'var a=e.target.closest("a[href]");'
            'if(!a)return;'
            'var h=a.getAttribute("href")||"";'
            'if(h.startsWith("#"))return;'
            'e.preventDefault();'
            '},true);'
            '})();'
            '</script>'
        )
        _last_body = _mqc_raw.rfind('</body>')
        if _last_body >= 0:
            _mqc_raw = _mqc_raw[:_last_body] + _nav_guard + _mqc_raw[_last_body:]
        else:
            _mqc_raw += _nav_guard
        # Embed MultiQC as data:text/html;base64,... directly in the iframe src attribute.
        # This is the only approach that avoids all four known issues in file:// context:
        #   1. srcdoc >~1MB silently loads parent page (Chrome/Edge bug)
        #   2. <script type="application/json"> closed early by </script> inside MultiQC
        #   3. blob:null iframe #anchor clicks resolve to parent file URL (Edge bug)
        #   4. document.write() context confusion makes MultiQC JS execute in parent frame
        # Base64 alphabet (A-Za-z0-9+/=) has no HTML-special chars, so the attribute is safe.
        # Anchor links (href="#section") inside the data: iframe trigger same-document
        # fragment navigation — they scroll within the iframe, never navigate the parent.
        _mqc_b64 = base64.b64encode(_mqc_raw.encode('utf-8')).decode('ascii')
        multiqc_section = (
            '<details id="qc-section" open style="margin:16px 0 24px 0">\n'
            '  <summary style="cursor:pointer;padding:10px 16px;background:#f0f9ff;'
            'border:1px solid #bae6fd;border-radius:8px;font-size:15px;font-weight:600;'
            'color:#0369a1;list-style:none;display:flex;align-items:center;gap:8px">\n'
            '    <span>▼</span> \U0001f4ca QC Report (MultiQC)\n'
            '    <span id="qc-collapse-lbl" data-en="(click to collapse)" '
            'style="font-size:12px;font-weight:400;color:#64748b;margin-left:4px">'
            '（點擊折疊）</span>\n'
            '  </summary>\n'
            '  <div style="margin-top:8px">\n'
            '    <iframe id="qc-iframe"\n'
            '      src="data:text/html;base64,' + _mqc_b64 + '"\n'
            '      style="width:100%;height:850px;border:1px solid #e2e8f0;border-radius:6px;background:#fff"\n'
            '      title="MultiQC Report"></iframe>\n'
            '  </div>\n'
            '</details>\n'
            '<script>\n'
            '(function(){\n'
            '  var det=document.getElementById("qc-section");\n'
            '  if(det){\n'
            '    det.addEventListener("toggle",function(){\n'
            '      var sp=det.querySelector("summary span:first-child");\n'
            '      var lb=document.getElementById("qc-collapse-lbl");\n'
            '      if(sp) sp.textContent=this.open?"▼":"▶";\n'
            '      if(lb){\n'
            '        var isEn=typeof _LANG!=="undefined"&&_LANG==="en";\n'
            '        var openZh="（點擊折疊）",closeZh="（點擊展開）";\n'
            '        var openEn="(click to collapse)",closeEn="(click to expand)";\n'
            '        lb.dataset.en=this.open?openEn:closeEn;\n'
            '        lb.dataset.zh=this.open?openZh:closeZh;\n'
            '        lb.textContent=this.open?(isEn?openEn:openZh):(isEn?closeEn:closeZh);\n'
            '      }\n'
            '      if(this.open){\n'
            '        var fr=document.getElementById("qc-iframe");\n'
            '        function _r(){\n'
            '          try{\n'
            '            var hc=fr.contentWindow.Highcharts;\n'
            '            if(hc&&hc.charts) hc.charts.forEach(function(c){if(c)c.reflow();});\n'
            '          }catch(e){}\n'
            '        }\n'
            '        setTimeout(_r,300); setTimeout(_r,800);\n'
            '      }\n'
            '    });\n'
            '  }\n'
            '})();\n'
            '</script>'
        )
    else:
        multiqc_section = ""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>circRNA Analysis Report – {project_id}</title>
  {_STYLE}
  {_SCRIPT}
  <script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
</head>
<body>

<div class="print-bar no-print">
  <!-- CircDEX brand -->
  <div class="cd-rpt-brand">
    <svg width="52" height="52" viewBox="0 0 38 38" fill="none" aria-hidden="true">
      <circle cx="19" cy="19" r="15" stroke="rgba(0,180,198,.18)" stroke-width="1.5"/>
      <path d="M19 4 A15 15 0 1 1 7.3 28" stroke="#00B4C6" stroke-width="2.6"
            stroke-linecap="round" fill="none"/>
      <circle cx="19" cy="4"   r="2.8" fill="#00B4C6"/>
      <circle cx="7.3" cy="28" r="1.8" fill="#00B4C6" opacity=".55"/>
    </svg>
    <div class="cd-rpt-meta">
      <div class="cd-rpt-wm"><span class="cd-rpt-circ">Circ</span><span class="cd-rpt-dex">DEX</span></div>
      <div class="cd-rpt-sub">From reads to circRNA biomarkers</div>
      <div class="cd-rpt-chips">
        <span class="cd-rpt-chip">Dual-tool consensus</span>
        <span class="cd-rpt-chip">Differential expression</span>
        <span class="cd-rpt-chip">6D biomarker ranking</span>
      </div>
    </div>
  </div>
  <!-- Nav right -->
  <div class="cd-rpt-nav">
    <span class="cd-rpt-proj">{project_id}</span>
    <span class="method-tag" style="font-size:12px">{de_method}</span>
    <span style="font-size:11px;color:rgba(255,255,255,.35);padding:0 4px">{datetime.now().strftime('%Y-%m-%d')}</span>
    <button class="print-btn" data-en="🖨 Print / PDF" onclick="window.print()">🖨 列印 / PDF</button>
    <div style="display:flex;gap:3px;margin-left:4px">
      <button class="lang-btn-rpt active" data-lang="zh" onclick="switchReportLang('zh')">中文</button>
      <button class="lang-btn-rpt" data-lang="en" onclick="switchReportLang('en')">EN</button>
    </div>
  </div>
</div>
  {"<p style='font-size:12px;color:#888;margin-top:20px'>Interaction data pre-fetched for "
   + str(n_ixn) + " circRNAs (" + str(n_ixn_mirna) + " with miRNA data, from CircInteractome). "
   + "Click any circ_position to view exon diagram, miRNA sponge sites, and RBP binding sites.</p>" if n_ixn > 0
   else "<p style='font-size:12px;color:#aaa;margin-top:20px'>Click any circ_position to view exon diagram "
        "(interaction data not yet fetched — run predict_interactions rule).</p>"}
  <p style="font-size:11px;color:#888">
    <button onclick="showCircDetail('{list(interactions.keys())[0] if interactions else ''}')"
      style="background:#2c6fad;color:#fff;border:none;border-radius:4px;padding:3px 10px;cursor:pointer;font-size:11px">
      &#128269; Test: Open first circRNA detail
    </button>
  </p>

  {sample_html}

  {multiqc_section}

  <!-- ① Summary stat boxes + DE method switcher -->
  <h2 data-en="Summary">摘要</h2>
  {_msw_html}
  <div>
    <div class="stat-box"><div class="num">{n_sample}</div><div class="lbl">Samples</div></div>
    <div class="stat-box"><div class="num">{n_total}</div><div class="lbl">Total circRNAs</div></div>
    <div class="stat-box"><div class="num" id="stat-n-sig">{n_sig}</div><div class="lbl">Significant ({sig_label}, |log2FC|&gt;{lfc})</div></div>
    <div class="stat-box"><div class="num" id="stat-n-up">{n_up}</div><div class="lbl">Up-regulated</div></div>
    <div class="stat-box"><div class="num" id="stat-n-dn">{n_dn}</div><div class="lbl">Down-regulated</div></div>
  </div>

  <!-- ② PCA — sample quality and clustering -->
  <h2 data-en="PCA">PCA（樣本群聚分析）</h2>
  {pca_html}

  <!-- ③ Volcano — overall DE landscape -->
  <h2 data-en="Volcano Plot">Volcano Plot（全局差異表現）</h2>
  <p style="font-size:12px;color:#888">&#9711; Heatmap top {heatmap_top_n} up + {heatmap_top_n} down markers: use the toggle button in the chart to show/hide.</p>
  {volcano_html}

  <!-- ── Clustering Heatmap (hierarchical row clustering, all sig circRNAs) ── -->
  <div id="clust-heatmap-section" style="margin-top:32px">
    <h2 style="font-size:18px;font-weight:600;color:#2c3e50;
               padding:8px 0;border-bottom:2px solid #e0e8f0"
        data-en="Clustering Heatmap (all significant DE circRNAs, hierarchical row clustering)">
      聚類熱圖（全部顯著 DE circRNA，階層式 row 聚類）
    </h2>
    <p style="font-size:12px;color:#888;margin:6px 0 2px"
       data-en="Shows all significant DE circRNAs (primary method: {de_method}) ordered by hierarchical clustering (Ward linkage). Use Plotly zoom tools or click-drag to inspect regions of interest. Double-click to reset zoom.">
      顯示全部顯著 DE circRNA（主方法：{de_method}），依階層聚類（Ward linkage）排列。可用 Plotly 工具列或拖曳放大感興趣區域，雙擊還原。
    </p>
    <div id="clust-heatmap-plot" style="width:100%"></div>
  </div>

  <!-- ④ Top DE tables — up / down circRNAs -->
  <div id="de-tables-section">
  <h2 id="de-tables-heading">Top Differentially Expressed circRNAs ({sig_label}, |log2FC| &gt; {lfc})</h2>
  <p style="font-size:12px;color:#666">&#128204; Click a <strong>circ_position</strong> to view exon diagram, miRNA sponge sites, and RBP binding sites.</p>
  {_de_split_tables(top_table, tumor_label=tumor_label, normal_label=normal_label, interactions=interactions)}
  </div>

  <!-- ⑤ Type I / II proportion -->
  {type_html}

  <!-- ⑥ 3-method Venn — method agreement -->
  {"<h2 data-en='3-Method DE Venn Diagram'>三方法 DE 結果 Venn Diagram</h2><p style='font-size:13px;color:#555' data-en='Compares significant DE circRNAs across three methods (edgeR FSJ offset, DESeq2, limma-voom) at the same threshold.'>比較三種方法（edgeR FSJ offset、DESeq2、limma-voom）在相同閾值下的顯著 DE circRNA 交集。</p>" + venn_html if venn_html else ""}

  <!-- ⑦ Biomarker candidates — final ranked list -->
  <div id="biomarker-section">
  {biomarker_html}
  </div>

  <!-- ⑧ Isoform switching -->
  <div id="isoform-section">
  {isoform_html}
  <p style="font-size:11px;color:#999;margin:-4px 0 8px" data-en="&#8505; Isoform switching is based on IUI; not affected by DE method switching.">&#8505; Isoform switching 依據 IUI 計算，不受 DE 方法切換影響。</p>
  </div>

{_modal_html}
{_modal_js}
</body>
</html>
"""
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(output_file).write_text(html, encoding="utf-8")
    print(f"[OK] Report written → {output_file}")


# ── Snakemake entry point ────────────────────────────────────────────────────

if "snakemake" in dir():
    build_report(
        de_file        = snakemake.input.de,                    # type: ignore[name-defined]
        matrix_file    = snakemake.input.matrix,                # type: ignore[name-defined]
        volcano_pdf    = snakemake.input.volcano,               # type: ignore[name-defined]
        heatmap_pdf    = snakemake.input.heatmap,               # type: ignore[name-defined]
        pca_pdf        = snakemake.input.pca,                   # type: ignore[name-defined]
        output_file    = snakemake.output[0],                   # type: ignore[name-defined]
        project_id     = snakemake.params.project_id,           # type: ignore[name-defined]
        fdr            = float(snakemake.params.fdr),           # type: ignore[name-defined]
        lfc            = float(snakemake.params.lfc),           # type: ignore[name-defined]
        de_method      = str(snakemake.params.de_method),       # type: ignore[name-defined]
        biomarker_file = snakemake.input.biomarkers,            # type: ignore[name-defined]
        switching_file = snakemake.input.switching,             # type: ignore[name-defined]
        groups_file    = getattr(snakemake.input, "groups", None),  # type: ignore[name-defined]
        de_sig_by      = str(getattr(snakemake.params, "de_sig_by", "auto")),   # type: ignore[name-defined]
        tumor_label    = str(snakemake.params.tumor_label),     # type: ignore[name-defined]
        normal_label   = str(snakemake.params.normal_label),    # type: ignore[name-defined]
        isoform_file   = getattr(snakemake.input, "isoform_groups", None),  # type: ignore[name-defined]
        circbase_file  = getattr(snakemake.input, "circbase_annot", None),  # type: ignore[name-defined]
        heatmap_top_n      = int(getattr(snakemake.params, "heatmap_top_n", 10)),  # type: ignore[name-defined]
        interactions_file  = getattr(snakemake.input, "interactions", None),       # type: ignore[name-defined]
        multiqc_file       = getattr(snakemake.input, "multiqc", None),            # type: ignore[name-defined]
        de_files           = {m: getattr(snakemake.input, a, None)                 # type: ignore[name-defined]
                              for m, a in [("edgeR_ciriquant","de_edger"),
                                           ("deseq2","de_deseq"),
                                           ("limma","de_limma")]
                              if getattr(snakemake.input, a, None)},               # type: ignore[name-defined]
        study_title        = str(getattr(snakemake.params, "study_title", "")),    # type: ignore[name-defined]
    )
