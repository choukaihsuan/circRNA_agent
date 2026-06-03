"""
generate_report.py – Build a self-contained HTML summary report.

Called as a Snakemake script.  Embeds PDF plots as base64-encoded images
and includes top DE results as a table.
"""

import base64
import math
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


def _df_to_html(df: pd.DataFrame, max_rows: int = 50) -> str:
    return (
        _fmt_floats(df.head(max_rows))
          .to_html(index=False, classes="table", border=0, na_rep="—")
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


def _venn_3_svg(sig_sets: dict, lfc: float) -> str:
    """Generate a 3-circle (or 2-circle) Venn SVG for DE method overlap."""
    methods = [k for k in ("edgeR_ciriquant", "deseq2", "limma") if k in sig_sets]
    if len(methods) < 2:
        return ""
    sets = [sig_sets[m] for m in methods]
    # 7-region computation
    A, B = sets[0], sets[1]
    C = sets[2] if len(sets) > 2 else set()
    ab = A & B; ac = A & C; bc = B & C; abc = A & B & C
    counts = {
        "A":   len(A - B - C),  "B":  len(B - A - C),  "C":   len(C - A - B),
        "AB":  len(ab - C),     "AC": len(ac - B),      "BC":  len(bc - A),
        "ABC": len(abc),
    }
    total_union = len(A | B | C) if len(sets) > 2 else len(A | B)
    # SVG layout (3-circle)
    W, H = 460, 295
    cx = [150, 310, 230]; cy = [128, 128, 218]; r = 90
    colors = ["rgba(214,39,40,0.20)", "rgba(44,160,44,0.20)", "rgba(44,119,214,0.20)"]
    strokes = ["#d62728", "#2CA02C", "#2c6fad"]
    labels  = {"edgeR_ciriquant": "edgeR (FSJ offset)", "deseq2": "DESeq2", "limma": "limma-voom"}
    # Label anchor positions (outside each circle)
    lpos = [(cx[0]-r-2, cy[0]-r+5), (cx[1]+r+2, cy[1]-r+5), (cx[2], cy[2]+r+14)]
    # Region count positions
    rpos = {
        "A":  (cx[0]-46, cy[0]),     "B":  (cx[1]+46, cy[1]),    "C":  (cx[2], cy[2]+46),
        "AB": (230, cy[0]-20),        "AC": (165, cy[0]+62),      "BC": (295, cy[1]+62),
        "ABC":(230, 178),
    }
    svg_parts = [f'<svg width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
                 f'xmlns="http://www.w3.org/2000/svg" style="font-family:sans-serif;max-width:480px">']
    # Circles
    for i in range(min(len(methods), 3)):
        svg_parts.append(f'<circle cx="{cx[i]}" cy="{cy[i]}" r="{r}" '
                         f'fill="{colors[i]}" stroke="{strokes[i]}" stroke-width="2.5"/>')
    # Method labels
    for i, m in enumerate(methods[:3]):
        lx, ly = lpos[i]
        anc = "end" if i == 0 else ("start" if i == 1 else "middle")
        svg_parts.append(f'<text x="{lx}" y="{ly}" font-size="11" font-weight="bold" '
                         f'fill="{strokes[i]}" text-anchor="{anc}">'
                         f'{labels.get(m, m)}</text>')
        svg_parts.append(f'<text x="{lx}" y="{ly+14}" font-size="10" '
                         f'fill="{strokes[i]}" text-anchor="{anc}">'
                         f'(n={len(sig_sets[m])})</text>')
    # Region counts
    for region, (rx, ry) in rpos.items():
        if len(methods) < 3 and region in ("C", "AC", "BC", "ABC"):
            continue
        v = counts[region]
        if v == 0:
            continue
        svg_parts.append(f'<text x="{rx}" y="{ry}" font-size="13" font-weight="bold" '
                         f'text-anchor="middle" dominant-baseline="middle" fill="#333">{v}</text>')
    svg_parts.append('</svg>')
    note = (f'<p style="font-size:12px;color:#666;margin:4px 0 0">数字 = 各區域顯著 circRNA 數量'
            f'（閾值：|log₂FC| &gt; {lfc}）；三方法聯集共 {total_union} 個。</p>')
    return f'<div style="text-align:center">{"".join(svg_parts)}</div>{note}'


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
                up_p = de[de["log2FC"] > 0].dropna(subset=[p_col]).sort_values(p_col).head(hm_pool)
                dn_p = de[de["log2FC"] < 0].dropna(subset=[p_col]).sort_values(p_col).head(hm_pool)
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
    return {"stats": {"n_sig": n_sig, "n_up": n_up, "n_dn": n_dn},
            "volcano": vol_rows, "heatmap": hm_data}


def _load_interactions(json_file: Optional[str]) -> dict:
    """Load pre-fetched ENCORI interactions JSON, return {} on any error."""
    if not json_file:
        return {}
    try:
        import json as _json
        return _json.loads(Path(json_file).read_text(encoding="utf-8"))
    except Exception:
        return {}


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
  .table th { background: #2c6fad; color: #fff; padding: 6px 10px; text-align: left; }
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
  .itable th { background:#f0f5ff; color:#335; padding:5px 10px; border:1px solid #dde; text-align:left; }
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
  .print-bar  { position:sticky; top:0; z-index:99; background:#2c6fad; color:#fff;
                padding:8px 24px; display:flex; align-items:center; gap:12px;
                flex-wrap:wrap; }
  .print-bar span { font-size:13px; font-weight:bold; flex:1; }
  .print-btn  { padding:5px 14px; border:2px solid #fff; border-radius:4px;
                background:transparent; color:#fff; font-size:13px;
                font-weight:bold; cursor:pointer; }
  .print-btn:hover { background:rgba(255,255,255,.15); }

  /* DE method switcher */
  .msw-btn { padding:5px 14px; border:1px solid #2c6fad; border-radius:4px;
             background:white; color:#2c6fad; cursor:pointer; font-size:13px;
             margin-right:6px; transition:all 0.15s; }
  .msw-btn.active { background:#2c6fad; color:white; }
  .msw-btn:hover:not(.active) { background:#eaf3ff; }

  @media print {
    .print-bar, .dl-btn, .no-print { display:none !important; }
    body { max-width:100%; margin:0; padding:0 12px; }
    .table th { background:#2c6fad !important; -webkit-print-color-adjust:exact;
                print-color-adjust:exact; }
    .table { page-break-inside:avoid; font-size:11px; }
    h1, h2 { color:#2c6fad !important; -webkit-print-color-adjust:exact;
             print-color-adjust:exact; }
  }
</style>
"""

_SCRIPT = """
<script>
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


def _plot_isoform_usage(sig: pd.DataFrame, top_n: int = 10) -> str:
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
                x           = [f"{gene}<br>Normal", f"{gene}<br>Tumor"],
                y           = [row.get("iui_normal", 0), row.get("iui_tumor", 0)],
                marker_color= colors[int(i) % len(colors)],
                legendgroup = gene,
            ))

    fig.update_layout(
        barmode       = "stack",
        title         = "Isoform Usage Index — Tumor vs Normal (significant switching only)",
        yaxis_title   = "Isoform Usage Index (IUI)",
        xaxis_tickangle = -30,
        height        = 480,
        legend_title  = "Gene | BSJ coords",
        plot_bgcolor  = "white",
        paper_bgcolor = "white",
    )
    return fig.to_html(include_plotlyjs="cdn", full_html=False)


def _isoform_section(switching_file: Optional[str],
                     isoform_file:  Optional[str] = None,
                     circbase_file: Optional[str] = None) -> str:
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

    plot_html = _plot_isoform_usage(sig)

    show_cols = [c for c in
                 ["gene_name", "circ_id", "exon_structure", "circbase_id",
                  "iui_normal", "iui_tumor", "delta_iui", "p_value"]
                 if c in sig.columns]
    sort_keys = [c for c in ["gene_name", "p_value"] if c in sig.columns]
    table_html = (
        _dl_wrap(_df_to_html(sig.sort_values(sort_keys)[show_cols], max_rows=40),
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


def _biomarker_section(biomarker_file: Optional[str],
                       interactions: Optional[dict] = None) -> str:
    if not biomarker_file or not Path(biomarker_file).exists():
        return ""
    try:
        bm = pd.read_csv(biomarker_file, sep="\t")
    except Exception:
        return ""
    if bm.empty:
        return ""
    show_cols = [c for c in ["rank", "circ_id", "log2FC", "n_mirna", "n_rbp",
                              "biomarker_score", "in_circbase", "circbase_id",
                              "circbase_gene", "Type"]
                 if c in bm.columns]
    has_interactions = "n_mirna" in bm.columns and (bm["n_mirna"].max() + bm["n_rbp"].max()) > 0
    if has_interactions:
        score_desc = "Score = mean of: −log₁₀(p-value), |log₂FC|, confidence score, circBase bonus, #miRNA binders, #RBP binders (each normalised 0–1; 6D)."
    else:
        score_desc = "Score = mean of: −log₁₀(p-value), |log₂FC|, confidence score, circBase known bonus (each normalised 0–1; 4D)."

    disp = _fmt_floats(bm[show_cols].head(30).copy())
    if "circ_id" in disp.columns and interactions is not None:
        def _link(v: str) -> str:
            tip = "in interaction data" if str(v) in interactions else "no interaction data pre-fetched"
            return (f'<a class="circ-link" onclick="showCircDetail(\'{v}\')" '
                    f'title="{tip}">{v}</a>')
        disp["circ_id"] = disp["circ_id"].astype(str).apply(_link)

    table_html = _dl_wrap(
        disp.to_html(index=False, classes="table", border=0, na_rep="—", escape=False),
        "tbl_biomarker", "biomarker_candidates.csv"
    )

    return f"""
  <h2>Biomarker Candidates (top {min(len(bm), 30)} by composite score)</h2>
  <p style="font-size:13px;color:#555;margin-bottom:8px;">{score_desc}</p>
  <p style="font-size:12px;color:#666">&#128204; Click a <strong>circ_id</strong> to view exon diagram, miRNA and RBP binding sites.</p>
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
    return fig.to_html(include_plotlyjs="cdn", full_html=False, div_id=div_id)


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
    if groups_file and Path(groups_file).exists():
        try:
            grp = pd.read_csv(groups_file)
            cond_map = dict(zip(grp["srr_id"], grp["condition"]))
            normal_cols = [c for c in log_sub.columns if cond_map.get(c, "") == normal_label]
        except Exception:
            pass
    ref_cols = normal_cols if normal_cols else log_sub.columns.tolist()
    row_mean = log_sub[ref_cols].mean(axis=1)       # center on normal mean
    row_std  = log_sub.std(axis=1).clip(lower=0.1)  # scale on all-sample SD
    z = log_sub.sub(row_mean, axis=0).div(row_std, axis=0)

    n_up = sum(i in up_ids for i in avail)
    n_dn = sum(i in dn_ids for i in avail)
    y_labels = [_label(i) for i in avail]

    fig = go.Figure(go.Heatmap(
        z=z.values.tolist(), x=z.columns.tolist(), y=y_labels,
        colorscale=[[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']], zmid=0,
        colorbar=dict(title="z-score<br>(normal-<br>centered)"),
        hovertemplate="<b>%{y}</b><br>%{x}<br>z-score: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=f"Top {n_up} up + {n_dn} down DE circRNAs — Heatmap (normal-centered z-score)",
                   font_size=14),
        yaxis=dict(tickfont=dict(size=8), autorange="reversed"),
        height=max(420, len(avail) * 22 + 120),
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(t=60, l=300),
    )
    return fig.to_html(include_plotlyjs="cdn", full_html=False, div_id="main-heatmap-plot")


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
    color_map  = {"tumor": "#d62728", "normal": "#2CA02C", "unknown": "#888"}
    fallback   = ["#2563eb", "#e07b39", "#16a34a", "#9333ea", "#dc2626"]
    for i, c in enumerate(sorted(set(conditions))):
        color_map.setdefault(c, fallback[i % len(fallback)])

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
    return fig.to_html(include_plotlyjs="cdn", full_html=False)


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
                _up_pool = de[de["log2FC"] > 0].dropna(subset=[_pc_hm]).sort_values(_pc_hm).head(_HM_POOL)
                _dn_pool = de[de["log2FC"] < 0].dropna(subset=[_pc_hm]).sort_values(_pc_hm).head(_HM_POOL)
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
        except Exception as _exc:
            import sys as _sys
            print(f"[report] DE file for {_mkey} failed: {_exc}", file=_sys.stderr)

    # Ensure primary method is in venn sets
    if de_method not in sig_sets_venn:
        sig_sets_venn[de_method] = _sig_ids_from_de(de, p_col, sig_thr, lfc)

    all_de_methods_js = _json.dumps(all_de_data, ensure_ascii=False)
    venn_html = _venn_3_svg(sig_sets_venn, lfc) if len(sig_sets_venn) >= 2 else ""

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
        '<span style="font-size:13px;color:#555;margin-right:6px">分析方法：</span>'
        + "".join(_msw_btns)
        + '<span style="font-size:11px;color:#999;margin-left:8px">切換後 Volcano / Heatmap / 統計即時更新</span>'
        + '</div>'
    ) if len(_msw_btns) > 1 else ""

    type_html      = _type_section(sig)
    biomarker_html = _biomarker_section(biomarker_file, interactions=interactions)
    isoform_html   = _isoform_section(switching_file,
                                       isoform_file=isoform_file,
                                       circbase_file=circbase_file)

    # Interactive Plotly charts; fall back to static PDF embeds when unavailable
    p_volcano = _plotly_volcano(de, fdr, lfc, de_method, p_col=p_col, sig_thr=sig_thr,
                                p_label=_sig_label_str, heatmap_ids=heatmap_ids)
    p_heatmap = _plotly_heatmap(de, matrix, top_n=heatmap_top_n, p_col=p_col,
                                groups_file=groups_file, normal_label=normal_label)
    p_pca     = _plotly_pca(matrix, groups_file)
    volcano_html = p_volcano if p_volcano else _embed_pdf(volcano_pdf)
    heatmap_html = p_heatmap if p_heatmap else _embed_pdf(heatmap_pdf)
    pca_html     = p_pca     if p_pca     else _embed_pdf(pca_pdf)

    n_ixn = len(interactions)
    n_ixn_mirna = sum(len(v.get("mirna", [])) > 0 for v in interactions.values())

    _modal_js = f"""
<script>
const CIRC_DATA         = {interactions_js};
const VOLCANO_DATA      = {volcano_data_js};
const FULL_HEATMAP_DATA  = {full_heatmap_data_js};
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
      [info.gene_name, info.strand, info.region, info.exon_span].filter(Boolean).join(' &nbsp;|&nbsp; ');
    const circEl = document.getElementById('cm-circle-wrap');
    circEl.innerHTML = '';
    _drawCircleRNA(circId, circEl);
    document.getElementById('cm-mirna').innerHTML = _buildInteractionTable(
      d.mirna||[], ['_priority','miRNAName','siteType','circ_pos','_seq_logo','clipExpNum','cellType','source','in_circ'],
      ['Priority','miRNA','Site Type','Chr Position','Binding Seq','CLIP Exp.','Cell Type','Source','In circRNA'], circId, 'mirna');
    document.getElementById('cm-rbp').innerHTML = _buildInteractionTable(
      d.rbp||[], ['_priority','RBPName','bindingSites','circ_pos','_seq_logo','location','clipExpNum','cellType','source','in_circ'],
      ['Priority','RBP','Sites','Chr Position','Binding Seq','Location','CLIP Exp.','Cell Type','Source','In circRNA'], circId, 'rbp');
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
  const blob = new Blob([lines.join('\n')], {{type:'text/csv;charset=utf-8;'}});
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
  const totalLen = parseInt(info.spliced_length) || 0;
  const exonBds  = info.exon_boundaries || [];

  const W=480, H=440, cx=240, cy=210;
  const ROUT=148, RIN=115, MI_OUT=176, MI_IN=153, RBP_OUT=113, RBP_IN=90;

  const MI_COLORS  = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628',
                      '#c4007a','#0077b6','#c96a00','#2d6a4f','#6d4c91','#b5470b',
                      '#1d7a3a','#7b2d8b','#0059a6','#8b6900','#c1121f','#3d7a00',
                      '#1d91c0','#6a3d9a'];
  const RBP_COLORS = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e',
                      '#e6ab02','#a6761d','#333333','#1f78b4','#b2df8a','#fb8072','#80b1d3'];

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

  // ── Proportional exon boundaries ──
  if(exonBds.length>0&&totalLen>0){{
    exonBds.forEach(eb=>{{
      // divider line at exon junction
      if(eb.cum_start>0){{
        const a=posToAngle(eb.cum_start);
        const [lx1,ly1]=polar(RIN-2,a),[lx2,ly2]=polar(ROUT+2,a);
        svg+=`<line x1="${{lx1.toFixed(1)}}" y1="${{ly1.toFixed(1)}}" x2="${{lx2.toFixed(1)}}" y2="${{ly2.toFixed(1)}}" stroke="#999" stroke-width="1.5"/>`;
      }}
      // exon label on the middle ring (between RIN and ROUT)
      const arcFrac=(eb.cum_end-eb.cum_start)/totalLen;
      const arcPx=arcFrac*2*Math.PI*(RIN+ROUT)/2;
      if(arcPx>30){{  // only label if enough space
        const aMid=posToAngle((eb.cum_start+eb.cum_end)/2);
        const [lx,ly]=polar((RIN+ROUT)/2,aMid);
        const isLow=(aMid>0&&aMid<Math.PI);
        const tdeg=isLow?(aMid*180/Math.PI-90):(aMid*180/Math.PI+90);
        const exonLabel=eb.label.replace(/^e(\\d+)$/,'exon $1');
        svg+=`<text transform="translate(${{lx.toFixed(1)}},${{ly.toFixed(1)}}) rotate(${{tdeg.toFixed(1)}})" text-anchor="middle" dominant-baseline="central" font-size="11" fill="#444" font-weight="600">${{exonLabel}}</text>`;
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
  const miMap={{}};let miIdx=0;
  const miLegend=[];
  const miBadgeAngs={{}};
  const miArcGroups={{}};
  if(totalLen>0){{
    mirnaList.forEach(item=>{{
      const name=item.miRNAName||'';
      const m=(item.circ_pos||'').match(/(\\d+)[–-](\\d+)/);
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
      miArcGroups[name].push({{d:arcPath(MI_OUT,MI_IN,a1,a2),title:`${{name}} · ${{item.siteType||''}} · ${{item.circ_pos}}`}});
      if(!(name in miBadgeAngs)) miBadgeAngs[name]=normA((a1+a2)/2);
    }});
    // Draw arc groups (each name wrapped in <g> for show/hide)
    Object.entries(miArcGroups).forEach(([name,arcs])=>{{
      const e=miMap[name];
      svg+=`<g id="_miarc_${{e.num}}" class="_miarc">`;
      arcs.forEach(a=>svg+=`<path d="${{a.d}}" fill="${{e.color}}" opacity="0.85"><title>${{a.title}}</title></path>`);
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
  }}

  // ── RBP arcs (inner ring) — draw arcs first, then de-overlapped badges ──
  const rbpMap={{}};let rbpIdx=0;
  const rbpLegend=[];
  const rbpBadgeAngs={{}};
  const rbpArcGroups={{}};
  if(totalLen>0){{
    rbpList.forEach(item=>{{
      if(item.location!=='internal')return;
      const name=item.RBPName||'';
      if(!(name in rbpMap)){{
        rbpIdx++;
        rbpMap[name]={{num:rbpIdx,color:RBP_COLORS[(rbpIdx-1)%RBP_COLORS.length]}};
        rbpLegend.push({{num:rbpIdx,name,color:rbpMap[name].color,ns:(item.sites||[]).length,src:item.source||''}});
        rbpArcGroups[name]=[];
      }}
      const e=rbpMap[name];
      const sites=item.sites||[];
      if(sites.length>0){{
        sites.forEach(s=>{{
          const a1=posToAngle(s.circ_start),a2=posToAngle(Math.max(s.circ_end,s.circ_start+1));
          rbpArcGroups[name].push({{d:arcPath(RBP_OUT,RBP_IN,a1,a2),title:`${{name}} · ${{s.circ_pos}}`}});
          if(!(name in rbpBadgeAngs)) rbpBadgeAngs[name]=normA((a1+a2)/2);
        }});
      }} else {{
        if(!(name in rbpBadgeAngs)){{
          const fallbackAng=normA(-Math.PI/2+(rbpIdx-1)/(Math.max(rbpList.filter(r=>r.location==='internal').length,1))*2*Math.PI);
          rbpBadgeAngs[name]=fallbackAng;
        }}
      }}
    }});
    // Draw RBP arc groups
    Object.entries(rbpArcGroups).forEach(([name,arcs])=>{{
      const e=rbpMap[name];
      svg+=`<g id="_rbparc_${{e.num}}" class="_rbparc">`;
      arcs.forEach(a=>svg+=`<path d="${{a.d}}" fill="${{e.color}}" opacity="0.8"><title>${{a.title}}</title></path>`);
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
  }}

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
  const gname=info.gene_name?`${{info.gene_name}}${{info.strand?' ('+info.strand+')':''}}`:'';
  const eline=info.exon_span?`exon ${{info.exon_span}}`:(info.region||'');
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
      ${{toggleBtn('全顯示','_toggleAll(\\'mi\\',true,this.closest(\\'div\\'))')}}
      ${{toggleBtn('全隱藏','_toggleAll(\\'mi\\',false,this.closest(\\'div\\'))')}}
      <div style="display:flex;flex-wrap:wrap;gap:3px 8px;margin-top:4px">`;
    miLegend.forEach(e=>{{
      legHtml+=`<span id="_leg_mi_${{e.num}}" onclick="_toggleBadge('mi',${{e.num}},this)"
        title="點擊顯示/隱藏" style="display:inline-flex;align-items:center;gap:2px;cursor:pointer;padding:1px 3px;border-radius:3px;border:1px solid #eee">
        ${{badge(e.color,e.num)}}${{e.name}}${{e.src?srcBadge(e.src):''}}</span>`;
    }});
    legHtml+='</div></div>';
  }}
  if(rbpLegend.length>0){{
    legHtml+=`<div style="margin-top:6px;font-size:10px">
      <span style="font-weight:bold">RBP (inner ring):</span>
      ${{toggleBtn('全顯示','_toggleAll(\\'rbp\\',true,this.closest(\\'div\\'))')}}
      ${{toggleBtn('全隱藏','_toggleAll(\\'rbp\\',false,this.closest(\\'div\\'))')}}
      <div style="display:flex;flex-wrap:wrap;gap:3px 8px;margin-top:4px">`;
    rbpLegend.forEach(e=>{{
      legHtml+=`<span id="_leg_rbp_${{e.num}}" onclick="_toggleBadge('rbp',${{e.num}},this)"
        title="點擊顯示/隱藏" style="display:inline-flex;align-items:center;gap:2px;cursor:pointer;padding:1px 3px;border-radius:3px;border:1px solid #eee">
        ${{badge(e.color,e.num)}}${{e.name}} (${{e.ns}} sites)${{e.src?srcBadge(e.src):''}}</span>`;
    }});
    legHtml+='</div></div>';
  }}
  container.innerHTML=`<div style="text-align:center">${{svg}}</div>`+legHtml;
}}

// ── Mini heatmap ──────────────────────────────────────────────────────────────
function _buildMiniHeatmap(circId) {{
  const el=document.getElementById('cm-heatmap');
  if(!el)return;
  if(typeof Plotly==='undefined'){{el.innerHTML='<p class="no-data">Plotly not available.</p>';return;}}
  const _hmd=_HEATMAP_DATA_CACHE||FULL_HEATMAP_DATA;
  if(!_hmd||!_hmd.rows||Object.keys(_hmd.rows).length===0){{
    el.innerHTML='<p class="no-data">Heatmap data not available.</p>';return;
  }}
  const _inTop=_hmd.rows[circId]!=null;
  const allIds=Object.keys(_hmd.rows);
  const tIdx=allIds.indexOf(circId);
  const samps=_hmd.samples||[];
  const zMatrix=allIds.map(id=>_hmd.rows[id].z);
  const yLabels=allIds.map(id=>_hmd.rows[id].label||id);
  const conds=_hmd.conditions||{{}};
  const TUMOR_COL='#d62728', NORMAL_COL='#2c6fad';
  // merge consecutive samples of same condition into one labelled bar
  const grps=[];let cur=null;
  samps.forEach((s,i)=>{{
    const c=conds[s]||'';
    const col=c==='{tumor_label}'?TUMOR_COL:c==='{normal_label}'?NORMAL_COL:'#888';
    if(!cur||cur.c!==c){{cur={{c,col,s:i,e:i}};grps.push(cur);}}
    else cur.e=i;
  }});
  const groupShapes=grps.map(g=>{{
    return{{type:'rect',xref:'x',yref:'paper',
      x0:g.s-0.45,x1:g.e+0.45,y0:1.02,y1:1.055,
      fillcolor:g.col,line:{{width:0}}}};
  }});
  const groupAnno=grps.map(g=>{{
    return{{xref:'x',yref:'paper',x:(g.s+g.e)/2,y:1.0375,
      yanchor:'middle',xanchor:'center',
      text:g.c,showarrow:false,
      font:{{size:9,color:'white',family:'sans-serif'}}}};
  }});
  const hlShape=tIdx>=0?[{{
    type:'rect',x0:-0.5,x1:samps.length-0.5,y0:tIdx-0.5,y1:tIdx+0.5,
    line:{{color:'#ff8c00',width:2.5}},fillcolor:'rgba(255,140,0,0.07)'
  }}]:[];
  const titleText=_inTop?`${{circId}} ← highlighted`:`Top DE Heatmap (${{circId}} not in top set)`;
  Plotly.newPlot(el,[{{
    type:'heatmap',z:zMatrix,x:samps,y:yLabels,
    colorscale:[[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']],
    zmid:0,showscale:true,
    colorbar:{{title:'z-score',titlefont:{{size:10}}}},
    hovertemplate:'<b>%{{y}}</b><br>%{{x}}<br>z: %{{z:.2f}}<extra></extra>',
  }}],{{
    height:Math.max(300,allIds.length*20+120),
    margin:{{t:65,b:60,l:180,r:80}},
    title:{{text:titleText,font:{{size:11}}}},
    xaxis:{{tickangle:-35,tickfont:{{size:9}}}},
    yaxis:{{tickfont:{{size:9}},autorange:'reversed'}},
    shapes:[...groupShapes,...hlShape],
    annotations:groupAnno,
    plot_bgcolor:'white',paper_bgcolor:'white',
  }},{{responsive:true,displayModeBar:false}});
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
    if(!chrom||!v||v==='—') return v;
    const pm=String(v).match(/(\\d+)[–-](\\d+)/);
    if(!pm) return v;
    const a=chromStart+parseInt(pm[1])-1;
    const b=chromStart+parseInt(pm[2])-1;
    return chrom+':'+a.toLocaleString()+'-'+b.toLocaleString();
  }};

  const priTitle='Priority score (click to sort):\\n'
    +(tableType==='mirna'
      ?'8mer=+4, 7mer-m8=+3, 7mer-1A=+2, 6mer=+1\\nCLIP exp>0=+3, ENCORI=+2, in_circ=+1'
      :'bindingSites log2×2 (max+3), internal=+2\\nCLIP exp>0=+3, ENCORI=+2, in_circ=+1');
  let html='<table class="itable"><thead><tr>';
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
      if(k==='circ_pos') v=_absPos(v);
      if(k==='_seq_logo') {{
        const rawPos=String(r.circ_pos||'');
        const pm=rawPos.match(/(\\d+)[–-](\\d+)/);
        if(pm&&chrom) {{
          const s0=chromStart+parseInt(pm[1])-1;
          const e0=chromStart+parseInt(pm[2]);
          const len=e0-s0;
          if(len>80)
            v=`<span style="color:#aaa;font-size:10px">${{len}}bp</span>`;
          else
            v=`<span class="_seqCell" data-chrom="${{chrom}}" data-s="${{s0}}" data-e="${{e0}}"
                 style="color:#aaa;font-size:10px;font-family:monospace">⟳</span>`;
        }} else {{ v='—'; }}
      }}
      if(k==='source') v=_srcBadge(String(r[k]||''));
      if(k==='in_circ'){{
        const ic=r[k];
        if(ic===undefined||ic===true||ic==='true')
          v='<span style="color:#2ca02c;font-weight:bold">✓</span>';
        else
          v='<span style="color:#aaa">✗</span>';
      }}
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
}}

function _toggleAll(type, show, container) {{
  const wrap = document.getElementById('cm-circle-wrap');
  if (!wrap) return;
  wrap.querySelectorAll(`._${{type}}badge`).forEach(el => el.style.display = show ? '' : 'none');
  wrap.querySelectorAll(`._${{type}}arc`).forEach(el   => el.style.display = show ? '' : 'none');
  if (container) {{
    container.querySelectorAll(`[id^="_leg_${{type}}_"]`).forEach(chip => {{
      chip.style.opacity        = show ? '1' : '0.3';
      chip.style.textDecoration = show ? '' : 'line-through';
    }});
  }}
}}

// ── Main heatmap dynamic update ───────────────────────────────────────────────
function updateMainHeatmap() {{
  const _hmData=_HEATMAP_DATA_CACHE||FULL_HEATMAP_DATA;
  if(!_hmData||typeof Plotly==='undefined')return;
  const n=Math.max(1,Math.min(50,parseInt(document.getElementById('heatmap-n-input').value)||10));
  document.getElementById('heatmap-n-input').value=n;
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
  const TUMOR_COL='#d62728',NORMAL_COL='#2c6fad';
  const grps=[];let cur=null;
  samps.forEach((s,i)=>{{
    const c=conds[s]||'';
    const col=c==='{tumor_label}'?TUMOR_COL:c==='{normal_label}'?NORMAL_COL:'#888';
    if(!cur||cur.c!==c){{cur={{c,col,s:i,e:i}};grps.push(cur);}}else cur.e=i;
  }});
  const groupShapes=grps.map(g=>{{return{{type:'rect',xref:'x',yref:'paper',
    x0:g.s-0.45,x1:g.e+0.45,y0:1.02,y1:1.055,fillcolor:g.col,line:{{width:0}}}}}});
  const groupAnno=grps.map(g=>{{return{{xref:'x',yref:'paper',x:(g.s+g.e)/2,y:1.0375,
    yanchor:'middle',xanchor:'center',text:g.c,showarrow:false,
    font:{{size:9,color:'white',family:'sans-serif'}}}}}});
  Plotly.react('main-heatmap-plot',[{{
    type:'heatmap',z:zMatrix,x:samps,y:yLabels,
    colorscale:[[0,'#2ca02c'],[0.5,'white'],[1,'#d62728']],
    zmid:0,colorbar:{{title:'z-score<br>(normal-<br>centered)'}},
    hovertemplate:'<b>%{{y}}</b><br>%{{x}}<br>z-score: %{{z:.2f}}<extra></extra>',
  }}],{{
    title:{{text:`Top ${{nUp}} up + ${{nDn}} down DE circRNAs — Heatmap (normal-centered z-score)`,font:{{size:14}}}},
    yaxis:{{tickfont:{{size:8}},autorange:'reversed'}},
    height:Math.max(420,validIds.length*22+120),
    plot_bgcolor:'white',paper_bgcolor:'white',
    margin:{{t:60,l:300}},
    shapes:groupShapes,
    annotations:groupAnno,
  }});
  const titleEl=document.getElementById('heatmap-section-title');
  if(titleEl)titleEl.textContent=`Heatmap (top ${{nUp}} up + ${{nDn}} down DE circRNAs)`;
  const maxUp=(_hmData.up_order||[]).length;
  const maxDn=(_hmData.dn_order||[]).length;
  const statusEl=document.getElementById('heatmap-status');
  if(statusEl)statusEl.textContent=`（已顯示 ${{validIds.length}} 筆；資料庫最多 ${{maxUp}} up + ${{maxDn}} down）`;
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
      xaxis_title:'log₂ Fold Change (Tumor / Normal)',
      yaxis_title:'−log₁₀('+_SIG_LABEL+')',
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
  if(md.heatmap){{
    _HEATMAP_DATA_CACHE=md.heatmap;
    updateMainHeatmap();
  }}
}}

document.addEventListener('keydown',e=>{{if(e.key==='Escape')closeCircModal();}});
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
  <span>circRNA Analysis Report — {project_id}</span>
  <button class="print-btn" onclick="window.print()">🖨 列印 / 存為 PDF</button>
</div>

  <h1>circRNA Analysis Report</h1>
  <p><strong>Project:</strong> {project_id} &nbsp;&nbsp;
     <strong>Method:</strong> <span class="method-tag">{de_method}</span> &nbsp;&nbsp;
     <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}
     {"&nbsp;&nbsp;<a href='multiqc_report.html' target='_blank' style='background:#2c6fad;color:white;padding:3px 10px;border-radius:4px;text-decoration:none;font-size:13px'>&#128202; QC Report (MultiQC)</a>" if multiqc_file and __import__('os').path.exists(multiqc_file) else ""}
  </p>
  {"<p style='font-size:12px;color:#888'>Interaction data pre-fetched for "
   + str(n_ixn) + " circRNAs (" + str(n_ixn_mirna) + " with miRNA data, from CircInteractome). "
   + "Click any circ_position to view exon diagram, miRNA sponge sites, and RBP binding sites.</p>" if n_ixn > 0
   else "<p style='font-size:12px;color:#aaa'>Click any circ_position to view exon diagram "
        "(interaction data not yet fetched — run predict_interactions rule).</p>"}

  <h2>Summary</h2>
  {_msw_html}
  <div>
    <div class="stat-box"><div class="num">{n_sample}</div><div class="lbl">Samples</div></div>
    <div class="stat-box"><div class="num">{n_total}</div><div class="lbl">Total circRNAs</div></div>
    <div class="stat-box"><div class="num" id="stat-n-sig">{n_sig}</div><div class="lbl">Significant ({sig_label}, |log2FC|&gt;{lfc})</div></div>
    <div class="stat-box"><div class="num" id="stat-n-up">{n_up}</div><div class="lbl">Up-regulated</div></div>
    <div class="stat-box"><div class="num" id="stat-n-dn">{n_dn}</div><div class="lbl">Down-regulated</div></div>
  </div>

  {type_html}

  {biomarker_html}

  {isoform_html}

  {"<h2>三方法 DE 結果 Venn Diagram</h2><p style='font-size:13px;color:#555'>比較三種方法（edgeR FSJ offset、DESeq2、limma-voom）在相同閾值下的顯著 DE circRNA 交集。</p>" + venn_html if venn_html else ""}

  <h2>Top Differentially Expressed circRNAs ({sig_label}, |log2FC| &gt; {lfc})</h2>
  <p style="font-size:12px;color:#666">&#128204; Click a <strong>circ_position</strong> to view exon diagram, miRNA sponge sites, and RBP binding sites.</p>
  {_de_split_tables(top_table, tumor_label=tumor_label, normal_label=normal_label, interactions=interactions)}

  <h2>Volcano Plot</h2>
  <p style="font-size:12px;color:#888">&#9711; Heatmap top {heatmap_top_n} up + {heatmap_top_n} down markers: use the toggle button in the chart to show/hide.</p>
  {volcano_html}

  <h2>PCA</h2>
  {pca_html}

  <h2 id="heatmap-section-title">Heatmap (top {heatmap_top_n} up + {heatmap_top_n} down DE circRNAs)</h2>
  <div style="display:flex;align-items:center;gap:10px;margin:8px 0 12px;background:#f4f8ff;padding:10px 16px;border-radius:6px;border:1px solid #d0e4f7;flex-wrap:wrap">
    <span style="font-size:13px">每方向顯示 top</span>
    <input type="number" id="heatmap-n-input" value="{heatmap_top_n}" min="1" max="50"
           style="width:60px;padding:3px 6px;border:1px solid #bbb;border-radius:4px;font-size:13px"
           onkeydown="if(event.key==='Enter')updateMainHeatmap()">
    <span style="font-size:13px">up + down</span>
    <button onclick="updateMainHeatmap()"
            style="background:#2c6fad;color:white;border:none;border-radius:4px;padding:5px 16px;cursor:pointer;font-size:13px">更新 Heatmap</button>
    <span id="heatmap-status" style="font-size:12px;color:#888"></span>
  </div>
  {heatmap_html}

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
    )
