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


def _df_to_html(df: pd.DataFrame, max_rows: int = 50) -> str:
    return (
        df.head(max_rows)
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


def _de_split_tables(sig: pd.DataFrame, tumor_label: str = "tumor",
                     normal_label: str = "normal") -> str:
    """Return two HTML tables: up-regulated and down-regulated in tumor."""
    if sig.empty or "log2FC" not in sig.columns:
        return _df_to_html(sig)
    up   = sig[sig["log2FC"] > 0].copy()
    down = sig[sig["log2FC"] < 0].copy()
    html_parts = []
    if not up.empty:
        html_parts.append(
            f'<h3 style="color:#d62728">&#8593; Up-regulated in {tumor_label} '
            f'(log₂FC &gt; 0) — {len(up)} circRNAs</h3>'
        )
        html_parts.append(_df_to_html(up))
    if not down.empty:
        html_parts.append(
            f'<h3 style="color:#1f77b4">&#8595; Down-regulated in {tumor_label} / '
            f'Up-regulated in {normal_label} (log₂FC &lt; 0) — {len(down)} circRNAs</h3>'
        )
        html_parts.append(_df_to_html(down))
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
</style>
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


def _isoform_section(switching_file: Optional[str]) -> str:
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

    plot_html = _plot_isoform_usage(sig)

    show_cols = [c for c in
                 ["gene_name", "circ_id", "iui_normal", "iui_tumor",
                  "delta_iui", "p_value", "padj_global"]
                 if c in sig.columns]
    table_html = (
        _df_to_html(sig.sort_values("padj_global")[show_cols], max_rows=20)
        if show_cols else ""
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


def _biomarker_section(biomarker_file: Optional[str]) -> str:
    if not biomarker_file or not Path(biomarker_file).exists():
        return ""
    try:
        bm = pd.read_csv(biomarker_file, sep="\t")
    except Exception:
        return ""
    if bm.empty:
        return ""
    show_cols = [c for c in ["rank", "circ_id", "log2FC", "padj", "biomarker_score",
                              "in_circbase", "circbase_id", "circbase_gene", "Type"]
                 if c in bm.columns]
    return f"""
  <h2>Biomarker Candidates (top {min(len(bm), 30)} by composite score)</h2>
  <p style="font-size:13px;color:#555;margin-bottom:8px;">
    Score = mean of: −log₁₀(padj), |log₂FC|, confidence score, circBase known bonus (each normalised 0–1).
  </p>
  {_df_to_html(bm[show_cols], max_rows=30)}
"""


def _plotly_volcano(de: pd.DataFrame, fdr: float, lfc: float, de_method: str,
                    use_pvalue: bool = False) -> str:
    """Interactive Plotly volcano; returns '' when plotly is unavailable or DE is empty."""
    if not _PLOTLY:
        return ""
    p_col = "pvalue" if use_pvalue else "padj"
    if p_col not in de.columns:
        p_col = "padj" if "padj" in de.columns else None
    if p_col is None:
        return ""
    df = de[de[p_col].notna() & de["log2FC"].notna()].copy()
    if df.empty:
        return ""

    df["nlp"] = df[p_col].clip(lower=1e-300).apply(lambda p: -math.log10(p))

    def _sig(row):
        if row[p_col] < fdr and row["log2FC"] > lfc:
            return "Up"
        if row[p_col] < fdr and row["log2FC"] < -lfc:
            return "Down"
        return "NS"

    df["sig"] = df.apply(_sig, axis=1)
    has_type = "Type" in df.columns
    p_label = "p-value (nominal)" if use_pvalue else "adjusted p-value"

    colors = {"Up": "#d62728", "Down": "#1f77b4", "NS": "rgba(150,150,150,0.35)"}
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

    fig = go.Figure(traces)
    fig.add_vline(x=lfc,   line_dash="dot", line_color="#aaa", line_width=1)
    fig.add_vline(x=-lfc,  line_dash="dot", line_color="#aaa", line_width=1)
    fig.add_hline(y=-math.log10(fdr), line_dash="dot", line_color="#aaa", line_width=1)
    fig.update_layout(
        title=dict(text=f"Volcano Plot [{de_method}]", font_size=14),
        xaxis_title="log₂ Fold Change (Tumor / Normal)",
        yaxis_title=f"−log₁₀({p_label})",
        height=500, plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(title="", orientation="h", y=1.02, x=0),
        margin=dict(t=60),
    )
    fig.update_xaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor="#f0f0f0", zeroline=False)
    return fig.to_html(include_plotlyjs="cdn", full_html=False)


def _plotly_heatmap(de: pd.DataFrame, matrix: pd.DataFrame, top_n: int = 50,
                    use_pvalue: bool = False) -> str:
    """Interactive Plotly heatmap of top DE circRNAs (log2 + z-score)."""
    if not _PLOTLY:
        return ""
    p_col = "pvalue" if (use_pvalue and "pvalue" in de.columns) else "padj"
    if "circ_id" not in de.columns or p_col not in de.columns:
        return ""

    top_ids = (
        de.dropna(subset=[p_col]).sort_values(p_col).head(top_n)["circ_id"].tolist()
    )
    avail = [i for i in top_ids if i in matrix.index]
    if len(avail) < 2:
        return ""

    sub = matrix.loc[avail].astype(float)
    log_sub = np.log2(sub + 1)
    row_mean = log_sub.mean(axis=1)
    row_std  = log_sub.std(axis=1).clip(lower=0.01)
    z = log_sub.sub(row_mean, axis=0).div(row_std, axis=0)

    fig = go.Figure(go.Heatmap(
        z=z.values.tolist(), x=z.columns.tolist(), y=z.index.tolist(),
        colorscale="RdBu_r", zmid=0,
        hovertemplate="<b>%{y}</b><br>%{x}<br>Z-score: %{z:.2f}<extra></extra>",
    ))
    fig.update_layout(
        title=dict(text=f"Top {len(avail)} DE circRNAs — Heatmap (log₂, z-scored)", font_size=14),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
        height=max(420, len(avail) * 14 + 120),
        plot_bgcolor="white", paper_bgcolor="white",
        margin=dict(t=60, l=180),
    )
    return fig.to_html(include_plotlyjs="cdn", full_html=False)


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
    color_map  = {"tumor": "#d62728", "normal": "#1f77b4", "unknown": "#888"}
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
    use_pvalue:     bool  = False,
    tumor_label:    str   = "tumor",
    normal_label:   str   = "normal",
    isoform_file:   Optional[str] = None,
    circbase_file:  Optional[str] = None,
) -> None:
    de     = pd.read_csv(de_file, sep="\t")
    matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)

    # Normalise column names — analysis.R always outputs log2FC
    if "log2FC" not in de.columns and "log2FoldChange" in de.columns:
        de = de.rename(columns={"log2FoldChange": "log2FC"})

    # Enrich with host gene, circBase ID, exon span
    de = _enrich_de(de, isoform_file, circbase_file)

    p_col = "pvalue" if (use_pvalue and "pvalue" in de.columns) else "padj"
    sig_mask = (de[p_col] < fdr) & (de["log2FC"].abs() > lfc) if p_col in de.columns else pd.Series(False, index=de.index)
    sig: pd.DataFrame = de.loc[sig_mask]

    n_total  = len(matrix)
    n_sig    = len(sig)
    n_up     = int((sig["log2FC"] > 0).sum()) if len(sig) else 0
    n_dn     = int((sig["log2FC"] < 0).sum()) if len(sig) else 0
    n_sample = matrix.shape[1]
    sig_label = f"p&lt;{fdr}" if use_pvalue else f"FDR&lt;{fdr}"

    # Rename circ_id to circ_position for display
    sig = sig.copy()
    if "circ_id" in sig.columns:
        sig = sig.rename(columns={"circ_id": "circ_position"})

    top_cols = [c for c in [
        "circ_position", "gene_name", "strand", "region", "exon_span", "circbase_id",
        "log2FC", "pvalue", "padj", "Type",
    ] if c in sig.columns]
    top_table = sig.sort_values(p_col)[top_cols] if top_cols else sig.head(20)

    type_html      = _type_section(sig)
    biomarker_html = _biomarker_section(biomarker_file)
    isoform_html   = _isoform_section(switching_file)

    # Interactive Plotly charts; fall back to static PDF embeds when unavailable
    p_volcano = _plotly_volcano(de, fdr, lfc, de_method, use_pvalue=use_pvalue)
    p_heatmap = _plotly_heatmap(de, matrix, use_pvalue=use_pvalue)
    p_pca     = _plotly_pca(matrix, groups_file)
    volcano_html = p_volcano if p_volcano else _embed_pdf(volcano_pdf)
    heatmap_html = p_heatmap if p_heatmap else _embed_pdf(heatmap_pdf)
    pca_html     = p_pca     if p_pca     else _embed_pdf(pca_pdf)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>circRNA Analysis Report – {project_id}</title>
  {_STYLE}
</head>
<body>
  <h1>circRNA Analysis Report</h1>
  <p><strong>Project:</strong> {project_id} &nbsp;&nbsp;
     <strong>Method:</strong> <span class="method-tag">{de_method}</span> &nbsp;&nbsp;
     <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

  <h2>Summary</h2>
  <div>
    <div class="stat-box"><div class="num">{n_sample}</div><div class="lbl">Samples</div></div>
    <div class="stat-box"><div class="num">{n_total}</div><div class="lbl">Total circRNAs</div></div>
    <div class="stat-box"><div class="num">{n_sig}</div><div class="lbl">Significant ({sig_label}, |log2FC|&gt;{lfc})</div></div>
    <div class="stat-box"><div class="num">{n_up}</div><div class="lbl">Up-regulated</div></div>
    <div class="stat-box"><div class="num">{n_dn}</div><div class="lbl">Down-regulated</div></div>
  </div>

  {type_html}

  {biomarker_html}

  {isoform_html}

  <h2>Top Differentially Expressed circRNAs ({sig_label}, |log2FC| &gt; {lfc})</h2>
  {_de_split_tables(top_table, tumor_label=tumor_label, normal_label=normal_label)}

  <h2>Volcano Plot</h2>
  {volcano_html}

  <h2>PCA</h2>
  {pca_html}

  <h2>Heatmap (top 50 DE circRNAs)</h2>
  {heatmap_html}

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
        use_pvalue     = bool(getattr(snakemake.params, "use_pvalue", False)),  # type: ignore[name-defined]
        tumor_label    = str(snakemake.params.tumor_label),     # type: ignore[name-defined]
        normal_label   = str(snakemake.params.normal_label),    # type: ignore[name-defined]
        isoform_file   = getattr(snakemake.input, "isoform_groups", None),  # type: ignore[name-defined]
        circbase_file  = getattr(snakemake.input, "circbase_annot", None),  # type: ignore[name-defined]
    )
