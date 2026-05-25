"""
generate_report.py – Build a self-contained HTML summary report.

Called as a Snakemake script.  Embeds PDF plots as base64-encoded images
and includes top DE results as a table.
"""

from __future__ import annotations

import base64
from datetime import datetime
from pathlib import Path

import pandas as pd


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
    n_type1 = int((sig["Type"] == "I").sum())
    n_type2 = int((sig["Type"] == "II").sum())
    total   = n_type1 + n_type2
    if total == 0:
        return ""
    pct1 = round(100 * n_type1 / total)
    pct2 = 100 - pct1
    return f"""
  <h2>Type I / II Classification <small style="font-size:0.7em;color:#888">(edgeR + FSJ offset)</small></h2>
  <p>
    <span class="badge badge-type1">Type I</span>&nbsp;circRNA-specific change (FSJ stable) &nbsp;|&nbsp;
    <span class="badge badge-type2">Type II</span>&nbsp;Gene-level change drives BSJ (FSJ also DE)
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


def _biomarker_section(biomarker_file: str | None) -> str:
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
    biomarker_file: str | None = None,
) -> None:
    de     = pd.read_csv(de_file, sep="\t")
    matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)

    # Normalise column names — analysis.R always outputs log2FC
    if "log2FC" not in de.columns and "log2FoldChange" in de.columns:
        de = de.rename(columns={"log2FoldChange": "log2FC"})

    sig_mask = (de["padj"] < fdr) & (de["log2FC"].abs() > lfc) if "padj" in de.columns else pd.Series(False, index=de.index)
    sig: pd.DataFrame = de.loc[sig_mask]

    n_total  = len(matrix)
    n_sig    = len(sig)
    n_up     = int((sig["log2FC"] > 0).sum()) if len(sig) else 0
    n_dn     = int((sig["log2FC"] < 0).sum()) if len(sig) else 0
    n_sample = matrix.shape[1]

    # Top table — include Type column when present
    top_cols = [c for c in ["circ_id", "log2FC", "pvalue", "padj", "Type"] if c in sig.columns]
    top_table = sig.sort_values("padj")[top_cols] if top_cols else sig.head(20)

    type_html = _type_section(sig)
    biomarker_html = _biomarker_section(biomarker_file)

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
    <div class="stat-box"><div class="num">{n_sig}</div><div class="lbl">Significant (FDR&lt;{fdr}, |log2FC|&gt;{lfc})</div></div>
    <div class="stat-box"><div class="num">{n_up}</div><div class="lbl">Up-regulated</div></div>
    <div class="stat-box"><div class="num">{n_dn}</div><div class="lbl">Down-regulated</div></div>
  </div>

  {type_html}

  {biomarker_html}

  <h2>Top Differentially Expressed circRNAs (FDR &lt; {fdr}, |log2FC| &gt; {lfc})</h2>
  {_df_to_html(top_table)}

  <h2>Volcano Plot</h2>
  {_embed_pdf(volcano_pdf)}

  <h2>PCA</h2>
  {_embed_pdf(pca_pdf)}

  <h2>Heatmap (top 50 DE circRNAs)</h2>
  {_embed_pdf(heatmap_pdf)}

</body>
</html>
"""
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    Path(output_file).write_text(html, encoding="utf-8")
    print(f"[OK] Report written → {output_file}")


# ── Snakemake entry point ────────────────────────────────────────────────────

if "snakemake" in dir():
    build_report(
        de_file        = snakemake.input.de,              # type: ignore[name-defined]
        matrix_file    = snakemake.input.matrix,          # type: ignore[name-defined]
        volcano_pdf    = snakemake.input.volcano,         # type: ignore[name-defined]
        heatmap_pdf    = snakemake.input.heatmap,         # type: ignore[name-defined]
        pca_pdf        = snakemake.input.pca,             # type: ignore[name-defined]
        output_file    = snakemake.output[0],             # type: ignore[name-defined]
        project_id     = snakemake.params.project_id,    # type: ignore[name-defined]
        fdr            = float(snakemake.params.fdr),    # type: ignore[name-defined]
        lfc            = float(snakemake.params.lfc),    # type: ignore[name-defined]
        de_method      = str(snakemake.params.de_method),  # type: ignore[name-defined]
        biomarker_file = snakemake.input.biomarkers,     # type: ignore[name-defined]
    )
