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
    """Encode a PDF as base64 for embedding in HTML via <embed>."""
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
</style>
"""


def build_report(
    de_file:      str,
    matrix_file:  str,
    volcano_pdf:  str,
    heatmap_pdf:  str,
    pca_pdf:      str,
    output_file:  str,
    project_id:   str,
    fdr:          float = 0.05,
    lfc:          float = 1.0,
) -> None:
    de     = pd.read_csv(de_file, sep="\t")
    matrix = pd.read_csv(matrix_file, sep="\t", index_col=0)

    bool_mask = (de["padj"] < fdr) & (de["log2FoldChange"].abs() > lfc)
    sig: pd.DataFrame = de.loc[bool_mask] if "padj" in de.columns else de.head(0)
    n_total  = len(matrix)
    n_sig    = len(sig)
    n_up     = int((sig["log2FoldChange"] > 0).sum()) if len(sig) else 0
    n_dn     = int((sig["log2FoldChange"] < 0).sum()) if len(sig) else 0
    n_sample = matrix.shape[1]

    top_table: pd.DataFrame = (
        sig.sort_values("padj")
           [["circ_id", "log2FoldChange", "pvalue", "padj"]]
           .rename(columns={"log2FoldChange": "log2FC"})
        if all(c in sig.columns for c in ["circ_id", "log2FoldChange", "pvalue", "padj"])
        else sig.head(20)
    )

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
     <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

  <h2>Summary</h2>
  <div>
    <div class="stat-box"><div class="num">{n_sample}</div><div class="lbl">Samples</div></div>
    <div class="stat-box"><div class="num">{n_total}</div><div class="lbl">Total circRNAs</div></div>
    <div class="stat-box"><div class="num">{n_sig}</div><div class="lbl">Significant (FDR&lt;{fdr}, |log2FC|&gt;{lfc})</div></div>
    <div class="stat-box"><div class="num">{n_up}</div><div class="lbl">Up-regulated</div></div>
    <div class="stat-box"><div class="num">{n_dn}</div><div class="lbl">Down-regulated</div></div>
  </div>

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
        de_file     = snakemake.input.de,           # type: ignore[name-defined]
        matrix_file = snakemake.input.matrix,       # type: ignore[name-defined]
        volcano_pdf = snakemake.input.volcano,      # type: ignore[name-defined]
        heatmap_pdf = snakemake.input.heatmap,      # type: ignore[name-defined]
        pca_pdf     = snakemake.input.pca,          # type: ignore[name-defined]
        output_file = snakemake.output[0],          # type: ignore[name-defined]
        project_id  = snakemake.params.project_id,  # type: ignore[name-defined]
        fdr         = float(snakemake.params.fdr),  # type: ignore[name-defined]
        lfc         = float(snakemake.params.lfc),  # type: ignore[name-defined]
    )
