"""
app.py – Streamlit web UI for the circRNA Agent.

Run: streamlit run app.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent / "scripts"))
from utils import load_config, save_config, generate_ciriquant_config, ensure_dirs

# ─────────────────────────────────────────────────────────────────────────────
st.set_page_config(page_title="circRNA Agent", page_icon="🧬", layout="wide")
st.title("🧬 circRNA Agent")
st.caption("End-to-end circRNA analysis from GEO/SRA to HTML report")

cfg = load_config("config.yaml")

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Settings")
    cfg["threads"] = st.slider("CPU cores", 1, 64, int(cfg.get("threads", 8)))
    st.markdown("---")
    st.markdown("**Significance thresholds**")
    cfg["de"]["fdr_cutoff"]    = st.number_input("FDR cut-off",    value=float(cfg["de"]["fdr_cutoff"]),    step=0.01)
    cfg["de"]["log2fc_cutoff"] = st.number_input("|log₂FC| cut-off", value=float(cfg["de"]["log2fc_cutoff"]), step=0.1)
    if st.button("💾 Save settings"):
        save_config(cfg)
        st.success("config.yaml updated")

# ── Tabs ──────────────────────────────────────────────────────────────────────
tab_input, tab_groups, tab_run, tab_results = st.tabs(
    ["📥 Input", "🏷️ Groups", "🚀 Run", "📊 Results"]
)

# ── Tab 1: Input ──────────────────────────────────────────────────────────────
with tab_input:
    st.subheader("Data source")
    source = st.radio("", ["GEO accession", "Upload RunInfo.csv"], horizontal=True)

    if source == "GEO accession":
        gse = st.text_input("GEO accession", placeholder="e.g. GSE113230")
        if st.button("📥 Fetch metadata"):
            if not gse.strip():
                st.error("Please enter a GSE accession.")
            else:
                with st.spinner(f"Querying SRA for {gse} …"):
                    try:
                        from download_geo import prepare_metadata
                        ensure_dirs()
                        df = prepare_metadata(gse_id=gse.strip(), output_file=cfg["metadata"])
                        cfg["project_id"] = gse.strip()
                        save_config(cfg)
                        st.success(f"✅ {len(df)} samples fetched")
                        st.dataframe(df)
                    except Exception as exc:
                        st.error(str(exc))
    else:
        uploaded = st.file_uploader("Upload RunInfo.csv", type="csv")
        if uploaded and st.button("📥 Load CSV"):
            tmp = Path("metadata/_uploaded_runinfo.csv")
            tmp.parent.mkdir(parents=True, exist_ok=True)
            tmp.write_bytes(uploaded.read())
            from download_geo import prepare_metadata
            ensure_dirs()
            df = prepare_metadata(runinfo_csv=str(tmp), output_file=cfg["metadata"])
            st.success(f"✅ {len(df)} samples loaded")
            st.dataframe(df)

    # Show current metadata
    if Path(cfg["metadata"]).exists():
        st.markdown("---")
        st.markdown("**Current metadata**")
        st.dataframe(pd.read_csv(cfg["metadata"]))

# ── Tab 2: Groups ─────────────────────────────────────────────────────────────
with tab_groups:
    st.subheader("Tumor / Normal assignment")

    meta_path   = Path(cfg["metadata"])
    groups_path = Path(cfg["groups"])

    if not meta_path.exists():
        st.info("Fetch metadata in the Input tab first.")
    else:
        meta = pd.read_csv(meta_path)

        if groups_path.exists():
            groups_df = pd.read_csv(groups_path)
        else:
            from prepare_metadata import auto_assign
            groups_df = auto_assign(meta)

        st.markdown("Edit the **condition** column (values: `tumor` or `normal`):")
        edited = st.data_editor(
            groups_df,
            column_config={"condition": st.column_config.SelectboxColumn(
                "Condition", options=["tumor", "normal", "unknown"])},
            use_container_width=True,
            num_rows="fixed",
        )
        if st.button("💾 Save groups"):
            groups_path.parent.mkdir(parents=True, exist_ok=True)
            edited.to_csv(groups_path, index=False)
            st.success(f"Saved → {groups_path}")

# ── Tab 3: Run ────────────────────────────────────────────────────────────────
with tab_run:
    st.subheader("Launch pipeline")

    col1, col2 = st.columns(2)
    dry_run = col1.checkbox("Dry run (show plan only)")
    col2.markdown(f"**Cores:** {cfg['threads']}")

    if st.button("🚀 Run Snakemake pipeline", type="primary"):
        # Regenerate CIRIquant config in case genome paths changed
        generate_ciriquant_config(cfg, output=cfg["ciriquant_config"])
        save_config(cfg)

        cmd = [
            sys.executable, "-m", "snakemake",
            "--snakefile", "workflow/Snakefile",
            "--cores", str(cfg["threads"]),
            "--rerun-incomplete",
            "--printshellcmds",
        ]
        if dry_run:
            cmd.append("--dry-run")

        st.code(" ".join(cmd), language="bash")
        with st.spinner("Running … (check your terminal for live logs)"):
            result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            st.success("✅ Pipeline finished successfully!")
        else:
            st.error("❌ Pipeline failed – see error below")

        if result.stdout:
            with st.expander("stdout"):
                st.code(result.stdout[-4000:])
        if result.stderr:
            with st.expander("stderr"):
                st.code(result.stderr[-4000:])

# ── Tab 4: Results ────────────────────────────────────────────────────────────
with tab_results:
    st.subheader("Analysis results")

    report = Path("results/report.html")
    de_tsv  = Path("results/de/de_results.tsv")
    matrix  = Path("results/circRNA/count_matrix.tsv")

    cols = st.columns(3)
    cols[0].metric("Report", "✅ Ready" if report.exists() else "⏳ Not yet")
    cols[1].metric("DE results", "✅ Ready" if de_tsv.exists()  else "⏳ Not yet")
    cols[2].metric("Count matrix", "✅ Ready" if matrix.exists() else "⏳ Not yet")

    if de_tsv.exists():
        de_df  = pd.read_csv(de_tsv, sep="\t")
        fdr    = cfg["de"]["fdr_cutoff"]
        lfc    = cfg["de"]["log2fc_cutoff"]
        sig_df = de_df[(de_df["padj"] < fdr) & (de_df["log2FoldChange"].abs() > lfc)] \
                 if all(c in de_df.columns for c in ["padj", "log2FoldChange"]) else de_df

        st.metric("Significant circRNAs", len(sig_df))
        st.markdown("**Top 20 DE circRNAs**")
        st.dataframe(sig_df.head(20))

        csv_bytes = sig_df.to_csv(index=False).encode()
        st.download_button("⬇️ Download DE results (CSV)", csv_bytes,
                           file_name="circRNA_DE.csv", mime="text/csv")

    for plot in ["results/plots/volcano.pdf", "results/plots/heatmap.pdf", "results/plots/pca.pdf"]:
        p = Path(plot)
        if p.exists():
            with open(p, "rb") as fh:
                st.download_button(
                    label=f"⬇️ {p.name}",
                    data=fh.read(),
                    file_name=p.name,
                    mime="application/pdf",
                )
