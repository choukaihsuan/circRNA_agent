# thesis cross-dataset numbers (19 vs 15 datasets), 2026-10-03
Inputs: `data/` (read-only copy from HPC server, `de_results.tsv` columns circ_id/log2FC/pvalue/Type, isoform_groups gene map, audit/per_dataset_summary.tsv).
- `python3 make_thesis_numbers.py` -> `../thesis_numbers_15ds.tsv`, `../thesis_dataset_table_19ds.tsv` (stdlib only)
- `python3 make_figures.py` (matplotlib, numpy) -> `figs/main_15ds/` (main) and `figs/appendix_19ds/`
Significance: pvalue<0.05 & |log2FC|>1 (edgeR_ciriquant main DE). 15-set = 19 minus GSE58135, GSE229705, GSE192410, GSE192849 (effective single-tool >=97.8% of consensus rows).
Not modified: the paper, web_ui cross-dataset page (still shows all datasets).
