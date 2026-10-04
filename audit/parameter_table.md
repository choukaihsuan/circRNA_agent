| 參數名 | 程式碼預設值 | 已發表執行使用值 | 狀態 | 作用 | 定義位置 |
|---|---|---|---|---|---|
| `de_method` | edgeR_ciriquant | varies: edgeR_ciriquant (19); limma (1) [config.de.method] | DIFF | analysis.R parameter; passed from workflow/Snakefile (DE_METHOD); analysis.R fallback deseq2 | scripts/analysis.R:71 |
| `de_sig_by` | pvalue | pvalue [config.de.de_sig_by] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:35 |
| `deseq2_size_factor_type` | poscounts | poscounts [config.de.deseq2_size_factor_type] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:58 |
| `fdr` | 0.05 | 0.05 [config.de.fdr_cutoff] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:31 |
| `filter_expr_min_count` | 5 | 5 [config.de.filter_expr_min_count] | same | analysis.R parameter; passed from workflow/rules/de.smk; analysis.R fallback 5L | scripts/analysis.R:50 |
| `filter_expr_min_total_count` | 15 | 15 [config.de.filter_expr_min_total_count] | same | analysis.R parameter; passed from workflow/rules/de.smk; analysis.R fallback 15L | scripts/analysis.R:53 |
| `fsj_concordance_lfc` | 0.0 | 0.0 [config.de.fsj_concordance_lfc] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:42 |
| `heatmap_top_n` | 10 | 10 [config.de.heatmap_top_n] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:44 |
| `lfc` | 1.0 | 1.0 [config.de.log2fc_cutoff] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:32 |
| `normal_label` | normal | normal [config.de.normal_label] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:34 |
| `paired_design_override` | auto | auto [config.de.paired_design_override] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:65 |
| `tumor_label` | tumor | tumor [config.de.tumor_label] | same | analysis.R parameter; passed from workflow/rules/de.smk | scripts/analysis.R:33 |
| `--adaptive` | False | not passed by pipeline (code default) | n/a | Enable adaptive fallback: if one tool detects far fewer circRNAs than another (< --adaptive-ratio), automatically lower min_tools to 1 to avoid near-zero recall | scripts/consensus_filter.py:346 |
| `--adaptive-ratio` | 0.1 | 0.1 [config.consensus.adaptive_ratio] | same | Threshold for adaptive fallback: if min(tool_counts) / max(tool_counts) < ratio, trigger single-tool mode (default: 0.1) | scripts/consensus_filter.py:350 |
| `--annot` | (required) | not passed by pipeline (code default) | n/a | circBase-annotated summary TSV | scripts/rank_biomarkers.py:523 |
| `--case-label` | tumor | not passed by pipeline (code default) | n/a | Label for the case/treatment group (default: tumor) | scripts/prepare_metadata.py:272 |
| `--circbase` | (not given) | not passed by pipeline (code default) | n/a | circbase_annotated.tsv for circBase ID lookup | scripts/predict_interactions.py:800 |
| `--circbase-file` | (not given) | not passed by pipeline (code default) | n/a | Local circBase hg19 file (auto-downloaded if absent) | scripts/annotate_circbase.py:208 |
| `--circexplorer2` | (not given) | not passed by pipeline (code default) | n/a | CIRCexplorer2 back_spliced_junction.bed or known_circ.txt (optional) | scripts/consensus_filter.py:327 |
| `--ciri2` | (not given) | not passed by pipeline (code default) | n/a | CIRI2 output text file (optional) | scripts/consensus_filter.py:325 |
| `--cirique` | (not given) | not passed by pipeline (code default) | n/a | CIRIquant GTF output (optional) | scripts/consensus_filter.py:324 |
| `--clip-exp-num` | 1 | 1 [config.de.encori_clip_exp_num] | same | ENCORI: min CLIP experiments (1=all, 2=moderate, 3=high) | scripts/predict_interactions.py:808 |
| `--cluster` | (not given) | not passed by pipeline (code default) | n/a | Snakemake --cluster submit command | scripts/agent.py:132 |
| `--components-json` | (not given) | not passed by pipeline (code default) | n/a | JSON dict: {component: {enabled, weight}} for sig/fc/conf/known/mirna/rbp. Default: all enabled, weight 1.0 (equal, this paper's setting). | scripts/rank_biomarkers.py:549 |
| `--control-label` | normal | not passed by pipeline (code default) | n/a | Label for the control group (default: normal) | scripts/prepare_metadata.py:274 |
| `--cores` | 8 | not passed by pipeline (code default) | n/a | CPU cores (default: 8) | scripts/agent.py:129 |
| `--count-matrix` | (required) | not passed by pipeline (code default) | n/a | BSJ count matrix TSV (row names: chrom:start|end) | scripts/assign_isoforms.py:283 |
| `--cross-method-min` | 3 | not passed by pipeline (code default) | n/a | Minimum number of DE methods (of those supplied) in which a candidate must be significant to count as cross-method (default 3 = all supplied). | scripts/rank_biomarkers.py:560 |
| `--datasets` | [] | not passed by pipeline (code default) | n/a | Dataset names to combine with --results-root | scripts/cross_dataset_biomarkers.py:199 |
| `--dcc` | (not given) | not passed by pipeline (code default) | n/a | DCC CircCoordinates file (optional) | scripts/consensus_filter.py:326 |
| `--de` | (required) | not passed by pipeline (code default) | n/a |  | scripts/predict_interactions.py:797 |
| `--de` | (required) | not passed by pipeline (code default) | n/a | DE results TSV (from analysis.R) | scripts/rank_biomarkers.py:522 |
| `--de-deseq2` | (not given) | not passed by pipeline (code default) | n/a | DESeq2 DE results TSV | scripts/predict_interactions.py:805 |
| `--de-deseq2` |  | not passed by pipeline (code default) | n/a | DESeq2 DE results TSV (optional; for n_sig_methods) | scripts/rank_biomarkers.py:527 |
| `--de-edger` | (not given) | not passed by pipeline (code default) | n/a | edgeR DE results TSV | scripts/predict_interactions.py:804 |
| `--de-edger` |  | not passed by pipeline (code default) | n/a | edgeR DE table (required unless edgeR is the primary, in which case --de is used) | scripts/rank_biomarkers.py:566 |
| `--de-limma` | (not given) | not passed by pipeline (code default) | n/a | limma DE results TSV | scripts/predict_interactions.py:806 |
| `--de-limma` |  | not passed by pipeline (code default) | n/a | limma DE results TSV (optional; for n_sig_methods) | scripts/rank_biomarkers.py:528 |
| `--de-sig-by` | auto | pvalue [config.de.de_sig_by] | same | Significance column: auto=cascade (q<0.2→p<0.05), pvalue/padj/qvalue=fixed | scripts/rank_biomarkers.py:531 |
| `--dry-run` | (not given) | not passed by pipeline (code default) | n/a | Show Snakemake plan without executing | scripts/agent.py:131 |
| `--event` | (required) | not passed by pipeline (code default) | n/a | 通知類型 | scripts/notify.py:219 |
| `--exclude` | GSE192410 | not passed by pipeline (code default) | n/a | Dataset to drop for the sensitivity analysis | scripts/run_permutation.py:253 |
| `--fc-cap` | 5.0 | not passed by pipeline (code default) | n/a | Cap for |log2FC| in fc_norm component (default: 5) | scripts/rank_biomarkers.py:539 |
| `--fdr` | 0.05 | 0.05 [config.de.fdr_cutoff] | same |  | scripts/rank_biomarkers.py:529 |
| `--filter-bed` | (not given) | not passed by pipeline (code default) | n/a | High-confidence BED files for filtering (optional) | scripts/merge_counts.py:104 |
| `--find-circ` | (not given) | not passed by pipeline (code default) | n/a | find_circ splice_sites.bed output (optional) | scripts/consensus_filter.py:329 |
| `--fully-evaluated-only` | (not given) | not passed by pipeline (code default) | n/a | Also write <output>_fully_evaluated.tsv (or --fully-evaluated-output path) containing only candidates with interaction_queried=True, re-ranked within that evidence-consistent subset. Does not replace the main --output file. | scripts/rank_biomarkers.py:541 |
| `--fully-evaluated-output` | (not given) | not passed by pipeline (code default) | n/a | Output path for the fully-evaluated-only TSV (default: <output> with _fully_evaluated inserted before the extension) | scripts/rank_biomarkers.py:546 |
| `--groups` | metadata/sample_groups.csv | not passed by pipeline (code default) | n/a |  | scripts/prepare_metadata.py:269 |
| `--gse` | (not given) | not passed by pipeline (code default) | n/a | GEO or SRA accession: GSE113230 / PRJNA808398 / SRP156355 | scripts/agent.py:126 |
| `--gse` | (not given) | not passed by pipeline (code default) | n/a | GEO accession (e.g. GSE113230) | scripts/download_geo.py:238 |
| `--gtf` | (required) | /home3/choukaihsuan/reference/hg19/genes.gtf [config.genome.gtf] | same | Genome annotation GTF | scripts/assign_isoforms.py:287 |
| `--gtf` | (not given) | /home3/choukaihsuan/reference/hg19/genes.gtf [config.genome.gtf] | same | Reference GTF for exon boundary lengths (optional) | scripts/predict_interactions.py:802 |
| `--gtfs` | (required) | not passed by pipeline (code default) | n/a | CIRIquant GTF files | scripts/merge_counts.py:100 |
| `--host` | 0.0.0.0 | not passed by pipeline (code default) | n/a |  | scripts/web_ui.py:2050 |
| `--inputs` | [] | not passed by pipeline (code default) | n/a | biomarker_candidates.tsv paths (name inferred from *_results/ path); or NAME=path pairs to set the name explicitly | scripts/cross_dataset_biomarkers.py:194 |
| `--interactions` |  | not passed by pipeline (code default) | n/a | interactions.json from predict_interactions (optional; adds miRNA/RBP dimensions) | scripts/rank_biomarkers.py:525 |
| `--interactive` | (not given) | not passed by pipeline (code default) | n/a | Manually assign tumor/normal labels | scripts/agent.py:130 |
| `--interactive` | (not given) | not passed by pipeline (code default) | n/a | Prompt user for each sample's condition | scripts/prepare_metadata.py:270 |
| `--iso` | (required) | not passed by pipeline (code default) | n/a |  | scripts/predict_interactions.py:798 |
| `--isoforms` |  | not passed by pipeline (code default) | n/a | isoform_groups.tsv (assign_isoforms.py): adds gene_name and region to the candidate file | scripts/rank_biomarkers.py:568 |
| `--jobs` | 50 | not passed by pipeline (code default) | n/a | Max concurrent cluster jobs | scripts/agent.py:133 |
| `--lfc` | 1.0 | 1.0 [config.de.log2fc_cutoff] | same |  | scripts/rank_biomarkers.py:530 |
| `--log` |  | not passed by pipeline (code default) | n/a | Pipeline log 檔路徑（failure 用） | scripts/notify.py:228 |
| `--loo-dir` |  | not passed by pipeline (code default) | n/a | Directory from scripts/loo_stability.R (NONE/ + one dir per dropped sample). If given, effect-size stability is tested for three-method candidates; otherwise effect_stable is NA. | scripts/rank_biomarkers.py:553 |
| `--max-junction-ratio` | 1.0 | 1.0 [config.consensus.max_junction_ratio] | same | Max BSJ/FSJ ratio for CIRIquant pseudo-circ QC (default: 1.0) | scripts/consensus_filter.py:339 |
| `--metadata` | metadata/library_info.csv | not passed by pipeline (code default) | n/a |  | scripts/prepare_metadata.py:268 |
| `--metadata-root` | metadata | not passed by pipeline (code default) | n/a | Root dir holding per-dataset subdirs (default: metadata) | scripts/validate_condition_labels.py:144 |
| `--min-bsj` | 2 | 2 [config.consensus.min_bsj_reads] | same | Minimum BSJ reads per tool before voting (default: 2) | scripts/consensus_filter.py:337 |
| `--min-confidence` | 0.0 | not passed by pipeline (code default) | n/a | Minimum confidence score threshold after voting (default: 0.0). Removes low-evidence circRNAs to improve precision. | scripts/consensus_filter.py:354 |
| `--min-datasets` | 2 | not passed by pipeline (code default) | n/a | Recurrence threshold to flag/report (default 2) | scripts/cross_dataset_biomarkers.py:205 |
| `--min-tools` | 2 | varies: 2 (17); 1 (3) [config.consensus.min_tools] | DIFF | Minimum tools that must agree (default: 2) | scripts/consensus_filter.py:333 |
| `--n-perm` | 10000 | not passed by pipeline (code default) | n/a |  | scripts/run_permutation.py:255 |
| `--no-run-params` | (not given) | not passed by pipeline (code default) | n/a | Do not write interaction coverage into <results>/run_params.yaml. Required for v2 runs so published results directories are not modified. | scripts/rank_biomarkers.py:570 |
| `--out` | (required) | not passed by pipeline (code default) | n/a | Output isoform groups TSV | scripts/assign_isoforms.py:288 |
| `--out` | (not given) | not passed by pipeline (code default) | n/a | Output recurrence TSV | scripts/cross_dataset_biomarkers.py:207 |
| `--out` | metadata/library_info.csv | not passed by pipeline (code default) | n/a | Output CSV path | scripts/download_geo.py:240 |
| `--out` | (required) | not passed by pipeline (code default) | n/a |  | scripts/predict_interactions.py:799 |
| `--out` | (not given) | not passed by pipeline (code default) | n/a | Optional TSV of every sample's verdict | scripts/validate_condition_labels.py:146 |
| `--out-dir` | (not static) | not passed by pipeline (code default) | n/a |  | scripts/run_permutation.py:252 |
| `--output` | (required) | not passed by pipeline (code default) | n/a | Output annotated TSV | scripts/annotate_circbase.py:206 |
| `--output` | (required) | not passed by pipeline (code default) | n/a | Output BED file | scripts/consensus_filter.py:331 |
| `--output` | (required) | not passed by pipeline (code default) | n/a | Output BSJ count matrix TSV | scripts/merge_counts.py:101 |
| `--output` | (required) | not passed by pipeline (code default) | n/a | Output ranked TSV | scripts/rank_biomarkers.py:524 |
| `--output-fsj` | (not given) | not passed by pipeline (code default) | n/a | Output FSJ count matrix TSV (optional) | scripts/merge_counts.py:102 |
| `--port` | 5000 | not passed by pipeline (code default) | n/a |  | scripts/web_ui.py:2051 |
| `--primary-method` | edgeR | not passed by pipeline (code default) | n/a | Which method's table is --de (the scored table). The cross-method flag and LOO gate are the same whichever method is primary. | scripts/rank_biomarkers.py:563 |
| `--profile` | (not given) | not passed by pipeline (code default) | n/a | Snakemake profile directory | scripts/agent.py:134 |
| `--program-num` | 2 | 2 [config.de.encori_program_num] | same | ENCORI: min prediction programs agreeing | scripts/predict_interactions.py:810 |
| `--project` | circRNA | not passed by pipeline (code default) | n/a | 專案 ID（如 GSE113230） | scripts/notify.py:222 |
| `--qc-bsj-threshold` | 5 | 5 [config.consensus.qc_bsj_threshold] | same | BSJ count threshold: only apply ratio QC when BSJ < this value (default: 5). High-BSJ circRNAs are exempt from ratio filtering. | scripts/consensus_filter.py:342 |
| `--report` |  | not passed by pipeline (code default) | n/a | report.html 路徑（success 用） | scripts/notify.py:224 |
| `--reset` | (not given) | not passed by pipeline (code default) | n/a | Remove existing groups file first | scripts/agent.py:143 |
| `--results-root` | (not given) | not passed by pipeline (code default) | n/a | Root holding {DATASET}_results/de/biomarker_candidates.tsv | scripts/cross_dataset_biomarkers.py:197 |
| `--rule` | unknown | not passed by pipeline (code default) | n/a | 失敗的 rule 名稱（failure 用） | scripts/notify.py:226 |
| `--runinfo` | (not given) | not passed by pipeline (code default) | n/a | Path to NCBI SRA RunInfo.csv | scripts/agent.py:128 |
| `--runinfo` | (not given) | not passed by pipeline (code default) | n/a | Path to NCBI SRA RunInfo.csv | scripts/download_geo.py:239 |
| `--seed` | 20260930 | not passed by pipeline (code default) | n/a |  | scripts/run_permutation.py:256 |
| `--sig-cap` | 10.0 | not passed by pipeline (code default) | n/a | Cap for -log10(p) in sig_norm component (default: 10) | scripts/rank_biomarkers.py:537 |
| `--slop` | 10 | 10 [config.consensus.slop] | same | Coordinate tolerance in bp (default: 10) | scripts/annotate_circbase.py:207 |
| `--slop` | 10 | 10 [config.consensus.slop] | same | Coordinate tolerance in bp (default: 10) | scripts/consensus_filter.py:335 |
| `--slop` | 10 | 10 [config.consensus.slop] | same | bp tolerance when matching novel circRNAs across datasets (default 10) | scripts/cross_dataset_biomarkers.py:201 |
| `--stability-min-ratio` | 0.5 | not passed by pipeline (code default) | n/a | A dropped-sample run fails if |log2FC_run| < ratio*|log2FC_full| (default 0.5, i.e. a drop of more than 50%%). Sign flips always fail. | scripts/rank_biomarkers.py:557 |
| `--stats` | {} | not passed by pipeline (code default) | n/a | 統計 JSON，如 '{"n_sig":10,"n_up":6,"n_down":4}' | scripts/notify.py:230 |
| `--status-json` | (not given) | not passed by pipeline (code default) | n/a | Write an explicit run status JSON (status, counts, parameters) here. | scripts/rank_biomarkers.py:573 |
| `--summary` | (required) | not passed by pipeline (code default) | n/a | One or more consensus_summary.tsv files from consensus_filter.py | scripts/annotate_circbase.py:204 |
| `--summary` | (required) | not passed by pipeline (code default) | n/a | Output summary TSV file | scripts/consensus_filter.py:332 |
| `--top-n` | 0 | 50 [config.de.interaction_top_n] | same | Only consider each dataset's top-N ranked candidates (0 = all, default) | scripts/cross_dataset_biomarkers.py:203 |
| `--top-n` | 50 | 50 [config.de.interaction_top_n] | same |  | scripts/predict_interactions.py:807 |
| `--use-pvalue` | (not given) | not passed by pipeline (code default) | n/a | (not static) | scripts/rank_biomarkers.py:536 |
| `config.circbase_file` |  | varies:  (11); /home3/choukaihsuan/circRNA_agent/reference/circbase_hg19.txt (9) | DIFF | leave empty for auto-download; or set local path | config.yaml:85 |
| `config.ciriquant_config` | config/ciriquant.yaml | config/ciriquant.yaml | same |  | config.yaml:1 |
| `config.consensus.adaptive` | True | True | same | default True matches what every published dataset actually ran with | config.yaml:31 |
| `config.consensus.adaptive_ratio` | 0.1 | 0.1 | same |  | config.yaml:34 |
| `config.consensus.max_junction_ratio` | 1.0 | 1.0 | same |  | config.yaml:29 |
| `config.consensus.min_bsj_reads` | 2 | 2 | same |  | config.yaml:26 |
| `config.consensus.min_confidence` | 0.0 | 0.0 | same | consensus_filter.py --min-confidence; previously never wired into circrna.smk | config.yaml:35 |
| `config.consensus.min_tools` | 2 | varies: 2 (17); 1 (3) | DIFF |  | config.yaml:27 |
| `config.consensus.qc_bsj_threshold` | 5 | 5 | same | only apply ratio QC when BSJ < this; high-BSJ circRNAs exempt | config.yaml:30 |
| `config.consensus.slop` | 10 | 10 | same |  | config.yaml:28 |
| `config.consensus.tools` | [ciriquant, dcc] | varies: ['ciriquant', 'dcc'] (17); ['ciriquant'] (3) | DIFF |  | config.yaml:11 |
| `config.de.biomarker.components.conf.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.conf.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.components.fc.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.fc.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.components.known.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.known.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.components.mirna.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.mirna.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.components.rbp.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.rbp.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.components.sig.enabled` | True | True | same |  | config.yaml:0 |
| `config.de.biomarker.components.sig.weight` | 1.0 | 1.0 | same |  | config.yaml:0 |
| `config.de.biomarker.fc_cap` | 5 | 5 | same | |log2FC| cap in rank_biomarkers.py fc_norm; was hardcoded | config.yaml:59 |
| `config.de.biomarker.sig_cap` | 10 | 10 | same | -log10(p) cap in rank_biomarkers.py sig_norm; was hardcoded | config.yaml:58 |
| `config.de.de_sig_by` | pvalue | pvalue | same | pvalue = nominal p<0.05 (paper §2.6/§2.7); auto = q<0.2 (Storey) then p<0.05 | config.yaml:43 |
| `config.de.delta_iui_cutoff` | 0.1 | 0.1 | same | minimum |ΔIUI| to call switching | config.yaml:46 |
| `config.de.deseq2_size_factor_type` | poscounts | poscounts | same | estimateSizeFactors(type=); was unconditionally hardcoded in analysis.R | config.yaml:55 |
| `config.de.encori_clip_exp_num` | 1 | 1 | same | min CLIP experiments (1=all, 2=moderate, 3=high confidence) | config.yaml:49 |
| `config.de.encori_program_num` | 2 | 2 | same | min prediction programs agreeing | config.yaml:50 |
| `config.de.fdr_cutoff` | 0.05 | 0.05 | same |  | config.yaml:41 |
| `config.de.filter_expr_min_count` | 5 | 5 | same | edgeR filterByExpr min.count; was hardcoded in analysis.R | config.yaml:53 |
| `config.de.filter_expr_min_total_count` | 15 | 15 | same | edgeR filterByExpr min.total.count (edgeR's own implicit default, now explicit) | config.yaml:54 |
| `config.de.fsj_concordance_lfc` | 0.0 | 0.0 | same | min |logFC_fsj| for Type II concordance; 0=direction-only (recommended) | config.yaml:44 |
| `config.de.heatmap_top_n` | 10 | 10 | same | top N up + N down circRNAs shown in heatmap | config.yaml:47 |
| `config.de.interaction_top_n` | 50 | 50 | same | top N circRNAs to query ENCORI for interactions | config.yaml:48 |
| `config.de.isoform_fdr_cutoff` | 0.1 | 0.1 | same | within-gene FDR for isoform switching (less strict than global) | config.yaml:45 |
| `config.de.log2fc_cutoff` | 1.0 | 1.0 | same |  | config.yaml:42 |
| `config.de.method` | edgeR_ciriquant | varies: edgeR_ciriquant (19); limma (1) | DIFF |  | config.yaml:40 |
| `config.de.normal_label` | normal | normal | same |  | config.yaml:51 |
| `config.de.paired_design_override` | auto | auto | same | auto = detect patient_id column; force_paired/force_unpaired to override | config.yaml:56 |
| `config.de.tumor_label` | tumor | tumor | same |  | config.yaml:52 |
| `config.download.ascp_speed` | 500m | 500m | same |  | config.yaml:77 |
| `config.download.max_concurrent` | 2 | 2 | same |  | config.yaml:73 |
| `config.download.retry` | 3 | 3 | same |  | config.yaml:74 |
| `config.download.sra_cache_dir` | /home/choukaihsuan/sra_cache | varies: /home3/choukaihsuan/GSE77509/sra_cache (2); /home3/choukaihsuan/GSE323364/sra_cache (2); /home3/choukaihsuan/GSE108735/sra_cache (1); /home3/choukaihsuan/GSE113230/sra_cache (1); /home3/choukaihsuan/GSE121842/sra_cache (1); /home3/choukaihsuan/GSE130078/sra_cache (1); /home3/choukaihsuan/GSE133998/sra_cache (1); /home3/choukaihsuan/GSE136569/sra_cache (1); /home3/choukaihsuan/GSE143797/sra_cache (1); /home3/choukaihsuan/GSE148036/sra_cache (1); /home3/choukaihsuan/GSE171011/sra_cache (1); /home3/choukaihsuan/GSE192410/sra_cache (1); /home3/choukaihsuan/GSE192849/sra_cache (1); /home3/choukaihsuan/GSE229705/sra_cache (1); /home3/choukaihsuan/GSE248612/sra_cache (1); /home3/choukaihsuan/sra_cache (1); /home3/choukaihsuan/GSE97239/sra_cache (1); /home3/choukaihsuan/PRJNA553289/sra_cache (1) | DIFF |  | config.yaml:75 |
| `config.download.tmp_dir` | /home/choukaihsuan/sra_tmp | varies: /home3/choukaihsuan/GSE323364/sra_tmp (4); /home3/choukaihsuan/GSE108735/sra_tmp (1); /home3/choukaihsuan/GSE113230/sra_tmp (1); /home3/choukaihsuan/GSE121842/sra_tmp (1); /home3/choukaihsuan/GSE130078/sra_tmp (1); /home3/choukaihsuan/GSE133998/sra_tmp (1); /home3/choukaihsuan/GSE136569/sra_tmp (1); /home3/choukaihsuan/GSE143797/sra_tmp (1); /home3/choukaihsuan/GSE148036/sra_tmp (1); /home3/choukaihsuan/GSE171011/sra_tmp (1); /home3/choukaihsuan/GSE192410/sra_tmp (1); /home3/choukaihsuan/GSE192849/sra_tmp (1); /home3/choukaihsuan/GSE229705/sra_tmp (1); /home3/choukaihsuan/GSE248612/sra_tmp (1); /home3/choukaihsuan/sra_tmp (1); /home3/choukaihsuan/GSE97239/sra_tmp (1); /home3/choukaihsuan/PRJNA553289/sra_tmp (1) | DIFF |  | config.yaml:76 |
| `config.genome.bwa_index` | /mnt/d/ref/hg19/bwa_index/hg19 | /home3/choukaihsuan/reference/hg19/hg19.fa | same |  | config.yaml:79 |
| `config.genome.fasta` | /mnt/d/ref/hg19/hg19.fa | /home3/choukaihsuan/reference/hg19/hg19.fa | same |  | config.yaml:80 |
| `config.genome.gtf` | /mnt/d/ref/hg19/hg19.gtf | /home3/choukaihsuan/reference/hg19/genes.gtf | same |  | config.yaml:81 |
| `config.genome.hisat2_index` | /mnt/d/ref/hg19/hisat2_index/hg19 | varies: /home3/choukaihsuan/reference/hg19/hisat2_index/hg19 (11); /home3/choukaihsuan/reference/hg19/hg19_hisat2_index (9) | DIFF |  | config.yaml:82 |
| `config.genome.species` | hg19 | hg19 | same |  | config.yaml:83 |
| `config.genome.star_index` | /mnt/d/ref/hg19/star_index | /home3/choukaihsuan/reference/hg19/star_index | same |  | config.yaml:84 |
| `config.groups` | metadata/sample_groups.csv | varies: metadata/GSE108735/sample_groups.csv (1); metadata/GSE113230/sample_groups.csv (1); metadata/GSE121842/sample_groups.csv (1); metadata/GSE130078/sample_groups.csv (1); metadata/GSE133998/sample_groups.csv (1); metadata/GSE136569/sample_groups.csv (1); metadata/GSE143797/sample_groups.csv (1); metadata/GSE148036/sample_groups.csv (1); metadata/GSE171011/sample_groups.csv (1); metadata/GSE192410/sample_groups.csv (1); metadata/GSE192849/sample_groups.csv (1); metadata/GSE221107/sample_groups.csv (1); metadata/GSE229705/sample_groups.csv (1); metadata/GSE248612/sample_groups.csv (1); metadata/GSE323364/sample_groups.csv (1); metadata/GSE58135/sample_groups.csv (1); metadata/GSE77509/sample_groups.csv (1); metadata/GSE97239/sample_groups.csv (1); metadata/PRJNA553289/sample_groups.csv (1); metadata/SRP156355/sample_groups.csv (1) | DIFF |  | config.yaml:86 |
| `config.log_dir` | logs | logs | same |  | config.yaml:87 |
| `config.metadata` | metadata/library_info.csv | varies: metadata/GSE108735/library_info.csv (1); metadata/GSE113230/library_info.csv (1); metadata/GSE121842/library_info.csv (1); metadata/GSE130078/library_info.csv (1); metadata/GSE133998/library_info.csv (1); metadata/GSE136569/library_info.csv (1); metadata/GSE143797/library_info.csv (1); metadata/GSE148036/library_info.csv (1); metadata/GSE171011/library_info.csv (1); metadata/GSE192410/library_info.csv (1); metadata/GSE192849/library_info.csv (1); metadata/GSE221107/library_info.csv (1); metadata/GSE229705/library_info.csv (1); metadata/GSE248612/library_info.csv (1); metadata/GSE323364/library_info.csv (1); metadata/GSE58135/library_info.csv (1); metadata/GSE77509/library_info.csv (1); metadata/GSE97239/library_info.csv (1); metadata/PRJNA553289/library_info.csv (1); metadata/SRP156355/library_info.csv (1) | DIFF |  | config.yaml:88 |
| `config.notify.email_to` | <redacted> | <redacted: personal address> | DIFF | key present in some published projects only | config.yaml:0 |
| `config.project_id` | GSE113230 | varies: GSE108735 (1); GSE113230 (1); GSE121842 (1); GSE130078 (1); GSE133998 (1); GSE136569 (1); GSE143797 (1); GSE148036 (1); GSE171011 (1); GSE192410 (1); GSE192849 (1); GSE221107 (1); GSE229705 (1); GSE248612 (1); GSE323364 (1); GSE58135 (1); GSE77509 (1); GSE97239 (1); PRJNA553289 (1); SRP156355 (1) | DIFF |  | config.yaml:89 |
| `config.raw_dir` | /mnt/d/circRNA_data/raw_fastq | varies: /home3/choukaihsuan/GSE108735/raw (1); /home3/choukaihsuan/GSE113230/raw (1); /home3/choukaihsuan/GSE121842/raw (1); /home3/choukaihsuan/GSE130078/raw (1); /home3/choukaihsuan/GSE133998/raw (1); /home3/choukaihsuan/GSE136569/raw (1); /home3/choukaihsuan/GSE143797/raw (1); /home3/choukaihsuan/GSE148036/raw (1); /home3/choukaihsuan/GSE171011/raw (1); /home3/choukaihsuan/GSE192410/raw (1); /home3/choukaihsuan/GSE192849/raw (1); /home3/choukaihsuan/GSE221107/raw (1); /home3/choukaihsuan/GSE229705/raw (1); /home3/choukaihsuan/GSE248612/raw (1); /home3/choukaihsuan/GSE323364/raw (1); /home3/choukaihsuan/GSE58135/raw (1); /home3/choukaihsuan/GSE77509/raw (1); /home3/choukaihsuan/GSE97239/raw (1); /home3/choukaihsuan/PRJNA553289/raw (1); /home3/choukaihsuan/SRP156355/raw (1) | DIFF |  | config.yaml:90 |
| `config.results_dir` | /mnt/d/circRNA_results/GSE113230 | varies: /home3/choukaihsuan/GSE108735_results (1); /home3/choukaihsuan/GSE113230_results (1); /home3/choukaihsuan/GSE121842_results (1); /home3/choukaihsuan/GSE130078_results (1); /home3/choukaihsuan/GSE133998_results (1); /home3/choukaihsuan/GSE136569_results (1); /home3/choukaihsuan/GSE143797_results (1); /home3/choukaihsuan/GSE148036_results (1); /home3/choukaihsuan/GSE171011_results (1); /home3/choukaihsuan/GSE192410_results (1); /home3/choukaihsuan/GSE192849_results (1); /home3/choukaihsuan/GSE221107_results (1); /home3/choukaihsuan/GSE229705_results (1); /home3/choukaihsuan/GSE248612_results (1); /home3/choukaihsuan/GSE323364_results (1); /home3/choukaihsuan/GSE58135_results (1); /home3/choukaihsuan/GSE77509_results (1); /home3/choukaihsuan/GSE97239_results (1); /home3/choukaihsuan/PRJNA553289_results (1); /home3/choukaihsuan/SRP156355_results (1) | DIFF |  | config.yaml:91 |
| `config.study_title` | (not in config.yaml) | varies:  (2); ncRNA-Seq profiling of Renal Cell Carcinoma tumor and adjacent normal kidney tissues (GSE108735) (1); High-throughput RNA sequencing on circular RNA profiles of human triple-negative breast cancer and adjacent normal tissues (1); Expression profiles of circular RNAs in human colorectal cancer based on RNA deep sequencing (Wang et al. 2019) breast cancer and adjacent normal tissues (1); A Novel LincRNA HERES Epigenetically Regulates Wnt Signaling Pathway via Interaction with EZH2 in Esophageal Squamous Cell Carcinoma (1); Transcriptomic sequencing of breast cancer (1); circNFIB suppresses lymphatic metastasis of pancreatic cancer breast cancer and adjacent normal tissues (1); Circular RNA expression in Nasopharyngeal Carcinoma (NPC) tumor vs. nasopharyngitis normal tissues (Ding et al. 2020) (1); Lung transcriptional profiling of pulmonary tuberculosis, adenocarcinoma and sarcoidosis patients using RNA-Seq (1); Next Generation Sequencing of Papillary Thyroid Cancer tissue sample and adjacent normal tissue breast cancer and adjacent normal tissues (1); Circular RNA expression profiling in ovarian cancer tumor and adjacent normal tissues (GSE192410) (1); Circular RNA expression profiling in node-positive primary breast cancer tumor and adjacent normal tissue (1); RNA Sequencing of six pairs of gastric cancer (1); Coordinating catalytic and non-canonical functions of EZH2 sensitizes tumors to CAR-T cell therapy (1); Breast Cancer RNA-seq (1); RNA-seq analysis of hepatocellular carcinoma with portal vein tumor thrombosis (1); Circular RNA expression profiling in bladder cancer tumor and adjacent normal tissues (GSE97239) (1); RNA-seq profiling of small cell lung cancer tumor and adjacent normal lung tissues (Nanjing Chest Hospital) (1); RNA-seq of early-stage invasive ductal carcinoma (IDC) of breast, paired tumor and adjacent normal (1) | DIFF | key present in some published projects only | config.yaml:0 |
| `config.threads` | 8 | varies: 12 (9); 18 (6); 24 (3); 8 (2) | DIFF | -T; was hardcoded in circrna.smk. benchmark/config_benchmark.yaml pins this | config.yaml:16 |
| `config.tools.ciriquant.anchor` | 5 | 5 | same | -a/--anchor: min anchor length for junction alignment | config.yaml:14 |
| `config.tools.ciriquant.no_gene` | False | False | same | --no-gene: skip stringtie gene abundance (speed only; no effect on circRNA calls) | config.yaml:13 |
| `config.tools.dcc.count_threshold` | 2 | 2 | same | -Nr arg 1 | config.yaml:18 |
| `config.tools.dcc.end_tol` | 5 | 5 | same | -E | config.yaml:22 |
| `config.tools.dcc.filter_mito` | True | True | same | -M | config.yaml:20 |
| `config.tools.dcc.filter_recommended` | True | True | same | -F | config.yaml:21 |
| `config.tools.dcc.max_length` | 1000000 | 1000000 | same | -m | config.yaml:23 |
| `config.tools.dcc.min_length` | 30 | 30 | same | -n | config.yaml:24 |
| `config.tools.dcc.replicate_threshold` | 1 | 1 | same | -Nr arg 2 | config.yaml:19 |
| `config.tools.dcc.threads` | 4 | 4 | same | -T; was hardcoded in circrna.smk. benchmark/config_benchmark.yaml pins this | config.yaml:16 |
| `config.trimmed_dir` | /mnt/d/circRNA_data/trimmed | varies: /home3/choukaihsuan/GSE108735/trimmed (1); /home3/choukaihsuan/GSE113230/trimmed (1); /home3/choukaihsuan/GSE121842/trimmed (1); /home3/choukaihsuan/GSE130078/trimmed (1); /home3/choukaihsuan/GSE133998/trimmed (1); /home3/choukaihsuan/GSE136569/trimmed (1); /home3/choukaihsuan/GSE143797/trimmed (1); /home3/choukaihsuan/GSE148036/trimmed (1); /home3/choukaihsuan/GSE171011/trimmed (1); /home3/choukaihsuan/GSE192410/trimmed (1); /home3/choukaihsuan/GSE192849/trimmed (1); /home3/choukaihsuan/GSE221107/trimmed (1); /home3/choukaihsuan/GSE229705/trimmed (1); /home3/choukaihsuan/GSE248612/trimmed (1); /home3/choukaihsuan/GSE323364/trimmed (1); /home3/choukaihsuan/GSE58135/trimmed (1); /home3/choukaihsuan/GSE77509/trimmed (1); /home3/choukaihsuan/GSE97239/trimmed (1); /home3/choukaihsuan/PRJNA553289/trimmed (1); /home3/choukaihsuan/SRP156355/trimmed (1) | DIFF |  | config.yaml:93 |
| `config.use_container` | False | False | same | set true when running with --use-singularity | config.yaml:2 |
