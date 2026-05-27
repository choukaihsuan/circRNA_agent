# analysis.R – circRNA differential expression analysis
#
# Supported methods (config de.method):
#   edgeR_ciriquant  – edgeR GLM + FSJ offset (replicates CIRI_DE_replicate)
#   deseq2           – DESeq2 RLE normalisation
#   limma            – limma-voom

suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

# ── I/O from Snakemake ────────────────────────────────────────────────────────
count_file   <- snakemake@input[["matrix"]]
groups_file  <- snakemake@input[["groups"]]
de_out       <- snakemake@output[["de"]]
volcano_out  <- snakemake@output[["volcano"]]
heatmap_out  <- snakemake@output[["heatmap"]]
pca_out      <- snakemake@output[["pca"]]
fdr_cutoff   <- as.numeric(snakemake@params[["fdr"]])
lfc_cutoff   <- as.numeric(snakemake@params[["lfc"]])
tumor_label  <- snakemake@params[["tumor_label"]]
normal_label <- snakemake@params[["normal_label"]]
use_pvalue   <- isTRUE(tryCatch(snakemake@params[["use_pvalue"]], error = function(e) FALSE))

# Backward-compatible: older DAG may not pass de_method or fsj_matrix
de_method <- tryCatch(snakemake@params[["de_method"]], error = function(e) NULL)
fsj_file  <- tryCatch(snakemake@input[["fsj_matrix"]], error = function(e) NULL)

# Fall back to deseq2 if edgeR_ciriquant is requested but FSJ matrix is unavailable
if (is.null(de_method)) de_method <- "deseq2"
if (de_method == "edgeR_ciriquant" && (is.null(fsj_file) || !file.exists(fsj_file))) {
  message("[WARN] fsj_count_matrix not available; falling back to deseq2")
  de_method <- "deseq2"
}

dir.create(dirname(de_out),      recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(volcano_out), recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
counts <- read.table(count_file, sep = "\t", header = TRUE, row.names = 1,
                     check.names = FALSE)
groups <- read.csv(groups_file)

common_samples <- intersect(colnames(counts), groups$srr_id)
if (length(common_samples) < 2)
  stop("Fewer than 2 samples match between count matrix and groups file")

counts <- counts[, common_samples, drop = FALSE]
groups <- groups[match(common_samples, groups$srr_id), ]

condition <- factor(groups$condition, levels = c(normal_label, tumor_label))

counts <- data.matrix(round(counts))
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

message(sprintf("[DE] method=%s  samples=%d  circRNAs=%d", de_method,
                ncol(counts), nrow(counts)))


# ══════════════════════════════════════════════════════════════════════════════
# edgeR_ciriquant — replicates CIRI_DE_replicate
#   BSJ counts tested via edgeR GLM, with log(FSJ lib size) as offset.
#   This models the BSJ / FSJ ratio change, controlling for host-gene expression.
#   Type I  = circRNA-specific change (FSJ stable)
#   Type II = gene-level change drives BSJ change (FSJ also DE)
# ══════════════════════════════════════════════════════════════════════════════
if (de_method == "edgeR_ciriquant") {
  suppressPackageStartupMessages(library(edgeR))

  fsj <- read.table(fsj_file, sep = "\t", header = TRUE, row.names = 1,
                    check.names = FALSE)
  fsj <- data.matrix(round(fsj[, common_samples, drop = FALSE]))
  storage.mode(fsj) <- "integer"

  # Align rows: only circRNAs present in both matrices
  common_circs <- intersect(rownames(counts), rownames(fsj))
  if (length(common_circs) == 0)
    stop("No circRNAs in common between BSJ and FSJ matrices")
  bsj <- counts[common_circs, ]
  fsj <- fsj[common_circs, ]

  # Min-count filter: require >= 1 read in at least min_samples samples
  min_samp <- min(table(condition))
  keep <- rowSums(bsj >= 1) >= min_samp
  bsj  <- bsj[keep, ]
  fsj  <- fsj[keep, ]
  message(sprintf("[edgeR_v2] %d circRNAs after min-count filter", sum(keep)))

  # ── Per-locus FSJ offset ──────────────────────────────────────────────────
  # Each circRNA uses its own FSJ counts as offset, so the GLM tests the
  # BSJ/FSJ ratio change rather than absolute BSJ abundance.
  fsj_pseudo <- fsj + 0.5
  y_fsj_norm <- DGEList(counts = fsj, group = condition)
  y_fsj_norm <- calcNormFactors(y_fsj_norm, method = "TMM")
  fsj_scaled <- t(t(fsj_pseudo) * y_fsj_norm$samples$norm.factors)
  fsj_cpm    <- t(t(fsj_scaled) / colSums(fsj_scaled) * 1e6)
  offset_mat <- log(fsj_cpm + 0.5)   # n_circ × n_sample

  # ── BSJ GLM with per-locus FSJ offset ────────────────────────────────────
  design <- model.matrix(~ condition)
  dge    <- DGEList(counts = bsj, group = condition)
  dge$offset <- offset_mat
  dge    <- estimateDisp(dge, design, robust = TRUE)
  fit    <- glmQLFit(dge, design, robust = TRUE)
  qlf    <- glmQLFTest(fit, coef = 2)
  res_bsj <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))

  # ── Independent FSJ test (host gene expression) ──────────────────────────
  fsj_safe <- fsj + 1L
  dge_fsj  <- DGEList(counts = fsj_safe, group = condition)
  dge_fsj  <- calcNormFactors(dge_fsj, method = "TMM")
  dge_fsj  <- estimateDisp(dge_fsj, design, robust = TRUE)
  fit_fsj  <- glmQLFit(dge_fsj, design, robust = TRUE)
  qlf_fsj  <- glmQLFTest(fit_fsj, coef = 2)
  res_fsj  <- as.data.frame(topTags(qlf_fsj, n = Inf, sort.by = "none"))

  # ── Merge BSJ + FSJ results ───────────────────────────────────────────────
  res_bsj$circ_id <- rownames(res_bsj)
  res_fsj$circ_id <- rownames(res_fsj)
  # Include PValue from both tables so FSJ nominal p-value is available.
  # With PValue in both tables it becomes PValue_bsj / PValue_fsj after merge.
  res <- merge(
    res_bsj[, c("circ_id", "logFC", "logCPM", "PValue", "FDR")],
    res_fsj[, c("circ_id", "logFC", "PValue", "FDR")],
    by = "circ_id", all.x = TRUE, suffixes = c("_bsj", "_fsj")
  )

  # ── Type I / II / III classification ─────────────────────────────────────
  # Type_I   : BSJ/FSJ ratio changes; FSJ stable (circRNA-specific switching)
  # Type_II  : BSJ/FSJ ratio changes; FSJ also DE in same direction (gene-level)
  # Type_III : FSJ DE but BSJ/FSJ ratio not significant (host gene only)
  # NS       : neither significant
  bsj_sig_col <- if (use_pvalue) "PValue_bsj" else "FDR_bsj"
  fsj_sig_col <- if (use_pvalue) "PValue_fsj" else "FDR_fsj"
  res$sig_bsj <- res[[bsj_sig_col]] < fdr_cutoff & abs(res$logFC_bsj) >= lfc_cutoff
  res$sig_fsj <- !is.na(res[[fsj_sig_col]]) & res[[fsj_sig_col]] < fdr_cutoff
  # concordant requires same direction AND FSJ |logFC| >= 0.5 to avoid
  # noise when both tests are marginally significant
  # concordant: FSJ changes in the same direction as BSJ.
  # |logFC_fsj| threshold removed — FSJ significance (fsj_sig_col) already
  # filters noise; requiring an additional FC floor was too conservative.
  res$concordant <- with(res,
    !is.na(logFC_fsj) &
    sign(logFC_bsj) == sign(logFC_fsj)
  )
  res$Type <- dplyr::case_when(
    res$sig_bsj & !res$sig_fsj                   ~ "Type_I",
    res$sig_bsj & res$sig_fsj & res$concordant   ~ "Type_II",
    res$sig_bsj & res$sig_fsj & !res$concordant  ~ "Type_I",
    !res$sig_bsj & res$sig_fsj                   ~ "Type_III",
    TRUE                                          ~ "NS"
  )

  # ── Circular Splicing Index (CSI = BSJ / (BSJ + FSJ + 1)) ────────────────
  tumor_idx  <- which(condition == tumor_label)
  normal_idx <- which(condition == normal_label)
  csi_fn <- function(idx)
    rowMeans(bsj[, idx, drop = FALSE] /
             (bsj[, idx, drop = FALSE] + fsj[, idx, drop = FALSE] + 1))
  csi_df <- data.frame(
    circ_id    = rownames(bsj),
    csi_tumor  = csi_fn(tumor_idx),
    csi_normal = csi_fn(normal_idx)
  )
  csi_df$delta_csi <- csi_df$csi_tumor - csi_df$csi_normal
  res <- merge(res, csi_df, by = "circ_id", all.x = TRUE)

  # ── Assemble final res_df ─────────────────────────────────────────────────
  res_df <- res %>%
    rename(log2FC = logFC_bsj, pvalue = PValue_bsj, padj = FDR_bsj,
           pvalue_fsj = PValue_fsj) %>%
    select(circ_id, log2FC, pvalue, padj, Type,
           logFC_fsj, pvalue_fsj, FDR_fsj, delta_csi, csi_tumor, csi_normal, logCPM) %>%
    arrange(padj)

  log_cpm <- cpm(dge, log = TRUE, offset = offset_mat)

} else if (de_method == "deseq2") {

  # ══════════════════════════════════════════════════════════════════════════
  # DESeq2 — RLE normalisation, Wald test
  # ══════════════════════════════════════════════════════════════════════════
  suppressPackageStartupMessages(library(DESeq2))

  col_data <- data.frame(row.names = common_samples, condition = condition)
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data,
                                 design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", tumor_label, normal_label))
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("circ_id") %>%
    rename(log2FC = log2FoldChange) %>%
    arrange(padj)

  vsd     <- vst(dds, blind = FALSE)
  log_cpm <- assay(vsd)

} else if (de_method == "limma") {

  # ══════════════════════════════════════════════════════════════════════════
  # limma-voom
  # ══════════════════════════════════════════════════════════════════════════
  suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
  })

  dge    <- DGEList(counts = counts, group = condition)
  dge    <- calcNormFactors(dge, method = "TMM")
  design <- model.matrix(~ condition)
  v      <- voom(dge, design, plot = FALSE)
  fit    <- lmFit(v, design)
  fit    <- eBayes(fit)
  res    <- topTable(fit, coef = 2, n = Inf, sort.by = "P")
  res_df <- res %>%
    tibble::rownames_to_column("circ_id") %>%
    rename(log2FC = logFC, pvalue = P.Value, padj = adj.P.Val) %>%
    arrange(padj)

  log_cpm <- v$E

} else {
  stop(paste("Unknown de_method:", de_method,
             "— must be edgeR_ciriquant, deseq2, or limma"))
}

write.table(res_df, de_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("[OK] DE results (", nrow(res_df), " circRNAs) → ", de_out)


# ── Shared plot helpers ───────────────────────────────────────────────────────
sig_col <- c(Up = "#d62728", Down = "#1f77b4", NS = "grey70")

plot_sig_col <- if (use_pvalue) "pvalue" else "padj"
plot_y_label <- if (use_pvalue) expression(-log[10]~"(p-value, nominal)") else expression(-log[10]~"(adjusted p-value)")

plot_df <- res_df %>%
  filter(!is.na(.data[[plot_sig_col]])) %>%
  mutate(sig = case_when(
    .data[[plot_sig_col]] < fdr_cutoff & log2FC >  lfc_cutoff ~ "Up",
    .data[[plot_sig_col]] < fdr_cutoff & log2FC < -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  ))

# ── Volcano ───────────────────────────────────────────────────────────────────
pdf(volcano_out, width = 7, height = 6)
ggplot(plot_df, aes(x = log2FC, y = -log10(.data[[plot_sig_col]]), colour = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = sig_col) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed",
             colour = "grey40") +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
             colour = "grey40") +
  labs(
    title  = paste0("Volcano (", tumor_label, " vs ", normal_label,
                    ")  [", de_method, "]"),
    x      = expression(log[2]~"Fold Change"),
    y      = plot_y_label,
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")
dev.off()
message("[OK] Volcano → ", volcano_out)

# ── PCA ───────────────────────────────────────────────────────────────────────
pca_res  <- prcomp(t(log_cpm), scale. = FALSE)
var_pct  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2))
pca_data <- data.frame(
  PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
  condition = condition,
  name = common_samples
)

pdf(pca_out, width = 6, height = 5)
ggplot(pca_data, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title  = paste0("PCA  [", de_method, "]"),
    x      = paste0("PC1: ", var_pct[1], "%"),
    y      = paste0("PC2: ", var_pct[2], "%"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 13)
dev.off()
message("[OK] PCA → ", pca_out)

# ── Heatmap (top 50 by padj) ─────────────────────────────────────────────────
top_ids <- res_df %>%
  filter(!is.na(padj)) %>%
  slice_min(padj, n = 50) %>%
  pull(circ_id)

mat <- log_cpm[intersect(top_ids, rownames(log_cpm)), , drop = FALSE]
mat <- mat - rowMeans(mat)

ann_col <- data.frame(condition = condition, row.names = common_samples)
ann_col_colours <- list(
  condition = setNames(c("#d62728", "#1f77b4"), c(tumor_label, normal_label))
)

pdf(heatmap_out, width = 10, height = 10)
pheatmap(
  mat,
  annotation_col    = ann_col,
  annotation_colors = ann_col_colours,
  color             = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  show_rownames     = nrow(mat) <= 60,
  fontsize_row      = 7,
  border_color      = NA,
  main              = paste0("Top ", nrow(mat), " DE circRNAs  [", de_method, "]")
)
dev.off()
message("[OK] Heatmap → ", heatmap_out)
