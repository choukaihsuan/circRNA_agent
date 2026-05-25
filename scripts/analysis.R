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

counts <- round(counts)
storage.mode(counts) <- "integer"
counts <- counts[rowSums(counts) > 0, ]

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
  fsj <- fsj[, common_samples, drop = FALSE]
  fsj <- round(fsj)
  storage.mode(fsj) <- "integer"

  # Align rows: only circRNAs present in both matrices
  common_circs <- intersect(rownames(counts), rownames(fsj))
  if (length(common_circs) == 0)
    stop("No circRNAs in common between BSJ and FSJ matrices")
  bsj <- counts[common_circs, ]
  fsj <- fsj[common_circs, ]

  # ── TMM normalisation on FSJ to get library-size offset ──────────────────
  fsj_safe <- fsj + 1L                  # avoid log(0) for zero-FSJ loci
  fsj_dge  <- DGEList(counts = fsj_safe, group = condition)
  fsj_dge  <- calcNormFactors(fsj_dge, method = "TMM")
  # offset[i,j] = log( effective FSJ library size for sample j )
  eff_lib <- fsj_dge$samples$lib.size * fsj_dge$samples$norm.factors
  offset_mat <- matrix(log(eff_lib), nrow = nrow(bsj), ncol = ncol(bsj), byrow = TRUE)

  # ── BSJ GLM with FSJ offset ───────────────────────────────────────────────
  dge    <- DGEList(counts = bsj, group = condition)
  design <- model.matrix(~ condition)
  dge    <- estimateDisp(dge, design = design, offset = offset_mat)
  fit    <- glmQLFit(dge, design = design, offset = offset_mat)
  qlf    <- glmQLFTest(fit, coef = 2)
  res    <- as.data.frame(topTags(qlf, n = Inf, sort.by = "PValue"))

  # ── Type classification via separate FSJ test ────────────────────────────
  fsj_dge2 <- calcNormFactors(DGEList(counts = fsj_safe, group = condition), method = "TMM")
  fsj_dge2 <- estimateDisp(fsj_dge2, design = design)
  fsj_fit  <- glmQLFit(fsj_dge2, design = design)
  fsj_qlf  <- glmQLFTest(fsj_fit, coef = 2)
  fsj_res  <- as.data.frame(topTags(fsj_qlf, n = Inf))

  # Type I  : ratio changes, FSJ stable  (circRNA-specific)
  # Type II : FSJ also significantly changes in same direction (gene-level)
  type_vec <- ifelse(
    rownames(res) %in% rownames(fsj_res[fsj_res$FDR < 0.05, ]) &
      sign(res$logFC) == sign(fsj_res[rownames(res), "logFC"]),
    "II", "I"
  )

  res_df <- res %>%
    tibble::rownames_to_column("circ_id") %>%
    mutate(Type = type_vec) %>%
    rename(log2FC = logFC, pvalue = PValue, padj = FDR) %>%
    select(circ_id, log2FC, pvalue, padj, Type, everything()) %>%
    arrange(padj)

  # Log-normalised counts for plots (log2 CPM with FSJ offset)
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

plot_df <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(sig = case_when(
    padj < fdr_cutoff & log2FC >  lfc_cutoff ~ "Up",
    padj < fdr_cutoff & log2FC < -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  ))

# ── Volcano ───────────────────────────────────────────────────────────────────
pdf(volcano_out, width = 7, height = 6)
ggplot(plot_df, aes(x = log2FC, y = -log10(padj), colour = sig)) +
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
    y      = expression(-log[10]~"(adjusted p-value)"),
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
