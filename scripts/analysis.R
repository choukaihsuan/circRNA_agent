# analysis.R – DESeq2-based differential expression for circRNAs
# Called as a Snakemake Rscript; receives inputs/outputs via the snakemake object.

suppressPackageStartupMessages({
  library(DESeq2)
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

dir.create(dirname(de_out),      recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(volcano_out), recursive = TRUE, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
counts <- read.table(count_file, sep = "\t", header = TRUE, row.names = 1,
                     check.names = FALSE)
groups <- read.csv(groups_file)

# Keep only samples that appear in groups table
common_samples <- intersect(colnames(counts), groups$srr_id)
if (length(common_samples) < 2) stop("Fewer than 2 samples match between count matrix and groups file")

counts <- counts[, common_samples, drop = FALSE]
groups <- groups[groups$srr_id %in% common_samples, ]
groups <- groups[match(colnames(counts), groups$srr_id), ]

# Round to integers (DESeq2 requires integer counts)
counts <- round(counts)
storage.mode(counts) <- "integer"

# Remove circRNAs with zero counts in all samples
counts <- counts[rowSums(counts) > 0, ]

# ── DESeq2 ────────────────────────────────────────────────────────────────────
col_data <- data.frame(
  row.names = groups$srr_id,
  condition = factor(groups$condition, levels = c(normal_label, tumor_label))
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = col_data,
  design    = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", tumor_label, normal_label))
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("circ_id") %>%
  arrange(padj)

write.table(res_df, de_out, sep = "\t", quote = FALSE, row.names = FALSE)
message("[OK] DE results → ", de_out)

# ── Variance-stabilised counts for plots ─────────────────────────────────────
vsd <- vst(dds, blind = FALSE)

# ── Volcano plot ──────────────────────────────────────────────────────────────
plot_df <- res_df %>%
  mutate(
    sig = case_when(
      padj < fdr_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
      padj < fdr_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(padj))

pdf(volcano_out, width = 7, height = 6)
ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(Up = "#d62728", Down = "#1f77b4", NS = "grey70")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(fdr_cutoff),           linetype = "dashed", colour = "grey40") +
  labs(
    title  = paste0("Volcano Plot (", tumor_label, " vs ", normal_label, ")"),
    x      = expression(log[2]~"Fold Change"),
    y      = expression(-log[10]~"(adjusted p-value)"),
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")
dev.off()
message("[OK] Volcano plot → ", volcano_out)

# ── PCA ───────────────────────────────────────────────────────────────────────
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
var_pct  <- round(100 * attr(pca_data, "percentVar"))

pdf(pca_out, width = 6, height = 5)
ggplot(pca_data, aes(x = PC1, y = PC2, colour = condition, label = name)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title  = "PCA – VST-normalised counts",
    x      = paste0("PC1: ", var_pct[1], "% variance"),
    y      = paste0("PC2: ", var_pct[2], "% variance"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 13)
dev.off()
message("[OK] PCA → ", pca_out)

# ── Heatmap (top 50 DE circRNAs) ──────────────────────────────────────────────
top50 <- res_df %>%
  filter(!is.na(padj)) %>%
  slice_min(padj, n = 50) %>%
  pull(circ_id)

mat <- assay(vsd)[top50, , drop = FALSE]
mat <- mat - rowMeans(mat)   # centre for visual clarity

ann_col <- data.frame(
  condition = col_data$condition,
  row.names = rownames(col_data)
)
ann_colours <- list(
  condition = setNames(
    c("#d62728", "#1f77b4"),
    c(tumor_label, normal_label)
  )
)

pdf(heatmap_out, width = 10, height = 10)
pheatmap(
  mat,
  annotation_col  = ann_col,
  annotation_colors = ann_colours,
  color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  show_rownames   = nrow(mat) <= 60,
  fontsize_row    = 7,
  border_color    = NA,
  main            = paste0("Top ", nrow(mat), " DE circRNAs (VST, centred)")
)
dev.off()
message("[OK] Heatmap → ", heatmap_out)
