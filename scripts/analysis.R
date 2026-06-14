# analysis.R – circRNA differential expression analysis
#
# Always runs ALL THREE methods (edgeR_ciriquant, DESeq2, limma-voom).
# Outputs:
#   de_results.tsv                      – primary method (config de.method)
#   de_results_edgeR_ciriquant.tsv      – edgeR/FSJ offset  (if output declared)
#   de_results_deseq2.tsv               – DESeq2            (if output declared)
#   de_results_limma.tsv                – limma-voom        (if output declared)
#   volcano.pdf / heatmap.pdf / pca.pdf – primary method only

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

# Optional per-method outputs (may not be declared in all rules)
de_edger_out <- tryCatch(snakemake@output[["de_edger"]], error = function(e) NULL)
de_deseq_out <- tryCatch(snakemake@output[["de_deseq"]], error = function(e) NULL)
de_limma_out <- tryCatch(snakemake@output[["de_limma"]], error = function(e) NULL)

fdr_cutoff   <- as.numeric(snakemake@params[["fdr"]])
lfc_cutoff   <- as.numeric(snakemake@params[["lfc"]])
tumor_label  <- snakemake@params[["tumor_label"]]
normal_label <- snakemake@params[["normal_label"]]
de_sig_by    <- tryCatch(snakemake@params[["de_sig_by"]], error = function(e) "auto")
if (!is.character(de_sig_by) || !de_sig_by %in% c("auto", "pvalue", "padj", "qvalue"))
  de_sig_by <- if (isTRUE(de_sig_by)) "pvalue" else "auto"
# Minimum |logFC_fsj| required for concordance (Type II classification).
# 0.0 = direction match only; 0.5 = old strict setting.
# CIRIquant FSJ counts are junction-specific (not whole-gene CPM), so
# fold-changes tend to be small even when the host gene genuinely changes.
fsj_lfc_thr  <- tryCatch(as.numeric(snakemake@params[["fsj_concordance_lfc"]]),
                          error = function(e) 0.0)
heatmap_top_n_raw <- tryCatch(snakemake@params[["heatmap_top_n"]], error = function(e) 10L)
heatmap_top_n <- if (length(heatmap_top_n_raw) == 0 || is.null(heatmap_top_n_raw)) 10L else as.integer(heatmap_top_n_raw)
circbase_annot_file <- tryCatch(snakemake@input[["circbase_annot"]], error = function(e) NULL)

# Resolve primary method
de_method <- tryCatch(snakemake@params[["de_method"]], error = function(e) "deseq2")
if (is.null(de_method)) de_method <- "deseq2"
fsj_file  <- tryCatch(snakemake@input[["fsj_matrix"]], error = function(e) NULL)

# Fallback: if edgeR requested but FSJ unavailable, warn (still try, guard inside block)
if (de_method == "edgeR_ciriquant" && (is.null(fsj_file) || !file.exists(fsj_file))) {
  message("[WARN] fsj_count_matrix not available; primary method changed to deseq2")
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

has_patient <- "patient_id" %in% colnames(groups)
if (has_patient) {
  patient <- factor(groups$patient_id)
  design  <- model.matrix(~ patient + condition)
  message("[DE] Paired design: ~ patient + condition")
} else {
  design <- model.matrix(~ condition)
  message("[DE] Unpaired design: ~ condition")
}
# condition coefficient is always the last column
cond_coef <- ncol(design)

message(sprintf("[DE] primary=%s  samples=%d  circRNAs=%d",
                de_method, ncol(counts), nrow(counts)))


# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

# Cascade significance: auto → Storey q (if available + finds hits) → nominal p
do_cascade <- function(res_df, de_sig_by, fdr_cutoff) {
  eff_col <- "pvalue"; eff_thr <- fdr_cutoff
  if (de_sig_by == "auto") {
    q_result <- tryCatch({
      suppressPackageStartupMessages(library(qvalue))
      pv_vec <- setNames(res_df$pvalue, res_df$circ_id)
      pv_use <- pv_vec[!is.na(pv_vec) & pv_vec > 0 & pv_vec <= 1]
      q_out  <- qvalue(pv_use)
      list(qvalues = setNames(q_out$qvalues, names(pv_use)), pi0 = q_out$pi0)
    }, error = function(e) {
      message("[cascade] qvalue failed: ", conditionMessage(e)); NULL
    })
    if (!is.null(q_result)) {
      res_df$qvalue <- q_result$qvalues[res_df$circ_id]
      n_sig <- sum(!is.na(res_df$qvalue) & res_df$qvalue < 0.2)
      message(sprintf("[cascade] q < 0.2: %d significant (pi0=%.3f)", n_sig, q_result$pi0))
      if (n_sig > 0) { eff_col <- "qvalue"; eff_thr <- 0.2 } else {
        message("[cascade] q < 0.2 found nothing — fallback to nominal p")
      }
    }
  } else if (de_sig_by == "pvalue") {
    eff_col <- "pvalue"; eff_thr <- fdr_cutoff
  } else if (de_sig_by == "qvalue") {
    eff_col <- "qvalue"; eff_thr <- 0.2
  } else {
    eff_col <- "padj"; eff_thr <- fdr_cutoff
  }
  message(sprintf("[cascade] significance: %s < %.3g", eff_col, eff_thr))
  list(res_df = res_df, eff_col = eff_col, eff_thr = eff_thr)
}

# Type I / II / III classification (edgeR_ciriquant only; NA for all others)
add_type_col <- function(res_df, method_name, eff_col, eff_thr, fdr_cutoff, lfc_cutoff,
                         fsj_lfc_thr = 0.0) {
  if (method_name == "edgeR_ciriquant" && "logFC_fsj" %in% colnames(res_df)) {
    fsj_sig_col <- if (eff_col == "padj") "FDR_fsj" else "pvalue_fsj"
    sig_bsj  <- !is.na(res_df[[eff_col]]) & res_df[[eff_col]] < eff_thr & abs(res_df$log2FC) >= lfc_cutoff
    # sig_fsj: FSJ independently significant AND fold-change >= fsj_lfc_thr.
    # Note: direction concordance (sign check) is NOT used here because the BSJ/FSJ
    # ratio test (offset model) creates an inherent anti-correlation between log2FC
    # direction and logFC_fsj direction — checking sign equality produces near-zero
    # Type II. Instead, Type II = BOTH ratio-test AND independent FSJ-test significant,
    # regardless of direction; the offset already removes pure-FSJ-driven ratio changes.
    sig_fsj  <- !is.na(res_df[[fsj_sig_col]]) &
                res_df[[fsj_sig_col]] < fdr_cutoff &
                (!is.na(res_df[["logFC_fsj"]]) & abs(res_df[["logFC_fsj"]]) >= fsj_lfc_thr)
    res_df$Type <- dplyr::case_when(
      sig_bsj & sig_fsj    ~ "Type_II",
      sig_bsj & !sig_fsj   ~ "Type_I",
      !sig_bsj & sig_fsj   ~ "Type_III",
      TRUE                 ~ "NS"
    )
  } else {
    res_df$Type <- NA_character_
  }
  res_df
}


# ══════════════════════════════════════════════════════════════════════════════
# Method 1 – edgeR_ciriquant
# ══════════════════════════════════════════════════════════════════════════════
result_edger <- tryCatch({
  if (is.null(fsj_file) || !file.exists(fsj_file))
    stop("FSJ matrix not available")
  suppressPackageStartupMessages(library(edgeR))

  fsj <- read.table(fsj_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
  fsj <- data.matrix(round(fsj[, common_samples, drop = FALSE]))
  storage.mode(fsj) <- "integer"

  common_circs <- intersect(rownames(counts), rownames(fsj))
  if (length(common_circs) == 0) stop("No circRNAs in common between BSJ and FSJ matrices")
  bsj <- counts[common_circs, ]; fsj <- fsj[common_circs, ]

  keep <- filterByExpr(bsj, group = condition, min.count = 5)
  bsj  <- bsj[keep, ]; fsj <- fsj[keep, ]
  message(sprintf("[edgeR_ciriquant] %d circRNAs after filterByExpr", sum(keep)))

  # Per-locus FSJ offset
  fsj_pseudo <- fsj + 0.5
  y_fsj_norm <- DGEList(counts = fsj, group = condition)
  y_fsj_norm <- calcNormFactors(y_fsj_norm, method = "TMM")
  fsj_scaled <- t(t(fsj_pseudo) * y_fsj_norm$samples$norm.factors)
  fsj_cpm    <- t(t(fsj_scaled) / colSums(fsj_scaled) * 1e6)
  offset_mat <- log(fsj_cpm + 0.5)

  # BSJ GLM
  dge    <- DGEList(counts = bsj, group = condition)
  dge$offset <- offset_mat
  dge    <- estimateDisp(dge, design, robust = TRUE)
  fit    <- glmQLFit(dge, design, robust = TRUE)
  qlf    <- glmQLFTest(fit, coef = cond_coef)
  res_bsj <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))

  # Independent FSJ test
  fsj_safe <- fsj + 1L
  dge_fsj  <- DGEList(counts = fsj_safe, group = condition)
  dge_fsj  <- calcNormFactors(dge_fsj, method = "TMM")
  dge_fsj  <- estimateDisp(dge_fsj, design, robust = TRUE)
  fit_fsj  <- glmQLFit(dge_fsj, design, robust = TRUE)
  qlf_fsj  <- glmQLFTest(fit_fsj, coef = cond_coef)
  res_fsj  <- as.data.frame(topTags(qlf_fsj, n = Inf, sort.by = "none"))

  res_bsj$circ_id <- rownames(res_bsj)
  res_fsj$circ_id <- rownames(res_fsj)
  res <- merge(
    res_bsj[, c("circ_id", "logFC", "logCPM", "PValue", "FDR")],
    res_fsj[, c("circ_id", "logFC", "PValue", "FDR")],
    by = "circ_id", all.x = TRUE, suffixes = c("_bsj", "_fsj")
  )

  # CSI
  tumor_idx  <- which(condition == tumor_label)
  normal_idx <- which(condition == normal_label)
  csi_fn <- function(idx)
    rowMeans(bsj[, idx, drop = FALSE] / (bsj[, idx, drop = FALSE] + fsj[, idx, drop = FALSE] + 1))
  csi_case_col    <- paste0("csi_", tumor_label)
  csi_control_col <- paste0("csi_", normal_label)
  csi_df <- data.frame(circ_id = rownames(bsj), stringsAsFactors = FALSE)
  csi_df[[csi_case_col]]    <- csi_fn(tumor_idx)
  csi_df[[csi_control_col]] <- csi_fn(normal_idx)
  csi_df$delta_csi <- csi_df[[csi_case_col]] - csi_df[[csi_control_col]]
  res <- merge(res, csi_df, by = "circ_id", all.x = TRUE)

  res_df <- res %>%
    rename(log2FC = logFC_bsj, pvalue = PValue_bsj, padj = FDR_bsj,
           pvalue_fsj = PValue_fsj) %>%
    select(circ_id, log2FC, pvalue, padj,
           logFC_fsj, pvalue_fsj, FDR_fsj, delta_csi,
           all_of(c(csi_case_col, csi_control_col)), logCPM) %>%
    arrange(padj)

  log_cpm <- cpm(dge, log = TRUE, offset = offset_mat)
  list(res_df = res_df, log_cpm = log_cpm, success = TRUE)
}, error = function(e) {
  message("[WARN] edgeR_ciriquant: ", conditionMessage(e))
  list(success = FALSE)
})


# ══════════════════════════════════════════════════════════════════════════════
# Method 2 – DESeq2
# ══════════════════════════════════════════════════════════════════════════════
result_deseq <- tryCatch({
  suppressPackageStartupMessages(library(DESeq2))
  col_data <- data.frame(row.names = common_samples, condition = condition)
  if (has_patient) col_data$patient <- patient
  deseq_design <- if (has_patient) ~ patient + condition else ~ condition
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = deseq_design)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", tumor_label, normal_label))
  res_df <- as.data.frame(res)
  names(res_df)[names(res_df) == "log2FoldChange"] <- "log2FC"
  res_df <- res_df %>% tibble::rownames_to_column("circ_id") %>% arrange(padj)
  vsd <- tryCatch(
    vst(dds, blind = FALSE),
    error = function(e) varianceStabilizingTransformation(dds, blind = FALSE)
  )
  list(res_df = res_df, log_cpm = assay(vsd), success = TRUE)
}, error = function(e) {
  message("[WARN] DESeq2: ", conditionMessage(e))
  list(success = FALSE)
})


# ══════════════════════════════════════════════════════════════════════════════
# Method 3 – limma-voom
# ══════════════════════════════════════════════════════════════════════════════
result_limma <- tryCatch({
  suppressPackageStartupMessages({ library(limma); library(edgeR) })
  dge    <- DGEList(counts = counts, group = condition)
  dge    <- calcNormFactors(dge, method = "TMM")
  v      <- voom(dge, design, plot = FALSE)
  fit    <- lmFit(v, design)
  fit    <- eBayes(fit)
  res    <- topTable(fit, coef = cond_coef, n = Inf, sort.by = "P")
  res_df <- res %>% tibble::rownames_to_column("circ_id")
  names(res_df)[names(res_df) == "logFC"]     <- "log2FC"
  names(res_df)[names(res_df) == "P.Value"]   <- "pvalue"
  names(res_df)[names(res_df) == "adj.P.Val"] <- "padj"
  res_df <- res_df %>% arrange(padj)
  list(res_df = res_df, log_cpm = v$E, success = TRUE)
}, error = function(e) {
  message("[WARN] limma: ", conditionMessage(e))
  list(success = FALSE)
})


# ══════════════════════════════════════════════════════════════════════════════
# Apply cascade + Type I/II and write method-specific TSVs
# ══════════════════════════════════════════════════════════════════════════════
method_map <- list(
  edgeR_ciriquant = list(result = result_edger, out = de_edger_out),
  deseq2          = list(result = result_deseq, out = de_deseq_out),
  limma           = list(result = result_limma, out = de_limma_out)
)
final_results <- list()

for (mname in names(method_map)) {
  entry <- method_map[[mname]]
  r     <- entry$result
  if (!isTRUE(r$success)) {
    if (!is.null(entry$out)) {
      empty_df <- data.frame(circ_id = character(0), log2FC = numeric(0),
                             pvalue = numeric(0), padj = numeric(0))
      dir.create(dirname(entry$out), recursive = TRUE, showWarnings = FALSE)
      write.table(empty_df, entry$out, sep = "\t", quote = FALSE, row.names = FALSE)
      message(sprintf("[WARN] %s failed — wrote empty placeholder to %s", mname, entry$out))
    }
    next
  }

  casc     <- do_cascade(r$res_df, de_sig_by, fdr_cutoff)
  final_df <- add_type_col(casc$res_df, mname, casc$eff_col, casc$eff_thr,
                            fdr_cutoff, lfc_cutoff, fsj_lfc_thr)
  final_results[[mname]] <- list(
    res_df  = final_df,
    log_cpm = r$log_cpm,
    eff_col = casc$eff_col,
    eff_thr = casc$eff_thr
  )

  # Write per-method TSV only when the output path is declared
  if (!is.null(entry$out)) {
    dir.create(dirname(entry$out), recursive = TRUE, showWarnings = FALSE)
    write.table(final_df, entry$out, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[OK] %s → %s (%d circRNAs)", mname, entry$out, nrow(final_df)))
  }
}

# ── Primary de_results.tsv ────────────────────────────────────────────────────
primary_method <- de_method
if (is.null(final_results[[primary_method]])) {
  for (fb in c("deseq2", "limma", "edgeR_ciriquant")) {
    if (!is.null(final_results[[fb]])) { primary_method <- fb; break }
  }
}
if (is.null(final_results[[primary_method]]))
  stop("All three DE methods failed — cannot write de_results.tsv")

primary <- final_results[[primary_method]]
res_df  <- primary$res_df
log_cpm <- primary$log_cpm
eff_col <- primary$eff_col
eff_thr <- primary$eff_thr

write.table(res_df, de_out, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[OK] Primary (%s) → %s (%d circRNAs)", primary_method, de_out, nrow(res_df)))


# ══════════════════════════════════════════════════════════════════════════════
# Shared plot helpers  (primary method only)
# ══════════════════════════════════════════════════════════════════════════════
sig_col <- c(Up = "#d62728", Down = "#2CA02C", NS = "grey70")

plot_y_label <- switch(eff_col,
  "pvalue"  = expression(-log[10]~"(p-value, nominal)"),
  "qvalue"  = expression(-log[10]~"(Storey q-value)"),
  expression(-log[10]~"(adjusted p-value)")
)

plot_df <- res_df %>%
  filter(!is.na(.data[[eff_col]])) %>%
  mutate(sig = case_when(
    .data[[eff_col]] < eff_thr & log2FC >  lfc_cutoff ~ "Up",
    .data[[eff_col]] < eff_thr & log2FC < -lfc_cutoff ~ "Down",
    TRUE ~ "NS"
  ))

# Heatmap top IDs (used for volcano annotation + heatmap)
heat_up_ids <- res_df %>%
  filter(!is.na(.data[[eff_col]]), log2FC > 0) %>%
  arrange(.data[[eff_col]]) %>%
  head(heatmap_top_n) %>% pull(circ_id)
heat_dn_ids <- res_df %>%
  filter(!is.na(.data[[eff_col]]), log2FC < 0) %>%
  arrange(.data[[eff_col]]) %>%
  head(heatmap_top_n) %>% pull(circ_id)
heat_ids    <- c(heat_up_ids, heat_dn_ids)
heat_plot_df <- plot_df %>% filter(circ_id %in% heat_ids)


# ── Volcano ───────────────────────────────────────────────────────────────────
p_volc <- ggplot(plot_df, aes(x = log2FC, y = -log10(.data[[eff_col]]),
                               colour = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = sig_col) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed",
             colour = "grey40") +
  geom_hline(yintercept = -log10(eff_thr), linetype = "dashed",
             colour = "grey40") +
  labs(
    title  = paste0("Volcano (", tumor_label, " vs ", normal_label,
                    ")  [", primary_method, "]"),
    x      = expression(log[2]~"Fold Change"),
    y      = plot_y_label,
    colour = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

if (nrow(heat_plot_df) > 0) {
  p_volc <- p_volc +
    geom_point(data = heat_plot_df,
               aes(x = log2FC, y = -log10(.data[[eff_col]])),
               shape = 21, size = 4, fill = NA,
               colour = "black", stroke = 1.2, inherit.aes = FALSE) +
    ggrepel::geom_text_repel(
               data = heat_plot_df,
               aes(x = log2FC, y = -log10(.data[[eff_col]]), label = circ_id),
               size = 2.5, colour = "black", max.overlaps = 20, inherit.aes = FALSE)
}

pdf(volcano_out, width = 7, height = 6)
print(p_volc)
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
  scale_colour_manual(values = setNames(c("#d62728", "#2CA02C"), c(tumor_label, normal_label))) +
  labs(
    title  = paste0("PCA  [", primary_method, "]"),
    x      = paste0("PC1: ", var_pct[1], "%"),
    y      = paste0("PC2: ", var_pct[2], "%"),
    colour = "Condition"
  ) +
  theme_bw(base_size = 13)
dev.off()
message("[OK] PCA → ", pca_out)


# ── Heatmap ───────────────────────────────────────────────────────────────────
# Build row labels: circBase ID for known circRNAs, raw circ_id otherwise
cb_label_map <- setNames(character(0), character(0))
if (!is.null(circbase_annot_file) && file.exists(circbase_annot_file)) {
  cb_df <- tryCatch(
    read.table(circbase_annot_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cb_df) && all(c("circ_id", "circbase_id", "in_circbase") %in% colnames(cb_df))) {
    cb_known     <- cb_df[cb_df$in_circbase == 1, ]
    cb_label_map <- setNames(cb_known$circbase_id, cb_known$circ_id)
  }
}
make_heatmap_label <- function(cid) {
  if (cid %in% names(cb_label_map)) cb_label_map[[cid]] else cid
}

mat <- log_cpm[intersect(heat_ids, rownames(log_cpm)), , drop = FALSE]
rownames(mat) <- sapply(rownames(mat), make_heatmap_label)

# Normal-centered z-score
normal_cols <- common_samples[condition == normal_label]
normal_mat  <- mat[, normal_cols, drop = FALSE]
row_mean_n  <- rowMeans(normal_mat)
row_sd_all  <- apply(mat, 1, sd)
row_sd_all[row_sd_all < 0.1] <- 0.1
mat <- sweep(sweep(mat, 1, row_mean_n, "-"), 1, row_sd_all, "/")

ann_col <- data.frame(condition = condition, row.names = common_samples)
ann_col_colours <- list(
  condition = setNames(c("#d62728", "#2CA02C"), c(tumor_label, normal_label))
)

n_up_heat <- sum(heat_up_ids %in% rownames(log_cpm))
n_dn_heat <- sum(heat_dn_ids %in% rownames(log_cpm))
pdf(heatmap_out, width = 10, height = max(6, nrow(mat) * 0.4 + 3))
if (nrow(mat) >= 2) {
  pheatmap(
    mat,
    annotation_col    = ann_col,
    annotation_colors = ann_col_colours,
    color             = colorRampPalette(c("#2ca02c", "white", "#d62728"))(100),
    show_rownames     = TRUE,
    fontsize_row      = 8,
    border_color      = NA,
    main              = paste0("Top ", n_up_heat, " up + ", n_dn_heat,
                               " down DE circRNAs  [", primary_method, "]")
  )
} else {
  plot.new()
  text(0.5, 0.5, paste0("Too few DE circRNAs for heatmap (n=", nrow(mat), ")"), cex = 1.2)
}
dev.off()
message("[OK] Heatmap → ", heatmap_out)
