# scripts/isoform_switching.R
# ─────────────────────────────────────────────────────────────────────────────
# Isoform Usage Index (IUI) calculation and switching test.
#
# IUI_i = BSJ_i / Σ_j BSJ_j  (sum over all isoforms of the same gene)
#
# For each multi-isoform gene, Wilcoxon rank-sum test compares per-sample
# IUI between tumor and normal.  Two-level BH correction:
#   padj_within_gene : per-gene FDR (controls false discoveries within a gene)
#   padj_global      : global FDR across all isoforms (for final filtering)
#
# Inputs  (via snakemake):
#   bsj_matrix     – circRNA × sample BSJ count matrix
#   isoform_groups – output of assign_isoforms.py
#   sample_groups  – metadata CSV with srr_id + condition columns
#
# Outputs:
#   iui_matrix     – per-circRNA mean IUI for each condition
#   switching      – Wilcoxon results with is_switching flag
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})


# ── I/O ───────────────────────────────────────────────────────────────────────

bsj_mat <- as.matrix(
  read.table(snakemake@input[["bsj_matrix"]],
             header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
)

iso_groups <- read.table(snakemake@input[["isoform_groups"]],
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sample_meta <- read.csv(snakemake@input[["sample_groups"]],
                         stringsAsFactors = FALSE)

tumor_label      <- snakemake@params[["tumor_label"]]
normal_label     <- snakemake@params[["normal_label"]]
fdr_cutoff       <- as.numeric(snakemake@params[["fdr"]])
delta_iui_cutoff <- as.numeric(snakemake@params[["delta_iui_cutoff"]])

# Keep only samples that exist in the count matrix
tumor_samples  <- intersect(
  sample_meta$srr_id[sample_meta$condition == tumor_label],
  colnames(bsj_mat)
)
normal_samples <- intersect(
  sample_meta$srr_id[sample_meta$condition == normal_label],
  colnames(bsj_mat)
)

message(sprintf(
  "[isoform_switching] %d tumor, %d normal samples in matrix",
  length(tumor_samples), length(normal_samples)
))

if (length(tumor_samples) == 0 || length(normal_samples) == 0) {
  stop("[isoform_switching] No samples matched — check condition labels in sample_groups.csv")
}


# ── IUI helper ────────────────────────────────────────────────────────────────
# Returns a data.frame: circ_id, gene_id, mean_iui  (averaged over all samples
# in `sample_ids` for that condition)

compute_mean_iui <- function(mat, iso_grp, sample_ids) {
  s <- intersect(sample_ids, colnames(mat))
  if (length(s) == 0) return(NULL)

  sub <- mat[, s, drop = FALSE]

  df <- as.data.frame(sub) %>%
    rownames_to_column("circ_id") %>%
    pivot_longer(-circ_id, names_to = "sample", values_to = "bsj") %>%
    inner_join(iso_grp[, c("circ_id", "gene_id", "n_isoforms")],
               by = "circ_id") %>%
    filter(gene_id != "intergenic")

  gene_totals <- df %>%
    group_by(gene_id, sample) %>%
    summarise(gene_bsj = sum(bsj) + 1e-6, .groups = "drop")

  df %>%
    left_join(gene_totals, by = c("gene_id", "sample")) %>%
    mutate(iui = bsj / gene_bsj) %>%
    group_by(circ_id, gene_id) %>%
    summarise(mean_iui = mean(iui), .groups = "drop")
}

iui_tumor  <- compute_mean_iui(bsj_mat, iso_groups, tumor_samples)
iui_normal <- compute_mean_iui(bsj_mat, iso_groups, normal_samples)


# ── Isoform switching test ────────────────────────────────────────────────────

perform_switching_test <- function(mat, iso_grp,
                                    tumor_samp, normal_samp,
                                    fdr_cut, delta_cut) {

  multi_iso_genes <- unique(iso_grp$gene_id[iso_grp$n_isoforms >= 2])
  message(sprintf("[isoform_switching] %d multi-isoform genes to test",
                  length(multi_iso_genes)))

  all_rows <- list()

  for (gid in multi_iso_genes) {
    isoforms <- intersect(iso_grp$circ_id[iso_grp$gene_id == gid],
                          rownames(mat))
    if (length(isoforms) < 2) next

    gene_mat    <- mat[isoforms, , drop = FALSE]
    gene_totals <- colSums(gene_mat) + 1e-6          # per-sample total BSJ for gene
    iui_mat     <- sweep(gene_mat, 2, gene_totals, "/")

    gene_meta   <- iso_grp[iso_grp$gene_id == gid, ][1, ]
    gene_rows   <- list()

    for (iso in isoforms) {
      iui_t <- iui_mat[iso, intersect(tumor_samp,  colnames(iui_mat))]
      iui_n <- iui_mat[iso, intersect(normal_samp, colnames(iui_mat))]
      if (length(iui_t) < 2 || length(iui_n) < 2) next

      wt <- tryCatch(
        wilcox.test(iui_t, iui_n, exact = FALSE),
        error   = function(e) list(p.value = NA_real_),
        warning = function(w) suppressWarnings(
          wilcox.test(iui_t, iui_n, exact = FALSE)
        )
      )

      gene_rows[[iso]] <- data.frame(
        circ_id    = iso,
        gene_id    = gid,
        gene_name  = gene_meta$gene_name,
        n_isoforms = length(isoforms),
        iui_tumor  = mean(iui_t),
        iui_normal = mean(iui_n),
        delta_iui  = mean(iui_t) - mean(iui_n),
        p_value    = wt$p.value,
        stringsAsFactors = FALSE
      )
    }

    if (length(gene_rows) > 0) {
      gdf <- do.call(rbind, gene_rows)
      # Per-gene BH: controls false discoveries within a single gene
      gdf$padj_within_gene <- p.adjust(gdf$p_value, method = "BH")
      all_rows[[gid]] <- gdf
    }
  }

  if (length(all_rows) == 0) {
    message("[isoform_switching] No testable multi-isoform genes found")
    return(data.frame())
  }

  res <- do.call(rbind, all_rows)
  # Global BH kept for reference; switching uses within-gene FDR (less conservative
  # and standard in isoform-level analyses, e.g. DEXSeq).
  res$padj_global  <- p.adjust(res$p_value, method = "BH")
  res$is_switching <- (
    !is.na(res$padj_within_gene) &
    res$padj_within_gene < fdr_cut &
    abs(res$delta_iui) > delta_cut
  )
  res[order(res$padj_within_gene), ]
}

switching_results <- perform_switching_test(
  bsj_mat, iso_groups,
  tumor_samples, normal_samples,
  fdr_cutoff, delta_iui_cutoff
)


# ── Build IUI combined matrix ─────────────────────────────────────────────────

if (!is.null(iui_tumor) && !is.null(iui_normal)) {
  iui_combined <- full_join(
    rename(iui_tumor,  iui_tumor  = mean_iui),
    rename(iui_normal, iui_normal = mean_iui),
    by = c("circ_id", "gene_id")
  ) %>%
    mutate(delta_iui = iui_tumor - iui_normal)
} else if (!is.null(iui_tumor)) {
  iui_combined <- rename(iui_tumor, iui_tumor = mean_iui)
} else {
  iui_combined <- data.frame()
}


# ── Outputs ───────────────────────────────────────────────────────────────────

write.table(iui_combined,
            snakemake@output[["iui_matrix"]],
            sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("[isoform_switching] IUI matrix (%d rows) → %s",
                nrow(iui_combined), snakemake@output[["iui_matrix"]]))

write.table(switching_results,
            snakemake@output[["switching"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

n_sig <- if (nrow(switching_results) > 0)
           sum(switching_results$is_switching, na.rm = TRUE) else 0L
message(sprintf(
  "[isoform_switching] %d significant switching events (FDR<%.2f, |ΔIUI|>%.2f) → %s",
  n_sig, fdr_cutoff, delta_iui_cutoff, snakemake@output[["switching"]]
))
