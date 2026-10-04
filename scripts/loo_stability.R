# Leave-one-sample-out (LOO) runs of the three DE chains, for effect-size stability.
#
# Mirrors the model setup of scripts/analysis.R exactly:
#   - design: ~ patient + condition if the groups file has a patient_id column
#             (the same auto rule as analysis.R), otherwise ~ condition
#   - edgeR_ciriquant: filterByExpr(group = condition), FSJ TMM offset, glmQLF
#   - DESeq2: poscounts size factors, results(contrast = tumor vs normal)
#   - limma-voom: TMM, voom, eBayes, topTable on the condition coefficient
# Run once for NONE (full data) and once per dropped sample. Errors are NOT caught:
# any failure stops the script with a non-zero exit status, so no partial output is
# written silently. Validate the NONE runs against the published de_results*.tsv
# before using any LOO output (see audit/fsj_three_case/validate_loo_vs_published.py).
#
# usage: Rscript loo_stability.R <results_dir> <groups_csv> <out_dir>
#                                [tumor_label=tumor] [normal_label=normal]
# writes <out_dir>/<drop>/{edger,deseq2,limma}.tsv  (circ_id, log2FC, pvalue)
# and   <out_dir>/run_info.tsv                       (drop, n_samples, design, rows per method)

suppressPackageStartupMessages({ library(edgeR); library(DESeq2); library(limma) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: Rscript loo_stability.R <results_dir> <groups_csv> <out_dir> [tumor] [normal]")
results_dir <- args[1]; groups_file <- args[2]; out_dir <- args[3]
tumor_label  <- if (length(args) >= 4) args[4] else "tumor"
normal_label <- if (length(args) >= 5) args[5] else "normal"

counts_all <- read.table(file.path(results_dir, "circRNA", "count_matrix.tsv"), sep = "\t",
                         header = TRUE, row.names = 1, check.names = FALSE)
fsj_all <- read.table(file.path(results_dir, "circRNA", "fsj_count_matrix.tsv"), sep = "\t",
                      header = TRUE, row.names = 1, check.names = FALSE)
groups_full <- read.csv(groups_file)
common_all <- intersect(colnames(counts_all), groups_full$srr_id)
if (length(common_all) < 2) stop("fewer than 2 samples shared by count matrix and groups file")
if (!all(common_all %in% colnames(fsj_all))) stop("FSJ matrix lacks samples present in the BSJ matrix")

has_patient_full <- "patient_id" %in% colnames(groups_full)
drops <- c("NONE", common_all)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
info <- list()

run_one <- function(drop) {
  common <- common_all
  if (drop != "NONE") common <- setdiff(common_all, drop)
  if (length(common) < 4) stop(sprintf("drop=%s leaves %d samples; need >= 4", drop, length(common)))
  groups <- groups_full[match(common, groups_full$srr_id), , drop = FALSE]
  condition <- factor(groups$condition, levels = c(normal_label, tumor_label))
  if (any(is.na(condition))) stop("condition labels other than tumor/normal present")
  if (any(table(condition) < 2)) stop(sprintf("drop=%s leaves a group with < 2 samples", drop))

  counts <- data.matrix(round(counts_all[, common, drop = FALSE])); storage.mode(counts) <- "integer"
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  fsj <- data.matrix(round(fsj_all[, common, drop = FALSE])); storage.mode(fsj) <- "integer"

  if (has_patient_full) {
    patient <- factor(groups$patient_id)
    design <- model.matrix(~ patient + condition); deseq_design <- ~ patient + condition
  } else {
    design <- model.matrix(~ condition); deseq_design <- ~ condition
  }
  cond_coef <- ncol(design)  # condition is always the last column, as in analysis.R
  design_name <- if (has_patient_full) "paired(~patient+condition)" else "unpaired(~condition)"

  od <- file.path(out_dir, drop); dir.create(od, showWarnings = FALSE)

  # ---- edgeR_ciriquant (analysis.R lines ~200-219) ----
  common_circs <- intersect(rownames(counts), rownames(fsj))
  bsj <- counts[common_circs, ]; fsj_c <- fsj[common_circs, ]
  keep <- filterByExpr(bsj, group = condition, min.count = 5, min.total.count = 15)
  bsj <- bsj[keep, ]; fsj_c <- fsj_c[keep, ]
  fsj_pseudo <- fsj_c + 0.5
  y_fsj <- calcNormFactors(DGEList(counts = fsj_c, group = condition), method = "TMM")
  fsj_scaled <- t(t(fsj_pseudo) * y_fsj$samples$norm.factors)
  fsj_cpm <- t(t(fsj_scaled) / colSums(fsj_scaled) * 1e6)
  offset_mat <- log(fsj_cpm + 0.5)
  dge <- DGEList(counts = bsj, group = condition); dge$offset <- offset_mat
  dge <- estimateDisp(dge, design, robust = TRUE)
  qlf <- glmQLFTest(glmQLFit(dge, design, robust = TRUE), coef = cond_coef)
  r_e <- as.data.frame(topTags(qlf, n = Inf, sort.by = "none"))
  e <- data.frame(circ_id = rownames(r_e), log2FC = r_e$logFC, pvalue = r_e$PValue)
  write.table(e, file.path(od, "edger.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- DESeq2 (analysis.R lines ~268-280; poscounts size factors) ----
  col_data <- data.frame(row.names = common, condition = condition)
  if (has_patient_full) col_data$patient <- patient
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = deseq_design)
  dds <- estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds)
  rd <- as.data.frame(results(dds, contrast = c("condition", tumor_label, normal_label)))
  d2 <- data.frame(circ_id = rownames(rd), log2FC = rd$log2FoldChange, pvalue = rd$pvalue)
  write.table(d2, file.path(od, "deseq2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # ---- limma-voom (analysis.R lines ~298-305) ----
  dgl <- calcNormFactors(DGEList(counts = counts, group = condition), method = "TMM")
  v <- voom(dgl, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  rl <- topTable(fit, coef = cond_coef, n = Inf, sort.by = "P")
  l <- data.frame(circ_id = rownames(rl), log2FC = rl$logFC, pvalue = rl$P.Value)
  write.table(l, file.path(od, "limma.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  data.frame(drop = drop, n_samples = length(common), design = design_name,
             edger_rows = nrow(e), deseq2_rows = nrow(d2), limma_rows = nrow(l))
}

for (dr in drops) info[[dr]] <- run_one(dr)
info <- do.call(rbind, info)
write.table(info, file.path(out_dir, "run_info.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[loo] %d runs written to %s (%s)\n", nrow(info), out_dir,
            if (has_patient_full) "paired" else "unpaired"))
