#!/usr/bin/env Rscript
#
# run_de_timed.R -- benchmark-only wrapper to time a DE analysis.R run.
#
# Builds the same Snakemake S4 mock object used by the documented manual
# rerun pattern (CLAUDE.md "Server 端獨立重跑腳本"), then sources the
# *unmodified* scripts/analysis.R. This exists purely so the benchmark's
# DE step can be wrapped in `/usr/bin/time -v` at the shell level without
# touching analysis.R itself (which is live-read by the production
# `de_analysis` Snakemake rule for datasets that may be running concurrently).
#
# All arguments are key=value pairs (order-independent).
#
# Usage:
#   Rscript run_de_timed.R \
#     analysis_r=/path/to/scripts/analysis.R \
#     matrix=... fsj_matrix=... groups=... circbase_annot=... \
#     de_method=edgeR_ciriquant fdr=0.05 lfc=1.0 \
#     tumor_label=tumor normal_label=normal de_sig_by=auto \
#     heatmap_top_n=10 fsj_concordance_lfc=0.0 \
#     out_de=... out_de_edger=... out_de_deseq=... out_de_limma=... \
#     out_volcano=... out_heatmap=... out_pca=... \
#     log=...

args_raw <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args_raw, "=", fixed = TRUE)
args <- setNames(vapply(kv, `[`, character(1), 2), vapply(kv, `[`, character(1), 1))

req <- function(name, default = NULL) {
  if (!is.na(args[name]) && nzchar(args[name])) return(args[[name]])
  if (!is.null(default)) return(default)
  stop(sprintf("missing required argument: %s", name))
}

setClass("Snakemake", representation(input = "list", output = "list",
                                      params = "list", log = "list"))

# fsj_matrix is optional: omitting it (or passing fsj_matrix=) makes
# analysis.R fall back to deseq2/limma-only, matching how the existing
# de_deseq2_baseline benchmark rule simulates a tool with no BSJ/FSJ
# ratio capability (see CLAUDE.md's documented analysis.R fallback).
fsj_raw <- args["fsj_matrix"]
fsj_val <- if (!is.na(fsj_raw) && nzchar(fsj_raw)) fsj_raw else NULL

snakemake <- new(
  "Snakemake",
  input = list(
    matrix         = req("matrix"),
    fsj_matrix     = fsj_val,
    groups         = req("groups"),
    circbase_annot = req("circbase_annot")
  ),
  output = list(
    de       = req("out_de"),
    de_edger = req("out_de_edger"),
    de_deseq = req("out_de_deseq"),
    de_limma = req("out_de_limma"),
    volcano  = req("out_volcano"),
    heatmap  = req("out_heatmap"),
    pca      = req("out_pca")
  ),
  params = list(
    de_method           = req("de_method"),
    fdr                 = as.numeric(req("fdr", "0.05")),
    lfc                 = as.numeric(req("lfc", "1.0")),
    tumor_label         = req("tumor_label", "tumor"),
    normal_label        = req("normal_label", "normal"),
    de_sig_by           = req("de_sig_by", "auto"),
    heatmap_top_n       = as.numeric(req("heatmap_top_n", "10")),
    fsj_concordance_lfc = as.numeric(req("fsj_concordance_lfc", "0.0"))
  ),
  log = list(req("log", "/dev/stderr"))
)

analysis_r <- req("analysis_r")
source(analysis_r)
