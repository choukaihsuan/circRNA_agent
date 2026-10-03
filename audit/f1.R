suppressPackageStartupMessages(library(edgeR))
H <- "/home3/choukaihsuan/"; M <- "/home3/choukaihsuan/circRNA_agent/metadata/"
load_ds <- function(d){
  b <- read.table(paste0(H,d,"_results/circRNA/count_matrix.tsv"),sep="\t",header=TRUE,row.names=1,check.names=FALSE)
  f <- read.table(paste0(H,d,"_results/circRNA/fsj_count_matrix.tsv"),sep="\t",header=TRUE,row.names=1,check.names=FALSE)
  g <- read.csv(paste0(M,d,"/sample_groups.csv"))
  cs <- intersect(colnames(b),g$srr_id); g <- g[match(cs,g$srr_id),]
  cond <- factor(g$condition,levels=c("normal","tumor"))
  b <- data.matrix(round(b[,cs,drop=FALSE])); f <- data.matrix(round(f[,cs,drop=FALSE]))
  b <- b[rowSums(b)>0,,drop=FALSE]
  cc <- intersect(rownames(b),rownames(f)); list(b=b[cc,],f=f[cc,],cond=cond,samples=cs)
}
cat("== filterByExpr applied to BSJ vs to FSJ (same params as analysis.R:200-202)\n")
for (d in c("GSE108735","GSE133998","GSE113230","GSE130078","GSE148036")) {
  x <- load_ds(d)
  kb <- filterByExpr(x$b, group=x$cond, min.count=5, min.total.count=15)
  kf <- filterByExpr(x$f, group=x$cond, min.count=5, min.total.count=15)
  cat(sprintf("%s circs=%d pass_BSJ=%d (%.2f%%) pass_if_FSJ=%d (%.2f%%) both=%d | median per-sample BSJ sum=%d FSJ sum=%d | FSJ==0 cells=%.1f%% BSJ==0 cells=%.1f%%\n",
     d,nrow(x$b),sum(kb),100*mean(kb),sum(kf),100*mean(kf),sum(kb&kf),as.integer(median(colSums(x$b))),as.integer(median(colSums(x$f))),100*mean(x$f==0),100*mean(x$b==0)))
}
