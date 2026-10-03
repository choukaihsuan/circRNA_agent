suppressPackageStartupMessages(library(edgeR))
d <- "GSE133998"; H <- "/home3/choukaihsuan/"
b <- read.table(paste0(H,d,"_results/circRNA/count_matrix.tsv"),sep="\t",header=TRUE,row.names=1,check.names=FALSE)
f <- read.table(paste0(H,d,"_results/circRNA/fsj_count_matrix.tsv"),sep="\t",header=TRUE,row.names=1,check.names=FALSE)
g <- read.csv(paste0("/home3/choukaihsuan/circRNA_agent/metadata/",d,"/sample_groups.csv"))
cs <- intersect(colnames(b),g$srr_id); g <- g[match(cs,g$srr_id),]
condition <- factor(g$condition,levels=c("normal","tumor"))
counts <- data.matrix(round(b[,cs])); storage.mode(counts) <- "integer"; counts <- counts[rowSums(counts)>0,,drop=FALSE]
fsj <- data.matrix(round(f[,cs])); storage.mode(fsj) <- "integer"
cc <- intersect(rownames(counts),rownames(fsj)); bsj <- counts[cc,]; fsj <- fsj[cc,]
keep <- filterByExpr(bsj,group=condition,min.count=5,min.total.count=15); bsj <- bsj[keep,]; fsj <- fsj[keep,]
design <- model.matrix(~condition); cond_coef <- ncol(design)
# ---- identical to analysis.R:207-212
fsj_pseudo <- fsj + 0.5
y_fsj_norm <- DGEList(counts=fsj, group=condition)
y_fsj_norm <- calcNormFactors(y_fsj_norm, method="TMM")
fsj_scaled <- t(t(fsj_pseudo) * y_fsj_norm$samples$norm.factors)
fsj_cpm <- t(t(fsj_scaled)/colSums(fsj_scaled)*1e6)
offset_mat <- log(fsj_cpm + 0.5)
dge <- DGEList(counts=bsj, group=condition); dge$offset <- offset_mat
dge <- estimateDisp(dge, design, robust=TRUE); fit <- glmQLFit(dge, design, robust=TRUE)
qlf <- glmQLFTest(fit, coef=cond_coef); res <- as.data.frame(topTags(qlf,n=Inf,sort.by="none"))
# compare with shipped de_results.tsv
de <- read.delim(paste0(H,d,"_results/de/de_results.tsv"), check.names=FALSE); rownames(de) <- de$circ_id
common <- intersect(rownames(res), rownames(de))
cat(sprintf("tested here=%d, in de_results.tsv=%d, common=%d\n",nrow(res),nrow(de),length(common)))
cat(sprintf("max |log2FC diff|=%.3g ; max |log10(p) diff|=%.3g\n", max(abs(res[common,"logFC"]-de[common,"log2FC"])), max(abs(log10(res[common,"PValue"])-log10(de[common,"pvalue"])))))
cat("TMM norm.factors (FSJ):\n"); print(setNames(round(y_fsj_norm$samples$norm.factors,4),cs))
cat("FSJ library size over tested circs (raw):\n"); print(setNames(colSums(fsj),cs))
cat("scaled lib (sum fsj_scaled):\n"); print(setNames(round(colSums(fsj_scaled),1),cs))
out <- NULL
for (id in c("chr13:33091994|33101669","chr8:141856359|141900868","chr10:7839010|7844817")) {
  if (!id %in% rownames(bsj)) next
  o <- data.frame(circ_id=id, sample=cs, condition=as.character(condition), BSJ=as.integer(bsj[id,]), FSJ=as.integer(fsj[id,]),
     TMM_norm_factor=round(y_fsj_norm$samples$norm.factors,4), FSJ_pseudo=fsj_pseudo[id,], FSJ_scaled=round(fsj_scaled[id,],4),
     FSJ_CPM_pseudo_TMM=round(fsj_cpm[id,],3), offset_ln_CPMplus0.5=round(offset_mat[id,],4), BSJ_over_expFSJoffset=round(bsj[id,]/exp(offset_mat[id,]),5), stringsAsFactors=FALSE)
  out <- rbind(out,o)
}
out$file_log2FC <- de[out$circ_id,"log2FC"]; out$file_pvalue <- de[out$circ_id,"pvalue"]; out$file_logFC_fsj <- de[out$circ_id,"logFC_fsj"]; out$file_pvalue_fsj <- de[out$circ_id,"pvalue_fsj"]; out$file_Type <- de[out$circ_id,"Type"]
out$refit_log2FC <- res[out$circ_id,"logFC"]; out$refit_pvalue <- res[out$circ_id,"PValue"]
write.table(out,"/home/choukaihsuan/audit_tmp/fsj_worked_example.tsv",sep="\t",quote=FALSE,row.names=FALSE)
print(out[out$circ_id=="chr13:33091994|33101669", c("sample","condition","BSJ","FSJ","FSJ_CPM_pseudo_TMM","offset_ln_CPMplus0.5")])
print(unique(out[,c("circ_id","file_log2FC","file_pvalue","file_logFC_fsj","file_pvalue_fsj","file_Type","refit_log2FC","refit_pvalue")]))
cat("getOffset uses y$offset when present:\n"); print(body(edgeR::getOffset))
