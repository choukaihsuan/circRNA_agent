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

fitfun <- function(drop=character(0)){
  keepS <- !(cs %in% drop)
  bs <- bsj[,keepS,drop=FALSE]; fs <- fsj[,keepS,drop=FALSE]; cd <- droplevels(condition[keepS])
  des <- model.matrix(~cd)
  fp <- fs+0.5; y <- calcNormFactors(DGEList(counts=fs,group=cd),method="TMM")
  sc <- t(t(fp)*y$samples$norm.factors); cp <- t(t(sc)/colSums(sc)*1e6); off <- log(cp+0.5)
  dg <- DGEList(counts=bs,group=cd); dg$offset <- off
  dg <- estimateDisp(dg,des,robust=TRUE); ft <- glmQLFit(dg,des,robust=TRUE)
  as.data.frame(topTags(glmQLFTest(ft,coef=ncol(des)),n=Inf,sort.by="none"))
}
full <- fitfun()
chk <- list(c("chr13:33091994|33101669","SRR11600338"), c("chr8:141856359|141900868","SRR11600331,SRR11600337"), c("chr10:7839010|7844817","SRR11600335"))
for (x in chk) {
  dr <- strsplit(x[2],",")[[1]]; r <- fitfun(dr)
  cat(sprintf("%s | all 12: log2FC=%.3f p=%.4g | drop %s: log2FC=%.3f p=%.4g\n", x[1], full[x[1],"logFC"], full[x[1],"PValue"], x[2], r[x[1],"logFC"], r[x[1],"PValue"]))
}
# BSJ-only crude (CPM of BSJ library over tested set)
lib <- colSums(bsj)
for (x in chk) { id <- x[1]; cpm_b <- bsj[id,]/lib*1e6
  cat(sprintf("%s | BSJ-only crude log2(mean tumor CPM / mean normal CPM)=%.3f\n", id, log2(mean(cpm_b[condition=="tumor"])/mean(cpm_b[condition=="normal"])))) }
