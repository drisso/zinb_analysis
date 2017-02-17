library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(rARPACK)

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

data <- read.csv("GSE59114_C57BL6_GEO_all.csv", skip=1, stringsAsFactors = FALSE)
genenames <- data[,1]
data <- data[,3:1430]

filter <- rowSums(data>10)>=10
raw <- data[filter,] ## data in log2(TPM+1)

detection_rate <- colSums(data>0)
coverage <- colSums(data)
qc <- cbind(detection_rate, coverage)
head(qc)

tmp <- strsplit(colnames(raw), "_")
age <- as.factor(sapply(tmp, function(x) x[1]))
stage <- as.factor(sapply(tmp, function(x) x[2]))
table(stage, age)
cross <- as.factor(paste0(stage, "_", age))
# ST: short-term, LT: long-term, MPP: multi-potent progenitors

fastpca <- function(expr, scale=FALSE) {
  svd_raw <- svds(scale(t(expr), center=TRUE, scale=scale), k=3, nu=3, nv=0)
  pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:3])
  return(pc_raw)
}

pc_tpm <- fastpca(raw)

plot(pc_tpm, pch=19, col=col1[cross], main="PCA")
pairs(pc_tpm, pch=19, col=col1[cross], main="PCA")
cor(pc_tpm[,1], detection_rate)
cor(pc_tpm[,2], detection_rate)

##zinb (workaround -- on rounded TPM)
rawcount <- round(2^raw - 1)

zinb <- zinbFit(rawcount, ncores = 7, K = 2, epsilon=1e3)
plot(zinb@W, pch=19, col=col1[cross], main="ZINB")
plot(zinb@W, pch=19, col=col1[age], main="ZINB")
plot(zinb@W, pch=19, col=col1[stage], main="ZINB")

cor(zinb@W[,1], detection_rate)
cor(zinb@W[,2], detection_rate)

library(cluster)
pamres <- pam(zinb@W, k=4)
cl <- pamres$clustering
cl[cl==2] <- 1
cl[cl==3] <- 1
cl[cl==4] <- 2

plot(zinb@W, pch=19, col=col1[cl], main="ZINB")
plot(pc_tpm, pch=19, col=col1[cl], main="PCA")

library(stringr)
genenames<- gsub("'", '', genenames)
geneupper <- str_to_upper(genenames)

library(limma)
design <- model.matrix(~as.factor(cl))
fit <- lmFit(raw, design)
fit <- eBayes(fit)
interesting_genes <- genenames[as.numeric(rownames(topTable(fit, coef=2, n=100)))]

library(RCurl)
url <- "http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"
f <- getBinaryURL(url)
load(rawConnection(f))

library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
bm <- getBM(attributes=c('entrezgene', 'mgi_symbol'), mart = mart)
bm <- bm[match(genenames, bm[,2]),]
bm <- na.omit(bm)
msigdb_idx <- ids2indices(Mm.c2, bm[,1])
msigdb_idx <- msigdb_idx[sapply(msigdb_idx, length)>=10]
keep <- grepl("^BIOCARTA_", names(msigdb_idx)) | grepl("^KEGG_", names(msigdb_idx)) | grepl("^REACTOME_", names(msigdb_idx)) | grepl("^PID_", names(msigdb_idx)) | grepl("^SA_", names(msigdb_idx)) | grepl("^ST_", names(msigdb_idx))

msigdb_idx <- msigdb_idx[keep]

pvals <- sapply(seq_along(msigdb_idx), function(i) {
  geneset_genes <- bm[msigdb_idx[[i]],1]
  a <- sum(interesting_genes %in% geneset_genes)
  b <- length(interesting_genes) - a
  c <- sum(setdiff(bm[,1], interesting_genes) %in% geneset_genes)
  d <- length(setdiff(bm[,1], interesting_genes)) - c
  tab <- matrix(c(a, b, c, d), ncol=2)
  return(fisher.test(tab)$p.value)
})
names(pvals) <- names(msigdb_idx)

adjp <- p.adjust(pvals, method = "BH")
print(head(sort(pvals)))

