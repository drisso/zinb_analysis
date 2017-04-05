## The following code is from https://github.com/MarioniLab/PlateEffects2016/blob/master/reference/ESpresso.R

all.counts <- read.table("ESpresso/counttable_es.csv", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 704)))
serum <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(all.counts))
batch <- sub("ola_mES_[^_]+_([^_]+)_.*", "\\1", colnames(all.counts))
targets <- data.frame(Serum=serum, Batch=batch)

# Only using data from two batches.
keep <- targets$Batch %in% c("2", "3")
all.counts <- all.counts[,keep]
targets <- targets[keep,]
targets$Plate <- as.integer(factor(paste0(targets$Serum, targets$Batch)))
targets[] <- lapply(targets, factor)
targets$Serum <- factor(targets$Serum, c("lif", "2i", "a2i"))

# Removing spike-ins.
is.mouse <- grepl("^ENSMUSG", rownames(all.counts))
all.counts <- all.counts[is.mouse,]

## end of copied code

library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(rARPACK)
library(edgeR)
library(digest)
library(scone)
library(matrixStats)

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

filter <- rowSums(all.counts>10)>=10
table(filter)
raw <- all.counts[filter,]

detection_rate <- colSums(all.counts>0)
coverage <- colSums(all.counts)
qc <- cbind(detection_rate, coverage)
head(qc)

fastpca <- function(expr, scale=FALSE) {
  svd_raw <- svds(scale(t(expr), center=TRUE, scale=scale), k=3, nu=3, nv=0)
  pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:3])
  return(pc_raw)
}

totalcount = function (ei)
{
  sums = colSums(ei)
  eo = t(t(ei)*mean(sums)/sums)
  return(eo)
}

tc <- totalcount(raw)
fq <- as.matrix(FQT_FN(raw))

vars <- rowVars(log1p(fq))
names(vars) <- rownames(fq)
vars <- sort(vars, decreasing = TRUE)
vargenes <- names(vars)[1:1000]
pc_cpm <- fastpca(log1p(tc[vargenes,]))

plot(pc_cpm, pch=19, col=col1[targets$Serum], main="PCA, by Serum")
plot(pc_cpm, pch=19, col=col2[targets$Batch], main="PCA, by Batch")
plot(pc_cpm, pch=19, col=col2[targets$Plate], main="PCA, by Plate")

cor(pc_cpm[,1], detection_rate)
cor(pc_cpm[,2], detection_rate)

library(zinb)
zinb <- zinbFit(as.matrix(raw[vargenes,]), ncores = 7, K = 2, epsilon=1e3)
plot(zinb@W, pch=19, col=col1[targets$Serum], main="ZINB, by Serum")
plot(zinb@W, pch=19, col=col2[targets$Batch], main="ZINB, by Batch")
plot(zinb@W, pch=19, col=col2[targets$Plate], main="ZINB, by Plate")

cor(zinb@W[,1], detection_rate)
cor(zinb@W[,2], detection_rate)

zinb_batch <- zinbFit(as.matrix(raw[vargenes,]), ncores = 7, K = 2, epsilon=1e3,
                      X=model.matrix(~Batch, data=targets))
plot(zinb_batch@W, pch=19, col=col1[targets$Serum], main="ZINB, by Serum")
plot(zinb_batch@W, pch=19, col=col2[targets$Batch], main="ZINB, by Batch")
plot(zinb_batch@W, pch=19, col=col2[targets$Plate], main="ZINB, by Plate")

cor(zinb_batch@W[,1], detection_rate)
cor(zinb_batch@W[,2], detection_rate)

zinb_bio <- zinbFit(as.matrix(raw[vargenes,]), ncores = 7, K = 2, epsilon=1e3,
                      X=model.matrix(~Serum, data=targets))
plot(zinb_bio@W, pch=19, col=col1[targets$Serum], main="ZINB, by Serum")
plot(zinb_bio@W, pch=19, col=col2[targets$Batch], main="ZINB, by Batch")
plot(zinb_bio@W, pch=19, col=col2[targets$Plate], main="ZINB, by Plate")

cor(zinb_bio@W[,1], detection_rate)
cor(zinb_bio@W[,2], detection_rate)
