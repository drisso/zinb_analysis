library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(rARPACK)
library(edgeR)
library(digest)

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

data <- read.table("retina/GSM1626793_P14Retina_1.digital_expression.txt.gz",
                   header=TRUE, stringsAsFactors=FALSE, row.names=1)
cluster_lab <- read.table("retina/retina_clusteridentities.txt",
                          stringsAsFactors = FALSE)
cl <- cluster_lab[grep("^r1_", cluster_lab[,1]),2]
names(cl) <- str_sub(cluster_lab[grep("^r1_", cluster_lab[,1]),1], 4)
cl <- as.factor(cl)
table(cl)
stopifnot(all(names(cl) %in% colnames(data)))

counts <- data[, names(cl)]

filter <- rowSums(counts>5)>=5
table(filter)
raw <- counts[filter,]

detection_rate <- colSums(counts>0)
coverage <- colSums(counts)
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

pc_cpm <- fastpca(log1p(tc))

plot(pc_cpm, pch=19, col=col2[cl], main="PCA")
cor(pc_cpm[,1], detection_rate)
cor(pc_cpm[,2], detection_rate)

zinb <- zinbFit(as.matrix(raw), ncores = 7, K = 2, epsilon=1e3)
plot(zinb@W, pch=19, col=col1[day], main="ZINB")
cor(zinb@W[,1], detection_rate)
cor(zinb@W[,2], detection_rate)
