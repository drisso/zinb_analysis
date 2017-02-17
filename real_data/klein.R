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

esd0 <- read.csv("Klein/GSM1599494_ES_d0_main.csv.bz2", header=FALSE,
                 stringsAsFactors=FALSE, row.names=1)
esd2 <- read.csv("Klein/GSM1599497_ES_d2_LIFminus.csv.bz2", header=FALSE,
                 stringsAsFactors=FALSE, row.names=1)
esd4 <- read.csv("Klein/GSM1599498_ES_d4_LIFminus.csv.bz2", header=FALSE,
                 stringsAsFactors=FALSE, row.names=1)
counts <- cbind(esd0, esd2, esd4)

day <- c(rep("d0", ncol(esd0)), rep("d2", ncol(esd2)), rep("d4", ncol(esd4)))
day <- as.factor(day)
table(day)

filter <- rowSums(counts>10)>=10
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

plot(pc_cpm, pch=19, col=col1[day], main="PCA")
cor(pc_cpm[,1], detection_rate)
cor(pc_cpm[,2], detection_rate)

zinb <- zinbFit(as.matrix(raw), ncores = 7, K = 2, epsilon=1e3)
plot(zinb@W, pch=19, col=col1[day], main="ZINB")
cor(zinb@W[,1], detection_rate)
cor(zinb@W[,2], detection_rate)

wrapRzifa <- function(Y, block = T){
  # wrapper R function for ZIFA.
  # md5 hashing and temporary files are used not to re-run zifa
  # if it has already be run on this computer.
  d = digest(Y, "md5")
  tmp = paste0(tempdir(), '/', d)
  write.csv(Y, paste0(tmp, '.csv'))

  if (!file.exists(paste0(tmp, '_zifa.csv'))){
    print('run ZIFA')
    bb = ifelse(block, '-b ', '')
    cmd = sprintf('python run_zifa.py %s%s.csv %s_zifa.csv', bb, tmp, tmp)
    system(cmd)
  }
  read.csv(sprintf("%s_zifa.csv", tmp), header=FALSE)
}

zifa_raw <- wrapRzifa(log1p(tc))
plot(zifa_raw, pch=19, col=col1[day], main="ZIFA")
cor(zifa_raw[,1], detection_rate)
cor(zifa_raw[,2], detection_rate)

data.frame(Dim1=pc_cpm[,1], Dim2=pc_cpm[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=detection_rate)) +
  scale_colour_gradient(low="blue", high="green") + geom_point()

data.frame(Dim1=zifa_raw[,1], Dim2=zifa_raw[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=detection_rate)) +
  scale_colour_gradient(low="blue", high="green") + geom_point()

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=detection_rate)) +
  scale_colour_gradient(low="blue", high="green") + geom_point()

cors <- sapply(seq_len(nrow(tc)), function(i) {
  cor(pc_cpm[,1], tc[i,])
})
names(cors) <- rownames(tc)
cor_tc <- names(sort(cors, decreasing = TRUE)[1:100])

cors <- sapply(seq_len(nrow(tc)), function(i) {
  cor(zinb@W[,1], tc[i,])
})
names(cors) <- rownames(tc)
cor_zinb <- names(sort(cors, decreasing = TRUE)[1:100])

cors <- sapply(seq_len(nrow(tc)), function(i) {
  cor(zifa_raw[,1], tc[i,])
})
names(cors) <- rownames(tc)
cor_zifa <- names(sort(cors, decreasing = TRUE)[1:100])

library(RCurl)
url <- "http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"
f <- getBinaryURL(url)
load(rawConnection(f))

library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
bm <- getBM(attributes=c('entrezgene', 'mgi_symbol'), mart = mart)
bm <- bm[match(rownames(tc), bm[,2]),]
bm <- na.omit(bm)
msigdb_idx <- ids2indices(Mm.c2, bm[,1])
msigdb_idx <- msigdb_idx[sapply(msigdb_idx, length)>=10]
keep <- grepl("^BIOCARTA_", names(msigdb_idx)) | grepl("^KEGG_", names(msigdb_idx)) | grepl("^REACTOME_", names(msigdb_idx)) | grepl("^PID_", names(msigdb_idx)) | grepl("^SA_", names(msigdb_idx)) | grepl("^ST_", names(msigdb_idx))

msigdb_idx <- msigdb_idx[keep]

pvals <- sapply(seq_along(msigdb_idx), function(i) {
  geneset_genes <- bm[msigdb_idx[[i]],2]
  a <- sum(cor_tc %in% geneset_genes)
  b <- length(cor_tc) - a
  c <- sum(setdiff(bm[,1], cor_tc) %in% geneset_genes)
  d <- length(setdiff(bm[,1], cor_tc)) - c
  tab <- matrix(c(a, b, c, d), ncol=2)
  return(fisher.test(tab)$p.value)
})
names(pvals) <- names(msigdb_idx)

adjp <- p.adjust(pvals, method = "BH")
print(head(sort(adjp)))

pvals <- sapply(seq_along(msigdb_idx), function(i) {
  geneset_genes <- bm[msigdb_idx[[i]],2]
  a <- sum(cor_zifa %in% geneset_genes)
  b <- length(cor_zifa) - a
  c <- sum(setdiff(bm[,1], cor_zifa) %in% geneset_genes)
  d <- length(setdiff(bm[,1], cor_zifa)) - c
  tab <- matrix(c(a, b, c, d), ncol=2)
  return(fisher.test(tab)$p.value)
})
names(pvals) <- names(msigdb_idx)

adjp <- p.adjust(pvals, method = "BH")
print(head(sort(adjp)))

pvals <- sapply(seq_along(msigdb_idx), function(i) {
  geneset_genes <- bm[msigdb_idx[[i]],2]
  a <- sum(cor_zinb %in% geneset_genes)
  b <- length(cor_zinb) - a
  c <- sum(setdiff(bm[,1], cor_zinb) %in% geneset_genes)
  d <- length(setdiff(bm[,1], cor_zinb)) - c
  tab <- matrix(c(a, b, c, d), ncol=2)
  return(fisher.test(tab)$p.value)
})
names(pvals) <- names(msigdb_idx)

adjp <- p.adjust(pvals, method = "BH")
print(head(sort(adjp)))

