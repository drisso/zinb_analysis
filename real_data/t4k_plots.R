library(cowplot)
library(RColorBrewer)
library(magrittr)
library(zinbwave)
library(scater)
library(cellrangerRkit)

load("tenx_t4k/zinb_res.rda")
W <- getW(zinb_res)

pal <- clusterExperiment::bigPalette
pal2 <- c(brewer.pal(8, "Reds")[c(3, 6)], brewer.pal(8, "Purples")[c(3, 4, 6, 8)],
          brewer.pal(8, "Blues")[c(3, 5, 7)])

## data
pipestance_path <- "tenx_t4k"
pbmc <- load_cellranger_matrix(pipestance_path)
use_genes <- get_nonzero_genes(pbmc)

load("~/Google Drive/SCONE_COPY/t_4k/scone_update.rda")
keep <- which(colnames(pbmc) %in% colnames(scone_obj))

dense <- as.matrix(exprs(pbmc[use_genes, keep]))

vars <- rowVars(log1p(dense))
names(vars) <- rownames(dense)
vars <- sort(vars, decreasing = TRUE)
vargenes <- names(vars)[1:1000]

dense <- dense[vargenes,]

clusters <- read.csv("tenx_t4k/clusters.csv")
cl <- clusters[,2]
names(cl) <- clusters[,1]
cl <- cl[keep]

sceset <- newSCESet(countData = exprs(pbmc[use_genes, keep]))

sceset <- calculateQCMetrics(sceset)
qc <- pData(sceset)[,c(2, 5, 7, 12:15)]

ribo_idx <- grep("^RPL", fData(pbmc)[,2])
mito_idx <- grep("^MT-", fData(pbmc)[,2])
ribo_pct <- colSums(as.matrix(exprs(pbmc[ribo_idx,keep])))/colSums(as.matrix(exprs(pbmc[,keep]))) * 100
mito_pct <- colSums(as.matrix(exprs(pbmc[mito_idx,keep])))/colSums(as.matrix(exprs(pbmc[,keep]))) * 100

qc <- cbind(qc, pct_ribo = ribo_pct, pct_mito = mito_pct)

## PCA
library(rARPACK)
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

tc <- totalcount(as.matrix(exprs(pbmc[use_genes, keep])))
pc_tc <- fastpca(log1p(tc[rownames(dense),]))

## ZIFA (TODO!)

## Figure
data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2], OriginalCluster=as.factor(cl)) %>%
  ggplot(aes(Dim1, Dim2, color=OriginalCluster)) + geom_point() +
  scale_color_manual(values = pal) -> panel1_pca

# data.frame(Dim1=zifa_fq[,1], Dim2=zifa_fq[,2]) %>%
#   ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
#   scale_color_brewer(palette="Set1") -> panel1_zifa

data.frame(Dim1=W[,1], Dim2=W[,2], OriginalCluster=as.factor(cl)) %>%
  ggplot(aes(Dim1, Dim2, colour=OriginalCluster)) + geom_point()  +
  scale_color_manual(values = pal)  -> panel1_zinb

# p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
#           panel1_zifa + theme(legend.position = "none"),
#           panel1_zinb + theme(legend.position = "none"),
#           labels=c("a", "c", "e"), align = "h", ncol=3)
#
# legend <- get_legend(panel1_pca)
# upper <- plot_grid(p1, legend, rel_widths = c(3, .6))

cors <- lapply(1:3, function(i) abs(cor(pc_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(colnames(qc), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=pal2) + ylim(0, 1) -> panel2_pca

## repeat for ZIFA

cors <- lapply(1:3, function(i) abs(cor(W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(colnames(qc), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=pal2) + ylim(0, 1) -> panel2_zinb

# p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
#                 panel2_zifa + theme(legend.position = "none"),
#                 panel2_zinb + theme(legend.position = "none"),
#                 labels=c("b", "d", "f"), align = "h", ncol=3)
#
# legend2 <- get_legend(panel2_pca)
# lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

# fig1bis <- plot_grid(upper, lower, ncol=1, nrow=2)
#
# save_plot("patel_fig1bis.pdf", fig1bis,
#           ncol = 3,
#           nrow = 3,
#           base_aspect_ratio = 1.3
# )

ggsave("10x_t4k_cor2.png", panel2_pca)

if(run_zinb) {

  library(BiocParallel)
  library(doParallel)
  registerDoParallel(6)
  register(DoparParam())

  Xmod <- model.matrix(~pct_ribo + pct_mito, data = qc)
  system.time(zinb_res <- zinbFit(dense, K=3, X=Xmod))
  save(zinb_res, file="tenx_t4k/zinb_ribomito.rda")
} else {
  load("tenx_t4k/zinb_ribomito.rda")
}
W <- getW(zinb_res)

data.frame(Dim1=W[,1], Dim2=W[,2], OriginalCluster=as.factor(cl)) %>%
  ggplot(aes(Dim1, Dim2, colour=OriginalCluster)) + geom_point()  +
  scale_color_manual(values = pal)

cors <- lapply(1:3, function(i) abs(cor(W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(colnames(qc), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=pal2) + ylim(0, 1) -> panel3_zinb

ggsave("10x_t4k_cor3.png", panel2_zinb)
ggsave("10x_t4k_cor4.png", panel3_zinb)
