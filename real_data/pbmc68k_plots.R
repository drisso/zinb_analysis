rm(list=ls())
library(cellrangerRkit)
library(zinbwave)
library(matrixStats)

analysis <- FALSE

if(analysis) {
  set.seed(1241241)
  pipestance_path <- "tenx_pbmc"
  if(!file.exists(paste0(pipestance_path, "/outs"))) {
    download_sample(sample_name="fresh_68k_pbmc_donor_a",sample_dir=pipestance_path,
                    host="http://cf.10xgenomics.com/samples/cell-exp/1.1.0/")
  }
  pbmc <- load_cellranger_matrix(pipestance_path)
  dim(pbmc)

  ## filtering
  use_genes <- get_nonzero_genes(pbmc)
  dense <- as.matrix(exprs(pbmc[use_genes,]))
  vars <- rowVars(log1p(dense))
  names(vars) <- rownames(dense)
  vars <- sort(vars, decreasing = TRUE)
  vargenes <- names(vars)[1:1000]

  dense <- dense[vargenes,]
  dim(dense)

  ## normalization and PCA
  library(rARPACK)
  fastpca <- function(expr, scale=FALSE) {
    k <- 50
    svd_raw <- svds(scale(t(expr), center=TRUE, scale=scale), k=k, nu=k, nv=0)
    pc_raw <- svd_raw$u %*% diag(svd_raw$d[1:k])
    return(pc_raw)
  }

  norm10x = function (ei)
  {
    sums = colSums(ei)
    eo = t(t(ei)*median(sums)/sums)
    return(eo)
  }

  tc <- norm10x(as.matrix(exprs(pbmc[use_genes,])))
  pc_tc <- fastpca(log1p(tc[rownames(dense),]), scale=TRUE)

  ## t-sne
  library(Rtsne)
  tsne_data <- Rtsne(pc_tc)

  load("tenx_pbmc/zinb_10.rda")
  tsne_zinb10 <- Rtsne(getW(zinb_10), pca = FALSE)

  ## k-means
  library(cluster)
  km <- kmeans(pc_tc, centers=10)

  ## ZINB-WaVE Clustering
  library(Seurat)
  W <- getW(zinb_10)
  pbmc.data <- Read10X("tenx_pbmc/outs/filtered_gene_bc_matrices/hg19/")
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 1, min.genes = 1,
                             project = "10X_PBMC")

  ## Build SNN
  k.param = 10
  k.scale = 10
  n.cells = NCOL(pbmc.data)
  data.use <- W
  my.knn <- FNN::get.knn(as.matrix(data.use), k = min(k.scale * k.param, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
  nn.large <- my.knn$nn.index

  w <- Seurat:::CalcSNNSparse(cell.names = colnames(pbmc.data),
                              k.param = k.param,
                              nn.large = nn.large,
                              nn.ranked = nn.ranked,
                              prune.SNN = 1/15,
                              print.output = FALSE)

  pbmc@snn <- w

  ## Run modularity clustering
  resolution <- 0.6
  pbmc <- Seurat:::RunModularityClustering(object = pbmc, SNN = w,
                                           modularity = 1, resolution = resolution,
                                           algorithm = 1, n.start = 100,
                                           n.iter = 10, random.seed = 0,
                                           print.output = FALSE, temp.file.location = NULL)
  pbmc <- Seurat:::GroupSingletons(pbmc, pbmc@snn)
  name <- paste("res.", resolution, sep = "")
  pbmc <- StashIdent(pbmc, name)

  z_cl <- pbmc@ident
  names(z_cl) <- paste0(names(pbmc@ident), "-1")

  table(z_cl)

  ## PCA clustering
  pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 1, min.genes = 1,
                             project = "10X_PBMC")

  ## Build SNN
  data.use <- pc_tc[,1:10]
  my.knn <- FNN::get.knn(as.matrix(data.use), k = min(k.scale * k.param, n.cells - 1))
  nn.ranked <- cbind(1:n.cells, my.knn$nn.index[, 1:(k.param-1)])
  nn.large <- my.knn$nn.index

  w <- Seurat:::CalcSNNSparse(cell.names = colnames(pbmc.data),
                              k.param = k.param,
                              nn.large = nn.large,
                              nn.ranked = nn.ranked,
                              prune.SNN = 1/15,
                              print.output = FALSE)

  pbmc@snn <- w

  ## Run modularity clustering
  resolution <- 0.2
  pbmc <- Seurat:::RunModularityClustering(object = pbmc, SNN = w,
                                           modularity = 1, resolution = resolution,
                                           algorithm = 1, n.start = 100,
                                           n.iter = 10, random.seed = 0,
                                           print.output = FALSE, temp.file.location = NULL)
  pbmc <- Seurat:::GroupSingletons(pbmc, pbmc@snn)
  name <- paste("res.", resolution, sep = "")
  pbmc <- StashIdent(pbmc, name)

  pca_cl <- pbmc@ident
  names(pca_cl) <- paste0(names(pbmc@ident), "-1")

  save(pc_tc, tc, km, tsne_data, dense, vargenes, tsne_zinb10,
       z_cl, pca_cl, file="tenx_pbmc/kmeans_res.rda")
}

if(analysis) {
  ## sequential clustering
  library(clusterExperiment)
  load("tenx_pbmc/zinb_10.rda")
  seW <- SummarizedExperiment(t(getW(zinb_10)))

  set.seed(4364236)
  cs <- clusterSingle(seW, subsample=FALSE, sequential=TRUE,
                      seqArgs = list(k0 = 15, beta = 0.95),
                      mainClusterArgs = list(clusterFunction = "kmeans"))
  cs
  cs_cl <- factor(primaryClusterNamed(cs), levels=c(-1, 1:17))

  sePCA <- SummarizedExperiment(t(pc_tc[,1:10]))
  set.seed(35235)
  cs2 <- clusterSingle(sePCA, subsample=FALSE, sequential=TRUE,
                      seqArgs = list(k0 = 15, beta = 0.9),
                      mainClusterArgs = list(clusterFunction = "kmeans"))
  cs2
  cs2_cl <- factor(primaryClusterNamed(cs2), levels=c(-1, 1:12))

  sePCA <- SummarizedExperiment(t(pc_tc[,1:50]))
  set.seed(352)
  cs3 <- clusterSingle(sePCA, subsample=FALSE, sequential=TRUE,
                       seqArgs = list(k0 = 15, beta = 0.8),
                       mainClusterArgs = list(clusterFunction = "kmeans"))
  cs3
  cs3_cl <- factor(primaryClusterNamed(cs3), levels=c(-1, 1:22))

  save(cs_cl, cs2_cl, cs3_cl, file="tenx_pbmc/seq_cluster.rda")
}

## color palette
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"))

## t-sne
load("tenx_pbmc/kmeans_res.rda")
load("tenx_pbmc/zinb_res.rda")
load("tenx_pbmc/zinb_10.rda")

## 10X genomics clusters
clusters <- read.csv("tenx_pbmc/clusters.csv")
cl <- clusters[,2]
names(cl) <- clusters[,1]
clus <- rep(NA, NCOL(tc))
names(clus) <- colnames(tc)
clus[names(cl)] <- cl
table(clus, z_cl)

plot(tsne_data$Y, pch=20, col=pal[clus])
plot(tsne_data$Y, pch=20, col=pal[z_cl])
plot(tsne_data$Y, pch=20, col=pal[pca_cl])

plot(getW(zinb_res), col=pal[z_cl], pch=20, xlab="W1", ylab="W2")
legend("topright", levels(z_cl), fill=pal, cex=.5)

plot(tsne_zinb10$Y, col=pal[z_cl], pch=20)
legend("topright", levels(z_cl), fill=pal, cex=.5)

## zinb-wave
plot(getW(zinb_res), col=pal[clus], pch=20)

pdf("pbmc68k_W2.pdf")
plot(getW(zinb_res), col=pal[z_cl], pch=20, xlab="W1", ylab="W2")
dev.off()

pdf("pbmc68k_W10_tsne.pdf")
plot(tsne_zinb10$Y, col=pal[z_cl], pch=20)
dev.off()

pdf("pbmc68k_W10_tsne_legend.pdf")
plot(tsne_zinb10$Y, col=pal[z_cl], pch=20)
legend("topright", levels(z_cl), fill=pal, cex=.5)
dev.off()

plot(getW(zinb_res), col=pal[pca_cl], pch=20)
plot(tsne_zinb10$Y, col=pal[pca_cl], pch=20)

## PCA
plot(pc_tc[,1:2], col=pal[clus], pch=20)
plot(pc_tc[,1:2], col=pal[z_cl], pch=20)
plot(pc_tc[,1:2], col=pal[pca_cl], pch=20)

## Let's look at Marker genes in z_cl clustering
table(z_cl)
library(limma)
pipestance_path <- "tenx_pbmc"
pbmc <- load_cellranger_matrix(pipestance_path)

y <- log1p(tc[rownames(dense),])
rownames(y) <- fData(pbmc)[rownames(y),2]
design <- model.matrix(~z_cl - 1)
fit <- lmFit(y, design)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(z_cl18 -
                                   (z_cl0 + z_cl1 + z_cl2 + z_cl3 + z_cl4 + z_cl5 +
                                    z_cl6 + z_cl7 + z_cl8 + z_cl9)/10, levels=design)

contrast.matrix <- makeContrasts((z_cl9 + z_cl10)/2 -
                                   (z_cl13 + z_cl11)/2, levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

top <- topTable(fit2, number=50)
head(top[order(top[,1], decreasing = TRUE),], 20)

## Cl0, 1: CD4+ T-cells (31%)
## Cl2,3,4,6: CD8+ T-cells (40%)
## Cl5: Activated cytotoxic T-cells (9%)
## Cl7: NK cells (6%)
## Cl8: B-cells (5%)
## Cl9: CD14+ Monocytes (CD14, S100A8, S100A9, LYZ) (3%)
## Cl10: CD16+ monocytes (FCGR3A/CD16, AIF1, FTL, LST1, SERPINA1)  (2.5%)
## Cl11: CD1C+ Dendritic Cells (CD1C, LYZ, HLA-DRA, HLA-DPA1, HLA-DPB1) (1%)
## Cl12: ? (0.5%)
## Cl13: pDC (plasmacytoid DC; GZMB, SERPINF1, ITM2C) (0.5%)
## Cl14: Enriched for mytochondrial genes (0.3%)
## Cl15: Megakatyocytes (PF4) (0.2%)
## Cl16: ? (0.2%)
## Cl17, Cl18: Enriched for ribosomal proteins (0.1%)

## 80% t-cells, 6% NK cells, 5% b-cells, 7% myeloid cells
## overall recapitulate results, but with finer details
## especially we are able to identify sub-pop of dendritic cells


df %>%
  ggplot(aes(x = Dim1, y = Dim2, color = CD1C)) +
  geom_point() + scale_color_continuous(low = "grey", high = "red")

##Let's look at plots with marker genes
boxplot(y["CD3D",]~z_cl)
boxplot(y["CD8A",]~z_cl)
boxplot(y["NKG7",]~z_cl)
boxplot(y["S100A8",]~z_cl)

tc_sym <- tc
rownames(tc_sym) <- fData(pbmc)[rownames(tc),2]

cl_anno <- as.character(z_cl)
cl_anno[z_cl == 0] <- "CD4+ T-cells (1)"
cl_anno[z_cl == 1] <- "CD4+ T-cells (2)"
cl_anno[z_cl == 2] <- "CD8+ T-cells (1)"
cl_anno[z_cl == 3] <- "CD8+ T-cells (2)"
cl_anno[z_cl == 4] <- "CD8+ T-cells (3)"
cl_anno[z_cl == 6] <- "CD8+ T-cells (4)"
cl_anno[z_cl == 5] <- "Activated T-cells"
cl_anno[z_cl == 7] <- "NK cells"
cl_anno[z_cl == 8] <- "B-cells"
cl_anno[z_cl == 9] <- "CD14+ monocytes"
cl_anno[z_cl == 10] <- "CD16+ monocytes"
cl_anno[z_cl == 11] <- "CD1C+ DCs"
cl_anno[z_cl == 13] <- "pDCs"
cl_anno[z_cl == 15] <- "Megakaryocytes"
cl_anno[z_cl == 12] <- "Uncharacterized (1)"
cl_anno[z_cl == 14] <- "Uncharacterized (2)"
cl_anno[z_cl == 16] <- "Uncharacterized (3)"
cl_anno[z_cl == 17] <- "Uncharacterized (4)"
cl_anno[z_cl == 18] <- "Uncharacterized (5)"

library(magrittr)
library(cowplot)

z_cl2 <- factor(paste0("Cluster", z_cl), levels=paste0("Cluster", 0:18))
df1 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                 Cluster=z_cl2)
df1 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb

save_plot("pbmc68k_tsne_zinb10.png", fig_tsne_zinb)

df2 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                  Cluster=cl_anno)
df2 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb_anno

save_plot("pbmc68k_tsne_zinb10_anno.png", fig_tsne_zinb_anno, base_aspect_ratio = 1.6)

## PCA only finds 10 clusters

## Cl0,2: CD8+ T-cells
## Cl1: CD4+ T-cells
## Cl3: activated cytotoxic T-cells
## Cl4: NK cells
## Cl5: B-cells
## Cl6: Myeloid cells
## Cl7: Dendritic cells
## Cl8: pDC?
## Cl9: Megakaryocyte

## PCA0: Cl0,2,6
## PCA1: Cl1,4,12
## PCA2: Cl3,14
## PCA3: Cl5
## PCA4: Cl7
## PCA5: Cl8,16
## PCA6: Cl9,11
## PCA7: Cl10
## PCA8: Cl12,13
## PCA9: Cl15

pca_anno <- as.character(pca_cl)
pca_anno[pca_cl == 0] <- "CD4+/CD8+ T-cells (1)"
pca_anno[pca_cl == 2] <- "CD4+/CD8+ T-cells (2)"
pca_anno[pca_cl == 1] <- "CD4+/CD8+ T-cells (3)"
pca_anno[pca_cl == 3] <- "Activated T-cells"
pca_anno[pca_cl == 4] <- "NK cells"
pca_anno[pca_cl == 5] <- "B-cells"
pca_anno[pca_cl == 6] <- "Myeloid cells"
pca_anno[pca_cl == 7] <- "Dendritic cells (1)"
pca_anno[pca_cl == 8] <- "Dendritic cells (2)"
pca_anno[pca_cl == 9] <- "Megakaryocytes"

df3 <- data.frame(Dim1=tsne_data$Y[,1], Dim2=tsne_data$Y[,2],
                 Cluster=paste0("Cl", pca_cl))

df3 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_pca

save_plot("pbmc68k_tsne_pca.png", fig_tsne_pca)

df4 <- data.frame(Dim1=tsne_data$Y[,1], Dim2=tsne_data$Y[,2],
                  Cluster=pca_anno)

df4 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_pca_anno

save_plot("pbmc68k_tsne_pca_anno.png", fig_tsne_pca_anno, base_aspect_ratio = 1.6)

design <- model.matrix(~pca_cl - 1)
fit <- lmFit(y, design)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(pca_cl8 -
                                   (pca_cl6 + pca_cl7)/2, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

top <- topTable(fit2, number=50)
head(top[order(top[,1], decreasing = TRUE),], n=10)

table(z_cl, cs_cl)

plot(tsne_zinb10$Y, col=pal[z_cl], pch=20)
plot(tsne_zinb10$Y, col=c("white", pal)[cs_cl], pch=20)
