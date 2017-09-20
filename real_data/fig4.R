## color palette
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"))

## t-sne
load("tenx_pbmc/kmeans_res.rda")
load("tenx_pbmc/zinb_res.rda")
load("tenx_pbmc/zinb_10.rda")
load("tenx_pbmc/seq_cluster.rda")

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
pca_cl2 <- factor(paste0("Cluster", pca_cl), levels=paste0("Cluster", 0:9))

df1 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                  W1=getW(zinb_res)[,1], W2=getW(zinb_res)[,2],
                  PC1=pc_tc[,1], PC2=pc_tc[,2],
                  Cluster=z_cl2)
df1 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb

save_plot("zinb_tsne_pbmc.png", fig_tsne_zinb, base_aspect_ratio = 1.5, base_height = 8)

df1 %>%
  ggplot(aes(x = W1, y = W2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_W

df1 %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_pca

df2 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                  W1=getW(zinb_res)[,1], W2=getW(zinb_res)[,2],
                  PC1=pc_tc[,1], PC2=pc_tc[,2],
                  Cluster=cl_anno)

df2 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb_anno

df2 %>%
  ggplot(aes(x = W1, y = W2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_W_anno

df2 %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_pca_anno

sup_fig <- plot_grid(fig_W + theme(legend.position = "none"),
          fig_pca, nrow = 1, ncol=2, rel_widths = c(1, 1.5), labels = "auto")
save_plot("Sup_pbmc.png", sup_fig, base_aspect_ratio = 1.5, base_height = 8)

df3 <- data.frame(Dim1=tsne_zinb10$Y[cs_cl != "-1",1], Dim2=tsne_zinb10$Y[cs_cl != "-1",2],
                  Cluster=droplevels(cs_cl[cs_cl != "-1"]))

df3 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb_seq
save_plot("Sup_pbmc_seq.png", fig_tsne_zinb_seq, base_aspect_ratio = 1.5, base_height = 8)

## heatmap of marker genes
pipestance_path <- "tenx_pbmc"
pbmc <- cellrangerRkit::load_cellranger_matrix(pipestance_path)
# use_genes <- cellrangerRkit::get_nonzero_genes(pbmc)
# epbmc <- as.matrix(exprs(pbmc[use_genes,]))

rownames(dense) <- fData(pbmc)[rownames(dense), 2]
rownames(tc) <- fData(pbmc)[rownames(tc), 2]

se <- SummarizedExperiment(dense)
zw <- zinbwave(se, fitted_model = zinb_10)

norm <- assay(zw, "normalizedValues")
means <- t(apply(norm, 1, tapply, z_cl, mean))
means2 <- t(apply(tc, 1, tapply, z_cl, mean))

genes <- c("CD3D", "CD8A", "NKG7", "CD79A", "CD14", #"S100A9", "LYZ",
           "FCGR3A", #"AIF1", "FTL", "LST1",
           "CD1C", #"HLA.DR", "GZMB",
           "SERPINF1", #"ITM2C",
           "PF4")
gg <- genes[which(genes %in% rownames(means))]
gg2 <- genes[which(genes %in% rownames(means2))]

library(pheatmap)
cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
pdf("pbmc68k_heatmap.pdf")
pheatmap(log1p(means2)[gg2,], cluster_cols = FALSE, cluster_rows = FALSE, scale="row", color = cols)
dev.off()

## Silhouette width (ZINB)
W <- getW(zinb_10)
centers <- apply(W, 2, tapply, z_cl, mean)
f <- function(i) {
  d <- sapply(levels(z_cl), function(c) {
    sqrt(sum((W[i,] - centers[c,])^2))
  })
  ai <- d[as.character(z_cl[i])]
  bi <- min(d[names(d) != as.character(z_cl[i])])
  s <- (bi - ai)/max(ai, bi)
  return(as.numeric(s))
}

s <- sapply(seq_len(nrow(W)), f)
mean(s)
tapply(s, z_cl, mean)

## Silhouette width (PCA)
library(Seurat)
pbmc.data <- Read10X("tenx_pbmc/outs/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 1, min.genes = 1,
                           project = "10X_PBMC")

## Build SNN
k.param = 10
k.scale = 10
n.cells = NCOL(pbmc.data)
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
resolution <- 0.6
pbmc <- Seurat:::RunModularityClustering(object = pbmc, SNN = w,
                                         modularity = 1, resolution = resolution,
                                         algorithm = 1, n.start = 100,
                                         n.iter = 10, random.seed = 0,
                                         print.output = FALSE, temp.file.location = NULL)
pbmc <- Seurat:::GroupSingletons(pbmc, pbmc@snn)
name <- paste("res.", resolution, sep = "")
pbmc <- StashIdent(pbmc, name)

pca2_cl <- pbmc@ident
names(pca2_cl) <- paste0(names(pbmc@ident), "-1")

table(pca2_cl, pca_cl)
table(pca2_cl, z_cl)

P <- pc_tc[,1:10]
centers <- apply(P, 2, tapply, pca_cl, mean)
g <- function(i) {
  d <- sapply(levels(pca_cl), function(c) {
    sqrt(sum((P[i,] - centers[c,])^2))
  })
  ai <- d[as.character(pca_cl[i])]
  bi <- min(d[names(d) != as.character(pca_cl[i])])
  s <- (bi - ai)/max(ai, bi)
  return(as.numeric(s))
}

sp <- sapply(seq_len(nrow(P)), g)
mean(sp)
tapply(sp, pca_cl, mean)

means3 <- t(apply(tc, 1, tapply, pca_cl, mean))
pheatmap(log1p(means3)[gg2,], cluster_cols = FALSE, cluster_rows = FALSE, scale="row", color = cols)

means4 <- t(apply(tc, 1, tapply, pca2_cl, mean))
pheatmap(log1p(means4)[gg2,], cluster_cols = FALSE, cluster_rows = FALSE, scale="row", color = cols)
