library(cowplot)
library(RColorBrewer)
library(magrittr)
library(zinbwave)
library(cellrangerRkit)
library(Seurat)

## Clustering with ZINB-WaVE
load("tenx_pbmc4k/zinb_eps.rda")
W <- getW(zinb_res3)

pbmc.data <- Read10X("tenx_pbmc4k/outs/filtered_gene_bc_matrices/GRCh38/")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 1, min.genes = 1,
                           project = "10X_PBMC")

## Build SNN
k.param = 10
k.scale = 10
n.cells = NCOL(pbmc.data)
data.use <- getW(zinb_res3)
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

z_cl <- pbmc@ident
names(z_cl) <- paste0(names(pbmc@ident), "-1")

table(z_cl)

pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"))

plot(W, pch=19, col=pal[z_cl])

#############
## PCA clustering
load("tenx_pbmc4k/pca_res.rda")
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 1, min.genes = 1,
                           project = "10X_PBMC")

## Build SNN
data.use <- pc_tc[,1:50]
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

table(pca_cl)
table(pca_cl, z_cl)

plot(pc_tc[,1:2], pch=19, col=pal[pca_cl])
plot(W, pch=19, col=pal[pca_cl])

### Original clustering
clusters <- read.csv("tenx_pbmc4k/clusters.csv")
cl <- clusters[,2]
names(cl) <- clusters[,1]

plot(pc_tc[,1:2], pch=19, col=pal[cl])
plot(W, pch=19, col=pal[cl])
table(cl, z_cl)

## redo with clusterExperiment?
