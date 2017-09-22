library(clusterExperiment)
library(slingshot)
library(zinbwave)
library(RColorBrewer)

data_dir <- "/Users/dar2062/git/bioc2017singlecell/data/"
load(sprintf('%sse_after_zinbwave.rda', data_dir))
load(sprintf('%sceObj_after_RSEC.rda', data_dir))
load(sprintf('%sseObj_after_RSEC.rda', data_dir))

palDF <- ceObj@clusterLegend[[1]]
pal <- palDF[, "color"]
names(pal) <- palDF[, "name"]
pal["-1"] = "transparent"

our_cl <- primaryClusterNamed(ceObj)
cl <- our_cl[!our_cl %in% c("-1", "c4")]
pal <- pal[!names(pal) %in% c("-1", "c4")]
pal <- brewer.pal(6, "Set2")
names(pal) <- paste0("c", c(1:3, 5:7))

## ZINB-Wave
W <- colData(se)[, grepl("^W", colnames(colData(se)))]
W <- as.matrix(W)

X <- W[!our_cl %in% c("-1", "c4"), ]
X <- cmdscale(dist(X), k = 3)

plot(X, pch=19, col=pal[cl])

lineages <- getLineages(X, clusterLabels = cl, start.clus = "c1")

pdf("oe_slingshot_zinbwave.pdf")
pairs(lineages, type="lineages", col = pal[cl], show.constraints = TRUE, pch=19)
dev.off()

pdf("oe_slingshot_zinbwave_2d.pdf")
plot(X[,c(1, 2)], pch=19, col=pal[cl], xlab="Dim1", ylab="Dim2", cex.lab=1.5, cex.axis=1.5)
lines(lineages, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()

## PCA
load("oeHBCdiff_finalClusterObject.Rda")
cl2 <- factor(primaryClusterNamed(cmobj2), levels=paste0("m", c(1:5, 7:12, 14, 15)))
pal2 <- c("#1B9E77", "antiquewhite2", "cyan", "#E7298A",
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2",
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")

pc <- prcomp(t(log1p(assay(cmobj2))))
X2 <- pc$x[,1:5]
lineages2 <- getLineages(X2, clusterLabels = cl2, start.clus = "m1", end.clus = c("m12", "m4", "m15"))

pdf("oe_slingshot_pca.pdf")
pairs(lineages2, type="lineages", col = pal2[cl2], show.constraints = TRUE, pch=19)
dev.off()
lineages2

pdf("oe_slingshot_pca_2d.pdf")
plot(X2[,c(1, 2)], pch=19, col=pal2[cl2], xlab="Dim1", ylab="Dim2", cex.lab=1.5, cex.axis=1.5)
lines(lineages2, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()

## no constraints
lineages3 <- getLineages(X2, clusterLabels = cl2, start.clus = "m1")

pdf("oe_slingshot_pca_noc.pdf")
pairs(lineages3, type="lineages", col = pal2[cl2], show.constraints = TRUE, pch=19)
dev.off()
lineages3

pdf("oe_slingshot_pca_2d_noc.pdf")
plot(X2[,c(1, 2)], pch=19, col=pal2[cl2], xlab="Dim1", ylab="Dim2", cex.lab=1.5, cex.axis=1.5)
lines(lineages3, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()

## no constraints -- zinbwave clusters
names(cl) <- rownames(X)
names(cl2) <- rownames(X2)
idx <- intersect(names(cl), names(cl2))
cl3 <- cl[idx]
X3 <- X2[idx,]
lineages4 <- getLineages(X3, clusterLabels = cl3, start.clus = "c1")

pdf("oe_slingshot_pca_noc_zinbcl.pdf")
pairs(lineages4, type="lineages", col = pal[cl3], show.constraints = TRUE, pch=19)
dev.off()
lineages4

pdf("oe_slingshot_pca_2d_noc_zinbcl.pdf")
plot(X3[,c(1, 2)], pch=19, col=pal[cl3], xlab="Dim1", ylab="Dim2")
lines(lineages4, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()

## 3-d
X4 <- pc$x[,1:3]
lineages5 <- getLineages(X4, clusterLabels = cl2, start.clus = "m1", end.clus = c("m12", "m4", "m15"))

pdf("oe_slingshot_pca3.pdf")
pairs(lineages5, type="lineages", col = pal2[cl2], show.constraints = TRUE, pch=19)
dev.off()
lineages5

pdf("oe_slingshot_pca3_2d.pdf")
plot(X4[,c(1, 2)], pch=19, col=pal2[cl2], xlab="Dim1", ylab="Dim2")
lines(lineages5, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()

## ZINB-WaVE (5 components)
X5 <- W[!our_cl %in% c("-1", "c4"), ]
X5 <- cmdscale(dist(X5), k = 5)

plot(X5, pch=19, col=pal[cl])

lineages6 <- getLineages(X5, clusterLabels = cl, start.clus = "c1")

pdf("oe_slingshot_zinbwave5.pdf")
pairs(lineages6, type="lineages", col = pal[cl], show.constraints = TRUE, pch=19)
dev.off()
lineages6

pdf("oe_slingshot_zinbwave5_2d.pdf")
plot(X[,c(1, 2)], pch=19, col=pal[cl], xlab="Dim1", ylab="Dim2")
lines(lineages6, type="lineages", show.constraints = TRUE, pch=19, dims = 1:2)
dev.off()
