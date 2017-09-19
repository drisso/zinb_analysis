library(clusterExperiment)
library(slingshot)
library(zinbwave)

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


## PCA
their_cl <- factor(colData(seObj)$publishedClusters)
cl2 <- droplevels(their_cl[!their_cl %in% c("-2")])
pal2 <- c("#1B9E77", "antiquewhite2", "cyan", "#E7298A",
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2",
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")

X2 <- W[!their_cl %in% c("-2"),]
X2 <- cmdscale(dist(X2), k = 3)
plot(X2, pch=19, col=pal2[cl2])

lineages2 <- getLineages(X2, clusterLabels = cl2, start.clus = "1", end.clus = "")
pairs(lineages2, type="lineages", col = pal2[cl2], show.constraints = TRUE, pch=19)
lineages2
