library(scone)
library(slingshot)
data_dir <- "/Users/dar2062/git/bioc2017singlecell/data/"
load(paste0(data_dir, "GSE95601_oeHBCdiff_Cufflinks_eSet_reduced.Rda"))

# Count matrix
E <- assayData(Cufflinks_eSet)$counts_table

# Remove undetected genes
E <- na.omit(E)
E <- E[rowSums(E)>0,]

# Remove ERCC and CreER genes
cre <- E["CreER",]
ercc <- E[grep("^ERCC-", rownames(E)),]
E <- E[grep("^ERCC-", rownames(E), invert = TRUE), ]
E <- E[-which(rownames(E)=="CreER"), ]

# Extract QC metrics
qc <- as.matrix(protocolData(Cufflinks_eSet)@data)[,c(1:5, 10:18)]
qc <- cbind(qc, CreER = cre, ERCC_reads = colSums(ercc))

# Extract metadata
batch <- droplevels(pData(Cufflinks_eSet)$MD_c1_run_id)
bio <- droplevels(pData(Cufflinks_eSet)$MD_expt_condition)
clusterLabels <- read.table(paste0(data_dir, "oeHBCdiff_clusterLabels.txt"),
                            sep = "\t", stringsAsFactors = FALSE)
m <- match(colnames(E), clusterLabels[, 1])

# Create metadata data.frame
metadata <- data.frame("Experiment" = bio,
                       "Batch" = batch,
                       "publishedClusters" = clusterLabels[m,2],
                       qc)

# Symbol for cells not assigned to a lineage in original data
metadata$publishedClusters[is.na(metadata$publishedClusters)] <- -2

se <- SummarizedExperiment(assays = list(counts = E),
                           colData = metadata)
se

# QC-metric-based sample-filtering
data("housekeeping")
hk = rownames(se)[toupper(rownames(se)) %in% housekeeping$V1]

mfilt <- metric_sample_filter(assay(se),
                              nreads = colData(se)$NREADS,
                              ralign = colData(se)$RALIGN,
                              pos_controls = rownames(se) %in% hk,
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)

mfilt <- !apply(simplify2array(mfilt[!is.na(mfilt)]), 1, any)
se <- se[, mfilt]
dim(se)

publishedClusters <- colData(core)[, "publishedClusters"]
col_clus <- c("transparent", "#1B9E77", "antiquewhite2", "cyan", "#E7298A",
              "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2",
              "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")
names(col_clus) <- sort(unique(publishedClusters))
table(publishedClusters)

cl2 <- as.factor(publishedClusters[!publishedClusters %in% c("-2")])
pal2 <- c("#1B9E77", "antiquewhite2", "cyan", "#E7298A",
  "#A6CEE3", "#666666", "#E6AB02", "#FFED6F", "darkorchid2",
  "#B3DE69", "#FF7F00", "#A6761D", "#1F78B4")

# Normalization
fq <- FQ_FN(assay(se))
pc_fq <- prcomp(t(log1p(fq)))
X_pca <- pc_fq$x[!publishedClusters %in% c("-2"),1:3]
plot(X_pca, pch=19, col=pal2[cl2])
lineages_pca <- getLineages(X_pca, clusterLabels = cl2, start.clus = "1",
                            end.clus = c("12", "7", "15"))
lineages_pca
pairs(lineages_pca, type="lineages", col = pal2[cl2], show.constraints = TRUE, pch=19)
