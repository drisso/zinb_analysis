library(cowplot)
library(RColorBrewer)
library(magrittr)
library(ggplot2)

load("espresso_covariates.rda")

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

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

detection_rate <- colSums(all.counts>0)
coverage <- colSums(all.counts)

library(scater)
rownames(targets) <- colnames(all.counts)
all.counts = as.matrix(all.counts)
sceset <- SingleCellExperiment(assays = list(counts = all.counts),
                               colData = targets)

keep_feature <- rowSums(assay(sceset) > 0) > 0
sceset <- sceset[keep_feature,]

sceset <- calculateQCMetrics(sceset)
qc <- colData(sceset)[,c(4, 6, 8:11)]
pct_dropout <- colMeans(assay(sceset) == 0)
qc$pct_dropout <- pct_dropout
qc = as.matrix(qc)

filter <- rowSums(all.counts>10)>=10
raw <- all.counts[filter,]

colMerged <- col1[targets$Serum]
colBatch <- col2[targets$Batch]
level1 <- targets$Serum
level2 <- targets$Batch

condition <- level1
batch <- level2

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2], batch = batch) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel1_zinb

data.frame(Dim1=zinb_batch@W[,1], Dim2=zinb_batch@W[,2], batch = batch) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel2_zinb

library(cluster)

methods_sub <- list("ZINB-WaVE"=zinb@W, "ZINB-Batch"=zinb_batch@W)
sil_cond <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  s_cond <- silhouette(as.numeric(condition), d)
  ss_cond <-  summary(s_cond)
  return(ss_cond$clus.avg.widths)
})

bars <- data.frame(AverageSilhouette=unlist(sil_cond),
                   Method=rep(names(methods_sub), each=nlevels(condition)),
                   Cluster=factor(rep(levels(condition), length(methods_sub)),
                                  levels=levels(condition)))

bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil

sil_batch <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  s_cond <- silhouette(as.numeric(batch), d)
  ss_cond <-  summary(s_cond)
  return(ss_cond$clus.avg.widths)
})

bars <- data.frame(AverageSilhouette=unlist(sil_batch),
                   Method=rep(names(methods_sub), each=nlevels(batch)),
                   Cluster=factor(rep(paste0("Batch", levels(batch)), length(methods_sub)),
                                  levels=paste0("Batch", levels(batch))))

bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil2

load("patel_covariates.rda")

counts <- read.table("patel/glioblastoma_raw_rnaseq_SCandbulk_counts.txt", header=TRUE, stringsAsFactors = FALSE, row.names=NULL)

gene_symbols <- counts[,1]
ensembl_ids <- counts[,2]
sample_names <- colnames(counts)[-(1:2)]

all.counts <- counts[,-(1:2)]
rownames(all.counts) <- ensembl_ids

metadata <- read.table("patel/SraRunTable.txt", sep='\t', stringsAsFactors = FALSE, header=TRUE, row.names=5, na.strings = "<not provided>")
metadata <- metadata[sample_names,]

# select only single-cell samples from patients
keep <- which(grepl("^Single cell", metadata$source_name_s) &
                !is.na(metadata$patient_id_s) &
                !is.na(metadata$subtype_s))
metadata <- metadata[keep,]

all.counts <- all.counts[,keep]

stopifnot(all(rownames(metadata)==colnames(all.counts)))

level1 <- as.factor(metadata$patient_id_s[!is.na(metadata$subtype_s)])

col3 <- col2[-(1:2)]

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch)) + geom_point()  +
  scale_color_manual(values = col3) -> patel_zinb

data.frame(Dim1=zinb_cdr@W[,1], Dim2=zinb_cdr@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch)) + geom_point()  +
  scale_color_manual(values = col3) -> patel_zinb2

fig2_tris <- plot_grid(panel1_zinb + theme(legend.position = "none"),
                       panel2_zinb,
                       sil, sil2,
                       patel_zinb + theme(legend.position = "none"),
                       patel_zinb2,
                       labels=c("a", "b", "c", "d", "e", "f"),
                       ncol=2, nrow=3, rel_widths = c(1, 1.25, 1, 1, 1, 1.25))
fig2_tris

save_plot("espresso_fig5.pdf", fig2_tris,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1
)

table(batch)
detection_rate <- colSums(all.counts>0)

df <- data.frame(detection_rate = detection_rate, batch = batch,
                 PC1 = pc_combat_fq[,1], PC2 = pc_combat_fq[,2],
                 W1 = zinb_batch@W[,1], W2 = zinb_batch@W[,2])

ggplot(df, aes(x = batch, y = detection_rate, fill = batch)) +
  geom_boxplot() + scale_fill_manual(values = col3) -> panel_c

ggplot(df, aes(PC1, PC2, colour=batch)) + geom_point()  +
  scale_color_manual(values = col3) -> panel_a

ggplot(df, aes(W1, W2, colour=batch)) + geom_point()  +
  scale_color_manual(values = col3) -> panel_b

fig_sup <- plot_grid(panel_a + theme(legend.position = "none"),
                     panel_b + theme(legend.position = "none"),
                     panel_c + theme(axis.text.x = element_text(angle = 30, hjust = 1))
                     , ncol=3, nrow=1, rel_widths = c(1, 1, 1.25), labels="auto")

fig_sup

save_plot("patel_combat.pdf", fig_sup,
          base_aspect_ratio = 3
)
