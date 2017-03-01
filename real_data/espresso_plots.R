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
sceset <- newSCESet(countData = all.counts, phenoData = AnnotatedDataFrame(targets))

keep_feature <- rowSums(exprs(sceset) > 0) > 0
sceset <- sceset[keep_feature,]

sceset <- calculateQCMetrics(sceset)
qc <- pData(sceset)[,c(4, 7, 10, 14:17)]

filter <- rowSums(all.counts>10)>=10
raw <- all.counts[filter,]

colMerged <- col1[targets$Serum]
colBatch <- col2[targets$Batch]
level1 <- targets$Serum
level2 <- targets$Batch

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> panel1_pca

data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> panel1_zifa

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point()  +
  scale_color_brewer(palette="Set1")  -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
          panel1_zifa + theme(legend.position = "none"),
          panel1_zinb + theme(legend.position = "none"),
          labels=c("A", "C", "E"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, .6))

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow") -> panel2_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow") -> panel2_zifa
data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2)) + geom_point(aes(color=detection_rate)) + scale_colour_gradient(low="blue", high="yellow") -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("B", "D", "F"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, .6))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("epresso_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

cors <- lapply(1:2, function(i) abs(cor(pc_tc[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_pca

cors <- lapply(1:2, function(i) abs(cor(zifa_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zifa

cors <- lapply(1:2, function(i) abs(cor(zinb@W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("B", "D", "F"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1bis <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("espresso_fig1bis.pdf", fig1bis,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)


library(cluster)

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB"))

sil_cl <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(level1), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_cl, Method=names(methods), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Type)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none") -> sil

sil2 <- plot_grid(sil, NULL, NULL, ncol=3, nrow=1, labels="G")
fig1_tris <- plot_grid(upper, lower, sil2, ncol=1, nrow=3)
fig1_tris

save_plot("espresso_fig1tris.pdf", fig1_tris,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(level1), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB"), each=nlevels(level1)),
                   Cluster=factor(rep(levels(level1), length(methods_sub)),
                                  levels=levels(level1)))

library(dplyr)
bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("A", "C", "E"), align = "h", ncol=3)

upper <- plot_grid(p1, sil, labels=c("", "G"), rel_widths = c(3, 1))

fig1_4 <- plot_grid(upper, lower, ncol=1, nrow=2)
fig1_4

save_plot("espresso_fig1_v4.pdf", fig1_4,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

save_plot("espresso_supp_sil.pdf", sil)

condition <- level1
batch <- level2

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel1_zinb

data.frame(Dim1=zinb_batch@W[,1], Dim2=zinb_batch@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel2_zinb

fig2 <- plot_grid(panel1_zinb + theme(legend.position = "none"),
                  panel2_zinb,
                  labels=c("A", "B"), ncol=2, nrow=1, rel_widths = c(1, 1.25))
fig2

save_plot("espresso_fig2.pdf", fig2,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.3
)

methods_sub <- list(ZINB=zinb@W, ZINB_BATCH=zinb_batch@W)
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
                   Cluster=factor(rep(levels(condition), length(methods_sub)),
                                  levels=levels(condition)))

bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil
