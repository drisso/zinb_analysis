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
          labels=c("a", "c", "e"), align = "h", ncol=3)
p1

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
                labels=c("b", "d", "f"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, .6))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)
fig1

save_plot("epresso_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

cors <- lapply(1:2, function(i) abs(cor(pc_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 0.8) -> panel2_pca

cors <- lapply(1:2, function(i) abs(cor(zifa_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, .8) -> panel2_zifa

cors <- lapply(1:2, function(i) {
  yy <- tapply(getW(zinb)[,i], level1, identity)
  apply(qc, 2, function(x) {
    xx <- tapply(x, level1, identity)
    mean(sapply(seq_along(xx), function(i) abs(cor(xx[[i]], yy[[i]]))))
  })
})
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 2),
                   Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, .8) -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("b", "d", "f"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1bis <- plot_grid(upper, lower, ncol=1, nrow=2)
fig1bis

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
                    "ZINB-Wave")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-Wave"))

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
                   Method=rep(c("PCA", "ZIFA", "ZINB-Wave"), each=nlevels(level1)),
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
                labels=c("a", "c", "e"), align = "h", ncol=3)

upper <- plot_grid(p1, sil, labels=c("", "g"), rel_widths = c(3, 1))

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

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2], batch = batch) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel1_zinb

data.frame(Dim1=zinb_batch@W[,1], Dim2=zinb_batch@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel2_zinb

fig2 <- plot_grid(panel1_zinb + theme(legend.position = "none"),
                  panel2_zinb,
                  labels=c("a", "b"), ncol=2, nrow=1, rel_widths = c(1, 1.25))
fig2

save_plot("espresso_fig2.pdf", fig2,
          ncol = 2,
          nrow = 1,
          base_aspect_ratio = 1.3
)

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

fig2_bis <- plot_grid(panel1_zinb + theme(legend.position = "none"),
                  panel2_zinb,
                  sil, sil2,
                  labels=c("a", "b", "c", "d"), ncol=2, nrow=2, rel_widths = c(1, 1.25, 1, 1))
fig2_bis

save_plot("espresso_fig2bis.pdf", fig2_bis,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3
)

# Compare with combat
methods_sub <- list("ZINB-WaVE"=zinb@W,
                    "ZINB-Batch"=zinb_batch@W,
                    "PCA ComBat RAW" = pc_combat_raw,
                    "PCA ComBat TC" = pc_combat_tc,
                    "PCA ComBat TMM" = pc_combat_tmm,
                    "PCA ComBat FQ" = pc_combat_fq,
                    "PCA RAW" = pc_raw,
                    "PCA TC" = pc_tc,
                    "PCA TMM" = pc_tmm,
                    "PCA FQ" = pc_fq)

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


data.frame(Dim1=pc_combat_fq[,1], Dim2=pc_combat_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel_pc_combat_fq

data.frame(Dim1=pc_fq[,1], Dim2=pc_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=batch, shape=condition)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel_pc_fq

supp_combat <- plot_grid(panel_pc_fq + theme(legend.position = "none"),
                         panel_pc_combat_fq,
                         sil, sil2,
                      labels=c("a", "b", "c", "d"), ncol=2, nrow=2, rel_widths = c(1, 1.25, 1, 1))
supp_combat

save_plot("espresso_supp_combat.pdf", supp_combat,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3
)


# ## goodness-of-fit
# library(matrixStats)
# totalcount = function (ei)
# {
#   sums = colSums(ei)
#   eo = t(t(ei)*mean(sums)/sums)
#   return(eo)
# }
#
# #-------------------- Functions
# myExp <- function(x, eps=1)
#   exp(x)-eps
#
# myLog <- function(x, eps=1)
#   log(x+eps)
#
# invLogit <- function(x)
#   exp(x)/(1 + exp(x))
#
# #----- MD-plot
# MD <- function(x, y, log=FALSE, pts=NULL, pch=20, col=2, smooth=TRUE, main="", ...)
# {
#   if(log)
#   {
#     m <- (myLog(x)+myLog(y))/2
#     d <- myLog(y)-myLog(x)
#   }
#   if(!log)
#   {
#     m <- (x+y)/2
#     d <- y-x
#   }
#   if(smooth)
#     smoothScatter(m,d,xlab="M",ylab="D",main=main,...)
#   if(!smooth)
#     plot(m,d,xlab="M",ylab="D",main=main,...)
#   lines(lowess(d ~ m), col=2, lwd=2)
#   abline(h=0)
#
#   if(!is.null(pts))
#     points(m[pts],d[pts],pch=pch,col=col)
# }
#
#
# tc <- totalcount(raw)
#
# vars <- rowVars(log1p(tc))
# names(vars) <- rownames(tc)
# vars <- sort(vars, decreasing = TRUE)
# vargenes <- names(vars)[1:1000]
#
# core <- raw[vargenes,]
#
# obs_mean <- rowMeans(log1p(core[,level1=="2i"]))
# obs_prop <- rowMeans(core[,level1=="2i"]>0)
#
# zinb_mu = getLogMu(zinb)[level1=="2i",]
# zinb_pi = getPi(zinb)[level1=="2i",]
#
# zinb_mean <- colMeans((1 - zinb_pi) * zinb_mu)
# zinb_prop <- colMeans(zinb_pi + (1 - zinb_pi) * (1 + getPhi(zinb)[1] * zinb_mu)^(1/getPhi(zinb)[1]))
#
# MD(obs_mean, zinb_mean, log=FALSE)
# MD(obs_pi, zinb_pi, log=FALSE)


# PCA and ZIFA for raw, TC, TMM, FQ
data.frame(Dim1=pc_raw[,1], Dim2=pc_raw[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_raw
data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_tc
data.frame(Dim1=pc_tmm[,1], Dim2=pc_tmm[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_tmm
data.frame(Dim1=pc_fq[,1], Dim2=pc_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> pca_fq

fig_pca <- plot_grid(pca_raw, pca_tc, pca_tmm, pca_fq, labels=c("a", "b", "c", "d"))
save_plot("espresso_supp_pca.pdf", fig_pca,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3)

data.frame(Dim1=zifa_raw[,1], Dim2=zifa_raw[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> zifa_raw
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> zifa_tc
data.frame(Dim1=zifa_tmm[,1], Dim2=zifa_tmm[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> zifa_tmm
data.frame(Dim1=zifa_fq[,1], Dim2=zifa_fq[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set1") -> zifa_fq

fig_zifa <- plot_grid(zifa_raw, zifa_tc, zifa_tmm, zifa_fq, labels=c("a", "b", "c", "d"))
save_plot("espresso_supp_zifa.pdf", fig_zifa,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3)

## new figure batch
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

