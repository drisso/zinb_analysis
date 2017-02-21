library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
load("allen_covariates_1000.rda")

data("allen")
allen_core <- allen[grep("^ERCC-", rownames(allen), invert = TRUE),
                    which(colData(allen)$Core.Type=="Core")]

filter <- rowSums(assay(allen_core)>10)>=10
raw <- assay(allen_core)[filter,]

detection_rate <- colSums(raw>0)
coverage <- colSums(raw)
qc <- cbind(detection_rate, as.data.frame(colData(allen_core)[,1:15]))

layer <- as.factor(colData(allen_core)$driver_1_s)
cluster <- as.factor(colData(allen_core)$Primary.Type)
cluster2 <- as.character(cluster)
cluster2[cluster2 %in% c("Pvalb Tacr3", "Sst Myh8")] <- "Other"
cluster2[grep("^L4", cluster2)] <- "L4"
cluster2[grep("^L6", cluster2)] <- "L6a"
cluster2[grep("^L5a", cluster2)] <- "L5a"
cluster2[cluster2 == "L5 Ucma"] <- "L5a"
cluster2[grep("^L5b", cluster2)] <- "L5b"
cluster2 <- as.factor(cluster2)

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))
collayer <- col1[layer]

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=layer)) + geom_point() -> panel1_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=layer)) + geom_point() -> panel1_zifa
data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=layer)) + geom_point() -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("A", "B", "C"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, .6))

cors <- lapply(1:2, function(i) abs(cor(pc_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, .6) -> panel2_pca

cors <- lapply(1:2, function(i) abs(cor(zifa_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, .6) -> panel2_zifa

cors <- lapply(1:2, function(i) abs(cor(zinb@W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors, QC=rep(stringr::str_to_lower(colnames(qc)), 2), Dimension=as.factor(rep(1:2, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) + geom_bar(stat="identity", position='dodge') + scale_fill_manual(values=col2) + ylim(0, .6) -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("D", "E", "F"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("allen_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

## Alternatively

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster2)) + geom_point() + scale_color_manual(values=col2) -> panel1_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster2)) + geom_point() + scale_color_manual(values=col2) -> panel1_zifa
data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster2)) + geom_point() + scale_color_manual(values=col2) -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("A", "B", "C"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, 1))

fig1_bis <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("allen_fig1bis.pdf", fig1_bis,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

# ### t-sne
#
# library(Rtsne)
# tsne_pca <- Rtsne(pc_tc[,1:2], pca=FALSE)
# plot(tsne_pca$Y, col=col2[cluster], pch=19)
# tsne_zifa <- Rtsne(zifa_tc[,1:2], pca=FALSE)
# plot(tsne_zifa$Y, col=col2[cluster], pch=19)
# tsne_zinb <- Rtsne(zinb@W[,1:2], pca=FALSE)
# plot(tsne_zinb$Y, col=col2[cluster], pch=19)

### Silhouette

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
  ss <- silhouette(as.numeric(cluster), d)
  mean(ss[,3])
})

sil_lay <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(cluster2), d)
  mean(ss[,3])
})

pdf("allen_sil.pdf")
barplot(sil_lay, col=col1[met_type], horiz = TRUE, names.arg = names(methods), las=2, main="by layer")
dev.off()

bars <- data.frame(AverageSilhouette=sil_lay, Method=names(methods), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Type)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none") -> sil

sil2 <- plot_grid(sil, NULL, NULL, ncol=3, nrow=1, labels="G")
fig1_tris <- plot_grid(upper, lower, sil2, ncol=1, nrow=3)
fig1_tris

save_plot("allen_fig1tris.pdf", fig1_tris,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)
