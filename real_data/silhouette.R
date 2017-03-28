library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
library(cluster)
library(dplyr)

col1 <- brewer.pal(9, "Set1")
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))
mycol <- c(brewer.pal(11,"RdYlGn")[c(8:11, 1:4)], brewer.pal(11,"RdYlBu")[8:11])
colsil <- mycol[c(5:12, 2, 3)]
## Allen

load("allen_covariates_1000.rda")

data("allen")
allen_core <- allen[grep("^ERCC-", rownames(allen), invert = TRUE),
                    which(colData(allen)$Core.Type=="Core" &
                            !(colData(allen)$Primary.Type %in% c("Pvalb Tacr3", "Sst Myh8")))]
cluster2 <- colData(allen_core)$Primary.Type
cluster2[grep("^L4", cluster2)] <- "L4"
cluster2[grep("^L6", cluster2)] <- "L6a"
cluster2[grep("^L5a", cluster2)] <- "L5a"
cluster2[cluster2 == "L5 Ucma"] <- "L5"
cluster2[grep("^L5b", cluster2)] <- "L5b"
cluster2 <- as.factor(cluster2)

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB-WaVE")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-WaVE"))

sil_lay <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(cluster2), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_lay,
                   Method=factor(names(methods),
                                 levels=c("PCA_RAW", "PCA_TC", "PCA_TMM", "PCA_FQ",
                                          "ZIFA_RAW", "ZIFA_TC", "ZIFA_TMM", "ZIFA_FQ",
                                          "ZINB-WaVE")),
                   Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Method)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=colsil) + coord_flip() +
  theme(legend.position = "none") -> sil_allen

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(cluster2), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB-WaVE"), each=nlevels(cluster2)),
                   Cluster=rep(levels(cluster2), length(methods_sub)))

bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + coord_flip() +
  theme(legend.position = "none") -> sil_allen_pc

## Zeisel

load("zeisel_covariates.rda")

data <- read.table("expression_mRNA_17-Aug-2014.txt", sep='\t', stringsAsFactors = FALSE, comment.char = '%')

tissue <- as.factor(as.matrix(data)[1,-(1:2)])
table(tissue)
group <- as.factor(as.matrix(data)[2,-(1:2)])
table(tissue, group)
nmolecule <- as.numeric(as.matrix(data)[3,-(1:2)])
well <- as.factor(as.matrix(data)[4,-(1:2)])
sex <- as.factor(as.matrix(data)[5,-(1:2)])
age <- as.numeric(as.matrix(data)[6,-(1:2)])
table(tissue, sex)
table(tissue, age)
diameter <- as.numeric(as.matrix(data)[7,-(1:2)])
cell_id <- as.matrix(data)[8,-(1:2)]
batch <- as.factor(sapply(strsplit(cell_id, "_"), function(x) x[1]))
position <- as.factor(sapply(strsplit(cell_id, "_"), function(x) x[2]))

level1 <- as.factor(as.matrix(data)[9,-(1:2)])

methods <- list(pc_raw[,1:3], pc_tc[,1:3], pc_tmm[,1:3], pc_fq[,1:3],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB-WaVE")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-WaVE"))

sil_cl <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(level1), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_cl,
                   Method=factor(names(methods),
                                 levels=c("PCA_RAW", "PCA_TC", "PCA_TMM", "PCA_FQ",
                                          "ZIFA_RAW", "ZIFA_TC", "ZIFA_TMM", "ZIFA_FQ",
                                          "ZINB-WaVE")), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Method)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=colsil) + coord_flip() +
  theme(legend.position = "none") -> sil_zeisel

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(level1), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB-WaVE"), each=nlevels(level1)),
                   Cluster=rep(levels(level1), length(methods_sub)))

library(dplyr)
bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil_zeisel_pc

## Espresso

load("espresso_covariates.rda")

all.counts <- read.table("ESpresso/counttable_es.csv", header=TRUE, row.names=1, colClasses=c("character", rep("integer", 704)))
serum <- sub("ola_mES_([^_]+)_.*", "\\1", colnames(all.counts))
batch <- sub("ola_mES_[^_]+_([^_]+)_.*", "\\1", colnames(all.counts))
targets <- data.frame(Serum=serum, Batch=batch)

# Only using data from two batches.
keep <- targets$Batch %in% c("2", "3")
targets <- targets[keep,]
level1 <- targets$Serum

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W, zinb_batch@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB-WaVE", "ZINB-Batch")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-WaVE", "ZINB-WaVE"))

sil_cl <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(level1), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_cl,
                   Method=factor(names(methods),
                                 levels=c("PCA_RAW", "PCA_TC", "PCA_TMM", "PCA_FQ",
                                          "ZIFA_RAW", "ZIFA_TC", "ZIFA_TMM", "ZIFA_FQ",
                                          "ZINB-WaVE", "ZINB-Batch")), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Method)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=colsil) + coord_flip() +
  theme(legend.position = "none") -> sil_espresso

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(level1), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB-WaVE"), each=nlevels(level1)),
                   Cluster=factor(rep(levels(level1), length(methods_sub)),
                                  levels=levels(level1)))

library(dplyr)
bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1[c(2, 3, 1)]) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil_espresso_pc

## Patel

load("patel_covariates.rda")

counts <- read.table("Patel/glioblastoma_raw_rnaseq_SCandbulk_counts.txt", header=TRUE, stringsAsFactors = FALSE, row.names=NULL)

gene_symbols <- counts[,1]
ensembl_ids <- counts[,2]
sample_names <- colnames(counts)[-(1:2)]

metadata <- read.table("patel/SraRunTable.txt", sep='\t', stringsAsFactors = FALSE, header=TRUE, row.names=5, na.strings = "<not provided>")
metadata <- metadata[sample_names,]

# select only single-cell samples from patients
keep <- which(grepl("^Single cell", metadata$source_name_s) &
                !is.na(metadata$patient_id_s) &
                !is.na(metadata$subtype_s))
metadata <- metadata[keep,]
level1 <- as.factor(metadata$patient_id_s)

methods <- list(pc_raw[,1:2], pc_tc[,1:2], pc_tmm[,1:2], pc_fq[,1:2],
                zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
                zinb@W)
names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
                    paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
                    "ZINB-WaVE")
met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB-WaVE"))

sil_cl <- sapply(seq_along(methods), function(i) {
  d <- dist(methods[[i]])
  ss <- silhouette(as.numeric(level1), d)
  mean(ss[,3])
})

bars <- data.frame(AverageSilhouette=sil_cl,
                   Method=factor(names(methods),
                                 levels=c("PCA_RAW", "PCA_TC", "PCA_TMM", "PCA_FQ",
                                          "ZIFA_RAW", "ZIFA_TC", "ZIFA_TMM", "ZIFA_FQ",
                                          "ZINB-WaVE")), Type=met_type)

bars %>%
  ggplot(aes(Method, AverageSilhouette, group=Type, fill=Method)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=colsil) + coord_flip() +
  theme(legend.position = "none") -> sil_patel

methods_sub <- methods[c(2, 6, 9)]
sil_pc <- lapply(seq_along(methods_sub), function(i) {
  d <- dist(methods_sub[[i]])
  ss <- silhouette(as.numeric(level1), d)
  sss <-  summary(ss)
  sss$clus.avg.widths
})

bars <- data.frame(AverageSilhouette=unlist(sil_pc),
                   Method=rep(c("PCA", "ZIFA", "ZINB-WaVE"), each=nlevels(level1)),
                   Cluster=rep(levels(level1), length(methods_sub)))

library(dplyr)
bars %>%
  dplyr::mutate(ClusterByMethod = paste0(Cluster, " ", Method)) %>%
  ggplot(aes(ClusterByMethod, AverageSilhouette, fill=Cluster)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col1) + coord_flip() +
  theme(legend.position = "none", axis.text = element_text(size=8)) -> sil_patel_pc

# ## Olfactory
# load("olfactory.rda")
#
# load("Expt4c2b_filtdata.Rda")
# load("E4c2b_slingshot_wsforkelly.RData")
#
# methods <- list(pc_raw[,1:3], pc_tc[,1:3], pc_tmm[,1:3], pc_fq[,1:3],
#                 zifa_raw, zifa_tc, zifa_tmm, zifa_fq,
#                 zinb_Vall@W)
# names(methods) <- c(paste0("PCA_", c("RAW", "TC", "TMM", "FQ")),
#                     paste0("ZIFA_", c("RAW", "TC", "TMM", "FQ")),
#                     "ZINB")
# met_type <- as.factor(c(rep(c("PCA", "ZIFA"), each=4), "ZINB"))
#
# sil_cl <- sapply(seq_along(methods), function(i) {
#   d <- dist(methods[[i]])
#   ss <- silhouette(as.numeric(clus.labels), d)
#   mean(ss[,3])
# })
#
# bars <- data.frame(AverageSilhouette=sil_cl, Method=names(methods), Type=met_type)
#
# bars %>%
#   ggplot(aes(Method, AverageSilhouette, group=Type, fill=Type)) +
#   geom_bar(stat="identity", position='dodge') +
#   scale_fill_manual(values=col1) + coord_flip() +
#   theme(legend.position = "none") -> sil_olfactory

sil <- plot_grid(sil_allen, sil_zeisel,
                 sil_patel, sil_espresso,
                 labels=c("a", "b", "c", "d"))

save_plot("silhouette.pdf",
          sil,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3
)

sil_pc <- plot_grid(sil_allen_pc, sil_zeisel_pc,
                 sil_patel_pc, sil_espresso_pc,
                 labels=c("a", "b", "c", "d"))

save_plot("silhouette_pc.pdf",
          sil_pc,
          ncol = 2,
          nrow = 2,
          base_aspect_ratio = 1.3
)
