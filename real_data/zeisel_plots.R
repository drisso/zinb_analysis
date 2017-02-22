library(cowplot)
library(RColorBrewer)
library(magrittr)
library(ggplot2)

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
level2 <- as.factor(as.matrix(data)[10,-(1:2)])

counts <- as.matrix(data[12:NROW(data),-(1:2)])
counts <- matrix(as.numeric(counts), ncol=ncol(counts), nrow=nrow(counts))
rownames(counts) <- data[12:NROW(data),1]
colnames(counts) <- data[8, -(1:2)]
filter <- rowSums(counts>10)>=10
raw <- counts[filter,]

detection_rate <- colSums(raw>0)
coverage <- colSums(raw)

colPal1 <- colorRampPalette(c("blue", "yellow"))(nlevels(batch))
colBatch <- colPal1[batch]
col1 <- brewer.pal(9, "Set1")
col2 <- brewer.pal(8, "Set2")
colMerged <- col2[level1]

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set2") -> panel1_pca

data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() +
  scale_color_brewer(palette="Set2") -> panel1_zifa

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point()  +
  scale_color_brewer(palette="Set2")  -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
          panel1_zifa + theme(legend.position = "none"),
          panel1_zinb + theme(legend.position = "none"),
          labels=c("A", "B", "C"), align = "h", ncol=3)

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
                labels=c("D", "E", "F"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, .6))

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("zeisel_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)

qc <- data.frame(age, coverage, detection_rate, diameter, nmolecule, well)
qc <- apply(qc, 2, as.numeric)

cors <- lapply(1:3, function(i) abs(cor(pc_tc[,i], qc, method="spearman")))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_pca

cors <- lapply(1:3, function(i) abs(cor(zifa_tc[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zifa

cors <- lapply(1:3, function(i) abs(cor(zinb@W[,i], qc)))
cors <- unlist(cors)
bars <- data.frame(AbsoluteCorrelation=cors,
                   QC=rep(stringr::str_to_lower(colnames(qc)), 3),
                   Dimension=as.factor(rep(1:3, each=ncol(qc))))

bars %>%
  ggplot(aes(Dimension, AbsoluteCorrelation, group=QC, fill=QC)) +
  geom_bar(stat="identity", position='dodge') +
  scale_fill_manual(values=col2) + ylim(0, 1) -> panel2_zinb

p2 <- plot_grid(panel2_pca + theme(legend.position = "none"),
                panel2_zifa + theme(legend.position = "none"),
                panel2_zinb + theme(legend.position = "none"),
                labels=c("D", "E", "F"), align = "h", ncol=3)

legend2 <- get_legend(panel2_pca)
lower <- plot_grid(p2, legend2, rel_widths = c(3, 1))

fig1bis <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("zeisel_fig1bis.pdf", fig1bis,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)


library(cluster)

methods <- list(pc_raw[,1:3], pc_tc[,1:3], pc_tmm[,1:3], pc_fq[,1:3],
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

save_plot("zeisel_fig1tris.pdf", fig1_tris,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)
