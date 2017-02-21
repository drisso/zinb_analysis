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
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() -> panel1_pca

data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point() -> panel1_zifa

data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=level1)) + geom_point()  -> panel1_zinb

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

library(cluster)
# d <- dist(zinb@W)
# s <- silhouette(as.numeric(level1), d)
# mean(s[,3])
# tapply(s[,3], s[,1], mean)
# plot(s, border=col2[sort(as.numeric(level1), decreasing = TRUE)], main="ZINB")
#
# d <- dist(pc_tc[,1:2])
# s <- silhouette(as.numeric(level1), d)
# mean(s[,3])
# tapply(s[,3], s[,1], mean)
# plot(s, border=col2[sort(as.numeric(level1), decreasing = TRUE)], main="PCA")
#
# d <- dist(zifa_tc)
# s <- silhouette(as.numeric(level1), d)
# mean(s[,3])
# tapply(s[,3], s[,1], mean)
# plot(s, border=col2[sort(as.numeric(level1), decreasing = TRUE)], main="ZIFA")

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

pdf("zeisel_sil.pdf")
barplot(sil_cl, col=col1[met_type], horiz = TRUE, names.arg = names(methods), las=2, main="Average silhouette width")
dev.off()

cors <- sapply(seq_along(methods), function(i) {
  max(abs(cor(methods[[i]][,1], detection_rate)), abs(cor(methods[[i]][,2], detection_rate)))
})

pdf("zeisel_sil.pdf")
barplot(cors, col=col1[met_type], horiz = TRUE, names.arg = names(methods), las=2, main="Absolute correlation with Detection Rate")
dev.off()
