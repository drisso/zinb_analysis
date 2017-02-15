library(cowplot)
library(scRNAseq)
load("allen_covariates.rda")

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

load("allen_covariates_1000.rda")

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
upper <- plot_grid(p1, legend, rel_widths = c(3, 1))

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

plot_grid(upper, lower, ncol=1, nrow=2)

## Alternatively (too crowded -- move upper to supplementary?)

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster)) + geom_point() + scale_color_manual(values=col2) -> panel1_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster)) + geom_point() + scale_color_manual(values=col2) -> panel1_zifa
data.frame(Dim1=zinb@W[,1], Dim2=zinb@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, shape=layer, color=cluster)) + geom_point() + scale_color_manual(values=col2) -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("A", "B", "C"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, 1))

plot_grid(upper, lower, ncol=1, nrow=2)

### Silhouette

library(cluster)
d <- dist(zinb@W)
ss <- silhouette(as.numeric(cluster), d)
plot(ss, border=col2[sort(as.numeric(cluster), decreasing = TRUE)], main="ZINB")
pamres <- pam(d, k=4)
s <- silhouette(pamres, d)
mean(s[,3])
tapply(s[,3], s[,1], mean)
plot(s, border=collayer[as.numeric(rownames(s))], main="ZINB")

d <- dist(pc_tc[,1:2])
ss <- silhouette(as.numeric(cluster), d)
plot(ss, border=col2[sort(as.numeric(cluster), decreasing = TRUE)], main="PCA")
mean(ss[,3])
hist(ss[,3])
pamres <- pam(d, k=4)
s <- silhouette(pamres, d)
mean(s[,3])
tapply(s[,3], s[,1], mean)
plot(s, border=collayer[as.numeric(rownames(s))], main="PCA")

d <- dist(zifa_tc)
ss <- silhouette(as.numeric(cluster), d)
plot(ss, border=col2[sort(as.numeric(cluster), decreasing = TRUE)], main="ZIFA")
mean(ss[,3])
hist(ss[,3])
pamres <- pam(d, k=4)
s <- silhouette(pamres, d)
mean(s[,3])
tapply(s[,3], s[,1], mean)
plot(s, border=collayer[as.numeric(rownames(s))], main="ZIFA")
