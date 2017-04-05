library(cowplot)
library(scRNAseq)
library(RColorBrewer)
library(magrittr)
library(ggplot2)
load("olfactory.rda")

load("Expt4c2b_filtdata.Rda")
load("E4c2b_slingshot_wsforkelly.RData")

colMerged <- cc_rev[clus.labels]

qc <- qc[names(clus.labels),]
col2 <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Set3"), brewer.pal(8, "Set1"))

data.frame(Dim1=pc_tc[,1], Dim2=pc_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=clus.labels)) + geom_point() +
  scale_color_manual(values=cc_rev)-> panel1_pca
data.frame(Dim1=zifa_tc[,1], Dim2=zifa_tc[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=clus.labels)) + geom_point() +
  scale_color_manual(values=cc_rev)  -> panel1_zifa
data.frame(Dim1=zinb_Vall@W[,1], Dim2=zinb_Vall@W[,2]) %>%
  ggplot(aes(Dim1, Dim2, colour=clus.labels)) + geom_point() +
  scale_color_manual(values=cc_rev) -> panel1_zinb

p1 <- plot_grid(panel1_pca + theme(legend.position = "none"),
                panel1_zifa + theme(legend.position = "none"),
                panel1_zinb + theme(legend.position = "none"),
                labels=c("A", "B", "C"), align = "h", ncol=3)

legend <- get_legend(panel1_pca)
upper <- plot_grid(p1, legend, rel_widths = c(3, .6))

cors <- lapply(1:3, function(i) abs(cor(pc_tc[,i], qc)))
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

cors <- lapply(1:3, function(i) abs(cor(zinb_Vall@W[,i], qc)))
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

fig1 <- plot_grid(upper, lower, ncol=1, nrow=2)

save_plot("olfactory_fig1.pdf", fig1,
          ncol = 3,
          nrow = 3,
          base_aspect_ratio = 1.3
)
