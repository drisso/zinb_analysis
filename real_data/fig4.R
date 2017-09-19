## color palette
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"))

## t-sne
load("tenx_pbmc/kmeans_res.rda")
load("tenx_pbmc/zinb_res.rda")
load("tenx_pbmc/zinb_10.rda")
load("tenx_pbmc/seq_cluster.rda")

cl_anno <- as.character(z_cl)
cl_anno[z_cl == 0] <- "CD4+ T-cells (1)"
cl_anno[z_cl == 1] <- "CD4+ T-cells (2)"
cl_anno[z_cl == 2] <- "CD8+ T-cells (1)"
cl_anno[z_cl == 3] <- "CD8+ T-cells (2)"
cl_anno[z_cl == 4] <- "CD8+ T-cells (3)"
cl_anno[z_cl == 6] <- "CD8+ T-cells (4)"
cl_anno[z_cl == 5] <- "Activated T-cells"
cl_anno[z_cl == 7] <- "NK cells"
cl_anno[z_cl == 8] <- "B-cells"
cl_anno[z_cl == 9] <- "CD14+ monocytes"
cl_anno[z_cl == 10] <- "CD16+ monocytes"
cl_anno[z_cl == 11] <- "CD1C+ DCs"
cl_anno[z_cl == 13] <- "pDCs"
cl_anno[z_cl == 15] <- "Megakaryocytes"
cl_anno[z_cl == 12] <- "Uncharacterized (1)"
cl_anno[z_cl == 14] <- "Uncharacterized (2)"
cl_anno[z_cl == 16] <- "Uncharacterized (3)"
cl_anno[z_cl == 17] <- "Uncharacterized (4)"
cl_anno[z_cl == 18] <- "Uncharacterized (5)"

library(magrittr)
library(cowplot)

z_cl2 <- factor(paste0("Cluster", z_cl), levels=paste0("Cluster", 0:18))
pca_cl2 <- factor(paste0("Cluster", pca_cl), levels=paste0("Cluster", 0:9))

df1 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                  W1=getW(zinb_res)[,1], W2=getW(zinb_res)[,2],
                  PC1=pc_tc[,1], PC2=pc_tc[,2],
                  Cluster=z_cl2)
df1 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb

df1 %>%
  ggplot(aes(x = W1, y = W2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_W

df1 %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_pca

df2 <- data.frame(Dim1=tsne_zinb10$Y[,1], Dim2=tsne_zinb10$Y[,2],
                  W1=getW(zinb_res)[,1], W2=getW(zinb_res)[,2],
                  PC1=pc_tc[,1], PC2=pc_tc[,2],
                  Cluster=cl_anno)

df2 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb_anno

df2 %>%
  ggplot(aes(x = W1, y = W2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_W_anno

df2 %>%
  ggplot(aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_pca_anno

sup_fig <- plot_grid(fig_W + theme(legend.position = "none"),
          fig_pca, nrow = 1, ncol=2, rel_widths = c(1, 1.5), labels = "auto")
save_plot("Sup_pbmc.png", sup_fig, base_aspect_ratio = 1.5, base_height = 8)

df3 <- data.frame(Dim1=tsne_zinb10$Y[cs_cl != "-1",1], Dim2=tsne_zinb10$Y[cs_cl != "-1",2],
                  Cluster=droplevels(cs_cl[cs_cl != "-1"]))

df3 %>%
  ggplot(aes(x = Dim1, y = Dim2, color = Cluster)) +
  geom_point(size = 0.5) + scale_color_manual(values=pal) -> fig_tsne_zinb_seq
save_plot("Sup_pbmc_seq.png", fig_tsne_zinb_seq, base_aspect_ratio = 1.5, base_height = 8)
