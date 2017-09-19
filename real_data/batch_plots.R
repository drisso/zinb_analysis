rm(list=ls())
library(zinbwave)
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(9, "Set3"))

analysis <- FALSE

if(analysis) {
  ## COMBAT
  library(sva)
  library(rARPACK)
  library(Rtsne)
  combat_res <- ComBat(log(as.matrix(tc) + .1), batch)
  pc_combat  <- fastpca(combat_res[rownames(dense),])
  set.seed(22002)
  tsne_combat <- Rtsne(pc_combat[,1:50], pca=FALSE)

  save(pc_combat, tsne_combat, file="pbmc_combat.rda")
}

zinb_df <- read.table("zinb_tenx_combined_clusters.txt")
pca_df <- read.table("pca_tenx_combined_clusters.txt")
load("zinb_tenx_combined.rda")
load("pbmc_combined.rda")
load("pbmc_combat.rda")

plot(tsne_zinb$Y, pch=19, col=pal[batch])
plot(tsne_data$Y, pch=19, col=pal[batch])
plot(tsne_combat$Y, pch=19, col=pal[batch])

plot(tsne_combat$Y, pch=19, col=pal[zinb_df$paste0..Cluster...z_cl.])
plot(tsne_zinb$Y, pch=19, col=pal[zinb_df$paste0..Cluster...z_cl.])
plot(tsne_data$Y, pch=19, col=pal[zinb_df$paste0..Cluster...z_cl.])

plot(tsne_combat$Y, pch=19, col=pal[pca_df$paste0..PCA...pca_cl.])
plot(tsne_zinb$Y, pch=19, col=pal[pca_df$paste0..PCA...pca_cl.])
plot(tsne_data$Y, pch=19, col=pal[pca_df$paste0..PCA...pca_cl.])

## predict batch from data and look at prediction error
library(MASS)
fit_zinb <- lda(x = zinb_df[,1:10], grouping = batch)
pred_zinb = predict(fit_zinb)
sum(pred_zinb$class != batch) / length(batch)

fit_pca <- lda(x = pca_df[,1:10], grouping = batch)
pred_pca = predict(fit_pca)
sum(pred_pca$class != batch) / length(batch)

fit_combat <- lda(x = pc_combat[,1:10], grouping = batch)
pred_combat = predict(fit_combat)
sum(pred_combat$class != batch) / length(batch)

## use detection rate as a covariate instead of the batch

