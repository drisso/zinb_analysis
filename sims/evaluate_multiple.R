setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims')
library(cluster)
library(zinb)
library(EDASeq)

# input data
load("k2_Xintercept_Vintercept.rda")
load('k2_Xintercept_Vintercept_fitted.rda')
# fittedSim is a list with order k (in 1:4), Xintercept (T or F),
# Vintercept (T or F) and commondispersion (T or F), simulated dataset.
# For example fittedSim[[3]][[1]][[1]][[2]][[10]] was fitted with k = 3,
# Xintercept = T, Vintercept = T, and commondispersion = F
# for simulated dataset 10. See code in run_sim.R

# useful functions
eval_cor <- function(dtrue, dest) {
  corr <- sapply(seq_len(NCOL(dtrue)), function(i) cor(dtrue[,i], dest[,i]))
  return(corr)
}

eval_diff <- function(dtrue, dest) {
  diff <- sapply(seq_len(NCOL(dtrue)), function(i) dest[,i] - dtrue[,i])
  return(as.vector(diff))
}

eval_sil <- function(labels, dtrue, dest) {
  sest <- silhouette(labels, dest)
  return(sest)
}

# useful vars
true_W <- sim_obj@W
dtrue <- as.matrix(dist(true_W))
pamtrue <- pam(dtrue, k=3)$clustering
strue <- silhouette(pamtrue, dtrue)

fit = lapply(1:4, function(i) fittedSim[[i]][[1]][[1]][[1]][[1]])
dest = lapply(1:4, function(i) as.matrix(dist(fit[[i]]@W)))
corr = lapply(1:4, function(i) eval_cor(dtrue, dest[[i]]))
sil = lapply(1:4, function(i) eval_sil(pamtrue, dtrue, dest[[i]]))

par(mfrow=c(2,2))
lapply(1:4, function(i){
  hist(corr[[i]], breaks=20, xlim=c(0, 1), main=paste("K = ",i),
       xlab=paste0("Mean correlation: ", round(mean(corr[[i]]), 2)))
})

lapply(1:4, function(i){
  plot(sil[[i]], main=paste("K = ",i), col=1:3)
})

par(mfrow = c(1,1))
plot(true_W, pch=19, col=pamtrue, main="True W", xlab="W1", ylab="W2")

par(mfrow = c(1,2))
plot(fit[[1]]@W, rep(0, NROW(true_W)), pch=19, col=pamtrue, main="K=1")
plot(fit[[2]]@W, pch=19, col=pamtrue, main="K=2")
pairs(fit[[3]]@W, pch=19, col=pamtrue, main="K=3")
pairs(fit[[4]]@W, pch=19, col=pamtrue, main="K=4")


## comments: there is no guarantee on the order of the W's
eval_data <- function(fittedSim, k = 1:4, Xintercept = 1, Vintercept = 1,
                      commondisp = 1, simData = 1) {
  
  fit = lapply(1:4, function(i){
    fittedSim[[i]][[Xintercept]][[Vintercept]][[commondisp]][[simData]]
  })
  
  fq <- betweenLaneNormalization(t(sim_data[[simData]]$counts), which="full")
  pca <- prcomp(t(log1p(fq)))
  
  dest = lapply(1:4, function(i) as.matrix(dist(fit[[i]]@W)))
  diff = lapply(1:4, function(i) eval_diff(dtrue, dest[[i]]))
  corr = lapply(1:4, function(i) eval_cor(dtrue, dest[[i]]))
  sil = lapply(1:4, function(i) eval_sil(pamtrue, dtrue, dest[[i]]))

  dest[[5]] <- as.matrix(dist(pca$x[,1:2]))
  diff[[5]] <- eval_diff(dtrue, dest[[5]])
  corr[[5]] <- eval_cor(dtrue, dest[[5]])
  sil[[5]] <- eval_sil(pamtrue, dtrue, dest[[5]])
  
  retval <- list(diff[[1]], diff[[2]], diff[[3]], diff[[4]], diff[[5]],
                 corr[[1]], corr[[2]], corr[[3]], corr[[4]], corr[[5]],
                 sil[[1]][,3], sil[[2]][,3], sil[[3]][,3], sil[[4]][,3],
                 sil[[5]][,3])
  return(retval)
}

res <- mclapply(1:50, function(i){
  eval_data(fittedSim, simData = i)
})

diff <- sapply(1:4, function(i) rowMeans(sapply(res, function(x) x[[i]])))
boxplot(diff, main="Difference between true and estimated sample distances in W", border=c(1, 2, 1, 1), xlab="k", ylab="Estimated distance - true distance")
abline(h=0, lty=2)

diff <- sapply(1:5, function(i) rowMeans(sapply(res, function(x) x[[i]])))
colnames(diff) <- c(paste0("zinb k=", 1:4), "FQ PCA (k=2)")
boxplot(diff, main="Difference between true and estimated sample distances in W", border=c(1, 2, 1, 1), xlab="method", ylab="Estimated distance - true distance")
abline(h=0, lty=2)

corr <- sapply(6:9, function(i) rowMeans(sapply(res, function(x) x[[i]])))
boxplot(corr, main="Correlation between true and estimated sample distances in W", border=c(1, 2, 1, 1), xlab="k", ylab="Correlation")

corr <- sapply(6:10, function(i) rowMeans(sapply(res, function(x) x[[i]])))
colnames(corr) <- c(paste0("zinb k=", 1:4), "FQ PCA (k=2)")
boxplot(corr, main="Correlation between true and estimated sample distances in W", border=c(1, 2, 1, 1), xlab="method", ylab="Correlation")

sil <- sapply(11:14, function(i) rowMeans(sapply(res, function(x) x[[i]])))
boxplot(sil, main="Silhouette width of the PAM clustering based on true W", border=c(1, 2, 1, 1),  xlab="k", ylab="Silhouette width")
boxplot(apply(sil, 2, function(x) x - strue[,3]), main="Silhouette width of the PAM clustering based on true W", border=c(1, 2, 1, 1),  xlab="k", ylab="Silhouette in estimated W - Silhouette in true W")
abline(h=0, lty=2)

sil <- sapply(11:15, function(i) rowMeans(sapply(res, function(x) x[[i]])))
colnames(sil) <- c(paste0("zinb k=", 1:4), "FQ PCA (k=2)")
boxplot(sil, main="Silhouette width of the PAM clustering based on true W", border=c(1, 2, 1, 1),  xlab="method", ylab="Silhouette width")
boxplot(apply(sil, 2, function(x) x - strue[,3]), main="Silhouette width of the PAM clustering based on true W", border=c(1, 2, 1, 1),  xlab="k", ylab="Silhouette in estimated W - Silhouette in true W")
abline(h=0, lty=2)
dev.off()
