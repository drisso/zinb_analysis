library(cluster)
library(zinb)
library(EDASeq)

# input data
load("k2_Xintercept_Vintercept.rda")
ncores <- 2
outprefix <- "k2_Xintercept_Vintercept_out"

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

fit1 <- zinbFit(t(sim_data[[1]]$counts), ncores=ncores, K=1)
fit2 <- zinbFit(t(sim_data[[1]]$counts), ncores=ncores, K=2)
fit3 <- zinbFit(t(sim_data[[1]]$counts), ncores=ncores, K=3)
fit4 <- zinbFit(t(sim_data[[1]]$counts), ncores=ncores, K=4)

dest1 <- as.matrix(dist(fit1@W))
dest2 <- as.matrix(dist(fit2@W))
dest3 <- as.matrix(dist(fit3@W))
dest4 <- as.matrix(dist(fit4@W))

corr1 <- eval_cor(dtrue, dest1)
corr2 <- eval_cor(dtrue, dest2)
corr3 <- eval_cor(dtrue, dest3)
corr4 <- eval_cor(dtrue, dest4)

sil1 <- eval_sil(pamtrue, dtrue, dest1)
sil2 <- eval_sil(pamtrue, dtrue, dest2)
sil3 <- eval_sil(pamtrue, dtrue, dest3)
sil4 <- eval_sil(pamtrue, dtrue, dest4)

pdf(paste0(outprefix, "_figure_one_sim.pdf"))
hist(corr1, breaks=20, xlim=c(0, 1), main="K=1", xlab=paste0("Mean correlation: ", round(mean(corr1), 2)))
hist(corr2, breaks=20, xlim=c(0, 1), main="K=2", xlab=paste0("Mean correlation: ", round(mean(corr2), 2)))
hist(corr3, breaks=20, xlim=c(0, 1), main="K=3", xlab=paste0("Mean correlation: ", round(mean(corr3), 2)))
hist(corr4, breaks=20, xlim=c(0, 1), main="K=4", xlab=paste0("Mean correlation: ", round(mean(corr4), 2)))
plot(sil1, main="K=1", col=1:3)
plot(sil2, main="K=2", col=1:3)
plot(sil3, main="K=3", col=1:3)
plot(sil4, main="K=4", col=1:3)

plot(true_W, pch=19, col=pamtrue, main="True W", xlab="W1", ylab="W2")
plot(fit1@W, rep(0, NROW(true_W)), pch=19, col=pamtrue, main="K=1")
plot(fit2@W, pch=19, col=pamtrue, main="K=2")
pairs(fit3@W, pch=19, col=pamtrue, main="K=3")
pairs(fit4@W, pch=19, col=pamtrue, main="K=4")
dev.off()

## comments: there is no guarantee on the order of the W's
eval_data <- function(x) {
  fit2 <- zinbFit(t(x$counts), ncores=1, K=2)
  fit1 <- zinbFit(t(x$counts), ncores=1, K=1)
  fit3 <- zinbFit(t(x$counts), ncores=1, K=3)
  fit4 <- zinbFit(t(x$counts), ncores=1, K=4)

  fq <- betweenLaneNormalization(t(x$counts), which="full")
  pca <- prcomp(t(log1p(fq)))

  dest1 <- as.matrix(dist(fit1@W))
  dest2 <- as.matrix(dist(fit2@W))
  dest3 <- as.matrix(dist(fit3@W))
  dest4 <- as.matrix(dist(fit4@W))
  destPC <- as.matrix(dist(pca$x[,1:2]))

  diff1 <- eval_diff(dtrue, dest1)
  diff2 <- eval_diff(dtrue, dest2)
  diff3 <- eval_diff(dtrue, dest3)
  diff4 <- eval_diff(dtrue, dest4)
  diffPC <- eval_diff(dtrue, destPC)

  corr1 <- eval_cor(dtrue, dest1)
  corr2 <- eval_cor(dtrue, dest2)
  corr3 <- eval_cor(dtrue, dest3)
  corr4 <- eval_cor(dtrue, dest4)
  corrPC <- eval_cor(dtrue, destPC)

  sil1 <- eval_sil(pamtrue, dtrue, dest1)
  sil2 <- eval_sil(pamtrue, dtrue, dest2)
  sil3 <- eval_sil(pamtrue, dtrue, dest3)
  sil4 <- eval_sil(pamtrue, dtrue, dest4)
  silPC <- eval_sil(pamtrue, dtrue, destPC)

  retval <- list(diff1, diff2, diff3, diff4, diffPC,
                 corr1, corr2, corr3, corr4, corrPC,
                 sil1[,3], sil2[,3], sil3[,3], sil4[,3], silPC[,3])
  return(retval)
}

res <- mclapply(sim_data, eval_data, mc.cores = ncores)

pdf(paste0(outprefix, "_figure_average.pdf"))
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
