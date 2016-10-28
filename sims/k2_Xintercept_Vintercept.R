# Use the allen data to simulate with realistic parameters
library(scRNAseq)
library(EDASeq)
library(zinb)
library(matrixStats)

set.seed(1829)
data("allen")

# select only core cells and endogenous genes
allen_core <- allen[grep("^ERCC-", rownames(allen), invert = TRUE),
                    which(colData(allen)$Core.Type=="Core")]

# filter out lowly expressed genes
filter <- apply(assay(allen_core)>10, 1, sum)>=10

# randomly select 100 samples
sampleidx <- sample(colnames(allen_core), 100)

allen_core <- allen_core[filter, sampleidx]
core <- assay(allen_core)
bio <- as.factor(colData(allen_core)$driver_1_s)

# select only 1000 most variable genes
vars <- rowVars(log1p(core))
names(vars) <- rownames(core)
vars <- sort(vars, decreasing = TRUE)
core <- core[names(vars)[1:1000],]

zinb_true <- zinbFit(core, ncores = 3, K = 2)

true_W <- zinb_true@W
plot(true_W, pch=19, col=bio)

mod <- model.matrix(~true_W)

zinb_pars <- zinbFit(core, ncores = 3, X=mod, which_X_mu=1:3, which_X_pi=1:3, commondispersion=TRUE)

true_alpha_mu <- zinb_pars@beta_mu[-1,]
true_alpha_pi <- zinb_pars@beta_pi[-1,]
true_beta_mu <- zinb_pars@beta_mu[1,,drop=FALSE]
true_beta_pi <- zinb_pars@beta_pi[1,,drop=FALSE]
true_gamma_mu <- zinb_pars@gamma_mu
true_gamma_pi <- zinb_pars@gamma_pi

sim_obj <- zinbModel(W=true_W, alpha_mu=true_alpha_mu, alpha_pi=true_alpha_pi, beta_mu=true_beta_mu, beta_pi=true_beta_pi, zeta = zinb_pars@zeta, gamma_mu=true_gamma_mu, gamma_pi=true_gamma_pi)

B <- 50
sim_data <- lapply(seq_len(B), function(i) zinbSim(sim_obj, seed=i))
save(sim_obj, sim_data, file="k2_Xintercept_Vintercept.rda")

# tmp <- zinbSim(sim_obj, seed=i)
# fit <- zinbFit(t(tmp$counts), ncores=3, K=2)
# plot(fit@W, pch=19, col=bio)
# pca <- prcomp(log1p(tmp$counts), center=TRUE, scale=TRUE)
# plot(pca$x, pch=19, col=bio)
#
# dtrue <- as.matrix(dist(true_W))
# dest <- as.matrix(dist(fit@W))
# dpca <- as.matrix(dist(pca$x[,1:2]))
# hist(sapply(1:100, function(i) cor(dtrue[,i], dest[,i])), main="ZINB")
# hist(sapply(1:100, function(i) cor(dtrue[,i], dpca[,i])), main="PCA")
#
# library(cluster)
# pamtrue <- pam(dtrue, k=3)
# strue <- silhouette(pamtrue$clustering, dtrue)
# plot(strue, col=1:3)
# sest <- silhouette(pamtrue$clustering, dest)
# plot(sest, col=1:3)
# spca <- silhouette(pamtrue$clustering, dpca)
# plot(spca, col=1:3)
#
# tapply(strue[,3], strue[,1], mean)
# tapply(sest[,3], sest[,1], mean)
# tapply(spca[,3], spca[,1], mean)
