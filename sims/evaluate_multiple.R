# setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims/datasets/')
# Comment: I commented out the line above to make it work in all local copies
# The setwd() should be taken care of by the .Rproj file

library(cluster)
library(zinb)
library(EDASeq)

# input data
load("sim_allen5_fitted.rda")
load("sim_allen5.rda")
# fittedSim is a list with order k (in 1:4), Xintercept (T or F),
# Vintercept (T or F) and commondispersion (T or F), simulated dataset.
# For example fittedSim[[3]][[1]][[2]][[10]] was fitted with k = 3,
# Vintercept = T, and commondispersion = F for simulated dataset 10.
# See code in run_sim.R

## Let's look at bias and MSE

## levels: k (4), Vint (2), commondisp (2), nsim (50)
true_W <- simModel@W
K <- 1:4
Vintercept <- commondispersion <- 1:2


###########
# BIAS
###########
biasList <- lapply(K, function(k){
  tmp <- lapply(Vintercept, function(Vint){
    tmp <- lapply(commondispersion, function(commondisp){
      gamma_mu <- rowMeans(sapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_mu[1,]
      }))
      
      beta_mu <- rowMeans(sapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_mu[1,]
      }))
      
      gamma_pi <- rowMeans(sapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_pi[1,]
      }))
      
      beta_pi <- rowMeans(sapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_pi[1,]
      }))
      
      theta <- rowMeans(sapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@zeta
      }))
      
      tmp <- lapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_mu
      })
      walpha_mu <- Reduce("+", tmp)/length(tmp)
      
      tmp <- lapply(seq_along(simData), function(i) {
        fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_pi
      })
      walpha_pi <- Reduce("+", tmp)/length(tmp)
      
      mu <- simModel@X %*% beta_mu + t(simModel@V %*% gamma_mu) + walpha_mu
      pi <- simModel@X %*% beta_pi + t(simModel@V %*% gamma_pi) + walpha_pi
      
      return(list(gamma_mu=gamma_mu - simModel@gamma_mu[1,],
                  gamma_pi=gamma_pi - simModel@gamma_pi[1,],
                  beta_mu=beta_mu - simModel@beta_mu[1,],
                  beta_pi=beta_pi - simModel@beta_pi[1,],
                  theta=theta - simModel@zeta,
                  walpha_mu=as.vector(walpha_mu - simModel@W %*% simModel@alpha_mu),
                  walpha_pi=as.vector(walpha_pi) - simModel@W %*% simModel@alpha_pi,
                  mu=as.vector(mu - getLogMu(simModel)),
                  pi=as.vector(pi - getLogitPi(simModel))))
    })
    names(tmp) <- c("common_disp", "genewise_disp")
    return(tmp)
  })
  names(tmp) <- c("V", "noV")
  return(tmp)
})
names(biasList) <- paste0("K", K)
plotbias <- unlist(unlist(unlist(biasList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
bias = lapply(seq_along(plotbias), function(i){
  nn = strsplit(names(plotbias)[i], '\\.')[[1]]
  data.frame(bias = as.vector(plotbias[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
})
bias = data.frame(do.call(rbind, bias), stringsAsFactors = F)
save(bias, file = "sim_allen5_bias.rda")


###########
# VARIANCE
###########
n = length(simData)
varianceList <- lapply(K, function(k){
  tmp <- lapply(Vintercept, function(Vint){
    tmp <- lapply(commondispersion, function(commondisp){
      gamma_mu <- (1/n-1) * rowSums(sapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_mu[1,] - simModel@gamma_mu[1,])^2
      }))
      
      beta_mu <- (1/n-1) * rowSums(sapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_mu[1,] - simModel@beta_mu[1,])^2
      }))
      
      gamma_pi <- (1/n-1) * rowSums(sapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_pi[1,] - simModel@gamma_pi[1,])^2
      }))
      
      beta_pi <- (1/n-1) * rowSums(sapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_pi[1,] - simModel@beta_pi[1,])^2
      }))
      
      theta <- (1/n-1) * rowSums(sapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@zeta - simModel@zeta)^2
      }))
      
      tmp <- lapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_mu -
           simModel@W %*% simModel@alpha_mu)^2
      })
      walpha_mu <- Reduce("+", tmp)/(n-1)
      
      tmp <- lapply(seq_along(simData), function(i) {
        (fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_pi -
           simModel@W %*% simModel@alpha_pi)^2
      })
      walpha_pi <- Reduce("+", tmp)/(n-1)
      
      tmp <- lapply(seq_along(simData), function(i) {
        (simModel@X %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_mu +
           t(simModel@V %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_mu) +
           fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_mu -
           getLogMu(simModel))^2
      })
      mu <- Reduce("+", tmp)/(n-1)
      
      tmp <- lapply(seq_along(simData), function(i) {
        (simModel@X %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_pi +
           t(simModel@V %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_pi) +
           fittedSim[[k]][[Vint]][[commondisp]][[i]]@W %*% fittedSim[[k]][[Vint]][[commondisp]][[i]]@alpha_pi -
           getPi(simModel))^2
      })
      pi <- Reduce("+", tmp)/(n-1)
      
       
      return(list(gamma_mu=gamma_mu,
                  gamma_pi=gamma_pi,
                  beta_mu=beta_mu,
                  beta_pi=beta_pi,
                  theta=theta,
                  walpha_mu=as.vector(walpha_mu),
                  walpha_pi=as.vector(walpha_pi),
                  mu=as.vector(mu),
                  pi=as.vector(pi)))
    })
    names(tmp) <- c("common_disp", "genewise_disp")
    return(tmp)
  })
  names(tmp) <- c("V", "noV")
  return(tmp)
})
names(varianceList) <- paste0("K", K)
plotvariance <- unlist(unlist(unlist(varianceList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
variance = lapply(seq_along(plotvariance), function(i){
  nn = strsplit(names(plotvariance)[i], '\\.')[[1]]
  data.frame(variance = as.vector(plotvariance[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
})
variance = data.frame(do.call(rbind, variance), stringsAsFactors = F)
save(variance, file = "sim_allen5_variance.rda")


##########
# FIT
##########
fit <- lapply(1:4, function(i) fittedSim[[i]][[1]][[1]][[1]])
save(fit, file="sim_allen5_first.rda")

############################
# dim. red. and Silhouette
############################
# useful functions
eval_cor <- function(dtrue, dest) {
  corr <- sapply(seq_len(NCOL(dtrue)), function(i) cor(dtrue[,i], dest[,i]))
  return(corr)
}

eval_diff <- function(dtrue, dest) {
  diff <- sapply(seq_len(NCOL(dtrue)), function(i) dest[,i] - dtrue[,i])
  return(as.vector(diff))
}

eval_sil <- function(labels, dest) {
  sest <- silhouette(labels, dest)
  return(sest)
}

# useful vars
true_W <- simModel@W
dtrue <- as.matrix(dist(true_W))
biotrue <- as.numeric(factor(bio))
strue <- silhouette(biotrue, dtrue)

fit = lapply(1:4, function(i) fittedSim[[i]][[1]][[1]][[1]])
dest = lapply(1:4, function(i) as.matrix(dist(fit[[i]]@W)))
corr = lapply(1:4, function(i) eval_cor(dtrue, dest[[i]]))
sil = lapply(1:4, function(i) eval_sil(biotrue, dest[[i]]))

par(mfrow=c(2,2))
lapply(1:4, function(i){
  hist(corr[[i]], breaks=20, xlim=c(0, 1), main=paste("K = ",i),
       xlab=paste0("Mean correlation: ", round(mean(corr[[i]]), 2)))
})

lapply(1:4, function(i){
  plot(sil[[i]], main=paste("K = ",i), col=1:3)
})


## comments: there is no guarantee on the order of the W's
eval_data <- function(fittedSim, simData, k = 1:4, Vintercept = 1,
                      commondisp = 2, sim = 1) {

  fit = lapply(1:4, function(i){
    fittedSim[[i]][[Vintercept]][[commondisp]][[sim]]
  })

  fq <- betweenLaneNormalization(t(simData[[sim]]$counts), which="full")
  pca <- prcomp(t(log1p(fq)))

  dest = lapply(1:4, function(i) as.matrix(dist(fit[[i]]@W)))
  diff = lapply(1:4, function(i) eval_diff(dtrue, dest[[i]]))
  corr = lapply(1:4, function(i) eval_cor(dtrue, dest[[i]]))
  sil = lapply(1:4, function(i) eval_sil(biotrue, dest[[i]]))

  dest[[5]] <- as.matrix(dist(pca$x[,1:2]))
  diff[[5]] <- eval_diff(dtrue, dest[[5]])
  corr[[5]] <- eval_cor(dtrue, dest[[5]])
  sil[[5]] <- eval_sil(biotrue, dest[[5]])

  retval <- list(diff[[1]], diff[[2]], diff[[3]], diff[[4]], diff[[5]],
                 corr[[1]], corr[[2]], corr[[3]], corr[[4]], corr[[5]],
                 sil[[1]][,3], sil[[2]][,3], sil[[3]][,3], sil[[4]][,3],
                 sil[[5]][,3])
  return(retval)
}

load('sim_allen25.rda')
load('sim_allen25_fitted.rda')
true_W <- simModel@W
dtrue <- as.matrix(dist(true_W))
biotrue <- as.numeric(factor(bio))
strue <- silhouette(biotrue, dtrue)
res <- mclapply(1:10, function(i){
  eval_data(fittedSim, simData, sim = i)
})
save(res, file="sim_allen25_dist.rda")

load('sim_allen5.rda')
load('sim_allen5_fitted.rda')
true_W <- simModel@W
dtrue <- as.matrix(dist(true_W))
biotrue <- as.numeric(factor(bio))
strue <- silhouette(biotrue, dtrue)
res <- mclapply(1:10, function(i){
  eval_data(fittedSim, simData, sim = i)
})
save(res, file="sim_allen5_dist.rda")

load('sim_allen75.rda')
load('sim_allen75_fitted.rda')
true_W <- simModel@W
dtrue <- as.matrix(dist(true_W))
biotrue <- as.numeric(factor(bio))
strue <- silhouette(biotrue, dtrue)
res <- mclapply(1:10, function(i){
  eval_data(fittedSim, simData, sim = i)
i})
save(res, file="sim_allen75_dist.rda")


###############
# Model fit
###############
computeExp <- function(zinbModel, model = 'zinb'){
  if (model == 'zinb'){
    (1 - t(getPi(zinbModel))) * t(getMu(zinbModel))
  }else{
    t(getMu(zinbModel))
  }
  
}

computeVar <- function(zinbModel, model = 'zinb'){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zz))
  if (model == 'zinb'){
    (1 - pi) * mu * (1 + mu*(phi + pi))
  }else{
    mu * (1 + mu*phi)
  }
}

computeP0 <- function(zinbModel, model = 'zinb'){
  mu = t(getMu(zinbModel))
  pi = t(getPi(zinbModel))
  phi = exp(-getZeta(zinbModel))
  if (model == 'zinb'){
    pi + (1 - pi) * (1 + phi * mu) ^ (-1/phi)
  }else{
    (1 + phi * mu) ^ (-1/phi)
  }
}

plotMD <- function(x, y, xlim = c(0,10), ylim = c(-5, 5),
                   main = 'ZINB: MD-plot estimated vs. observed mean count, log scale'){
  mm = .5*(x + y)
  dd = x - y
  smoothScatter(mm, dd, xlim = xlim, ylim = ylim, xlab = 'Mean', ylab = 'Diff', main = main)
  abline(h = 0, col = 'gray')
  fit = loess(dd ~ mm)
  xpred = seq(0, 10, .1)
  pred = predict(fit, xpred)
  lines(xpred, pred, col = 'red', type = 'l', lwd=2)
}

zz = fittedSim[[2]][[1]][[2]][[1]]
zinbMu = log(rowMeans(computeExp(zz, model = 'zinb')))
obsMu = log(rowMeans(originalCounts))
nbMu = log(rowMeans(computeExp(zz, model = 'nb')))
zinbPi = rowMeans(computeP0(zz, model = 'zinb'))
obsPi = rowMeans(originalCounts == 0)
nbPi = rowMeans(computeP0(zz, model = 'nb'))

par(mfrow=c(1,2))
plotMD(obsMu, zinbMu, xlim = c(0, 60), ylim = c(-100, 10),
       main = 'ZINB: MD-plot estimated vs. observed mean count, log scale')
plotMD(obsMu, zinbMu, xlim = c(2, 9), ylim = c(-4, 2),
       main = 'ZINB: MD-plot estimated vs. observed mean count, log scale')

plotMD(obsMu, nbMu, xlim = c(0, 60), ylim = c(-110, 2),
       main = 'NB: MD-plot estimated vs. observed mean count, log scale')
plotMD(obsMu, nbMu, xlim = c(2, 9), ylim = c(-4, 2),
       main = 'NB: MD-plot estimated vs. observed mean count, log scale')

plotMD(obsPi, zinbPi, xlim = c(0, 1), ylim = c(-.2, .2),
       main = 'ZINB: MD-plot estimated vs. observed zero probability')
plotMD(obsPi, nbPi, xlim = c(0, 1), ylim = c(-.2, 1), 
       main = 'NB: MD-plot estimated vs. observed zero probability')

smoothScatter(obsMu, rowMeans(originalCounts == 0), xlab = 'log average count',
              ylab = 'proportion of zeros', main = 'Zero probability versus Mean')
fit = loess(zinbPi ~ zinbMu)
xpred = seq(1,10,.1)
pred = predict(fit, xpred)
lines(xpred, pred, col = 'red', type = 'l', lwd=2)
fit = loess(nbPi ~ nbMu)
xpred = seq(1,10,.1)
pred = predict(fit, xpred)
lines(xpred, pred, col = 'green', type = 'l', lwd=2)
fit = loess(rowMeans(originalCounts == 0) ~ obsMu)
xpred = seq(1,10,.1)
pred = predict(fit, xpred)
lines(xpred, pred, col = 'blue', type = 'l', lwd=2)

smoothScatter(rowMeans(originalCounts == 0), exp(-zz@zeta),
              xlab = 'Observed zero probability', ylab = 'Estimated dispersion',
              main = 'ZINB: dispersion versus zero probability')
fit = loess(exp(-zz@zeta) ~ rowMeans(originalCounts == 0))
xpred = seq(0, 1, .05)
pred = predict(fit, xpred)
lines(xpred, pred, col = 'red', type = 'l', lwd=2)










