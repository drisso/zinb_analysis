# setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims')
# Comment: I commented out the line above to make it work in all local copies
# The setwd() should be taken care of by the .Rproj file

library(cluster)
library(zinb)
library(EDASeq)
library(digest) #to use zifa R wrapper

# input data
# load("sim_allen25_fitted.rda")
# load("sim_allen25.rda")
# fittedSim is a list with order k (in 1:4), Xintercept (T or F),
# Vintercept (T or F) and commondispersion (T or F), simulated dataset.
# For example fittedSim[[3]][[1]][[2]][[10]] was fitted with k = 3,
# Vintercept = T, and commondispersion = F for simulated dataset 10.
# See code in run_sim.R

## Let's look at bias and MSE

## levels: k (4), Vint (2), commondisp (2), nsim (10)

computeBiasList <- function(prefix, K = 1:4, Vintercept = 1:2,
                            commondispersion = 1:2){
  load(paste0(prefix, '.rda'))
  load(paste0(prefix, '_fitted.rda'))
  
  biasList = lapply(K, function(k){
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
  biasList
}

computeVarianceList <- function(prefix, K = 1:4, Vintercept = 1:2,
                                commondispersion = 1:2){
  load(paste0(prefix, '.rda'))
  load(paste0(prefix, '_fitted.rda'))
  n = length(simData)
  varianceList <- lapply(K, function(k){
    tmp <- lapply(Vintercept, function(Vint){
      tmp <- lapply(commondispersion, function(commondisp){
        gamma_mu <- (1/(n-1)) * rowSums(sapply(seq_along(simData), function(i) {
          (fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_mu[1,] - simModel@gamma_mu[1,])^2
        }))
        
        beta_mu <- (1/(n-1)) * rowSums(sapply(seq_along(simData), function(i) {
          (fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_mu[1,] - simModel@beta_mu[1,])^2
        }))
        
        gamma_pi <- (1/(n-1)) * rowSums(sapply(seq_along(simData), function(i) {
          (fittedSim[[k]][[Vint]][[commondisp]][[i]]@gamma_pi[1,] - simModel@gamma_pi[1,])^2
        }))
        
        beta_pi <- (1/(n-1)) * rowSums(sapply(seq_along(simData), function(i) {
          (fittedSim[[k]][[Vint]][[commondisp]][[i]]@beta_pi[1,] - simModel@beta_pi[1,])^2
        }))
        
        theta <- (1/(n-1)) * rowSums(sapply(seq_along(simData), function(i) {
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
  varianceList
}

############
# BIAS
############
for (k in c(75)){
  biasList = computeBiasList(paste0('./datasets/sim_allen', k))
  plotbias <- unlist(unlist(unlist(biasList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
  bias = lapply(seq_along(plotbias), function(i){
    nn = strsplit(names(plotbias)[i], '\\.')[[1]]
    data.frame(bias = as.vector(plotbias[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
  })
  bias = data.frame(do.call(rbind, bias), stringsAsFactors = F)
  save(bias, file = sprintf("./datasets/sim_allen%s_bias.rda", k))
}

###########
# VARIANCE
###########
for (k in c(25, 5)){
  varianceList = computeVarianceList(paste0('./datasets/sim_allen', k))
  plotvariance <- unlist(unlist(unlist(varianceList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
  variance = lapply(seq_along(plotvariance), function(i){
    nn = strsplit(names(plotvariance)[i], '\\.')[[1]]
    data.frame(variance = as.vector(plotvariance[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
  })
  variance = data.frame(do.call(rbind, variance), stringsAsFactors = F)
  save(variance, file = sprintf("./datasets/sim_allen%s_variance.rda", k))
}

##############
# ZIFA
##############
wrapRzifa <- function(Y, block = F){
  # wrapper R function for ZIFA.
  # md5 hashing and temporary files are used not to re-run zifa 
  # if it has already be run on this computer.
  d = digest(Y, "md5")
  tmp = paste0(tempdir(), '/', d)
  write.csv(Y, paste0(tmp, '.csv'))
  
  if (!file.exists(paste0(tmp, '_zifa.csv'))){
    print('run ZIFA')
    bb = ifelse(block, '-b ', '')
    cmd = sprintf('python ../real_data/run_zifa.py %s%s.csv %s_zifa.csv', bb, tmp, tmp)
    system(cmd)
  }
  read.csv(sprintf("%s_zifa.csv", tmp), header=FALSE)
}

for (k in c(25, 5, 75)){
  load(sprintf('./datasets/sim_allen%s.rda', k))
  zifa = lapply(1:10, function(i){
    Y = log1p(t(simData[[i]]$counts))
    wrapRzifa(Y)
  })
  save(zifa, file = sprintf('./datasets/sim_allen%s_zifa.rda', k))
}


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

## comments: there is no guarantee on the order of the W's
eval_data <- function(fittedSim, simData, simModel, zifa,
                      k = 1:4, Vintercept = 1, commondisp = 2, sim = 1) {

  fit = lapply(k, function(i){
    fittedSim[[i]][[Vintercept]][[commondisp]][[sim]]
  })

  fq <- betweenLaneNormalization(t(simData[[sim]]$counts), which="full")
  pca <- prcomp(t(log1p(fq)))
  
  true_W <- simModel@W
  dtrue <- as.matrix(dist(true_W))
  biotrue <- as.numeric(factor(bio))

  dest = lapply(k, function(i) as.matrix(dist(fit[[i]]@W)))
  diff = lapply(k, function(i) eval_diff(dtrue, dest[[i]]))
  corr = lapply(k, function(i) eval_cor(dtrue, dest[[i]]))
  sil = lapply(k, function(i) eval_sil(biotrue, dest[[i]]))

  dest[[5]] <- as.matrix(dist(pca$x[,1:2]))
  diff[[5]] <- eval_diff(dtrue, dest[[5]])
  corr[[5]] <- eval_cor(dtrue, dest[[5]])
  sil[[5]] <- eval_sil(biotrue, dest[[5]])
  
  dest[[6]] <- dtrue
  diff[[6]] <- eval_diff(dtrue, dest[[6]])
  corr[[6]] <- eval_cor(dtrue, dest[[6]])
  sil[[6]] <- eval_sil(biotrue, dest[[6]])
  
  dest[[7]] <- as.matrix(dist(zifa[[sim]]))
  diff[[7]] <- eval_diff(dtrue, dest[[7]])
  corr[[7]] <- eval_cor(dtrue, dest[[7]])
  sil[[7]] <- eval_sil(biotrue, dest[[7]])

  retval <- list(diff[[1]], diff[[2]], diff[[3]], diff[[4]], diff[[5]],diff[[6]], diff[[7]],
                 corr[[1]], corr[[2]], corr[[3]], corr[[4]], corr[[5]],corr[[6]],corr[[7]],
                 sil[[1]][,3], sil[[2]][,3], sil[[3]][,3], sil[[4]][,3],
                 sil[[5]][,3], sil[[6]][,3], sil[[7]][,3])
  return(retval)
}

for (k in c(25, 5, 75)){
  load(sprintf('./datasets/sim_allen%s.rda', k))
  load(sprintf('./datasets/sim_allen%s_fitted.rda', k))
  load(sprintf('./datasets/sim_allen%s_zifa.rda', k))
  res <- mclapply(1:10, function(i){
    eval_data(fittedSim, simData, simModel, zifa, sim = i, k = 1:2)
  })
  save(res, file = sprintf("./datasets/sim_allen%s_dist.rda", k))
}


