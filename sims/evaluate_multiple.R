# setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims')
# Comment: I commented out the line above to make it work in all local copies
# The setwd() should be taken care of by the .Rproj file

library(cluster)
library(zinb)
library(EDASeq)
library(digest) #to use zifa R wrapper
library(edgeR)
library(DESeq2)

# input data
# load("simAllen25_fitted.rda")
# load("simAllen25.rda")
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
ds = 'Allen'
for (k in c(25, 50, 75)){
  print(k)
  biasList = computeBiasList(paste0('./datasets/sim', ds, k))
  plotbias <- unlist(unlist(unlist(biasList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
  bias = lapply(seq_along(plotbias), function(i){
    nn = strsplit(names(plotbias)[i], '\\.')[[1]]
    data.frame(bias = as.vector(plotbias[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
  })
  bias = data.frame(do.call(rbind, bias), stringsAsFactors = F)
  save(bias, file = sprintf("./datasets/sim%s%s_bias.rda", ds, k))
}

###########
# VARIANCE
###########
for (k in c(25)){
  varianceList = computeVarianceList(paste0('./datasets/sim',ds, k))
  plotvariance <- unlist(unlist(unlist(varianceList, recursive=FALSE), recursive=FALSE), recursive=FALSE)
  variance = lapply(seq_along(plotvariance), function(i){
    nn = strsplit(names(plotvariance)[i], '\\.')[[1]]
    data.frame(variance = as.vector(plotvariance[[i]]), K=nn[1], V=nn[2], disp=nn[3], param=nn[4])
  })
  variance = data.frame(do.call(rbind, variance), stringsAsFactors = F)
  save(variance, file = sprintf("./datasets/sim%s%s_variance.rda", ds, k))
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

for (k in c(50, 75)){
  load(sprintf('./datasets/sim%s%s.rda', ds, k))
  zifa = lapply(1:10, function(i){
    Y = t(simData[[i]]$counts)
    Y = Y[rowSums(Y) != 0, ]
    Y = log1p(Y)
    wrapRzifa(Y)
  })
  save(zifa, file = sprintf('./datasets/sim%s%s_zifa.rda',ds, k))
}

# 75, allen, c(1:2, 2, 4:9, 9)
for (k in c(75)){
  load(sprintf('./datasets/sim%s%s.rda', ds,k))
  zifaFQ = lapply(c(1:2, 2, 4:9, 9), function(i){
    print(i)
    fq <- betweenLaneNormalization(t(simData[[i]]$counts), which="full")
    Y = log1p(fq)
    Y = Y[rowSums(Y) != 0, ]
    wrapRzifa(Y)
  })
  save(zifaFQ, file = sprintf('./datasets/sim%s%s_zifaFQ.rda', ds,k))
}

for (k in c(25, 50, 75)){
  load(sprintf('./datasets/sim%s%s.rda', ds,k))
  zifaTMM = lapply(1:10, function(i){
    counts = t(simData[[i]]$counts)
    y = DGEList(counts)
    y = calcNormFactors(y, method="TMM")
    tmm = t(t(counts) / (y$samples$lib.size * y$samples$norm.factors))
    Y = log1p(tmm)
    Y = Y[rowSums(Y) != 0, ]
    wrapRzifa(Y)
  })
  save(zifaTMM, file = sprintf('./datasets/sim%s%s_zifaTMM.rda', ds,k))
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
eval_data <- function(fittedSim, simData, simModel, zifa, zifaFQ, zifaTMM,
                      k = 1:4, Vintercept = 1, commondisp = 2, sim = 1) {
  # ini
  dest = corr = sil = list()
  
  ## TRUE
  true_W <- simModel@W
  dtrue <- as.matrix(dist(true_W))
  biotrue <- as.numeric(factor(bio))
  dest[[1]] <- dtrue
  #diff[[1]] <- eval_diff(dtrue, dest[[1]])
  corr[[1]] <- eval_cor(dtrue, dest[[1]])
  sil[[1]] <- eval_sil(biotrue, dest[[1]])
  
  ## ZINB
  fit = lapply(k, function(i){
    fittedSim[[i]][[Vintercept]][[commondisp]][[sim]]
  })
  dest[2:(length(k)+1)] = lapply(k, function(i) as.matrix(dist(fit[[i]]@W)))
  corr[2:(length(k)+1)] = lapply(k, function(i) eval_cor(dtrue, dest[[i+1]]))
  sil[2:(length(k)+1)] = lapply(k, function(i) eval_sil(biotrue, dest[[i+1]]))
  
  ## PCA
  m = length(dest) + 1
  counts = t(simData[[sim]]$counts)
  pca <- prcomp(log1p(t(counts)))
  dest[[m]] <- as.matrix(dist(pca$x[,1:2]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## PCA TC
  m = m + 1
  mult = sum(counts) / (ncol(counts) * nrow(counts))
  fact = colSums(counts)
  tc = mult * (t(counts) / fact)
  pcatc <- prcomp(log1p(tc))
  #plot(pcatmm$x,col=bio)
  dest[[m]] <- as.matrix(dist(pcatc$x[,1:2]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## PCA tmm normalized counts (edgeR)
  m = m + 1
  y = DGEList(counts)
  y = calcNormFactors(y, method="TMM")
  tmm <- t(counts) / (y$samples$lib.size * y$samples$norm.factors)
  pcatmm <- prcomp(log1p(tmm))
  #plot(pcatmm$x,col=bio)
  dest[[m]] <- as.matrix(dist(pcatmm$x[,1:2]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## PCA DESeq2
  #if (sum(rowSums(counts == 0) > 0) < 1000){
  #  m = m + 1
  #  se = SummarizedExperiment(list(counts), colData = DataFrame(rep(1, ncol(counts))))
  #  dds = DESeqDataSet(se, design = ~ 1 )
  #  dds = estimateSizeFactors(dds)
  #  pcadeseq2 <- prcomp(log1p(t(counts(dds, normalized=TRUE))))
  #  dest[[m]] <- as.matrix(dist(pcadeseq2$x[,1:2]))
  #  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  #  sil[[m]] <- eval_sil(biotrue, dest[[m]]) 
  #}
  
  ## PCA FQ
  m = m + 1
  fq <- betweenLaneNormalization(counts, which="full")
  pcafq <- prcomp(t(log1p(fq)))
  dest[[m]] <- as.matrix(dist(pcafq$x[,1:2]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## ZIFA
  m = m + 1
  dest[[m]] <- as.matrix(dist(zifa[[sim]]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## ZIFA TMM
  m = m + 1
  dest[[m]] <- as.matrix(dist(zifaTMM[[sim]]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])
  
  ## ZIFA FQ
  m = m + 1
  dest[[m]] <- as.matrix(dist(zifaFQ[[sim]]))
  corr[[m]] <- eval_cor(dtrue, dest[[m]])
  sil[[m]] <- eval_sil(biotrue, dest[[m]])

  retval <- c(corr, lapply(1:m, function(i) sil[[i]][, 3]))
  return(retval)
}

ds = 'Allen'
for (k in c(75)){
  print(k)
  load(sprintf('./datasets/sim%s%s.rda',ds, k))
  load(sprintf('./datasets/sim%s%s_fitted.rda',ds, k))
  load(sprintf('./datasets/sim%s%s_zifa.rda', ds, k))
  load(sprintf('./datasets/sim%s%s_zifaFQ.rda', ds, k))
  load(sprintf('./datasets/sim%s%s_zifaTMM.rda', ds, k))
  res <- mclapply(1:10, function(i){
    eval_data(fittedSim, simData, simModel, zifa, zifaFQ, zifaTMM,
              sim = i, k = 1:2)
  })
  save(res, file = sprintf("./datasets/sim%s%s_dist.rda", ds, k))
}


