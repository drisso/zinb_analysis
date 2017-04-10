simulateW <- function(zinb, ncells = 100, nclust = 3, ratioSSW_SSB = 1, colIni = 1){
  par(mfrow = c(2,2))
  # zinbW
  xlim = c(min(zinb@W[,1]) - 1, max(zinb@W[,1]) + 1) 
  ylim = c(min(zinb@W[,2]) - 1, max(zinb@W[,2]) + 1)
  plot(zinb@W, col = colIni, xlim = xlim, ylim = ylim,
       main = 'zinb fitted W\ncolor = brain area')
  
  # mclustW
  mclustW = Mclust(zinb@W, G = nclust)
  plot(mclustW$data, main = 'multivar gauss fit\nmclust K = 3', xlab = 'W1', ylab = 'W2', 
       xlim = xlim, ylim = ylim, col = mclustW$classification)
  
  # multivar gaussian
  clust = sample(mclustW$classification, ncells, replace = TRUE)
  stopifnot(length(unique(clust)) == nclust)
  # a = b = 1
  simW1 = lapply(clust, function(i){
    mvrnorm(n = 1, mu = mclustW$parameters$mean[, i], 
            Sigma = mclustW$parameters$variance$sigma[,, i])
  })
  simW1 = do.call(rbind, simW1)
  plot(simW1, col = clust,
       main = paste0('multivar gauss sim\nncells=', ncells, ', scaleRatioSSW_SSB = 1'),
       xlab = 'W1', ylab = 'W2', xlim = xlim, ylim = ylim)
  
  simW2 = computeNewW(simW1, clust, ratioSSW_SSB)
  plot(simW2, col = clust,
       main = paste0('multivar gauss sim\nncells=', ncells, ', scaleRatioSSW_SSB =' , ratioSSW_SSB),
       xlab = 'W1', ylab = 'W2', xlim = xlim, ylim = ylim)
  par(mfrow = c(1, 1))
  
  return(list(simW = simW2, bio = clust))
}

computeNewW = function(W, labels, ratioSSW_SSB = 1){
  Vtot = apply(W, 2, var)
  W_bar = apply(W, 2, mean)
  N = nrow(W)
  
  cc = table(labels)
  ccNames = names(cc)
  nk = as.vector(cc)
  W_bar_c = sapply(as.numeric(ccNames), function(i){
    apply(W[labels == i, ], 2, mean)
  })
  colnames(W_bar_c) = ccNames
  SS_between = colSums(nk*t((W_bar_c - W_bar)^2))
  
  SS_within = rowSums(sapply(seq_len(N), function(i){
    (W[i, ] - W_bar_c[, as.character(labels[i])])^2
  }))

  b2 = ratioSSW_SSB
  a = sqrt( ( (N-1) * Vtot ) / ( SS_between + b2 * SS_within ) )
  W_start = sapply(seq_len(N), function(i){
    (1 - a) * W_bar +  W_bar_c[, as.character(labels[i])] * a * (1 - sqrt(b2)) +  a * sqrt(b2) * W[i,]
  })
  t(W_start)
}

simulateGamma <- function(zinb, ncells = 100, gammapiOffset = 0, colIni = 1,
                          colSim = 1){
  # gamma zinb
  gamma = data.frame(gammaMu = zinb@gamma_mu[1, ],
                     gammaPi = zinb@gamma_pi[1, ])
  # mclustW
  mclustGamma = Mclust(gamma, G = 1)
  
  # multivar gaussian
  simGamma = mvrnorm(n = ncells,
                     mu = mclustGamma$parameters$mean[,1] + c(0, gammapiOffset), 
                     Sigma = mclustGamma$parameters$variance$sigma[,, 1])
  
  par(mfrow = c(1,2))
  xlim = c(min(c(gamma[,1], simGamma[,1])) - .5,
           max(c(gamma[,1], simGamma[,1])) + .5) 
  ylim = c(min(c(gamma[,2], simGamma[,2])) - .5,
           max(c(gamma[,2], simGamma[,2])) + .5) 
  plot(gamma[,1], gamma[,2], col = colIni,
       xlim = xlim, ylim = ylim, xlab = 'gamma_mu', ylab = 'gamma_pi', 
       main = 'zinb fitted Gamma\ncolor = brain area')
  plot(simGamma, col = colSim,
       main = paste0('bivar gauss sim\nncells=', ncells, 
                     ', gammaPi offset = ', gammapiOffset),
       xlab = 'gamma_mu', ylab = 'gamma_pi', xlim = xlim, ylim = ylim)
  par(mfrow = c(1, 1))
  
  return(simGamma)
}

zinbSimWrapper <- function(core, colIni, ncells = 100, ngenes = 1000, nclust = 3, 
                           ratioSSW_SSB = 1, gammapiOffset = 0, B = 1, ncores = NULL,
                           fileName = 'zinbSim.rda'){
  # sample ngenes 
  set.seed(9128)
  if (ngenes > nrow(core)) repl = T else repl = F
  core = core[sample(1:nrow(core), ngenes, repl = repl),]
  
  # fit zinb (if you already fitted zinb, it is cached)
  d = digest(core, "md5")
  tmp = paste0(tempdir(), '/', d)
  fileZinb = sprintf("%s_zinb.rda", tmp)
  if (!file.exists(fileZinb)){
    print('run ZINB')
    if (is.null(ncores)) ncores = max(1, detectCores() - 1)
    zinb <- zinbFit(core, ncores = ncores, K = 2, commondispersion = FALSE,
                    epsilon = ngenes)
    save(zinb, file = fileZinb)
  }else{
    load(fileZinb)
  }
  
  # sim W
  w = simulateW(zinb, ncells, nclust, ratioSSW_SSB, colIni)
  simW = w$simW
  bio = w$bio
  
  # sim gamma
  simGamma = simulateGamma(zinb, ncells, gammapiOffset,
                           colIni = colIni, colSim = bio)
  
  # sim model
  simModel = zinbModel(W=simW, gamma_mu = matrix(simGamma[,1], nrow = 1),
                  gamma_pi = matrix(simGamma[,2], nrow = 1),
                  alpha_mu=zinb@alpha_mu, alpha_pi=zinb@alpha_pi,
                  beta_mu=zinb@beta_mu, beta_pi=zinb@beta_pi, zeta = zinb@zeta)
  
  # sim data
  if (B == 1){
    simData = zinbSim(simModel, seed = 1, no_cores=ncores)
  } else{
    simData = lapply(seq_len(B), function(j){
      zinbSim(simModel, seed = j, no_cores=ncores)
    })
  }
  
  save(bio, simModel, simData, file = fileName)
}

###############################################################################
library(scRNAseq)
library(zinbwave)
library(mclust)
library(digest)
library(RColorBrewer)
library(MASS)

##########
## ALLEN
##########
data("allen")
cols = brewer.pal(8, "Set1")
prefilter = allen[grep("^ERCC-", rownames(allen), invert = TRUE),
                   which(colData(allen)$Core.Type=="Core")]
filterGenes = apply(assay(prefilter) > 5, 1, sum) >= 5
postfilter = prefilter[filterGenes, ]
bioIni =  as.factor(colData(postfilter)$driver_1_s)
core = assay(postfilter)
colIni = cols[bioIni]


###########################################################################
## Correlation and Silhouette plot
## Figures 6, S13, S14 
for (nc in c(100, 1000, 10000)){
  for (b2 in c(1, 5, 50)){
    for (offs in c(0, 2, 5)){
      ff = sprintf('fig6ad-S13-S14/simAllen_nc%s_ratio%s_offs%s.rda', nc, b2, offs)
      zinbSimWrapper(core = core, colIni = colIni, ncells = nc, nclust = 3, 
                     ratioSSW_SSB = b2, gammapiOffset = offs, B = 10, 
                     fileName = ff, ncores = 2)
    }
  }
}



##########
## ZEISEL
##########
data <- read.table("../datasets/expression_mRNA_17-Aug-2014.txt", sep='\t',
                   stringsAsFactors = FALSE, comment.char = '%')
counts <- as.matrix(data[12:NROW(data),-(1:2)])
counts <- matrix(as.numeric(counts), ncol=ncol(counts), nrow=nrow(counts))
rownames(counts) <- data[12:NROW(data),1]
colnames(counts) <- data[8, -(1:2)]
level1 <- as.factor(as.matrix(data)[9,-(1:2)])
set.seed(21986)
filter = sample(1:ncol(counts), 2000, replace = F)
counts = counts[, filter]
level1 = droplevels(level1[filter])
filterGenes = apply(counts > 5, 1, sum) >= 5
counts <- counts[filterGenes, ]
col <- brewer.pal(8, "Set2")
colIni <- col[level1]

###################################################################
## Correlation and Silhouette plots
# Figures 6, S13, S14 
###################################
for (nc in c(100, 1000, 10000)){
  for (b2 in c(1, 5, 10)){
    for (offs in c(-1.5, 0.5, 2)){
      ff = sprintf('fig6ad-S13-S14/simZeisel_nc%s_ratio%s_offs%s.rda', nc, b2, offs)
      zinbSimWrapper(core = counts, colIni = colIni, ncells = nc, nclust = 3, 
                     ratioSSW_SSB = b2, gammapiOffset = offs, B = 10, 
                     fileName = ff, ncores = 2)
    }
  }
}


#######################################################################
# Mean-Difference plot
# Figure S12
######################
nc = 1000
b2 = 1
offs = 2
ff = sprintf('figS12/simZeisel_nc%s_ratio%s_offs%s.rda', nc, b2, offs)

## simulate
zinbSimWrapper(core = counts, colIni = colIni, ncells = nc, nclust = 3, 
               ratioSSW_SSB = b2, gammapiOffset = offs, B = 1, ncores = 2,
               fileName = ff)

# remove genes with only zeros and keep track of removed genes for mean diff
load(ff)
counts = t(simData$counts)
keep = rowSums(counts) > 0
simData$counts = simData$counts[,keep] 
save(bio, simModel, simData, keep, file = ff)

# fit
fittedSim = zinbFit(t(simData$counts), K = 2, commondispersion = FALSE, ncores = 2)
save(fittedSim, file = gsub('.rda', '_fitted.rda', ff))



#######################################################################
# bias_mse_allParam (figures 5, S10, S11) + bias_mse_ncells (figure S9) + cpuTime (figure S15)
################
b2 = 1
offs = 2
for (nc in c(50, 100, 500, 1000, 5000, 10000)){
  ff = sprintf('fig5-S10-S11-S15-S9/simZeisel_nc%s_ratio%s_offs%s.rda', nc, b2, offs)
  zinbSimWrapper(core = counts, colIni = colIni, ncells = nc, nclust = 3, 
                 ratioSSW_SSB = b2, gammapiOffset = offs, B = 10,
                 fileName = ff, ncores = 2)
  
  # remove genes with only zeros in at least one of the B simulated ds
  load(ff)
  keep = lapply(1:length(simData), function(i){
    counts = t(simData[[i]]$counts)
    rowSums(counts) != 0
  })
  keep = Reduce('+', keep) == 10
  for (i in 1:length(simData)){
    simData[[i]]$counts = simData[[i]]$counts[,keep] 
  }
  save(bio, simModel, simData, keep, file = ff)
}



