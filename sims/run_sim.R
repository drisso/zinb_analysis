# Use the allen data to simulate with realistic parameters
library(zinb)
#setwd('~/Documents/BRAIN/gitrepo/zinb_analysis/sims/datasets/')

makeZinbFit <- function(Xintercept = T, Vintercept = T, K = 2,
                        commondispersion = T, ngenes = 1000, ncells = 100){
  # X intercept or not
  Xint = as.numeric(Xintercept)
  X = matrix(Xint, ncol = 1, nrow = ncells)
  
  # V intercept or not
  Vint = as.numeric(Vintercept)
  V = matrix(Vint, ncol = 1, nrow = ngenes)
  
  buildZinb <- function(Y, ...){
    zinbFit(Y, X = X, V = V, K = K, commondispersion = commondispersion, ...)
  }
  return(buildZinb)
}

K = 1:4
Vintercept = c(TRUE, FALSE)
commondispersion = c(TRUE, FALSE)
ncores = 10

load("sim_allen5.rda")
fittedSim = lapply(K, function(k){
  lapply(Vintercept, function(Vint){
    lapply(commondispersion, function(commondisp){
      mclapply(1:length(simData), function(i){
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) != 0, ] 
        myZinbFit = makeZinbFit(Xintercept = TRUE, Vintercept = Vint,
                                K = k, commondispersion = commondisp,
                                ngenes = nrow(counts), 
                                ncells = ncol(counts))
        myZinbFit(counts)
      },mc.cores =  ncores)
    })
  })
})
save(fittedSim, file = 'sim_allen5_fitted.rda')
