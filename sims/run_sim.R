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
ncores = max(1, detectCores() - 2)

load("sim_allen_lowZero.rda")
load("sim_allen_highZero.rda")
fittedSim = lapply(K, function(k){
  lapply(Vintercept, function(Vint){
    lapply(commondispersion, function(commondisp){
      myZinbFit = makeZinbFit(Xintercept = TRUE,
                              Vintercept = Vint,
                              K = k,
                              commondispersion = commondisp)
      mclapply(1:length(simData), function(i){
        myZinbFit(t(simData[[i]]$counts))
      },mc.cores = ncores)
    })
  })
})
save(fittedSim, file = 'sim_allen_lowZero_fitted.rda')
fittedSim = lapply(K, function(k){
  lapply(Vintercept, function(Vint){
    lapply(commondispersion, function(commondisp){
      myZinbFit = makeZinbFit(Xintercept = TRUE,
                              Vintercept = Vint,
                              K = k,
                              commondispersion = commondisp)
      mclapply(1:length(simData_more0), function(i){
        myZinbFit(t(simData_more0[[i]]$counts))
      },mc.cores = ncores)
    })
  })
})
save(fittedSim_more0, file = 'sim_allen_highZero_fitted.rda')


