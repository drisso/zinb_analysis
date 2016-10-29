# Use the allen data to simulate with realistic parameters
library(zinb)

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

# load simulated datasets
load("k2_Xintercept_Vintercept.rda")

# 32 different sets of parameters 
K = 1:4
Vintercept = c(TRUE, FALSE)
Xintercept = c(TRUE, FALSE)
commondispersion = c(TRUE, FALSE)
ncores = max(1, detectCores() - 2)

fittedSim = lapply(K, function(k){
  lapply(Xintercept, function(Xint){
    lapply(Vintercept, function(Vint){
      lapply(commondispersion, function(commondisp){
        myZinbFit = makeZinbFit(Xintercept = Xint,
                                Vintercept = Vint,
                                K = k,
                                commondispersion = commondisp)
        mclapply(1:length(sim_data), function(i){
          myZinbFit(t(sim_data[[i]]$counts))
        },mc.cores = ncores)
      })
    })
  })
})
save(fittedSim, file = 'k2_Xintercept_Vintercept_fitted.rda')


