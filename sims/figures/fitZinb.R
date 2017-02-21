library(zinb)

fileName = sprintf('simZeisel_10000_a1_offs3.5_seed9128')
load(paste0(fileName, ".rda"))
fittedSim = lapply(1:length(simData), function(i){
  counts = t(simData[[i]]$counts)
  counts = counts[rowSums(counts) != 0, ]
  ngenes = nrow(counts)
  ncells = ncol(counts)
  X = matrix(1, ncol = 1, nrow = ncells)
  V = matrix(1, ncol = 1, nrow = ngenes)
  zinbFit(counts, X = X, V = V, K = k, commondispersion = FALSE,
          epsilon = ngenes, ncores = 4)
})
out = paste0(fileName, '_fitted.rda')
save(fittedSim, file = out)


fileName = sprintf('simZeisel_10000_a1_offs3.5_seed9128')
load(paste0(fileName, ".rda"))
fittedSim = list()
for (k in 1:4){
  for (i in 1:10){
    counts = t(simData[[i]]$counts)
    counts = counts[rowSums(counts) != 0, ]
    ngenes = nrow(counts)
    ncells = ncol(counts)
    X = matrix(1, ncol = 1, nrow = ncells)
    V = matrix(1, ncol = 1, nrow = ngenes)
    fit = zinbFit(counts, X = X, V = V, K = k, commondispersion = FALSE,
                  epsilon = ngenes, ncores = 32)
    fittedSim[[k]][[i]] = fit
    rm(fit)
    rm(counts)
  }
}
out = paste0(fileName, '_fitted.rda')
save(fittedSim, file = out)








