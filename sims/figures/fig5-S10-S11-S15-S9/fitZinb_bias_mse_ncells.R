library(zinbwave)

ncores = c(1, 5, 5, 10, 20, 20)
nc = c(50, 100, 500, 1000, 5000, 10000)
ds = 'Zeisel'
b2 = 1
offs = 2
eps = 1000

lapply(1:length(nc), function(i){
  pp = sprintf('fig5-S10-S11-S15-S9/simZeisel_nc%s_ratio1_offs2', nc[i])
  load(paste0(pp, ".rda"))
  fittedSim = mclapply(1:length(simData), function(j){
    counts = simData[[j]]$counts
    counts = counts[rowSums(counts) != 0, ]
    zinbFit(counts, K = 2, commondispersion = FALSE,
            epsilon = eps, ncores = ncores[i])
  })
  out = paste0(pp, '_fitted.rda')
  save(fittedSim, file = out)
})

