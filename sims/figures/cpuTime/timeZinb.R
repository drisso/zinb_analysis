library(zinb)

cpuTime = lapply(c(50, 100, 500, 1000, 5000, 10000), function(nc){
  fileName = sprintf('simZeisel_%s_a1_offs3.5_seed9128', nc)
  load(paste0(fileName, ".rda"))
  tt = lapply(1:10, function(j){
    counts = t(simData[[j]]$counts)
    counts = counts[rowSums(counts) != 0, ] 
    ngenes = nrow(counts)
    ncells = ncol(counts)
    X = matrix(1, ncol = 1, nrow = ncells)
    V = matrix(1, ncol = 1, nrow = ngenes)
    system.time(zinbFit(counts, X = X, V = V, K = 2, commondispersion = FALSE,
                        epsilon = ngenes, ncores = 4))
  })
  tt
})
save(cpuTime, file = 'cpuTime.rda')
