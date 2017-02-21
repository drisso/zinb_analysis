library(zinb)

for (ncells in c(100, 1000, 10000)){
  fileName = sprintf('simLun_%s', ncells)
  load(paste0(fileName, '.rda'))
  fittedSim = lapply(1:4, function(k){
    lapply(1:length(simData), function(i){
      counts = simData[[i]]
      counts = counts[rowSums(counts) != 0, ]
      ngenes = nrow(counts)
      zinbFit(counts, K = k, commondispersion = FALSE,
              epsilon = ngenes, ncores = 4)
    })
  })
  out = paste0(fileName, '_fitted.rda')
  save(fittedSim, file = out)
}








