library(zinbwave)

ncores = 2

for (ncells in c(100, 1000, 10000)){
  for (add in c('_ziadd0', '_ziadd0.33', '_ziadd0.67')){
    fileName = sprintf('fig6e-g/simLun_%s%s', ncells, add)
    load(paste0(fileName, '.rda'))
    fittedSim = lapply(1:length(simData), function(i){
      counts = simData[[i]]$counts
      counts = counts[rowSums(counts) != 0, ]
      zinbFit(counts, K = 2, ncores = ncores, commondispersion = FALSE)
    })
    out = paste0(fileName, '_fitted.rda')
    save(fittedSim, file = out)
  }
}

