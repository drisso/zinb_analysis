library(zinbwave)
library(doParallel)
library(BiocParallel)

NCORES = 2
registerDoParallel(NCORES)
register(DoparParam())

for (ncells in c(100, 1000, 10000)){
  for (add in c('_ziadd0', '_ziadd0.33', '_ziadd0.67')){
    fileName = sprintf('fig6e-g/simLun_%s%s', ncells, add)
    load(paste0(fileName, '.rda'))
    fittedSim = lapply(1:length(simData), function(i){
      counts = simData[[i]]$counts
      counts = counts[rowSums(counts) != 0, ]
      zinbFit(counts, K = 2, commondispersion = FALSE)
    })
    out = paste0(fileName, '_fitted.rda')
    save(fittedSim, file = out)
  }
}

