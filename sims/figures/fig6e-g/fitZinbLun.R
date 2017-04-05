library(zinbwave)


for (ncells in c(100)){
  for (add in c('', '_ziadd0.33', '_ziadd0.67')){
    fileName = sprintf('simLun_%s%s', ncells, add)
    load(paste0(fileName, '.rda'))
    fittedSim = lapply(c(2), function(k){
      lapply(1:length(simData), function(i){
        if (add == ''){
          counts = simData[[i]]
        }else{
          counts = simData[[i]]$counts
        }
        counts = counts[rowSums(counts) != 0, ]
        ngenes = nrow(counts)
        zinbFit(counts, K = k, commondispersion = FALSE, ncores = 2)
      })
    })
    out = paste0(fileName, '_fitted.rda')
    save(fittedSim, file = out)
  }
}








