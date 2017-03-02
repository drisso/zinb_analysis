library(zinb)

nc = 100
for (b2 in c(1, 5, 10)){
  for (offs in c(-1.5,0.5,2)){
    ff = sprintf('simZeisel_nc%s_ratio%s_offs%s', nc, b2, offs)
    load(paste0(ff, '.rda'))
    fittedSim = lapply(1:4, function(k){
      lapply(1:length(simData), function(i){
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) != 0, ]
        ngenes = nrow(counts)
        zinbFit(counts, K = k, commondispersion = FALSE,
                epsilon = ngenes, ncores = 5)
      })
    })
    out = paste0(ff, '_fitted.rda')
    save(fittedSim, file = out)
  }
}
