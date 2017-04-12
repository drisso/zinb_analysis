library(zinbwave)

# if want to do it all at once, loop over ncells = c(100, 1000, 1000)
# and datasets = c('Allen', 'Zeisel')
nc = 100
ds = 'Zeisel'
ncores = 2

for (b2 in c(1, 5, 10)){
  for (offs in c(-1.5,0.5,2)){
    ff = sprintf('fig6ad-S13-S14/sim%s_nc%s_ratio%s_offs%s', ds, nc, b2, offs)
    load(paste0(ff, '.rda'))
    fittedSim = lapply(1:4, function(k){
      lapply(1:length(simData), function(i){
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) != 0, ]
        ngenes = nrow(counts)
        zinbFit(counts, K = k, commondispersion = FALSE,
                epsilon = ngenes, ncores = ncores)
      })
    })
    out = paste0(ff, '_fitted.rda')
    save(fittedSim, file = out)
  }
}
