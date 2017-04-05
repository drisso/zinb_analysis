library(zinbwave)

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
ds = args[1]
nc = args[2]
b2 = strsplit(args[3], '-offs')[[1]][1]
offs = strsplit(args[3], '-offs')[[1]][2]

ff = sprintf('sim%s_nc%s_ratio%s_offs%s', ds, nc, b2, offs)
load(paste0(ff, '.rda'))
fittedSim = lapply(1:4, function(k){
  print(k)
  aa = lapply(1:2, function(i){
    print(i)
    set.seed(73839)
    counts = t(simData[[i]]$counts)
    counts = counts[rowSums(counts) != 0, ]
    ngenes = nrow(counts)
    zinbFit(counts, K = k, commondispersion = FALSE, verbose = T, epsilon = ngenes, ncores = 8)
  })
  save(aa, file = sprintf('%s_k%s.rda', ff, k))
  aa
})
out = paste0(ff, '_fitted.rda')
save(fittedSim, file = out)


