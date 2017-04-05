library(zinbwave)

makeZinbFit <- function(Xintercept = T, Vintercept = T, K = 2,
                        commondispersion = T, ngenes = 1000, ncells = 100){
  # X intercept or not
  Xint = as.numeric(Xintercept)
  X = matrix(Xint, ncol = 1, nrow = ncells)
  
  # V intercept or not
  Vint = as.numeric(Vintercept)
  V = matrix(Vint, ncol = 1, nrow = ngenes)
  
  buildZinb <- function(Y, ...){
    zinbFit(Y, X = X, V = V, K = K, commondispersion = commondispersion,
            epsilon = ngenes, ...)
  }
  return(buildZinb)
}

K = 1:4
Vintercept = c(TRUE, FALSE)
commondispersion = c(TRUE, FALSE)
ncores = 10
ds = 'Zeisel'
nc = 1000
b2 = 1
offs = 2

pp = sprintf('sim%s_nc%s_ratio%s_offs%s', ds, nc, b2, offs)
load(paste0(pp,".rda"))
fittedSim = lapply(K, function(k){
  lapply(Vintercept, function(Vint){
    lapply(commondispersion, function(commondisp){
      mclapply(1:length(simData), function(i){
        counts = t(simData[[i]]$counts)
        counts = counts[rowSums(counts) !=0, ]
        #zeros = (rowSums(counts) == 0)
        #if (sum(zeros) > 0){
        #  mm =  matrix(0, ncol = ncol(counts), nrow = sum(zeros)) 
        #  mm = sapply(1:nrow(mm), function(i){
        #    zz = rep(0, ncol(mm))
        #    zz[sample(length(zz), 1)] = 1
        #    zz
        #  })
        #  counts[zeros, ] = t(mm)
        #}
        myZinbFit = makeZinbFit(Xintercept = TRUE, Vintercept = Vint,
                                K = k, commondispersion = commondisp,
                                ngenes = nrow(counts), 
                                ncells = ncol(counts))
        myZinbFit(counts, ncores = ncores)
      },mc.cores =  ncores)
    })
  })
})
out = paste0(pp, '_fittedAll.rda')
save(fittedSim, file = out)

