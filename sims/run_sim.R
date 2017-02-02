library(zinb)

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
ncores = 20

ds = 'Allen'
nc = 100
for (aa in c(1, .85)){
  for (offs in c(-3.5, 0, 3.5)){
    pp = sprintf('simAllen_%s_a%s_offs%s_seed9128', nc, aa, offs)
    load(paste0(pp,".rda"))
    fittedSim = lapply(K, function(k){
      lapply(Vintercept, function(Vint){
        lapply(commondispersion, function(commondisp){
          mclapply(1:length(simData), function(i){
            counts = t(simData[[i]]$counts)
            zeros = (rowSums(counts) == 0)
            if (sum(zeros) > 0){
              mm =  matrix(0, ncol = ncol(counts), nrow = sum(zeros)) 
              mm[sample(length(mm), sum(zeros))] = 1
              counts[zeros, ] = mm
            }
            myZinbFit = makeZinbFit(Xintercept = TRUE, Vintercept = Vint,
                                    K = k, commondispersion = commondisp,
                                    ngenes = nrow(counts), 
                                    ncells = ncol(counts))
            myZinbFit(counts, ncores = 2)
          },mc.cores =  ncores)
        })
      })
    })
    out = paste0(pp, '_fittedAll.rda')
    save(fittedSim, file = out)
  }
}

