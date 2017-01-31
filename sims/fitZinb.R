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
    zinbFit(Y, X = X, V = V, K = K, commondispersion = commondispersion, epsilon = ngenes, ...)
  }
  return(buildZinb)
}

K = 1:4
Vint = TRUE
commondisp = FALSE
ncores = 40

for (nc in c(100, 1000, 10000)){
  for (aa in c(1, .85)){
    for (offs in c(-3.5, 0, 3.5)){
      fileName = sprintf('simAllen_%s_a%s_offs%s_seed9128', nc, aa, offs)
      load(paste0(fileName,".rda"))
      fittedSim = lapply(K, function(k){
            mclapply(1:length(simData), function(i){
              counts = t(simData[[i]]$counts)
              counts = counts[rowSums(counts) != 0, ] 
              myZinbFit = makeZinbFit(Xintercept = TRUE, Vintercept = Vint,
                                      K = k, commondispersion = commondisp,
                                      ngenes = nrow(counts), 
                                      ncells = ncol(counts))
              myZinbFit(counts, ncores = 4)
            }, mc.cores =  ncores)
          })
      out = paste0(fileName, '_fitted.rda')
      save(fittedSim, file = out)
    }
  }
}

