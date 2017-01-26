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
    zinbFit(Y, X = X, V = V, K = K, commondispersion = commondispersion, ...)
  }
  return(buildZinb)
}

K = 1:4
Vintercept = c(TRUE, FALSE)
commondispersion = c(TRUE, FALSE)
ncores = 20

ds = 'Allen'
for (j in 1){
  for (l in c(25, 45, 75)){
    fileName = sprintf('sim%s_var%s_z%s', ds, j, l)
    load(paste0(fileName,".rda"))
    fittedSim = lapply(K, function(k){
      lapply(Vintercept, function(Vint){
        lapply(commondispersion, function(commondisp){
          mclapply(1:length(simData), function(i){
            counts = t(simData[[i]]$counts)
            zeros = (rowSums(counts) == 0)
            counts = counts[!zeros, ] 
            myZinbFit = makeZinbFit(Xintercept = TRUE, Vintercept = Vint,
                                    K = k, commondispersion = commondisp,
                                    ngenes = nrow(counts), 
                                    ncells = ncol(counts))
            ff = myZinbFit(counts, ncores = 2)
            
            
            ff
          },mc.cores =  ncores)
        })
      })
    })
    out = paste0(fileName, '_fitted.rda')
    save(fittedSim, file = out)
  }
}

