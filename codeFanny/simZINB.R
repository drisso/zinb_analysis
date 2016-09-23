# ref: "Pooling across cells to normalize single-cell 
# RNA sequencing data with many zero counts"
# Aaron T. L. Lun, Karsten Bach and John C. Marioni 

######################
# FUNCTIONS
######################
simulateZINB <- function(j = 100, s = 2, n = rep(10, 2), nbDisp = .1,
                         gammaShape = 2, gammaRate = 2, fc = rep(5, 2),
                         nDE = rep(10, 2), upDEprop = rep(.5, 2), 
                         cellBiasSD = .25, ZIprop = .5, ncores = NULL){
  
  stopifnot(upDEprop >= 0 | upDEprop <= 1)
  stopifnot(ZIprop >= 0 | ZIprop <= 1)
  
  #### mean NB
  # cell-specific bias for cell j
  log2theta = lapply(n, function(x) rnorm(x, 0, cellBiasSD))
  theta = lapply(log2theta, function(x) 2^x)
  
  # fold change
  lambdaBaseline = rgamma(j, shape = gammaShape, rate = gammaRate)
  fcList = lapply(1:s, function(x){
    fc = rep(1, j)
    if (nDE[x] != 0){
      de = sample(j, nDE[x])
      fc[de] = 0 # DE downregulated
      propUp = round(nDE[x] * upDEprop[x])
      if (propUp > 0){
        upDe = sample(de, propUp)
        fc[upDe] = FC[x] # DE upregulated
      }
    }
    list(fc = fc, upDe = upDe, downDE = de[!de %in% upDe])
  })
  lambda = lapply(lapply(fcList, '[[', 1), function(x) x * lambdaBaseline)
  NBmean = lapply(1:s, function(x) lambda[[x]] %*% t(theta[[x]]))
  NBmean = do.call(cbind, NBmean)
  
  #### counts
  if (is.null(ncores)) ncores = max(detectCores() - 2, 1)
  Y = mclapply(1:sum(n), function(y){
    sapply(1:j, function(x){
      rnbinom(1, mu = NBmean[x, y], size = nbDisp)  
    })
  }, mc.cores = ncores)
  Y = do.call(cbind, Y)
  
  #### ZI 
  sumZeroBefore = sum(Y == 0)
  stochasticZero = sample(length(Y), round(ZIprop * prod(dim(Y))))
  Y[stochasticZero] = 0
  stopifnot(sumZeroBefore <= sum(Y == 0))
  list(counts = Y, stochasticZero = stochasticZero,
       upDE = lapply(fcList, '[[', 2), downDE = lapply(fcList, '[[', 3))
}


#######################
# SIMULATE
#######################
J            = 10000 # number of genes
S            = 2 # number of subpopulations
N            = rep(250, S) #number of cells for each subpopulation
NB_DISP      = .5
SHAPE        = 2
RATE         = 2
FC           = rep(5, S) #fold change
DE           = rep(100, S) # number of DE genes
UP_DE        = rep(.2, S) # proportion of DE genes that are upregulated
CELL_BIAS_SD = .25 #gaussian
ZI           = .5 # prop of genes with technical zero

Y = simulateZINB(j = J, s = S, n = N, nbDisp = NB_DISP, gammaShape = SHAPE,
                 gammaRate = RATE, fc = FC, nDE = DE, upDEprop = UP_DE,
                 cellBiasSD = CELL_BIAS_SD, ZI = .5, ncores = 4)

#####################
# PLOTS
#####################
color = c(rep('red', 250), rep('blue', 250))

genezero = apply(Y$counts, 2, function(x) length(x[x == 0])/length(x))
par(lty = 0)
barplot(genezero, col = color, main = 'Perc. of genes with zero count', space=0)
par(lty = 1)

ngeneDetected = apply(Y$counts, 1, function(x) length(x[x != 0]))
plot(table(ngeneDetected), main= 'Number of cells in which a gene is detected')
plot(density(ngeneDetected), main = 'Number of cells in which a gene is detected')
abline(v = mean(ngeneDetected), col = 'red')

aveCounts = apply(Y$counts, 1, mean)
h = hist(aveCounts, plot = F)
h$counts = h$counts/sum(h$counts)
plot(h, xlim = c(0, 10), ylim = c(0, 1))
par(new = T)
plot(density(rgamma(J, SHAPE, RATE)), xlim = c(0,10), ylim = c(0, 1))


########################
# FIT
########################
library(pscl)
library(RUVSeq)

#norm
set <- newSeqExpressionSet(Y$counts,
                           phenoData = data.frame(group = c(rep(c(0,1), each=250))))
fq <- betweenLaneNormalization(set, which = "full")
plotPCA(fq, col = factor(pData(set)$group), labels = F, main = 'After FQ')

#zinb
dd = as.data.frame(cbind(t(Y$counts), layer = factor(pData(set)$group)))
zinb = zeroinfl(layer ~ . | ., data = dd)













