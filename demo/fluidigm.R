#########################################################################################################
################    FLUIDIGM data set
#########################################################################################################
#install scRNAseq from local repository
#install.packages("/Users/svetlanagribkova/GitHub/scRNAseq",repos=NULL,type="source")
library(scRNAseq)
data("fluidigm")
fluidigm
#object consisting of several data sets
assayNames(fluidigm) #"counts" and "fpkm"
dim(fluidigm)
#?assays gives information on mode of accession of data
#print assays(...) to pick assays and then $... to pick the data I need

names(colData(fluidigm))
length(colData(fluidigm)$Biological_Condition)
sum(colData(fluidigm)$Coverage_Type=="High")#65 cells
high.cov.cells=assays(fluidigm)$counts[,which(colData(fluidigm)$Coverage_Type=="High")]
colData(fluidigm)$Coverage_Type=="High"
colData(fluidigm)$Cluster2
    
ncol(high.cov.cells)#65
nrow(high.cov.cells)#26262

######## SIZE FACTORS ESTIMATION  ############
#to estimate size factors, find genes never 0
all1=(apply(fluidigm2==0,2,sum)==0) 
sum(all1) #81
#estimate size factors
sizeFgenes=fluidigm2[,all1==1]
#calculate geometric means of lines
geom.means.lines=apply(sizeFgenes^(1/ncol(sizeFgenes)),2,prod)
#divide each line by its geom mean
for(i in 1:nrow(sizeFgenes)){
    sizeFgenes[i,]=sizeFgenes[i,]/geom.means.lines[i]
}
#compute size factors
size.factors=apply(sizeFgenes,1,median)
###############################################
#find genes to be filtered out
#filter.out=apply(high.cov.cells>10,1,sum)<20 4 groups
filter.out=apply(high.cov.cells>10,1,sum)<10
sum(!filter.out) #6995
#exclude genes with too low expressions
fluidigm2=high.cov.cells[!filter.out,]
#fluidigm2bis=t(fluidigm2)/matrix(rep(size.factors,J),ncol=J)

#initialization with PCA
#gene.exp.zeroinf=t(fluidigm2)[,500:800]
gene.exp.zeroinf=t(fluidigm2)[,1:500]
#n=nrow(gene.exp.zeroinf)
#J=ncol(gene.exp.zeroinf)
PCA.init=prcomp(log(1+gene.exp.zeroinf),center=TRUE,scale.=TRUE)
plot(PCA.init$x[,1:2])
plot(PCA.init$x[,1:2],col=colData(fluidigm)$Cluster2[colData(fluidigm)$Coverage_Type=="High"])

zinbres=zinb(gene.exp.zeroinf)

#for(l in 1:J){
#theta0[l]=fitdistr(gene.exp.zeroinf[,l],densfun="negative binomial")$estimate[1]
#}

############## EWING DATA
test=read.table("/Users/svetlanagribkova/GitHub/zinb/demo/rawcounts.txt",header=F)
test=test[2:nrow(test),3:98]
exp.mat=matrix(as.numeric(t(test)),nrow=nrow(t(test)))
exp.mat=exp.mat[,apply(exp.mat!=0,2,sum)>20]
dim(exp.mat)
exp.mat.pca=prcomp(exp.mat,scale.=TRUE,center=TRUE)
n=nrow(exp.mat)
J=ncol(exp.mat)
