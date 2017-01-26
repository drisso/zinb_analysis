library(biomaRt)
library(scRNAseq)
library(matrixStats)

#### core
data("allen")
prefilter <- allen[grep("^ERCC-", rownames(allen), invert = TRUE),
                   which(colData(allen)$Core.Type=="Core")]
filterGenes = apply(assay(prefilter) > 5, 1, sum) >= 5
postfilter <- prefilter[filterGenes, ]
core <- assay(postfilter)
vars <- rowVars(log1p(core))
names(vars) <- rownames(core)
vars <- sort(vars, decreasing = TRUE)
core <- core[names(vars)[1:1000],]

### length
mart = useMart("ensembl")
mart = useDataset("hsapiens_gene_ensembl", mart = mart)
attrs = c("hgnc_symbol", "entrezgene")
bm = getBM(attributes=attrs, mart = mart)
geneName = toupper(rownames(core))
bm = bm[match(geneName, bm[,1]),]
bm = na.omit(bm)
gene_info = getGeneLengthAndGCContent(as.character(bm[,2]), "hg38", mode="org.db")
rownames(gene_info) = bm[,1]
gene_info = na.omit(gene_info)
head(geneInfo, 2)