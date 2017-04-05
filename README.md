# Data analysis and simulations for the ZINB-WaVE paper

This repository is designed to allow interested people to reproduce the results and figures of the paper:

Risso D, Perraudeau F, Gribkova S, Dudoit S, Vert JP. ZINB-WaVE: A general and flexible method for signal extraction from single-cell RNA-seq data.

## Dependencies

In order to be able to run the code in this repo, it is required to have `R` (>=3.3) and `python` (>=2.7).

### R packages

- [zinbwave (R package)](https://github.com/drisso/zinbwave)
- cluster
- matrixStats
- magrittr
- RColorBrewer
- ggplot2
- reshape
- dplyr
- knitr
- rmarkdown
- mclust
- cowplot
- rARPACK
- Rtsne
- parallel
- digest


### Bioconductor packages

- EDASeq
- biomaRt
- scRNAseq
- SummarizedExperiment
- edgeR
- scran
- scater
- scone
- DESeq2


### python packages

- [ZIFA (python package)](https://github.com/epierson9/ZIFA)

## Getting started

### Real data

For each of the real datasets analyzed in the paper, there are a `.Rmd` file and a `.R` file in the `real_data` folder, e.g.,
for the Patel data, the files are [patel_covariates.Rmd](https://github.com/drisso/zinb_analysis/blob/master/real_data/patel_covariates.Rmd)
and [patel_plots.R](https://github.com/drisso/zinb_analysis/blob/master/real_data/patel_plots.R).

One needs to compile the .Rmd file first. This will have two effects: (i) it will create an HTML report with useful analyses of
the dataset; and (ii) it will create a `.rda` file with the results of `zinbwave`, `pca`, and `zifa`. Once this file is generated,
one can use the `.R` file to generate the dataset-specific plots found in the paper.

To generate the plots related to silhouette width, one needs to source the [silhouette.R](https://github.com/drisso/zinb_analysis/blob/master/real_data/silhouette.R)
file.

### Simulations

TODO
