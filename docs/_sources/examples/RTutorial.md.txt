# Promimity metrics in R

This is an optional notebook that enables users to run proximity metrics that are only in R.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("WGCNA")
BiocManager::install("reticulate",force=TRUE)
BiocManager::install("devtools",force=TRUE)
BiocManager::install("LoomExperiment",force=TRUE)
BiocManager::install("SingleCellExperiment",force=TRUE)
BiocManager::install("Seurat",force=TRUE, dependencies=TRUE)

library(devtools)
devtools::install_github("skinnider/dismay")
library(devtools)
devtools::install_github("cellgeni/sceasy",force = TRUE)

library('Seurat')
library('WGCNA')
library('dismay')

```

Options for proximity metrics in R are:

- phi_s
- rho_p
- weighted_rank
- zi-kendall

From the implementation provided by the dismay package (https://github.com/skinnider/dismay). 

As dismay package computes metrics as similarities, they are converted to dissimilarities when merged with python metrics. 

```
library(Seurat)
library(reticulate)
# Note here is where we suggest you use the conda env with scproximite installed
use_condaenv('scproximite')
library(sceasy)

data_dir = '../data/R_datasets/'

datasets = c("Discrete_Rare_Processed","Discrete_Abundant_Processed",
             "Continuous_Rare_Processed","Continuous_Abundant_Processed")

metrics = c('kendall') # 'phi_s', 'rho_p', 'weighted_rank','zi_kendall'

# Convert the datsets to rds format
for (data in datasets){
  name = paste(data_dir,data,sep="")
  sceasy::convertFormat(paste(name,".h5ad", sep=""), from="anndata", to="seurat", outFile=paste(name,".rds", sep=""))
}

```

Computes a distance matrix for each metric listed on each dataset provided, and outputs them to the data directory as 
.csv files for integration into the anndata object in python.


```
for (data in datasets){
  dat <- readRDS(paste(data_dir,data,".rds",sep =""))
  x <- as.matrix(dat@assays$RNA@data)
  for (m in metrics){
    out = dismay(x, m)
    nm = paste(data,m,'.csv', sep="_")
    write.csv(out, file=paste(data_dir, nm, sep=""))
  }
}
```