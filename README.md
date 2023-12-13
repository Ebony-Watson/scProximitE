# Similarity Metrics at High Dimensionality - testing for rare cell types
This package is designed for evaluating the performance of various proximity metrics (including distance, similarity, dissimilarity, correlation etc. metrics) with respect to quantifying cell-cell similarity in scRNA-seq datasets. The study for which the package was originally created and the performance of the metrics included in the package with respect to various dataset-specific properties of scRNA-seq data is available at https://doi.org/10.1093/bib/bbac387.

If relevant, please cite this package using the paper citation:
Ebony Rose Watson, Ariane Mora, Atefeh Taherian Fard, Jessica Cara Mar, How does the structure of data impact cell–cell similarity? Evaluating how structural properties influence the performance of proximity metrics in single cell RNA-seq data, Briefings in Bioinformatics, Volume 23, Issue 6, November 2022, bbac387, https://doi.org/10.1093/bib/bbac387


## Docs
Documentation and reproducibility are available at:

https://ebony-watson.github.io/scProximitE

## Install
```
pip install scproximite
```

Note: scproximite was developed using Python 3.8, of you have any issues we recommend using conda and creating a new
environment before installing:
```
conda create --name scproximite python=3.8
```
```
conda activate scproximite
```
```
pip install scproximite
```

## Run tutorials

1. Get tutorial data from zeonodo: https://zenodo.org/record/6443267 (DOI: 10.5281/zenodo.6443266)
2. Add to the `data/framework` folder
3. Run `jupyter notebook` in the `tutorials` folder

You should now be able to run the tutorial notebooks. Note if you don't have `R` installed you won't be able to 
run the notebook that uses `R` metrics: `Proximity_Metrics_R.ipynb`.
