# Similarity Metrics at High Dimensionality - testing for rare cell types

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
