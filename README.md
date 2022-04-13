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

## Datasets

###  Cellsius
A benchmark dataset of ~ 12,000 single-cell transcriptomes from eight human cell lines. The eight human cell lines were individually profiled by bulk RNA-seq, and mixed in four batches containing mixtures of two or three cell lines each for scRNA-seq profiling.

###### Batch1: IMR90 and HCT116 (50/50)
- IMR90 is a fibroblast cell line, isolated from fetal lung. Female.
- HCT116 is from human colon carcinoma with epithelium-like morphology. Male.

###### Batch2: A549 and Ramos (50/50)
- A549 is from human lung carcinoma, cell type is epithelial. Male.
- Ramos cells are from Burkitt’s lymphoma. They are lymphoblasts with B-cell characteristics. Ramos cells are very small (7-10um), so we usually find that they have fewer detected features and lower total count than other cell lines. Male.
 
###### Batch3: HEK293 and H1437 (50/50)
- HEK293 is a cell line form human embryonic kidney cells. Female.
- H1437 is from lung adenocarcinoma (i.e. origin is epithelial / glandular). Male.

###### DA234 (Batch 4): Jurkat, K562, Ramos (40% Jurkat, 55% K562 and 5% Ramos)
- Jurkat is a T-cell lymphoblast cell line. Male.
- K562 is a lymphoblast cell line wih granulocyte/erythrocyte/monocyte characteristics (fairly undifferentiated). Female.

#### Cell-type annotation:
Correlation of the single-cell to bulk expression profiles was used for cell type assignment, & Single cells were assigned to the cell type correlating most with their expression profile. Cells were excluded if their z-score correlation < 0.2, or if they correlated strongly with more than one bulk expression profile (likely doublets).

#### Subsets

| Cell-type  | Complete| Subset 1 | Subset 2|
| ------------- | ------------- | ------------- | ------------- |
| HCT116  | 1743  | 1400  | 1600  |
| HEK293  | 2002  | 1600  | 2000  |
| IMR90  | 1039  | 500  | 100 |
| A549  | 1320  | 400  | 80  |
| Ramos  | 1892  | 350  | 125  |
| H1437  | 1116  | 270  | 3  |
| K562  | 1606  | 380  | 70  |
| Jurkat  | 962  | 100  | 6  |

Datasets are pre-annotated with cell_idx, Batch, cell_line, cell_cycle_phase, gene names etc. and a range of QC metrics (would not necessarily trust).
Data is downloaded as an R data object, and were subsequently processed in R. Then convereted from seurat to anndata object using SCEasy.

##### Final datasets are located in RDM under code/DimensionalityReduction_Aim2/data/Cellsius/:
- Cellsius_Complete_Raw(sceasy).h5ad (Full dataset of all 8 cell lines, only pre-cursor filtering)
- Cellsius_Subset1_Raw(sceasy).h5ad (Subset 1 dataset of all 8 cell lines, only pre-cursor filtering)
- Cellsius_Subset2_Raw(sceasy).h5ad (Subset 2 dataset of all 8 cell lines, only pre-cursor filtering)
- subset1_sce_cleaned(SCEeasy).h5ad (Subset 1 dataset of all 8 cell lines, pre-cursor + some additional filtering)
- subset2_sce_cleaned(sceasy).h5ad (Subset 2 dataset of all 8 cell lines, pre-cursor + some additional filtering)

*None of the datasets have been normalised/transformed/scaled*

##### Filtering:
Precursor (done by authors prior to uploading data publically):
- ≥ 10.5 genes per cell [log2]
- ≥ 12.0 total UMIs / cell [log2]
- ≥ 10% mitochondrial genes

Additional:
- Outliers
- ≥ 3 counts in at least 1 cell

Sourced from:https://zenodo.org/record/3238275#.YWYVKBx_VhE

Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1739-7 
