***********
scProximitE
***********

The accurate identification of cell sub-populations is paramount to the quality of downstream analyses and overall
interpretations of single-cell RNA-seq (scRNA-seq) datasets but remains an ongoing challenge for the field.
The quality of single-cell clustering has been demonstrated to be largely dependent on the proximity metric
used to generate cell-to-cell distances provided as input for the clustering algorithm. Accordingly, proximity
metrics have been benchmarked for scRNA-seq clustering previously, with results averaged across numerous datasets
to identify a single metric for recommendation. However, the metric identified as ‘best-performing’ varies between
studies, and the performance of proximity metrics has been found differ significantly between individual datasets.
These results suggest that the unique structural properties of scRNA-seq datasets impacts the performance of proximity
metrics, and as such identification of a single best-performing metric may be impossible. To investigate this,
we developed a framework for the in-depth evaluation of 17 proximity metrics performance with respect to core
structural properties of scRNA-seq data, including sparsity, dimensionality, cell population distribution and rarity.
We find that clustering performance can be improved substantially by the selection of an appropriate proximity metric
and neighbourhood size for the structural properties of a given dataset, as well as suitable processing and dimensionality
reduction approaches. Furthermore, popular metrics such as Euclidean and Manhattan distance performed poorly in
comparison to several lessor applied metrics, suggesting the default metric of choice for many scRNA-seq methods
should be re-evaluated. Our findings highlight the critical nature of tailoring scRNA-seq analyses pipelines to
the dataset at hand and provide practical guidance for researchers looking to optimise cell similarity search
for the structural properties of their own data.

.. figure:: _static/Fig2.png
   :width: 1000
   :align: center


Datasets
========

Cellsius
--------
A benchmark dataset of ~ 12,000 single-cell transcriptomes from eight human cell lines. The eight human cell lines were individually profiled by bulk RNA-seq, and mixed in four batches containing mixtures of two or three cell lines each for scRNA-seq profiling.

Batch1: IMR90 and HCT116 (50/50)
--------------------------------
- IMR90 is a fibroblast cell line, isolated from fetal lung. Female.
- HCT116 is from human colon carcinoma with epithelium-like morphology. Male.

Batch2: A549 and Ramos (50/50)
--------------------------------
- A549 is from human lung carcinoma, cell type is epithelial. Male.
- Ramos cells are from Burkitt’s lymphoma. They are lymphoblasts with B-cell characteristics. Ramos cells are very small (7-10um), so we usually find that they have fewer detected features and lower total count than other cell lines. Male.

Batch3: HEK293 and H1437 (50/50)
--------------------------------
- HEK293 is a cell line form human embryonic kidney cells. Female.
- H1437 is from lung adenocarcinoma (i.e. origin is epithelial / glandular). Male.

DA234 (Batch 4): Jurkat, K562, Ramos (40% Jurkat, 55% K562 and 5% Ramos)
------------------------------------------------------------------------
- Jurkat is a T-cell lymphoblast cell line. Male.
- K562 is a lymphoblast cell line wih granulocyte/erythrocyte/monocyte characteristics (fairly undifferentiated). Female.

Cell-type annotation:
---------------------
Correlation of the single-cell to bulk expression profiles was used for cell type assignment, & Single cells were assigned to the cell type correlating most with their expression profile. Cells were excluded if their z-score correlation < 0.2, or if they correlated strongly with more than one bulk expression profile (likely doublets).

Subsets
-------

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

Final datasets are located in RDM under code/DimensionalityReduction_Aim2/data/Cellsius/:
-----------------------------------------------------------------------------------------

- Cellsius_Complete_Raw(sceasy).h5ad (Full dataset of all 8 cell lines, only pre-cursor filtering)
- Cellsius_Subset1_Raw(sceasy).h5ad (Subset 1 dataset of all 8 cell lines, only pre-cursor filtering)
- Cellsius_Subset2_Raw(sceasy).h5ad (Subset 2 dataset of all 8 cell lines, only pre-cursor filtering)
- subset1_sce_cleaned(SCEeasy).h5ad (Subset 1 dataset of all 8 cell lines, pre-cursor + some additional filtering)
- subset2_sce_cleaned(sceasy).h5ad (Subset 2 dataset of all 8 cell lines, pre-cursor + some additional filtering)

*None of the datasets have been normalised/transformed/scaled*

Filtering:
----------

Precursor (done by authors prior to uploading data publically):
- ≥ 10.5 genes per cell [log2]
- ≥ 12.0 total UMIs / cell [log2]
- ≥ 10% mitochondrial genes

Additional:
- Outliers
- ≥ 3 counts in at least 1 cell

Sourced from:https://zenodo.org/record/3238275#.YWYVKBx_VhE

Paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1739-7

Tang et al. (2019)
------------------

single-cell RNA-seq (Drop-seq) of 400 HCA2 (Human foreskin) fibroblasts at population doubling 38 (young group (PD38)), population doubling 48 (middle-age group (PD48)), and population doubling 71 (senescent group (RS)), as well as population doubling 38 cells irradiated at 50Gy & left for 20 days (IS). Resulting in a total of 1200 cells. Senescence in IS and RS cells was confirmed with Beta-Gal staining & EdU assay.

Datasets are annotated with 'orig.ident' (PD38, PD48, RS, IR), and a few QC metrics.
Data is downloaded as seperate count matrices for each population, and were loaded into R & merged into a Seurat object, followed by filtering for the 'Filtered' files, or nothing for the 'UNFiltered' files. The Seurat object was then convereted from seurat to anndata object using SCEasy.

#### Final data is available in RDM under code/DimensionalityReduction_Aim2/data/Tang_2019/:
- Tang_All_Filtered(sceasy).h5ad (This contains all four groups, has undergone filtering)
- Tang_All_UNFiltered(sceasy).h5ad (This contains all four groups)
- Tang_PD38_PD48_PD71_Filtered_(sceasy).h5ad (This only contains the time-based groups, has undergone filtering)
- Tang_PD38_PD48_PD71_UNFiltered_(sceasy).h5ad (This only contains the time-based groups)
- GSM3384106_LowPDCtrl.dge.txt.gz (Raw txt file containing count matrix for PD38)
- GSM3384107_LowPD50Gy.dge.txt.gz (Raw txt file containing count matrix for IR)
- GSM3384108_HighPDCtrl.dge.txt.gz (Raw txt file containing count matrix for PD48)
- GSM3384109_senescence.dge.txt.gz (Raw txt file containing count matrix for RS)

*None of the datasets have been normalised/transformed/scaled*

Filtering:
----------
Filtered cells for:
- ≥ 200 genes expressed per cell
- ≥ 10% mitochondrial genes

Filtered genes for:
- Expression in ≥ 10% of all cells.


Sourced from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119807

Paper: https://link.springer.com/article/10.1007/s13238-018-0591-y
