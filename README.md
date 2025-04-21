# BrCa Project

This repository contains the analysis code for the project.

## GenerateAnndata1.ipynb

This script creates the Anndata objects from QuPath output from cell segmentation. We first corrected artifacts that exist on certain channels of the image stacks. The resulting data from the corrected images are then incorporated into the Anndata objects we created in this script.

![Picture1](https://github.com/user-attachments/assets/40a7273b-bc83-4dd2-9b8a-b636998794c8)
![Picture2](https://github.com/user-attachments/assets/e31f5ebb-419d-4d1b-ae85-f18bd79521c5)

This script then performs QC on the samples, filter out cells with total fluorescence levels lower than 3 standard deviations below the mean and filter out cells with sizes that are 3 standard deviations higher than the mean. Lastly, we filter out cells with DAPI level less than a variable ammount to remove any cell segmentation artifacts.

![Unknown](https://github.com/user-attachments/assets/fc743039-1890-4d15-a43a-953df531be50)
![Unknown-1](https://github.com/user-attachments/assets/99707bfc-cd2d-48ec-aae5-50d91e10aa05)

## RunSciMap.py

This script runs the SCIMAP gating tool that performs cell annotation. To run the tool use the following command:

```bash
python RunSciMap.py INPUT.h5ad WORKFLOW.csv OUTPUT.h5ad
```

I ran SCIMAP separately on each sample alone which tend to return more uniform results and higher fidelity cell calling.

## RunHarmony.py

This script runs Harmony batch correction and stores the new PCA values in adata.obsm['X_pca_harmony'].

```bash
python RunHarmony.py INPUT.h5ad OUTPUT.h5ad
```

## DimReduc.py

This script runs n-neighbors, UMAP, and TSNE.

```bash
python DimReduc.py 
```

## PhenotypeAnalysis.ipynb

This script examines the annotation quality by examining the spatial maps of each cell type called.

<img width="987" alt="Screenshot 2025-04-21 at 2 21 52 PM" src="https://github.com/user-attachments/assets/ae08227d-e44c-4453-a46e-587dc8473054" />

This script examines samples individually and assigns the cell types back into the raw Anndata objects.

## CombineAnndatas.ipynb

This script takes all the separate Anndata objects and combines them together. This script generates the spatial mapping of phenotypes

![Unknown](https://github.com/user-attachments/assets/0576c5f4-528b-4c93-841f-d464f1daf0f5)

The UMAP/TSNE representation of the entire dataset

![Unknown-1](https://github.com/user-attachments/assets/f96f03fc-12ea-4fb6-8a84-2c34d876703c)

The cell proportion plot

![Unknown-2](https://github.com/user-attachments/assets/63fa3856-13e6-4087-bdd4-55a22e696436)

And the dot plot to examine cell type genetic expression

![Unknown-3](https://github.com/user-attachments/assets/1a9f55ab-d129-45d2-990e-07f4981ae276)

## CN_init.py

Run this script on the combined Anndata object to receive an elbow plot that shows the optimum number of neighborhoods for the dataset

```bash
python CN_init.py 
```

## CN.py

After determining the number of neighborhoods, run this script to conduct cell neighborhooding analysis.

```bash
python CN.py 
```

This script also generates the heatmap for which each cellular neighborhood is defined

![pdf](https://github.com/user-attachments/assets/d9872b09-d687-4c01-a59a-d122d3f8a470)

## PostCN.ipynb

This script generates the cellular neighborhood map

![Unknown](https://github.com/user-attachments/assets/815bbdd7-5eb2-43ec-8227-e5cbae41a4a0)

And the cellular neighborhood proportions

![Unknown-1](https://github.com/user-attachments/assets/9339662b-c34d-4d87-a87c-5a0880ba97c3)

## SpatialContext.ipynb

This script conducts the spatial context analysis on the cellular neighborhoods. In this script, I grouped the samples into their respective racial categories and generated the spatial context interaction plots for each racial group.

![Unknown](https://github.com/user-attachments/assets/f1bd6ac8-d5f8-4d82-81e1-1266e819401e)










