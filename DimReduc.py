import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import rapids_singlecell as rsc

adata = sc.read('adata.h5ad')
rsc.get.anndata_to_GPU(adata)

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
rsc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca')
rsc.tl.umap(adata)
rsc.tl.tsne(adata)

adata.write_h5ad('adata.h5ad')