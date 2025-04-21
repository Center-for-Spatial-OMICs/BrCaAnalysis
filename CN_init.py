# import spacec first
import spacec as sp

# silencing warnings
import warnings
warnings.filterwarnings('ignore')

#import standard packages
import os
import scanpy as sc

sc.settings.set_figure_params(dpi=80, facecolor='white')

output_dir = '/mnt/scratch2/Luke/BrCaSCIMAP/'
adata_path = './adata.h5ad'
adata = sc.read(adata_path)

# compute for CNs
# tune k and n_neighborhoods to obtain the best result
adata = sp.tl.neighborhood_analysis(
    adata, 
    unique_region = "batch", 
    cluster_col = "phenotype", 
    X = 'imagecol', Y = 'imagerow',
    k = 20, # k nearest neighbors
    n_neighborhoods = 20, #number of CNs
    elbow = True)