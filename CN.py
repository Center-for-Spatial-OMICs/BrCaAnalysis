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
n = 9
adata = sc.read(adata_path)

# compute for CNs
# tune k and n_neighborhoods to obtain the best result
adata = sp.tl.neighborhood_analysis(
    adata, 
    unique_region = "batch", 
    cluster_col = "phenotype", 
    X = 'imagecol', Y = 'imagerow',
    k = 20, # k nearest neighbors
    n_neighborhoods = n, #number of CNs
    elbow = False)

print(adata)

# to better visualize the cellular neighborhood (CN), we choose a color palette
# but if you set palette = None in the following function, it will randomly generate a palette for you
cn_palette = {
    0: '#829868',
    1: '#3C5FD7',
    2: '#44CB63',
    3: '#FDA9AA',
    4: '#E623B1',
    5: '#204F89',
    6: '#F28E2B',
    7: '#8C564B',  # Brownish tone
    8: '#17BECF',  # Cyan/Teal
}

# save the palette in the adata
adata.uns[f'CN_k20_n{n}_colors'] = cn_palette.values()

# plot CN to see what cell types are enriched per CN so that we can annotate them better
sp.pl.cn_exp_heatmap(
    adata, # anndata
    cluster_col = "phenotype", # cell type column
    cn_col = f"CN_k20_n{n}", # CN column
    palette=cn_palette, # color palette for CN
    savefig = True, # save the figure
    output_dir = output_dir, # output directory
    rand_seed = 1 # random seed for reproducibility
)

# Convert dict_values to a list
adata.uns[f'CN_k20_n{n}_colors'] = list(adata.uns[f'CN_k20_n{n}_colors'])

# Save the AnnData object
adata.write_h5ad('adataN.h5ad')