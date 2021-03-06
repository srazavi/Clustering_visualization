import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import os

sc.settings.verbosity = 3             				  # verbosity: errors (0), warnings (1), info (2), hints (3)
results_folder = 'write'
adata = sc.read_csv("dataset_cleaned.csv", first_column_names = True)

nn_number_list = [10, 12] 				  			  # a list of numbers of nearest neighbors to consider for clustering
resolution_list = [0.001, 0.010, 0.015, 0.02]		  # a list of resolutions to consider for cluster annotation

def main():
	# common pre-processing steps
	sc.pp.scale(adata)                      
	sc.tl.pca(adata, svd_solver='auto')
	# the next steps vary per hyper-parameter
	for n_neighbors in nn_number_list:
		#clustering
		nn_key_added = str(n_neighbors)+'_nn'
		umap_obsm_key = 'X_umap_'+nn_key_added
		sc.pp.neighbors(adata, method='umap',n_neighbors=n_neighbors, 
			n_pcs=20, key_added=nn_key_added)
		adata_copy = sc.tl.umap(adata, min_dist=1, spread=1, \
								neighbors_key=nn_key_added, copy=True)
		adata.obsm[umap_obsm_key] = adata_copy.obsm['X_umap']
		#cluster annotation
		for res in resolution_list:
			leiden_key_added = 'leiden_res_'+str(res)+'_'+nn_key_added
			sc.tl.leiden(adata, resolution=res, key_added=leiden_key_added, 
				neighbors_key=nn_key_added)

	os.mkdir('write') if not os.path.exists('write') else None
	adata.write(results_folder+'/results.h5ad')

if __name__== "__main__":
	main()