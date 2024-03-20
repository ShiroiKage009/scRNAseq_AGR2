# This file does the automatic lableing for colon using whatever dataset we designate

# Importing the libraries
import scanpy as sc
import umap
import anndata as ad
import decoupler as dc
import numpy as np
import pandas

#%% reading files
duo_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_colon.h5ad")

#%% preprocessing

sc.pp.filter_cells(duo_unfilt,
                   min_counts = 2000)

sc.pp.filter_cells(duo_unfilt,
                   min_genes= 700) # best practices says 700 but we're keeping it at this for now

sc.pp.filter_genes(duo_unfilt,
                   min_cells = 50)

duo_mt_filt = duo_unfilt[duo_unfilt.obs['pct_counts_mt'] > 11]

sc.pp.normalize_total(duo_mt_filt, target_sum=1e4)

sc.pp.log1p(duo_mt_filt)

sc.pp.highly_variable_genes(duo_mt_filt, n_top_genes=2000)

duo_mt_filt_var = duo_mt_filt[:, duo_mt_filt.var.highly_variable]

#%% Dim red by PCA and tSNE

sc.pp.scale(duo_mt_filt_var, max_value = 10)
sc.tl.pca(duo_mt_filt_var, svd_solver = 'arpack')
sc.tl.tsne(duo_mt_filt_var)

#%% Dim UMAP

sc.pp.neighbors(duo_mt_filt_var)
# sc.tl.leiden(duo_mt_filt_var,resolution = 0.8)
sc.tl.leiden(duo_mt_filt_var)
sc.tl.umap(duo_mt_filt_var)

#%% Plotting

sc.set_figure_params(dpi = 600)
sc.pl.umap(duo_mt_filt_var)

#%%

sc.pl.umap(duo_mt_filt_var, color=["leiden"], ncols=2, size = 5)
sc.tl.rank_genes_groups(duo_mt_filt_var, groupby = 'leiden', method = 'wilcoxon')
sc.pl.rank_genes_groups(duo_mt_filt_var, groupby = 'leiden')

#%%
markers = dc.get_resource('PanglaoDB')
# Filter by canonical_marker and human
markers2 = markers[(markers['human']=='True')&(markers['canonical_marker']=='True')]

# Remove duplicated entries
markers2 = markers2[~markers2.duplicated(['cell_type', 'genesymbol'])]

print("end of cell")

#%%
dc.run_ora(mat = duo_mt_filt_var,
           net = markers2,
           source = 'cell_type',
           target = 'genesymbol',
           verbose = 1,
           use_raw = 0,
           min_n = 3
           )

duo_mt_filt_var.obsm['ora_estimate']

#%%
acts = dc.get_acts(duo_mt_filt_var, obsm_key = 'ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e


#%%
sc.pl.umap(acts, color=['Fibroblasts', 'leiden'])
# sc.pl.umap(acts, color=['cell_type']) #I don't think there is a cell type by this point. This would work for the previously-labeled data

#%%
test_ranks = dc.rank_sources_groups(acts, groupby = 'leiden', reference = 'rest', method = 't-test_overestim_var')
test_ranks

n_ctypes = 3
ctypes_dict = test_ranks.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

sc.pl.matrixplot(acts, ctypes_dict, 'leiden', dendrogram=True, standard_scale='var', colorbar_title='Z-scaled scores', cmap='RdBu_r')

sc.pl.violin(acts, keys=['Enterocytes', 'Goblet cells', 'Enteroendocrine cells', 'Fibroblasts'], groupby='leiden')
sc.pl.violin(acts, keys=['T cells', 'B cells', 'Platelets', 'Monocytes', 'NK cells'], groupby='leiden')

annotation_dict = test_ranks.groupby('group').head(1).set_index('group')['names'].to_dict()
annotation_dict

# Add cell type column based on annotation
acts.obs['cell_type'] = [annotation_dict[clust] for clust in acts.obs['leiden']]

# Visualize
sc.pl.umap(acts, color='cell_type')
#%%

sc.pl.violin(acts, keys=['Pericytes', 'Goblet cells', 'Enteroendocrine cells', 'Podocytes', 'Ductal cells'], groupby='leiden')
sc.pl.violin(acts, keys=['T cells', 'B cells', 'Platelets', 'Monocytes', 'NK cells'], groupby='leiden')
