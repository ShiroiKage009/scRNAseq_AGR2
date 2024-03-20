# Manual labeling antrum

# Importing the libraries
import scanpy as sc
import pandas as pd
import numpy as np

#%% Changing settings
sc.settings.verbosity = 3
results_path = 'test_results.h5ad'
sc.set_figure_params(dpi = 600)

#%% reading files
ant_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_antrum.h5ad")

#%% preprocessing
sc.pp.filter_cells(ant_unfilt, min_counts = 2000)
sc.pp.filter_cells(ant_unfilt, min_genes= 700) # best practices says 700 but we're keeping it at this for now
sc.pp.filter_genes(ant_unfilt, min_cells = 50)
ant_mt_filt = ant_unfilt[ant_unfilt.obs['pct_counts_mt'] < 10]

sc.pp.normalize_total(ant_mt_filt, target_sum=1e4)
sc.pp.log1p(ant_mt_filt)
sc.pp.highly_variable_genes(ant_mt_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
sc.pl.highly_variable_genes(ant_mt_filt)

ant_mt_filt.raw = ant_mt_filt

ant_mt_filt = ant_mt_filt[:, ant_mt_filt.var.highly_variable]

sc.pp.regress_out(ant_mt_filt, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(ant_mt_filt, max_value=10)

#%% PCA and tSNE
sc.tl.pca(ant_mt_filt, svd_solver='arpack')
sc.tl.tsne(ant_mt_filt)

sc.pl.pca(ant_mt_filt, color='MUC6')
sc.pl.tsne(ant_mt_filt, color='MUC6')

sc.pl.pca_variance_ratio(ant_mt_filt, log=True)

#%% UMAPping the UMAP

sc.pp.neighbors(ant_mt_filt)
sc.tl.leiden(ant_mt_filt, resolution = 0.8)
sc.tl.paga(ant_mt_filt)
sc.pl.paga(ant_mt_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(ant_mt_filt, init_pos='paga')

sc.tl.umap(ant_mt_filt)
sc.pl.umap(ant_mt_filt, color = ['leiden'])

#%% Gene expression ranking per cluster

sc.tl.rank_genes_groups(ant_mt_filt, 'leiden', method='t-test')
sc.pl.rank_genes_groups(ant_mt_filt, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(ant_mt_filt, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(ant_mt_filt, n_genes=25, sharey=False)

ant_mt_filt.write(results_path)

pd.DataFrame(ant_mt_filt.uns['rank_genes_groups']['names']).head(20)

result = ant_mt_filt.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd_df = pd.DataFrame(
             {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(20)

pd_df_df = pd.DataFrame(pd_df)
pd_df_df.to_csv('markers.csv')


#%% Intermediate testing for subsetting


if 'LGR5' in ant_mt_filt.var['highly_variable']:
    print("true")
else:
    print ("false")

# Get the index of highly variable genes
highly_variable_genes_index = ant_mt_filt.var['highly_variable']

# Subset the AnnData object based on highly variable genes
subset_ant_mt_filt = ant_mt_filt[:, highly_variable_genes_index]

# Now subset_ant_mt_filt contains only the highly variable genes.
ant_mt_filt_subset = ant_mt_filt[ant_mt_filt[:, 'LGR5'].X > 4]

sc.pl.umap(ant_mt_filt_subset, color = ['LGR5', 'MKI67', 'MUC5AC', 'TFF3', 'MUC6'], use_raw = 1)
sc.pl.scatter(ant_mt_filt_subset, x = 'MKI67', y = 'LGR5')

#%% Testing cell

#sc.pl.rank_genes_groups_violin(ant_mt_filt, groups='17', n_genes=12)

testing_df = pd.DataFrame(ant_mt_filt.uns['rank_genes_groups']['names']).head(20)

sc.pl.umap(ant_mt_filt, color = ['DEF6', 'LYZ', 'CD3D', 'CD19', 'NCAM1'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['MUC6', 'MUC5AC', 'ANPEP', 'LYZ', 'LGR5', 'MKI67', 'MUC2'], use_raw = 1)

sc.pl.umap(ant_mt_filt, color = ['BMPR1A', 'LGR5', 'PTEN', 'SMAD4', 'MKI67'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['LGR5'], use_raw = 0)

sc.pl.umap(ant_mt_filt, color = ['GKN2', 'MUC5AC', 'TFF2', 'MKI67'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['GKN2', 'MUC5AC', 'TFF2', 'MUC6', 'FOXQ1'], use_raw = 1)



#%%
sc.pl.scatter(ant_mt_filt, x='MUC5AC', y='GKN2', color = 'MKI67')
sc.pl.umap(ant_mt_filt, color = ['leiden'], use_raw = 0)

sc.pl.umap(ant_mt_filt, color = ['ANPEP', 'MUC2', 'MUC6', 'MUC5AC', 'TFF3', 'LGR5'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['ANPEP', 'MUC2','TFF3', 'MUC6', 'MUC5AC', 'TFF2', 'GKN2', 'MUC5B', 'LGR5', 'MKI67', 'leiden'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['MUC6', 'leiden'], use_raw = 0)
sc.pl.umap(ant_mt_filt, color = ['CD19', 'CD4', 'CD8A', 'CD3D', 'NCAM1'], use_raw = 0)

sc.pl.umap(ant_mt_filt, color = ['AXIN2', 'MKI67', 'LGR5'], use_raw = 1)
sc.pl.umap(ant_mt_filt, color = ['TNFRSF19'], use_raw = 1)

sc.pl.scatter(ant_mt_filt, x = 'MUC5AC', y = 'MUC6', color = 'MKI67')

sc.pl.umap(ant_mt_filt, color = ['LRIG1'], use_raw = 1)

#%%

# import scFates as scf

# new_data = ant_mt_filt
# sc.pp.scale(new_data)
# sc.pp.pca(new_data)
# scf.tl.curve(new_data,Nodes=30,use_rep="X_pca",ndims_rep=2,)

# scf.tl.root(new_data,'LGR5')
# scf.tl.pseudotime(new_data)
# sc.pl.pca(new_data,color="t")
# scf.pl.trajectory(new_data,basis="umap",arrows=True,arrow_offset=3)

