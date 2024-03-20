# Manual labeling antrum

# Importing the libraries
import scanpy as sc
import pandas as pd

#%% Changing settings
sc.settings.verbosity = 3
results_path = 'test_results.h5ad'
sc.set_figure_params(dpi = 600)

#%% reading files
col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")

#%% preprocessing

sc.pp.filter_cells(col_org_unfilt, min_counts = 2000)

sc.pp.filter_cells(col_org_unfilt, min_genes= 700) # best practices says 700 but we're keeping it at this for now

sc.pp.filter_genes(col_org_unfilt, min_cells = 50)

col_org_mt_filt = col_org_unfilt[col_org_unfilt.obs['pct_counts_mt'] > 10]
sc.pp.normalize_total(col_org_mt_filt, target_sum=1e4)
sc.pp.log1p(col_org_mt_filt)
sc.pp.highly_variable_genes(col_org_mt_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
sc.pl.highly_variable_genes(col_org_mt_filt)

col_org_mt_filt.raw = col_org_mt_filt

col_org_mt_filt = col_org_mt_filt[:, col_org_mt_filt.var.highly_variable]

sc.pp.regress_out(col_org_mt_filt, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(col_org_mt_filt, max_value=10)

#%% PCA and tSNE
sc.tl.pca(col_org_mt_filt, svd_solver='arpack')
sc.tl.tsne(col_org_mt_filt)

sc.pl.pca(col_org_mt_filt, color='MUC6')
sc.pl.tsne(col_org_mt_filt, color='MUC6')

sc.pl.pca_variance_ratio(col_org_mt_filt, log=True)

#%% UMAPping the UMAP

sc.pp.neighbors(col_org_mt_filt)
sc.tl.leiden(col_org_mt_filt, resolution = 0.8)
sc.tl.paga(col_org_mt_filt)
sc.pl.paga(col_org_mt_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(col_org_mt_filt, init_pos='paga')

sc.tl.umap(col_org_mt_filt)
sc.pl.umap(col_org_mt_filt, color = ['leiden'])

#%% Gene expression ranking per cluster

sc.tl.rank_genes_groups(col_org_mt_filt, 'leiden', method='t-test')
sc.pl.rank_genes_groups(col_org_mt_filt, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(col_org_mt_filt, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(col_org_mt_filt, n_genes=25, sharey=False)

col_org_mt_filt.write(results_path)

pd.DataFrame(col_org_mt_filt.uns['rank_genes_groups']['names']).head(20)

result = col_org_mt_filt.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd_df = pd.DataFrame(
             {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(20)

pd_df_df = pd.DataFrame(pd_df)
pd_df_df.to_csv('markers.csv')


#%% Marker definitions


#%% Testing cell

#sc.pl.rank_genes_groups_violin(col_org_mt_filt, groups='17', n_genes=12)

testing_df = pd.DataFrame(col_org_mt_filt.uns['rank_genes_groups']['names']).head(20)

sc.pl.umap(col_org_mt_filt, color = ['MUC6', 'MUC5AC', 'ANPEP', 'LYZ', 'LGR5', 'MKI67', 'MUC2', 'XBP1', 'AGR2', 'DEFA1'], use_raw = 1)

sc.pl.umap(col_org_mt_filt, color = ['AXIN2', 'BMPR1A', 'LGR5', 'PTEN', 'SMAD4', 'MKI67', 'LGR4'], use_raw = 1)

sc.pl.umap(col_org_mt_filt, color = ['AXIN2', 'BMPR1A', 'LGR5', 'PTEN', 'SMAD4', 'MKI67', 'LGR4'], use_raw = 1)

sc.pl.umap(col_org_mt_filt, color = ['TFF3', 'MUC2', 'MUC5B', 'MKI67', 'ANPEP', 'LGR5'], use_raw = 1)

sc.pl.umap(col_org_mt_filt, color = ['AGR2'], use_raw = 1)





