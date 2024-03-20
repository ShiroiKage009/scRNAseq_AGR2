# Manual labeling duodenum

# Importing the libraries
import scanpy as sc
import pandas as pd

#%% Changing settings
sc.settings.verbosity = 3
results_path = 'test_results.h5ad'
sc.set_figure_params(dpi = 600)

#%% reading files
duo_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_duodenum.h5ad")

#%% preprocessing

sc.pp.filter_cells(duo_unfilt, min_counts = 2000)

sc.pp.filter_cells(duo_unfilt, min_genes= 700) # best practices says 700 but we're keeping it at this for now

sc.pp.filter_genes(duo_unfilt, min_cells = 50)

duo_mt_filt = duo_unfilt[duo_unfilt.obs['pct_counts_mt'] > 10]
sc.pp.normalize_total(duo_mt_filt, target_sum=1e4)
sc.pp.log1p(duo_mt_filt)
sc.pp.highly_variable_genes(duo_mt_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
sc.pl.highly_variable_genes(duo_mt_filt)

duo_mt_filt.raw = duo_mt_filt

duo_mt_filt = duo_mt_filt[:, duo_mt_filt.var.highly_variable]

sc.pp.regress_out(duo_mt_filt, ['total_counts', 'pct_counts_mt'])

sc.pp.scale(duo_mt_filt, max_value=10)

#%% PCA and tSNE
sc.tl.pca(duo_mt_filt, svd_solver='arpack')
sc.tl.tsne(duo_mt_filt)

sc.pl.pca(duo_mt_filt, color='MUC6')
sc.pl.tsne(duo_mt_filt, color='MUC6')

sc.pl.pca_variance_ratio(duo_mt_filt, log=True)

#%% UMAPping the UMAP

sc.pp.neighbors(duo_mt_filt)
sc.tl.leiden(duo_mt_filt, resolution = 0.8)
sc.tl.paga(duo_mt_filt)
sc.pl.paga(duo_mt_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(duo_mt_filt, init_pos='paga')

sc.tl.umap(duo_mt_filt)
sc.pl.umap(duo_mt_filt, color = ['leiden'])
sc.pl.umap(duo_mt_filt, color = ['CDX2', 'MUC2', 'LYZ', 'ANPEP', 'LGR5', 'MUC6'])
sc.pl.umap(duo_mt_filt, color = ['MUC2', 'LYZ', 'LGR5', 'MUC6'], use_raw = False)

sc.pl.umap(duo_mt_filt, color = ['MUC6', 'ANPEP', 'LYZ', 'LGR5', 'MKI67', 'MUC2', 'XBP1', 'AGR2', 'DEFA1'], use_raw = 1)

sc.pl.scatter(duo_mt_filt, x = 'MUC2', y = 'AGR2', color = 'XBP1')

sc.pl.scatter(duo_mt_filt, x = 'LYZ', y = 'MUC6')

#%% Gene expression ranking per cluster

sc.tl.rank_genes_groups(duo_mt_filt, 'leiden', method='t-test')
sc.pl.rank_genes_groups(duo_mt_filt, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(duo_mt_filt, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(duo_mt_filt, n_genes=25, sharey=False)

duo_mt_filt.write(results_path)

pd.DataFrame(duo_mt_filt.uns['rank_genes_groups']['names']).head(20)

result = duo_mt_filt.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd_df = pd.DataFrame(
             {group + '_' + key[:1]: result[key][group]
             for group in groups for key in ['names', 'pvals']}).head(20)

pd_df_df = pd.DataFrame(pd_df)
pd_df_df.to_csv('markers.csv')


#%% Marker definitions


#%% Testing cell

sc.pl.rank_genes_groups_violin(duo_mt_filt, groups='17', n_genes=12)

testing_df = pd.DataFrame(duo_mt_filt.uns['rank_genes_groups']['names']).head(20)

sc.pl.umap(duo_mt_filt, color = ['MUC2', 'TFF3'], use_raw = 0)

sc.pl.umap(duo_mt_filt, color = ['MUC2', 'TFF3'], use_raw = 1)

sc.pl.scatter(duo_mt_filt, x='MUC6', y='LYZ')
