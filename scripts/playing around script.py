# Import packages
import scanpy as sc
import pandas as pd
import anndata as ad

#%% function definitions
def process_for_UMAP(adata, normed = 0):
        sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
        sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
        sc.pp.filter_genes(adata, min_cells = 10) # Filter genes based on the minimum number of cells expressing it
        adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
        if not normed:
            adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 50] # Filter on the cells with fewer than 10% mitochondrial reads
        else:
            adata_filt = adata_prefilt
        sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
        sc.pp.log1p(adata_filt) # Log scaling
        sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
        adata_filt.raw = adata_filt # Store the raw files in its own layer
        adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
        sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
        sc.pp.scale(adata_filt, max_value = 10) # Scale the data
        sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
        sc.tl.tsne(adata_filt)
        sc.pp.neighbors(adata_filt)
        sc.tl.leiden(adata_filt, resolution = 0.8)
        sc.tl.paga(adata_filt)
        sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
        sc.tl.umap(adata_filt, init_pos='paga')
        sc.tl.umap(adata_filt)
        sc.pl.umap(adata_filt, color = ['leiden'])
        return adata_filt

def process_until_norm(adata, cells):
        sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
        sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
        sc.pp.filter_genes(adata, min_cells = cells) # Filter genes based on the minimum number of cells expressing it
        adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 10] # Filter on the cells with fewer than 10% mitochondrial reads
        sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
        sc.pp.log1p(adata_filt) # Log scaling
        sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
        adata_filt.raw = adata_filt # Store the raw files in its own layer
        return adata_filt
    

def isolate_cells_by_gene(data, gene, threshold):
    # Now subset_ant_mt_filt contains only the highly variable genes
    data_subset = data[data[:, gene].X > threshold]
    
    return data_subset


    
#%% Read files
col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")
ant_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_antrum.h5ad")
duo_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_duodenum.h5ad")
col_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_colon.h5ad")

#%% Separating the thing
ant_patient = ant_unfilt[ant_unfilt.obs['Patient'] == 'P26']
ant_control = ant_unfilt[ant_unfilt.obs['Patient'] == 'GI6253']
duo_patient = duo_unfilt[duo_unfilt.obs['Patient'] == 'P26']


#%% process everything
col_org_filt = process_for_UMAP(col_org_unfilt)
ant_filt = process_for_UMAP(ant_patient)
ant_filt_cont = process_for_UMAP(ant_control)
duo_filt = process_for_UMAP(duo_unfilt)
col_filt = process_for_UMAP(col_unfilt)

#%%

sc.pl.umap(ant_filt, color = ['LRIG1'], use_raw = 1)

#%%

sc.pl.umap(ant_filt, color = ['leiden'])
sc.pl.umap(ant_filt, color = ['MUC5AC', 'TFF2', 'MUC6', 'LYZ', 'LGR5', 'MKI67'])
sc.pl.umap(ant_filt, color = ['ANPEP', 'MUC2', 'TFF3', 'LYZ', 'LGR5', 'MKI67'])
sc.pl.umap(ant_filt, color = ['CD3D', 'CD4', 'CD8A', 'CD19', 'NCAM1'])
sc.pl.umap(ant_filt, color = ['LGR5', 'MKI67', 'SMAD2', 'SMAD3', 'LRIG1'])
sc.pl.umap(ant_filt, color = ['KRT8', 'EPCAM'])

#%%
normed_ant = process_until_norm(ant_patient, cells = 25)
normed_ant_sub_MKI = normed_ant[(normed_ant[:, 'MKI67'].X > 0.5)]
normed_ant_sub_LGR5 = normed_ant[normed_ant[:, 'LGR5'].X > 0.5]
normed_ant_sub = ad.concat([normed_ant_sub_MKI, normed_ant_sub_LGR5], join='outer')
normed_ant_sub.obs['source'] = 'antrum'
normed_ant_sub_processed = process_for_UMAP(normed_ant_sub, normed = 1)

# sc.pl.umap(normed_ant_sub_processed, color = ['leiden', 'LGR5', 'MKI67', 'ANPEP', 'TFF2', 'TFF3'], use_raw = 1)
# sc.pl.tsne(normed_ant_sub_processed, color = ['leiden', 'LGR5', 'MKI67', 'ANPEP', 'TFF2', 'TFF3'], use_raw = 1)

#%%
normed_duo = process_until_norm(duo_patient, cells = 25)
normed_duo_sub_MKI = normed_duo[(normed_duo[:, 'MKI67'].X > 0.5)]
normed_duo_sub_LGR5 = normed_duo[normed_duo[:, 'LGR5'].X > 0.5]
normed_duo_sub = ad.concat([normed_duo_sub_MKI, normed_duo_sub_LGR5], join='outer')
normed_duo_sub.obs['source'] = 'duodenum'
normed_duo_sub_processed = process_for_UMAP(normed_duo_sub, normed = 1)

# sc.pl.umap(normed_duo_sub_processed, color = ['leiden', 'LGR5', 'MKI67', 'ANPEP', 'TFF3'], use_raw = 1)
# sc.pl.tsne(normed_duo_sub_processed, color = ['leiden', 'LGR5', 'MKI67', 'ANPEP', 'TFF2', 'TFF3'], use_raw = 1)

#%%

test_conc = ad.concat([normed_ant_sub, normed_duo_sub], join = 'outer')

processing_temp = process_for_UMAP(test_conc)
sc.pl.umap(processing_temp, color = ['leiden', 'LGR5', 'source'], use_raw = 1)


#%%

temp1 = normed_ant_sub_MKI
temp1.obs['source'] = 'antrum'

temp2 = normed_duo_sub_MKI
temp2.obs['source'] = 'duodenum'

temp_combine = ad.concat([temp1, temp2], join = 'outer')

temp_proc = process_for_UMAP(temp_combine)
sc.pl.umap(temp_proc, color = ['leiden', 'source', 'ANPEP'], use_raw = 1)
sc.pl.scatter(temp_proc, x = 'MKI67', y = 'ANPEP', color = 'source')

#%%

sc.tl.rank_genes_groups(temp_proc, 'leiden', method='t-test')
sc.pl.rank_genes_groups(temp_proc, n_genes=25, sharey=False)

sc.tl.rank_genes_groups(temp_proc, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(temp_proc, n_genes=25, sharey=False)