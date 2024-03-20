# Import packages
import scanpy as sc
import anndata as ad

#%% function definitions
def process_for_UMAP(adata, normed = 0, leiden_res = 0.8):
    print("Filtering ...")
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = 10) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 50] # Filter on the cells with fewer than 10% mitochondrial reads
    else:
        adata_filt = adata_prefilt
    print("Normalizing ...")
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print("Log scaling ...")
    sc.pp.log1p(adata_filt) # Log scaling
    print("Finding variable genes ...")
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print("Saving raw data ...")
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    print("Filtering on variable genes ...")
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    print("Regressing ...")
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print("Scaling ...")
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print("Calculating PCA ...")
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("Calculating tSNE ...")
    sc.tl.tsne(adata_filt)
    print("Calculating neighbors ...")
    sc.pp.neighbors(adata_filt)
    print("Calculating Leiden ...")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("Calculating PAGA ...")
    sc.tl.paga(adata_filt)
    print("Plotting PAGA ...")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("Calculating UMAP init_pos = paga")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("Calculating UMAP ...")
    sc.tl.umap(adata_filt)
    print("Plotting UMAP ...")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def recalc_UMAP(adata_filt, leiden_res = 0.8):
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("Calculating tSNE ...")
    sc.tl.tsne(adata_filt)
    print("Calculating neighbors ...")
    sc.pp.neighbors(adata_filt)
    print("Calculating Leiden ...")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("Calculating PAGA ...")
    sc.tl.paga(adata_filt)
    print("Plotting PAGA ...")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("Calculating UMAP init_pos = paga")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("Calculating UMAP ...")
    sc.tl.umap(adata_filt)
    print("Plotting UMAP ...")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def process_until_norm(adata, cells):
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = cells) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 10] # Filter on the cells with fewer than 10% mitochondrial reads
    sc.pp.normalize_total(adata_filt, target_sum = 1e4) # Normalize
    sc.pp.log1p(adata_filt) # Log scaling
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    return adata_filt
    

def isolate_cells_by_gene(data, gene, threshold):
    # Now subset_ant_mt_filt contains only the highly variable genes
    data_subset = data[data[:, gene].X > threshold]
    
    return data_subset

#%% Environment settings and misc variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)
plot_LGRMK = ['LGR5', 'MKI67', 'leiden', 'Localization']
global_res = 0.5


#%% Read the files
col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")
ant_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_antrum.h5ad")
duo_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_duodenum.h5ad")
col_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_colon.h5ad")

#%% Tagging and labeling
combined = ad.concat([ant_unfilt, duo_unfilt, col_unfilt], join = 'outer')
combined_no_col = ad.concat([ant_unfilt, duo_unfilt], join = 'outer')
combined.obs['Localization'] = combined.obs['Site'].astype(str) + ' ' + combined.obs['Patient'].astype(str)
combined_no_col.obs['Localization'] = combined_no_col.obs['Site'].astype(str) + ' ' + combined_no_col.obs['Patient'].astype(str)

#%% Processing the UMAP
combined_processed = process_for_UMAP(combined, leiden_res = global_res)
restrict_process = process_for_UMAP(combined_no_col, leiden_res = global_res)
ant_pat_processed = process_for_UMAP(ant_unfilt, leiden_res = global_res)

#%% picking out the proliferative cells and LGR5 stem cells then recalculating the UMAP

# Separating cells
combined_processed_LGR5 = isolate_cells_by_gene(data = combined_processed, gene = 'LGR5', threshold = 0.5)
combined_processed_MKI67 = isolate_cells_by_gene(data = combined_processed, gene = 'MKI67', threshold = 1.5)

# Separating cells no colon
restrict_process_LGR5 = isolate_cells_by_gene(data = restrict_process, gene = 'LGR5', threshold = 0.5)
restrict_process_MKI67 = isolate_cells_by_gene(data = restrict_process, gene = 'MKI67', threshold = 1.5)

# Recalculating the UMAP
reprocessed_LGR5 = recalc_UMAP(combined_processed_LGR5, leiden_res = 0.1)
reprocessed_MKI67 = recalc_UMAP(combined_processed_MKI67, leiden_res = 0.1)

# Recalculating the UMAP no colon
restrict_process_LGR5 = recalc_UMAP(restrict_process_LGR5, leiden_res = 0.1)
restrict_process_MKI67 = recalc_UMAP(restrict_process_MKI67, leiden_res = 0.1)

#%% Plotting

sc.pl.umap(combined_processed, color = ['Localization'], size = 0.9)
sc.pl.umap(combined_processed, color = ['Site'])
sc.pl.umap(combined_processed, color = ['Patient'])
sc.pl.umap(combined_processed, color = ['EPCAM', 'KRT8', 'CD3D'])
sc.pl.umap(combined_processed, color = ['LGR5', 'MKI67', 'KRT8', 'EPCAM', 'CD3D', 'CD19', 'NCAM1'])

sc.pl.umap(reprocessed_LGR5, color = ['Localization'], title = 'LGR5 Localization')
sc.pl.umap(reprocessed_LGR5, color = ['LGR5', 'MKI67', 'Localization'])
sc.pl.umap(reprocessed_MKI67, color = ['LGR5', 'MKI67', 'Localization'])

sc.pl.umap(reprocessed_LGR5, color = 'leiden', size = 70)

sc.pl.umap(reprocessed_MKI67, color = ['Localization'], title = 'MKI67 Localization', size = 40)
sc.pl.umap(reprocessed_MKI67, color = ['EPCAM', 'KRT8', 'CD3D', 'CD19', 'NCAM1'])
sc.pl.umap(reprocessed_MKI67, color = ['MKI67'])



subset_dat = reprocessed_LGR5[reprocessed_LGR5.obs['Localization'] == 'Antrum P26' ]
sc.pl.umap(subset_dat, color = plot_LGRMK)

#%%
sc.pl.umap(reprocessed_MKI67, color = plot_LGRMK, size = 70)
sc.pl.umap(reprocessed_LGR5, color = plot_LGRMK, size = 70)
sc.pl.umap(restrict_process_MKI67, color = plot_LGRMK, size = 70)
sc.pl.umap(restrict_process_LGR5, color = plot_LGRMK, size = 70)

sc.pl.umap(reprocessed_MKI67, color = 'Localization', size = 70)
sc.pl.umap(reprocessed_LGR5, color = 'Localization', size = 70)
sc.pl.umap(restrict_process_MKI67, color = ['leiden', 'Localization'], size = 70)
sc.pl.umap(restrict_process_LGR5, color = ['leiden', 'Localization'], size = 70)

sc.pl.umap(restrict_process, color = 'Localization', size = 4)


#%% Checking the LGR5 proliferation in the patient antrum alone
ant_pat = ant_unfilt[ant_unfilt.obs['Patient'] == 'P26']
ant_pat_processed = process_for_UMAP(ant_pat)
sc.pl.umap(ant_pat_processed, color = ['LGR5', 'MKI67', 'leiden'])

#%%
temp_test = restrict_process_MKI67[restrict_process_MKI67.obs['Localization'] == 'Antrum P26']

sc.pl.scatter(restrict_process_MKI67[restrict_process_MKI67.obs['Patient'] == 'P26'], y = 'MKI67', x = 'leiden', color = 'Localization')
sc.pl.scatter(restrict_process_MKI67, y = 'MKI67', x = 'leiden', color = 'Localization')

sc.pl.violin(restrict_process_MKI67, keys = 'MKI67', groupby = 'Localization', rotation = 45)
sc.pl.violin(temp_test, keys = 'MKI67', groupby = 'leiden', rotation = 0)

#%%

sc.pl.umap(ant_pat_processed, color = ['LGR5', 'MKI67'])
