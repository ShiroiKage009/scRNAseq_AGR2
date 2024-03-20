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
# col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")
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

sc.pl.umap(adata = combined_processed, color = 'leiden', legend_loc = 'on data')

#%% Picking out epithelial leiden clusters. The clusters were identified manually.

ep1 = combined_processed[combined_processed.obs['leiden'] == '27']
ep2 = combined_processed[combined_processed.obs['leiden'] == '26']
ep3 = combined_processed[combined_processed.obs['leiden'] == '21']
ep4 = combined_processed[combined_processed.obs['leiden'] == '0']
ep5 = combined_processed[combined_processed.obs['leiden'] == '7']
ep6 = combined_processed[combined_processed.obs['leiden'] == '29']
ep7 = combined_processed[combined_processed.obs['leiden'] == '17']
ep8 = combined_processed[combined_processed.obs['leiden'] == '18']
ep9 = combined_processed[combined_processed.obs['leiden'] == '24']
ep10 = combined_processed[combined_processed.obs['leiden'] == '14']
ep11 = combined_processed[combined_processed.obs['leiden'] == '16']
ep12 = combined_processed[combined_processed.obs['leiden'] == '25']
ep13 = combined_processed[combined_processed.obs['leiden'] == '5']
ep14 = combined_processed[combined_processed.obs['leiden'] == '8']
ep15 = combined_processed[combined_processed.obs['leiden'] == '1']
ep16 = combined_processed[combined_processed.obs['leiden'] == '9']
ep17 = combined_processed[combined_processed.obs['leiden'] == '15']

epithelial_filt = ad.concat(adatas = [ep1, ep2, ep3, ep4, ep5, ep6, ep7, ep8, ep9, ep10, ep11, ep12, ep13, ep14, ep15, ep16, ep17], join = 'outer')
sc.pl.umap(epithelial_filt, color = ['LGR5', 'MKI67', 'EPCAM', 'KRT8', 'leiden'])
epithelial_filt_proc = recalc_UMAP(epithelial_filt, leiden_res = global_res)
epithelial_filt_proc_2 = recalc_UMAP(epithelial_filt_proc, leiden_res = 0.1)

sc.pl.umap(epithelial_filt_proc, color = ['LGR5', 'MKI67', 'EPCAM', 'KRT8', 'leiden', 'Site'])

test_ant = combined_processed[combined_processed.obs['Site'] == 'Antrum']
test_ant.obs['Localization'] = test_ant.obs['Site'].astype(str) + ' ' + test_ant.obs['Patient'].astype(str)
antLGR5 = isolate_cells_by_gene(data=test_ant, gene = 'LGR5', threshold = 0.5)
sc.pl.umap(recalc_UMAP(antLGR5, leiden_res = 0.5), color = ['LGR5', 'Site', 'EPCAM', 'Localization'])
test_ant_2 = test_ant[test_ant.obs['Localization'] == 'Antrum P26']


test_iso_2 = isolate_cells_by_gene(data = epithelial_filt_proc, gene = 'LGR5', threshold = 0.1)
test_iso_2.obs['Localization'] = test_iso_2.obs['Site'].astype(str) + ' ' + test_iso_2.obs['Patient'].astype(str)

sc.pl.umap(adata = test_iso_2, color = ['LGR5', 'Site', 'leiden', 'Localization'])



#%% picking out the proliferative cells and LGR5 stem cells then recalculating the UMAP

# Separating cells
combined_processed_LGR5 = isolate_cells_by_gene(data = epithelial_filt_proc, gene = 'LGR5', threshold = 0.5)
combined_processed_MKI67 = isolate_cells_by_gene(data = combined_processed, gene = 'MKI67', threshold = 1.5)

# Separating cells no colon
restrict_process_LGR5 = isolate_cells_by_gene(data = restrict_process, gene = 'LGR5', threshold = 0.5)
restrict_process_MKI67 = isolate_cells_by_gene(data = restrict_process, gene = 'MKI67', threshold = 1.5)

# Recalculating the UMAP
reprocessed_LGR5 = recalc_UMAP(combined_processed_LGR5, leiden_res = global_res)
reprocessed_MKI67 = recalc_UMAP(combined_processed_MKI67, leiden_res = global_res)

# Recalculating the UMAP no colon
restrict_process_LGR5 = recalc_UMAP(restrict_process_LGR5, leiden_res = global_res)
restrict_process_MKI67 = recalc_UMAP(restrict_process_MKI67, leiden_res = global_res)

#%%

sc.pl.umap(epithelial_filt_proc, color = 'leiden')
epith_recalc_test = recalc_UMAP(adata_filt = epithelial_filt_proc, leiden_res = 0.01)
sc.pl.umap(adata = epith_recalc_test, color = 'leiden')
