# Import packages
import scanpy as sc
import anndata as ad
import pandas as pd

#%% function definitions
def process_for_UMAP(data, normed = 0, leiden_res = 0.8):
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = 10) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 50] # Filter on the cells with fewer than 10% mitochondrial reads
    else:
        adata_filt = adata_prefilt
    print("################# Normalizing ... #################")
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print("################# Log scaling ... #################")
    sc.pp.log1p(adata_filt) # Log scaling
    print("################# Finding variable genes ... #################")
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print("################# Saving raw data ... #################")
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    print("################# Filtering on variable genes ... #################")
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    print("################# Regressing ... #################")
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print("################# Scaling ... #################")
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print("################# Calculating PCA ... #################")
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("################# Calculating tSNE ... #################")
    sc.tl.tsne(adata_filt)
    print("################# Calculating neighbors ... #################")
    sc.pp.neighbors(adata_filt)
    print("################# Calculating Leiden ... #################")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("################# Calculating PAGA ... #################")
    sc.tl.paga(adata_filt)
    print("################# Plotting PAGA ... #################")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("################# Calculating UMAP init_pos = paga #################")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("################# Calculating UMAP ... #################")
    sc.tl.umap(adata_filt)
    print("#################Plotting UMAP ... #################")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def recalc_UMAP(data_filt, leiden_res = 0.8):
    adata_filt = data_filt
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("################# Calculating tSNE ... #################")
    sc.tl.tsne(adata_filt)
    print("################# Calculating neighbors ... #################")
    sc.pp.neighbors(adata_filt)
    print("################# Calculating Leiden ... #################")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("################# Calculating PAGA ... #################")
    sc.tl.paga(adata_filt)
    print("################# Plotting PAGA ... #################")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("################# Calculating UMAP init_pos = paga#################")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("################# Calculating UMAP ... #################")
    sc.tl.umap(adata_filt)
    print("################# Plotting UMAP ... #################")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def process_until_norm(data, cells):
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = 10) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 50] # Filter on the cells with fewer than 10% mitochondrial reads
    print("################# Normalizing ... #################")
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print("################# Log scaling ... #################")
    sc.pp.log1p(adata_filt) # Log scaling
    print("#################Finding variable genes ... #################")
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print("################# Saving raw data ... #################")
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    return adata_filt
    

def isolate_cells_by_gene(data, gene, threshold):
    # Now subset_ant_mt_filt contains only the highly variable genes
    data_subset = data[data[:, gene].X > threshold]
    
    return data_subset

# This function filters the leiden clusters that are positivefor the gene you specify
# It assumes that you already did the differential expression analysis. 
# diff is boolean specifying if differential expresion is already done
# threshold is the threshold of expression
def filter_clusters_by_gene(data, gene, threshold = 0.5):
    # Load your AnnData object
    adata = data
    sc.tl.rank_genes_groups(adata, groupby='leiden')
    # Extract the DataFrame for the differential expression results
    de_results = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    # Define a threshold for significant expression (adjust as needed)
    expression_threshold = threshold
    # Find clusters with significant EPCAM expression
    significant_clusters = []
    for cluster in de_results.columns:
        epcam_gene = de_results[cluster].str.contains('EPCAM')
        epcam_expression = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][epcam_gene]
        if any(epcam_expression >= expression_threshold):
            significant_clusters.append(cluster)
    # Subset the data to include only cells from the significant clusters
    adata_subset = adata[adata.obs['leiden'].isin(significant_clusters)].copy()
    return adata_subset

#%% Environment settings and misc variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

inspect_stem = ['LGR5', 'MKI67', 'leiden', 'Localization']
global_res = 0.5


#%% Read the files
col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")
ant_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_antrum.h5ad")
duo_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_duodenum.h5ad")
col_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2_unfilt_colon.h5ad")

#%% Tagging and labeling
combined = ad.concat([ant_unfilt, duo_unfilt, col_unfilt], join = 'outer')
combined_nocol = ad.concat([ant_unfilt, duo_unfilt], join = 'outer')
combined.obs['Localization'] = combined.obs['Site'].astype(str) + ' ' + combined.obs['Patient'].astype(str)
combined_nocol = ad.concat([ant_unfilt, duo_unfilt], join = 'outer')
combined_nocol.obs['Localization'] = combined_nocol.obs['Site'].astype(str) + ' ' + combined_nocol.obs['Patient'].astype(str)
ant_unfilt.obs['Localization'] = ant_unfilt.obs['Site'].astype(str) + ' ' + ant_unfilt.obs['Patient'].astype(str)

#%% Initial processing the UMAP
combined_proc = process_for_UMAP(combined, leiden_res = global_res)
combined_nocol_proc = process_for_UMAP(combined_nocol, leiden_res = global_res)
antrum_proc = process_for_UMAP(ant_unfilt, leiden_res = global_res)

#%% Write files to save then load in the next script
combined_proc.write_h5ad(filename = 'C:/anaconda_envs/spyder/project data cache/testing integration with separation and the stem cells part 2/saved files/combined_proc.h5ad')
combined_nocol_proc.write_h5ad(filename = 'C:/anaconda_envs/spyder/project data cache/testing integration with separation and the stem cells part 2/saved files/combined_nocol_proc.h5ad')
antrum_proc.write_h5ad(filename = 'C:/anaconda_envs/spyder/project data cache/testing integration with separation and the stem cells part 2/saved files/antrum_proc.h5ad')

#%% Gating on epithelial cells in the fully-combined file
combined_epithelium = filter_clusters_by_gene(data = combined_proc, gene = 'EPCAM')
antrum_epithelium = filter_clusters_by_gene(data = antrum_proc, gene = 'EPCAM')
nocol_epithelium = filter_clusters_by_gene(data = combined_nocol_proc, gene = 'EPCAM')

#%%

combined_ep_LGR5 = isolate_cells_by_gene(data = combined_epithelium, gene = 'LGR5', threshold = 0.5)
combined_ep_MKI67 = isolate_cells_by_gene(data = combined_epithelium, gene = 'MKI67', threshold = 0.5)

antrum_ep_LGR5 = isolate_cells_by_gene(data = antrum_epithelium, gene = 'LGR5', threshold = 0.5)
antrum_ep_MKI67 = isolate_cells_by_gene(data = antrum_epithelium, gene = 'MKI67', threshold = 0.5)

nocol_ep_LGR5 = isolate_cells_by_gene(data = nocol_epithelium, gene = 'LGR5', threshold = 0.5)
nocol_ep_MKI67 = isolate_cells_by_gene(data = nocol_epithelium, gene = 'MKI67', threshold = 0.5)


#%% Test cell


