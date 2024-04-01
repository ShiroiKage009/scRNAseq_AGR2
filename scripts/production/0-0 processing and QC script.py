# This script is meant for the processing of different files into a format that's ready for filtering by expression
# and downstream analyses

# Import packages
import scanpy as sc
import anndata as ad
import pandas as pd

#%% function definitions
def filter_cells_for_UMAP(data, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 50, max_genes = 8000, normed = 0): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filtering based on percentage of mitochondrial genes
    else:
        adata_filt = adata_prefilt
    return adata_filt    

def process_for_UMAP(data, normed = 0, leiden_res = 0.8, filtering = 1, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 50, max_genes = 8000): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    if filtering:
        adata_filt = filter_cells_for_UMAP(data = adata, min_ct = min_ct, min_gen = min_gen, min_cell = min_cell, max_genes = max_genes, mt_pct = mt_pct)
    else:
        adata_filt = adata       
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
    sc.tl.tsne(adata_filt) # Calculate tsne
    print("################# Calculating neighbors ... #################")
    sc.pp.neighbors(adata_filt) # Calculate neighbors
    print("################# Calculating Leiden ... #################")
    sc.tl.leiden(adata_filt, resolution = leiden_res) # Calculate Leiden clusters
    print("################# Calculating PAGA ... #################")
    sc.tl.paga(adata_filt) # Calculate PAGA
    print("################# Plotting PAGA ... #################")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("################# Calculating UMAP init_pos = paga #################")
    sc.tl.umap(adata_filt, init_pos='paga') # Plot PAGA
    print("################# Calculating UMAP ... #################")
    sc.tl.umap(adata_filt) # Calculate UMAP
    print("#################Plotting UMAP ... #################")
    sc.pl.umap(adata_filt, color = ['leiden']) # Plot UMAP and show Leiden clusters
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################


def recalc_UMAP(data_filt, leiden_res = 0.8):
    adata_filt = data_filt
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("################# Calculating tSNE ... #################")
    sc.tl.tsne(adata_filt) # Calculate tsne
    print("################# Calculating neighbors ... #################")
    sc.pp.neighbors(adata_filt) # Calculate neighbors
    print("################# Calculating Leiden ... #################")
    sc.tl.leiden(adata_filt, resolution = leiden_res) # Calculate Leiden clusters
    print("################# Calculating PAGA ... #################")
    sc.tl.paga(adata_filt) # Calculate PAGA
    print("################# Plotting PAGA ... #################")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("################# Calculating UMAP init_pos = paga#################")
    sc.tl.umap(adata_filt, init_pos='paga') # Calculate PAGA
    print("################# Calculating UMAP ... #################")
    sc.tl.umap(adata_filt) # Calculate UMAP
    print("################# Plotting UMAP ... #################")
    sc.pl.umap(adata_filt, color = ['leiden']) # Plot UMAP and show Leiden clusters
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################


def process_until_norm(data, cells, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 50, max_genes = 8000): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
    adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filter on the cells with fewer than 10% mitochondrial reads
    print("################# Normalizing ... #################")
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print("################# Log scaling ... #################")
    sc.pp.log1p(adata_filt) # Log scaling
    print("#################Finding variable genes ... #################")
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print("################# Saving raw data ... #################")
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################
    

def isolate_cells_by_gene(data, gene, threshold):
    # Now subset_ant_mt_filt contains only the highly variable genes
    data_subset = data[data[:, gene].X > threshold]
    return data_subset
#######################################################
################## FUNCTION DEF END ###################
#######################################################

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
        gene_presence = de_results[cluster].str.contains(gene)
        gene_expression = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][gene_presence]
        if any(gene_expression >= expression_threshold):
            significant_clusters.append(cluster)
    # Subset the data to include only cells from the significant clusters
    adata_subset = adata[adata.obs['leiden'].isin(significant_clusters)].copy()
    return adata_subset

#%% Environment settings and misc variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

# MIK67 = Ki67, TNSFRSF19 = TROY
inspect_stem = ['LGR5', 'MKI67', 'TNFRSF19', 'BMI1', 'LRIG1', 'leiden', 'Localization']
global_res = 0.5


#%% Read the files
col_org_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")
ant_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_antrum.h5ad")
duo_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_duodenum.h5ad")
col_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_colon.h5ad")

#%% Tagging and labeling
combined = ad.concat([ant_unfilt, duo_unfilt, col_unfilt], join = 'outer')
combined_unique = combined[~combined.obs['cellbarcode'].duplicated(keep = 'first')].copy()
combined_nocol = ad.concat([ant_unfilt, duo_unfilt], join = 'outer')
combined.obs['Localization'] = combined.obs['Site'].astype(str) + ' ' + combined.obs['Patient'].astype(str)
combined_nocol = ad.concat([ant_unfilt, duo_unfilt], join = 'outer')
combined_nocol.obs['Localization'] = combined_nocol.obs['Site'].astype(str) + ' ' + combined_nocol.obs['Patient'].astype(str)
ant_unfilt.obs['Localization'] = ant_unfilt.obs['Site'].astype(str) + ' ' + ant_unfilt.obs['Patient'].astype(str)
combined_control = combined_unique[combined_unique.obs['Patient'] == 'GI6253', :] # Using the unique object seems to work for filtering while usign the combined object directly doesn't. No idea why.
combined_control.obs['Localization'] = combined_control.obs['Site'].astype(str) + ' ' + combined_control.obs['Patient'].astype(str)


#%% Initial processing the UMAP
combined_proc = process_for_UMAP(combined, leiden_res = global_res)
combined_nocol_proc = process_for_UMAP(combined_nocol, leiden_res = global_res)
antrum_proc = process_for_UMAP(ant_unfilt, leiden_res = global_res)
combined_control_proc = process_for_UMAP(combined_control, leiden_res = global_res)

#%% Write files to save then load in the next script
combined_proc.write_h5ad(filename = 'C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/combined_proc.h5ad')
combined_nocol_proc.write_h5ad(filename = 'C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/combined_nocol_proc.h5ad')
antrum_proc.write_h5ad(filename = 'C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/antrum_proc.h5ad')
combined_control_proc.write_h5ad(filename = 'C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/combined_control_proc.h5ad')

#%% Printing the bookend
# =============================================================================
print("######################################################################")
print("####                        END OF SCRIPT                         ####")
print("######################################################################")
# =============================================================================

#%% Testing cell within the script


