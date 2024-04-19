# This script is meant for the processing of different files into a format that's ready for filtering by expression
# and downstream analyses

# Import packages
import scanpy as sc
import anndata as ad
import pandas as pd

#%% function definitions
# DEFAULT QC VALUES. Calibrated to Sarah Teichmann's paper "Cells of the human intestinal tract mapped across space and time." These QC values will apply by default for this entire script.
def filter_cells_for_UMAP(data, min_ct = 2000, min_gen = 200, min_cell = 3, mt_pct = 50, max_genes = 8000, normed = 0): 
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if max_genes > 0:
        adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
        
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filtering based on percentage of mitochondrial genes
    else:
        adata_filt = adata_prefilt
    return adata_filt    

def process_for_UMAP(data, normed = 0, leiden_res = 0.8, filtering = 1, min_ct = 2000, min_gen = 200, min_cell = 3, mt_pct = 50, max_genes = 8000): # DEFAULT QC VALUES
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


def process_until_norm(data, cells, min_ct = 2000, min_gen = 200, min_cell = 3, mt_pct = 50, max_genes = 8000): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    print("################# Filtering ... #################")
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if max_genes > 0:
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

def map_to_column(data, map_set, column = 'Localization'):
    data.obs[column + '_old'] = data.obs[column]
    data.obs[column] = data.obs[column].map(map_set)
    print(data.obs[column])
    print(data.obs[column + '_old'])
    return 'Mapping function done'

#%% Environment settings and misc variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

# MIK67 = Ki67, TNSFRSF19 = TROY
inspect_stem = ['LGR5', 'MKI67', 'TNFRSF19', 'BMI1', 'LRIG1', 'AGR2', 'HSPA5', 'leiden', 'Localization']
global_res = 0.5
#%%

ant_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_antrum.h5ad")

patient_ant = ant_unfilt[ant_unfilt.obs['Patient'] == 'P26']

patient_ant_filt = process_for_UMAP(patient_ant)
patient_epith = filter_clusters_by_gene(data = patient_ant_filt, gene = 'EPCAM', threshold = 0.5)
patient_epith_barcodes = patient_epith.obs_names.tolist()
gated_unproc = ant_unfilt[patient_epith_barcodes]
gated_proc = process_for_UMAP(data = gated_unproc, leiden_res = 0.5)

#%%
sc.pl.umap(gated_proc, color = ['TPSAB1', 'MUC5AC', 'MUC6', 'MUC2', 'ANPEP', 'LGR5', 'MKI67', 'AGR2', 'HSPA5', 'GHRL', 'leiden'])
sc.pl.umap(patient_ant_filt, color = ['EPCAM', 'THY1', 'leiden'])
sc.tl.rank_genes_groups(adata = gated_proc, groupby = 'leiden', method = 'wilcoxon')
sc.pl.rank_genes_groups(gated_proc, n_genes = 25)

clusters_mapping = {
    '0' : 'Eterocytes',
    '1' : 'GHRL enteroendocrine',
    '2' : 'LGR5+ stem cells',
    '3' : '3',
    '4' : 'Prolif 1',
    '5' : 'ECL cells',
    '6' : 'Gastric SMC and GMC',
    '7' : '7',
    '8' : '8'
        }

map_to_column(data = gated_proc, map_set = clusters_mapping, column = 'leiden')

#%%

sc.pl.umap(gated_proc, color = ['SERPINH1', 'MUC6', 'ARF4', 'leiden'])
