# This script performs some analysis using hte cells identified and filtered based on previous analyses but begins from their raw values.# This script is for poking around in the files filtered 
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
    print("################# Plotting UMAP ... #################")
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
def filter_clusters_by_gene(data, gene, threshold = 0.5, method = 't-test'):
    # Load your AnnData object
    adata = data
    sc.tl.rank_genes_groups(adata, groupby='leiden', method = method)
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


#%% Setting up environmental variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)
inspect_stem = ['LGR5', 'MKI67', 'TNFRSF19', 'BMI1', 'SMOC2', 'ASCL2', 'LRIG1', 'Patient']
global_res = 0.5
LGR5_threshold = 0.5
diff_exp_method = 'wilcoxon'

#%% Readin the files
# Writing the refiltered and processed files
combined_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/combined_unfilt_subset.h5ad')
antrum_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/antrum_unfilt_subset.h5ad')
combined_LGR5_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/combined_LGR5_unfilt_subset.h5ad')
antrum_LGR5_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/antrum_LGR5_unfilt_subset.h5ad')
nocol_LGR5_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/nocol_LGR5_unfilt_subset.h5ad')
nocol_refilt_subset = sc.read('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/refiltering from raw and reprocessing/nocol_unfilt_subset.h5ad')

#%%

processed_antrum = antrum_subset.copy()
processed_antrum = process_for_UMAP(processed_antrum, leiden_res = 0.5)
sc.pl.umap(processed_antrum, color = ['LGR5', 'ANPEP', 'MUC5AC', 'MUC2', 'MUC6', 'GKN2', 'MKI67', 'SMOC2', 'GHRL', 'TFF1', 'TFF2', 'GAST', 'CCKBR', 'CHGA', 'leiden', 'Patient'], size = 10)
sc.tl.rank_genes_groups(adata = processed_antrum, groupby = 'leiden')
sc.pl.rank_genes_groups(processed_antrum, n_genes = 25)

enterocytes = filter_clusters_by_gene(data = processed_antrum, gene = 'ANPEP', method = 'wilcoxon', threshold = 2)
LGR5 = filter_clusters_by_gene(data = processed_antrum, gene = 'LGR5', method = 'wilcoxon', threshold = 2)

sc.pl.umap(enterocytes, color = ['Patient', 'ANPEP', 'leiden'])
sc.pl.umap(LGR5, color = ['Patient', 'LGR5', 'leiden'])
sc.pl.umap(processed_antrum, color = inspect_stem)
sc.pl.scatter(adata = processed_antrum, x = 'LGR5', y = 'SMOC2', color = 'Patient')

entero_LGR5 = ad.concat(adatas = [LGR5, enterocytes], join = 'outer')

mapping_key = {
    '0' : '0',
    '1' : 'Foveolar2 (MUC5AClo/TFF1lo)',
    '2' : 'G Cells',
    '3' : 'Prolif 1',
    '4' : 'Foveolar3 (MUC5AC+/GKN2-)',
    '5' : 'Enterocytes',
    '6' : 'Base mucous gland cells (MUC6)',
    '7' : '7',
    '8' : 'Ghrelin Enteroendocrine cells',
    '9' : 'Foveolar (MUC5AC+/GKN2+)',
    '10' : '10',
    '11' : 'LGR5+ patient stem cells',
    '12' : 'LGR5+ control stem? cells',
    '13' : 'Prolif 2',
    '14' : 'ECL',
    '15' : '15',
    '16' : '16',
    '17' : '17',
    '18' : '18'
    }

map_to_column(data = entero_LGR5, map_set = mapping_key, column = 'leiden')
map_to_column(data = processed_antrum, map_set = mapping_key, column = 'leiden')
sc.pl.umap(entero_LGR5, color = 'leiden')

sc.pl.umap(entero_LGR5, color = ['LGR5', 'ANPEP', 'leiden', 'Patient', 'leiden_old'], size = 10)
sc.pp.highly_variable_genes(entero_LGR5, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
sc.tl.rank_genes_groups(adata = entero_LGR5, groupby = 'leiden', method = 'wilcoxon')
sc.pl.rank_genes_groups(entero_LGR5, n_genes = 25)

sc.pp.highly_variable_genes(processed_antrum, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
sc.tl.rank_genes_groups(adata = processed_antrum, groupby = 'leiden', method = 'wilcoxon')

#%%
sc.pl.umap(processed_antrum, color = ['SERPINH1', 'MUC6', 'ARF4', 'MUC5AC', 'TFE3', 'Patient', 'leiden'])

sc.pl.violin(adata = processed_antrum, keys = 'ARF4', groupby = 'Patient', rotation = 90)

#%%
# Extract the relevant arrays
gene_names = entero_LGR5.uns['rank_genes_groups']['names']
logfoldchanges = entero_LGR5.uns['rank_genes_groups']['logfoldchanges']
pvals_adj = entero_LGR5.uns['rank_genes_groups']['pvals_adj']


dataframes = {}
for group in ['Enterocytes', 'LGR5+ patient stem cells', 'LGR5+ control stem cells']:
    # Extract information for each group
    names = gene_names[group]
    lfc = logfoldchanges[group]
    pval_adj = pvals_adj[group]

    # Create DataFrame
    df = pd.DataFrame({
        'Gene Names': names,
        'Log Fold Change': lfc,
        'Adjusted P-Value': pval_adj
    }).set_index('Gene Names')

    dataframes[group] = df

# Now you have two DataFrames in the `dataframes` dict, one for each group
# Access them like this:
enterocytes_df = dataframes['Enterocytes']
pat_LGR5 = dataframes['LGR5+ patient stem cells']
cont_LGR5 = dataframes['LGR5+ control stem cells']

enterocytes_df.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/enterocytes_df.csv')
pat_LGR5.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/pat_LGR5.csv')
cont_LGR5.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/cont_LGR5.csv')

#%%
# Extract the relevant arrays
gene_names_2 = processed_antrum.uns['rank_genes_groups']['names']
logfoldchanges_2 = processed_antrum.uns['rank_genes_groups']['logfoldchanges']
pvals_adj_2 = processed_antrum.uns['rank_genes_groups']['pvals_adj']


dataframes_2 = {}
for group in ['Enterocytes', 'LGR5+ patient stem cells', 'LGR5+ control stem cells']:
    # Extract information for each group
    names = gene_names_2[group]
    lfc = logfoldchanges_2[group]
    pval_adj = pvals_adj_2[group]

    # Create DataFrame
    df = pd.DataFrame({
        'Gene Names': names,
        'Log Fold Change': lfc,
        'Adjusted P-Value': pval_adj
    }).set_index('Gene Names')

    dataframes_2[group] = df

# Now you have two DataFrames in the `dataframes_2` dict, one for each group
# Access them like this:
enterocytes_procant_df = dataframes_2['Enterocytes']
pat_LGR5_procant = dataframes_2['LGR5+ patient stem cells']
cont_LGR5_procant = dataframes_2['LGR5+ control stem cells']

enterocytes_procant_df.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/enterocytes_procant_df.csv')
pat_LGR5_procant.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/pat_LGR5_procant.csv')
cont_LGR5_procant.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/cont_LGR5_procant.csv')

