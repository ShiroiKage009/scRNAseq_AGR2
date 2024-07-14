import scanpy as sc
import pandas as pd
import time

start_time = time.time()

# DEFAULT QC VALUES. Calibrated to Sarah Teichmann's paper 'Cells of the human intestinal tract mapped across space and time.' These QC values will apply by default for this entire script.
def filter_cells_for_UMAP(data, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 60, normed = 0): 
    adata = data # This is to avoid writing into the file that's entered as an argument
    print('################# Filtering ... #################')
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['doublet_scores'] < 0.24]
   # if max_genes > 0:
   #     adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
        
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filtering based on percentage of mitochondrial genes
    else:
        adata_filt = adata_prefilt
    return adata_filt    

def process_for_UMAP(data, normed = 0, leiden_res = 0.8, filtering = 1, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 60): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    if filtering:
        adata_filt = filter_cells_for_UMAP(data = adata, min_ct = min_ct, min_gen = min_gen, min_cell = min_cell, mt_pct = mt_pct)
    else:
        adata_filt = adata       
    print('################# Normalizing ... #################')
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print('################# Log scaling ... #################')
    sc.pp.log1p(adata_filt) # Log scaling
    print('################# Finding variable genes ... #################')
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print('################# Saving raw data ... #################')
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    #print('################# Filtering on variable genes ... #################')
    #adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    print('################# Regressing ... #################')
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print('################# Scaling ... #################')
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print('################# Calculating tSNE ... #################')
    sc.tl.tsne(adata_filt) # Calculate tsne
    print('################# Calculating neighbors ... #################')
    sc.pp.neighbors(adata_filt, n_neighbors = 15) # Calculate neighbors
    print('################# Calculating Leiden ... #################')
    sc.tl.leiden(adata_filt, resolution = leiden_res) # Calculate Leiden clusters
    print('################# Calculating PAGA ... #################')
    sc.tl.paga(adata_filt) # Calculate PAGA
    print('################# Plotting PAGA ... #################')
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print('################# Calculating UMAP init_pos = paga #################')
    sc.tl.umap(adata_filt, init_pos='paga') # Plot PAGA
    print('################# Calculating UMAP ... #################')
    sc.tl.umap(adata_filt) # Calculate UMAP
    print('#################Plotting UMAP ... #################')
    sc.pl.umap(adata_filt, color = ['leiden']) # Plot UMAP and show Leiden clusters
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################

# diff is boolean specifying if differential expresion is already done
# threshold is the threshold of expression
def filter_clusters_by_gene(data, gene, threshold = 0.5):
    # Load your AnnData object
    adata = data
    # Extract the DataFrame for the differential expression results
    de_results = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    # Define a threshold for significant expression (adjust as needed)
    expression_threshold = threshold
    # Find clusters with significant gene expression
    significant_clusters = []
    for cluster in de_results.columns:
        gene_presence = de_results[cluster].str.contains(gene)
        gene_expression = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][gene_presence]
        if any(gene_expression >= expression_threshold):
            significant_clusters.append(cluster)
    # Subset the data to include only cells from the significant clusters
    adata_subset = adata[adata.obs['leiden'].isin(significant_clusters)].copy()
    return adata_subset
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

thresholded = 0.45

ann_antrum = sc.read_h5ad('C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/Data annotation/antrum/antrum_annotated_final1122.h5ad')
ant_unfilt = sc.read_h5ad('C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_antrum.h5ad')
ant_filt = filter_cells_for_UMAP(ant_unfilt)
ant_proc = process_for_UMAP(ant_unfilt)
sc.tl.rank_genes_groups(adata = ant_proc, groupby = 'leiden', metohd = 't-test')
ant_epith_clust = filter_clusters_by_gene(data = ant_proc, gene = 'EPCAM', threshold = thresholded)
ant_epith_cells = isolate_cells_by_gene(data= ant_proc, gene = 'EPCAM', threshold = thresholded)

end_time = time.time()
print ("Finished in", end_time - start_time, "seconds")
print ("Finished in", (end_time - start_time)/60, "minutes")

#%%
def filter_clusters_by_gene2(data, gene, threshold = 0.5):
    # Load your AnnData object
    adata = data
    # Extract the DataFrame for the differential expression results
    de_results = pd.DataFrame(adata.uns['rank_genes_groups']['names'])
    # Define a threshold for significant expression (adjust as needed)
    expression_threshold = threshold
    # Find clusters with significant gene expression
    significant_clusters = []
    for cluster in de_results.columns:
        gene_presence = de_results[cluster].str.contains(gene)
        gene_expression = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][gene_presence]
        if any(gene_expression >= expression_threshold):
            significant_clusters.append(cluster)
    # Subset the data to include only cells from the significant clusters
    adata_subset = adata[adata.obs['leiden'].isin(significant_clusters)].copy()
    return adata_subset

filted = filter_clusters_by_gene2(data = ant_proc, gene = 'EPCAM', threshold = 0.05)

#%%
print("Tom's thing")
ann_antrum.obs['SampleID'].value_counts()
#%%
print("Epithelial cluster isolation at threshold", thresholded)
ant_epith_clust.obs['Patient'].value_counts()
#%%
print("Epithelial cluster isolation at threshold", thresholded)
filted.obs['Patient'].value_counts()
