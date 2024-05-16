# testing published datasets

# Import packages
import scanpy as sc
import pandas as pd
import time

# Define functions
def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper

# DEFAULT QC VALUES. Calibrated to Sarah Teichmann's paper 'Cells of the human intestinal tract mapped across space and time.' These QC values will apply by default for this entire script.
def filter_cells_for_UMAP(data, min_gen = 500, min_cell = 3, mt_pct = 60, normed = 0, doublet_threshold = 0.24): #min_ct is removed. Consider restoring it later. You will find it commented out below.
    adata = data # This is to avoid writing into the file that's entered as an argument
    print('################# Filtering ... #################')
    #sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['doublet_scores'] < doublet_threshold]
   # if max_genes > 0:
   #     adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
        
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filtering based on percentage of mitochondrial genes
    else:
        adata_filt = adata_prefilt
    return adata_filt    

@time_it
def process_for_UMAP(data, normed = 0, leiden_res = 0.8, filtering = 1, min_gen = 500, min_cell = 3, mt_pct = 60, neighbors = 15): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    if filtering:
        adata_filt = filter_cells_for_UMAP(data = adata, min_gen = min_gen, min_cell = min_cell, mt_pct = mt_pct) #min ct is removed. Think about restoring it later.
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
    print('################# Filtering on variable genes ... #################')
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    #print('################# Regressing ... #################')
    #sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print('################# Scaling ... #################')
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    #print('################# Calculating tSNE ... #################')
    #sc.tl.tsne(adata_filt) # Calculate tsne
    print('################# Calculating neighbors ... #################')
    sc.pp.neighbors(adata_filt, n_neighbors = neighbors) # Calculate neighbors
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

@time_it
def prep_data(data, filtering = 1, min_ct = 2000, min_gen = 500, min_cell = 3, mt_pct = 60): # DEFAULT QC VALUES
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
    #print('################# Regressing ... #################')
    #sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print('################# Scaling ... #################')
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    #print('################# Calculating tSNE ... #################')
    #sc.tl.tsne(adata_filt) # Calculate tsne
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################

@time_it
def recalc_UMAP(data_filt, leiden_res = 0.8):
    adata_filt = data_filt
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print('################# Calculating tSNE ... #################')
    sc.tl.tsne(adata_filt) # Calculate tsne
    print('################# Calculating neighbors ... #################')
    sc.pp.neighbors(adata_filt) # Calculate neighbors
    print('################# Calculating Leiden ... #################')
    sc.tl.leiden(adata_filt, resolution = leiden_res) # Calculate Leiden clusters
    print('################# Calculating PAGA ... #################')
    sc.tl.paga(adata_filt) # Calculate PAGA
    print('################# Plotting PAGA ... #################')
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print('################# Calculating UMAP init_pos = paga#################')
    sc.tl.umap(adata_filt, init_pos='paga') # Calculate PAGA
    print('################# Calculating UMAP ... #################')
    sc.tl.umap(adata_filt) # Calculate UMAP
    print('################# Plotting UMAP ... #################')
    sc.pl.umap(adata_filt, color = ['leiden']) # Plot UMAP and show Leiden clusters
    return adata_filt
#######################################################
################## FUNCTION DEF END ###################
#######################################################


@time_it
def process_until_norm(data, cells, min_ct = 2000, min_gen = 200, min_cell = 3, mt_pct = 50, max_genes = 8000): # DEFAULT QC VALUES
    adata = data # This is to avoid writing into the file that's entered as an argument
    print('################# Filtering ... #################')
    sc.pp.filter_cells(adata, min_counts = min_ct) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= min_gen) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = min_cell) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if max_genes > 0:
        adata_prefilt = adata_prefilt[adata_prefilt.obs['n_genes_by_counts'] < max_genes]
    adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < mt_pct] # Filter on the cells with fewer than 10% mitochondrial reads
    print('################# Normalizing ... #################')
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print('################# Log scaling ... #################')
    sc.pp.log1p(adata_filt) # Log scaling
    print('#################Finding variable genes ... #################')
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print('################# Saving raw data ... #################')
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
def filter_clusters_by_gene(data, gene, threshold = 0.5, method = None):
    # Load your AnnData object
    adata = data
    print('################# Ranked gene expression ... #################')
    sc.tl.rank_genes_groups(adata, groupby='cell_type', method = method)
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
    adata_subset = adata[adata.obs['cell_type'].isin(significant_clusters)].copy()
    return adata_subset
#######################################################
################## FUNCTION DEF END ###################
#######################################################

#break Execution start begin with timing
start_time = time.time()
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

# MIK67 = Ki67, TNSFRSF19 = TROY
inspect_stem_symbol = ['cell_type', 'MKI67', 'LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2', 'Detailed_Cell_Type', 'Patient status']
inspect_stem_ensembl = ['cell_type', 'ENSG00000148773', 'ENSG00000139292', 'ENSG00000183734', 'ENSG00000026508', 'ENSG00000112562', 'ENSG00000117632', 'ENSG00000102837', 'ENSG00000168283', 'ENSG00000144749', 'ENSG00000180535', 'ENSG00000127863', 'ENSG00000181449', 'Detailed_Cell_Type', 'Patient_status']

inspect_stem_duo_symbol = ['cell_type', 'MKI67', 'LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'Detailed_Cell_Type']
inspect_stem_duo_ensembl = ['cell_type', 'ENSG00000148773', 'LGENSG00000139292R5', 'ENSG00000183734', 'ENSG00000026508', 'ENSG00000112562', 'ENSG00000117632', 'ENSG00000102837', 'ENSG00000168283', 'ENSG00000144749', 'Detailed_Cell_Type']

inspect_stem_gast_symbol = ['cell_type', 'MKI67', 'LGR5', 'TNFRSF19', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'Detailed_Cell_Type']
inspect_stem_gast_ensembl = ['cell_type', 'ENSG00000148773', 'ENSG00000139292', 'ENSG00000127863', 'ENSG00000183734', 'ENSG00000026508', 'ENSG00000112562', 'ENSG00000181449', 'ENSG00000117632', 'ENSG00000102837', 'ENSG00000168283', 'ENSG00000144749', 'Detailed_Cell_Type']


gastric_markers_symbol = ['cell_type', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1', 'Detailed_Cell_Type']
gastric_markers_ensembl = ['cell_type', 'ENSG00000148773', 'ENSG00000139292', 'ENSG00000127863', 'ENSG00000180535', 'ENSG00000183734', 'ENSG00000026508', 'ENSG00000112562', 'ENSG00000181449', 'ENSG00000117632', 'ENSG00000102837', 'ENSG00000168283', 'ENSG00000144749', 'ENSG00000184956', 'ENSG00000160181', 'ENSG00000215182', 'ENSG00000183607', 'ENSG00000160182', 'ENSG00000157017', 'ENSG00000161798', 'ENSG00000185499', 'Detailed_Cell_Type']


duodenal_markers = ['cell_type', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'LYZ', 'Detailed_Cell_Type']
duodenal_markers = ['cell_type', 'ENSG00000148773', 'ENSG00000139292', 'ENSG00000127863', 'ENSG00000180535', 'ENSG00000183734', 'ENSG00000026508', 'ENSG00000112562', 'ENSG00000181449', 'ENSG00000117632', 'ENSG00000102837', 'ENSG00000168283', 'ENSG00000144749', 'ENSG00000184956', 'ENSG00000198788', 'ENSG00000160180', 'ENSG00000157017', 'ENSG00000161798', 'ENSG00000185499', 'ENSG00000166825', 'ENSG00000090382', 'Detailed_Cell_Type']

global_res = 0.5

unif_data = sc.read_h5ad('C:/Work cache/data cache/testing/gim and eim unified umap of all file.h5ad')
unif_data.obsm['X_umap'] = unif_data.obsm['X_umap_MinDist_0.5_N_Neighbors_15']
sc.pl.umap(unif_data, color = 'Detailed_Cell_Type')
sc.pl.umap(unif_data, color = 'cell_type')

unif_GIM_healthy = unif_data[unif_data.obs['Patient_status'] != 'BE']

unif_GIM_healthy_prc = process_for_UMAP(unif_GIM_healthy, filtering = 0) # tSNE appears to be very hard to calculate for this dataset so I commented it out
GIM = unif_data[unif_data.obs['Patient_status'] == 'IM']
GIM_proc = process_for_UMAP(GIM, filtering = 0)
healthy = unif_data[unif_data.obs['Patient_status'] == 'Healthy']
healthy_proc = process_for_UMAP(healthy, filtering = 0)

GIM_epith = filter_clusters_by_gene(data = GIM_proc, gene = 'ENSG00000119888', method = 'wilcoxon')
healthy_epith = filter_clusters_by_gene(data = healthy_proc, gene = 'ENSG00000119888', method = 'wilcoxon')

GIM_epith_barcodes = GIM_epith.obs_names.tolist()
healthy_epith_barcodes = healthy_epith.obs_names.tolist()

GIM_epith_refilt = unif_data[GIM_epith_barcodes]
healthy_epith_refilt = unif_data[healthy_epith_barcodes]

GIM_epith_refilt_prepped = prep_data(GIM_epith_refilt, filtering = 0)
'ENSG00000139292' in GIM_epith_refilt_prepped.var_names

end_time = time.time()
print('done in', end_time - start_time, 'seconds')
print('done in', (end_time - start_time)/60, 'minutes')
#%%

GIM_LGR5 = isolate_cells_by_gene(data = GIM_epith_refilt, gene = 'ENSG00000139292', threshold = 0.1)
GIM_LGR5_barcodes = GIM_LGR5.obs_names.tolist()
GIM_LGR5_refilt = unif_data[GIM_LGR5_barcodes]
GIM_LGR5_reproc = process_for_UMAP(GIM_LGR5_refilt, filtering = 0)

sc.pl.umap(GIM_LGR5_reproc, color = ['leiden', 'Sample'])


healthy_LGR5 = isolate_cells_by_gene(data = healthy_epith_refilt, gene = 'ENSG00000139292', threshold = 0.1)
healthy_LGR5_barcodes = healthy_LGR5.obs_names.tolist()
healthy_LGR5_refilt = unif_data[healthy_LGR5_barcodes]
healthy_LGR5_reproc = process_for_UMAP(healthy_LGR5_refilt, filtering = 0)

sc.pl.umap(healthy_LGR5_reproc, color = ['leiden', 'Sample'])

LGR5_combo = healthy_LGR5_barcodes + GIM_LGR5_barcodes
combo_LGR5_refilt = unif_data[LGR5_combo]
combo_LGR5_reproc = process_for_UMAP(data = combo_LGR5_refilt, filtering = 0)

sc.pl.umap(combo_LGR5_refilt, color = inspect_stem_ensembl, title = inspect_stem_symbol, size = 90)
sc.pl.scatter(combo_LGR5_refilt, x = 'ENSG00000139292', y = 'ENSG00000102837', color = 'ENSG00000112562')
sc.pl.scatter(combo_LGR5_refilt, x = 'ENSG00000139292', y = 'ENSG00000112562', color = 'ENSG00000102837')
sc.pl.scatter(combo_LGR5_refilt, x = 'ENSG00000139292', y = 'ENSG00000102837', color = 'Patient_status')

sc.pl.umap(combo_LGR5_reproc, color = inspect_stem_ensembl, title = inspect_stem_symbol, size = 90)
sc.pl.scatter(combo_LGR5_reproc, x = 'ENSG00000139292', y = 'ENSG00000102837', color = 'ENSG00000112562')
sc.pl.scatter(combo_LGR5_reproc, x = 'ENSG00000139292', y = 'ENSG00000112562', color = 'ENSG00000102837')
sc.pl.scatter(combo_LGR5_reproc, x = 'ENSG00000139292', y = 'ENSG00000102837', color = 'Tissue_in_paper')

sc.pl.violin(combo_LGR5_reproc, keys = ['ENSG00000139292', 'ENSG00000112562', 'ENSG00000102837'], groupby = 'Tissue_in_paper')

# The object includes intestinal tissue as well as gastric tissue. You must filter down to only the gastric tissues.

#%%
sc.pl.umap(GIM_epith, color = inspect_stem_ensembl, title = inspect_stem_symbol)
sc.pl.scatter(GIM_epith, x = 'ENSG00000112562', y = 'ENSG00000102837', color = 'ENSG00000139292', size = 150)
sc.pl.scatter(GIM_epith, x = 'ENSG00000112562', y = 'ENSG00000102837', color = 'Sample', size = 150)

GIM_LGR5 = isolate_cells_by_gene(data = GIM_epith, gene = 'ENSG00000139292', threshold = 0.1)

#%%
sc.tl.rank_genes_groups(adata = unif_data, groupby = 'Detailed_Cell_Type', method = 'wilcoxon')
sc.pl.rank_genes_groups(adata = unif_data, n_genes = 20)
