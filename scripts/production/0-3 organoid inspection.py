# Testing of annotated organoid data
import scanpy as sc
import pandas as pd
import time

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper

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

@time_it
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
    print('################# Filtering on variable genes ... #################')
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
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
    print('################# Filtering on variable genes ... #################')
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    print('################# Regressing ... #################')
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print('################# Scaling ... #################')
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print('################# Calculating tSNE ... #################')
    sc.tl.tsne(adata_filt) # Calculate tsne
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
def filter_clusters_by_gene(data, gene, threshold = 0.5):
    # Load your AnnData object
    adata = data
    sc.tl.rank_genes_groups(adata, groupby='leiden')
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

#break Execution start begin with timing
#start_time = time.time()
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

inspect_stem = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'final']
gastric_markers = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1', 'final']
duodenal_markers = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'LYZ', 'final']
inspect_prolif = ['MKI67', 'PLK1', 'E2F1', 'FOXM1', 'MCM2', 'MCM7', 'BUB1', 'CCNE1', 'CCND1', 'CCNB1', 'TOP2A']

antrum_orgs = sc.read('C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/Data annotation/Antrum_organoids/ant_organoid_annotated.h5ad')
str(antrum_orgs)
sc.pl.umap(antrum_orgs, color = gastric_markers)
sc.pl.umap(antrum_orgs, color = duodenal_markers)
sc.pl.umap(antrum_orgs, color = inspect_stem)

org_LGR5 = isolate_cells_by_gene(data = antrum_orgs, gene = 'LGR5', threshold = 0.1)
org_MKI67 = isolate_cells_by_gene(data = antrum_orgs, gene = 'MKI67', threshold = 0.1)

cont_org = antrum_orgs[antrum_orgs.obs['sample_id'] == 'Ant_organoid_Ctrl']
AGR2_org = antrum_orgs[antrum_orgs.obs['sample_id'] == 'Ant_organoid_AGR2']

sc.pl.umap(cont_org, color = ['MKI67', 'LGR5', 'OLFM4', 'SMOC2'], title = 'Control')
sc.pl.umap(AGR2_org, color = ['MKI67', 'LGR5', 'OLFM4', 'SMOC2'], title = 'Patient')

sc.pl.umap(antrum_orgs, color = inspect_prolif)

cont_org_LGR5 = isolate_cells_by_gene(data = cont_org, gene = 'LGR5', threshold = 0.1)
AGR2_org_LGR5 = isolate_cells_by_gene(data = AGR2_org, gene = 'LGR5', threshold = 0.1)

cont_org_MKI67 = isolate_cells_by_gene(data = cont_org, gene = 'MKI67', threshold = 0.1)
AGR2_org_MKI67 = isolate_cells_by_gene(data = AGR2_org, gene = 'MKI67', threshold = 0.1)

sc.pl.violin(antrum_orgs, keys = ['LGR5', 'OLFM4', 'SMOC2', 'STMN1', 'MKI67'], groupby = 'sample_id', rotation = 45)
sc.pl.violin(org_LGR5, keys = ['LGR5', 'OLFM4', 'SMOC2', 'STMN1', 'MKI67'], groupby = 'sample_id', rotation = 45)
sc.pl.violin(org_MKI67, keys = ['LGR5', 'OLFM4', 'SMOC2', 'STMN1', 'MKI67'], groupby = 'sample_id', rotation = 45)

sc.pl.scatter(antrum_orgs, x = 'LGR5', y = 'OLFM4', color = 'SMOC2', groups = 'patient_id')
sc.pl.scatter(org_LGR5, x = 'LGR5', y = 'OLFM4', color = 'SMOC2', groups = 'patient_id')
sc.pl.scatter(org_MKI67, x = 'LGR5', y = 'OLFM4', color = 'SMOC2', groups = 'patient_id')

sc.pl.scatter(adata = cont_org_LGR5, x = 'LGR5', y = 'OLFM4', color = 'SMOC2')
sc.pl.scatter(adata = cont_org_LGR5, x = 'LGR5', y = 'SMOC2', color = 'OLFM4')
sc.pl.scatter(adata = cont_org_LGR5, x = 'MKI67', y = 'SMOC2', color = 'OLFM4')

sc.pl.scatter(adata = AGR2_org_LGR5, x = 'LGR5', y = 'OLFM4', color = 'SMOC2')
sc.pl.scatter(adata = AGR2_org_LGR5, x = 'LGR5', y = 'SMOC2', color = 'OLFM4')
sc.pl.scatter(adata = AGR2_org_LGR5, x = 'MKI67', y = 'SMOC2', color = 'OLFM4')

sc.pl.umap(cont_org, color = gastric_markers, title = 'cont_org_gastric_markers')
sc.pl.umap(cont_org, color = duodenal_markers, title = 'cont_org_duod_markers')
sc.pl.umap(cont_org, color = inspect_stem, title = 'cont_org_stem_markers')

sc.pl.umap(AGR2_org, color = gastric_markers, title = 'AGR2_org_gastric_markers')
sc.pl.umap(AGR2_org, color = duodenal_markers, title = 'AGR2_org_duod_markers')
sc.pl.umap(AGR2_org, color = inspect_stem, title = 'AGR2_org_stem_markers')

#%%

stem_cluster = antrum_orgs[antrum_orgs.obs['final'] == 'Proliferative3_LGR5hi']
sc.pl.umap(stem_cluster, color = gastric_markers)
stem_cluster2 = stem_cluster.copy()
sc.pp.neighbors(stem_cluster, n_neighbors = 50) # Calculate neighbors
sc.pl.umap(stem_cluster)
