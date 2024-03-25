# This script tries to isolate the LGR5 cells and get their barcodes from the annotated and QC'ed files.
# "filtering and clustering analysis.py" has a test cell with identifying cells from raw on my own then retrieving the 
# cells from the annotated files and reclustering.

#%% Library imports
import scanpy as sc
import pandas as pd

#%% function definitions
def process_for_UMAP(data, normed = 0, leiden_res = 0.8):
    adata = data.copy() # This is to avoid writing into the file that's entered as an argument
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

def recalc_UMAP(data, leiden_res = 0.8):
    adata_filt = data
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
    adata = data.copy() # This is to avoid writing into the file that's entered as an argument
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
    print("################# Finding variable genes ... #################")
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
    adata = data.copy()
    sc.tl.rank_genes_groups(adata, groupby='leiden', method = 'wilcoxon')
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

#%% Setting up environmental variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)
inspect_stem = ['LGR5', 'MKI67', 'TNFRSF19', 'BMI1', 'LRIG1', 'Patient']
global_res = 0.5
LGR5_threshold = 0.5
diff_exp_method = 'wilcoxon'

#%% Reading the files
antrum_annot = sc.read('C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/Data annotation/antrum/antrum_annotated_final1122.h5ad')

#%% Plotting for QC
sc.pl.umap(antrum_annot, color = inspect_stem)
antrum_annot_LGR5 = isolate_cells_by_gene(data = antrum_annot, gene = 'LGR5', threshold = 0.5)
sc.pl.umap(antrum_annot_LGR5, color = 'Patient', size = 120)
print(antrum_annot.obs['Patient'].value_counts())

#%% Isolate the LGR5 cells by barcode from the raw file
antrum_annot_LGR5_barcodes = antrum_annot_LGR5.obs_names.tolist()
ant_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_antrum.h5ad")
antrum_unfilt_LGR5 = ant_unfilt[antrum_annot_LGR5_barcodes].copy()

#%% Process and plot the isolated cells
antrum_LGR5_calc = process_for_UMAP(antrum_unfilt_LGR5, leiden_res = global_res)
sc.pl.umap(antrum_LGR5_calc, color = inspect_stem)

