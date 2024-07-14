# Currently, this file is being used to evaluate the LGR5 niche specifically for transformation alongside other stem cell niches
# Import packages
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import time

def time_it(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper

# DEFAULT QC VALUES
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
    print('################# Regressing ... #################')
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print('################# Scaling ... #################')
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print('################# Calculating PCA ... #################')
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print('################# Calculating tSNE ... #################')
    sc.tl.tsne(adata_filt) # Calculate tsne
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

def scatter_density_plot(data, background = 'grey', cmap = 'hsv', upxlim = None, loxlim = None, upylim = None, loylim = None, size = None, bw_adjust = 0.5):
    # Assuming 'cont_org' is your AnnData object and it has the necessary columns
    x = 'Gastric stemness'
    y = 'Intestinal stemness'
    
    # Extract the data from the AnnData object
    x_data = data.obs[x].values
    y_data = data.obs[y].values
    
    # Creating the density plot
    plt.figure()
    
    # Set the axes
    plt.xlim(loxlim, upxlim);
    plt.ylim(loylim, upylim);
    
    # Create the base scatter plot
    plt.scatter(x_data, y_data, c = 'gray', s = size, edgecolor='none')

    # Create the heatmap on top    
    sns.kdeplot(x = x_data, y = y_data, cmap = "hsv", shade=True, bw_adjust = bw_adjust)
    
    # Label the axes
    plt.xlabel(x)
    plt.ylabel(y)
    
    # Show the plot
    plt.show()
    
    return(plt)

#######################################################
################## FUNCTION DEF END ###################
#######################################################

# Function for mapping names to columns specifically
def map_to_column(data, map_set, column = 'Localization'):
    data.obs[column + '_old'] = data.obs[column]
    data.obs[column] = data.obs[column].map(map_set)
    print(data.obs[column])
    print(data.obs[column + '_old'])
    return 'Mapping function done'

#######################################################
################## FUNCTION DEF END ###################
#######################################################

#break Execution start begin with timing
start_time = time.time()
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

# MIK67 = Ki67, TNSFRSF19 = TROY
stem = ['LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2']
gastric_stem = ['LGR5', 'ASCL2', 'CD44', 'STMN1', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2']
inspect_stem = ['leiden', 'MKI67', 'LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2', 'leiden', 'Localization']
inspect_stem_duo = ['leiden', 'MKI67', 'LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'leiden', 'Localization']

inspect_stem_gast = ['leiden', 'MKI67', 'LGR5', 'TNFRSF19', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'leiden', 'Localization']

gastric_markers = ['leiden', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1', 'leiden', 'Localization']
duodenal_markers = ['leiden', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'LYZ', 'leiden', 'Localization']

global_res = 0.5


#Break
duo_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_duodenum.h5ad")
ant_unfilt = sc.read("C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/raw_data/agr2_unfilt_antrum.h5ad")
duo_unfilt.obs['Localization'] = duo_unfilt.obs['Site'].astype(str) + ' ' + duo_unfilt.obs['Patient'].astype(str)
ant_unfilt.obs['Localization'] = ant_unfilt.obs['Site'].astype(str) + ' ' + ant_unfilt.obs['Patient'].astype(str)
combo_unfilt = sc.concat(adatas = [duo_unfilt, ant_unfilt], join = 'outer')

ant_patient = ant_unfilt[ant_unfilt.obs['Patient'] == 'P26']
ant_cont = ant_unfilt[ant_unfilt.obs['Patient'] == 'GI6253']
duo_cont = duo_unfilt[duo_unfilt.obs['Patient'] == 'GI6253']
duo_pat = duo_unfilt[duo_unfilt.obs['Patient'] == 'P26']

ant_pat_proc = process_for_UMAP(data = ant_patient)
ant_cont_proc = process_for_UMAP(data = ant_cont)
duo_cont_proc = process_for_UMAP(data = duo_cont)
duo_pat_proc = process_for_UMAP(data = duo_pat)

ant_pat_proc_filt = filter_clusters_by_gene(data = ant_pat_proc, gene = 'EPCAM')
ant_cont_proc_filt = filter_clusters_by_gene(data = ant_cont_proc, gene = 'EPCAM')
duo_cont_proc_filt = filter_clusters_by_gene(data = duo_cont_proc, gene = 'EPCAM')
duo_pat_proc_filt = filter_clusters_by_gene(data = duo_pat_proc, gene = 'EPCAM')

pat_epith_barcode = ant_pat_proc_filt.obs_names.tolist()
cont_epith_barcode = ant_cont_proc_filt.obs_names.tolist()
duocont_epith_barcode = duo_cont_proc_filt.obs_names.tolist()
duopat_epith_barcode = duo_pat_proc_filt.obs_names.tolist()

ant_pat_refilt = ant_unfilt[pat_epith_barcode]
ant_cont_refilt = ant_unfilt[cont_epith_barcode]
duo_cont_refilt = duo_unfilt[duocont_epith_barcode]
duo_pat_refilt = duo_unfilt[duopat_epith_barcode]
combo_refilt = sc.concat(adatas = [ant_pat_refilt, ant_cont_refilt, duo_cont_refilt], join = 'outer')

pat_reproc = process_for_UMAP(data = ant_pat_refilt, leiden_res = 0.5)
cont_reproc = process_for_UMAP(data = ant_cont_refilt, leiden_res = 0.5)
duocont_reproc = process_for_UMAP(data = duo_cont_refilt, leiden_res = 0.5)
duopat_reproc = process_for_UMAP(data = duo_pat_refilt, leiden_res = 0.5)
combo_reproc = process_for_UMAP(data = combo_refilt, leiden_res = 0.8)

sc.tl.score_genes(adata = pat_reproc, gene_list = ['SMOC2', 'OLFM4'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = cont_reproc, gene_list = ['SMOC2', 'OLFM4'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = combo_reproc, gene_list = ['SMOC2', 'OLFM4'], score_name = 'Intestinal stemness')

sc.tl.score_genes(adata = pat_reproc, gene_list = stem, score_name = 'Stemness')
sc.tl.score_genes(adata = cont_reproc, gene_list = stem, score_name = 'Stemness')
sc.tl.score_genes(adata = combo_reproc, gene_list = stem, score_name = 'Stemness')

sc.tl.score_genes(adata = pat_reproc, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = cont_reproc, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = combo_reproc, gene_list = gastric_stem, score_name = 'Gastric stemness')

plotA = sc.pl.scatter(pat_reproc, x = 'Stemness', y = 'Intestinal stemness', color = 'LGR5', show = False)
plotB = sc.pl.scatter(cont_reproc, x = 'Stemness', y = 'Intestinal stemness', color = 'LGR5', show = False)
plotA.set_ylim(-1.25,4)
plotB.set_ylim(-1.25,4)
plotA
plotB

plotC = sc.pl.scatter(pat_reproc, x = 'Gastric stemness', y = 'Intestinal stemness', color = 'LGR5', show = False)
plotD = sc.pl.scatter(cont_reproc, x = 'Gastric stemness', y = 'Intestinal stemness', color = 'LGR5', show = False)
plotC.set_ylim(-1.25,4)
plotD.set_ylim(-1.25,4)
plotC
plotD

plotE = sc.pl.scatter(pat_reproc, x = 'LGR5', y = 'Intestinal stemness', color = 'Gastric stemness', show = False)
plotF = sc.pl.scatter(cont_reproc, x = 'LGR5', y = 'Intestinal stemness', color = 'Gastric stemness', show = False)
plotE.set_ylim(-1.25,4)
plotF.set_ylim(-1.25,4)
plotE
plotF

end_time = time.time()
print("Script executed in", end_time - start_time, "seconds")
print("Script executed in", (end_time - start_time)/60, "minutes")
#%%

sc.pl.umap(pat_reproc, color = ['leiden', 'MKI67', 'LGR5', 'BMI1', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'TNFRSF19', 'STMN1', 'OLFM4', 'BHLHA15', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1'], title = 'Patient Leiden', save = 'patient antrum.png')
sc.pl.umap(cont_reproc, color = ['leiden', 'MKI67', 'LGR5', 'BMI1', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'TNFRSF19', 'STMN1', 'OLFM4', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1'], title = 'Control Leiden', save = 'control antrum.png')
sc.pl.umap(duocont_reproc, color = ['leiden', 'MKI67', 'LGR5', 'BMI1', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'LYZ'], title = 'Control duo Leiden', save = 'control duodenum.png')
sc.pl.umap(duopat_reproc, color = ['leiden', 'MKI67', 'LGR5', 'BMI1', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'MUC1', 'ANPEP', 'LYZ'], title = 'Patient duo Leiden', save = 'patient duodenum.png')
sc.pl.umap(combo_reproc, color = ['leiden', 'MKI67', 'LGR5', 'BMI1', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'TNFRSF19', 'STMN1', 'OLFM4', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'MUC2', 'TFF3', 'MUC5B', 'LYZ', 'Localization'], title = 'Combo lieden', size = 10, save = 'combination.png')

sc.pl.umap(pat_reproc, color = inspect_stem, title = 'patient_ant_stem leiden', size = 20, save = 'patient_ant_stem.png')
sc.pl.umap(cont_reproc, color = inspect_stem_gast, title = 'control_ant_stem leiden', size = 20, save = 'control_ant_stem.png')
sc.pl.umap(duocont_reproc, color = inspect_stem_duo, title = 'control_duo_inspect_stem leiden', size = 20, save = 'control_duo_inspect_stem.png')
sc.pl.umap(duopat_reproc, color = inspect_stem_duo, title = 'patient_duo_inspect_stem leiden', size = 20, save = 'patient_duo_inspect_stem.png')


sc.tl.rank_genes_groups(adata = pat_reproc, groupby = 'leiden', method = 'wilcoxon')
sc.tl.rank_genes_groups(adata = cont_reproc, groupby = 'leiden', method = 'wilcoxon')
sc.tl.rank_genes_groups(adata = duocont_reproc, groupby = 'leiden', method = 'wilcoxon')
sc.tl.rank_genes_groups(adata = duopat_reproc, groupby = 'leiden', method = 'wilcoxon')
sc.tl.rank_genes_groups(adata = combo_reproc, groupby = 'leiden', method = 'wilcoxon')

sc.pl.rank_genes_groups(adata = cont_reproc, n_genes = 25, save = 'control antrum DGE.png')
sc.pl.rank_genes_groups(adata = pat_reproc, n_genes = 25, save = 'patient antrum DGE.png')
sc.pl.rank_genes_groups(adata = duocont_reproc, n_genes = 25, save = 'control duodenum DGE.png')
sc.pl.rank_genes_groups(adata = duopat_reproc, n_genes = 25, save = 'patient duodenum DGE.png')
sc.pl.rank_genes_groups(adata = combo_reproc, n_genes = 25, save = 'combination DGE.png')

end_time = time.time()
print("Script executed in", end_time - start_time, "seconds")
print("Script executed in", (end_time - start_time)/60, "minutes")

#%%
#sc.pl.rank_genes_groups_dotplot(adata = cont_reproc, n_genes = 10, groupby = 'leiden')
#sc.pl.rank_genes_groups_dotplot(adata = pat_reproc, n_genes = 10)

#%%
LGR5_cont_ant = isolate_cells_by_gene(data = cont_reproc, gene = 'LGR5', threshold = 0.1)
LGR5_pat_ant = isolate_cells_by_gene(data = pat_reproc, gene = 'LGR5', threshold = 0.1)
LGR5_cont_duo = isolate_cells_by_gene(data = duocont_reproc, gene = 'LGR5', threshold = 0.1)
LGR5_pat_duo = isolate_cells_by_gene(data = duopat_reproc, gene = 'LGR5', threshold = 0.1)

#%%
sc.pl.scatter(pat_reproc, x = 'LGR5', y = 'BMI1', color = 'MKI67')
sc.pl.umap(pat_reproc, color = 'GAPDH', use_raw = 1)

sc.pl.scatter(cont_reproc, x = 'LGR5', y = 'OLFM4', color = 'SMOC2')
sc.pl.scatter(pat_reproc, x = 'LGR5', y = 'OLFM4', color = 'SMOC2')

combonation = sc.concat(adatas = [cont_reproc, pat_reproc], join = 'outer')
crazy_combo = sc.concat(adatas = [cont_reproc, pat_reproc, duocont_reproc], join = 'outer')
combonation_LGR5 = isolate_cells_by_gene(data = combonation, gene = 'LGR5', threshold = 0.1)
sc.pl.violin(adata = combonation_LGR5, keys = 'LGR5', groupby = 'Patient')
sc.pl.violin(adata = combonation, keys = 'LGR5', groupby = 'Patient')

crazy_combo_LGR5 = isolate_cells_by_gene(data = crazy_combo, gene = 'LGR5', threshold = 0.1)
sc.pl.violin(adata = crazy_combo, keys = 'LGR5', groupby = 'Localization')
sc.pl.violin(adata = crazy_combo_LGR5, keys = 'LGR5', groupby = 'Localization')
sc.pl.violin(adata = crazy_combo, keys = 'OLFM4', groupby = 'Localization')

crazy_combo_ASCL2 = isolate_cells_by_gene(data = crazy_combo, gene = 'ASCL2', threshold = 0.2)
sc.pl.violin(adata = crazy_combo_ASCL2, keys = 'ASCL2', groupby = 'Localization')

sc.pl.scatter(pat_reproc, x = 'ASCL2', y = 'SMOC2', color = 'LGR5', title = 'LGR5 Patient Antrum')
sc.pl.scatter(cont_reproc, x = 'ASCL2', y = 'SMOC2', color = 'LGR5', title = 'LGR5 control antrum')
sc.pl.scatter(duocont_reproc, x = 'ASCL2', y = 'SMOC2', color = 'LGR5', title = 'LGR5 duodenum control')
sc.pl.scatter(crazy_combo, x = 'ASCL2', y = 'SMOC2', color = 'Localization')

sc.pl.scatter(pat_reproc, x = 'OLFM4', y = 'SMOC2', color = 'LGR5', title = 'LGR5 Patient Antrum')
sc.pl.scatter(cont_reproc, x = 'OLFM4', y = 'SMOC2', color = 'LGR5', title = 'LGR5 control antrum')
sc.pl.scatter(duocont_reproc, x = 'OLFM4', y = 'SMOC2', color = 'LGR5', title = 'LGR5 duodenum control')
sc.pl.scatter(crazy_combo, x = 'OLFM4', y = 'SMOC2', color = 'Localization')

combonation_SMOC2 = isolate_cells_by_gene(data = combonation, gene = 'SMOC2', threshold = 0.1)

#%% violin plots for comparison
LGR5contant = isolate_cells_by_gene(data = cont_reproc, gene = 'LGR5', threshold = 0.1)
LGR5patant = isolate_cells_by_gene(data = pat_reproc, gene = 'LGR5', threshold = 0.1)
#SMOC2contant = isolate_cells_by_gene(data = cont_reproc, gene = 'SMOC2', threshold = 0)
SMOC2patant = isolate_cells_by_gene(data = pat_reproc, gene = 'SMOC2', threshold = 0)

sc.pl.violin(adata = LGR5contant, keys = 'LGR5')
sc.pl.violin(adata = LGR5patant, keys = 'LGR5')

LGR5_combo = sc.concat(adatas = [LGR5contant, LGR5patant], join = 'outer')
sc.pl.violin(adata = LGR5_combo, keys = 'LGR5', groupby = 'Localization')
LGR5_barcode = LGR5_combo.obs_names.tolist()
LGR5_refilt = combo_unfilt[LGR5_barcode]
LGR5_reproc = process_for_UMAP(LGR5_refilt, leiden_res = 0.2)
sc.pl.umap(LGR5_reproc, color = ['Patient', 'leiden'])
sc.pl.violin(adata = LGR5_reproc, keys = 'LGR5', groupby = 'Localization')
sc.pl.violin(adata = combonation, keys = 'SMOC2', groupby = 'Localization')
sc.pl.violin(adata = combonation, keys = 'OLFM4', groupby = 'Localization')

plot1 = sc.pl.scatter(cont_reproc, x = 'LGR5', y = 'SMOC2', color = 'OLFM4', groups = 'Localization', show = False)
plot2 = sc.pl.scatter(pat_reproc, x = 'LGR5', y = 'SMOC2', color = 'OLFM4', groups = 'Localization', show = False)
plot1.set_ylim(0,3.5);
plot2.set_ylim(0,3.5);

plot1
plot2

#%%
sc.pl.umap(pat_reproc, color = inspect_stem)
sc.pl.umap(pat_reproc, color = duodenal_markers)
sc.pl.umap(pat_reproc, color = gastric_markers)

#%%
plt.hist(cont_reproc.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Control', color='blue', density = 1)
plt.hist(pat_reproc.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Patient', color='red', density = 1)

# Add titles and labels
plt.title('Antral Epithalial Intestinal Stemness')
plt.xlabel('Intestinal Stemness')
plt.ylabel('Count')
plt.xlim(-2, 3)
plt.legend(loc='upper right')

# Show plot
plt.show()
#%%
plt.hist(LGR5_cont_ant.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Control', color='blue')
plt.hist(LGR5_pat_ant.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Patient', color='red')

# Add titles and labels
plt.title('Antral Epithalial Intestinal Stemness')
plt.xlabel('Intestinal Stemness')
plt.ylabel('Count')
plt.xlim(-2, 3)
plt.legend(loc='upper right')

# Show plot
plt.show()


#%%
scatter_density_plot(data = cont_reproc, upxlim = 1, loxlim = -0.5, upylim = 4.5, loylim = -1.5, size = 10)
scatter_density_plot(data = pat_reproc, upxlim = 1, loxlim = -0.5, upylim = 4.5, loylim = -1.5, size = 10)
plot1 = sc.pl.scatter(pat_reproc, y = 'Intestinal stemness', x = 'Gastric stemness', color = 'LGR5', show = 0, size = 50)
plot1.set_xlim(-0.5, 1);
plot1.set_ylim(-1.5, 4.5);
plot1


#%%

# Extract the relevant arrays
gene_names = cont_reproc.uns['rank_genes_groups']['names']
logfoldchanges = cont_reproc.uns['rank_genes_groups']['logfoldchanges']
pvals_adj = cont_reproc.uns['rank_genes_groups']['pvals_adj']

# Since you have two groups, you can loop or manually index each group
# Assuming group names are '0' and '1', adjust accordingly

dataframes = {}
for group in ['0', '1', '2', '3', '4', '5', '6', '7']:
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
control_ant = dataframes['2']

control_ant.to_csv('C:/Work cache/py_projs/scRNAseq_AGR2/project data cache/testing integration with separation and the stem cells part 2/saved files/volcano plot/control_cluster_2.csv')

#%%

leiden_map = {
    '0' : 'Metaplastic antrum',
    '1' : 'Gastric antrum'}

map_to_column(data = LGR5_reproc, map_set = leiden_map, column = 'leiden')

#%%
sc.tl.rank_genes_groups(adata = LGR5_reproc, groupby = 'leiden', method = 'wilcoxon')
ranked_genes = LGR5_reproc.uns['rank_genes_groups']

# Extract the top 20 genes for the first group
top_n = 20
groups = ranked_genes['names'].dtype.names  # Get the cluster names

group1 = groups[0]  # Get the name of the first group
group2 = groups[1]  # Get the name of the second group

top_genes_group1 = ranked_genes['names'][group1][:top_n]
top_genes_group2 = ranked_genes['names'][group2][:top_n]

# Combine the top genes in the specified order
ordered_genes = list(top_genes_group2) + list(reversed(top_genes_group1))

# Reverse the order of groups on the Y axis
reversed_groups = list(reversed(LGR5_reproc.obs['leiden'].cat.categories))

# Plot the combined top genes as a dotplot with reversed Y-axis order
sc.pl.dotplot(LGR5_reproc, ordered_genes, groupby='leiden', standard_scale='var', categories_order=reversed_groups)



#%%
end_time = time.time()
print("Script executed in", end_time - start_time, "seconds")
print("Script executed in", (end_time - start_time)/60, "minutes")
