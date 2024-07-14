# Testing of annotated organoid data
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

#break Execution start begin with timing
#start_time = time.time()
#break Environemtnal settings and variable definitions:
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)

inspect_stem = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'final']
gastric_markers = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'TFF2', 'MUC5AC', 'GKN2', 'TFF1', 'GHRL', 'AQP5', 'MUC1', 'final']
duodenal_markers = ['final', 'MKI67', 'LGR5', 'TNFRSF19', 'BHLHA15', 'ASCL2', 'CD44', 'SMOC2', 'SOX2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'MUC6', 'MUC2', 'TFF3', 'GHRL', 'AQP5', 'MUC1', 'ANPEP', 'LYZ', 'final']
inspect_prolif = ['MKI67', 'PLK1', 'E2F1', 'FOXM1', 'MCM2', 'MCM7', 'BUB1', 'CCNE1', 'CCND1', 'CCNB1', 'TOP2A']
stem = ['LGR5', 'ASCL2', 'CD44', 'SMOC2', 'STMN1', 'OLFM4', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2']
gastric_stem = ['LGR5', 'ASCL2', 'CD44', 'STMN1', 'BMI1', 'LRIG1', 'BHLHA15', 'TNFRSF19', 'SOX2']

#break
# Read files
antrum_orgs = sc.read('C:/Work cache/Project sync/PhD/Research projects/AGR2 follow-up/Data cache/ssRNAseq/Aline/Data annotation/Antrum_organoids/ant_organoid_annotated.h5ad')
# Object inspection
#str(antrum_orgs)
#sc.pl.umap(antrum_orgs, color = gastric_markers)
#sc.pl.umap(antrum_orgs, color = duodenal_markers)
#sc.pl.umap(antrum_orgs, color = inspect_stem)

# Isolation of combined LGR5+ cells 
org_LGR5 = isolate_cells_by_gene(data = antrum_orgs, gene = 'LGR5', threshold = 0.1)
org_MKI67 = isolate_cells_by_gene(data = antrum_orgs, gene = 'MKI67', threshold = 0.1)

# Separation of samples
cont_org = antrum_orgs[antrum_orgs.obs['sample_id'] == 'Ant_organoid_Ctrl']
AGR2_org = antrum_orgs[antrum_orgs.obs['sample_id'] == 'Ant_organoid_AGR2']

# Calculation of gene scores for separated epithelium
sc.tl.score_genes(adata = cont_org, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = AGR2_org, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = cont_org, gene_list = ['OLFM4', 'SMOC2'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = AGR2_org, gene_list = ['OLFM4', 'SMOC2'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = cont_org, gene_list = stem, score_name = 'Stemness')
sc.tl.score_genes(adata = AGR2_org, gene_list = stem, score_name = 'Stemness')

# Plots for inspection
sc.pl.umap(cont_org, color = ['MKI67', 'LGR5', 'OLFM4', 'SMOC2'], title = 'Control')
sc.pl.umap(AGR2_org, color = ['MKI67', 'LGR5', 'OLFM4', 'SMOC2'], title = 'Patient')
sc.pl.umap(antrum_orgs, color = inspect_prolif)

# Isolation of LGR5 cells for each sample separately
cont_org_LGR5 = isolate_cells_by_gene(data = cont_org, gene = 'LGR5', threshold = 0.5)
AGR2_org_LGR5 = isolate_cells_by_gene(data = AGR2_org, gene = 'LGR5', threshold = 0.5)

# Calculation of gene scores for separated LGR5
sc.tl.score_genes(adata = cont_org_LGR5, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = AGR2_org_LGR5, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = cont_org_LGR5, gene_list = ['OLFM4', 'SMOC2'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = AGR2_org_LGR5, gene_list = ['OLFM4', 'SMOC2'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = cont_org_LGR5, gene_list = stem, score_name = 'Stemness')
sc.tl.score_genes(adata = AGR2_org_LGR5, gene_list = stem, score_name = 'Stemness')

# Isolation of proliferative cells using MKI67 from each sample
cont_org_MKI67 = isolate_cells_by_gene(data = cont_org, gene = 'MKI67', threshold = 0.1)
AGR2_org_MKI67 = isolate_cells_by_gene(data = AGR2_org, gene = 'MKI67', threshold = 0.1)

# Inspection violins
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

#%% Calculate gene scores for the integrated object now that I got what I need from it
# Calculation of gene scores for separated epithelium
sc.tl.score_genes(adata = antrum_orgs, gene_list = gastric_stem, score_name = 'Gastric stemness')
sc.tl.score_genes(adata = antrum_orgs, gene_list = ['OLFM4', 'SMOC2'], score_name = 'Intestinal stemness')
sc.tl.score_genes(adata = antrum_orgs, gene_list = stem, score_name = 'Stemness')

#%%
plot1 = sc.pl.scatter(cont_org, x = 'Gastric stemness', y = 'Intestinal stemness', size = 15, show = 0)
plot2 = sc.pl.scatter(AGR2_org, x = 'Gastric stemness', y = 'Intestinal stemness', size = 15, show = 0)
plot1.set_ylim(-1.8,3.8);
plot1.set_xlim(-1.5, 2.5);
plot2.set_ylim(-1.8,3.8);
plot2.set_xlim(-1.5, 2.5);
plot1
plot2

#%%
sc.pl.umap(cont_org, color = 'Stemness', vmax = 2, vmin = -1, size = 5)
sc.pl.umap(AGR2_org, color = 'Stemness', vmax = 2, vmin = -1, size = 5)
sc.pl.umap(cont_org, color = 'final', size = 5)
sc.pl.umap(AGR2_org, color = 'final', size = 5)
sc.pl.umap(antrum_orgs, color = 'Stemness', vmax = 2, vmin = -1, size = 5, title = 'iiiiiiiiiiiii')
sc.pl.umap(antrum_orgs, color = 'final', size = 5)

#%% Histogram for the epithelium
plt.hist(cont_org.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Control', color='blue', density = 1)
plt.hist(AGR2_org.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Patient', color='red', density = 1)

# Add titles and labels
plt.title('Antral Epithalial Intestinal Stemness')
plt.xlabel('Intestinal Stemness')
plt.ylabel('Count')
plt.xlim(-2, 3)
plt.legend(loc='upper right')

# Show plot
plt.show()

#%% Histogram for LGR5
plt.hist(cont_org_LGR5.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Control', color='blue', density = 1)
plt.hist(AGR2_org_LGR5.obs['Intestinal stemness'], bins=50, alpha=0.5, label='Patient', color='red', density = 1)

# Add titles and labels
plt.title('Antral LGR5+ Intestinal Stemness')
plt.xlabel('Intestinal Stemness')
plt.ylabel('Count')
plt.xlim(-2, 3)
plt.legend(loc='upper right')

# Show plot
plt.show()

#%% Boxplot
# Extract the "Intestinal stemness" data
intestinal_stemness_data1 = cont_org_LGR5.obs['Intestinal stemness']
intestinal_stemness_data2 = AGR2_org_LGR5.obs['Intestinal stemness']

# Create a DataFrame for plotting
data1 = pd.DataFrame({'Intestinal stemness': intestinal_stemness_data1, 'Dataset': 'Control'})
data2 = pd.DataFrame({'Intestinal stemness': intestinal_stemness_data2, 'Dataset': 'Patient'})
combined_data = pd.concat([data1, data2])

# Plot boxplots
plt.figure(figsize=(5, 10))
sns.boxplot(x='Dataset', y='Intestinal stemness', data=combined_data, palette=['blue', 'red'])

# Add titles and labels
plt.title('Distribution of Intestinal Stemness in Two Datasets')
plt.xlabel('Dataset')
plt.ylabel('Intestinal Stemness')

# Show plot
plt.show()
#%%
epith_export_cont = cont_org.obs[['Intestinal stemness']]
epith_export_cont.to_csv('cont epith intestinal stemness export.csv')
epith_export_AGR2 = AGR2_org.obs[['Intestinal stemness']]
epith_export_AGR2.to_csv('AGR2 epith intestinal stemness export.csv')

#%%
LGR5_export_cont = cont_org_LGR5.obs[['Intestinal stemness']]
LGR5_export_cont.to_csv('cont LGR5 intestinal stemness export.csv')
LGR5_export_AGR2 = AGR2_org_LGR5.obs[['Intestinal stemness']]
LGR5_export_AGR2.to_csv('AGR2 LGR5 intestinal stemness export.csv')

#%%
scatter_density_plot(data = cont_org, upxlim = 4, loxlim = -1.5, upylim = 4, loylim = -2, size = 10)
scatter_density_plot(data = AGR2_org, upxlim = 4, loxlim = -1.5, upylim = 4, loylim = -2, size = 10)

#%%
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np

# Fill in the anndata object
data = cont_org

# Extract the gene expression data for 'AGR2' and 'LGR5'
gene1 = 'Intestinal stemness'
gene2 = 'LGR5'

#gene1_data = data[:, gene1].X.flatten() # This is for when it's a gene
gene1_data = data.obs[gene1].values # This is for when it's a score

gene2_data = data[:, gene2].X.flatten()

# Calculate the correlation coefficient
correlation_coefficient = np.corrcoef(gene1_data, gene2_data)[0, 1]

# Create a scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(gene1_data, gene2_data, c='blue', s=15, edgecolor='none')

# Set the plot title
plt.title(f"Correlation between {gene1} and {gene2}\nCorrelation: {correlation_coefficient:.2f}")

# Label the axes
plt.xlabel(gene1)
plt.ylabel(gene2)

# Show the plot
plt.show()
