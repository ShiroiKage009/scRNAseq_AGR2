# Test to read shit

import scanpy as sc

adata = sc.read("C:/Work cache/data cache/testing/COMBAT-CITESeq-DATA.h5ad")

# Import packages
import scanpy as sc
import anndata as ad

#%% function definitions
def process_for_UMAP(adata, normed = 0, leiden_res = 0.8):
    print("Filtering ...")
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = 10) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    if not normed:
        adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 50] # Filter on the cells with fewer than 10% mitochondrial reads
    else:
        adata_filt = adata_prefilt
    print("Normalizing ...")
    sc.pp.normalize_total(adata_filt, target_sum=1e4) # Normalize
    print("Log scaling ...")
    sc.pp.log1p(adata_filt) # Log scaling
    print("Finding variable genes ...")
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    print("Saving raw data ...")
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    print("Filtering on variable genes ...")
    adata_filt = adata_filt[:, adata_filt.var.highly_variable] # Filter on genes that are highly variable
    print("Regressing ...")
    sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt']) # Regression. Not sure what that is.
    print("Scaling ...")
    sc.pp.scale(adata_filt, max_value = 10) # Scale the data
    print("Calculating PCA ...")
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("Calculating tSNE ...")
    sc.tl.tsne(adata_filt)
    print("Calculating neighbors ...")
    sc.pp.neighbors(adata_filt)
    print("Calculating Leiden ...")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("Calculating PAGA ...")
    sc.tl.paga(adata_filt)
    print("Plotting PAGA ...")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("Calculating UMAP init_pos = paga")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("Calculating UMAP ...")
    sc.tl.umap(adata_filt)
    print("Plotting UMAP ...")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def recalc_UMAP(adata_filt, leiden_res = 0.8):
    sc.tl.pca(adata_filt, svd_solver='arpack') # Compute PCA
    print("Calculating tSNE ...")
    sc.tl.tsne(adata_filt)
    print("Calculating neighbors ...")
    sc.pp.neighbors(adata_filt)
    print("Calculating Leiden ...")
    sc.tl.leiden(adata_filt, resolution = leiden_res)
    print("Calculating PAGA ...")
    sc.tl.paga(adata_filt)
    print("Plotting PAGA ...")
    sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
    print("Calculating UMAP init_pos = paga")
    sc.tl.umap(adata_filt, init_pos='paga')
    print("Calculating UMAP ...")
    sc.tl.umap(adata_filt)
    print("Plotting UMAP ...")
    sc.pl.umap(adata_filt, color = ['leiden'])
    return adata_filt

def process_until_norm(adata, cells):
    sc.pp.filter_cells(adata, min_counts = 2000) # Filter cells based on number of RNA reads
    sc.pp.filter_cells(adata, min_genes= 700) # Filter cells based on the number of recognized genes
    sc.pp.filter_genes(adata, min_cells = cells) # Filter genes based on the minimum number of cells expressing it
    adata_prefilt = adata[adata.obs['predicted_doublets'] == False]
    adata_filt = adata_prefilt[adata_prefilt.obs['pct_counts_mt'] < 10] # Filter on the cells with fewer than 10% mitochondrial reads
    sc.pp.normalize_total(adata_filt, target_sum = 1e4) # Normalize
    sc.pp.log1p(adata_filt) # Log scaling
    sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5) # Compute differentially expressed genes within the sample
    adata_filt.raw = adata_filt # Store the raw files in its own layer
    return adata_filt
    

def isolate_cells_by_gene(data, gene, threshold):
    # Now subset_ant_mt_filt contains only the highly variable genes
    data_subset = data[data[:, gene].X > threshold]
    
    return data_subset

#%% Environment settings and misc variables
sc.settings.verbosity = 3
sc.set_figure_params(dpi = 600)
plot_LGRMK = ['LGR5', 'MKI67', 'leiden', 'Localization']


#%% Read the files
test = sc.read("C:/Work cache/data cache/testing/COMBAT-CITESeq-DATA.h5ad")

#%%

process_for_UMAP(test)

recalc_UMAP(test)
