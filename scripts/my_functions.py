# Function definitions


import scanpy as sc

def process_for_UMAP(adata):
        sc.pp.filter_cells(adata, min_counts = 2000)
        sc.pp.filter_cells(adata, min_genes= 700) # best practices says 700 but we're keeping it at this for now
        sc.pp.filter_genes(adata, min_cells = 50)
        adata_filt = adata[adata.obs['pct_counts_mt'] > 10]
        sc.pp.normalize_total(adata_filt, target_sum=1e4)
        sc.pp.log1p(adata_filt)
        sc.pp.highly_variable_genes(adata_filt, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
        sc.pl.highly_variable_genes(adata_filt)
        adata_filt.raw = adata_filt
        adata_filt = adata_filt[:, adata_filt.var.highly_variable]
        sc.pp.regress_out(adata_filt, ['total_counts', 'pct_counts_mt'])
        sc.pp.scale(adata_filt, max_value=10)
        sc.tl.pca(adata_filt, svd_solver='arpack')
        sc.tl.tsne(adata_filt)
        sc.pl.pca_variance_ratio(adata_filt, log=True)
        sc.pp.neighbors(adata_filt)
        sc.tl.leiden(adata_filt, resolution = 0.8)
        sc.tl.paga(adata_filt)
        sc.pl.paga(adata_filt, plot = 1)  # remove `plot=False` if you want to see the coarse-grained graph
        sc.tl.umap(adata_filt, init_pos='paga')
        sc.tl.umap(adata_filt)
        sc.pl.umap(adata_filt, color = ['leiden'])
        return adata_filt
    

col_org_unfilt = sc.read("C:/Work cache/data cache/Aline data/Aline/raw_data/agr2colon_organoids_unfilt.h5ad")

col_org_filt = process_for_UMAP(col_org_unfilt)
