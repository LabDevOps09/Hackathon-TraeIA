import scanpy as sc

# Carregar o dataset .h5ad
adata = sc.read_h5ad("dados_exemplo.h5ad")

# Qualidade: filtragem de células e genes
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Normalização
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identificar genes variáveis
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# PCA + UMAP
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Clusterização
sc.tl.leiden(adata)
