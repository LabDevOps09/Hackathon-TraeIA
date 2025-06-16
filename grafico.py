sc.pl.umap(adata, color=['leiden'], save="_clusters.png")
sc.pl.rank_genes_groups_heatmap(adata, groupby='leiden', save="_heatmap.png")
