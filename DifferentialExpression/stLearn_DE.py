# Note: This commands following by stLearn_Clustering.py in CellTypeAnnotation folder
# You can replace adata by adata_A, adata_B, adata_C or adata_D

sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)